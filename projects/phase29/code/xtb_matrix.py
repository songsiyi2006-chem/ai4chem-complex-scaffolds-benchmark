"""Frozen 144-SP molecular xTB method/environment sensitivity matrix.

Reported xTB total-state differences at fixed nuclei, with independently equilibrated
ALPB states. No reaction labels, electrode potentials, barriers or regioselectivity.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime,timezone
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

ROOT=Path(__file__).resolve().parents[1]
HARTREE_EV=27.211386245988
BOUNDARY='EXECUTED_ASSUMED_MOLECULAR_STATES: fixed single cation geometry for both states; gas electronic or independently equilibrated ALPB total-state differences (including model solvation terms). No electron/electrode reference, thermal correction, reaction path, interface, concentration or observed protonation. Not redox potentials, experimental electron affinities, chemical selectivity or uncertainty intervals.'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def js(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf8')
def write_csv(path,rows):
    if not rows:return
    with path.open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def parse_output(stdout,stderr,returncode):
    energy=re.findall(r'TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh',stdout)
    joined=stdout+'\n'+stderr
    warnings=[line.strip() for line in joined.splitlines() if 'IEEE_' in line or 'WARNING' in line.upper()]
    scc='convergence criteria satisfied' in stdout
    normal='normal termination of xtb' in joined.lower()
    success=returncode==0 and bool(energy) and scc and normal
    return {'status':'success_with_warnings' if success and warnings else 'success' if success else 'failed_or_unverified',
            'scc_converged':scc,'normal_termination':normal,'energy_hartree':float(energy[-1]) if success else None,
            'numerical_warnings':warnings,'warning_count':len(warnings),'returncode':returncode}

def state_difference(cation,radical):
    for key in ('molecule_id','geometry_sha256','gfn','environment'):
        if cation[key]!=radical[key]:raise ValueError('Incomparable states: '+key)
    if (cation['charge'],cation['uhf'],radical['charge'],radical['uhf'])!=(1,0,0,1):raise ValueError('Wrong charge/spin')
    if not all(j['status'].startswith('success') and j['energy_hartree'] is not None for j in (cation,radical)):
        raise ValueError('Both successful state energies required')
    return (radical['energy_hartree']-cation['energy_hartree'])*HARTREE_EV

def ordering_reversals(a,b,tolerance=1e-6):
    common=sorted(set(a)&set(b));reversed_pairs=[];ties=[]
    for i,j in itertools.combinations(common,2):
        x=a[i]-a[j];y=b[i]-b[j]
        if abs(x)<=tolerance or abs(y)<=tolerance:ties.append([i,j])
        elif x*y<0:reversed_pairs.append([i,j])
    return {'n_common':len(common),'n_pairwise':len(common)*(len(common)-1)//2,
            'rank_reversals':len(reversed_pairs),'reversed_pairs':reversed_pairs,'ties_within_1e-6_eV':ties}

def prepare_geometry(item,out,cfg):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol=Chem.MolFromSmiles(item['neutral_mapped_smiles'])
    ringn=[a for a in mol.GetAtoms() if a.GetAtomMapNum()==1]
    if len(ringn)!=1 or ringn[0].GetAtomicNum()!=7:raise ValueError('Ring N identity missing')
    ringn[0].SetFormalCharge(1);ringn[0].SetNumExplicitHs(1)
    for atom in mol.GetAtoms():atom.SetAtomMapNum(0)
    Chem.SanitizeMol(mol);smiles=Chem.MolToSmiles(mol);mol=Chem.AddHs(mol)
    nuclear=sum(a.GetAtomicNum() for a in mol.GetAtoms())
    if Chem.GetFormalCharge(mol)!=1 or (nuclear-1)%2:raise ValueError('Unexpected charge/electron parity')
    xyz=out/'inputs'/f"{item['id']}.xyz";xyz.parent.mkdir(parents=True,exist_ok=True)
    status=None
    if item['geometry_policy']=='exact_original_preflight':
        src=ROOT/item['geometry_source']
        if sha(src)!=item['geometry_sha256']:raise ValueError('Original geometry hash changed')
        shutil.copyfile(src,xyz)
        elements=[line.split()[0] for line in xyz.read_text().splitlines()[2:] if line.strip()]
        if sorted(elements)!=sorted(a.GetSymbol() for a in mol.GetAtoms()):raise ValueError('Reused geometry formula mismatch')
        pre=json.loads((ROOT/'results/quantum_preflight/manifest.json').read_text(encoding='utf8'))
        original=next(r for r in pre['panel'] if r['id']==item['original_id'])
        if Chem.MolToSmiles(Chem.MolFromSmiles(original['cation_smiles']))!=smiles:raise ValueError('Reused structure identity mismatch')
        status=original['mmff_status']
    else:
        params=AllChem.ETKDGv3();params.randomSeed=item['geometry_seed'];params.numThreads=cfg['threads']
        if AllChem.EmbedMolecule(mol,params)!=0 or not AllChem.MMFFHasAllMoleculeParams(mol):raise ValueError('Geometry embedding/parameters failed')
        status=AllChem.MMFFOptimizeMolecule(mol,mmffVariant=cfg['MMFF_variant'],maxIters=cfg['MMFF_max_iters'])
        if status!=0:raise ValueError('Single MMFF geometry not converged')
        Chem.MolToXYZFile(mol,str(xyz))
    return xyz,{'id':item['id'],'name':item['name'],'cation_smiles':smiles,'geometry_policy':item['geometry_policy'],
                'geometry_sha256':sha(xyz),'atom_count':mol.GetNumAtoms(),'nuclear_charge_sum':nuclear,
                'cation_electrons':nuclear-1,'radical_electrons':nuclear,'mmff_status':status,
                'geometry_source':item.get('geometry_source','generated_ETKDGv3_MMFF94')}

def summarize(out,manifest):
    import numpy as np
    from scipy.stats import spearmanr
    pairs=[]
    for item in manifest['panel']:
        for method in manifest['config']['methods']:
            for environment in manifest['config']['environments']:
                jobs=[j for j in manifest['jobs'] if j['molecule_id']==item['id'] and j['gfn']==method and j['environment']==environment]
                a=next((j for j in jobs if j['charge']==1),None);b=next((j for j in jobs if j['charge']==0),None)
                if a and b and all(j['status'].startswith('success') for j in (a,b)):
                    pairs.append({'molecule_id':item['id'],'name':item['name'],'gfn':method,'environment':environment,
                                  'geometry_sha256':item['geometry_sha256'],'delta_model_eV':state_difference(a,b),
                                  'evidence_layer':'executed_assumed_molecular_states','any_numerical_warning':bool(a['warning_count'] or b['warning_count'])})
    write_csv(out/'state_differences.csv',pairs)
    lookup={(r['molecule_id'],r['gfn'],r['environment']):r['delta_model_eV'] for r in pairs}
    deltas=[]
    for item in manifest['panel']:
        for environment in manifest['config']['environments']:
            a=lookup.get((item['id'],1,environment));b=lookup.get((item['id'],2,environment))
            if a is not None and b is not None:deltas.append({'molecule_id':item['id'],'environment':environment,'GFN1_eV':a,'GFN2_eV':b,'GFN2_minus_GFN1_eV':b-a,'absolute_method_difference_eV':abs(b-a)})
    write_csv(out/'method_sensitivity.csv',deltas)
    shifts=[]
    for item in manifest['panel']:
        for method in (1,2):
            gas=lookup.get((item['id'],method,'gas'))
            for environment in ('acetonitrile','dmf'):
                value=lookup.get((item['id'],method,environment))
                if gas is not None and value is not None:shifts.append({'molecule_id':item['id'],'gfn':method,'environment':environment,'gas_difference_eV':gas,'ALPB_difference_eV':value,'ALPB_minus_gas_eV':value-gas})
    write_csv(out/'environment_sensitivity.csv',shifts)
    rankings=[];groups=[(g,e) for g in (1,2) for e in ('gas','acetonitrile','dmf')]
    for (ga,ea),(gb,eb) in itertools.combinations(groups,2):
        a={i:lookup[(i,ga,ea)] for i in [x['id'] for x in manifest['panel']] if (i,ga,ea) in lookup}
        b={i:lookup[(i,gb,eb)] for i in [x['id'] for x in manifest['panel']] if (i,gb,eb) in lookup}
        common=sorted(set(a)&set(b));item={'first':f'GFN{ga}_{ea}','second':f'GFN{gb}_{eb}',**ordering_reversals(a,b)}
        item['spearman_rho']=float(spearmanr([a[i] for i in common],[b[i] for i in common]).statistic) if len(common)>2 else None
        rankings.append(item)
    previous=json.loads((ROOT/'results/quantum_preflight/manifest.json').read_text(encoding='utf8'))
    reproduction=[]
    for pid,qid in [('P01','Q01'),('P04','Q02'),('P09','Q03')]:
        for method in (1,2):
            old=[j for j in previous['jobs'] if j['molecule_id']==qid and j['gfn']==method]
            oldvalue=(next(j['energy_hartree'] for j in old if j['charge']==0)-next(j['energy_hartree'] for j in old if j['charge']==1))*HARTREE_EV
            new=lookup.get((pid,method,'acetonitrile'))
            reproduction.append({'molecule_id':pid,'original_id':qid,'gfn':method,'previous_eV':oldvalue,'matrix_eV':new,'absolute_difference_eV':abs(new-oldvalue) if new is not None else None})
    summary={'boundary':BOUNDARY,'complete_pairs':len(pairs),'expected_pairs':72,
             'success_jobs':sum(j['status'].startswith('success') for j in manifest['jobs']),'total_jobs':len(manifest['jobs']),
             'jobs_with_warnings':sum(bool(j.get('warning_count')) for j in manifest['jobs']),
             'method_difference_ranges_eV':{e:{'min':min(r['absolute_method_difference_eV'] for r in deltas if r['environment']==e),'max':max(r['absolute_method_difference_eV'] for r in deltas if r['environment']==e)} for e in ('gas','acetonitrile','dmf') if any(r['environment']==e for r in deltas)},
             'solvent_shift_range_eV':{'min':min(r['ALPB_minus_gas_eV'] for r in shifts),'max':max(r['ALPB_minus_gas_eV'] for r in shifts)} if shifts else None,
             'rank_comparisons':rankings,'original_preflight_reproduction':reproduction,
             'rank_meaning':'Ordering of computed model-state differences, not chemistry success or site preference',
             'no_experimental_or_potential_validation':True}
    js(out/'analysis.json',summary)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(14,5),sharey=True)
    ids=[p['id'] for p in manifest['panel']]
    for ax,e in zip(axes,('gas','acetonitrile','dmf')):
        for method in (1,2):ax.plot(ids,[lookup.get((i,method,e),math.nan) for i in ids],marker='o',label=f'GFN{method}')
        ax.set_title(e);ax.tick_params(axis='x',rotation=90);ax.set_xlabel('Fixed proposed substrate ID');ax.legend()
    axes[0].set_ylabel('E_model(radical) - E_model(cation) / eV')
    fig.suptitle('Executed molecular sensitivity matrix; no electrochemical or site prediction')
    fig.tight_layout();fig.savefig(out/'method_environment_matrix.png',dpi=160);plt.close(fig)
    return summary

def run(args):
    import rdkit
    cfg=json.loads((ROOT/'configs/xtb_matrix.json').read_text(encoding='utf8'))
    if sha(ROOT/'configs/preregistered_panel.json')!=cfg['source_panel_sha256']:raise ValueError('Source panel changed')
    n=len(cfg['species'])*len(cfg['methods'])*len(cfg['environments'])*len(cfg['states'])
    if n!=cfg['planned_single_points'] or n!=144:raise ValueError('Frozen matrix must have 144 jobs')
    executable=shutil.which(args.xtb)
    if not executable:raise FileNotFoundError('Use existing xTB; no automatic install')
    out=Path(args.out).resolve()
    if out.exists():raise ValueError('Use a new output directory; existing evidence is preserved')
    out.mkdir(parents=True)
    env=dict(os.environ,OMP_NUM_THREADS=str(cfg['threads']),MKL_NUM_THREADS=str(cfg['threads']),OPENBLAS_NUM_THREADS=str(cfg['threads']))
    version=subprocess.run([executable,'--version'],capture_output=True,text=True,errors='replace',env=env,timeout=20)
    (out/'xtb_version.txt').write_text(version.stdout+version.stderr,encoding='utf8')
    manifest={'started_utc':datetime.now(timezone.utc).isoformat(),'config':cfg,'config_sha256':sha(ROOT/'configs/xtb_matrix.json'),
              'script_sha256':sha(__file__),'executable_sha256':sha(executable),'rdkit_version':rdkit.__version__,'python':sys.version,
              'boundary':BOUNDARY,'panel':[],'jobs':[],'geometry_failures':[]}
    for item in cfg['species']:
        try:xyz,meta=prepare_geometry(item,out,cfg)
        except Exception as exc:
            manifest['geometry_failures'].append({'id':item['id'],'error':repr(exc)});js(out/'manifest.json',manifest);continue
        manifest['panel'].append(meta)
        for method,environment,state in itertools.product(cfg['methods'],cfg['environments'],cfg['states']):
            jid=f"{item['id']}_gfn{method}_{environment}_{state['name']}"
            job={'id':jid,'molecule_id':item['id'],'gfn':method,'environment':environment,'state':state['name'],'charge':state['charge'],'uhf':state['uhf'],
                 'geometry_sha256':sha(xyz),'electron_count':meta['nuclear_charge_sum']-state['charge']}
            if (job['electron_count']-job['uhf'])%2:raise ValueError('Electron/spin parity mismatch')
            command=[executable,'input.xyz','--gfn',str(method),'--chrg',str(state['charge']),'--uhf',str(state['uhf']),'--sp','--acc',str(cfg['accuracy']),'--iterations',str(cfg['max_scc_iterations'])]
            if environment!='gas':command+=['--alpb',environment]
            dest=out/'raw'/jid;dest.mkdir(parents=True)
            start=time.perf_counter()
            with tempfile.TemporaryDirectory(prefix='p29_matrix_') as temp:
                shutil.copyfile(xyz,Path(temp)/'input.xyz')
                try:
                    result=subprocess.run(command,cwd=temp,env=env,capture_output=True,text=True,errors='replace',timeout=cfg['timeout_seconds'])
                    stdout,stderr=result.stdout,result.stderr;job.update(parse_output(stdout,stderr,result.returncode))
                except subprocess.TimeoutExpired as exc:
                    stdout=exc.stdout or '';stderr=exc.stderr or ''
                    if isinstance(stdout,bytes):stdout=stdout.decode(errors='replace')
                    if isinstance(stderr,bytes):stderr=stderr.decode(errors='replace')
                    job.update(status='timeout',energy_hartree=None,warning_count=0,numerical_warnings=[],scc_converged=False,normal_termination=False)
                (dest/'stdout.txt').write_text(stdout,encoding='utf8');(dest/'stderr.txt').write_text(stderr,encoding='utf8')
                shutil.copyfile(Path(temp)/'input.xyz',dest/'input.xyz')
            job.update(elapsed_seconds=time.perf_counter()-start,command=['xtb']+command[1:],stdout_sha256=sha(dest/'stdout.txt'),stderr_sha256=sha(dest/'stderr.txt'))
            js(dest/'job.json',job);manifest['jobs'].append(job)
            if len(manifest['jobs'])%12==0:print(f"{len(manifest['jobs'])}/{n} jobs; last {jid}: {job['status']}",flush=True)
        js(out/'manifest.json',manifest)
    manifest['completed_utc']=datetime.now(timezone.utc).isoformat();js(out/'manifest.json',manifest)
    analysis=summarize(out,manifest)
    files=[{'path':str(f.relative_to(out)).replace('\\','/'),'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(out.rglob('*')) if f.is_file()]
    js(out/'file_manifest.json',{'files':files,'excludes':'this manifest'})
    print(json.dumps({k:v for k,v in analysis.items() if k in ['success_jobs','total_jobs','jobs_with_warnings','method_difference_ranges_eV','solvent_shift_range_eV']},indent=2))
    return int(analysis['success_jobs']!=144)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--xtb',default='xtb');p.add_argument('--out',required=True)
    return run(p.parse_args())
if __name__=='__main__':raise SystemExit(main())
