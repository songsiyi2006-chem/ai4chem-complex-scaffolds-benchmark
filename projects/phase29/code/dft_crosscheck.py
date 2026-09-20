"""Bounded Psi4 gas-phase crosscheck at frozen molecular geometries.

Requires an existing Psi4 Python interpreter. No installation, geometry search,
surface model, implicit solvent, thermal correction or electrode conversion.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
HARTREE_EV = 27.211386245988


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf8")


def pair_difference(cation, radical):
    for key in ('molecule_id', 'geometry_sha256', 'method', 'basis', 'environment'):
        if cation[key] != radical[key]:
            raise ValueError('Incompatible charge-state pair: '+key)
    if (cation['charge'], cation['multiplicity'], radical['charge'], radical['multiplicity']) != (1,1,0,2):
        raise ValueError('Wrong charge/spin states')
    if cation['status'] != 'converged' or radical['status'] != 'converged':
        raise ValueError('Only converged states may be paired')
    return (radical['energy_hartree']-cation['energy_hartree'])*HARTREE_EV


def worker(spec_path):
    import numpy as np
    import psi4
    spec_path=Path(spec_path).resolve()
    job=json.loads(spec_path.read_text(encoding='utf8'))
    folder=spec_path.parent
    xyz=folder/'input.xyz'
    if sha(xyz)!=job['geometry_sha256']:
        raise ValueError('Geometry checksum mismatch')
    cfg=job['config']
    psi4.set_memory(cfg['memory'])
    psi4.set_num_threads(cfg['threads'])
    psi4.core.set_output_file(str(folder/'psi4.out'), False)
    result={key:value for key,value in job.items() if key!='config'}
    result.update(psi4_version=psi4.__version__,status='started',started_utc=datetime.now(timezone.utc).isoformat())
    try:
        lines=xyz.read_text(encoding='utf8').splitlines()
        coordinates='\n'.join(lines[2:2+int(lines[0])])
        molecule=psi4.geometry(f"{job['charge']} {job['multiplicity']}\n{coordinates}\nunits angstrom\nsymmetry c1\nno_reorient\nno_com")
        electrons=int(round(sum(molecule.Z(i) for i in range(molecule.natom()))))-job['charge']
        unpaired=job['multiplicity']-1
        if (electrons-unpaired)%2 or electrons<unpaired:
            raise ValueError('Electron/spin parity mismatch')
        reference='rks' if job['multiplicity']==1 else 'uks'
        options={key:cfg[key] for key in ('scf_type','e_convergence','d_convergence','maxiter','dft_radial_points','dft_spherical_points')}
        options.update(basis=job['basis'],reference=reference,guess='sad')
        psi4.set_options(options)
        with tempfile.TemporaryDirectory(prefix='phase29_psi4_') as scratch:
            psi4.core.IOManager.shared_object().set_default_path(scratch)
            energy,wfn=psi4.energy(job['method'],molecule=molecule,return_wfn=True)
            ca=wfn.Ca_subset('AO','OCC').np
            cb=wfn.Cb_subset('AO','OCC').np
            overlap=ca.T@wfn.S().np@cb
            nalpha,nbeta=wfn.nalpha(),wfn.nbeta()
            sz=(nalpha-nbeta)/2
            s2=sz*(sz+1)+nbeta-float(np.sum(overlap*overlap))
            assert nalpha+nbeta==electrons
            result.update(status='converged',energy_hartree=float(energy),reference=reference,
                          electrons=electrons,nalpha=nalpha,nbeta=nbeta,nbf=wfn.basisset().nbf(),
                          spin_squared=float(s2),target_spin_squared=unpaired/2*(unpaired/2+1),
                          spin_squared_formula='Sz*(Sz+1)+Nbeta-sum_ij|Ca_occ.T S Cb_occ|^2',
                          wavefunction_stability='not_checked',scf_options=options)
            psi4.core.clean()
    except Exception as exc:
        result.update(status='failed',error_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc())
    finally:
        psi4.core.close_outfile()
        result['finished_utc']=datetime.now(timezone.utc).isoformat()
        result['engine_output_sha256']=sha(folder/'psi4.out') if (folder/'psi4.out').exists() else None
        save(folder/'result.json',result)
    return 0 if result['status']=='converged' else 1


def analyze(out):
    jobs=[json.loads(path.read_text(encoding='utf8')) for path in sorted((out/'raw').glob('*/result.json'))]
    rows=[]
    for molecule,basis in sorted({(j['molecule_id'],j['basis']) for j in jobs}):
        states={j['state']:j for j in jobs if j['molecule_id']==molecule and j['basis']==basis}
        if set(states)!= {'cation','radical'} or any(j['status']!='converged' for j in states.values()):
            continue
        rows.append({'molecule_id':molecule,'method':states['cation']['method'],'basis':basis,
                     'environment':'gas','geometry_sha256':states['cation']['geometry_sha256'],
                     'delta_E_radical_minus_cation_eV':pair_difference(states['cation'],states['radical']),
                     'radical_spin_squared':states['radical']['spin_squared']})
    with (out/'energy_differences.csv').open('w',encoding='utf8',newline='') as handle:
        fields=['molecule_id','method','basis','environment','geometry_sha256','delta_E_radical_minus_cation_eV','radical_spin_squared']
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    summary={'evidence_type':'executed_gas_phase_PBE0_single_points_at_MMFF_geometry','jobs_attempted':len(jobs),
             'jobs_converged':sum(j['status']=='converged' for j in jobs),'pairs':rows,
             'limits':['DFT is a model comparator, not ground truth','one cation geometry per molecule; not DFT optimized',
                       'no diffuse-basis convergence study','basis check covers Q01 only','wavefunction stability not checked',
                       'no solvent/electrode/counterion/free-energy or reaction-barrier calculation'],
             'basis_change_eV':{}}
    for molecule in {r['molecule_id'] for r in rows}:
        values={r['basis']:r['delta_E_radical_minus_cation_eV'] for r in rows if r['molecule_id']==molecule}
        if set(values)=={'def2-svp','def2-tzvp'}:
            summary['basis_change_eV'][molecule]=values['def2-tzvp']-values['def2-svp']
    save(out/'summary.json',summary)
    return summary


def run(args):
    config_path=ROOT/'configs/dft_crosscheck.json'
    cfg=json.loads(config_path.read_text(encoding='utf8'))
    out=Path(args.out).resolve()
    out.mkdir(parents=True,exist_ok=False)
    script=Path(__file__).resolve()
    manifest={'started_utc':datetime.now(timezone.utc).isoformat(),'script_sha256':sha(script),
              'config_sha256':sha(config_path),'config':cfg,'sequential_jobs':True,'jobs':[]}
    save(out/'manifest.json',manifest)
    env=dict(os.environ,OMP_NUM_THREADS=str(cfg['threads']),MKL_NUM_THREADS=str(cfg['threads']),OPENBLAS_NUM_THREADS=str(cfg['threads']),PYTHONUTF8='1')
    interpreter=Path(args.python).resolve()
    # Conda Windows interpreters need their own DLL directory, not a different environment's DLLs.
    dll=interpreter.parent/'Library/bin'
    if dll.is_dir():
        env['PATH']=str(dll)+os.pathsep+str(interpreter.parent)+os.pathsep+env.get('PATH','')
    configurations=[(m,cfg['basis']) for m in cfg['molecules']]+[(m,cfg['basis_check']) for m in cfg['basis_check_molecules']]
    for molecule,basis in configurations:
        for state in cfg['states']:
            job_id=f"{molecule}_{basis}_{state['name']}"
            folder=out/'raw'/job_id;folder.mkdir(parents=True)
            source=ROOT/'results/quantum_preflight/inputs'/f'{molecule}.xyz'
            shutil.copy2(source,folder/'input.xyz')
            job={'id':job_id,'molecule_id':molecule,'method':cfg['method'],'basis':basis,'environment':'gas',
                 'state':state['name'],'charge':state['charge'],'multiplicity':state['multiplicity'],
                 'geometry_sha256':sha(source),'config':cfg}
            spec=folder/'input.json';save(spec,job)
            command=[str(interpreter),str(script),'--worker',str(spec)]
            started=time.perf_counter()
            try:
                result=subprocess.run(command,capture_output=True,text=True,encoding='utf8',errors='replace',env=env,
                                      cwd=folder,timeout=cfg['timeout_seconds_per_job'])
                stdout,stderr=result.stdout,result.stderr
                returncode=result.returncode
            except subprocess.TimeoutExpired as exc:
                stdout,stderr=exc.stdout or '',exc.stderr or ''
                if isinstance(stdout,bytes):stdout=stdout.decode('utf8',errors='replace')
                if isinstance(stderr,bytes):stderr=stderr.decode('utf8',errors='replace')
                returncode=None
            (folder/'stdout.txt').write_text(stdout,encoding='utf8')
            (folder/'stderr.txt').write_text(stderr,encoding='utf8')
            result_path=folder/'result.json'
            if not result_path.exists():
                failure={key:value for key,value in job.items() if key!='config'}
                failure['status']='timeout' if returncode is None else 'worker_crash'
                save(result_path,failure)
            result_record=json.loads(result_path.read_text(encoding='utf8'))
            manifest['jobs'].append({'id':job_id,'status':result_record['status'],'returncode':returncode,
                                     'elapsed_seconds':time.perf_counter()-started,'result_sha256':sha(result_path)})
            save(out/'manifest.json',manifest)
            print(job_id,result_record['status'],flush=True)
    summary=analyze(out)
    manifest['finished_utc']=datetime.now(timezone.utc).isoformat()
    save(out/'manifest.json',manifest)
    print(json.dumps(summary,indent=2))
    return 0 if summary['jobs_attempted']==summary['jobs_converged'] else 1


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python',default=sys.executable,help='Existing Psi4 Python interpreter')
    parser.add_argument('--out',type=Path,help='New output directory')
    parser.add_argument('--worker',type=Path,help=argparse.SUPPRESS)
    parser.add_argument('--analyze',type=Path,help='Reanalyze an existing completed output without electronic calculations')
    args=parser.parse_args()
    if args.worker:return worker(args.worker)
    if args.analyze:
        print(json.dumps(analyze(args.analyze),indent=2));return 0
    if args.out is None:parser.error('Provide --out NEW_DIRECTORY')
    return run(args)


if __name__=='__main__':
    raise SystemExit(main())
