"""Auditable structure calculations and SYNTHETIC methodology diagnostics.

No experimental selectivity, yield or electronic-structure model is inferred.
All kinetic rate constants are dimensioned hypothetical inputs, not fitted data.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import spearmanr
from rdkit import Chem, rdBase, DataStructs
from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolDescriptors, Draw, rdFingerprintGenerator

ROOT = Path(__file__).resolve().parents[1]
SEED = 29092026
FARADAY_C_MOL = 96485.33212

def save_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')

def csv_write(path, rows):
    if not rows:
        raise ValueError('CSV must declare rows; experimental data is a separate empty template')
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

def canonical(mol):
    copy = Chem.Mol(mol)
    for a in copy.GetAtoms(): a.SetAtomMapNum(0)
    return Chem.MolToSmiles(copy)

def panel_mol(item):
    mol = Chem.MolFromSmiles(item['mapped_smiles'])
    if mol is None: raise ValueError(item['id'])
    nxt = 7
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() == 0:
            atom.SetAtomMapNum(nxt); nxt += 1
    return mol

def carbonyl_partner():
    return Chem.MolFromSmiles('[C:100](=[O:101])([CH3:102])[c:103]1[cH:104][cH:105][c:106]([O:109][CH3:110])[cH:107][cH:108]1')

def proposed_product(substrate, site, partner=None):
    """Heavy-atom bookkeeping for C-H + C=O -> C-C(OH); not route validation."""
    idx = {a.GetAtomMapNum():a.GetIdx() for a in substrate.GetAtoms()}
    target = substrate.GetAtomWithIdx(idx[site])
    if target.GetTotalNumHs() != 1: raise ValueError(f'Site C{site} has no available C-H')
    partner = carbonyl_partner() if partner is None else partner
    combo = Chem.RWMol(Chem.CombineMols(substrate, partner))
    maps = {a.GetAtomMapNum():a.GetIdx() for a in combo.GetAtoms()}
    combo.GetAtomWithIdx(maps[site]).SetNumExplicitHs(0)
    combo.GetBondBetweenAtoms(maps[100], maps[101]).SetBondType(Chem.BondType.SINGLE)
    combo.GetAtomWithIdx(maps[101]).SetNumExplicitHs(1)
    combo.AddBond(maps[site], maps[100], Chem.BondType.SINGLE)
    product = combo.GetMol(); Chem.SanitizeMol(product)
    return product

def descriptors(mol):
    return dict(canonical_smiles=canonical(mol), formula=rdMolDescriptors.CalcMolFormula(mol),
                MW_g_mol=Descriptors.MolWt(mol), TPSA_A2=rdMolDescriptors.CalcTPSA(mol),
                calculated_logP=Crippen.MolLogP(mol), formal_charge=Chem.GetFormalCharge(mol),
                HBD=rdMolDescriptors.CalcNumHBD(mol), HBA=rdMolDescriptors.CalcNumHBA(mol),
                rotatable_bonds=rdMolDescriptors.CalcNumRotatableBonds(mol))

def calculate_structures(out, panel):
    desc, audit, confs, images = [], [], [], []
    sdf = Chem.SDWriter(str(out/'data/calculated/proposed_structures.sdf'))
    for item in panel['substrates']:
        mol = panel_mol(item); images.append(mol)
        entities = [('substrate', mol)] + [(f'proposed_C{s}', proposed_product(mol,s)) for s in (2,4,6)]
        for role, entity in entities:
            cid = f"{item['id']}_{role}"
            desc.append(dict(id=cid, family=item['family'], split=item['split'],
                             evidence_layer='calculated_descriptor_of_proposed_structure', **descriptors(entity)))
            mm = Chem.AddHs(entity)
            params = AllChem.ETKDGv3(); params.randomSeed=SEED; params.numThreads=2
            params.pruneRmsThresh=0.2
            ids = AllChem.EmbedMultipleConfs(mm, numConfs=4, params=params)
            if not ids: raise RuntimeError('Embedding failed: '+cid)
            opt = AllChem.MMFFOptimizeMoleculeConfs(mm, numThreads=2, maxIters=300, mmffVariant='MMFF94s')
            for cid_conf,(status, energy) in zip(ids,opt):
                confs.append(dict(id=cid, conformer=int(cid_conf), converged=status==0,
                                  MMFF94s_energy_kcal_mol=float(energy),
                                  meaning='geometry diagnostic; compare only conformers of same species'))
            eligible=[pair for pair in zip(ids,opt) if pair[1][0]==0]
            if not eligible: raise RuntimeError('No converged conformer: '+cid)
            best = min(eligible, key=lambda pair:pair[1][1])[0]
            mm.SetProp('_Name', cid); mm.SetProp('evidence_layer','proposed_identity_calculated_geometry')
            sdf.write(mm, confId=int(best))
        for site in (2,4,6):
            prod = proposed_product(mol,site)
            maps_before=sorted(a.GetAtomMapNum() for m in (mol,carbonyl_partner()) for a in m.GetAtoms())
            maps_after=sorted(a.GetAtomMapNum() for a in prod.GetAtoms())
            audit.append(dict(substrate_id=item['id'], site=f'C{site}',
                heavy_atom_map_conserved=maps_before==maps_after,
                mapping=f'{Chem.MolToSmiles(mol)}.{Chem.MolToSmiles(carbonyl_partner())}>>{Chem.MolToSmiles(prod)}',
                product_canonical_smiles=canonical(prod), product_identity_status='proposed_not_synthesized',
                stereochemistry='new carbinol stereocenter unassigned; no stereoselectivity claim',
                literature_scope_status=item['literature_scope_status'],
                exact_route_feasibility='unverified', charge=0, protonation='neutral free-base bookkeeping; in-cell speciation unknown'))
    sdf.close()
    csv_write(out/'data/calculated/descriptors.csv',desc)
    csv_write(out/'data/calculated/conformers.csv',confs)
    csv_write(out/'data/reaction_identity_audit.csv',audit)
    svg=Draw.MolsToGridImage(images,molsPerRow=3,subImgSize=(340,210),legends=[f"{p['id']} {p['name']}" for p in panel['substrates']],useSVG=True)
    (out/'figures/01_preregistered_structure_panel.svg').write_text(svg,encoding='utf-8')
    return desc,audit,confs

def kinetic_rhs(t,y,p):
    """Hours; S,P2,P4,P6,D are pyridine-equivalents; relay/theta dimensionless."""
    s,p2,p4,p6,d,relay,theta=y
    pulse = 1.0
    if p.get('pulse',False):
        pulse = 2.0 if (t % 1.0)<0.5 else 0.0
    flux=pulse*p['activity']
    # A surface channel competes with a relay-fed homogeneous channel.
    # Rates are deliberately hypothetical; no electrode potential calibration.
    r2=(p['k2_solution']*relay+p['k2_surface']*theta)*flux*s
    r4=(p['k4_solution']*relay+p['k4_surface']*theta)*flux*s
    r6=p['k6']*relay*flux*s
    l2=p['loss2']*relay*p2; l4=p['loss4']*relay*p4; l6=p['loss6']*relay*p6
    ls=p['side']*flux*s
    return [-r2-r4-r6-ls, r2-l2, r4-l4, r6-l6, ls+l2+l4+l6,
            (p['relay_setpoint']*pulse-relay)/p['tau_relay_h'],
            (p['surface']*relay/(1+relay)-theta)/p['tau_surface_h']]

def simulate(p,times=None):
    if times is None: times=np.linspace(0,9,181)
    initial=[1.,0.,0.,0.,0.,p['relay_setpoint'],p['surface']*p['relay_setpoint']/(1+p['relay_setpoint'])]
    # Segment at every pulse discontinuity so adaptive integration cannot skip pulses.
    knots=np.arange(0,9.0001,0.5) if p.get('pulse',False) else np.array([0.,9.])
    result=[]; state=initial
    for lo,hi in zip(knots[:-1],knots[1:]):
        sol=solve_ivp(lambda t,y:kinetic_rhs(t,y,p),(lo,hi),state,rtol=1e-8,atol=1e-10,dense_output=True,max_step=0.08)
        if not sol.success: raise RuntimeError(sol.message)
        chosen=times[(times>=lo)&(times<hi)]
        if len(chosen): result.extend(zip(chosen,sol.sol(chosen).T))
        state=sol.y[:,-1]
    result.append((9.,np.asarray(state)))
    return np.array([x[0] for x in result]),np.array([x[1] for x in result])

def selectivity(state):
    total=float(state[1]+state[2]+state[3]); denom=float(state[1]+state[2])
    return {'s2_C2_over_C2_plus_C4':None if denom<1e-8 else float(state[1]/denom),
            'total_three_isomer_yield':total,'C6_fraction_all_target':None if total<1e-8 else float(state[3]/total),
            'mass_balance':float(sum(state[:5]))}

def kinetic_diagnostics(out,config):
    base=config['kinetics']; variants={
        'surface_on':dict(base),
        'matched_support':dict(base,surface=0),
        'relay_removed':dict(base,relay_setpoint=0),
        'selective_C4_consumption':dict(base,surface=0,loss4=1.8),
        'equal_charge_pulse':dict(base,pulse=True)}
    rows=[]; terminal={}; trajectories={}
    for name,p in variants.items():
        times,y=simulate(p); trajectories[name]=(times,y)
        terminal[name]=selectivity(y[-1])
        for t,z in zip(times,y): rows.append(dict(scenario=name,time_h=float(t),S=z[0],P2=z[1],P4=z[2],P6=z[3],D=z[4],relay=z[5],theta=z[6],evidence_layer='synthetic_hypothetical_ODE'))
    csv_write(out/'data/synthetic/kinetic_trajectories.csv',rows)
    # Finite-difference identifiability: endpoint P2/P4 vs richer timed species observations.
    names=['k2_solution','k4_solution','k2_surface','k4_surface','loss2','loss4']
    def observations(p,rich,intervention=False):
        _,y=simulate(p)
        obs=y[[20,60,120,180],1:5].ravel() if rich else y[-1,1:3]
        if intervention:
            _,z=simulate(dict(p,surface=0))
            obs=np.concatenate([obs,z[[20,60,120,180],1:5].ravel()])
        return obs
    identification={}
    for label,rich,intervention in [('endpoint_only',False,False),('time_resolved',True,False),('time_resolved_plus_support_intervention',True,True)]:
        origin=observations(base,rich,intervention); jac=[]
        for name in names:
            upper=dict(base); lower=dict(base)
            upper[name]*=1.0001; lower[name]*=.9999
            jac.append((observations(upper,rich,intervention)-observations(lower,rich,intervention))/.0002)
        singular=np.linalg.svd(np.array(jac).T,compute_uv=False)
        identification[label]={'n_parameters':len(names),'n_observations':len(origin),'singular_values':singular.tolist(),'numerical_rank_relative_tolerance_1e-4':int(sum(singular>singular[0]*1e-4))}
    rng=np.random.default_rng(SEED); ensemble=[]
    for _ in range(config['uncertainty_draws']):
        p=dict(base)
        for key in names+['k6','side']: p[key]*=float(np.exp(rng.normal(0,0.65)))
        _,y=simulate(p); ensemble.append(selectivity(y[-1])['s2_C2_over_C2_plus_C4'])
    interval=np.quantile(ensemble,[.025,.5,.975]).tolist()
    payload={'evidence_layer':'synthetic_hypothetical_parameter_stress_test', 'terminal':terminal,
             'identifiability':identification,'s2_parameter_assumption_95pct_interval':interval,
             'uncertainty_distribution':'independent lognormal rate multipliers, log SD 0.65; not a posterior',
             'interpretation':'No chemical switch is predicted; rate assumptions alone span selectivity outcomes.'}
    save_json(out/'results/synthetic_kinetics.json',payload)
    return trajectories,payload,ensemble

def group_split(groups,seed=SEED):
    unique=np.unique(groups); rng=np.random.default_rng(seed); rng.shuffle(unique)
    n=len(unique); a=int(n*.55); b=int(n*.75)
    return [np.flatnonzero(np.isin(groups,x)) for x in (unique[:a],unique[a:b],unique[b:])]

def fit_ridge(x,y,lam=1.):
    center=x.mean(0); scale=x.std(0); scale[scale==0]=1
    design=np.column_stack([np.ones(len(x)),(x-center)/scale])
    penalty=np.diag([0.]+[lam]*x.shape[1])
    coef=np.linalg.solve(design.T@design+penalty,design.T@y)
    return lambda z:np.clip(np.column_stack([np.ones(len(z)),(z-center)/scale])@coef,0,1)

def conformal_radius(errors,alpha=.1):
    if not len(errors): raise ValueError('Calibration data required')
    order=math.ceil((len(errors)+1)*(1-alpha))
    if order>len(errors):return math.inf
    return float(np.sort(errors)[order-1])

def acquisition_score(site_probability,total_yield,impurity_risk,cost,information_value):
    """Proposed multiobjective expected useful product / cost + information.

    Inputs must come from calibrated independent evidence before deployment.
    Risk > 0.2 rejects a candidate; cost is joint reaction + analytical + separation.
    """
    if cost<=0:raise ValueError('Cost must be positive')
    if any(not 0<=v<=1 for v in (site_probability,total_yield,impurity_risk)):
        raise ValueError('Probabilities/fractions must be in [0, 1]')
    if impurity_risk>.2:return None
    return (site_probability*total_yield*(1-impurity_risk)+.15*information_value)/cost

def synthetic_learning(out):
    rng=np.random.default_rng(SEED)
    # Families, not regioisomers, are independent units. One family-level row only.
    groups=np.arange(120); x=rng.uniform(-1,1,(120,4))
    y=np.clip(.50+.22*x[:,0]-.15*x[:,1]+.10*x[:,0]*x[:,2]+rng.normal(0,.05,120),0,1)
    costs=rng.integers(1,5,120).astype(float)
    train,cal,test=group_split(groups); pred=fit_ridge(x[train],y[train])
    # Family-exchangeable split-conformal diagnostic for the constructed generator only.
    q=conformal_radius(np.abs(pred(x[cal])-y[cal])); yp=pred(x[test]); lo=np.maximum(0,yp-q); hi=np.minimum(1,yp+q)
    ood=np.any(np.abs(x[test])>np.max(np.abs(x[train]),axis=0)+.1,axis=1)
    rows=[dict(family=int(groups[i]),split=s,synthetic_site_fraction=float(y[i]),prediction=float(pred(x[[i]])[0]),cost_units=float(costs[i]),evidence_layer='synthetic_algorithm_check') for s,idx in [('train',train),('calibration',cal),('test',test)] for i in idx]
    csv_write(out/'data/synthetic/learning_rows.csv',rows)
    calibration=[]
    for j in range(5):
        mask=(yp>=j*.2)&(yp<(j+1)*.2 if j<4 else yp<=1)
        if mask.any(): calibration.append(dict(bin=j,count=int(mask.sum()),mean_prediction=float(yp[mask].mean()),mean_target=float(y[test][mask].mean())))
    # Budget matched cost-aware selection. Calibration/test never used for acquisitions.
    curves=[]
    seed_set=train[:6]; pool=train[6:]
    for policy in ['random','diversity_per_cost','uncertainty_proxy_per_cost']:
        selected=list(seed_set); remaining=list(pool); cumulative=float(costs[seed_set].sum())
        random_order=np.random.default_rng(SEED).permutation(len(pool)).tolist()
        for budget in [20,35,50,65,80,100,130]:
            while remaining:
                affordable=[i for i in remaining if cumulative+costs[i]<=budget]
                if not affordable: break
                if policy=='random': chosen=min(affordable,key=lambda i:random_order[list(pool).index(i)])
                else:
                    dist={i:float(np.min(np.linalg.norm(x[selected]-x[i],axis=1))) for i in affordable}
                    if policy=='diversity_per_cost': chosen=max(affordable,key=lambda i:dist[i]/costs[i])
                    else:
                        model=fit_ridge(x[selected],y[selected]); chosen=max(affordable,key=lambda i:(dist[i]+.1*model(x[[i]])[0])/costs[i])
                selected.append(chosen); remaining.remove(chosen); cumulative+=costs[chosen]
            model=fit_ridge(x[selected],y[selected])
            curves.append(dict(policy=policy,budget=float(budget),spent=float(cumulative),n_labels=len(selected),test_MAE=float(np.mean(np.abs(model(x[test])-y[test]))),evidence_layer='synthetic_single_seed_not_experiment_savings'))
    csv_write(out/'data/synthetic/budget_learning_curve.csv',curves)
    summary={'evidence_layer':'synthetic_algorithm_check','n_families':120,'split_sizes':list(map(len,[train,cal,test])),
             'ridge_test_MAE':float(np.mean(np.abs(yp-y[test]))),'constant_baseline_test_MAE':float(np.mean(np.abs(y[train].mean()-y[test]))),
             'conformal_nominal_coverage':.9,'conformal_observed_coverage':float(np.mean((lo<=y[test])&(y[test]<=hi))),
             'conformal_radius':q,'OOD_rejected':int(ood.sum()),'calibration_bins':calibration,
             'deployment_status':'ABSTAIN_ALL_REAL_REACTIONS_NO_EXPERIMENTAL_LABELS',
             'clinical_or_experimental_savings_claim':False}
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, WhiteKernel
    from sklearn.tree import DecisionTreeRegressor
    tree=DecisionTreeRegressor(max_depth=3,min_samples_leaf=6,random_state=SEED).fit(x[train],y[train])
    gp=make_pipeline(StandardScaler(),GaussianProcessRegressor(kernel=RBF(2.0)+WhiteKernel(.01),optimizer=None,normalize_y=True)).fit(x[train],y[train])
    summary['tree_test_MAE']=float(np.mean(np.abs(tree.predict(x[test])-y[test])))
    summary['GP_test_MAE']=float(np.mean(np.abs(gp.predict(x[test])-y[test])))
    summary['hyperparameters']='fixed a priori; no test optimization'
    summary['proposed_deployment_acquisition_function']='probability*total_yield*(1-impurity_risk)/joint_cost + 0.15*information/joint_cost; reject risk>0.2; NOT evaluated in this simulation'
    summary['actual_evaluated_acquisition_rule']='uncertainty_proxy_per_cost=(nearest_feature_distance+0.1*ridge_predicted_synthetic_site_fraction)/cost; diversity_per_cost=nearest_feature_distance/cost'
    save_json(out/'results/synthetic_learning.json',summary)
    return summary,curves

def manufacturing_metrics(current_A,time_h,voltage_V,substrate_mol,yield_fraction,mw_g_mol,inputs_g,water_g,volume_L):
    if not(0<yield_fraction<=1) or min(current_A,time_h,voltage_V,substrate_mol,mw_g_mol,volume_L)<=0:
        raise ValueError('Positive dimensions and usable yield required')
    mass=substrate_mol*yield_fraction*mw_g_mol
    energy=current_A*voltage_V*time_h/1000
    return {'product_g':mass,'charge_C':current_A*time_h*3600,
            'energy_kWh_per_kg':energy/(mass/1000), 'STY_g_L_h':mass/volume_L/time_h,
            'PMI_without_water':inputs_g/mass,'PMI_including_water':(inputs_g+water_g)/mass,
            'FE_if_2e_percent':100*2*FARADAY_C_MOL*substrate_mol*yield_fraction/(current_A*time_h*3600)}

def manufacturing(out):
    # Only Q and conditional FE are literature-derived. Voltage/solvent masses below are scenarios.
    rows=[]
    for recovery in [.5,.7,.9]:
        for s2 in [.2,.5,.8]:
            m=manufacturing_metrics(.025,9,5,.0002,.81*s2*recovery,305.376,25,10,.010)
            rows.append(dict(scenario='ASSUMED_C2_target_cost_sensitivity',s2=s2,isolated_recovery=recovery,
                             assumed_total_assay_yield=.81,assumed_cell_voltage_V=5,assumed_inputs_g=25,
                             assumed_water_g=10,assumed_volume_L=.010,**m))
    csv_write(out/'data/synthetic/manufacturing_scenarios.csv',rows)
    q=.025*9*3600
    summary={'literature_NC_standard':{'substrate_mol':.0002,'isolated_major_yield':.81,'current_A':.025,'time_h':9,
            'charge_C_calculated':q,'charge_F_per_mol_substrate':q/FARADAY_C_MOL/.0002,
            'FE_percent_recalculated_if_two_electrons':100*2*FARADAY_C_MOL*.0002*.81/q,
            'reported_FE_percent':3.86,'electron_assumption':'two electrons per product; source SI S6 FE convention verified',
            'measured_cell_voltage_V':None,'published_own_experiment':False},
            'energy_or_route_advantage':'NOT_DETERMINED_voltage_full_inventory_and_matched_route_missing',
            'scenario_table_evidence':'hypothetical unit and recovery sensitivity only'}
    save_json(out/'results/manufacturing.json',summary)
    return summary,rows

def figures(out,trajectories,kin,ensemble,learn,curves,manufact):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':130})
    fig,ax=plt.subplots(figsize=(9,5))
    for name,(t,y) in trajectories.items():
        denominator=y[:,1]+y[:,2]; valid=denominator>1e-7
        ax.plot(t[valid],y[valid,1]/denominator[valid],label=name)
    ax.set(xlabel='Time (h)',ylabel='C2 / (C2 + C4)',ylim=(0,1),title='SYNTHETIC kinetic stress test: no measured site switch')
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(out/'figures/02_synthetic_site_competition.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4));ax.hist(ensemble,bins=20,color='#2a768e'); ax.axvline(.8,color='black',linestyle='--');ax.axvline(.2,color='black',linestyle='--')
    ax.set(xlabel='Final C2 / (C2 + C4)',ylabel='Hypothetical parameter draws',title='SYNTHETIC uncertainty: assumed rates, no fitted posterior');fig.tight_layout();fig.savefig(out/'figures/03_assumed_parameter_uncertainty.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,5)); bins=learn['calibration_bins'];ax.plot([0,1],[0,1],'--',color='grey')
    ax.scatter([b['mean_prediction'] for b in bins],[b['mean_target'] for b in bins],s=[20*b['count'] for b in bins]);ax.set(xlim=(0,1),ylim=(0,1),xlabel='Mean synthetic prediction',ylabel='Mean synthetic label',title='SYNTHETIC calibration; real deployment = ABSTAIN');fig.tight_layout();fig.savefig(out/'figures/04_synthetic_calibration.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4))
    for policy in sorted(set(x['policy'] for x in curves)):
        sel=[x for x in curves if x['policy']==policy];ax.plot([x['budget'] for x in sel],[x['test_MAE'] for x in sel],marker='o',label=policy)
    ax.set(xlabel='Same available labeling cost (arbitrary units)',ylabel='Frozen synthetic test MAE',title='SYNTHETIC single-seed workflow test; no real savings');ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'figures/05_synthetic_budget_learning.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4))
    for rec in [.5,.7,.9]:
        sel=[x for x in manufact if x['isolated_recovery']==rec];ax.plot([x['s2'] for x in sel],[x['energy_kWh_per_kg'] for x in sel],marker='o',label=f'Assumed recovery {rec}')
    ax.set(xlabel='Assumed desired-site fraction',ylabel='Electricity (kWh/kg isolated desired isomer)',title='SCENARIO ONLY: 5 V, 25 mA, 9 h; separation losses matter');ax.legend();fig.tight_layout();fig.savefig(out/'figures/06_assumed_separation_burden.png');plt.close(fig)

def run(out):
    for directory in ['data/calculated','data/synthetic','data/experimental','results','figures']:(out/directory).mkdir(parents=True,exist_ok=True)
    panel=json.loads((ROOT/'configs/preregistered_panel.json').read_text(encoding='utf-8'))
    config=json.loads((ROOT/'configs/compute.json').read_text(encoding='utf-8'))
    frozen={'seed':SEED,'threads_max':2,'ETKDGv3_candidates_per_species':4,'MMFF94s_max_iters':300,'kinetic_time_h':9}
    for key,expected in frozen.items():
        if config[key]!=expected:raise ValueError(f'Unsupported frozen configuration change: {key}; update implementation and preregistration first')
    desc,audit,confs=calculate_structures(out,panel)
    tr,kin,ensemble=kinetic_diagnostics(out,config);learn,curves=synthetic_learning(out);man,scenarios=manufacturing(out)
    from literature_audit import run_public_audit
    public_audit,electrical=run_public_audit(out)
    figures(out,tr,kin,ensemble,learn,curves,scenarios)
    manifest={'seed':SEED,'python':sys.version,'platform':platform.platform(),'rdkit':rdBase.rdkitVersion,
              'numpy':np.__version__,'descriptor_structures':len(desc),'mapped_proposed_reactions':len(audit),
              'conformers':len(confs),'conformer_converged':sum(r['converged'] for r in confs),
              'experimental_rows':0,'main_pipeline_electronic_structure_calculations':0,
              'optional_quantum_preflight':'separate results/quantum_preflight/manifest.json, not rerun by this pipeline','synthetic_ODE_scenarios':len(tr),
              'config_sha256':hashlib.sha256((ROOT/'configs/compute.json').read_bytes()).hexdigest(),
              'panel_sha256':hashlib.sha256((ROOT/'configs/preregistered_panel.json').read_bytes()).hexdigest()}
    save_json(out/'results/execution_manifest.json',manifest)
    print(json.dumps(manifest,indent=2))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir','--out',type=Path,help='Required new output directory; refuses an existing directory')
    parser.add_argument('--refresh-committed',action='store_true',help='Authoring only: explicitly refresh the committed phase29 evidence tree')
    parser.add_argument('--self-test',action='store_true',help='Run meaningful structure, kinetics, leakage and unit tests without full analysis')
    args=parser.parse_args()
    if args.self_test:
        import unittest
        suite=unittest.defaultTestLoader.discover(str(ROOT/'tests'))
        result=unittest.TextTestRunner(verbosity=2).run(suite)
        raise SystemExit(0 if result.wasSuccessful() else 1)
    if args.refresh_committed:
        if args.output_dir is not None:parser.error('--refresh-committed and --out are mutually exclusive')
        run(ROOT)
    else:
        if args.output_dir is None:parser.error('Provide --out NEW_DIRECTORY to run without overwriting evidence')
        target=args.output_dir.resolve()
        if target.exists():parser.error('Output directory already exists; choose a new directory')
        run(target)

if __name__=='__main__':main()
