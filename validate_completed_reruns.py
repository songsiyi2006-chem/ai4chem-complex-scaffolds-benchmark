"""Validate selected fresh outputs; never changes original numerical results."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np

def validate(attempt,phase):
    checks={}
    if phase==1:
        data=json.loads((attempt/'bench_results/benchmark_results.json').read_text(encoding='utf-8'))
        entries=data['results']
        checks['invalid_M09_rejected']=any(r['id']=='M09' and r['status']=='failed_parse' for r in entries)
        checks['ten_valid_entries']=len([r for r in entries if r['status']=='ok'])==10
        checks['all_valid_have_50_accepted']=all(r['conformers']['status']=='ok' and r['conformers']['n_accepted']==50 for r in entries if r['id']!='M09')
        limits='Finite converged force-field conformers only; no experimental or GNN accuracy validation.'
    elif phase==2:
        data=json.loads((attempt/'results_phase2/phase2_results.json').read_text(encoding='utf-8'))
        scans=data['stage2_torsion_scans']
        checks['all_stages_ok']=data['all_stages_ok'] is True
        checks['five_36_point_scans']=len(scans)==5 and all(r['status']=='ok' and r['n_points']==36 for r in scans.values())
        checks['angles_achieved']=all(r.get('achieved_max_dev_deg',float('inf'))<=5.1 for r in scans.values())
        checks['unrestrained_energies']=all(r.get('energy_includes_restraint') is False and np.isfinite(r['rel_energies']).all() for r in scans.values())
        md=data['stage3_md']
        checks['full_100000_step_MD']=md['production_steps']==100000 and md['n_frames']==200
        checks['finite_MD_trace']=all(all(np.isfinite(v) for v in row.values() if isinstance(v,(float,int))) for row in md['rows'])
        limits=f'200 ps implicit-solvent trajectory, not long-time convergence. Parameterization: {md["force_field"]}. Old reporter ETA counted equilibration incorrectly; display-only fix does not change trajectory.'
    elif phase==4:
        data=json.loads((attempt/'results_phase4/phase4_results.json').read_text(encoding='utf-8'))
        checks['all_stages_ok']=data.get('all_stages_ok') is True
        checks['NEB_converged']=data.get('stage2_neb',{}).get('converged') is True
        neb=data.get('stage2_neb',{})
        limits=(f"NEB residual force {neb.get('fmax_final_ev_A')} eV/A; target {neb.get('fmax_target')}. "
                'An unconverged path cannot establish a TS, spectrum or kinetic barrier. This check is not an independent Hessian/IRC validation.')
    elif phase==5:
        data=json.loads((attempt/'results_phase5/phase5_results.json').read_text(encoding='utf-8'))
        process=json.loads((attempt/'rerun_status.json').read_text(encoding='utf-8'))
        checks['execution_completed']=process.get('exit_code')==0 and data.get('execution_status')!='failed'
        checks['default_modules_present']=all(k in data for k in ['module_A','module_B','module_C'])
        checks['no_fatal_error']=not data.get('fatal_error')
        limits=(str(data.get('fatal_error','No fatal error recorded'))+
                '; constrained scan peaks are not validated stationary TSs. Default A/B/C only; D exploratory kinetics is not enabled.')
    elif phase==9:
        vals=json.loads((attempt/'results_phase9/ot2_validation_report.json').read_text(encoding='utf-8'))
        for key in ['round1','champion']:
            checks[key+'_official_simulation']=bool(vals[key].get('hardware_validation_passed',False))
        data=json.loads((attempt/'results_phase9/phase9_results.json').read_text(encoding='utf-8'))
        checks['50_synthetic_experiments']=data['module_9B']['n_experiments']==50
        checks['synthetic_origin_declared']=data.get('hardware_executed') is False
        limits='Protocol/software simulation only; approximate surrogate and batch EI, no wet-lab validation. HPLC errors must be reviewed separately.'
    elif phase==10:
        data=json.loads((attempt/'results_phase10/phase10_results.json').read_text(encoding='utf-8'))
        checks['full_110_training_episodes']=data['mode']=='full' and data['sac']['episodes']==110
        for key in ['tlog_open','tlog_ctrl']:
            trace=data[key]
            # JSON action log starts at 5 s; NPY state history also contains t=0.
            checks[key+'_complete_saved_trace']=len(trace)==720 and trace[0]['t']==5. and 3600.<=trace[-1]['t']<3600.02
            checks[key+'_finite']=all(all(np.isfinite(v) for v in row.values() if isinstance(v,(float,int))) for row in trace)
        checks['PAT_audit_finite']=np.isfinite(list(data['pat_audit'].values())).all().item()
        limits=('Synthetic plant/controller only, no hardware safety certification. '
                'Old integration floor caused endpoint overshoots of 10.4/2.1 ms; future-source fix bounds the final step. '
                'Saved-state rate checks found no recorded reaction timestep requirement below 1 ms; unsaved microsteps cannot be certified. '
                'Zero observed breaches in eight shielded validation cases is not a universal guarantee.')
    elif phase==18:
        path=attempt/'run_phase18_wholecell_metabolic_thermodynamics.py'
        spec=importlib.util.spec_from_file_location('rerun18_snapshot',path)
        mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
        net=mod.build_all()
        z=np.load(attempt/'results_phase18/phase18_tfba.npz')
        data=json.loads((attempt/'results_phase18/phase18_results.json').read_text(encoding='utf-8'))
        internal=[i for i,m in enumerate(net.mnames) if net.mets[m]['cls']!='ext']
        # External reservoir rows are open boundaries, not steady-state constraints.
        residual=float(np.max(np.abs(net.S[internal]@z['v'])))
        checks['stoichiometry_residual_below_1e-5']=residual<1e-5
        checks['all_fluxes_finite']=bool(np.isfinite(z['v']).all())
        checks['no_thermo_sign_violations']=data['tfba']['n_sign_violations']==0
        checks['no_closed_flux_cycle']=data['tfba']['zero_loop_certificate']['pass'] is True
        checks['enzyme_budget_satisfied']=data['tfba']['enzyme_budget_fraction']<=1+1e-6
        for key in ['moiety_drift_adenylate','moiety_drift_nad']:
            checks[key+'_below_1e-8']=data['dynamics'][key]<1e-8
        checks['all_dynamic_states_finite']=bool(np.isfinite(np.load(attempt/'results_phase18/phase18_dynamics.npz')['y']).all())
        limits=f'Assigned/curated network, not a validated organism. Stoichiometry residual={residual:.3g}; formation-energy benchmark RMS={data["curation"]["rms_kJ_mol"]:.3g} kJ/mol, not chemically accurate throughout.'
    else: raise ValueError('Implemented for1,2,4,5,9,10,18; other phases need phase-specific review')
    record=dict(phase=phase,checks=checks,passed=all(checks.values()),limitations=limits,
                validator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (attempt/'acceptance_review.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    status=attempt/'rerun_status.json'
    rec=json.loads(status.read_text(encoding='utf-8'))
    rec['scientific_acceptance']='specified numerical/software checks passed; see limitations' if record['passed'] else 'acceptance checks FAILED'
    status.write_text(json.dumps(rec,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(record,indent=2))
    return record['passed']

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('phase',type=int)
    ap.add_argument('attempt',type=Path)
    args=ap.parse_args()
    raise SystemExit(0 if validate(args.attempt,args.phase) else 1)
