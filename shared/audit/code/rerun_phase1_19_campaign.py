"""Isolated, logged Phase 1-19 reruns. Does not commit, push, or power off.

Each attempt snapshots current Python source into a fresh directory; historical
results/caches are not copied. Exit zero means process completion, NOT scientific
acceptance. Review phase outputs for explicit failure/fallback status.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT=Path(__file__).resolve().parent
PRIMARY=Path('C:/Users/HUIWEI/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe')
MOLECULAR=Path('C:/Users/HUIWEI/miniconda3/envs/phase2ff/python.exe')
DEPS=ROOT.parent/'phase24-deps'
EXTRA=ROOT.parent/'phasefix-test-deps'
SCRIPTS={1:'molecule_benchmark.py',2:'run_heavy_dynamics_benchmark.py'}
SCRIPTS.update({i:next(ROOT.glob(f'run_phase{i}_*.py')).name for i in range(3,20)})
ARGS={1:['--workers','1'],2:['--ff_env',str(MOLECULAR)],3:['--force_rerun'],
      4:['--force_rerun'],7:['--threads','2','--force'],8:['--threads','2']}

def save(path,data):
    temp=path.with_suffix('.tmp')
    temp.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    temp.replace(path)

def reconcile_exit(attempt,pid,code,evidence):
    """Record a separately observed exit after a queue parent was stopped.

    Call only with an exit code obtained from a retained child process handle.
    This does not validate scientific output or infer exit codes from PID absence.
    """
    path=attempt/'rerun_status.json'
    rec=json.loads(path.read_text(encoding='utf-8'))
    if rec.get('pid')!=pid or rec['status']!='running':
        raise ValueError('PID/status mismatch; refuse to overwrite an unrelated record')
    rec.update(status='process_completed' if code==0 else 'process_failed',exit_code=code,
               exit_reconciliation_evidence=evidence,exit_observed_at=time.strftime('%Y%m%dT%H%M%S'))
    rec['elapsed_s']=time.time()-time.mktime(time.strptime(rec['started'],'%Y%m%dT%H%M%S'))
    rec['elapsed_is_observation_upper_bound']=True
    save(path,rec)

def runtime_metadata(py,env):
    code="""import importlib.metadata as m, json, platform, sys
versions={}
for name in ['numpy','scipy','matplotlib','rdkit','torch','openmm','scikit-learn','ase','psi4','pyscf']:
    try: versions[name]=m.version(name)
    except m.PackageNotFoundError: versions[name]=None
print(json.dumps(dict(python=sys.version,executable=sys.executable,platform=platform.platform(),distributions=versions)))
"""
    try:
        probe=subprocess.run([str(py),'-c',code],env=env,capture_output=True,text=True,
                             encoding='utf-8',errors='replace',timeout=30)
        if probe.returncode==0: return json.loads(probe.stdout)
        return dict(error=probe.stderr[-1000:],exit_code=probe.returncode)
    except (OSError,subprocess.TimeoutExpired,json.JSONDecodeError) as exc:
        return dict(error=str(exc))

def run_one(phase,campaign,extra_args=(),profile=None,phase4_attempt=None):
    stamp=time.strftime('%Y%m%dT%H%M%S')
    dest=campaign/f'phase{phase:02d}'/stamp
    dest.mkdir(parents=True,exist_ok=False)
    sources={}
    for src in ROOT.glob('*.py'):
        if src.name.startswith('output_'): continue  # generated protocols are not input source
        shutil.copy2(src,dest/src.name)
        # Hash the executed snapshot, not a source that another worker may edit.
        sources[src.name]=hashlib.sha256((dest/src.name).read_bytes()).hexdigest()
    prerequisites=[]
    if phase in (6,7) and phase4_attempt:
        required=['reactant_3d.mol']
        if phase==6: required += ['product_3d.mol','images_idpp.xyz']
        for name in required:
            src=phase4_attempt.resolve()/'results_phase4'/name
            if not src.is_file(): raise FileNotFoundError(src)
            target=dest/'results_phase4'/name
            target.parent.mkdir(exist_ok=True)
            shutil.copy2(src,target)
            prerequisites.append(dict(source=str(src),destination=str(target.relative_to(dest)),
                                      sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
    use_primary=phase in [1,9,14,16,17,18]
    if profile: use_primary=profile=='primary'
    py=PRIMARY if use_primary else MOLECULAR
    env=os.environ.copy()
    env.update(PYTHONUNBUFFERED='1',PYTHONIOENCODING='utf-8',OMP_NUM_THREADS='2',
               MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',NUMEXPR_NUM_THREADS='2',OPENMM_CPU_THREADS='2')
    simulator=ROOT.parent/'phase24-ot2-env'/'Scripts'/'opentrons_simulate.exe'
    if simulator.is_file(): env['OPENTRONS_SIMULATE_EXE']=str(simulator)
    if use_primary:
        env['PYTHONPATH']=os.pathsep.join([str(DEPS),str(EXTRA)])
    else:
        env.pop('PYTHONPATH',None)
        env['PATH']=str(py.parent/'Library'/'bin')+os.pathsep+str(py.parent)+os.pathsep+env.get('PATH','')
    args=[str(py),'-X','faulthandler','-u',SCRIPTS[phase]]+ARGS.get(phase,[])+list(extra_args)
    scope='default/full attempt'
    if extra_args:
        scope='default/full attempt (platform override)' if list(extra_args) in (['--platform','CPU'],['--platform','Reference']) else 'custom-arguments attempt'
    rec=dict(phase=phase,started=stamp,source_sha256=sources,command=args,
             cwd=str(dest),status='running',scope=scope,
             prerequisite_artifacts=prerequisites,
             historical_cache_copied=False,scientific_acceptance='pending output review')
    rec['runtime']=runtime_metadata(py,env)
    rec['thread_limits']={k:env[k] for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS',
        'OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','OPENMM_CPU_THREADS')}
    save(dest/'rerun_status.json',rec)
    start=time.monotonic()
    print(f'PHASE {phase} START {dest}',flush=True)
    with (dest/'run.log').open('w',encoding='utf-8') as log:
        try:
            proc=subprocess.Popen(args,cwd=dest,env=env,stdout=log,stderr=subprocess.STDOUT)
            rec['pid']=proc.pid;save(dest/'rerun_status.json',rec)
            code=proc.wait()
        except OSError as exc:
            log.write(f'LAUNCH ERROR: {exc}\n')
            rec['launch_error']=str(exc)
            code=127
    rec.update(status='process_completed' if code==0 else 'process_failed',exit_code=code,
               elapsed_s=time.monotonic()-start)
    save(dest/'rerun_status.json',rec)
    print(f'PHASE {phase} EXIT {code} {rec["elapsed_s"]:.1f}s',flush=True)
    return rec

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phases')
    parser.add_argument('--campaign',type=Path,default=ROOT.parent/'phase1-19-rerun-20260910')
    parser.add_argument('--profile',choices=['primary','molecular'])
    parser.add_argument('--phase4-attempt',type=Path,help='Fresh Phase4 snapshot supplying geometry prerequisites for6/7')
    parser.add_argument('--reconcile-attempt',type=Path)
    parser.add_argument('--observed-pid',type=int)
    parser.add_argument('--observed-exit-code',type=int)
    parser.add_argument('--exit-evidence')
    args,extra=parser.parse_known_args()
    if args.reconcile_attempt:
        if args.observed_pid is None or args.observed_exit_code is None or not args.exit_evidence:
            parser.error('Reconciliation requires PID, directly observed exit code, and evidence description')
        reconcile_exit(args.reconcile_attempt,args.observed_pid,args.observed_exit_code,args.exit_evidence)
        raise SystemExit(0)
    if not args.phases: parser.error('--phases is required for computation')
    phases=[int(s) for s in args.phases.split(',')]
    if any(p not in SCRIPTS for p in phases):parser.error('Phases must be 1..19')
    for phase in phases: run_one(phase,args.campaign,extra,args.profile,args.phase4_attempt)
