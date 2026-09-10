"""Run cheap regressions in suitable isolated interpreters and retain evidence."""
import hashlib
import argparse
import json
import os
import re
from pathlib import Path
import subprocess
import time
from rerun_phase1_19_campaign import PRIMARY,MOLECULAR,DEPS,EXTRA

ROOT=Path(__file__).resolve().parent
GROUPS={
    'primary': ['test_phase1_19_campaign','test_phase1_5_audit','test_phase6_10_audit',
                'test_phase11_13_audit','test_phase14_17_audit','test_phase18_19_audit',
                'test_phase3_5_rerun','test_phase6_8_rerun','test_phase15_16_19_rerun','test_phase17_rerun','test_phase14_rerun'],
    'molecular': ['test_phase10_13_rerun', 'test_phase15_checkpoint', 'test_phase7_resume', 'test_phase5_plot',
                  'test_phase15_16_19_rerun.RerunTests.test_real_torch_training_checkpoint_roundtrip_and_no_overwrite'],
}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scientific-profile',choices=['primary','molecular'],default='primary',
                        help='Explicit interpreter for general scientific tests; recorded, never silent fallback')
    args=parser.parse_args()
    before={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py')}
    out=ROOT/'results_rerun_audit'/time.strftime('%Y%m%dT%H%M%S')
    out.mkdir(parents=True,exist_ok=False)
    records=[]
    for profile,modules in GROUPS.items():
        actual_profile=args.scientific_profile if profile=='primary' else 'molecular'
        py=PRIMARY if actual_profile=='primary' else MOLECULAR
        env=os.environ.copy()
        env.update(PYTHONIOENCODING='utf-8',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',OPENMM_CPU_THREADS='2')
        if actual_profile=='primary': env['PYTHONPATH']=os.pathsep.join([str(DEPS),str(EXTRA)])
        else:
            env.pop('PYTHONPATH',None)
            env['PATH']=str(py.parent/'Library/bin')+os.pathsep+str(py.parent)+os.pathsep+env.get('PATH','')
        cmd=[str(py),'-X','faulthandler','-m','unittest','-v']+modules
        started=time.monotonic()
        with (out/f'{profile}.log').open('w',encoding='utf-8') as log:
            proc=subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
        rec=dict(profile=profile,actual_runtime_profile=actual_profile,command=cmd,exit_code=proc.returncode,elapsed_s=time.monotonic()-started)
        log_text=(out/f'{profile}.log').read_text(encoding='utf-8',errors='replace')
        count=re.search(r'Ran (\d+) tests? in',log_text)
        skipped=re.search(r'(?:OK|FAILED) \([^\n)]*skipped=(\d+)',log_text)
        rec['tests_run']=int(count.group(1)) if count else None
        rec['skipped']=int(skipped.group(1)) if skipped else 0
        rec['passing_non_skipped']=(rec['tests_run']-rec['skipped']
            if proc.returncode==0 and count else None)
        records.append(rec)
        print(json.dumps(rec),flush=True)
    after={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py')}
    changed=sorted(k for k in before.keys() | after.keys() if before.get(k)!=after.get(k))
    record=dict(groups=records,source_sha256=after,source_sha256_before=before,
                sources_changed_during_tests=changed,
                final_source_coverage_claimed=not changed and all(r['exit_code']==0 for r in records))
    (out/'regression_status.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(f'Evidence: {out}')
    return int(any(r['exit_code']!=0 for r in records))

if __name__=='__main__': raise SystemExit(main())
