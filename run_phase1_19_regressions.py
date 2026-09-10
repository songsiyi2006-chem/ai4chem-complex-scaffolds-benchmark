"""Run cheap regressions in suitable isolated interpreters and retain evidence."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
from rerun_phase1_19_campaign import PRIMARY,MOLECULAR,DEPS,EXTRA

ROOT=Path(__file__).resolve().parent
GROUPS={
    'primary': ['test_phase1_19_campaign','test_phase1_5_audit','test_phase6_10_audit',
                'test_phase11_13_audit','test_phase14_17_audit','test_phase18_19_audit',
                'test_phase3_5_rerun','test_phase6_8_rerun','test_phase15_16_19_rerun','test_phase17_rerun','test_phase14_rerun'],
    'molecular': ['test_phase10_13_rerun'],
}

def main():
    before={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py')}
    out=ROOT/'results_rerun_audit'/time.strftime('%Y%m%dT%H%M%S')
    out.mkdir(parents=True,exist_ok=False)
    records=[]
    for profile,modules in GROUPS.items():
        py=PRIMARY if profile=='primary' else MOLECULAR
        env=os.environ.copy()
        env.update(PYTHONIOENCODING='utf-8',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',OPENMM_CPU_THREADS='2')
        if profile=='primary': env['PYTHONPATH']=os.pathsep.join([str(DEPS),str(EXTRA)])
        else:
            env.pop('PYTHONPATH',None)
            env['PATH']=str(py.parent/'Library/bin')+os.pathsep+str(py.parent)+os.pathsep+env.get('PATH','')
        cmd=[str(py),'-X','faulthandler','-m','unittest','-v']+modules
        started=time.monotonic()
        with (out/f'{profile}.log').open('w',encoding='utf-8') as log:
            proc=subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
        rec=dict(profile=profile,command=cmd,exit_code=proc.returncode,elapsed_s=time.monotonic()-started)
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
