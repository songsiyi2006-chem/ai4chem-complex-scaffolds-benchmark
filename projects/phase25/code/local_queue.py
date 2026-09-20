"""Finite local calculation queue; no scheduler, remote service, or auto-publish.

Can adopt an already running job by exact module/tag/cwd matching. Subsequent
jobs run serially, with caller-specified memory and threads. Results require a
separate human/agent review before publication or scientific acceptance.
"""
import argparse
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import subprocess
import time
import psutil


def now():return datetime.now(timezone.utc).isoformat()


def active_pilot(root,tag):
    found=[]
    for proc in psutil.process_iter(['pid','cmdline','cwd']):
        try:
            cmd=proc.info['cmdline'] or []
            if 'phase25_27.stationary_pilot' in cmd and tag in cmd and Path(proc.info['cwd']).resolve()==root:
                found.append(proc.info['pid'])
        except (psutil.Error,TypeError,OSError):continue
    return found


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);args=ap.parse_args()
    cfgpath=Path(args.config).resolve();cfg=json.loads(cfgpath.read_text(encoding='utf8'))
    root=Path(cfg['root']).resolve();out=cfgpath.parent;statusfile=out/'status.json'
    lock=out/'worker.lock'
    if lock.exists():
        oldpid=int(lock.read_text())
        if psutil.pid_exists(oldpid):raise RuntimeError('An existing worker may still own this queue')
        lock.unlink()
    with lock.open('x') as f:f.write(str(os.getpid()))
    status=dict(worker_pid=os.getpid(),started_utc=now(),state='RUNNING',jobs=[],auto_publish=False)
    def save():
        status['updated_utc']=now();tmp=statusfile.with_suffix('.tmp');tmp.write_text(json.dumps(status,indent=2),encoding='utf8');tmp.replace(statusfile)
    save();env=os.environ.copy();env['PYTHONUTF8']='1';env['PYTHONPATH']=cfg['dependency_path'];env['OMP_NUM_THREADS']=str(cfg['threads'])
    terminal={'OPTIMIZED_FREQUENCY_PENDING','FREQUENCIES_COMPUTED_REVIEW_REQUIRED','SINGLE_POINT_COMPLETED'}
    try:
        for item in cfg['jobs']:
            tag=item['tag'];entry=dict(tag=tag,started_utc=now(),state='CHECKING');status['jobs'].append(entry);save()
            while pids:=active_pilot(root,tag):
                entry.update(state='ADOPTED_EXISTING_PROCESS',pids=pids);save();time.sleep(15)
            ledger=root/'results_phase25_27/phase25/stationary_pilot'/tag/'result.json'
            previous=json.loads(ledger.read_text(encoding='utf8')) if ledger.exists() else {}
            if previous.get('status') in terminal:
                entry.update(state='FINISHED_EXISTING',result_status=previous['status'],finished_utc=now());save();continue
            inp=root/item['input']
            if not inp.exists():
                entry.update(state='BLOCKED_MISSING_INPUT',finished_utc=now());save();continue
            command=[cfg['psi4_python'],'-X','utf8','-u','-m','phase25_27.stationary_pilot','--input',item['input'],
              '--tag',tag,'--stage',item['stage'],'--method',item['method'],'--threads',str(cfg['threads']),'--memory-mb',str(cfg['memory_mb'])]
            if ledger.exists():command+=['--resume']
            entry.update(state='STARTING_LOCAL_JOB',command=command);save()
            with (out/(tag+'.log')).open('ab') as log:
                process=subprocess.Popen(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT)
                entry.update(state='RUNNING_LOCAL_JOB',pid=process.pid);save();entry['exit_code']=process.wait()
            result=json.loads(ledger.read_text(encoding='utf8')) if ledger.exists() else {}
            entry.update(state='FINISHED' if result.get('status') in terminal else 'FAILED_REVIEW_REQUIRED',result_status=result.get('status'),finished_utc=now());save()
        status.update(state='FINISHED_REVIEW_REQUIRED',finished_utc=now());save()
    except Exception as exc:
        status.update(state='WORKER_FAILED',error=str(exc));save();raise
    finally:
        if lock.exists() and lock.read_text()==str(os.getpid()):lock.unlink()


if __name__=='__main__':main()
