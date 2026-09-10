"""One-shot local calculation action for Windows Task Scheduler (no recurrence).

Prepare a job manifest, register it separately with the user's permission, then
execute --job. No credentials, elevation, shutdown, or unrelated processes.
"""
import argparse
import contextlib
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parent
PYTHON = Path('C:/Users/HUIWEI/miniconda3/envs/phase2ff/python.exe')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(data, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    temp.replace(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--prepare', choices=['probe', 'phase', 'resume7', 'regressions'])
    ap.add_argument('--phase', type=int, choices=range(1, 20))
    ap.add_argument('--extra', nargs=argparse.REMAINDER, default=[])
    ap.add_argument('--job', type=Path)
    args = ap.parse_args()
    if args.prepare:
        if args.prepare == 'phase' and args.phase is None:
            ap.error('--phase is required')
        stamp = time.strftime('%Y%m%dT%H%M%S')
        folder = ROOT/'results_rerun_audit'/'os_jobs'/f'{stamp}_{args.prepare}_{args.phase or 0}'
        folder.mkdir(parents=True, exist_ok=False)
        record = dict(kind=args.prepare, phase=args.phase, extra=args.extra,
                      worker_sha256=sha(__file__), root=str(ROOT),
                      task_name=f'GXNU_PhaseAudit_{stamp}_{args.prepare}_{args.phase or 0}',
                      source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                      scope='one-shot, same-user, no elevation or password, no startup trigger')
        save(folder/'job.json', record)
        print(str(folder/'job.json'))
        return 0
    if args.job is None:
        ap.error('--prepare or --job required')
    path = args.job.resolve()
    if not path.is_relative_to(ROOT/'results_rerun_audit'/'os_jobs'):
        raise ValueError('Job manifest outside this campaign')
    rec = json.loads(path.read_text(encoding='utf-8'))
    if rec['worker_sha256'] != sha(__file__) or Path(rec['root']) != ROOT:
        raise ValueError('Job source/workspace mismatch')
    status_path = path.parent/'status.json'
    if status_path.exists():
        raise FileExistsError('Job already dispatched; inspect its result before preparing another')
    env = os.environ
    env['PATH'] = str(PYTHON.parent/'Library/bin')+os.pathsep+str(PYTHON.parent)+os.pathsep+env.get('PATH', '')
    env.pop('PYTHONPATH', None)
    env.update({key: '2' for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS',
               'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'OPENMM_CPU_THREADS')})
    env.update(PYTHONIOENCODING='utf-8', PYTHONUNBUFFERED='1')
    status = dict(status='running', pid=os.getpid(), parent_pid=os.getppid(),
                  started=time.strftime('%Y-%m-%d %H:%M:%S'), job=rec,
                  interpreter=sys.executable)
    save(status_path, status)
    code = 1
    with (path.parent/'run.log').open('w', encoding='utf-8', buffering=1) as log:
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            try:
                if rec['kind'] == 'probe':
                    print('One-shot scheduler action started', flush=True)
                    time.sleep(2)
                    print('One-shot scheduler action completed', flush=True)
                    code = 0
                elif rec['kind'] == 'phase':
                    from rerun_phase1_19_campaign import run_one
                    from record_interrupted_campaign import process_state
                    campaign = ROOT.parent/'phase1-19-rerun-20260910'
                    for old in (campaign/f'phase{rec["phase"]:02d}').glob('*/rerun_status.json'):
                        previous = json.loads(old.read_text(encoding='utf-8'))
                        if previous.get('status') == 'running':
                            raise RuntimeError(f'Existing unresolved attempt: {old}')
                    result = run_one(rec['phase'], campaign, rec['extra'], profile='molecular')
                    code = result['exit_code']
                    status['attempt_result'] = result
                elif rec['kind'] == 'resume7':
                    from resume_phase7_rerun import main as resume7
                    code = resume7([str(ROOT.parent/'phase1-19-rerun-20260910/phase07/20260910T181009'), '--execute'])
                elif rec['kind'] == 'regressions':
                    code = subprocess.call([str(PYTHON), str(ROOT/'run_phase1_19_regressions.py'),
                                            '--scientific-profile', 'molecular'], cwd=ROOT,
                                           stdout=log, stderr=subprocess.STDOUT)
                else:
                    raise ValueError('Unknown job kind')
            except BaseException:
                traceback.print_exc()
                status['error'] = traceback.format_exc()
            finally:
                status.update(status='process_completed' if code == 0 else 'process_failed',
                              exit_code=code, finished=time.strftime('%Y-%m-%d %H:%M:%S'),
                              cleanup_due=True, scientific_acceptance='not established by process exit code')
                save(status_path, status)
    return code


if __name__ == '__main__':
    raise SystemExit(main())
