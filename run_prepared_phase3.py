"""Run a prepared, same-campaign Phase3 restart once after the existing queue.

No recurring trigger, old trajectory reuse, reduced sampling, code editing or
automatic Git publication. The queue dependency avoids competing memory gates.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

from prepare_phase3_restart import ROOT, SCRIPT, digest
from rerun_phase1_19_campaign import MOLECULAR, runtime_metadata, save
from resume_phase7_rerun import available_gib
from record_interrupted_campaign import process_state


def validate(attempt):
    attempt = Path(attempt).resolve()
    campaign = ROOT.parent / 'phase1-19-rerun-20260910'
    if attempt.parent != campaign / 'phase03':
        raise ValueError('Prepared attempt must be in this Phase3 campaign')
    record = json.loads((attempt / 'restart_provenance.json').read_text(encoding='utf-8'))
    if record['status'] != 'prepared_not_started':
        raise ValueError('Not an unused prepared attempt')
    if digest(attempt / SCRIPT) != record['restart_source_sha256']:
        raise ValueError('Prepared source hash mismatch')
    for name, expected in record['reused_file_sha256'].items():
        if Path(name).name != name or digest(attempt / 'results_phase3' / name) != expected:
            raise ValueError('Prepared input identity mismatch')
    if (attempt / 'rerun_status.json').exists():
        raise FileExistsError('Attempt already started; preserve it')
    return record


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('attempt', type=Path)
    ap.add_argument('--after-queue', type=Path, required=True)
    args = ap.parse_args()
    attempt = args.attempt.resolve()
    record = validate(attempt)
    dependency = args.after_queue.resolve()
    if not dependency.is_relative_to(ROOT / 'results_rerun_audit'):
        raise ValueError('Queue dependency outside this workspace')
    # Fail before claiming the attempt if the dependency cannot be read.
    json.loads(dependency.read_text(encoding='utf-8'))
    with (attempt / 'restart_launch_claim.json').open('x', encoding='utf-8') as fh:
        json.dump(dict(pid=os.getpid(), dependency=str(dependency)), fh)
    status_path = attempt / 'restart_launch_status.json'
    state = dict(status='waiting_for_existing_queue', pid=os.getpid(),
                 dependency=str(dependency), minimum_free_gib=3.5,
                 mode='one finite, one-shot calculation; no recurring trigger')
    save(status_path, state)
    try:
        while True:
            queue = json.loads(dependency.read_text(encoding='utf-8'))
            if queue['status'] == 'queue_finished':
                break
            if queue['status'] != 'running':
                raise RuntimeError('Existing queue did not finish normally; manual review required')
            if not process_state(queue['pid'])[0]:
                raise RuntimeError('Existing queue process is absent; completion is not inferred')
            time.sleep(5)
        state['status'] = 'waiting_for_memory'
        save(status_path, state)
        while available_gib() < 3.5:
            time.sleep(5)
        record = validate(attempt)
        for prior in attempt.parent.glob('*/rerun_status.json'):
            if json.loads(prior.read_text(encoding='utf-8')).get('status') == 'running':
                raise RuntimeError(f'Another unresolved Phase3 attempt: {prior}')
        env = os.environ.copy()
        env.pop('PYTHONPATH', None)
        env.pop('PYTHONHOME', None)
        env['PATH'] = str(MOLECULAR.parent / 'Library/bin') + os.pathsep + str(MOLECULAR.parent) + os.pathsep + env.get('PATH', '')
        env.update({k: '2' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                   'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'OPENMM_CPU_THREADS')})
        env.update(PYTHONUNBUFFERED='1', PYTHONIOENCODING='utf-8')
        command = [str(MOLECULAR), '-X', 'faulthandler', '-u', SCRIPT]
        run = dict(phase=3, started=time.strftime('%Y%m%dT%H%M%S'),
                   status='running', scope='default/full MD; same-campaign upstream reuse',
                   source_sha256={SCRIPT: record['restart_source_sha256']},
                   prerequisite_artifacts=record['reused_file_sha256'],
                   restart_provenance=record, historical_cache_copied=False,
                   scientific_acceptance='pending output review', command=command,
                   cwd=str(attempt), runtime=runtime_metadata(MOLECULAR, env))
        start = time.monotonic()
        save(attempt / 'rerun_status.json', run)
        with (attempt / 'run.log').open('x', encoding='utf-8') as log:
            try:
                process = subprocess.Popen(command, cwd=attempt, env=env,
                                           stdout=log, stderr=subprocess.STDOUT)
                run['pid'] = process.pid
                save(attempt / 'rerun_status.json', run)
                state.update(status='calculating', child_pid=process.pid)
                save(status_path, state)
                code = process.wait()
            except OSError as exc:
                run['launch_error'] = repr(exc)
                code = 127
        run.update(status='process_completed' if code == 0 else 'process_failed',
                   exit_code=code, elapsed_s=time.monotonic() - start)
        save(attempt / 'rerun_status.json', run)
        state.update(status=run['status'], exit_code=code,
                     scientific_acceptance='requires output review')
        return code
    except BaseException as exc:
        state.update(status='launcher_failed', error=repr(exc))
        raise
    finally:
        save(status_path, state)


if __name__ == '__main__':
    raise SystemExit(main())
