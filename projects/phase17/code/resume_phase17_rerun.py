"""Recompute corrected Phase17 downstream stages from this campaign's SCF data.

No historical calculation is imported; explicit source paths and hashes are
recorded. Figures stage includes bonding and avoids computing that stage twice.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from rerun_phase1_19_campaign import ROOT, PRIMARY, DEPS, EXTRA, save, runtime_metadata


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('atomic_attempt', type=Path)
    args = ap.parse_args()
    parent = args.atomic_attempt.resolve()
    atomic = json.loads((parent/'results_phase17/phase17_results.json').read_text(encoding='utf-8'))
    if not all(atomic['atoms'][k]['scf']['converged'] for k in ('Am3+', 'Eu3+', 'Am3+NR', 'Eu3+NR')):
        raise RuntimeError('Cannot resume unconverged atomic data')
    stamp = time.strftime('%Y%m%dT%H%M%S')
    dest = parent.parent/stamp
    dest.mkdir(exist_ok=False)
    hashes = {}
    for p in ROOT.glob('*.py'):
        if p.name.startswith('output_'):
            continue
        shutil.copy2(p, dest/p.name)
        hashes[p.name] = hashlib.sha256((dest/p.name).read_bytes()).hexdigest()
    prerequisites = []
    (dest/'results_phase17').mkdir()
    for name in ('phase17_results.json', 'phase17_radial_data.npz'):
        src = parent/'results_phase17'/name
        dst = dest/'results_phase17'/name
        shutil.copy2(src, dst)
        prerequisites.append(dict(source=str(src), destination=str(dst.relative_to(dest)),
                                  sha256=hashlib.sha256(dst.read_bytes()).hexdigest()))
    env = os.environ.copy()
    env.update(PYTHONUNBUFFERED='1', PYTHONIOENCODING='utf-8', OMP_NUM_THREADS='2',
               OPENBLAS_NUM_THREADS='2', MKL_NUM_THREADS='2', NUMEXPR_NUM_THREADS='2',
               PYTHONPATH=os.pathsep.join((str(DEPS), str(EXTRA))))
    commands = [[str(PRIMARY), '-X', 'faulthandler', '-u',
                 'run_phase17_relativistic_actinide_quantum.py', '--stage', stage]
                for stage in ('breit', 'multiplet', 'figures')]
    rec = dict(phase=17, started=stamp, status='running', cwd=str(dest),
               command=commands, source_sha256=hashes, prerequisite_artifacts=prerequisites,
               historical_cache_copied=False, fresh_atomic_parent=str(parent),
               scope='same-campaign full downstream recomputation; validated fresh atomic checkpoint',
               scientific_acceptance='pending output review', runtime=runtime_metadata(PRIMARY, env),
               stages=[])
    save(dest/'rerun_status.json', rec)
    t0 = time.monotonic()
    print(f'PHASE17 DOWNSTREAM START {dest}', flush=True)
    code = 0
    with (dest/'run.log').open('w', encoding='utf-8') as log:
        for cmd in commands:
            log.write(f'\nRECOMPUTE STAGE {cmd[-1]}\n'); log.flush()
            stage_start = time.monotonic()
            proc = subprocess.Popen(cmd, cwd=dest, env=env, stdout=log, stderr=subprocess.STDOUT)
            rec.update(pid=proc.pid, current_stage=cmd[-1]); save(dest/'rerun_status.json', rec)
            code = proc.wait()
            rec['stages'].append(dict(stage=cmd[-1], exit_code=code,
                                      elapsed_s=time.monotonic()-stage_start))
            save(dest/'rerun_status.json', rec)
            if code:
                break
    rec.update(status='process_completed' if code == 0 else 'process_failed',
               exit_code=code, elapsed_s=time.monotonic()-t0)
    save(dest/'rerun_status.json', rec)
    print(f'PHASE17 DOWNSTREAM EXIT {code} {rec["elapsed_s"]:.1f}s', flush=True)
    raise SystemExit(code)


if __name__ == '__main__':
    main()
