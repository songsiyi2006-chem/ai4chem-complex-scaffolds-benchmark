"""Re-render completed A/B/C outputs in a NEW directory; never mutate originals."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('attempt', type=Path)
    args = ap.parse_args()
    old = args.attempt.resolve()
    source = Path(__file__).resolve().parent
    data_path = old/'results_phase5/phase5_results.json'
    data = json.loads(data_path.read_text(encoding='utf-8'))
    if not all(data.get(k) for k in ('module_A', 'module_B', 'module_C')):
        raise ValueError('Require completed saved A/B/C data')
    out = old/('figure_recovery_'+time.strftime('%Y%m%dT%H%M%S'))
    out.mkdir(exist_ok=False)
    inputs = {str(data_path.relative_to(old)): sha(data_path)}
    for name in ('run_phase5_chemical_world_model.py', 'phase_audit.py'):
        shutil.copy2(source/name, out/name)
    shutil.copytree(old/'results_phase5', out/'results_phase5')
    for path in (old/'results_phase5').rglob('*'):
        if path.is_file():
            inputs[str(path.relative_to(old))] = sha(path)
    status = dict(status='running', scope='figures/reports only; original calculation unchanged',
                  source_sha256={p.name: sha(p) for p in out.glob('*.py')},
                  input_sha256=inputs, pid=os.getpid(), commands=[])
    record = out/'recovery_status.json'
    def save():
        record.write_text(json.dumps(status, indent=2)+'\n', encoding='utf-8')
    save()
    code = 1
    try:
        for stage in ('figures', 'reports'):
            cmd = [sys.executable, '-u', str(out/'run_phase5_chemical_world_model.py'), '--stage', stage]
            with (out/f'{stage}.log').open('w', encoding='utf-8') as log:
                code = subprocess.call(cmd, cwd=out, stdout=log, stderr=subprocess.STDOUT)
            status['commands'].append(dict(stage=stage, exit_code=code))
            if code:
                break
        if any(sha(old/rel) != digest for rel, digest in inputs.items()):
            raise RuntimeError('Original input changed during rendering')
        if code == 0:
            from PIL import Image
            figures = {}
            for p in (out/'figures_phase5').glob('*.png'):
                with Image.open(p) as img:
                    img.verify()
                with Image.open(p) as img:
                    if max(img.size) > 10000:
                        raise ValueError('Unexpected oversized figure')
                    figures[p.name] = dict(size=list(img.size), sha256=sha(p))
            if len(figures) != 3:
                raise ValueError('Missing figure')
            status['figures'] = figures
        status['status'] = 'completed' if code == 0 else 'failed'
        status['scientific_acceptance'] = 'not established: transition-state and kinetic-model limitations remain'
        return code
    finally:
        if status['status'] == 'running':
            status['status'] = 'failed'
        save()
        print(out, flush=True)


if __name__ == '__main__':
    raise SystemExit(main())
