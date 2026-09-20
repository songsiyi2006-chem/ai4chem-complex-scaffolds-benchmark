"""Resume audited Phase7 campaign; default is read-only planning, never --force.

Execution requires --execute and a main-agent allocated compute slot. Original
snapshot and rerun_status.json are immutable. Unknown caches fail closed.
"""
import argparse
import ctypes
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
from contextlib import contextmanager

SCRIPT = 'run_phase7_strong_correlation_wall.py'
SOURCE_HASH = '1141e04dc50ecdc55580416f11630187c85d09263d450570091e98bef8494e0d'
INPUT_HASH = 'e52cd7fbc8028bfbdf97841b04c1c87ec83b56c75097a894ec82e046e0db8c78'
SEEDS = {
    'R1.40': '1b6d6ad0b2a028ae954116d96678a394825fbebaa92af7a0d89aac50b713a0d5',
    'R1.50': '04fd7d66fad33a2ef9c3dcd1c63d6032f37a333dbd4a08767af60901f45eba37',
    'R1.60': '40882c5690183ebb6b4d1e869a51bc1515ef5889df6eda59058100e2b9cb9d12',
}
GRID = [1.4, 1.5, 1.6, 1.75, 1.9, 2.05, 2.2, 2.4, 2.6, 2.9, 3.2]
CUBES = [1.5, 2.2, 3.2]
BASES = ['def2-svp', '6-31g', 'sto-3g']
METHODS = ('rhf', 'uhf', 'rks', 'uks', 'uhf_triplet', 'uks_triplet')
CHEM = Path('C:/Users/HUIWEI/miniconda3/envs/phase2ff/python.exe')
QC = Path('C:/Users/HUIWEI/miniconda3/envs/phase7/python.exe')


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def checked_path(root, relative):
    path = (root / relative).resolve()
    require(path.is_relative_to(root.resolve()), f'Path escapes attempt: {relative}')
    return path


def checkpoint(stage, radius):
    return f'results_phase7/checkpoints/{stage}/R{radius:.2f}.json'


def valid_point(data, stage, radius, basis):
    require(data.get('R') == radius and data.get('tag') == f'R{radius:.2f}', 'Point identity mismatch')
    require(data.get('basis', '').lower() == basis.lower(), 'Basis mismatch')
    if stage == '7A':
        theories = data.get('theory', {})
        require(set(theories) == set(METHODS), 'Require all six methods')
        for name in METHODS:
            rec = theories[name]
            require(rec.get('converged') is True and finite(rec.get('energy_eh')), f'Invalid {name}')
            require(rec.get('basis', '').lower() == basis.lower(), f'{name} basis mismatch')
            require(rec.get('theory') == name, f'{name} identity mismatch')
            if name.startswith('u'):
                require(finite(rec.get('s2')) and finite(rec.get('ds2')), f'{name} spin missing')
                shift = 2 if name.endswith('triplet') else 0
                require(abs(rec['ds2'] - (rec['s2'] - shift)) < 1e-9, 'Spin convention mismatch')
    else:
        require(data.get('converged') is True, 'CASSCF unconverged')
        for key in ('e_casscf_eh', 'e_casscf_direct', 'noon_active_trace', 'y_luno', 'y_yamaguchi'):
            require(finite(data.get(key)), f'Missing/nonfinite {key}')
        noon = data.get('noon', [])
        require(len(noon) == 2 and all(finite(v) for v in noon), 'NOON missing/nonfinite')


def point_files(root, stage, radius, basis):
    rel = checkpoint(stage, radius)
    data = read(root / rel)
    valid_point(data, stage, radius, basis)
    files = {rel: digest(root / rel)}
    if stage == '7B' and radius in CUBES:
        names = data.get('cubes', [])
        require(len(names) >= 2, 'Required orbital cubes missing')
        for name in names:
            require(Path(name).name == name and name.endswith('.cube'), 'Invalid cube name')
            rel_cube = f'results_phase7/cubes/R{radius:.2f}/{name}'
            path = checked_path(root, rel_cube)
            require(path.stat().st_size > 0, 'Empty cube')
            files[rel_cube] = digest(path)
    return files


def verify_hashes(root, hashes):
    for rel, expected in hashes.items():
        require(digest(checked_path(root, rel)) == expected, f'Hash mismatch: {rel}')


def preflight(root):
    root = root.resolve()
    original = read(root / 'rerun_status.json')
    require(original.get('phase') == 7 and original.get('started') == '20260910T181009', 'Wrong campaign')
    require(original.get('status') == 'process_interrupted', 'Original attempt must remain interrupted')
    require(original['source_sha256'][SCRIPT] == SOURCE_HASH, 'Unexpected recorded source')
    verify_hashes(root, {SCRIPT: SOURCE_HASH, 'results_phase4/reactant_3d.mol': INPUT_HASH})
    inputs = {SCRIPT: SOURCE_HASH, 'results_phase4/reactant_3d.mol': INPUT_HASH,
              'rerun_status.json': digest(root / 'rerun_status.json')}
    scan_rel = 'results_phase7/checkpoints/G/scan.json'
    scan = read(root / scan_rel)['points']
    require([p['R'] for p in scan] == GRID, 'Default scan grid mismatch')
    inputs[scan_rel] = digest(root / scan_rel)
    for pt in scan:
        radius = pt['R']
        g_rel, ai_rel = checkpoint('G', radius), checkpoint('AI', radius)
        g, ai = read(root / g_rel), read(root / ai_rel)
        require(g == pt and g.get('converged') is True, 'Geometry/scan mismatch or unconverged')
        require(finite(g.get('bond_actual')) and abs(g['bond_actual'] - radius) < 1e-6, 'Constraint mismatch')
        require(ai.get('R') == radius and all(finite(ai.get(k)) for k in
                ('e_xtb_eV', 'e_mace_eV', 'e_ani_eV')), 'AI incomplete')
        xyz = checked_path(root, pt['xyz'])
        rows = xyz.read_text(encoding='utf-8').splitlines()
        require(int(rows[0]) == len(rows[2:]) and int(rows[0]) > 0, 'Malformed XYZ')
        require(all(len(row.split()) == 4 and all(finite(float(v)) for v in row.split()[1:])
                    for row in rows[2:]), 'Nonfinite XYZ')
        for rel in (g_rel, ai_rel, pt['xyz']):
            inputs[rel] = digest(checked_path(root, rel))
    state_path = root / 'continuation_status.json'
    prior = read(state_path) if state_path.exists() else None
    accepted = {}
    for tag, expected in SEEDS.items():
        radius = float(tag[1:])
        rel = checkpoint('7A', radius)
        verify_hashes(root, {rel: expected})
        accepted[rel] = {'basis': 'def2-svp', 'files': point_files(root, '7A', radius, 'def2-svp')}
    if prior:
        require(prior.get('tool') == 'resume_phase7_rerun' and prior.get('cwd') == str(root), 'Foreign continuation')
        require(prior['input_sha256'] == inputs, 'Continuation inputs changed')
        for rel, rec in prior['accepted'].items():
            verify_hashes(root, rec['files'])
            accepted[rel] = rec
    pending = []
    for stage in ('7A', '7B'):
        for i, radius in enumerate(GRID):
            rel = checkpoint(stage, radius)
            if rel in accepted:
                rec = accepted[rel]
                require(rec['basis'] in BASES, 'Unknown accepted basis')
                require(point_files(root, stage, radius, rec['basis']) == rec['files'], 'Accepted files changed')
            else:
                require(not (root / rel).exists(), f'Unregistered checkpoint preserved; review required: {rel}')
                pending.append((stage, i))
    return inputs, accepted, pending, prior


def save(path, data):
    temp = path.with_suffix('.json.tmp')
    temp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    temp.replace(path)


def available_gib():
    require(os.name == 'nt', 'Memory gate requires Windows')
    class MemoryStatus(ctypes.Structure):
        _fields_ = [('length', ctypes.c_ulong), ('load', ctypes.c_ulong)] + [
            (name, ctypes.c_ulonglong) for name in ('total_phys', 'avail_phys', 'total_page',
             'avail_page', 'total_virtual', 'avail_virtual', 'avail_extended')]
    status = MemoryStatus()
    status.length = ctypes.sizeof(status)
    require(ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)), 'Memory query failed')
    return status.avail_phys / 1024**3


def wait_memory(log):
    while True:
        free = available_gib()
        log.write(f'MEMORY available={free:.3f} GiB required=2.500 GiB\n'); log.flush()
        if free >= 2.5:
            return
        time.sleep(30)


@contextmanager
def execution_lock(root):
    # OS releases the byte lock after a crash; never delete another run's lock.
    import msvcrt
    with (root / 'phase7_continuation.lock').open('a+b') as handle:
        handle.seek(0, 2)
        if handle.tell() == 0:
            handle.write(b'0'); handle.flush()
        handle.seek(0)
        msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
        try:
            yield
        finally:
            handle.seek(0)
            msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)


def command(root, stage, index=None, basis='def2-svp', skip=0):
    py = QC if stage in ('7A', '7B') else CHEM
    cmd = [str(py), '-X', 'faulthandler', '-u', str(root / SCRIPT)]
    if stage == 'merge,fig':
        return cmd + ['--stage', stage, '--threads', '2']
    cmd += ['--worker', stage, '--point', str(index), '--basis', basis,
            '--mem', '2 GB', '--threads', '2']
    if stage == '7B':
        cmd += ['--skip-tiers', str(skip), '--cube-r', str(GRID[index] if GRID[index] in CUBES else -1.0)]
    return cmd


def launch(root, cmd, state, log, qc=True):
    verify_hashes(root, state['input_sha256'])
    for rec in state['accepted'].values():
        verify_hashes(root, rec['files'])
    if qc:
        wait_memory(log)
    env = os.environ.copy()
    prefix = Path(cmd[0]).parent
    env['PATH'] = os.pathsep.join([str(prefix / 'Library/bin'), str(prefix), env.get('PATH', '')])
    env.pop('PYTHONPATH', None)
    env.update({key: '2' for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                                    'OPENMM_CPU_THREADS', 'NUMEXPR_NUM_THREADS')})
    env.update(PYTHONUNBUFFERED='1', PYTHONIOENCODING='utf-8')
    record = {'command': cmd, 'started': time.time(), 'status': 'dispatching'}
    state['commands'].append(record)
    save(root / 'continuation_status.json', state)
    log.write('\nCOMMAND ' + json.dumps(cmd) + '\n'); log.flush()
    with subprocess.Popen(cmd, cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT) as proc:
        record.update(pid=proc.pid, status='running')
        save(root / 'continuation_status.json', state)
        rc = proc.wait()
    record.update(exit_code=rc, finished=time.time(), status='exited')
    save(root / 'continuation_status.json', state)
    return rc


def run_pending(root, state, pending, log):
    failed = []
    for stage, index in pending:
        radius = GRID[index]
        rel = checkpoint(stage, radius)
        ok = False
        for basis in BASES:
            for skip in ((0, 1) if stage == '7B' else (0,)):
                rc = launch(root, command(root, stage, index, basis, skip), state, log)
                try:
                    require(rc == 0, f'Worker exit {rc}')
                    files = point_files(root, stage, radius, basis)
                except (ValueError, OSError, KeyError, TypeError) as exc:
                    log.write(f'RETRY {stage} point={index} basis={basis} skip={skip}: {exc}\n')
                    # Preserve each failed checkpoint before the original worker overwrites it.
                    if (root / rel).exists():
                        evidence = root / state['log_directory'] / f'{stage}_{index}_{basis}_{skip}.json'
                        evidence.write_bytes((root / rel).read_bytes())
                    continue
                state['accepted'][rel] = {'basis': basis, 'files': files}
                save(root / 'continuation_status.json', state)
                ok = True
                break
            if ok:
                break
        if not ok:
            failed.append(rel)
    state['failed_points'] = failed
    # Run original merge/figure code even on scientific failure, preserving its diagnosis.
    rc = launch(root, command(root, 'merge,fig'), state, log, qc=False)
    master = read(root / 'results_phase7/phase7_results.json')
    return 0 if not failed and rc == 0 and master.get('all_stages_ok') is True else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('attempt', type=Path)
    ap.add_argument('--execute', action='store_true')
    args = ap.parse_args(argv)
    root = args.attempt.resolve()
    inputs, accepted, pending, prior = preflight(root)
    if not args.execute:
        print(json.dumps({'mode': 'read_only_plan', 'cwd': str(root), 'accepted': accepted,
                          'pending': pending, 'final_stage': 'merge,fig', 'input_sha256': inputs}, indent=2))
        return 0
    require(QC.is_file() and CHEM.is_file(), 'Required Python runtime missing')
    with execution_lock(root):
        inputs, accepted, pending, prior = preflight(root)
        stamp = time.strftime('%Y%m%dT%H%M%S') + f'_{os.getpid()}'
        log_dir = root / 'phase7_continuations' / stamp
        log_dir.mkdir(parents=True, exist_ok=False)
        if prior:
            save(log_dir / 'previous_continuation_status.json', prior)
        state = dict(tool='resume_phase7_rerun', phase=7, cwd=str(root), status='running',
                     started=time.time(), pid=os.getpid(), source_sha256={str(Path(__file__).resolve()): digest(__file__)},
                     input_sha256=inputs, accepted=accepted, commands=[], pending=pending,
                     log_directory=str(log_dir.relative_to(root)),
                     launch_command=[sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
                     protocol={'threads': 2, 'mem': '2 GB', 'minimum_free_gib': 2.5, 'basis_chain': BASES},
                     scientific_acceptance='pending output review', scope='same-campaign continuation')
        save(root / 'continuation_status.json', state)
        with (log_dir / 'run.log').open('w', encoding='utf-8') as log:
            try:
                code = run_pending(root, state, pending, log)
                state.update(status='process_completed' if code == 0 else 'process_failed', exit_code=code)
            except BaseException as exc:
                state.update(status='continuation_interrupted', error=repr(exc))
                raise
            finally:
                state['finished'] = time.time()
                save(root / 'continuation_status.json', state)
        return code


if __name__ == '__main__':
    raise SystemExit(main())
