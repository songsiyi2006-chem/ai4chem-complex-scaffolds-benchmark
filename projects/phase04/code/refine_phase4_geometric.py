"""Bounded alternative TS refinement on the SAME GFN2-xTB surface.

Requires the installed chem-ai4s extension's geomeTRIC. This is an independent
attempt, not a retroactive NEB convergence claim or a validated reaction rate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import time
from unittest.mock import patch
import numpy as np

EH_EV = 27.211386245988
BOHR_A = 0.52917721092


def atomic_units(energy_eV, forces_eV_A):
    return dict(energy=float(energy_eV / EH_EV),
                gradient=(-np.asarray(forces_eV_A) * BOHR_A / EH_EV).ravel())


def main():
    import geometric
    from geometric.engine import Engine
    from geometric.molecule import Molecule
    from geometric.optimize import run_optimizer
    from ase.io import read, write
    import run_phase4_reaction_mechanism as p
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('attempt', type=Path)
    ap.add_argument('--maxiter', type=int, default=200)
    ap.add_argument('--verify-from', type=Path,
                    help='Recheck saved converged geometry without repeating optimization')
    args = ap.parse_args()
    if not 1 <= args.maxiter <= 200:
        ap.error('--maxiter must be in 1..200; the fixed wall budget remains 240 seconds')
    attempt = args.attempt.resolve()
    source = attempt / 'results_phase4'
    metadata = json.loads((source / 'stage2.json').read_text())
    if metadata['fmax_target'] != .05 or not p.XTB_EXE:
        raise ValueError('Require original force target and installed xTB; no fallback')
    band = source / 'neb_final_path.xyz'
    band_hash = hashlib.sha256(band.read_bytes()).hexdigest()
    images = read(str(band), index=':')
    candidate = images[metadata['ts_image_index']]
    recovery = None
    if args.verify_from:
        recovery = args.verify_from.resolve()
        if recovery.parent != attempt:
            raise ValueError('Recovery must belong to this original attempt')
        old = json.loads((recovery / 'audit.json').read_text(encoding='utf-8'))
        if old['input_band_sha256'] != band_hash:
            raise ValueError('Recovery input band hash mismatch')
        if 'Converged! =D' not in (recovery / 'ts.log').read_text(encoding='utf-8', errors='replace'):
            raise ValueError('No recorded optimizer convergence')
        saved = read(str(recovery / 'refined_candidate.xyz'))
        if not np.array_equal(saved.numbers, candidate.numbers) or not np.isfinite(saved.positions).all():
            raise ValueError('Invalid saved atom identity or coordinates')
    out = attempt / ('tsopt_postaudit_' + time.strftime('%Y%m%dT%H%M%S'))
    out.mkdir(exist_ok=False)
    for src, name in ((Path(p.__file__), 'verified_source.py'),
                      (Path(__file__), Path(__file__).name)):
        shutil.copy2(src, out / name)
    write(str(out / 'initial.xyz'), candidate.copy())
    molecule = Molecule(str(out / 'initial.xyz'))
    engine = p.XTBWrap(charge=0, mult=1)
    numbers = candidate.get_atomic_numbers()
    start = time.monotonic()
    class XTBEngine(Engine):
        def calc_new(self, coords, dirname):
            energy, forces = engine._run(numbers, np.asarray(coords).reshape(-1, 3) * BOHR_A)
            result = atomic_units(energy, forces)
            if not np.isfinite(result['energy']) or not np.isfinite(result['gradient']).all():
                raise ValueError('Nonfinite xTB energy/gradient')
            return result
    adapter = XTBEngine(molecule)
    record = dict(status='running', method='geomeTRIC Cartesian TS / GFN2-xTB',
                  geometric_version=geometric.__version__, charge=0, multiplicity=1,
                  input_band_sha256=band_hash, initial_image=metadata['ts_image_index'],
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  scope='bounded_alternative_TS_refinement_not_default_NEB_rerun',
                  max_iterations=args.maxiter, gradient_force_target_eV_A=.05,
                  acceptance='not_established', irc_performed=False)
    if recovery:
        record.update(scope='saved_converged_geometry_gradient_and_hessian_reverification',
                      optimizer_repeated=False,
                      recovery_directory=str(recovery),
                      recovery_sha256={name: hashlib.sha256((recovery / name).read_bytes()).hexdigest()
                                       for name in ('audit.json', 'ts.log', 'refined_candidate.xyz')})
    path = out / 'audit.json'
    path.write_text(json.dumps(record, indent=2), encoding='utf-8')
    original_run = p.subprocess.run
    def bounded_run(*args, **kwargs):
        remaining = 240 - (time.monotonic() - start)
        if remaining <= 0:
            raise TimeoutError('Bounded TS allocation exhausted')
        kwargs['timeout'] = min(kwargs.get('timeout', remaining), remaining)
        return original_run(*args, **kwargs)
    try:
        with patch.object(p.subprocess, 'run', bounded_run):
            if recovery:
                final = saved.copy()
            else:
                result = run_optimizer(input=str(out / 'initial.xyz'), prefix=str(out / 'ts'),
                                       customengine=adapter, transition=True, coordsys='cart',
                                       hessian='first+last', maxiter=args.maxiter, trust=.01, tmax=.03)
                final = candidate.copy()
                final.positions = result.xyzs[-1]
            write(str(out / 'refined_candidate.xyz'), final)
            energy, forces = engine._run(numbers, final.positions)
            hessian = p.xtb_hessian(numbers, final.positions, charge=0)
        fmax = float(np.linalg.norm(forces, axis=1).max())
        frequencies = hessian['frequencies']
        (out / 'final_hessian.json').write_text(
            json.dumps(hessian, indent=2, default=lambda x: x.tolist(), allow_nan=False), encoding='utf-8')
        significant_imaginary = [f for f in frequencies if f < -20.]
        stationary_saddle = fmax <= .05 and len(significant_imaginary) == 1
        record.update(status='completed', optimizer_convergence_recorded=True, final_energy_eV=energy,
                      final_force_eV_A=fmax, significant_imaginary_cm=significant_imaginary,
                      stationary_saddle_checks_passed=stationary_saddle,
                      acceptance='stationary_saddle_only_reaction_connection_unverified' if stationary_saddle
                      else 'failed_stationary_saddle_checks', enzyme_or_rate_claims=False)
        if hashlib.sha256(band.read_bytes()).hexdigest() != band_hash:
            raise RuntimeError('Original band modified during optimization')
        record['input_hash_rechecked'] = True
    except BaseException as exc:
        record.update(status='failed', error=repr(exc), acceptance='not_established')
        raise
    finally:
        record.update(elapsed_s=time.monotonic() - start, energy_gradient_calls=engine.n_calls)
        path.write_text(json.dumps(record, indent=2), encoding='utf-8')
        print(json.dumps(dict(output=str(out), result=record), indent=2), flush=True)


if __name__ == '__main__':
    main()
