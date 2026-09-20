"""One fixed-band xTB force evaluation; no optimizer or TS acceptance."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import time
from unittest.mock import patch


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    from ase.io import read
    from ase.mep.neb import NEB
    import run_phase4_reaction_mechanism as p
    from diagnose_phase4_neb import replay_force_diagnostics
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('attempt', type=Path)
    args = ap.parse_args()
    attempt = args.attempt.resolve()
    source = attempt / 'results_phase4'
    inputs = [source / 'neb_final_path.xyz', source / 'stage2.json', attempt / 'rerun_status.json']
    hashes = {str(path): digest(path) for path in inputs}
    previous = json.loads(inputs[1].read_text())
    original = json.loads(inputs[2].read_text())
    if previous['fmax_target'] != .05 or original['phase'] != 4 or not p.XTB_EXE:
        raise ValueError('Require original Phase4 target and installed xTB')
    images = read(str(inputs[0]), index=':')
    if len(images) != previous['n_images'] or any(im.constraints or im.pbc.any() for im in images):
        raise ValueError('Unexpected band size or constrained/periodic input')
    engine = p.XTBWrap(charge=0, mult=1)
    for image in images:
        image.calc = p._PerAtomCalc(engine)
    neb = NEB(images, climb=True, k=.1, method='improvedtangent', parallel=False,
              remove_rotation_and_translation=False)
    out = attempt / ('force_postaudit_' + time.strftime('%Y%m%dT%H%M%S'))
    out.mkdir(exist_ok=False)
    for code in (Path(__file__), Path(p.__file__), Path(__file__).with_name('diagnose_phase4_neb.py')):
        shutil.copy2(code, out / code.name)
    meta = dict(status='running', scope='one_fixed_band_evaluation_no_optimizer',
                acceptance='no_TS_no_IRC', input_sha256=hashes,
                production_source_sha256=digest(Path(p.__file__)),
                xtb_executable=str(p.XTB_EXE), xtb_executable_sha256=digest(Path(p.XTB_EXE)),
                spring_eV_A2=.1, force_target_eV_A=.05)
    audit = out / 'audit.json'
    audit.write_text(json.dumps(meta, indent=2), encoding='utf-8')
    original_run = p.subprocess.run
    start = time.monotonic()
    def bounded_run(*args, **kwargs):
        remaining = 60 - (time.monotonic() - start)
        if remaining <= 0:
            raise TimeoutError('Fixed-band allocation exhausted')
        kwargs['timeout'] = min(kwargs.get('timeout', remaining), remaining)
        return original_run(*args, **kwargs)
    try:
        with patch.object(p.subprocess, 'run', bounded_run):
            projected = neb.get_forces()
        p._save_neb_force_diagnostics(out, neb, projected, engine, .05)
        path = out / 'neb_force_diagnostics.json'
        evidence = json.loads(path.read_text(encoding='utf-8'))
        all_raw = []
        for image in images:
            state = (image.get_atomic_numbers().tobytes(), image.get_positions().tobytes())
            if image.calc._state != state or image.calc._values is None:
                raise RuntimeError('Exact-state force cache missing')
            all_raw.append(image.calc._values[1].copy().tolist())
        evidence.update(all_image_raw_forces_eV_A=all_raw,
                        state='one_evaluation_of_saved_rounded_XYZ_no_geometry_updates',
                        probe_input_sha256=hashes)
        path.write_text(json.dumps(evidence, indent=2, allow_nan=False), encoding='utf-8')
        replay = replay_force_diagnostics(evidence)
        (out / 'force_replay.json').write_text(json.dumps(replay, indent=2), encoding='utf-8')
        if any(digest(path) != hashes[str(path)] for path in inputs):
            raise RuntimeError('Original inputs changed during probe')
        meta.update(status='completed', interpretation='evaluated_not_scientifically_accepted',
                    replay=replay, evidence_sha256=digest(out / 'neb_force_diagnostics.json'),
                    input_hashes_rechecked=True)
        print(json.dumps(dict(output=str(out), replay=replay), indent=2))
    except BaseException as exc:
        meta.update(status='failed', error=repr(exc))
        raise
    finally:
        meta.update(engine_calls=engine.n_calls, elapsed_s=time.monotonic() - start)
        audit.write_text(json.dumps(meta, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
