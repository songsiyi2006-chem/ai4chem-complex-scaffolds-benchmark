"""Prepare fresh MD from this campaign's completed upstream stages only.

No trajectory/velocity continuation is claimed. This command does not execute
MD or overwrite the old attempt. Review restart_provenance.json before launch.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import shutil
import time

ROOT = Path(__file__).resolve().parent
SCRIPT = 'run_phase3_complex_dynamics.py'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stages(path):
    wanted = {'stage1a_ingest', 'stage2_docking', 'stage3_ff_duality'}
    # Explicitly reviewed display-only change: the original log mislabeled
    # Sage as GAFF2. No expression, coordinate, energy or cached record changes.
    text = path.read_text(encoding='utf-8').replace('GAFF2 single point:', 'Sage 2.1 single point:')
    return {n.name: ast.dump(n) for n in ast.parse(text).body
            if isinstance(n, ast.FunctionDef) and n.name in wanted}


def prepare(old):
    old = old.resolve()
    state = json.loads((old/'rerun_status.json').read_text(encoding='utf-8'))
    if state.get('phase') != 3 or state.get('status') != 'process_interrupted':
        raise ValueError('Require an explicitly recorded interrupted Phase3 attempt')
    if state.get('source_sha256', {}).get(SCRIPT) != digest(old/SCRIPT):
        raise ValueError('Original executed source hash mismatch')
    if stages(old/SCRIPT) != stages(ROOT/SCRIPT):
        raise ValueError('Upstream stage implementations changed; review before reuse')
    source = old/'results_phase3'
    required = ['stage1a.json', 'stage2.json', 'stage3.json', '7RPZ_fixed_protein.pdb',
                'GDP_withH.pdb', '7RPZ_MG.pdb', 'T04_pose1.sdf', 'T04_pose1_Hrelaxed.sdf']
    if any(not (source/name).is_file() for name in required):
        raise ValueError('Required preparation evidence missing')
    # Explicit allowlist: never copy a partial DCD or stage4/5 result.
    names = required + ['7RPZ.pdb', '7RPZ_GDP.pdb', '7RPZ_protein_raw.pdb',
        'native_ligand.sdf', 'receptor_merged.pdb', 'receptor.pdbqt',
        'T04_docked_poses.pdbqt', 'T04_docking_input.pdbqt', 'T04_pose1.pdb']
    hashes = {name: digest(source/name) for name in names if (source/name).is_file()}
    dest = old.parent/('prepared_'+time.strftime('%Y%m%dT%H%M%S'))
    dest.mkdir(exist_ok=False)
    (dest/'results_phase3').mkdir()
    shutil.copy2(ROOT/SCRIPT, dest/SCRIPT)
    for name in hashes:
        shutil.copy2(source/name, dest/'results_phase3'/name)
        if digest(dest/'results_phase3'/name) != hashes[name]:
            raise RuntimeError('Copy verification failed')
    record = dict(status='prepared_not_started', reused_stages=[1, 2, 3],
                  reused_file_sha256=hashes, original_attempt=old.name,
                  original_source_sha256=digest(old/SCRIPT),
                  approved_display_only_difference='GAFF2 single point -> Sage 2.1 single point log label',
                  restart_source_sha256=digest(dest/SCRIPT),
                  scope='same-campaign upstream reuse; new full MD, not trajectory continuation',
                  command=['python', SCRIPT],
                  scientific_acceptance='not established; unchanged stage bodies do not certify upstream physics')
    (dest/'restart_provenance.json').write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8')
    print(dest)
    return dest


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('attempt', type=Path)
    prepare(ap.parse_args().attempt)
