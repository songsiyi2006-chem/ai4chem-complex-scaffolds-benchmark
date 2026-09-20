"""Independent artifact and optional official OT-2 simulator checks.

python validate_phase24_outputs.py
# In a separate environment with opentrons==8.8.2:
python validate_phase24_outputs.py --opentrons
"""
import argparse
import ast
import csv
import hashlib
import importlib.metadata
import json
import sqlite3
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value,indent=2,ensure_ascii=True,allow_nan=False)+'\n',
                    encoding='utf-8',newline='\n')


def check_outputs(root):
    res=root/'results_phase24'
    manifest=json.loads((res/'manifest.json').read_text(encoding='utf-8'))
    for name,expected in manifest['outputs'].items():
        assert digest(root/name)==expected, f'Output hash mismatch: {name}'
    assert digest(root/'run_ai_hts_platform_pilot.py')==manifest['script_sha256']
    metrics=json.loads((res/'metrics.json').read_text())
    data=json.loads((res/'measurements.json').read_text())
    assert len(data)==metrics['measurements']==312
    assert len(set(r['condition_id'] for r in data))==312
    assert all(r['provenance']=='simulated' for r in data)
    assert sum(r['qc_pass'] for r in data)==metrics['qc_pass_count']
    assert sum(r['qc_pass'] and r['yield_pct']>=90 and r['ee_abs_pct']>=95 for r in data)==metrics['joint_hits']
    rows=list(csv.DictReader((res/'conditions.csv').open(encoding='utf-8')))
    assert len(rows)==9216
    assert all(r['wet_lab_approved']=='False' for r in rows)
    db=sqlite3.connect(f'file:{(res/"master_measurements.sqlite").as_posix()}?mode=ro',uri=True)
    assert db.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
    assert not db.execute('PRAGMA foreign_key_check').fetchall()
    assert db.execute('SELECT COUNT(*) FROM measurements').fetchone()[0]==312
    for r in data:
        stored=db.execute('SELECT raw_sha256 FROM measurements WHERE round=? AND well=?',
                          (r['round'],r['well'])).fetchone()[0]
        assert digest(root/r['raw_path'])==stored
    db.close()
    past=set()
    for r in range(4):
        plate=list(csv.DictReader((res/f'plate_round{r}.csv').open(encoding='utf-8')))
        assert len(plate)==(24 if r==0 else 96)
        ids={int(p['condition_id']) for p in plate}
        assert not past & ids
        for p in plate:
            assert p['well'][0] in 'ABCDEFGH' and 1<=int(p['well'][1:])<=12
            assert float(p['total_uL'])==sum(float(p[k]) for k in
                ['substrate_uL','catalyst_ligand_premix_uL','base_uL','makeup_solvent_uL'])==80
        if r:
            audit=json.loads((res/f'acquisition_round{r}.json').read_text())
            assert set(audit['training_ids'])==past
            assert set(audit['selected_ids'])==ids
        past |= ids
    result=dict(artifacts_passed=True,measurements=312,unique_conditions=312,
                raw_hashes_verified=312,manifest_hashes_verified=len(manifest['outputs']),
                no_future_labels_in_training_ids=True,sqlite_integrity='ok',
                physical_chemistry_validated=False)
    save(res/'artifact_validation.json',result)
    return result


def check_opentrons(root):
    from opentrons import simulate
    outputs=[]
    for path in [root/'results_phase24'/f'ot2_round{r}.py' for r in range(1,4)]:
        with path.open(encoding='utf-8') as stream:
            log,_=simulate.simulate(stream)
        payloads=[x['payload'] for x in log]
        picks=[p for p in payloads if p.get('text','').startswith('Picking up tip')]
        counts={key:sum(key in str(p.get('instrument','')) for p in picks) for key in ['P20','P300']}
        assert counts=={'P20':36,'P300':12},counts
        assert [p['celsius'] for p in payloads if 'celsius' in p]==[25]
        assert any('Deactivating Temperature' in p.get('text','') for p in payloads)
        # Robot simulation returns all tip columns and liquid operations without errors.
        outputs.append(dict(protocol=str(path.relative_to(root)),sha256=digest(path),
                            commands=len(log),tip_column_pickups=counts,
                            physical_tip_count={k:8*v for k,v in counts.items()},
                            temperatures_C=[25],completed=True))
    final=root/'output_hts_96well_screening.py'
    assert digest(final)==digest(root/'results_phase24'/'ot2_round3.py')
    result=dict(simulator_version=importlib.metadata.version('opentrons'),api_level='2.15',
                results=outputs,final_protocol_sha256=digest(final),
                chemistry_branch_simulated=False,physical_hardware_tested=False,
                note='Official software simulation; default calibration only; WATER branch.')
    save(root/'results_phase24'/'opentrons_validation.json',result)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path('.'))
    parser.add_argument('--opentrons',action='store_true')
    args=parser.parse_args()
    result=(check_opentrons if args.opentrons else check_outputs)(args.root.resolve())
    print(json.dumps(result,ensure_ascii=True))
