"""Recheck generated artifacts, record software validation and hash the release.

python -m phase25_27.release_audit
Does not submit quantum jobs or upgrade any scientific acceptance criterion.
"""
import csv
from datetime import datetime,timezone
import hashlib
import io
import json
from pathlib import Path
import platform
import subprocess
import sys
import unittest
import numpy
import scipy
from rdkit import rdBase
from .build_inputs import OUT,ROOT,dump
from . import test_analysis
from . import test_pilot


def main():
    # Check execution-time hashes before making a new publication manifest.
    execution_files_checked=0
    for ledger in sorted((OUT/'phase25/xtb_pilot').glob('*/results.json')):
        for row in json.loads(ledger.read_text(encoding='utf8')):
            if row['status'] in ['STARTED','HESSIAN_RUNNING']:
                raise ValueError('Finish the finite xTB campaign before release')
            for name,digest in row.get('files',{}).items():
                path=ledger.parent/row['name']/name
                if hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
                    raise ValueError(f'Execution-time xTB hash mismatch: {path}')
                execution_files_checked+=1
    for ledger in sorted((OUT/'phase25/stationary_pilot').glob('*/result.json')):
        row=json.loads(ledger.read_text(encoding='utf8'))
        if row['status']=='STARTED':continue
        for name,digest in row.get('files',{}).items():
            path=ledger.parent/name
            if hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
                raise ValueError(f'Execution-time DFT hash mismatch: {path}')
            execution_files_checked+=1
    suite=unittest.defaultTestLoader.loadTestsFromModule(test_analysis)
    suite.addTests(unittest.defaultTestLoader.loadTestsFromModule(test_pilot))
    stream=io.StringIO(); result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    (OUT/'software_tests.log').write_text(stream.getvalue(),encoding='utf8')
    if not result.wasSuccessful():
        print(stream.getvalue());raise SystemExit(1)
    search=json.loads((OUT/'phase25'/'search_audit.json').read_text())
    requested=sum(x['requested_starts'] for x in search.values())
    converged=sum(x['converged'] for x in search.values())
    if requested!=144:raise ValueError('Unexpected search size')
    with (OUT/'phase25'/'conditions.csv').open(encoding='utf8') as f:conditions=list(csv.DictReader(f))
    split={k:sum(x['split']==k for x in conditions) for k in sorted({x['split'] for x in conditions})}
    source=json.loads((OUT/'literature_structures'/'manifest.json').read_text())
    for record in source:
        file=OUT/'literature_structures'/(record['name']+'.xyz')
        if hashlib.sha256(file.read_bytes()).hexdigest()!=record['sha256']:raise ValueError('Source structure hash mismatch')
    pilot=json.loads((OUT/'phase25'/'local_qm_pilot'/'results.json').read_text())
    successes=[r for r in pilot if r['status']=='COMPLETED_SINGLE_POINT_ONLY']
    for r in successes:
        source=OUT/'phase25'/'starting_structures'/f"I01_{r['E_Z']}.xyz"
        if hashlib.sha256(source.read_bytes()).hexdigest()!=r['input_sha256']:raise ValueError('Pilot starting geometry mismatch')
        output=OUT/'phase25'/'local_qm_pilot'/r['raw_output']
        if hashlib.sha256(output.read_bytes()).hexdigest()!=r['output_sha256']:raise ValueError('Pilot output hash mismatch')
        if 'Iterations' not in output.read_text(encoding='utf8'):raise ValueError('Unexpected raw output')
    summaries=[]
    for method in sorted({r['method'] for r in successes}):
        rows={r['E_Z']:r for r in successes if r['method']==method}
        if set(rows)=={'E','Z'}:
            summaries.append(dict(method=method,E_hartree=rows['E']['electronic_energy_hartree'],Z_hartree=rows['Z']['electronic_energy_hartree'],
               Z_minus_E_kcal_mol=(rows['Z']['electronic_energy_hartree']-rows['E']['electronic_energy_hartree'])*627.5094740631,
               interpretation='electronic single-point comparison of isolated MMFF geometries; not Gibbs energies or activation barriers'))
    dump(OUT/'phase25'/'local_qm_pilot'/'summary.json',dict(comparisons=summaries,
       total_successful_core_seconds=sum(r['allocated_core_seconds'] for r in successes),
       complete_cost_accounting=False,cost_note='failed logged jobs retained; initial interrupted and budget-stopped attempts lack complete core-second totals; not an AL cost benchmark'))
    validation=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),tests_run=result.testsRun,failures=len(result.failures),errors=len(result.errors),
      test_scope='analytic software fixtures + structural/split checks, not scientific validation',requested_forcefield_starts=requested,converged_forcefield_starts=converged,
      conditions=len(conditions),split_counts=split,successful_isolated_single_points=len(successes),
      quantum_production_paths=0,solvated_interface_trajectories=0,nonadiabatic_trajectories=0,
      python=sys.version,numpy=numpy.__version__,scipy=scipy.__version__,rdkit=rdBase.rdkitVersion,platform=platform.platform(),
      scientific_acceptance='NOT_PASSED',execution_file_hashes_verified=execution_files_checked)
    progress=OUT/'phase25/pilot_progress.json'
    if progress.exists():
        p=json.loads(progress.read_text(encoding='utf8'))
        validation.update(xtb_audited_pilot_jobs=p['xTB_audited_jobs'],xtb_preliminary_minima=p['xTB_accepted_preliminary_minima'],
          new_DFT_optimizations=sum(r['status']=='OPTIMIZED_FREQUENCY_PENDING' for r in p['DFT_jobs']),
          native_Hessian_snapshot_jobs=json.loads((OUT/'phase25/hessian_checkpoint_snapshot/manifest.json').read_text(encoding='utf8'))['completed_native_gradient_jobs'])
    dump(OUT/'software_validation.json',validation)
    files=[]
    omission_file=OUT/'publication_omissions.json'
    omitted=json.loads(omission_file.read_text())['prefixes'] if omission_file.exists() else []
    for path in sorted(OUT.rglob('*')):
        if path.is_file() and path.name!='file_manifest.json' and not any(path.relative_to(ROOT).as_posix().startswith(prefix) for prefix in omitted):
            files.append(dict(path=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    for path in sorted((ROOT/'phase25_27').glob('*.py'))+sorted(ROOT.glob('PHASE2[567]_REPORT_*.md'))+sorted(ROOT.glob('PHASE25_PROGRESS_*.md'))+[ROOT/'phase25_27'/'README.md',ROOT/'requirements_phase25_27.txt']:
        files.append(dict(path=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    dump(OUT/'file_manifest.json',dict(hash_algorithm='sha256',files=files))
    print('PASS',result.testsRun,'software/structure tests;',len(successes),'isolated single points; scientific acceptance NOT PASSED')


if __name__=='__main__':main()
