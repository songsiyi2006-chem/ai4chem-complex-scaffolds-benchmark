"""Correct a duplicated printout field while preserving original execution ledgers.

Raw engine outputs, frequencies, energies, timings and acceptance are unchanged.
Only run after the finite xTB campaigns have finished writing their ledgers.
"""
import json
import numpy as np
from .build_inputs import OUT
from .pilot_validation import frequency_block,xyz


def main():
    changed=[]
    for file in sorted((OUT/'phase25/xtb_pilot').glob('*/results.json')):
        original=file.read_bytes();rows=json.loads(original);edits=[]
        if any(r['status'] in ['STARTED','HESSIAN_RUNNING'] for r in rows):raise RuntimeError('An xTB ledger is still live')
        for row in rows:
            if 'projected_frequencies_cm1' in row:
                folder=file.parent/row['name'];symbols,_=xyz(folder/'xtbopt.xyz');count=3*len(symbols)
                old=row['projected_frequencies_cm1']
                modes=frequency_block((folder/'hessian.out').read_text(encoding='utf8'),len(symbols))
                if len(old)!=count:
                    if len(old)%count:raise ValueError('Unexpected frequency length')
                    for start in range(0,len(old),count):np.testing.assert_allclose(old[start+6:start+count],modes,atol=.011,rtol=0)
                    row['projected_frequencies_cm1']=old[:count]
                    row['frequency_parse_correction']=dict(original_array_length=len(old),corrected_length=count,
                      reason='Repeated identical native printouts concatenated by the first runner; authoritative raw values unchanged')
                    edits.append(row['name'])
            if row.get('calculation')=='single_point_without_geometry_changes' and row['optimization'] is not None:
                row['unused_optimization_option']=row['optimization'];row['optimization']=None;edits.append(row['name']+' metadata')
        if edits:
            archive=file.with_name('original_execution_ledger.json')
            if archive.exists():raise FileExistsError('Existing original ledger must not be overwritten')
            archive.write_bytes(original)
            file.write_bytes((json.dumps(rows,indent=2,allow_nan=False)+'\n').encode())
            changed.append(dict(campaign=file.parent.name,corrected_records=edits,original_ledger=archive.name))
    report=OUT/'phase25/xtb_pilot/ledger_corrections.json'
    if changed:report.write_bytes((json.dumps(changed,indent=2)+'\n').encode())
    print('Corrected duplicate frequency/unused-option fields in',len(changed),'ledgers; raw evidence retained')


if __name__=='__main__':main()
