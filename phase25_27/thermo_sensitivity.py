"""Reuse actual xTB Hessians for exact engine 50/100 cm^-1 rotor comparisons.

No new Hessian or fitted inertia approximation is used. Check 50 cm^-1 replay
against the original thermochemistry before accepting the 100 cm^-1 contrast.
"""
import argparse
import hashlib
import json
import os
import subprocess
from .build_inputs import OUT
from .pilot_validation import main as audit,last_number,HARTREE_KCAL


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--exe',required=True);args=ap.parse_args()
    audit();rows=json.loads((OUT/'phase25/xtb_pilot/validation.json').read_text())['records']
    dest=OUT/'phase25/thermo_sensitivity';dest.mkdir(exist_ok=True)
    env=os.environ.copy();env['OMP_NUM_THREADS']='1';records=[]
    for row in rows:
        if not row['accepted_local_minimum']:continue
        src=OUT/'phase25/xtb_pilot'/row['campaign']/row['name'];out=dest/row['campaign']/row['name'];out.mkdir(parents=True,exist_ok=True)
        geometry=src/'xtbopt.xyz';hessian=src/'hessian'
        meta=dict(campaign=row['campaign'],name=row['name'],geometry_sha256=hashlib.sha256(geometry.read_bytes()).hexdigest(),
          hessian_sha256=hashlib.sha256(hessian.read_bytes()).hexdigest(),temperature_K=298.15,
          algorithm='native xtb thermo; engine geometry-dependent rotational inertia retained',values=[])
        for cutoff in [50,100]:
            logfile=out/f'cutoff{cutoff}.out';command=[args.exe,'thermo','--sthr',str(cutoff),'--temp','298.15',str(geometry),str(hessian)]
            if not logfile.exists():
                with logfile.open('wb') as f:result=subprocess.run(command,env=env,cwd=out,stdout=f,stderr=subprocess.STDOUT,timeout=120)
                if result.returncode:raise RuntimeError('Thermochemistry replay failed')
            raw=logfile.read_text(encoding='utf8',errors='strict')
            if 'normal termination of xtb' not in raw:raise ValueError('Incomplete thermo output')
            corr=last_number(raw,r'G\(RRHO\) contrib\.\s+([-+\d.]+)\s+Eh')
            if corr is None:raise ValueError('Missing thermochemical correction')
            meta['values'].append(dict(cutoff_cm1=cutoff,G_correction_hartree=corr,output_sha256=hashlib.sha256(logfile.read_bytes()).hexdigest()))
        nativecorr=row['raw_engine_G_hartree']-row['energy_including_ALPB_hartree']
        discrepancy=meta['values'][0]['G_correction_hartree']-nativecorr
        meta.update(replay50_difference_hartree=discrepancy,replay50_pass=abs(discrepancy)<1e-6,
          G_shift_50_to100_kcal_mol=(meta['values'][1]['G_correction_hartree']-meta['values'][0]['G_correction_hartree'])*HARTREE_KCAL)
        records.append(meta)
    (dest/'results.json').write_bytes((json.dumps(records,indent=2)+'\n').encode())
    print('Native rotor-cutoff comparisons:',len(records),'50 cm^-1 replays passed:',sum(r['replay50_pass'] for r in records))


if __name__=='__main__':main()
