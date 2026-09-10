"""Optional real local single-point pilot on full isolated ketimine E/Z structures.

This does NOT compute catalytic barriers, Gibbs energies, or selectivity.
Run using a Psi4 environment from repo root. Raw outputs and failed attempts retained.
"""
import argparse
import hashlib
import json
import time
from datetime import datetime,timezone
from pathlib import Path
import psi4


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--threads',type=int,default=4)
    ap.add_argument('--methods',nargs='+',default=['pbe-d3bj','wb97x-d'])
    ap.add_argument('--run-id',default=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
    args=ap.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=root/'results_phase25_27'/'phase25'/'local_qm_pilot';out.mkdir(parents=True,exist_ok=True)
    psi4.set_memory('3 GB');psi4.set_num_threads(args.threads)
    scratch=root.parent/'qm_scratch';scratch.mkdir(exist_ok=True)
    psi4.core.IOManager.shared_object().set_default_path(str(scratch))
    ledger=out/'results.json'; records=json.loads(ledger.read_text()) if ledger.exists() else []
    for method in args.methods:
        for ez in ['E','Z']:
            source=root/'results_phase25_27'/'phase25'/'starting_structures'/f'I01_{ez}.xyz'
            xyz=source.read_text(); atomlines='\n'.join(xyz.splitlines()[2:])
            tag=args.run_id+'_'+method+'_'+ez; logfile=out/(tag+'.out')
            if logfile.exists(): raise FileExistsError('Use a unique --run-id to retain raw evidence')
            psi4.core.clean();psi4.core.clean_options();psi4.core.set_output_file(str(logfile),False)
            molecule=psi4.geometry('0 1\n'+atomlines+'\nunits angstrom\nsymmetry c1\nno_reorient\nno_com')
            psi4.set_options({'basis':'def2-svp','scf_type':'df','e_convergence':1e-8,'d_convergence':1e-8,'maxiter':150,'dft_radial_points':75,'dft_spherical_points':302})
            row=dict(run_id=args.run_id,raw_output=logfile.name,method=method,basis='def2-SVP',solvent='gas_phase',E_Z=ez,charge=0,multiplicity=1,psi4_version=psi4.__version__,threads=args.threads,
               calculation='single_point_on_MMFF_start',input_sha256=hashlib.sha256(xyz.encode()).hexdigest(),stationary_point_validated=False,production_usable=False)
            row.update(status='STARTED',started_utc=datetime.now(timezone.utc).isoformat())
            records.append(row); ledger.write_text(json.dumps(records,indent=2,allow_nan=False),encoding='utf8')
            t=time.monotonic()
            try:
                energy=psi4.energy(method,molecule=molecule)
                row.update(status='COMPLETED_SINGLE_POINT_ONLY',electronic_energy_hartree=float(energy))
            except Exception as e: row.update(status='FAILED',error=str(e),electronic_energy_hartree=None)
            row['wall_seconds']=time.monotonic()-t;row['allocated_core_seconds']=row['wall_seconds']*args.threads
            if logfile.exists(): row['output_sha256']=hashlib.sha256(logfile.read_bytes()).hexdigest()
            ledger.write_text(json.dumps(records,indent=2,allow_nan=False),encoding='utf8')
            print(tag,row['status'],round(row['wall_seconds'],2),flush=True)


if __name__=='__main__': main()
