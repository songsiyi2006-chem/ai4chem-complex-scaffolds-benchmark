"""Serial, checkpointed xTB optimization/Hessian campaign, never a DFT surrogate.

Specify a fresh tag and selected XYZ inputs; each raw run is isolated and retained.
Completed optimization requires exit code AND engine convergence marker.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time


def save(path,value):path.write_bytes((json.dumps(value,indent=2,allow_nan=False)+'\n').encode())


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--exe',required=True);ap.add_argument('--inputs',nargs='+',required=True)
    ap.add_argument('--tag',required=True);ap.add_argument('--solvent',default='toluene');ap.add_argument('--hessian',action='store_true')
    ap.add_argument('--threads',type=int,default=1);ap.add_argument('--opt',default='vtight')
    ap.add_argument('--single-point',action='store_true');args=ap.parse_args()
    if args.single_point and args.hessian:ap.error('Single-point and Hessian stages must be separate')
    root=Path(__file__).resolve().parents[1];out=root/'results_phase25_27/phase25/xtb_pilot'/args.tag
    out.mkdir(parents=True,exist_ok=False);exe=Path(args.exe).resolve()
    env=os.environ.copy();env['PATH']=str(exe.parent)+os.pathsep+env.get('PATH','')
    for key in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:env[key]=str(args.threads)
    env['OMP_STACKSIZE']='64M'
    version=subprocess.run([str(exe),'--version'],env=env,capture_output=True,text=True,timeout=30)
    (out/'engine_version.txt').write_bytes((version.stdout+version.stderr).encode())
    records=[];ledger=out/'results.json'
    for source in args.inputs:
        inp=(root/source).resolve();dest=out/inp.stem;dest.mkdir()
        data=inp.read_bytes();(dest/'input.xyz').write_bytes(data)
        row=dict(name=inp.stem,method='GFN2-xTB',solvation='ALPB('+args.solvent+')',
          charge=0,unpaired_electrons=0,threads=args.threads,optimization=args.opt,
          input_sha256=hashlib.sha256(data).hexdigest(),status='STARTED',
          started_utc=datetime.now(timezone.utc).isoformat(),production_catalysis_usable=False,
          stationary_point_validated=False)
        records.append(row);save(ledger,records);t=time.monotonic()
        base=['--gfn','2','--chrg','0','--uhf','0','--alpb',args.solvent,'--acc','0.2']
        try:
            command=[str(exe),'input.xyz',*base]+([] if args.single_point else ['--opt',args.opt])
            row['optimization_command']=['xtb',*command[1:]]
            with (dest/'opt.out').open('wb') as log:
                result=subprocess.run(command,cwd=dest,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=3600)
            raw=(dest/'opt.out').read_text(errors='replace');row['optimization_exit_code']=result.returncode
            converged=result.returncode==0 and 'GEOMETRY OPTIMIZATION CONVERGED' in raw and (dest/'xtbopt.xyz').exists()
            row['optimization_converged']=converged
            energies=re.findall(r'TOTAL ENERGY\s+([-+\d.]+)\s+Eh',raw)
            if energies:row['total_energy_hartree']=float(energies[-1])
            row['status']='OPTIMIZED_FREQUENCY_PENDING' if converged else 'OPTIMIZATION_FAILED'
            if args.single_point:
                row['status']='SINGLE_POINT_COMPLETED' if result.returncode==0 and energies else 'SINGLE_POINT_FAILED'
                row['calculation']='single_point_without_geometry_changes'
            if converged and args.hessian:
                row['status']='HESSIAN_RUNNING';save(ledger,records)
                cmd=[str(exe),'xtbopt.xyz',*base,'--hess','--strict']
                row['hessian_command']=['xtb',*cmd[1:]]
                with (dest/'hessian.out').open('wb') as log:
                    result=subprocess.run(cmd,cwd=dest,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=3600)
                raw=(dest/'hessian.out').read_text(errors='replace');row['hessian_exit_code']=result.returncode
                row['status']='HESSIAN_COMPUTED_REVIEW_REQUIRED' if result.returncode==0 and 'normal termination of xtb' in raw else 'HESSIAN_FAILED'
                # Independent parser in pilot_validation.py checks 3N count and
                # repeated blocks; retain the raw output as the primary record.
                row['frequency_validation']='Run python -m phase25_27.pilot_validation'
        except Exception as exc:row.update(status='FAILED',error=str(exc))
        row.update(wall_seconds=time.monotonic()-t,finished_utc=datetime.now(timezone.utc).isoformat())
        row['allocated_core_seconds']=row['wall_seconds']*args.threads
        row['files']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.iterdir() if p.is_file()}
        save(ledger,records);print(row['name'],row['status'],round(row['wall_seconds'],1),flush=True)


if __name__=='__main__':main()
