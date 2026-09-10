"""Real, serial Psi4 opt/freq jobs with explicit method and retained raw evidence.

Run in the Psi4 environment. Isolated gas-phase minima are not catalytic paths.
Each invocation creates a new directory; later freq jobs read an optimized XYZ.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
import traceback
import numpy as np
import psi4


def write_json(path, value):
    path.write_bytes((json.dumps(value, indent=2, allow_nan=False)+'\n').encode())


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',required=True)
    ap.add_argument('--tag',required=True)
    ap.add_argument('--method',default='pbe-d3bj')
    ap.add_argument('--stage',choices=['opt','freq','sp'],default='opt')
    ap.add_argument('--threads',type=int,default=2)
    ap.add_argument('--memory-mb',type=int,default=500)
    ap.add_argument('--maxiter',type=int,default=60)
    ap.add_argument('--resume',action='store_true')
    args=ap.parse_args()
    root=Path(__file__).resolve().parents[1]
    source=(root/args.input).resolve()
    out=root/'results_phase25_27'/'phase25'/'stationary_pilot'/args.tag
    if out.exists():
        if not args.resume:raise FileExistsError('Use a fresh tag or explicitly resume')
        old=json.loads((out/'result.json').read_text(encoding='utf8'))
        if (old['method'],old['stage'],old['input_sha256'])!=(args.method,args.stage,hashlib.sha256(source.read_bytes()).hexdigest()):
            raise ValueError('Resume input/method/stage mismatch')
        archive=out/'previous_attempts';archive.mkdir(exist_ok=True)
        (archive/(str(time.time_ns())+'.json')).write_bytes((out/'result.json').read_bytes())
    out.mkdir(parents=True,exist_ok=True)
    inp=source.read_bytes(); (out/'input.xyz').write_bytes(inp)
    molecule=psi4.geometry('0 1\n'+'\n'.join(inp.decode().splitlines()[2:])+'\nunits angstrom\nsymmetry c1\nno_reorient\nno_com')
    psi4.set_memory(f'{args.memory_mb} MB');psi4.set_num_threads(args.threads)
    scratch=root.parent/'qm_scratch'/args.tag;scratch.mkdir(parents=True,exist_ok=True)
    psi4.core.IOManager.shared_object().set_default_path(str(scratch))
    logfile=out/'psi4.out';psi4.core.set_output_file(str(logfile),False)
    psi4.set_options(dict(basis='def2-svp',scf_type='df',e_convergence=1e-9,d_convergence=1e-8,
      maxiter=150,dft_radial_points=75,dft_spherical_points=302,geom_maxiter=args.maxiter,
      g_convergence='gau_tight',t=298.15,p=101325.0))
    row=dict(tag=args.tag,stage=args.stage,method=args.method,basis='def2-SVP',charge=0,multiplicity=1,
      solvent='gas_phase',temperature_K=298.15,pressure_Pa=101325.0,grid=[75,302],
      psi4_version=psi4.__version__,input_sha256=hashlib.sha256(inp).hexdigest(),
      threads=args.threads,memory_MB=args.memory_mb,geometry_convergence='gau_tight',
      stationary_point_validated=False,production_catalysis_usable=False,status='STARTED',
      started_utc=datetime.now(timezone.utc).isoformat())
    ledger=out/'result.json';write_json(ledger,row);start=time.monotonic()
    try:
        if args.stage=='opt':
            energy,wfn=psi4.optimize(args.method,molecule=molecule,return_wfn=True)
            molecule.save_xyz_file(str(out/'optimized.xyz'),True)
            grad=wfn.gradient().np
            np.savetxt(out/'gradient_hartree_per_bohr.txt',grad)
            row.update(status='OPTIMIZED_FREQUENCY_PENDING',max_abs_gradient_hartree_per_bohr=float(np.max(np.abs(grad))),
              rms_gradient_hartree_per_bohr=float(np.sqrt(np.mean(grad**2))))
        elif args.stage=='freq':
            from .atomic_checkpoints import atomic_checkpoints
            with atomic_checkpoints(out/'atomic_jobs'):
                energy,wfn=psi4.frequency(args.method,molecule=molecule,return_wfn=True,dertype=1)
            omega=np.asarray(wfn.frequency_analysis['omega'].data,dtype=complex)
            # Preserve the imaginary component; never cast imaginary frequencies to real.
            row['frequencies_cm1']=[dict(real=float(x.real),imag=float(x.imag)) for x in omega]
            row['TRV_classification']=np.asarray(wfn.frequency_analysis['TRV'].data).tolist()
            np.savetxt(out/'hessian_hartree_per_bohr2.txt',wfn.hessian().np)
            row['thermochemistry_RRHO_hartree']={k:float(v) for k,v in psi4.core.variables().items()
              if any(s in k for s in ['GIBBS','ENTHALPY','ZERO K','ZPVE']) and isinstance(v,(float,int))}
            row.update(status='FREQUENCIES_COMPUTED_REVIEW_REQUIRED',frequency_derivatives='central finite difference of same-method gradients')
        else:
            energy,wfn=psi4.energy(args.method,molecule=molecule,return_wfn=True)
            row['status']='SINGLE_POINT_COMPLETED'
        row['electronic_energy_hartree']=float(energy)
    except Exception as exc:
        row.update(status='FAILED',error=str(exc),traceback=traceback.format_exc())
        molecule.save_xyz_file(str(out/('last_geometry_NOT_CONVERGED.xyz' if args.stage=='opt' else 'input_geometry_AFTER_FAILED_STAGE.xyz')),True)
    finally:
        row.update(wall_seconds=time.monotonic()-start,finished_utc=datetime.now(timezone.utc).isoformat())
        row['allocated_core_seconds']=row['wall_seconds']*args.threads
        psi4.core.close_outfile()
        row['files']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file() and p.name!='result.json'}
        write_json(ledger,row)
        print(args.tag,row['status'],round(row['wall_seconds'],2),flush=True)
    if row['status']=='FAILED':raise SystemExit(1)


if __name__=='__main__':main()
