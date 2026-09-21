"""Explicit PySCF adapters. Never silently replace the requested Hamiltonian.

Original omegaB97X-D is rejected by the current PySCF dispersion interface.
PBE0 here is an opt-in integration-test method, not fulfillment of that request.
NEVPT2 is a correction to a CAS reference, not to a DFT total energy.
"""
import argparse
import json
import math
from pathlib import Path
from common import EvidenceBlocked, write_json
from swarm.memory_ledger import digest

IMMUTABLE = ("atoms", "charge", "spin", "xc", "basis", "solvation", "calculation",
             "ncas", "nelecas", "active_orbitals")


def fingerprint(spec):
    return digest({key: spec.get(key) for key in IMMUTABLE})


def validate(spec):
    if spec.get("xc", "wb97x-d").lower().replace("_", "-") in {"wb97x-d", "ωb97x-d"}:
        raise EvidenceBlocked("Original omegaB97X-D unavailable in supported PySCF dispersion interface; no D3/D4 substitution")
    if spec.get("xc") != "pbe0" or not spec.get("explicit_alternative_authorized"):
        raise EvidenceBlocked("Only explicitly authorized PBE0 integration adapter is wired; production recipe remains blocked")
    if spec.get("solvation") != "vacuum":
        raise EvidenceBlocked("Solvation adapter not validated")
    if spec.get("basis") not in {"sto-3g", "def2-svp", "def2-tzvp"}:
        raise ValueError("unreviewed basis")
    if type(spec.get("charge")) is not int or type(spec.get("spin")) is not int or spec["spin"] < 0:
        raise ValueError("integer charge and spin=2S required")
    atoms = spec.get("atoms", [])
    if not 1 <= len(atoms) <= 200:
        raise ValueError("explicit atom list required, maximum 200")
    for symbol, xyz in atoms:
        if not isinstance(symbol, str) or len(xyz) != 3 or not all(math.isfinite(float(v)) for v in xyz):
            raise ValueError("invalid geometry in angstrom")
    if spec.get("calculation") not in {"single_point", "opt", "casscf_nevpt2"}:
        raise ValueError("unknown calculation")


def run(spec, numerical, out):
    validate(spec)
    try:
        from pyscf import gto, dft, scf, mcscf, mrpt, lib
    except ImportError as exc:
        raise EvidenceBlocked("PySCF not installed in this interpreter") from exc
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    lib.num_threads(1)
    mol = gto.M(atom=spec["atoms"], unit="Angstrom", basis=spec["basis"],
                charge=spec["charge"], spin=spec["spin"], max_memory=1024, verbose=4)
    if (mol.nelectron - mol.spin) % 2 or mol.spin > mol.nelectron:
        raise ValueError("electron number/spin parity mismatch")
    if spec["calculation"] == "casscf_nevpt2":
        reference = scf.RHF(mol) if mol.spin == 0 else scf.ROHF(mol)
    else:
        reference = dft.RKS(mol) if mol.spin == 0 else dft.UKS(mol)
        reference.xc = spec["xc"]
        reference.grids.level = 3
    reference.max_cycle = int(numerical.get("max_cycle", 120))
    reference.conv_tol = 1e-9
    reference.damp = float(numerical.get("damp", 0))
    reference.level_shift = float(numerical.get("level_shift", 0))
    if numerical.get("ediis", False):
        reference.diis = scf.EDIIS()
    reference.chkfile = str(out / "wavefunction.chk")
    energy = float(reference.kernel())
    result = dict(method_fingerprint=fingerprint(spec), converged=bool(reference.converged),
                  energy_hartree=energy, calculation=spec["calculation"],
                  evidence="EXECUTED_QM_NOT_BINDING_FREE_ENERGY", stability_verified=False)
    if not reference.converged:
        return result
    if spec["calculation"] == "opt":
        # Scanner convergence is required at every geometry step by PySCF wrapper.
        from pyscf.geomopt.geometric_solver import kernel
        geometry_converged, optimized = kernel(reference, maxsteps=80)
        reference = reference.reset(optimized)
        energy = float(reference.kernel())
        result.update(geometry_converged=bool(geometry_converged),
                      converged=bool(geometry_converged and reference.converged),
                      energy_hartree=energy, optimized_atoms=optimized.atom_coords(unit="Angstrom").tolist())
    elif spec["calculation"] == "casscf_nevpt2":
        ncas, nelecas = spec.get("ncas"), spec.get("nelecas")
        active = spec.get("active_orbitals")
        if type(ncas) is not int or not 2 <= ncas <= 12 or not isinstance(nelecas, list) or len(nelecas) != 2:
            raise ValueError("reviewed active space required, local ncas 2..12")
        if (any(type(n) is not int or n < 0 or n > ncas for n in nelecas)
                or sum(nelecas) > mol.nelectron or len(set(active or [])) != ncas):
            raise ValueError("inconsistent active space")
        cas = mcscf.CASSCF(reference, ncas, tuple(nelecas))
        mo = cas.sort_mo(active, base=0)
        cas.kernel(mo)
        result.update(converged=bool(cas.converged), cas_energy_hartree=float(cas.e_tot))
        if cas.converged:
            correction = float(mrpt.NEVPT(cas).kernel())
            result.update(nevpt2_correction_hartree=correction, energy_hartree=float(cas.e_tot + correction))
    if not math.isfinite(result["energy_hartree"]):
        raise ValueError("nonfinite electronic energy")
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--spec", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    request = json.loads(args.spec.read_text(encoding="utf-8"))
    write_json(args.out / "result.json", run(request["spec"], request["numerical"], args.out))
