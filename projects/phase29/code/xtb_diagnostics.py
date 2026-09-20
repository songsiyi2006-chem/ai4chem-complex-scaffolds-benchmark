"""Bounded, optional molecular xTB preflight; never a regioselectivity model.

The nuclei and solvent model are fixed between the pyridinium cation (+1, UHF 0)
and its one-electron-reduced neutral radical (0, UHF 1). ALPB independently
equilibrates each state's continuum: these are fixed-nuclei electronic energy
differences, NOT rigorous nonequilibrium-solvent vertical electron affinities,
electrode potentials, free energies, or activation barriers.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

HARTREE_EV = 27.211386245988
PANEL = [
    ("Q01", "pyridinium", "[nH+]1ccccc1"),
    ("Q02", "3-methoxypyridinium", "COc1ccc[nH+]c1"),
    ("Q03", "3-cyanopyridinium", "N#Cc1ccc[nH+]c1"),
]
BOUNDARY = (
    "EXECUTED_MOLECULAR_PREFLIGHT: hypothetical N-protonated solution species, "
    "not verified catalytic intermediates. Single MMFF94 geometry per species; "
    "no conformer ensemble, surface, electrolyte, counterion, explicit solvent, "
    "constant potential, thermal correction, or kinetic calculation. ALPB MeCN "
    "is independently equilibrated for each state. Method sensitivity is not "
    "experimental uncertainty. No redox-potential or C2/C4 prediction."
)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def fixed_nuclei_difference(cation, radical):
    for key in ("geometry_sha256", "gfn", "solvent"):
        if cation[key] != radical[key]:
            raise ValueError("Incompatible state pair: " + key)
    if (cation["charge"], cation["uhf"], radical["charge"], radical["uhf"]) != (1, 0, 0, 1):
        raise ValueError("Expected cation singlet and neutral radical doublet")
    return (radical["energy_hartree"] - cation["energy_hartree"]) * HARTREE_EV


def self_test():
    a = dict(geometry_sha256="same", gfn=1, solvent="acetonitrile", charge=1, uhf=0, energy_hartree=-10.)
    b = dict(a, charge=0, uhf=1, energy_hartree=-10.1)
    assert abs(fixed_nuclei_difference(a, b) + 0.1 * HARTREE_EV) < 1e-10
    rejected = 0
    for key, value in (("gfn", 2), ("geometry_sha256", "other"), ("solvent", "water"), ("uhf", 0)):
        try:
            fixed_nuclei_difference(a, dict(b, **{key: value}))
        except ValueError:
            rejected += 1
    assert rejected == 4
    print("xTB comparison invariants: 5 checks passed")


def run(args):
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem

    executable = shutil.which(args.xtb)
    if not executable:
        raise FileNotFoundError("Supply an existing xTB executable with --xtb; nothing is installed automatically")
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise ValueError("Use a new empty output directory to retain prior evidence")
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), MKL_NUM_THREADS=str(args.threads))
    version = subprocess.run([executable, "--version"], capture_output=True, text=True, errors="replace", env=env, timeout=20)
    (out / "xtb_version.txt").write_text(version.stdout + version.stderr, encoding="utf-8")
    manifest = {
        "started_utc": datetime.now(timezone.utc).isoformat(), "evidence_type": "executed_calculation",
        "boundary": BOUNDARY, "script_sha256": digest(__file__), "rdkit_version": rdkit.__version__,
        "python_version": sys.version, "executable_sha256": digest(executable),
        "threads": args.threads, "per_job_timeout_seconds": args.timeout, "seed": 20260920,
        "solvent": "acetonitrile", "geometry_method": "RDKit ETKDGv3 then MMFF94; one cation conformer",
        "panel": [], "jobs": [],
    }
    pairs = []
    for molecule_id, name, smiles in PANEL:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        if Chem.GetFormalCharge(mol) != 1:
            raise ValueError("Panel cation charge mismatch")
        ring_n = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 7 and a.GetIsAromatic()]
        assert len(ring_n) == 1
        electrons_cation = sum(a.GetAtomicNum() for a in mol.GetAtoms()) - 1
        assert electrons_cation % 2 == 0
        params = AllChem.ETKDGv3()
        params.randomSeed = 20260920
        params.numThreads = args.threads
        if AllChem.EmbedMolecule(mol, params) != 0 or not AllChem.MMFFHasAllMoleculeParams(mol):
            raise RuntimeError("Embedding or MMFF parameter failure: " + molecule_id)
        status = AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
        if status != 0:
            raise RuntimeError("MMFF did not converge: " + molecule_id)
        inputs = out / "inputs"
        inputs.mkdir(exist_ok=True)
        xyz = inputs / (molecule_id + ".xyz")
        Chem.MolToXYZFile(mol, str(xyz))
        xyz_hash = digest(xyz)
        manifest["panel"].append({"id": molecule_id, "name": name, "cation_smiles": smiles,
                                  "geometry_sha256": xyz_hash, "atom_count": mol.GetNumAtoms(),
                                  "cation_electron_count": electrons_cation, "mmff_status": status})
        for gfn in (1, 2):
            states = []
            for state, charge, uhf in (("cation", 1, 0), ("radical", 0, 1)):
                assert (sum(a.GetAtomicNum() for a in mol.GetAtoms()) - charge - uhf) % 2 == 0
                job_id = f"{molecule_id}_gfn{gfn}_{state}"
                job_out = out / "raw" / job_id
                job_out.mkdir(parents=True)
                job = dict(id=job_id, molecule_id=molecule_id, state=state, gfn=gfn, charge=charge, uhf=uhf,
                           solvent="acetonitrile", geometry_sha256=xyz_hash, status="not_started")
                command = [executable, "input.xyz", "--gfn", str(gfn), "--chrg", str(charge),
                           "--uhf", str(uhf), "--alpb", "acetonitrile", "--sp", "--acc", "0.1", "--iterations", "300"]
                # Only native text outputs and exact input are retained, not restart binaries.
                with tempfile.TemporaryDirectory(prefix="phase29_xtb_") as scratch:
                    shutil.copy2(xyz, Path(scratch) / "input.xyz")
                    start = time.perf_counter()
                    try:
                        result = subprocess.run(command, cwd=scratch, capture_output=True, text=True,
                                                errors="replace", env=env, timeout=args.timeout)
                        stdout, stderr = result.stdout, result.stderr
                        job["returncode"] = result.returncode
                        energies = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh", stdout)
                        normal = "normal termination of xtb" in (stdout + stderr).lower()
                        job["scc_converged"] = "convergence criteria satisfied" in stdout
                        job["numerical_warnings"] = [line.strip() for line in (stdout + stderr).splitlines()
                                                     if "IEEE_" in line or "WARNING" in line.upper()]
                        if result.returncode == 0 and energies and normal and job["scc_converged"]:
                            job.update(status="success", energy_hartree=float(energies[-1]))
                        else:
                            job["status"] = "failed_or_unverified_termination"
                    except subprocess.TimeoutExpired as exc:
                        stdout = exc.stdout or b""
                        stderr = exc.stderr or b""
                        stdout = stdout.decode(errors="replace") if isinstance(stdout, bytes) else stdout
                        stderr = stderr.decode(errors="replace") if isinstance(stderr, bytes) else stderr
                        job["status"] = "timeout"
                    job["elapsed_seconds"] = time.perf_counter() - start
                    (job_out / "stdout.txt").write_text(stdout, encoding="utf-8")
                    (job_out / "stderr.txt").write_text(stderr, encoding="utf-8")
                    for filename in ("charges", "wbo", "input.xyz"):
                        source = Path(scratch) / filename
                        if source.exists():
                            shutil.copy2(source, job_out / filename)
                job["command"] = ["xtb"] + command[1:]
                job["stdout_sha256"] = digest(job_out / "stdout.txt")
                job["stderr_sha256"] = digest(job_out / "stderr.txt")
                write_json(job_out / "job.json", job)
                manifest["jobs"].append(job)
                states.append(job)
                print(job_id, job["status"], flush=True)
            if all(s["status"] == "success" for s in states):
                pairs.append(dict(molecule_id=molecule_id, name=name, gfn=gfn,
                                  delta_e_eV=fixed_nuclei_difference(*states)))
    with (out / "fixed_nuclei_differences.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["molecule_id", "name", "gfn", "delta_e_eV"])
        writer.writeheader()
        writer.writerows(pairs)
    manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["successful_jobs"] = sum(j["status"] == "success" for j in manifest["jobs"])
    manifest["total_jobs"] = len(manifest["jobs"])
    manifest["jobs_with_numerical_warnings"] = sum(bool(j.get("numerical_warnings")) for j in manifest["jobs"])
    manifest["method_sensitivity_eV"] = {}
    for molecule_id, _, _ in PANEL:
        values = {p["gfn"]: p["delta_e_eV"] for p in pairs if p["molecule_id"] == molecule_id}
        if len(values) == 2:
            manifest["method_sensitivity_eV"][molecule_id] = abs(values[2] - values[1])
    write_json(out / "manifest.json", manifest)
    lines = ["# Molecular xTB preflight", "", BOUNDARY, "",
             "ΔE = E(neutral radical) − E(pyridinium cation), at identical nuclei within each method.",
             "Negative values indicate a lower total electronic energy under this model; the free electron reference,",
             "electrode reference and nonequilibrium solvent response needed for electrochemical interpretation are absent.", "",
             "| Molecule | GFN1 ΔE / eV | GFN2 ΔE / eV | Absolute method difference / eV |",
             "|---|---:|---:|---:|"]
    for molecule_id, name, _ in PANEL:
        values = {p["gfn"]: p["delta_e_eV"] for p in pairs if p["molecule_id"] == molecule_id}
        if len(values) == 2:
            lines.append(f"| {name} | {values[1]:.4f} | {values[2]:.4f} | {abs(values[2]-values[1]):.4f} |")
    lines += ["", "Raw input geometries, engine stdout/stderr, charge/spin settings and hashes are retained.",
              f"SCC convergence and normal termination: {manifest['successful_jobs']}/{manifest['total_jobs']} jobs.",
              f"Numerical warnings: {manifest['jobs_with_numerical_warnings']}/{manifest['total_jobs']} jobs. "
              "This Windows build reports IEEE floating-point flags despite SCC convergence; the warnings are retained "
              "in each job manifest and stderr. These results are flagged diagnostics, not validated electrochemistry.",
              "This is an optional diagnostic independent of the reaction-panel ranking and synthetic models.", ""]
    (out / "README.md").write_text("\n".join(lines), encoding="utf-8")
    return int(manifest["successful_jobs"] != manifest["total_jobs"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xtb", default="xtb")
    parser.add_argument("--out", default="work/phase29_quantum_preflight")
    parser.add_argument("--threads", type=int, choices=(1, 2), default=2)
    parser.add_argument("--timeout", type=int, default=90)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
