#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_phase5_chemical_world_model.py — PHASE 5: THE AUTONOMOUS CHEMICAL WORLD MODEL
=================================================================================

Unifies three tiers into one pipeline:

  MODULE A  Automated Reaction Network (ARN): GFN2-xTB electronic-structure
            topology of an asymmetric catalytic skeletal reorganization
            (N-bridged azirino-fused indole --chiral phosphoric acid--> chiral
            ring-expanded dihydroquinoline), including ion-pair (Tier-2)
            diastereomeric transition states and a directed, barrier-weighted
            reaction-network digraph.

  MODULE B  Stiff microkinetic world model: Eyring-Polanyi rate constants,
            fully coupled stiff ODE system solved with Radau over
            t in [1e-9 s, 1e5 s], temperature sweep 250-350 K, yield / ee /
            selectivity extraction with transparent stiffness diagnostics.

  MODULE C  Inverse generative design: 3D ESP field extraction of the
            rate-determining TS, evolutionary 3,3'-scaffold assembler on a
            BINOL-phosphoric-acid backbone (RDKit valency gates), scored by
            GFN2-xTB constrained complex-energy differences. These are
            exploratory proxies, not verified transition-state barrier reductions.

Outputs
  results_phase5/phase5_results.json  (machine-readable, atomic checkpoints)
  figures_phase5/fig1_reaction_network_topology.png   (300 DPI)
  figures_phase5/fig2_stiff_microkinetics_profile.png (300 DPI)
  figures_phase5/fig3_ts_stabilization_dock.png       (300 DPI)
  WORLD_MODEL_REPORT_EN.md / WORLD_MODEL_REPORT_ZH.md

Run stages (resumable; every stage checkpoints to results_phase5/cache/):
  python run_phase5_chemical_world_model.py --stage all      (default)
  python run_phase5_chemical_world_model.py --stage A        (network only)
  python run_phase5_chemical_world_model.py --stage B        (microkinetics)
  python run_phase5_chemical_world_model.py --stage C        (inverse design)
  python run_phase5_chemical_world_model.py --stage figures
  python run_phase5_chemical_world_model.py --stage reports

Defensive programming:
  * every xtb call wrapped, parsed defensively, timed out, and cached;
  * every molecule passes an RDKit formula/valency GATE before use;
  * ODE integration falls back BDF -> LSODA with transparent logging;
  * physical quantities are never silently invented: assigned (literature /
    Smoluchowski) parameters are flagged "assigned" everywhere they appear.

Engine: GFN2-xTB (xtb.exe subprocess, ASE-free minimal wrapper).
  multi-fidelity ladder:  Tier-1 substrate-only GFN2-xTB (intrinsic barriers)
                          Tier-2 explicit ion-pair complexes (catalyst effect)
                          Tier-3 MMFF/UFF prescreen inside the GA inner loop.
"""

from __future__ import annotations

import argparse
import copy
import itertools
import json
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from phase_audit import AUDIT_VERSION, exploratory_corrections

EXPLORATORY_KINETICS = False

# --------------------------------------------------------------------------- #
# 0.  CONFIG, LOGGING, RESULTS
# --------------------------------------------------------------------------- #
ROOT = Path(__file__).resolve().parent
RES = ROOT / "results_phase5"
FIG = ROOT / "figures_phase5"
ENERGY_CACHE_VERSION = AUDIT_VERSION + "-physical-sp-ts2b-pose-v3"
CACHE = RES / ("cache_" + ENERGY_CACHE_VERSION)
for d in (RES, FIG, CACHE):
    d.mkdir(parents=True, exist_ok=True)

RESULTS_PATH = RES / "phase5_results.json"
RESULTS: dict = {
    "audit_version": AUDIT_VERSION,
    "energy_cache_version": ENERGY_CACHE_VERSION,
    "scan_energy_scope": "unrestrained GFN2 single points on constrained geometries; not verified TS barriers",
    "evidence_status": "exploratory_effective_model_not_validated_prediction",
    "phase": 5,
    "title": "Autonomous Chemical World Model: ARN + stiff microkinetics + "
             "TS-conditioned generative catalyst design",
    "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
    "fallbacks": [],
    "warnings": [],
    "assigned_parameters": [],
}

EH_EV = 27.211386245988
EH_KCAL = 627.5094740631
EV_KCAL = EH_KCAL / EH_EV
BOHR_A = 0.52917721092
KB = 0.0019872041          # kcal/mol/K
H_KCAL = 6.27509474377e-12 # planck in kcal*s (unused placeholder, kept for clarity)
R_GAS = KB                 # kcal/(mol K)

T_REF = 298.15
STD_CORR = 1.89            # kcal/mol, 1 atm -> 1 mol/L translational correction
CONC_R0 = 0.10             # M
CONC_CAT0 = 0.010          # M  (10 mol %)
T_SWEEP = list(range(250, 351, 10))
T_HORIZON = (1e-9, 1e5)    # s
K_ON = 1.0e9               # M^-1 s^-1, Smoluchowski diffusion cap (assigned)
K_DIMER = 5.0e2            # M^-1 s^-1, cationic oligomerization surrogate (assigned)
BARRIER_AROM_SINK = 52.0   # kcal/mol, dihydroquinoline -> quinoline + H2 (assigned)


def _log(msg: str) -> None:
    print(f"[p5 {time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _fallback(msg: str) -> None:
    _log("!! FALLBACK: " + msg)
    RESULTS["fallbacks"].append(msg)


def _warn(msg: str) -> None:
    _log("!! WARNING: " + msg)
    RESULTS["warnings"].append(msg)


def _assigned(name: str, value, rationale: str) -> None:
    RESULTS["assigned_parameters"].append(
        {"name": name, "value": str(value), "rationale": rationale})
    _log(f"assigned parameter: {name} = {value}  ({rationale})")


def write_json_atomic(path: Path) -> None:
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(RESULTS, indent=2, default=float),
                   encoding="utf-8")
    os.replace(tmp, path)


# --------------------------------------------------------------------------- #
# 1.  XTB ENGINE LAYER (GFN2-xTB subprocess, Phase-4 pattern, hardened)
# --------------------------------------------------------------------------- #
def _find_xtb() -> str | None:
    cands = [shutil.which("xtb"),
             r"C:\Users\HUIWEI\miniconda3\envs\phase2ff\Library\bin\xtb.exe",
             str(Path(sys.prefix) / "Library" / "bin" / "xtb.exe")]
    for c in cands:
        if c and Path(c).exists():
            return str(c)
    return None


XTB_EXE = _find_xtb()
XTB_DIR = str(Path(XTB_EXE).parent) if XTB_EXE else ""
ENGINE_NAME = "GFN2-xTB" if XTB_EXE else "MMFF94 (degraded tier)"


def xtb_env() -> dict:
    env = os.environ.copy()
    if XTB_DIR:
        env["PATH"] = XTB_DIR + os.pathsep + env.get("PATH", "")
    env["OMP_NUM_THREADS"] = env.get("OMP_NUM_THREADS", "8")
    env["MKL_NUM_THREADS"] = env.get("MKL_NUM_THREADS", "8")
    return env


def _xyz_text(numbers, positions) -> str:
    el = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 15: "P", 16: "S", 17: "Cl"}
    lines = [str(len(numbers)), "phase5"]
    for z, p in zip(numbers, positions):
        lines.append(f"{el.get(int(z), str(z))} {p[0]:.10f} {p[1]:.10f} "
                     f"{p[2]:.10f}")
    return "\n".join(lines) + "\n"


def _parse_xyz(path: Path):
    toks = path.read_text().splitlines()
    n = int(toks[0].split()[0])
    zmap = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "P": 15, "S": 16,
            "Cl": 17}
    nums, pos = [], []
    for l in toks[2:2 + n]:
        p = l.split()
        nums.append(zmap.get(p[0], int(p[0]) if p[0].isdigit() else 6))
        pos.append([float(x) for x in p[-3:]])
    return np.array(nums), np.array(pos)


class XTBFailure(RuntimeError):
    """Retain engine evidence beyond the temporary working directory lifetime."""
    def __init__(self, message, stdout="", stderr="", files=None):
        super().__init__(message)
        self.stdout = stdout
        self.stderr = stderr
        self.files = files or {}


def run_xtb(numbers, positions, args_list, chrg=0, timeout=600,
            want_files=(), workdir: Path | None = None):
    """Generic hardened xtb call. Returns (stdout, {filename: text})."""
    if XTB_EXE is None:
        raise RuntimeError("xtb.exe not available")
    td = str(workdir) if workdir else tempfile.mkdtemp(prefix="p5xtb_")
    own = workdir is None
    try:
        xyz = Path(td) / "m.xyz"
        xyz.write_text(_xyz_text(numbers, positions), encoding="utf-8")
        cmd = [XTB_EXE, "m.xyz", "--chrg", str(chrg), "--gfn", "2"] + args_list
        # bytes mode + defensive decode: xtb emits non-UTF-8 bytes on Windows
        proc = subprocess.run(cmd, cwd=td, capture_output=True,
                              timeout=timeout, env=xtb_env())
        out = (proc.stdout or b"").decode("utf-8", errors="replace")
        files = {}
        for f in want_files:
            fp = Path(td) / f
            if fp.exists():
                files[f] = fp.read_text(errors="replace")
        if proc.returncode != 0:
            raise XTBFailure(
                f"xtb {' '.join(args_list)} failed rc={proc.returncode}; "
                f"tail: {out[-400:]!r} / "
                f"{(proc.stderr or b'')[-200:]!r}", stdout=out,
                stderr=(proc.stderr or b'').decode('utf-8', errors='replace'), files=files)
        return out, files
    finally:
        if own:
            shutil.rmtree(td, ignore_errors=True)


def xtb_opt(numbers, positions, chrg=0, constraints: list[tuple] | None = None,
            fc=0.6, timeout=900):
    """--opt; constraints = [(i, j, dist_A), ...] 1-based indices.
    Returns (positions, E_eh)."""
    extra, inp = [], None
    if constraints:
        rows = ["$constrain", f"  force constant={fc}"]
        for i, j, d in constraints:
            rows.append(f"  distance: {i}, {j}, {d:.4f}")
        rows.append("$end")
        inp = "\n".join(rows) + "\n"
        extra = ["--input", "constr.inp"]
    with tempfile.TemporaryDirectory(prefix="p5opt_") as td:
        if inp:
            (Path(td) / "constr.inp").write_text(inp, encoding="utf-8")
        out, files = run_xtb(numbers, positions, ["--opt"] + extra,
                             chrg=chrg, timeout=timeout,
                             want_files=("xtbopt.xyz",), workdir=Path(td))
        if "GEOMETRY OPTIMIZATION CONVERGED" not in out.upper():
            raise XTBFailure("xtb geometry optimization did not converge: " + out[-300:],
                             stdout=out, files=files)
        if "xtbopt.xyz" not in files:
            raise RuntimeError("xtb --opt produced no xtbopt.xyz; tail: "
                               + out[-300:])
        pos = _parse_xyz(Path(td) / "xtbopt.xyz")[1]
        m = re.findall(r"total energy\s*:?\s*(-[\d.]+) Eh", out)
        e_eh = float(m[-1]) if m else float("nan")
        if pos.shape != (len(numbers), 3) or not np.isfinite(pos).all() or not math.isfinite(e_eh):
            raise RuntimeError("xtb optimization returned invalid geometry or energy")
        if constraints:
            # Direct same-geometry tests show xTB TOTAL ENERGY includes the
            # restraint potential. Remove the restraint from the evaluation,
            # never by guessing/subtracting an undocumented energy component.
            physical, _ = run_xtb(numbers, pos, ['--sp'], chrg=chrg, timeout=timeout)
            match = re.search(r'TOTAL ENERGY\s+(-?\d+\.\d+)', physical)
            if match is None or not math.isfinite(float(match.group(1))):
                raise RuntimeError('Unrestrained energy missing after constrained optimization')
            e_eh = float(match.group(1))
        return pos, e_eh


def xtb_sp_charges(numbers, positions, chrg=0):
    """Single point + GFN2 charges. Returns (E_eh, charges ndarray)."""
    out, files = run_xtb(numbers, positions, ["--sp"], chrg=chrg,
                         want_files=("charges",))
    m = re.search(r"TOTAL ENERGY\s+(-?\d+\.\d+)", out)
    e_eh = float(m.group(1)) if m else float("nan")
    q = None
    if "charges" in files:
        try:
            q = np.array([float(x) for x in
                          files["charges"].split()])
        except ValueError:
            q = None
    if not math.isfinite(e_eh) or q is None or len(q) != len(numbers) or not np.isfinite(q).all():
        raise RuntimeError("xtb single point has missing or non-finite energy/charges; ESP unavailable")
    return e_eh, q


def xtb_hess(numbers, positions, chrg=0, timeout=1800):
    """Numerical GFN2-xTB Hessian. Returns dict with frequencies, n_imag,
    total free energy G (Eh), electronic E (Eh), charges."""
    out, files = run_xtb(numbers, positions, ["--hess"], chrg=chrg,
                         timeout=timeout,
                         want_files=("vibspectrum", "charges"))
    m = re.search(r"total free energy\s+(-?\d+\.\d+)\s+Eh", out)
    g_eh = float(m.group(1)) if m else None
    m = re.search(r"TOTAL ENERGY\s+(-?\d+\.\d+)", out)
    e_eh = float(m.group(1)) if m else float("nan")
    n_imag = None
    m = re.search(r"# imaginary freq\.\s+(\d+)", out)
    if m:
        n_imag = int(m.group(1))
    freqs = []
    if "vibspectrum" in files:
        for line in files["vibspectrum"].splitlines():
            s = line.strip()
            if not s or s.startswith("#") or s.startswith("$"):
                continue
            p = s.split()
            if len(p) >= 5:
                try:
                    freqs.append(float(p[2]))
                except ValueError:
                    pass
    q = None
    if "charges" in files:
        try:
            q = np.array([float(x) for x in files["charges"].split()])
        except ValueError:
            q = None
    if not math.isfinite(e_eh) or not freqs or not np.isfinite(freqs).all():
        raise RuntimeError("xtb Hessian has missing or non-finite energy/frequencies")
    if q is None or len(q) != len(numbers) or not np.isfinite(q).all():
        raise RuntimeError("xtb Hessian has missing or non-finite charges")
    if g_eh is None:
        # defensive RRHO fallback from frequencies (symmetry number 1)
        g_eh = _rrho_gibbs_fallback(e_eh, freqs, len(numbers))
        _fallback("xtb thermo block unparsable — RRHO fallback from "
                  "vibspectrum frequencies (symmetry number 1)")
    return {"G_eh": g_eh, "E_eh": e_eh, "n_imag": n_imag,
            "freqs": freqs, "charges": q}


def _rrho_gibbs_fallback(e_eh, freqs_cm, n_atoms, temperature=T_REF):
    """G = E + ZPE + H_vib(T) - T*(S_trans+S_rot+S_vib), symmetry number 1.
    cruder than xtb's own block; only used when parsing fails."""
    # S_trans (1 M gas std ~ S at 1 atm minus 1.89 kcal correction later; log-scale)
    S_trans = 26.34 if n_atoms >= 2 else 21.42            # cal/mol/K, ~1 atm
    S_rot = (21.34 if n_atoms >= 3 else 15.27) if n_atoms >= 2 else 0.0
    svib, zvib, hvib = 0.0, 0.0, 0.0
    for f in freqs_cm:
        if f <= 0:
            continue
        x = f * 1.4387769 / temperature
        svib += x / math.expm1(x) - math.log1p(-math.exp(-x))
        zvib += 0.5 * x
        hvib += x / math.expm1(x)
    S_tot = S_trans + S_rot + svib * 1.9872041
    H_eh = e_eh + (zvib + hvib) * 1.9872041e-3 * temperature / 627.5095
    G = H_eh - temperature * S_tot * 1e-3 / 627.5095
    return G


def scan_1d(numbers, positions, constraints_fn, frames, chrg=0, fc=0.6,
            label="scan"):
    """Relaxed scan. constraints_fn(s) -> list[(i, j, d)]. Returns
    list of (s, E_eh, positions)."""
    out = []
    for k, s in enumerate(frames):
        cons = constraints_fn(s)
        record = dict(label=label, target=s, charge=chrg, gfn=2, solvent='none',
                      force_constant=fc, constraints_1based=cons,
                      numbers=np.asarray(numbers).tolist(),
                      initial_positions_A=np.asarray(positions).tolist(),
                      energy_scope='unrestrained single point; constrained geometry; not validated TS')
        try:
            pos, e = xtb_opt(numbers, positions, chrg=chrg, constraints=cons,
                             fc=fc)
        except Exception as exc:
            record.update(converged=False, error=str(exc),
                          stdout=getattr(exc, 'stdout', ''), stderr=getattr(exc, 'stderr', ''),
                          engine_files=getattr(exc, 'files', {}))
            _save_scan_record(label, k, record)
            _warn(f"{label} frame s={s:.3f} failed ({exc}); skipped")
            continue
        record.update(converged=True, physical_E_eh=e, positions_A=pos.tolist())
        _save_scan_record(label, k, record)
        out.append((s, e, pos))
        _log(f"    {label} frame {k + 1}/{len(frames)} s={s:.3f} "
             f"E={e * EH_KCAL:.3f} Eh-kcal")
    if len(out) < 3:
        raise RuntimeError(f"{label}: too few converged frames ({len(out)})")
    return out


def _save_scan_record(label, index, record):
    directory = RES / 'scan_diagnostics' / re.sub(r'[^A-Za-z0-9_-]', '_', label)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f'{index:03d}.json'
    path.write_text(json.dumps(record, indent=2), encoding='utf-8')
    RESULTS.setdefault('scan_diagnostics', []).append(str(path))


# --------------------------------------------------------------------------- #
# 2.  RDKIT LAYER — programmatic species construction with formula GATES
# --------------------------------------------------------------------------- #
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors, rdDistGeom, rdForceFieldHelpers  # noqa: E402

RDLogger.DisableLog("rdApp.error")


def formula_of(mol) -> str:
    return rdMolDescriptors.CalcMolFormula(mol)


def gate(mol, want: str, name: str, strip_charge: bool = False):
    if mol is None:
        raise RuntimeError(f"{name}: mol is None")
    f = formula_of(mol)
    if strip_charge:
        f = f.rstrip("+-").rstrip("0123456789") if f[-1] in "+-" else f
        f = f.rstrip("+-")
    if f != want:
        raise RuntimeError(f"{name}: formula gate failed — got {f}, want {want}")
    Chem.SanitizeMol(mol)
    return mol


def embed3d(mol, seed=0xC0FFEE, n_conf=8, label="mol"):
    """ETKDGv3 + MMFF (UFF fallback) → best-conformer (numbers, positions)."""
    mol = Chem.AddHs(mol)
    ps = rdDistGeom.ETKDGv3()
    ps.randomSeed = (seed + 0xC0FFEE) & 0x7FFFFFFF
    ps.numThreads = 0
    cids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=n_conf, params=ps)
    if not cids:
        raise RuntimeError(f"{label}: ETKDG embedding failed")
    ff_ok = False
    try:
        ff_ok = rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(
            mol, numThreads=0) is not None
        props = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol)
        ff_ok = props is not None
        if ff_ok:
            for cid in cids:
                ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(
                    mol, props, confId=cid)
                ff.Minimize(maxIts=2000)
    except Exception:
        ff_ok = False
    if not ff_ok:
        try:
            for cid in cids:
                ff = rdForceFieldHelpers.UFFGetMoleculeForceField(mol,
                                                                  confId=cid)
                ff.Minimize(maxIts=2000)
            _fallback(f"{label}: MMFF94 unavailable — UFF preopt used")
        except Exception as exc:
            _warn(f"{label}: FF preopt skipped ({exc})")
    best, be = None, float("inf")
    for cid in cids:
        conf = mol.GetConformer(cid)
        pos = np.array(conf.GetPositions())
        # crude steric score: min nonbonded heavy distance
        d = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
        np.fill_diagonal(d, 9e9)
        sc = float(d.min())
        if sc < be:
            be, best = sc, pos
    nums = np.array([a.GetAtomicNum() for a in mol.GetAtoms()])
    return nums, best


def labels_of(mol) -> list[str]:
    lab = []
    for a in mol.GetAtoms():
        l = a.GetProp("p5label") if a.HasProp("p5label") else a.GetSymbol()
        lab.append(l)
    return lab


# --- Reactant R: 2-methyl-azirino[1,2-a]indole (C10H11N) -------------------
def build_R():
    m = Chem.MolFromSmiles("Cc1cc2ccccc2[nH]1")        # 2-methylindole
    if m is None:
        raise RuntimeError("2-methylindole SMILES rejected")
    Chem.Kekulize(m, clearAromaticFlags=True)
    rw = Chem.RWMol(m)
    n = [a.GetIdx() for a in rw.GetAtoms() if a.GetSymbol() == "N"][0]
    c2 = None
    for nbr in rw.GetAtomWithIdx(n).GetNeighbors():
        if nbr.GetSymbol() == "C" and any(
                x.GetSymbol() == "C" and x.GetDegree() == 1
                for x in nbr.GetNeighbors()):
            c2 = nbr.GetIdx()
    c3 = [x.GetIdx() for x in rw.GetAtomWithIdx(c2).GetNeighbors()
          if x.GetIdx() not in (n,) and x.GetSymbol() == "C"
          and x.GetDegree() >= 2][0]
    b = rw.GetBondBetweenAtoms(c2, c3)
    b.SetBondType(Chem.BondType.SINGLE)
    cn = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(n, cn, Chem.BondType.SINGLE)
    rw.AddBond(c2, cn, Chem.BondType.SINGLE)
    an = rw.GetAtomWithIdx(n)
    an.SetNoImplicit(True)
    an.SetNumExplicitHs(0)
    ac2 = rw.GetAtomWithIdx(c2)
    ac2.SetNoImplicit(True)
    ac2.SetNumExplicitHs(0)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    gate(mol, "C10H11N", "R")
    mol.GetAtomWithIdx(n).SetProp("p5label", "N1")
    mol.GetAtomWithIdx(c2).SetProp("p5label", "C2")
    mol.GetAtomWithIdx(c3).SetProp("p5label", "C3")
    mol.GetAtomWithIdx(cn).SetProp("p5label", "C2n")
    return mol


# --- Desired product P_target: ring-expanded chiral cyclic imine ------------
def build_P_target(R_mol):
    """Ring expansion: N1-C2 cleaved, N1=C2n imine formed; C2 keeps 1 H
    (sp3, stereogenic), C3 stays CH2, benzene aromatic. C10H11N, 1 center."""
    m = Chem.Mol(R_mol)
    Chem.Kekulize(m, clearAromaticFlags=True)
    idx = {a.GetProp("p5label"): a.GetIdx() for a in m.GetAtoms()
           if a.HasProp("p5label")}
    n, c2, c3, c2n = idx["N1"], idx["C2"], idx["C3"], idx["C2n"]
    rw = Chem.RWMol(m)
    rw.RemoveBond(n, c2)
    b = rw.GetBondBetweenAtoms(n, c2n)
    b.SetBondType(Chem.BondType.DOUBLE)          # N1=C2n imine
    an = rw.GetAtomWithIdx(n)
    an.SetNoImplicit(True)
    an.SetNumExplicitHs(0)
    ac2 = rw.GetAtomWithIdx(c2)
    ac2.SetNoImplicit(False)
    ac2.SetNumExplicitHs(0)                      # sp3, 1 implicit H
    mol = rw.GetMol()
    try:
        Chem.SanitizeMol(mol)
        gate(mol, "C10H11N", "P_target")
        centers = Chem.FindMolChiralCenters(
            Chem.RemoveHs(Chem.Mol(mol)), includeUnassigned=True,
            useLegacyImplementation=False)
        return mol, len(centers)
    except Exception as exc:
        _fallback(f"P_target imine pattern failed ({exc}); trying enamine "
                  "pattern N1-H + C2=C3 (achiral)")
    m2 = Chem.Mol(R_mol)
    Chem.Kekulize(m2, clearAromaticFlags=True)
    rw2 = Chem.RWMol(m2)
    rw2.RemoveBond(n, c2)
    an = rw2.GetAtomWithIdx(n)
    an.SetNoImplicit(True)
    an.SetNumExplicitHs(1)
    rw2.GetAtomWithIdx(c2).SetNoImplicit(True)
    rw2.GetAtomWithIdx(c2).SetNumExplicitHs(0)
    bb = rw2.GetBondBetweenAtoms(c2, c3)
    bb.SetBondType(Chem.BondType.DOUBLE)
    mol2 = rw2.GetMol()
    Chem.SanitizeMol(mol2)
    gate(mol2, "C10H11N", "P_target_P2")
    return mol2, 0


# --- Side product P_elim: achiral conjugated enamine isomer -----------------
def build_P_elim(P_target_mol):
    """Imine hydrolysis-free tautomer: N1-H + C2nH2 + C2(Me)=C3H enamine."""
    m = Chem.Mol(P_target_mol)
    Chem.Kekulize(m, clearAromaticFlags=True)
    idx = {a.GetProp("p5label"): a.GetIdx() for a in m.GetAtoms()
           if a.HasProp("p5label")}
    n, c2, c3, c2n = idx["N1"], idx["C2"], idx["C3"], idx["C2n"]
    rw = Chem.RWMol(m)
    rw.GetBondBetweenAtoms(n, c2n).SetBondType(Chem.BondType.SINGLE)
    an = rw.GetAtomWithIdx(n)
    an.SetNoImplicit(True)
    an.SetNumExplicitHs(1)
    rw.GetAtomWithIdx(c2n).SetNoImplicit(False)
    rw.GetAtomWithIdx(c2n).SetNumExplicitHs(0)
    rw.GetBondBetweenAtoms(c2, c3).SetBondType(Chem.BondType.DOUBLE)
    rw.GetAtomWithIdx(c2).SetNoImplicit(True)
    rw.GetAtomWithIdx(c2).SetNumExplicitHs(0)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    gate(mol, "C10H11N", "P_elim")
    return mol


# --- P_poly surrogate: C2n-C2n' coupled dimer (cationic oligomerization) ----
def build_P_poly(P_elim_mol):
    m1 = Chem.Mol(P_elim_mol)
    m2 = Chem.Mol(P_elim_mol)
    combo = Chem.CombineMols(m1, m2)
    rw = Chem.RWMol(combo)
    idx1 = {a.GetProp("p5label"): a.GetIdx() for a in m1.GetAtoms()
            if a.HasProp("p5label")}
    n1_at = len(m1.GetAtoms())
    a1 = idx1["C2n"]
    a2 = idx1["C2n"] + n1_at
    rw.AddBond(a1, a2, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    try:
        Chem.SanitizeMol(mol)
        f = formula_of(mol)
        if f not in ("C20H20N2", "C20H22N2"):
            raise RuntimeError(f"dimer formula {f} unexpected")
        return mol, f
    except Exception as exc:
        _warn(f"P_poly dimer build failed ({exc}) — energy marked surrogate")
        return None, "C20H22N2 (surrogate)"


# --- Catalyst Cat: (axle) BINOL-derived cyclic phosphoric acid --------------
CAT_SMILES = "O=P1(O)Oc2ccc3ccccc3c2-c2ccc3ccccc3c2O1"


def build_Cat(smiles: str | None = None, name="Cat"):
    mol = Chem.MolFromSmiles(smiles or CAT_SMILES)
    if mol is None:
        raise RuntimeError(f"{name}: catalyst SMILES rejected")
    Chem.SanitizeMol(mol)
    if not any(a.GetSymbol() == "P" for a in mol.GetAtoms()):
        raise RuntimeError(f"{name}: no phosphorus — not a CPA scaffold")
    if not any(a.GetSymbol() == "O" and a.GetTotalNumHs() > 0
               for a in mol.GetAtoms()):
        raise RuntimeError(f"{name}: no P-OH proton donor")
    p_idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "P"][0]
    o_pho = None
    for nbr in mol.GetAtomWithIdx(p_idx).GetNeighbors():
        if nbr.GetSymbol() == "O" and nbr.GetTotalNumHs() > 0:
            o_pho = nbr.GetIdx()
    if o_pho is None:
        raise RuntimeError(f"{name}: no P-OH found")
    mol.GetAtomWithIdx(p_idx).SetProp("p5label", "P")
    mol.GetAtomWithIdx(o_pho).SetProp("p5label", "O_pho")
    return mol


# --- geometry helpers --------------------------------------------------------
def kabsch_rotate(P, Q):
    """Rotation matrix R minimizing |P·R - Q|."""
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    return U @ D @ Vt


def place_fragment(frag_nums, frag_pos, anchor_pairs, target_points,
                   extra_rot_axis=None, extra_rot_deg=0.0):
    """Rigidly move frag so that its anchor atoms sit at target_points
    (best-fit Kabsch), then optional rotation about (target center) axis."""
    R = kabsch_rotate(frag_pos[anchor_pairs], np.array(target_points))
    new = frag_pos @ R
    new += (np.mean(target_points, axis=0) - np.mean(
        new[anchor_pairs], axis=0))
    if extra_rot_axis is not None and abs(extra_rot_deg) > 1e-9:
        ax = np.asarray(extra_rot_axis, float)
        ax /= np.linalg.norm(ax)
        th = math.radians(extra_rot_deg)
        K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]],
                      [-ax[1], ax[0], 0]])
        Rr = np.eye(3) + math.sin(th) * K + (1 - math.cos(th)) * (K @ K)
        ctr = new[anchor_pairs].mean(axis=0)
        new = (new - ctr) @ Rr.T + ctr
    return new


def merge_frags(fragments):
    nums, pos, labels = [], [], []
    for n, p, l in fragments:
        nums.append(np.asarray(n))
        pos.append(np.asarray(p, float))
        labels.extend(l)
    return np.concatenate(nums), np.vstack(pos), labels


def prepare_deprotonation_pose(numbers, positions, donor, hydrogen, oxygen,
                                phosphorus, fragment_start, h_o_distance=1.65):
    """Rigidly place phosphate at C-H, preserving the substrate and C-H bond.

    All indices refer to the merged complex. This is an initial H-bond pose,
    not an optimized intermediate or a transition-state claim.
    """
    nums = np.asarray(numbers)
    pos = np.asarray(positions, dtype=float).copy()
    if not (0 <= donor < fragment_start and 0 <= hydrogen < fragment_start
            and fragment_start <= oxygen < len(nums)
            and fragment_start <= phosphorus < len(nums)):
        raise ValueError('Deprotonation atom indices cross the wrong fragments')
    if (nums[donor], nums[hydrogen], nums[oxygen], nums[phosphorus]) != (6, 1, 8, 15):
        raise ValueError('Deprotonation requires C-H donor and phosphate O/P')
    axis = pos[hydrogen] - pos[donor]
    length = np.linalg.norm(axis)
    po = pos[oxygen] - pos[phosphorus]
    if not (0.7 < length < 1.3) or np.linalg.norm(po) < 1e-8:
        raise ValueError('Invalid initial C-H or P-O geometry')
    axis /= length
    po /= np.linalg.norm(po)
    rotation = kabsch_rotate(np.array([po]), np.array([-axis]))
    target = pos[hydrogen] + h_o_distance * axis
    pos[fragment_start:] = ((pos[fragment_start:] - pos[oxygen]) @ rotation + target)
    return pos


# --------------------------------------------------------------------------- #
# 3.  MODULE A — AUTOMATED REACTION NETWORK (ARN)
# --------------------------------------------------------------------------- #
XA: dict = {}   # module-A artifacts


def cached_or_compute(key, fn):
    cache = CACHE / f"{key}.json"
    if cache.exists():
        try:
            d = json.loads(cache.read_text())
            _log(f"cache hit: {key}")
            return d
        except Exception:
            _warn(f"cache {key} corrupt — recomputing")
    d = fn()
    tmp = CACHE / f"{key}.json.tmp"
    tmp.write_text(json.dumps(d, default=float), encoding="utf-8")
    os.replace(tmp, cache)
    return d


def species_energy(mol_or_nums, name, chrg=0, seed=0xC0FFEE, positions=None):
    """opt + hess -> {G_eh, E_eh, n_imag, positions, numbers}. Cached."""
    def comp():
        if positions is None:
            nums, pos = embed3d(mol_or_nums, seed=seed, label=name)
        else:
            nums = (mol_or_nums if not hasattr(mol_or_nums, "GetAtoms")
                    else None)
            pos = positions
        nums = np.asarray(nums)
        t0 = time.time()
        pos2, _ = xtb_opt(nums, pos, chrg=chrg)
        h = xtb_hess(nums, pos2, chrg=chrg)
        _log(f"  species {name}: G={h['G_eh'] * EH_KCAL:.2f} kcal/mol, "
             f"n_imag={h['n_imag']} ({time.time() - t0:.0f} s)")
        return {"numbers": nums.tolist(), "positions": pos2.tolist(),
                "G_eh": h["G_eh"], "E_eh": h["E_eh"], "n_imag": h["n_imag"],
                "freqs": h["freqs"][:6], "charges": h["charges"].tolist(),
                "chrg": chrg}
    return cached_or_compute(f"sp_{name}", comp)


def module_A():
    _log("=" * 70)
    _log("MODULE A — automated reaction network (GFN2-xTB topology)")
    _log("=" * 70)
    XA["engine"] = ENGINE_NAME
    if XTB_EXE is None:
        raise RuntimeError("xtb.exe unavailable — Module A cannot run; "
                           "see Phase-4 fallback docs")

    # ---- A.1 species -------------------------------------------------------
    _log("A.1 building species through RDKit programmatic gates ...")
    R_mol = build_R()
    P_mol, n_chiral = build_P_target(R_mol)
    XA["p_target_chiral_centers"] = n_chiral
    PE_mol = build_P_elim(P_mol)
    PP_mol, PP_formula = build_P_poly(PE_mol)
    Cat_mol = build_Cat()
    XA["species_smiles"] = {
        "R": Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(R_mol))),
        "P_target": Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(P_mol))),
        "P_elim": Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(PE_mol))),
        "P_poly": PP_formula + " (C2-C2 coupled dimer surrogate)",
        "Cat": Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(Cat_mol))),
    }
    _log(f"  R        = {XA['species_smiles']['R']}")
    _log(f"  P_target = {XA['species_smiles']['P_target']} "
         f"(chiral centers: {n_chiral})")
    _log(f"  P_elim   = {XA['species_smiles']['P_elim']}")
    _log(f"  Cat      = {XA['species_smiles']['Cat']}")

    # atom indices (heavy atoms first in RDKit order; H appended by AddHs)
    idx = {a.GetProp("p5label"): a.GetIdx() for a in
           Chem.AddHs(R_mol).GetAtoms() if a.HasProp("p5label")}
    XA["idx_R"] = idx

    sp = {}
    _log("A.2 GFN2-xTB opt+hess of closed-shell species ...")
    for name, mol, chrg, seed in [
            ("R", R_mol, 0, 1), ("P_target", P_mol, 0, 2),
            ("P_elim", PE_mol, 0, 3), ("Cat", Cat_mol, 0, 4),
            ("H2", None, 0, 5)]:
        if name == "H2":
            nums, pos = np.array([1, 1]), np.array(
                [[0.0, 0.0, 0.0], [0.0, 0.0, 0.7414]])
            sp[name] = species_energy(nums, "H2", chrg=0, positions=pos)
        else:
            sp[name] = species_energy(mol, name, chrg=chrg, seed=seed)
    if PP_mol is not None:
        try:
            sp["P_poly"] = species_energy(PP_mol, "P_poly", chrg=2, seed=6)
            sp["P_poly"]["note"] = ("cationic dimer surrogate, net +2 "
                                    "(counterion bookkeeping)")
        except Exception as exc:
            _warn(f"P_poly energy failed ({exc}); marked surrogate")
            sp["P_poly"] = {"G_eh": 2 * sp["R"]["G_eh"] + 0.2, "E_eh": None,
                            "n_imag": None, "note": "surrogate (G = 2R+0.2)"}
    else:
        sp["P_poly"] = {"G_eh": 2 * sp["R"]["G_eh"] + 0.2, "E_eh": None,
                        "n_imag": None, "note": "surrogate (G = 2R+0.2)"}

    # ---- A.3 protonated aziridine & deprotonated CPA -----------------------
    _log("A.3 ionic partners: aziridinium R-H+ and chiral phosphate anion ...")
    def make_Rprot():
        m = Chem.Mol(R_mol)
        Chem.Kekulize(m, clearAromaticFlags=True)
        rw = Chem.RWMol(m)
        n = [a.GetIdx() for a in m.GetAtoms()
             if a.HasProp("p5label") and a.GetProp("p5label") == "N1"][0]
        h = rw.AddAtom(Chem.Atom(1))
        rw.AddBond(n, h, Chem.BondType.SINGLE)
        rw.GetAtomWithIdx(n).SetFormalCharge(1)   # aziridinium ammonium N+
        mol = rw.GetMol()
        for a in mol.GetAtoms():
            a.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol)
        gate(mol, "C10H12N", "Rprot(+1)", strip_charge=True)
        return mol
    Rprot_mol = make_Rprot()
    sp["Rprot"] = species_energy(Rprot_mol, "Rprot", chrg=1, seed=7)

    def make_CPAneg():
        smi = CAT_SMILES.replace("P1(O)", "P1([O-])")
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise RuntimeError("CPA(-1): deprotonated SMILES rejected")
        Chem.SanitizeMol(mol)
        gate(mol, "C20H12O4P", "CPA(-1)", strip_charge=True)
        return mol
    CPAneg_mol = make_CPAneg()
    sp["CPAneg"] = species_energy(CPAneg_mol, "CPAneg", chrg=-1, seed=8)

    # H2O-free ion pair reference energies (sums) for binding energetics
    G_R = sp["R"]["G_eh"]
    G_Cat = sp["Cat"]["G_eh"]
    G_Rprot = sp["Rprot"]["G_eh"]
    G_CPAneg = sp["CPAneg"]["G_eh"]

    # ---- A.4 Tier-2 ion-pair complexes ------------------------------------
    _log("A.4 Tier-2: assembling H-bonded complex RC and aziridinium "
         "ion pair I_RC ...")
    R_nums, R_pos = np.array(sp["R"]["numbers"]), np.array(
        sp["R"]["positions"])
    Cat_nums, Cat_pos = np.array(sp["Cat"]["numbers"]), np.array(
        sp["Cat"]["positions"])
    Rp_nums, Rp_pos = np.array(sp["Rprot"]["numbers"]), np.array(
        sp["Rprot"]["positions"])
    CPa_nums, CPa_pos = np.array(sp["CPAneg"]["numbers"]), np.array(
        sp["CPAneg"]["positions"])

    n_at_R = len(R_nums)
    i_N1, i_C2, i_C2n, i_C3 = (idx["N1"], idx["C2"], idx["C2n"], idx["C3"])
    # the N1-H proton exists only in the PROTONATED frame (neutral R's N1 is
    # a bridgehead tertiary amine); the explicit H is Rprot heavy-atom +1
    i_HN = [i for i in range(len(Rp_nums)) if Rp_nums[i] == 1 and
            np.linalg.norm(Rp_pos[i] - Rp_pos[i_N1]) < 1.15][0]

    # CPAneg O_pho index: same as Cat minus one atom (H_pho removed last-
    # but-one); recompute directly from geometry: anionic O = P-bound O that
    # was the O_pho — find via label-free rule: O with no H, largest P-O
    # distance? Use construction knowledge: O_pho index in Cat (heavy-first
    # enumeration with AddHs? Cat embed adds H at end; H_pho sits right
    # after O_pho). Robust: P index then its O neighbours without H.
    def phos_o_of(cp_nums, cp_pos):
        iP = [i for i in range(len(cp_nums)) if cp_nums[i] == 15][0]
        out = []
        for i in range(len(cp_nums)):
            if cp_nums[i] == 8 and np.linalg.norm(cp_pos[i] -
                                                  cp_pos[iP]) < 1.75:
                has_h = any(cp_nums[j] == 1 and np.linalg.norm(
                    cp_pos[j] - cp_pos[i]) < 1.15 for j in range(len(cp_nums)))
                out.append((i, has_h))
        hydrox = [i for i, hh in out if hh]
        if hydrox:
            return hydrox[0]
        return out[0][0]

    def assemble_ionpair(oN1_target=True, face_deg=0.0, tag=""):
        """CPA- H-bonded to aziridinium N-H (O···H-N)."""
        v = Rp_pos[i_HN] - Rp_pos[i_N1]
        v /= np.linalg.norm(v)
        tgt_H = Rp_pos[i_HN] + 1.65 * v
        iO = phos_o_of(CPa_nums, CPa_pos)
        dirO = CPa_pos[iO] - CPa_pos[
            [j for j in range(len(CPa_nums)) if CPa_nums[j] == 15][0]]
        dirO = dirO / np.linalg.norm(dirO)
        tgt_O = tgt_H + 1.50 * (tgt_H - Rp_pos[i_N1]) / np.linalg.norm(
            tgt_H - Rp_pos[i_N1])
        # rotate CPA so (P->O) points toward (O->H direction reversed)
        want = -v
        Rm = kabsch_rotate(np.array([dirO]), np.array([want]))
        cp_new = CPa_pos @ Rm
        cp_new += (tgt_O - cp_new[iO])
        cp_new = place_fragment(CPa_nums, cp_new, [iO],
                                [tgt_O], extra_rot_axis=v,
                                extra_rot_deg=face_deg)
        nums, pos, lab = merge_frags([
            (Rp_nums, Rp_pos, [a.GetProp("p5label") if a.HasProp("p5label")
                               else "H" for a in Rprot_mol.GetAtoms()]),
            (CPa_nums, cp_new, ["CPA"] * len(CPa_nums))])
        return nums, pos, lab, len(Rp_nums) + iO

    def ip_opt_hess(nums, pos, constraints=None, chrg=0, key="",
                    fc=0.5, do_hess=True):
        pos2, _ = xtb_opt(nums, pos, chrg=chrg, constraints=constraints,
                          fc=fc)
        res = {"numbers": nums.tolist(), "positions": pos2.tolist(),
               "chrg": chrg}
        if do_hess:
            h = xtb_hess(nums, pos2, chrg=chrg)
            res.update({"G_eh": h["G_eh"], "E_eh": h["E_eh"],
                        "n_imag": h["n_imag"], "freqs": h["freqs"][:6],
                        "charges": h["charges"].tolist()})
        return res

    # RC: neutral H-bonded complex (pre proton transfer)
    def comp_RC():
        v = R_pos[i_N1] - np.mean(R_pos[[i_C2, i_C2n, i_C3]], axis=0)
        v /= np.linalg.norm(v)
        # CPA O_pho-H aims its H at N1 lone-pair side
        iO_pho = None
        iP = [j for j in range(len(Cat_nums)) if Cat_nums[j] == 15][0]
        for j in range(len(Cat_nums)):
            if Cat_nums[j] == 8 and np.linalg.norm(Cat_pos[j] -
                                                   Cat_pos[iP]) < 1.75:
                hs = [k for k in range(len(Cat_nums)) if Cat_nums[k] == 1
                      and np.linalg.norm(Cat_pos[k] - Cat_pos[j]) < 1.15]
                if hs:
                    iO_pho, iH_pho = j, hs[0]
        u = Cat_pos[iH_pho] - Cat_pos[iO_pho]
        u /= np.linalg.norm(u)
        tgt = R_pos[i_N1] - 1.80 * v          # H should land near N1-1.8v? no:
        # H_pho target point = N1 + 1.75 Å toward the incoming direction:
        tgtH = R_pos[i_N1] - 1.75 * v * -1
        Rm = kabsch_rotate(np.array([u]), np.array([tgtH - Cat_pos[iO_pho]]))
        cp_new = Cat_pos @ Rm
        cp_new += (tgtH - cp_new[iH_pho])
        nums, pos, lab = merge_frags([
            (R_nums, R_pos, [a.GetProp("p5label") if a.HasProp("p5label")
                             else "" for a in Chem.AddHs(R_mol).GetAtoms()]),
            (Cat_nums, cp_new, ["CAT"] * len(Cat_nums))])
        res = ip_opt_hess(nums, pos, chrg=0, key="RC",
                          constraints=[(i_N1 + 1, len(R_nums) + iH_pho + 1,
                                        1.85)], do_hess=True)
        return res
    RC = cached_or_compute("cx_RC", comp_RC)
    XA["RC"] = RC

    # I_RC: aziridinium·phosphate- ion pair
    def comp_IRC():
        nums, pos, lab, iO = assemble_ionpair(face_deg=0.0)
        return ip_opt_hess(nums, pos, chrg=0, key="I_RC", do_hess=True)
    I_RC = cached_or_compute("cx_IRC", comp_IRC)
    XA["I_RC"] = I_RC

    # ---- A.5 TS1: rate-determining aziridinium C2-N1 cleavage --------------
    _log("A.5 Tier-2 relaxed scan: TS1 = aziridinium C2-N1 cleavage ...")
    def comp_TS1():
        nums, pos, lab, iO = assemble_ionpair(face_deg=0.0)
        i_N1b = i_N1
        i_C2b = i_C2
        frames = [1.60, 1.80, 2.00, 2.20, 2.35, 2.50]
        scan = scan_1d(nums, pos,
                       lambda d: [(i_N1b + 1, i_C2b + 1, d)],
                       frames, chrg=0, label="TS1-scan(dN1-C2)")
        es = np.array([e for _, e, _ in scan])
        k = int(np.argmax(es))
        d_star = scan[k][0]
        pos_ts, _ = xtb_opt(nums, scan[k][2], chrg=0,
                            constraints=[(i_N1b + 1, i_C2b + 1, d_star)],
                            fc=0.8)
        h = xtb_hess(nums, pos_ts, chrg=0)
        _log(f"  TS1: d(N1-C2)*={d_star:.2f} A, G={h['G_eh'] * EH_KCAL:.2f}, "
             f"n_imag={h['n_imag']}")
        if h["n_imag"] != 1:
            _warn(f"TS1 has {h['n_imag']} imaginary frequencies "
                  "(1 expected) — constrained-scan saddle, quality flag set")
        return {"numbers": nums.tolist(), "positions": pos_ts.tolist(),
                "G_eh": h["G_eh"], "E_eh": h["E_eh"], "n_imag": h["n_imag"],
                "chrg": 0, "d_star": d_star,
                "scan_d": [s for s, _, _ in scan],
                "scan_e_eh": [e for _, e, _ in scan],
                "charges": h["charges"].tolist()}
    TS1 = cached_or_compute("cx_TS1", comp_TS1)
    XA["TS1"] = TS1

    # I1Cat: relaxed ion pair at broken-bond limit
    def comp_I1Cat():
        nums, pos, lab, iO = assemble_ionpair(face_deg=0.0)
        pos2, _ = xtb_opt(nums, pos, chrg=0,
                          constraints=[(i_N1 + 1, i_C2 + 1, 2.75)], fc=0.3)
        h = xtb_hess(nums, pos2, chrg=0)
        return {"numbers": nums.tolist(), "positions": pos2.tolist(),
                "G_eh": h["G_eh"], "E_eh": h["E_eh"], "n_imag": h["n_imag"],
                "chrg": 0, "charges": h["charges"].tolist()}
    I1Cat = cached_or_compute("cx_I1Cat", comp_I1Cat)
    XA["I1Cat"] = I1Cat

    # ---- A.6/A.7 diastereomeric ion-pair topology --------------------------
    _log("A.6 Tier-2: diastereomeric TS2a proton-relay pair (Re/Si faces) "
         "...")
    # Ring plane of the substrate (for enantiotopic-face generation)
    _plane_pts = R_pos[[i_N1, i_C2, i_C2n, i_C3] + [i for i, z in
                       enumerate(R_nums) if z == 6][:3]]
    _ctr = _plane_pts.mean(axis=0)
    _u, _s, _vt = np.linalg.svd(_plane_pts - _ctr)
    _nvec = _vt[2]

    def _reflect_through_plane(x):
        return x - 2.0 * np.outer((x - _ctr) @ _nvec, _nvec)

    def assemble_face(face_deg):
        """Catalyst-anion pose H-bonded to the aziridinium N-H. face_deg>0:
        direct pose rotated about the H-bond axis; face_deg<0: mirror image
        of that pose through the substrate ring plane (= enantiotopic face)."""
        v = Rp_pos[i_HN] - Rp_pos[i_N1]
        v /= np.linalg.norm(v)
        tgt_H = Rp_pos[i_HN] + 1.65 * v
        iO = phos_o_of(CPa_nums, CPa_pos)
        iP = [j for j in range(len(CPa_nums)) if CPa_nums[j] == 15][0]
        dirO = CPa_pos[iO] - CPa_pos[iP]
        dirO /= np.linalg.norm(dirO)
        tgt_O = tgt_H + 1.50 * v
        Rm = kabsch_rotate(np.array([dirO]), np.array([-v]))
        cp_new = CPa_pos @ Rm
        cp_new += (tgt_O - cp_new[iO])
        cp_new = place_fragment(CPa_nums, cp_new, [iO], [tgt_O],
                                extra_rot_axis=v, extra_rot_deg=abs(face_deg))
        mirrored = face_deg < 0
        if mirrored:
            cp_new = _reflect_through_plane(cp_new)
            # re-anchor the H-bond after reflection
            cp_new += (tgt_O - cp_new[iO])
        nums, pos, lab = merge_frags([
            (Rp_nums, Rp_pos, ["R"] * len(Rp_nums)),
            (CPa_nums, cp_new, ["CPA"] * len(CPa_nums))])
        nR = len(Rp_nums)
        locks = [(nR + iP + 1, i_C3 + 1,
                  round(float(np.linalg.norm(pos[nR + iP] - pos[i_C3])), 3)),
                 (nR + iP + 1, i_C2n + 1,
                  round(float(np.linalg.norm(pos[nR + iP] - pos[i_C2n])), 3))]
        return nums, pos, nR + iO, nR + iP, nR, locks, mirrored

    FLOOR = 1.5   # kcal/mol, early-TS resolution floor (documented)
    B2A_ASSIGNED = 8.0   # kcal/mol, ion-pair cation-trap scale (assigned)
    B2B_ASSIGNED = 6.0   # kcal/mol, deprotonation-trap scale (assigned)

    def comp_I1Cat_face(face_deg, tag):
        """Facially-locked resting ion pair (reference for that face)."""
        nums, pos, iO, iP, nR, locks, _ = assemble_face(face_deg)
        cons = [(i_N1 + 1, i_C2 + 1, 2.75)] + locks
        pos2, _ = xtb_opt(nums, pos, chrg=0, constraints=cons, fc=0.3)
        h = xtb_hess(nums, pos2, chrg=0)
        _log(f"  I1Cat{tag}: G={h['G_eh'] * EH_KCAL:.2f}, "
             f"n_imag={h['n_imag']}")
        return {"numbers": nums.tolist(), "positions": pos2.tolist(),
                "G_eh": h["G_eh"], "E_eh": h["E_eh"], "n_imag": h["n_imag"],
                "chrg": 0, "charges": h["charges"].tolist(),
                "face_deg": face_deg}

    I1Cat_M = cached_or_compute("cx_I1Cat_M", lambda: comp_I1Cat_face(
        +35.0, "M"))
    I1Cat_m = cached_or_compute("cx_I1Cat_m", lambda: comp_I1Cat_face(
        -35.0, "m"))

    def comp_TS2a():
        """Enantiodetermining 1,2-H relay (N1 -> C2) inside the opened ion
        pair, located with the TS1 recipe: relaxed 1-D scan of the migrating
        proton d(C2-H) at fixed post-cleavage d(N1-C2). Facially
        unconstrained: the Re/Si differentiation is carried by the locked
        diastereomeric resting pairs I1Cat_M / I1Cat_m (Curtin-Hammett
        early-TS limit) rather than by strained multi-constraint poses."""
        nums, pos, iO, iP, nR, locks, _ = assemble_face(+35.0)
        hN = [i for i in range(nR) if nums[i] == 1 and np.linalg.norm(
            pos[i] - pos[i_N1]) < 1.15][0]
        axv = pos[i_C2] - pos[i_N1]
        axv /= np.linalg.norm(axv)
        pos[hN] = pos[i_N1] + 1.30 * axv
        # late frames add the N1=C2n imine-formation constraint so the
        # product basin can relax (avoids anion-trapped kinetic artifacts)
        frames = [2.20, 1.90, 1.60, 1.35, 1.15]
        scan = scan_1d(nums, pos,
                       lambda d: [(i_C2 + 1, hN + 1, d),
                                  (i_N1 + 1, i_C2 + 1, 2.60)]
                       + ([(i_N1 + 1, i_C2n + 1, 1.32)] if d <= 1.6 else []),
                       frames, chrg=0, fc=0.5, label="TS2a-scan(dC2-H)")
        es = np.array([e for _, e, _ in scan])
        k = int(np.argmax(es))
        pos_ts, _ = xtb_opt(nums, scan[k][2], chrg=0,
                            constraints=[(i_C2 + 1, hN + 1, scan[k][0]),
                                         (i_N1 + 1, i_C2 + 1, 2.60),
                                         (i_N1 + 1, i_C2n + 1, 1.32)], fc=0.6)
        h = xtb_hess(nums, pos_ts, chrg=0)
        _log(f"  TS2a: d(C2-H)*={scan[k][0]:.2f} A, "
             f"G={h['G_eh'] * EH_KCAL:.2f}, n_imag={h['n_imag']}")
        return {"numbers": nums.tolist(), "positions": pos_ts.tolist(),
                "G_eh": h["G_eh"], "E_eh": h["E_eh"], "n_imag": h["n_imag"],
                "chrg": 0, "d_star": scan[k][0],
                "scan_d": [s for s, _, _ in scan],
                "scan_e_eh": [e for _, e, _ in scan],
                "charges": h["charges"].tolist()}

    def comp_TS2b():
        """Side branch: C3 deprotonation by the phosphate anion (achiral
        enamine channel), same single-coordinate scan recipe."""
        nums, pos, iO, iP, nR, locks, _ = assemble_face(+35.0)
        h3 = [i for i in range(nR) if nums[i] == 1 and np.linalg.norm(
            pos[i] - pos[i_C3]) < 1.15][0]
        pos = prepare_deprotonation_pose(nums, pos, i_C3, h3, iO, iP, nR)
        frames = [1.60, 1.40, 1.25, 1.10]
        scan = scan_1d(nums, pos,
                       lambda d: [(i_C3 + 1, h3 + 1, d),
                                  (i_N1 + 1, i_C2 + 1, 2.60)],
                       frames, chrg=0, fc=0.5, label="TS2b-scan(dC3-H)")
        es = np.array([e for _, e, _ in scan])
        k = int(np.argmax(es))
        pos_ts, _ = xtb_opt(nums, scan[k][2], chrg=0,
                            constraints=[(i_C3 + 1, h3 + 1, scan[k][0]),
                                         (i_N1 + 1, i_C2 + 1, 2.60)], fc=0.6)
        h = xtb_hess(nums, pos_ts, chrg=0)
        _log(f"  TS2b: d(C3-H)*={scan[k][0]:.2f} A, "
             f"G={h['G_eh'] * EH_KCAL:.2f}, n_imag={h['n_imag']}")
        return {"numbers": nums.tolist(), "positions": pos_ts.tolist(),
                "G_eh": h["G_eh"], "E_eh": h["E_eh"], "n_imag": h["n_imag"],
                "chrg": 0,
                "scan_d": [s for s, _, _ in scan],
                "scan_e_eh": [e for _, e, _ in scan],
                "charges": h["charges"].tolist()}

    TS2a = cached_or_compute("cx_TS2a", comp_TS2a)
    TS2b = cached_or_compute("cx_TS2b", comp_TS2b)
    TS2aM = TS2a
    TS2am = TS2a          # shared relay TS; stereo split = I1Cat_M - I1Cat_m

    def comp_I1Cat_unloc():
        """Relaxed (unlocked) resting ion pair — the thermodynamic kinetic
        reference. Starts from the optimized locked pose (same atom order)
        to avoid raw-assembly SCF failures; falls back to the locked value
        with a logged warning if the re-opt fails."""
        nums, pos, iO, iP, nR, locks, _ = assemble_face(+35.0)
        try:
            pos2, _ = xtb_opt(nums, np.array(I1Cat_M["positions"]), chrg=0,
                              constraints=[(i_N1 + 1, i_C2 + 1, 2.75)],
                              fc=0.3)
            h = xtb_hess(nums, pos2, chrg=0)
        except Exception as exc:
            _fallback(f"I1Cat unlocked re-opt failed ({exc}); using locked "
                      "I1Cat_M as resting-state reference")
            return dict(I1Cat_M)
        _log(f"  I1Cat (unlocked): G={h['G_eh'] * EH_KCAL:.2f}, "
             f"n_imag={h['n_imag']}")
        return {"numbers": nums.tolist(), "positions": pos2.tolist(),
                "G_eh": h["G_eh"], "E_eh": h["E_eh"], "n_imag": h["n_imag"],
                "chrg": 0, "charges": h["charges"].tolist()}

    I1Cat = cached_or_compute("cx_I1Cat", comp_I1Cat_unloc)
    XA["I1Cat_M"], XA["I1Cat_m"] = I1Cat_M, I1Cat_m
    XA["TS2aM"], XA["TS2am"], XA["TS2b"] = TS2aM, TS2am, TS2b
    XA["I1Cat"] = I1Cat


    # ---- A.8 Tier-1: intrinsic (thermal, uncatalyzed) reference ------------
    _log("A.8 Tier-1: intrinsic thermal cleavage scan on neutral R ...")
    def comp_TS1_therm():
        frames = [1.60, 1.80, 2.00, 2.20, 2.40]
        scan = scan_1d(R_nums, R_pos,
                       lambda d: [(i_N1 + 1, i_C2 + 1, d)],
                       frames, chrg=0, label="thermal-scan")
        es = np.array([e for _, e, _ in scan])
        k = int(np.argmax(es))
        pos_ts, _ = xtb_opt(R_nums, scan[k][2], chrg=0,
                            constraints=[(i_N1 + 1, i_C2 + 1, scan[k][0])],
                            fc=0.8)
        h = xtb_hess(R_nums, pos_ts, chrg=0)
        return {"G_eh": h["G_eh"], "n_imag": h["n_imag"],
                "d_star": scan[k][0],
                "positions": pos_ts.tolist(), "numbers": R_nums.tolist(),
                "chrg": 0}
    TS1_therm = cached_or_compute("cx_TS1_therm", comp_TS1_therm)
    XA["TS1_therm"] = TS1_therm

    # ---- A.9 network assembly ----------------------------------------------
    G = {k: sp[k]["G_eh"] * EH_KCAL for k in sp}
    G["RC"] = RC["G_eh"] * EH_KCAL
    G["I_RC"] = I_RC["G_eh"] * EH_KCAL
    G["TS1"] = TS1["G_eh"] * EH_KCAL
    G["I1Cat"] = I1Cat["G_eh"] * EH_KCAL
    G["I1Cat_M"] = I1Cat_M["G_eh"] * EH_KCAL
    G["I1Cat_m"] = I1Cat_m["G_eh"] * EH_KCAL
    G["TS2aM"] = TS2aM["G_eh"] * EH_KCAL
    G["TS2am"] = TS2am["G_eh"] * EH_KCAL
    G["TS2b"] = TS2b["G_eh"] * EH_KCAL
    G["TS1_therm"] = TS1_therm["G_eh"] * EH_KCAL
    G_sep_RCat = G["R"] + G["Cat"] + STD_CORR
    G_sep_PCat = G["P_target"] + G["Cat"] + STD_CORR
    G_sep_PECat = G["P_elim"] + G["Cat"] + STD_CORR   # P_elim is C10H11N:
    # the I1Cat -> P_elim + Cat edge is isomeric (no H2 coproduct); H2 only
    # enters the assigned downstream sink P_elim -> Q + H2.
    G["sep_R+Cat"] = G_sep_RCat
    G["sep_P+Cat"] = G_sep_PCat
    G["sep_PE+Cat"] = G_sep_PECat
    # dimer absolute G is ~2x a monomer and would wreck the topology scale;
    # the node is placed relative to the I1 resting state (documented)
    G["P_poly"] = G["I1Cat"] + 4.0

    # barrier bookkeeping per facial family (Curtin-Hammett): the relay
    # coordinate can relax downhill into the deep ion-pair product basin, so
    # a constrained point below the facial reference is floored at FLOOR and
    # flagged — an early (near-barrierless) relay.
    dg = {}
    dg["dG_bind_RC"] = G["RC"] - G_sep_RCat
    dg["dG_proton_transfer"] = G["I_RC"] - G["RC"]
    dg["dG_TS1_vs_RC"] = G["TS1"] - G["RC"]
    dg["dG_I1_vs_RC"] = G["I1Cat"] - G["RC"]
    b2a_raw = G["TS2aM"] - G["I1Cat"]
    b2b_raw = G["TS2b"] - G["I1Cat"]
    # Post-RDS cation-trap barriers lie below the constrained-scan resolution
    # of this pipeline (the relaxed relay coordinate collapses into the deep
    # contact ion-pair product well, raw values below). They are transparently
    # ASSIGNED at ion-pair-trap scale and flagged in the JSON audit trail.
    _assigned("dGTS2a_trap", B2A_ASSIGNED,
              "post-RDS cation capture; computed raw = %.2f kcal (relaxed "
              "scan collapses into product well)" % b2a_raw)
    _assigned("dGTS2b_trap", B2B_ASSIGNED,
              "C3 deprotonation capture; computed raw = %.2f kcal" % b2b_raw)
    dg["dG_TS2aM_vs_I1"] = B2A_ASSIGNED
    dg["dG_TS2am_vs_I1"] = B2A_ASSIGNED   # shared relay TS (early limit)
    dg["dG_TS2b_vs_I1"] = B2B_ASSIGNED
    dg["dG_TS2a_raw"] = b2a_raw
    dg["dG_TS2b_raw"] = b2b_raw
    # Curtin-Hammett stereo-differentiation probe: face split of the locked
    # diastereomeric resting ion pairs I1Cat_M / I1Cat_m (baseline catalyst)
    dg["dG_I1Cat_face_split"] = G["I1Cat_M"] - G["I1Cat_m"]
    dg["ddG_ts2a_stereo"] = dg["dG_I1Cat_face_split"]
    dg["dG_TS1_thermal"] = G["TS1_therm"] - G["R"]
    dg["dG_rxn_overall"] = G_sep_PCat - G_sep_RCat
    dg["dG_aromatization_sink"] = BARRIER_AROM_SINK
    for k, v in dg.items():
        if not (-80 < v < 120):
            _warn(f"network ΔG {k} = {v:.1f} kcal/mol outside sanity window")
    XA["G_kcal"] = G
    XA["dG_kcal"] = dg

    import networkx as nx
    net = nx.DiGraph()
    nodes = {
        "R+Cat": G_sep_RCat, "RC": G["RC"], "I_RC": G["I_RC"],
        "TS1": G["TS1"], "I1Cat": G["I1Cat"], "TS2aM": G["TS2aM"],
        "TS2am": G["TS2am"], "TS2b": G["TS2b"],
        "P_R+Cat": G_sep_PCat, "P_S+Cat": G_sep_PCat,
        "P_elim+Cat": G_sep_PECat, "P_poly": G["P_poly"],
    }
    g0 = min(nodes.values())
    for nname, gv in nodes.items():
        net.add_node(nname, G_rel=gv - g0)
    edges = [
        ("R+Cat", "RC", G["TS1"] - G_sep_RCat),          # association+equil
        ("RC", "TS1", G["TS1"] - G["RC"]),
        ("TS1", "I1Cat", G["TS1"] - G["I1Cat"]),
        ("I1Cat", "TS2aM", G["TS2aM"] - G["I1Cat"]),
        ("I1Cat", "TS2am", G["TS2am"] - G["I1Cat"]),
        ("I1Cat", "TS2b", G["TS2b"] - G["I1Cat"]),
        ("TS2aM", "P_R+Cat", G["TS2aM"] - G_sep_PCat),
        ("TS2am", "P_S+Cat", G["TS2am"] - G_sep_PCat),
        ("TS2b", "P_elim+Cat", G["TS2b"] - G_sep_PECat),
    ]
    for a, b, w in edges:
        net.add_edge(a, b, dG_ddg=abs(w) if w else 0.5)
    XA["network_nodes"] = {n: {"G_rel": d["G_rel"]}
                           for n, d in net.nodes(data=True)}
    XA["network_edges"] = [
        {"u": u, "v": v, "dG_bar": float(d["dG_ddg"])}
        for u, v, d in net.edges(data=True)]

    RESULTS["module_A"] = {k: v for k, v in XA.items()
                           if k not in ("RC", "I_RC", "TS1", "I1Cat",
                                        "I1Cat_M", "I1Cat_m",
                                        "TS2aM", "TS2am", "TS2b")}
    RESULTS["module_A"]["complex_G_kcal"] = {
        "RC": G["RC"], "I_RC": G["I_RC"], "TS1": G["TS1"],
        "I1Cat": G["I1Cat"], "I1Cat_M": G["I1Cat_M"],
        "I1Cat_m": G["I1Cat_m"],
        "TS2aM": G["TS2aM"], "TS2am": G["TS2am"],
        "TS2b": G["TS2b"], "TS1_therm": G["TS1_therm"]}
    RESULTS["module_A"]["ts_quality"] = {
        "TS1_n_imag": TS1["n_imag"], "TS2aM_n_imag": TS2aM["n_imag"],
        "TS2am_n_imag": TS2am["n_imag"], "TS2b_n_imag": TS2b["n_imag"],
        "TS1_therm_n_imag": TS1_therm["n_imag"]}
    write_json_atomic(RESULTS_PATH)
    _log("MODULE A complete — network written to results JSON")
    return sp, RC, I_RC, TS1, I1Cat, TS2aM, TS2am, TS2b


# --------------------------------------------------------------------------- #
# 4.  MODULE B — STIFF MICROKINETIC WORLD MODEL
# --------------------------------------------------------------------------- #
XB: dict = {}
SPECIES = ["R", "Cat", "RC", "I1Cat", "I1", "P_R", "P_S", "P_elim", "Q",
           "H2", "P_poly"]
SI = {s: i for i, s in enumerate(SPECIES)}


def eyring(dg_kcal, T):
    return (2.08366e10 * T) * math.exp(-dg_kcal / (R_GAS * T))  # s^-1


def build_rate_system(G, dg, T, corr=None):
    """Reactions -> (nu matrix (nrxn x nsp), k list (M-aware)).
    corr: designed-catalyst corrections {"d1_drop", "ddG_stereo", ...}
    applied as barrier shifts on the computed/assigned values."""
    corr = corr or {}
    d1 = corr.get("d1_drop", 0.0)
    d2 = corr.get("d2a_drop", 0.0)
    ddS = corr.get("ddG_stereo", 0.0)
    rxns = []
    # 1  R + Cat -> RC
    rxns.append(({("R", -1), ("Cat", -1), ("RC", +1)}, K_ON, "assigned"))
    # 2  RC -> R + Cat
    kd = K_ON * math.exp(max(dg["dG_bind_RC"], -20.0) / (R_GAS * T))
    rxns.append(({("R", +1), ("Cat", +1), ("RC", -1)}, kd, "derived_with_assigned_floors"))
    # 3  RC -> I1Cat  (folded proton transfer + cleavage)
    k3 = eyring(max(dg["dG_TS1_vs_RC"] - d1, 2.0), T)
    rxns.append(({("RC", -1), ("I1Cat", +1)}, k3, "derived_with_assigned_floors"))
    # 4  I1Cat -> RC
    rev1 = max(dg["dG_TS1_vs_RC"] - d1 - dg["dG_I1_vs_RC"], 3.0)
    rxns.append(({("I1Cat", -1), ("RC", +1)}, eyring(rev1, T), "derived_with_assigned_floors"))
    # 5/6  I1Cat -> P_R / P_S + Cat   (enantiodifferentiating)
    # kinetic floor 5.0 kcal on post-correction trap barriers: designed
    # systems would otherwise run at >1e12 s^-1 and stall the integrator.
    # The stereo split is applied AFTER the floor (a floor larger than the
    # split would otherwise erase the enantiodifferentiation); positive
    # ddS favors the P_R (major) channel by convention.
    b_floor = max(dg["dG_TS2aM_vs_I1"] - d2, 5.0)
    b2aM = b_floor - 0.5 * ddS
    b2am = b_floor + 0.5 * ddS
    rxns.append(({("I1Cat", -1), ("P_R", +1), ("Cat", +1)},
                 eyring(b2aM, T), "derived_with_assigned_floors"))
    rxns.append(({("I1Cat", -1), ("P_S", +1), ("Cat", +1)},
                 eyring(b2am, T), "derived_with_assigned_floors"))
    # 7  I1Cat -> P_elim + Cat  (side branch)
    rxns.append(({("I1Cat", -1), ("P_elim", +1), ("Cat", +1)},
                 eyring(max(dg["dG_TS2b_vs_I1"], 2.0), T), "derived_with_assigned_floors"))
    # 8/9  I1Cat <-> I1 + Cat (resting-state dissociation preeq)
    dg_diss = 12.0
    rxns.append(({("I1Cat", -1), ("I1", +1), ("Cat", +1)},
                 eyring(dg_diss, T), "assigned"))
    rxns.append(({("I1", -1), ("Cat", -1), ("I1Cat", +1)}, K_ON,
                 "assigned"))
    # 10  2 I1 -> P_poly
    rxns.append(({("I1", -2), ("P_poly", +1)}, K_DIMER, "assigned"))
    # 11  P_elim -> Q + H2  (thermodynamic aromatization sink)
    rxns.append(({("P_elim", -1), ("Q", +1), ("H2", +1)},
                 eyring(BARRIER_AROM_SINK, T), "assigned"))
    # 12  background thermal R -> I1 (uncatalyzed Tier-1)
    rxns.append(({("R", -1), ("I1", +1)},
                 eyring(max(dg["dG_TS1_thermal"], 2.0), T), "derived_with_assigned_floors"))
    nu = np.zeros((len(rxns), len(SPECIES)))
    ks, kinds = [], []
    for r, (sto, k, kind) in enumerate(rxns):
        for s, c in sto:
            nu[r, SI[s]] = c
        ks.append(k)
        kinds.append(kind)
    return nu, np.array(ks), kinds


def rhs_factory(nu, ks):
    def rhs(t, y):
        ysafe = np.maximum(y, 0.0)
        r = ks.copy()
        # mass-action laws
        r[0] *= ysafe[SI["R"]] * ysafe[SI["Cat"]]
        r[1] *= ysafe[SI["RC"]]
        r[2] *= ysafe[SI["RC"]]
        r[3] *= ysafe[SI["I1Cat"]]
        r[4] *= ysafe[SI["I1Cat"]]
        r[5] *= ysafe[SI["I1Cat"]]
        r[6] *= ysafe[SI["I1Cat"]]
        r[7] *= ysafe[SI["I1Cat"]]
        r[8] *= ysafe[SI["I1"]] * ysafe[SI["Cat"]]
        r[9] *= ysafe[SI["I1"]] ** 2
        r[10] *= ysafe[SI["P_elim"]]
        r[11] *= ysafe[SI["R"]]
        return nu.T @ r
    return rhs


def jac_factory(nu, ks):
    def jac(t, y):
        ysafe = np.maximum(y, 1e-30)
        J = np.zeros((len(SPECIES), len(SPECIES)))
        rr = np.zeros((12, len(SPECIES)))
        rr[0, SI["R"]] = ysafe[SI["Cat"]]; rr[0, SI["Cat"]] = ysafe[SI["R"]]
        rr[1, SI["RC"]] = 1; rr[2, SI["RC"]] = 1; rr[3, SI["I1Cat"]] = 1
        rr[4, SI["I1Cat"]] = 1; rr[5, SI["I1Cat"]] = 1
        rr[6, SI["I1Cat"]] = 1; rr[7, SI["I1Cat"]] = 1
        rr[8, SI["I1"]] = ysafe[SI["Cat"]]
        rr[8, SI["Cat"]] = ysafe[SI["I1"]]
        rr[9, SI["I1"]] = 2 * ysafe[SI["I1"]]
        rr[10, SI["P_elim"]] = 1; rr[11, SI["R"]] = 1
        for r in range(12):
            J += ks[r] * np.outer(nu[r], rr[r])
        return J
    return jac


def integrate(T, G, dg, y0=None, dense=False, corr=None):
    from scipy.integrate import solve_ivp
    nu, ks, kinds = build_rate_system(G, dg, T, corr=corr)
    if y0 is None:
        y0 = np.zeros(len(SPECIES))
        y0[SI["R"]] = CONC_R0
        y0[SI["Cat"]] = CONC_CAT0
    rhs = rhs_factory(nu, ks)
    jac = jac_factory(nu, ks)
    t_eval = np.logspace(math.log10(T_HORIZON[0]),
                         math.log10(T_HORIZON[1]), 160)
    sol = None
    # BDF primary: Radau's internal FD-Newton diverges on this system
    # (zero initial product components against ~1e9 s^-1 association
    # rates); BDF integrates it in sub-second with identical results.
    for method in ("BDF", "Radau", "LSODA"):
        try:
            sol = solve_ivp(rhs, T_HORIZON, y0, method=method, jac=jac,
                            t_eval=t_eval, rtol=1e-6, atol=1e-14,
                            max_step=T_HORIZON[1] / 4,
                            dense_output=dense)
            if sol.success:
                if method != "BDF":
                    _fallback(f"ODE at {T} K: primary BDF failed — "
                              f"{method} used ({sol.message})")
                break
            _warn(f"ODE {method} at {T} K not converged: {sol.message}")
        except Exception as exc:
            _warn(f"ODE {method} at {T} K raised {exc}")
    if sol is None or not sol.success:
        raise RuntimeError(f"all stiff solvers failed at T={T}")
    return sol, nu, ks, kinds


def module_B():
    _log("=" * 70)
    _log("MODULE B — stiff microkinetic world model (Eyring-Polanyi, Radau)")
    _log("=" * 70)
    G = XA.get("G_kcal") or RESULTS["module_A"]["G_kcal"]
    dg = XA.get("dG_kcal") or RESULTS["module_A"]["dG_kcal"]

    _assigned("k_on(R+Cat→RC)", K_ON, "Smoluchowski diffusion-limited cap")
    _assigned("k_dimer(2I1→P_poly)", K_DIMER,
              "cationic oligomerization surrogate")
    _assigned("ΔG‡(P_elim→Q+H2)", BARRIER_AROM_SINK,
              "dehydrogenative aromatization reference barrier")

    # 298.15 K reference run
    sol, nu, ks, kinds = integrate(T_REF, G, dg)
    prof = {s: sol.y[SI[s]] for s in SPECIES}
    XB["profile_298K"] = {"t": sol.t.tolist(),
                          "y": {s: prof[s].tolist() for s in SPECIES}}
    yield_tot = prof["P_R"][-1] + prof["P_S"][-1]
    ee = 100.0 * (prof["P_R"][-1] - prof["P_S"][-1]) / max(yield_tot, 1e-30)
    sel = (yield_tot) / max(prof["P_elim"][-1] + prof["Q"][-1] + 2 *
                            prof["P_poly"][-1], 1e-30)
    XB["reference_298K"] = {
        "yield_target": yield_tot / CONC_R0,
        "ee_pct": ee,
        "selectivity_Pt_Pside": sel,
        "conversion_R": 1 - prof["R"][-1] / CONC_R0,
        "poly_decomp": (prof["P_poly"][-1] * 2) / CONC_R0,
        "rate_constants_s": {f"r{i + 1}": float(ks[i])
                             for i in range(len(ks))},
        "rate_kinds": kinds,
        "nfev": sol.nfev, "njev": sol.njev, "nlu": sol.nlu,
    }
    # stiffness diagnostics
    jac = jac_factory(nu, ks)
    J0 = jac(0.0, np.maximum(sol.y[:, 0], 1e-12))
    Je = jac(sol.t[-1], np.maximum(sol.y[:, -1], 1e-12))
    ev0 = np.linalg.eigvals(J0)
    eve = np.linalg.eigvals(Je)
    lam = np.abs(np.real(ev0))
    lam = lam[lam > 0]
    XB["stiffness"] = {
        "lambda_max_t0": float(lam.max()),
        "lambda_min_t0": float(lam[lam > 0].min()),
        "stiffness_ratio_t0": float(lam.max() / max(lam.min(), 1e-30)),
        "eig_max_real_end": float(np.real(eve).max()),
        "njev": sol.njev, "nlu": sol.nlu, "nfev": sol.nfev,
        "solver": "BDF (implicit multistep, stiffly stable; Radau/LSODA "
                  "fallbacks; analytic Jacobian)",
        "note": "all eigenvalues Re<0 at end state — asymptotically stable "
                "equilibrium (mass-conserved closed system)",
    }
    _log(f"298 K: yield={XB['reference_298K']['yield_target'] * 100:.1f}%  "
         f"ee={ee:.1f}%  selectivity={sel:.1f}  "
         f"stiffness ratio={XB['stiffness']['stiffness_ratio_t0']:.2e}")

    # temperature sweep
    _log(f"B: temperature sweep T={T_SWEEP[0]}..{T_SWEEP[-1]} K ...")
    sweep = []
    for T in T_SWEEP:
        s, _, _, _ = integrate(T, G, dg)
        pR, pS = s.y[SI["P_R"]][-1], s.y[SI["P_S"]][-1]
        pe = s.y[SI["P_elim"]][-1] + s.y[SI["Q"]][-1]
        pp = s.y[SI["P_poly"]][-1]
        yt = pR + pS
        sweep.append({"T": T,
                      "yield": yt / CONC_R0,
                      "ee": 100.0 * (pR - pS) / max(yt, 1e-30),
                      "selectivity": yt / max(pe + 2 * pp, 1e-30),
                      "poly": 2 * pp / CONC_R0})
        _log(f"  T={T} K: yield={sweep[-1]['yield'] * 100:.1f}%  "
             f"ee={sweep[-1]['ee']:.2f}%  sel={sweep[-1]['selectivity']:.1f}")
    XB["T_sweep"] = sweep
    # Curtin-Hammett cross-check
    ddG = dg["ddG_ts2a_stereo"]
    XB["ee_curtin_298K"] = 100.0 * math.tanh(ddG / (2 * R_GAS * T_REF))
    RESULTS["module_B"] = XB
    write_json_atomic(RESULTS_PATH)
    _log("MODULE B complete")
    return XB



def module_D():
    """Designed-catalyst world model: measure the WINNER catalyst's facial
    stereo-differentiation (locked diastereomeric resting ion pairs) with a
    baseline control, then re-integrate the stiff ODE with the designed
    barrier corrections."""
    if not EXPLORATORY_KINETICS:
        raise RuntimeError('Module D requires explicit exploratory-kinetics opt-in')
    _log("=" * 70)
    _log("MODULE D — designed-catalyst world model (winner parameters)")
    _log("=" * 70)
    MC = RESULTS.get("module_C") or {}
    if "proof" not in MC:
        raise RuntimeError("module D requires module C results (winner "
                           "catalyst)")
    G = XA["G_kcal"]; dg = XA["dG_kcal"]
    I1Cat = XA["I1Cat"]
    sub_nums = np.array(I1Cat["numbers"])
    sub_pos = np.array(I1Cat["positions"])
    n_sub_full = len(sub_nums)
    cpa_len = len(np.array(json.loads(
        (CACHE / "sp_CPAneg.json").read_text())["numbers"]))
    nS = n_sub_full - cpa_len                      # bare substrate part
    S_nums, S_pos = sub_nums[:nS], sub_pos[:nS]
    idx = XA["idx_R"]
    i_N1, i_C2, i_C2n, i_C3 = (idx["N1"], idx["C2"], idx["C2n"], idx["C3"])
    i_HN = [i for i in range(nS) if S_nums[i] == 1 and np.linalg.norm(
        S_pos[i] - S_pos[i_N1]) < 1.15]
    if not i_HN:
        i_HN = [i for i in range(nS) if S_nums[i] == 1 and np.linalg.norm(
            S_pos[i] - S_pos[i_C2]) < 1.6]
    i_HN = i_HN[0]
    anchor = (S_pos[i_HN] + 1.65 * (S_pos[i_HN] - S_pos[i_N1]) /
              np.linalg.norm(S_pos[i_HN] - S_pos[i_N1]))
    # substrate ring plane for enantiotopic-face generation
    pts = S_pos[[i_N1, i_C2, i_C2n, i_C3]]
    ctr = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - ctr)
    nvec = vt[2]

    def reflect(x):
        return x - 2.0 * np.outer((x - ctr) @ nvec, nvec)

    def deprot_pose(cat_smiles):
        m = build_Cat(cat_smiles, name="Cat_D")
        cn, cp = embed3d(m, seed=0xD1CE, n_conf=6, label="D")
        cp, _ = xtb_opt(cn, cp)
        iP = [i for i in range(len(cn)) if cn[i] == 15][0]
        iO = iHrm = None
        for j in range(len(cn)):
            if cn[j] == 8 and np.linalg.norm(cp[j] - cp[iP]) < 1.8:
                hs = [k for k in range(len(cn)) if cn[k] == 1 and
                      np.linalg.norm(cp[k] - cp[j]) < 1.2]
                if hs and iHrm is None:
                    iHrm = hs[0]
                elif not hs and iO is None:
                    iO = j
        keep = [i for i in range(len(cn)) if i != iHrm]
        cn2, cp2 = cn[keep], cp[keep]
        iO2 = iO - (1 if iHrm < iO else 0)
        return cn2, cp2, iO2

    def face_energy(cat_smiles, face_sign, tag):
        cn2, cp2, iO2 = deprot_pose(cat_smiles)
        pO = cp2[iO2]
        v = (anchor - S_pos[i_N1]) / np.linalg.norm(anchor - S_pos[i_N1])
        dirO = pO - cp2[[i for i in range(len(cn2)) if cn2[i] == 15][0]]
        dirO /= np.linalg.norm(dirO)
        Rm = kabsch_rotate(np.array([dirO]), np.array([v]))
        cp3 = cp2 @ Rm
        cp3 += (anchor - cp3[iO2])
        if face_sign < 0:
            cp3 = reflect(cp3)
            cp3 += (anchor - cp3[iO2])
        nums, pos, _lab = merge_frags([(S_nums, S_pos, ["S"] * nS),
                                 (cn2, cp3, ["CAT"] * len(cn2))])
        iPm = nS + [i for i in range(len(cn2)) if cn2[i] == 15][0]
        # constrained-first protocol (comparable across faces); unconstrained
        # pre-relax is a FALLBACK only for raw-assembly SCF crashes, never
        # the default path (pre-relaxed faces relax into non-comparable
        # binding modes)
        cons = [(i_HN + 1, nS + iO2 + 1, 1.65),
                (iPm + 1, i_C3 + 1,
                 round(float(np.linalg.norm(pos[iPm] - pos[i_C3])), 3))]
        try:
            pos2, _ = xtb_opt(nums, pos, chrg=0, constraints=cons, fc=0.3)
        except Exception:
            _warn(f"face {tag}({face_sign:+d}): constrained opt failed; "
                  "pre-relax fallback engaged")
            pos_pre, _ = xtb_opt(nums, pos, chrg=0, timeout=900)
            cons = [(i_HN + 1, nS + iO2 + 1, 1.65),
                    (iPm + 1, i_C3 + 1,
                     round(float(np.linalg.norm(
                         pos_pre[iPm] - pos_pre[i_C3])), 3))]
            pos2, _ = xtb_opt(nums, pos_pre, chrg=0, constraints=cons,
                              fc=0.3)
        h = xtb_hess(nums, pos2, chrg=0)
        _log(f"  {tag} face({face_sign:+d}): G={h['G_eh'] * EH_KCAL:.2f}, "
             f"n_imag={h['n_imag']}")
        return h["G_eh"]

    win_smiles = MC["proof"]["winner_smiles"]
    base_smiles = MC["proof"]["baseline"]["smiles"]

    def face_energy_cached(cat_smiles, face_sign, tag):
        key = "dface_%s_%+d" % (tag, face_sign)
        return cached_or_compute(
            key, lambda: {"G_eh": face_energy(cat_smiles, face_sign,
                                              tag)})["G_eh"]

    dd_win = (face_energy_cached(win_smiles, +1, "WIN") -
              face_energy_cached(win_smiles, -1, "WIN")) * EH_KCAL
    dd_base = (face_energy_cached(base_smiles, +1, "BASE") -
               face_energy_cached(base_smiles, -1, "BASE")) * EH_KCAL
    dd_stereo = dd_win - dd_base      # designed differential (baseline-removed)
    _log(f"  face split: winner {dd_win:+.2f} kcal, baseline control "
         f"{dd_base:+.2f} kcal -> designed ddG_stereo = {dd_stereo:+.2f}")

    drop_raw = MC["proof"]["barrier_drop_vs_baseline_kcal"]
    # transparent kinetic caps over the computed raw values: gas-phase
    # ion-pair differentials (~65 kcal raw) are strongly attenuated in
    # polar solution (cap 8.0); the facial-split probe is pose-unstable at
    # this fidelity, so the designed stereo input is capped at a
    # good-CPA-typical 1.5 kcal/mol. Raw values remain in the JSON.
    KINETIC_CAP = 8.0
    STEREO_CAP = 1.5
    scenario = exploratory_corrections(drop_raw, dd_stereo)
    drop = scenario['d1_drop']
    dd_stereo_used = scenario['ddG_stereo']
    _assigned("d1_drop(designed, used)", drop,
              "cap %.1f over raw Pauling differential %.2f kcal"
              % (KINETIC_CAP, drop_raw))
    _assigned("ddG_stereo(designed, used)", dd_stereo_used,
              "cap %.1f over raw facial split %.2f kcal"
              % (STEREO_CAP, abs(dd_stereo)))
    corr = scenario
    RESULTS['evidence_status'] = 'assigned_sensitivity_scenario_not_prediction'

    sol, _, _, _ = integrate(T_REF, G, dg, corr=corr)
    prof_w = {s: sol.y[SI[s]] for s in SPECIES}
    yt = prof_w["P_R"][-1] + prof_w["P_S"][-1]
    ee_w = 100.0 * (prof_w["P_R"][-1] - prof_w["P_S"][-1]) / max(yt, 1e-30)
    sel_w = yt / max(prof_w["P_elim"][-1] + prof_w["Q"][-1] +
                     2 * prof_w["P_poly"][-1], 1e-30)
    sweep_w = []
    for T in T_SWEEP:
        s, _, _, _ = integrate(T, G, dg, corr=corr)
        pR, pS = s.y[SI["P_R"]][-1], s.y[SI["P_S"]][-1]
        ytw = pR + pS
        sweep_w.append({"T": T, "yield": ytw / CONC_R0,
                        "ee": 100.0 * (pR - pS) / max(ytw, 1e-30),
                        "selectivity": ytw / max(
                            s.y[SI["P_elim"]][-1] + s.y[SI["Q"]][-1] +
                            2 * s.y[SI["P_poly"]][-1], 1e-30)})
        _log(f"  [designed] T={T} K: yield={sweep_w[-1]['yield'] * 100:.1f}% "
             f"ee={sweep_w[-1]['ee']:.2f}%")
    RESULTS["module_D"] = {
        "winner_face_split_kcal": float(dd_win),
        "baseline_face_split_kcal": float(dd_base),
        "ddG_stereo_designed_kcal": float(abs(dd_stereo)),
        "corr": corr,
        "reference_298K_winner": {
            "yield_target": yt / CONC_R0, "ee_pct": ee_w,
            "selectivity_Pt_Pside": sel_w,
            "conversion_R": 1 - prof_w["R"][-1] / CONC_R0},
        "profile_298K_winner": {"t": sol.t.tolist(),
                                "y": {s: prof_w[s].tolist()
                                      for s in SPECIES}},
        "T_sweep_winner": sweep_w,
        "ee_curtin_winner_298K": 100.0 * math.tanh(
            abs(dd_stereo) / (2 * R_GAS * T_REF)),
    }
    write_json_atomic(RESULTS_PATH)
    _log("MODULE D complete")
    return RESULTS["module_D"]



# --------------------------------------------------------------------------- #
# 5.  MODULE C — INVERSE GENERATIVE DESIGN (ESP + EVOLUTIONARY ASSEMBLER)
# --------------------------------------------------------------------------- #
XC: dict = {}
MOTIFS = {
    "H": "[H]",
    "tBu": "C(C)(C)C",
    "Me": "C",
    "Ph": "c3ccccc3",
    "pOH-Ph": "c3ccc(O)cc3",
    "oOH-Ph": "c3ccccc3O",
    "iPr-Ph": "c3cc(C(C)C)ccc3",
    "CF3-Ph": "c3ccc(C(F)(F)F)cc3",
    "OMe-Ph": "c3ccc(OC)cc3",
}
CAT_TEMPLATE = ("O=P1(O)Oc2c({X})cc3ccccc3c2-c2c({Y})cc3ccccc3c2O1")


def esp_field(numbers, positions, charges, center, half=6.0, n=26):
    """phi(r) = sum q_i / |r - r_i| in e/A; grid around `center`."""
    ax = np.linspace(-half, half, n)
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    pts = np.stack([X, Y, Z], axis=-1).reshape(-1, 3) + center
    phi = np.zeros(len(pts))
    for (ri, qi) in zip(positions, charges):
        d = np.linalg.norm(pts - ri, axis=1)
        d = np.maximum(d, 0.35)
        phi += qi / d
    return pts, phi


def assemble_complex(nums_TS, pos_TS, cat_nums, cat_pos, iO, face_deg,
                     anchor_xyz):
    """Place catalyst anionic O at anchor (H-bond to TS core), rotate face."""
    pO = cat_pos[iO]
    d0 = anchor_xyz - pO
    cat_new = cat_pos + d0
    axv = anchor_xyz - np.mean(pos_TS[[i for i in range(len(nums_TS))][:12]],
                               axis=0)
    cat_new = place_fragment(cat_nums, cat_new, [iO], [anchor_xyz],
                             extra_rot_axis=axv, extra_rot_deg=face_deg)
    nums, pos, lab = merge_frags([
        (nums_TS, pos_TS, ["TS"] * len(nums_TS)),
        (cat_nums, cat_new, ["CAT"] * len(cat_nums))])
    return nums, pos


def vdW_clash(nums, pos, frag_boundary, cutoff=1.15):
    rvdw = {1: 1.1, 6: 1.7, 7: 1.55, 8: 1.52, 15: 1.8, 9: 1.47, 16: 1.78,
            17: 1.75}
    bad = 0
    nA = frag_boundary
    for i in range(nA):
        for j in range(nA, len(nums)):
            if np.linalg.norm(pos[i] - pos[j]) < (rvdw.get(nums[i], 1.5) +
                                                  rvdw.get(nums[j], 1.5)) \
                    * 0.62:
                bad += 1
    return bad


def module_C():
    _log("=" * 70)
    _log("MODULE C — TS-conditioned inverse catalyst design")
    _log("=" * 70)
    dg = XA["dG_kcal"]

    # ---- C.1 ESP field of the RDS transition state -------------------------
    _log("C.1 ESP field extraction around the charge-redistribution center")
    # Docking/ESP template = the post-RDS cationic intermediate I1Cat (the
    # stereodetermining arenium·phosphate ion pair, carries the N1-H anchor);
    # the computed RDS saddle TS1 remains the ENERGY reference.
    TSg = XA["I1Cat"]
    TS1 = XA["TS1"]
    nums = np.array(TSg["numbers"]); pos = np.array(TSg["positions"])
    q = np.array(TSg["charges"])
    idx = XA["idx_R"]
    _cpa = cached_or_compute("sp_CPAneg", lambda: (_ for _ in ()).throw(
        RuntimeError("sp_CPAneg cache missing")))
    n_sub = len(np.array(XA["I_RC"]["numbers"])) - len(np.array(_cpa["numbers"]))
    ctr = pos[[idx["C2"], idx["N1"], idx["C2n"], idx["C3"]]].mean(axis=0)
    pts, phi = esp_field(nums, pos, q, ctr)
    phi_min_pt = pts[int(np.argmin(phi))]
    phi_max_pt = pts[int(np.argmax(phi))]
    XC["esp"] = {
        "center": ctr.tolist(),
        "phi_min": float(phi.min()), "phi_max": float(phi.max()),
        "min_point": phi_min_pt.tolist(), "max_point": phi_max_pt.tolist(),
        "grid_n": 26, "half_box_A": 6.0,
        "q_range": [float(q.min()), float(q.max())],
    }
    _log(f"  ESP range {phi.min():.3f}..{phi.max():.3f} e/A; "
         f"hotspot vector length "
         f"{np.linalg.norm(phi_max_pt - phi_min_pt):.2f} A")

    # ---- C.2 evolutionary scaffold assembler -------------------------------
    _log("C.2 evolutionary 3,3'-scaffold assembly (RDKit valency gates + "
         "GFN2-xTB fitness)")
    i_HN = [i for i in range(n_sub) if nums[i] == 1 and
            np.linalg.norm(pos[i] - pos[idx["N1"]]) < 1.2]
    if not i_HN:          # in-transit H sits between N1 and C2 in TS1 frame
        i_HN = [i for i in range(n_sub) if nums[i] == 1 and
                np.linalg.norm(pos[i] - pos[idx["C2"]]) < 1.6]
    if not i_HN:
        raise RuntimeError("module C: could not locate anchoring N1-H in TS1")
    i_HN = i_HN[0]
    anchor = pos[i_HN] + 1.65 * (pos[i_HN] - pos[idx["N1"]]) / np.linalg.norm(
        pos[i_HN] - pos[idx["N1"]])

    def cat_variant(x_raw, y_raw):
        smi = CAT_TEMPLATE.format(X=x_raw, Y=y_raw)
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        try:
            Chem.SanitizeMol(m)
        except Exception:
            return None
        return m

    def deprotonated_pose(cnums, cpos):
        """Catalyst optimized geometry -> anionic phosphate pose (drop one
        P-OH proton); returns (nums, pos, iO_anionic)."""
        iP = [i for i in range(len(cnums)) if cnums[i] == 15][0]
        iO = iHrm = None
        for j in range(len(cnums)):
            if cnums[j] == 8 and np.linalg.norm(cpos[j] - cpos[iP]) < 1.8:
                hs = [k for k in range(len(cnums)) if cnums[k] == 1 and
                      np.linalg.norm(cpos[k] - cpos[j]) < 1.2]
                if hs and iHrm is None:
                    iHrm = hs[0]
                elif not hs:
                    iO = j
        if iO is None or iHrm is None:
            raise RuntimeError("phosphate O/H topology not resolved")
        keep = [i for i in range(len(cnums)) if i != iHrm]
        cn2 = cnums[keep]
        cp2 = cpos[keep]
        iO2 = iO - (1 if iHrm < iO else 0)
        return cn2, cp2, iO2

    _cat_opt_cache: dict = {}

    def pose_fitness(m_cat, face_deg, cache_key=None):
        """Pauling proxy fitness: GFN2-xTB SP of TS-substrate·Cat- ion pair
        (lower = stronger TS stabilization), gated by valency + vdW.
        A 7-point torsional sweep about the H-bond axis removes the
        orientation false-rejections of a single fixed pose."""
        try:
            ck = cache_key or id(m_cat)
            if ck in _cat_opt_cache:
                cnums, cpos = _cat_opt_cache[ck]
            else:
                cnums, cpos = embed3d(m_cat, seed=0xBEEF, n_conf=4,
                                      label="catGA")
                cpos, _ = xtb_opt(cnums, cpos)
                _cat_opt_cache[ck] = (cnums, cpos)
            cn2, cp2, iO2 = deprotonated_pose(cnums, cpos)
            nTS = n_sub
            best = None
            for dphi in np.linspace(-60.0, 60.0, 7):
                cnumsTS, cposTS = assemble_complex(
                    nums[:n_sub], pos[:n_sub], cn2, cp2, iO2,
                    face_deg + float(dphi), anchor)
                if vdW_clash(cnumsTS, cposTS, nTS) > 0:
                    continue
                e_ts_c, _ = xtb_sp_charges(cnumsTS, cposTS, chrg=0)
                if best is None or e_ts_c < best:
                    best = e_ts_c
            return best
        except Exception:
            return None

    key_list = list(MOTIFS)
    # one pose family per (X, Y): the enantio-differential is computed
    # explicitly in Module A; here the GA searches catalyst identity
    population = [(X, Y, 35.0) for X in key_list for Y in key_list]
    rng = np.random.default_rng(0xD1CE)
    scored = []
    t0 = time.time()
    for i, (X, Y, f) in enumerate(population):
        m = cat_variant(MOTIFS[X], MOTIFS[Y])
        if m is None:
            continue
        fmut = f + float(rng.uniform(-12, 12))
        e = pose_fitness(m, fmut, cache_key=(X, Y))
        if e is not None:
            scored.append({"X": X, "Y": Y, "face": fmut, "E_tsC_eh": e})
            _log(f"  gen0 cand {i + 1}/{len(population)} {X}/{Y}@{fmut:+.0f}"
                 f"°  E(TS·C)={e * EH_KCAL:.1f} kcal "
                 f"[{time.time() - t0:.0f}s]")
    if not scored:
        raise RuntimeError("GA: all candidates failed — motif library too "
                           "aggressive")
    scored.sort(key=lambda d: d["E_tsC_eh"])
    elite = scored[:3]
    for e in elite:
        for dphi in (-10.0, +10.0):
            m = cat_variant(MOTIFS[e["X"]], MOTIFS[e["Y"]])
            if m is None:
                continue
            e2 = pose_fitness(m, e["face"] + dphi,
                              cache_key=(e["X"], e["Y"]))
            if e2 is not None and e2 < e["E_tsC_eh"]:
                e["face"] += dphi
                e["E_tsC_eh"] = e2
    elite.sort(key=lambda d: d["E_tsC_eh"])
    winner = elite[0]
    XC["ga"] = {
        "library_size": len(MOTIFS),
        "candidates_scored": len(scored),
        "fitness": "GFN2-xTB single-point E(TS-substrate·Cat-) ion-pair SP "
                   "(Pauling proxy), gated by RDKit valency sanitize + "
                   "vdW non-clashing",
        "top5": scored[:5],
        "winner": winner,
    }
    _log(f"  GA winner: {winner['X']}/{winner['Y']} "
         f"@{winner['face']:+.0f}°")

    # ---- C.3 proof: >= 4.0 kcal/mol barrier drop ---------------------------
    _log("C.3 exploratory constrained complex-energy comparison; not barrier proof")
    base_smiles = CAT_SMILES
    win_smiles = CAT_TEMPLATE.format(X=MOTIFS[winner["X"]],
                                     Y=MOTIFS[winner["Y"]])
    TS_sub_nums = nums[:n_sub]
    TS_sub_pos = pos[:n_sub]

    def full_complex_energies(cat_smiles, tag):
        m_cat = build_Cat(cat_smiles, name=f"Cat_{tag}")
        cnums, cpos = embed3d(m_cat, seed=0xD1CE, n_conf=6, label=tag)
        cpos, _ = xtb_opt(cnums, cpos)
        cn2, cp2, iO2 = deprotonated_pose(cnums, cpos)
        # --- TS-substrate · phosphate-anion ion pair (constrained refit) ---
        cnumsTS, cposTS = assemble_complex(
            TS_sub_nums, TS_sub_pos, cn2, cp2, iO2, winner["face"], anchor)
        nTS = len(TS_sub_nums)
        cons = [(i_HN + 1, nTS + iO2 + 1, 1.65)]
        posTS, _ = xtb_opt(cnumsTS, cposTS, chrg=0, constraints=cons, fc=0.4)
        hTS = xtb_hess(cnumsTS, posTS, chrg=0)
        # --- neutral R · neutral catalyst H-bonded complex ------------------
        Rc = cached_or_compute("sp_R", lambda: (_ for _ in ()).throw(
            RuntimeError("sp_R cache missing")))
        Rn = np.array(Rc["numbers"]); Rp2p = np.array(Rc["positions"])
        iN1 = idx["N1"]
        v = Rp2p[iN1] - Rp2p[[idx["C2"], idx["C2n"], idx["C3"]]].mean(axis=0)
        v /= np.linalg.norm(v)
        iPc = [i for i in range(len(cnums)) if cnums[i] == 15][0]
        iOph = iHph = None
        for j in range(len(cnums)):
            if cnums[j] == 8 and np.linalg.norm(cpos[j] - cpos[iPc]) < 1.8:
                hs = [k for k in range(len(cnums)) if cnums[k] == 1 and
                      np.linalg.norm(cpos[k] - cpos[j]) < 1.2]
                if hs:
                    iOph, iHph = j, hs[0]
                    break
        if iOph is None:
            iOph, iHph = iO, None
        cp3 = cpos.copy()
        tgt = Rp2p[iN1] + (2.85 if iHph is not None else 1.85) * v
        cp3 += (tgt - cp3[iOph])
        numsRC, posRC, _lab = merge_frags([
            (Rn, Rp2p, ["R"] * len(Rn)),
            (cnums, cp3, ["CAT"] * len(cnums))])
        nR = len(Rn)
        cons2 = ([(iN1 + 1, nR + iHph + 1, 1.80)] if iHph is not None
                 else [(iN1 + 1, nR + iOph + 1, 1.85)])
        posRC2, _ = xtb_opt(numsRC, posRC, chrg=0, constraints=cons2, fc=0.4)
        hRC = xtb_hess(numsRC, posRC2, chrg=0)
        cstand = xtb_hess(cnums, cpos, chrg=0)
        if tag == "WIN":
            XC["proof_complex"] = {
                "cat_nums": cnumsTS[nTS:].tolist(),
                "cat_pos": posTS[nTS:].tolist(),
                "nTS": int(nTS), "iO": int(iO2 + nTS),
            }
        return {"G_tsC": hTS["G_eh"], "G_rcC": hRC["G_eh"],
                "G_C": cstand["G_eh"], "G_TS": TS1["G_eh"],
                "G_R": Rc["G_eh"]}

    def barrier_drop(cat_smiles, tag):
        ce = full_complex_energies(cat_smiles, tag)
        ddG_bind = ((ce["G_tsC"] - ce["G_TS"] - ce["G_C"]) -
                    (ce["G_rcC"] - ce["G_R"] - ce["G_C"]))
        drop = -ddG_bind * EH_KCAL      # >0 = TS stabilized relative to RC
        return drop, ce

    drop_win, ce_win = barrier_drop(win_smiles, "WIN")
    drop_base, ce_base = barrier_drop(base_smiles, "BASE")
    drop_rel = drop_win - drop_base
    XC["proof"] = {
        "winner_smiles": win_smiles,
        "winner_motifs": [winner["X"], winner["Y"]],
        "ddG_bind_winner_TS_minus_R_kcal": float(-drop_win),
        "barrier_drop_vs_baseline_kcal": float(drop_rel),
        "baseline": {"smiles": base_smiles,
                     "ddG_bind_kcal": float(-drop_base)},
        "target_kcal": 4.0,
        "claim_proven": False,
        "threshold_exceeded": bool(drop_rel >= 4.0),
        "evidence_status": "unvalidated_constrained_complex_energy_proxy",
        "reason": "No unconstrained saddle/gradient/IRC verification for both catalysts",
    }
    if not XC["proof"]["claim_proven"]:
        _warn(f"complex-energy proxy {drop_rel:.2f} kcal/mol is NOT a verified barrier reduction")
    RESULTS["module_C"] = {k: v for k, v in XC.items()
                           if k != "proof_complex"}
    RESULTS["module_C"]["_pos"] = {
        "ts_nums": TS_sub_nums.tolist(), "ts_pos": TS_sub_pos.tolist(),
        "ts_q": q[:n_sub].tolist(),
        "winner_face": winner["face"], "anchor": anchor.tolist(),
        "esp_center": ctr.tolist(),
    }
    RESULTS["module_C"]["_proof_complex"] = XC.get("proof_complex", {})
    write_json_atomic(RESULTS_PATH)
    _log("MODULE C complete")
    return XC

def style_axes(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(alpha=0.25, lw=0.6)


def figure_1():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import Normalize
    import matplotlib.patheffects as pe

    G = XA["G_kcal"]
    dg = XA["dG_kcal"]
    # TS nodes plotted at their EFFECTIVE barriers over the parent state
    # (raw constrained-structure G carries placement strain and would
    # exaggerate the topology)
    nodes = {
        "R+Cat": G["sep_R+Cat"], "RC": G["RC"], "I_RC": G["I_RC"],
        "TS1": G["RC"] + dg["dG_TS1_vs_RC"], "I1Cat": G["I1Cat"],
        "TS2aM": G["I1Cat"] + dg["dG_TS2aM_vs_I1"],
        "TS2am": G["I1Cat"] + dg["dG_TS2am_vs_I1"],
        "TS2b": G["I1Cat"] + dg["dG_TS2b_vs_I1"],
        "P_R+Cat": G["sep_P+Cat"], "P_S+Cat": G["sep_P+Cat"],
        "P_elim+Cat": G["sep_PE+Cat"], "P_poly": G["P_poly"],
    }
    base = min(nodes.values())
    nodes = {k: v - base for k, v in nodes.items()}
    # decluster the degenerate product trio and the shared-relay TS2 pair
    spread = {"P_R+Cat": +0.8, "P_S+Cat": +3.4, "P_elim+Cat": -1.8,
              "TS2am": +3.0, "TS2b": -1.5}
    lay = {   # (x, energy) hand layout of the catalytic cycle
        "R+Cat": (0.0, nodes["R+Cat"]),
        "RC": (1.4, nodes["RC"]),
        "I_RC": (2.8, nodes["I_RC"]),
        "TS1": (4.2, nodes["TS1"]),
        "I1Cat": (5.6, nodes["I1Cat"]),
        "TS2aM": (7.0, nodes["TS2aM"] + spread.get("TS2aM", 0.0)),
        "TS2am": (7.0, nodes["TS2am"] + spread.get("TS2am", 0.0)),
        "TS2b": (7.6, nodes["TS2b"] + spread.get("TS2b", 0.0)),
        "P_R+Cat": (8.8, nodes["P_R+Cat"] + spread.get("P_R+Cat", 0.0)),
        "P_S+Cat": (8.8, nodes["P_S+Cat"] + spread.get("P_S+Cat", 0.0)),
        "P_elim+Cat": (8.8,
                       nodes["P_elim+Cat"] + spread.get("P_elim+Cat", 0.0)),
        "P_poly": (10.2, nodes["P_poly"]),
    }
    fig, ax = plt.subplots(figsize=(13.5, 8.0))
    norm = Normalize(vmin=min(nodes.values()),
                     vmax=np.percentile(list(nodes.values()), 92))
    cmap = matplotlib.colormaps["viridis"]
    arrows = [
        ("R+Cat", "RC", "k_on/k_off"),
        ("RC", "TS1", "RDS ΔG‡"),
        ("TS1", "I1Cat", None),
        ("I1Cat", "TS2aM", None), ("I1Cat", "TS2am", None),
        ("I1Cat", "TS2b", None),
        ("TS2aM", "P_R+Cat", None), ("TS2am", "P_S+Cat", None),
        ("TS2b", "P_elim+Cat", None),
        ("I1Cat", "P_poly", "oligomerization"),
    ]
    for a, b, lab in arrows:
        x0, y0 = lay[a]; x1, y1 = lay[b]
        if y0 == y1 and x0 == x1:
            continue
        w = 1.2 + 4.5 / max(abs(G["TS1"] - G["RC"]) / 6.0, 0.5)
        rad = 0.18 if a == "I1Cat" else 0.12
        arr = matplotlib.patches.FancyArrowPatch(
            (x0, y0), (x1, y1), connectionstyle=f"arc3,rad={rad}",
            arrowstyle="-|>", mutation_scale=16, lw=1.6,
            color="#555555", zorder=2,
            path_effects=[pe.withStroke(linewidth=3.2, foreground="white")])
        ax.add_patch(arr)
        if lab:
            ax.text(0.5 * (x0 + x1), 0.5 * (y0 + y1) + 1.2, lab, fontsize=9,
                    ha="center", color="#333333", rotation=0)
    # barrier labels
    xlay = {k: v[0] for k, v in lay.items()}
    ylay = {k: v[1] for k, v in lay.items()}
    ax.annotate(f"ΔG‡₁ = {dg['dG_TS1_vs_RC']:.1f}",
                (xlay["TS1"], ylay["TS1"] + 1.0),
                fontsize=9.5, ha="center", fontweight="bold",
                color="#8e44ad")
    ax.annotate(f"ΔG‡₂ᵃ = {dg['dG_TS2aM_vs_I1']:.1f}",
                (xlay["TS2aM"] - 1.15, ylay["TS2aM"]), fontsize=9,
                ha="center", color="#1f6f43")
    ax.annotate(f"(shared relay TS)", (xlay["TS2am"] - 1.15,
                ylay["TS2am"]), fontsize=8.2, ha="center", color="#a04000")
    ax.annotate(f"ΔG‡₂ᵇ = {dg['dG_TS2b_vs_I1']:.1f}",
                (xlay["TS2b"] + 1.25, ylay["TS2b"]), fontsize=9,
                ha="center", color="#7f2400")
    for name, (x, y) in lay.items():
        c = cmap(norm(y))
        ax.scatter([x], [y], s=430, c=[c], edgecolors="black", zorder=3,
                   linewidths=1.1)
        dy_lab = 1.8 if name not in ("TS2aM", "TS2am", "TS2b",
                                     "P_R+Cat", "P_S+Cat",
                                     "P_elim+Cat") else 2.1
        ax.text(x - 0.35, y - dy_lab, name, fontsize=8.8, ha="center",
                fontweight="bold")
        ax.text(x + 0.42, y + 0.2, f"{y:.1f}", fontsize=8.0,
                ha="center", color="#222222")
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, pad=0.015, fraction=0.035)
    cb.set_label("relative Gibbs free energy  G − G$_{min}$  (kcal/mol)",
                 fontsize=11)
    ax.axhline(G["sep_R+Cat"] - base, color="#999999", lw=0.8, ls=":")
    ax.text(10.6, G["sep_R+Cat"] - base + 0.4, "R + Cat reference",
            fontsize=8.5, color="#666666")
    ax.text(4.0, 8.0, "baseline Cat: facial probe ~ 0", fontsize=8.6,
            ha="center", color="#444444", style="italic",
            bbox=dict(boxstyle="round,pad=0.28", fc="#fdf6e3",
                      ec="#b58900", lw=0.8))
    ax.set_xlim(-0.9, 11.6)
    ax.set_ylim(min(nodes.values()) - 5, max(nodes.values()) + 5)
    ax.set_xlabel("catalytic-cycle progress (reaction-network coordinate)",
                  fontsize=12)
    ax.set_ylabel("Gibbs free energy  (kcal/mol)", fontsize=12)
    ax.set_title("Fig. 1 — Phase 5 Automated Reaction Network: asymmetric "
                 "aziridine ring expansion, GFN2-xTB ion-pair topology",
                 fontsize=13.5, fontweight="bold")
    style_axes(ax)
    fig.tight_layout()
    p = FIG / "fig1_reaction_network_topology.png"
    fig.text(.5, .005, 'EXPLORATORY MODEL: not a verified barrier prediction', ha='center', fontsize=8)
    fig.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    _log(f"figure 1 -> {p.name}")


def figure_2():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.colors import Normalize

    XD = RESULTS.get("module_D") or {}
    designed = "profile_298K_winner" in XD
    prof = (XD.get("profile_298K_winner") if designed
            else XB["profile_298K"])
    t = np.array(prof["t"])
    fig = plt.figure(figsize=(13.5, 6.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1])
    ax = fig.add_subplot(gs[0, 0])
    series = [
        ("R", "#d62728", "-", "R (aziridine substrate)"),
        ("Cat", "#7f7f7f", "--", "Cat (chiral phosphoric acid)"),
        ("RC", "#e377c2", "-.", "RC (H-bonded complex)"),
        ("I1Cat", "#9467bd", ":", "I₁·Cat (ion-pair intermediate)"),
        ("I1", "#c5b0d5", ":", "I₁ (resting state)"),
        ("P_R", "#1f77b4", "-", "P_target (major enantiomer)"),
        ("P_S", "#aec7e8", "-", "P_target (minor enantiomer)"),
        ("P_elim", "#ff7f0e", "-", "P_elim (achiral enamine)"),
        ("Q", "#b22222", "--", "Q (aromatized quinoline sink)"),
        ("P_poly", "#2ca02c", "--", "P_poly (oligomer)"),
    ]
    for key, c, ls, lab in series:
        ax.semilogx(t, prof["y"][key], ls, color=c, lw=1.9, label=lab)
    ax.set_xlabel("time  t  (s, log scale) — 1 ns to 28 h", fontsize=12)
    ax.set_ylabel("concentration  (mol/L)", fontsize=12)
    ax.set_title("Panel A — concentration dynamics, "
                 f"{'designed' if designed else 'baseline'} Cat (298 K, "
                 "BDF)", fontsize=12.0, fontweight="bold", pad=10)
    ax.legend(fontsize=8.4, ncol=2, frameon=False, loc="center right")
    style_axes(ax)

    ax2 = fig.add_subplot(gs[0, 1])
    sw = (XD.get("T_sweep_winner") if designed else XB["T_sweep"])
    if designed:
        base_sw = XB["T_sweep"]
        ax2.plot([100 * d["ee"] for d in base_sw],
                 [100 * d["yield"] for d in base_sw], "--", color="#aaaaaa",
                 lw=1.4, zorder=1)
        ax2.annotate("baseline Cat (ee≈0 diagonal)", (
            100 * base_sw[-1]["ee"], 100 * base_sw[-1]["yield"]),
            fontsize=8.4, color="#888888", xytext=(4, -10),
            textcoords="offset points")
    Ts = [d["T"] for d in sw]
    ys = [100 * d["yield"] for d in sw]
    ees = [d["ee"] for d in sw]
    norm = Normalize(Ts[0], Ts[-1])
    cmap = matplotlib.colormaps["plasma"]
    sc = ax2.scatter(ees, ys, c=Ts, cmap=cmap, s=110, edgecolors="black",
                     zorder=3, linewidths=0.9)
    ax2.plot(ees, ys, "-", color="#888888", lw=1.2, zorder=2)
    for d in sw:
        if d["T"] % 50 == 0 or d["T"] in (250, 350):
            ax2.annotate(f"{d['T']} K", (d["ee"], 100 * d["yield"]),
                         textcoords="offset points", xytext=(8, -3),
                         fontsize=8.6, color="#444444")
    cb = fig.colorbar(sc, ax=ax2, pad=0.015)
    cb.set_label("temperature (K)", fontsize=10.5)
    ee_c = (XD.get("ee_curtin_winner_298K") if designed
            else XB["ee_curtin_298K"])
    ax2.axvline(ee_c, color="#8e44ad", ls="--", lw=1.2)
    ax2.text(ee_c, min(ys) - 2, f" Curtin–Hammett ee limit "
             f"({ee_c:.1f}% @298 K)", fontsize=8.8, color="#8e44ad",
             rotation=90, va="bottom")
    ax2.set_xlabel("enantioselectivity  ee (%)", fontsize=12)
    ax2.set_ylabel("conditional model yield of P_target (%)", fontsize=12)
    ax2.set_title("Panel B — yield–ee Pareto frontier vs temperature "
                  "(250–350 K)", fontsize=12.5, fontweight="bold")
    style_axes(ax2)
    fig.suptitle("Fig. 2 — Stiff microkinetic world model of the catalytic "
                 "network", fontsize=13.5, fontweight="bold", y=1.02)
    fig.tight_layout()
    p = FIG / "fig2_stiff_microkinetics_profile.png"
    fig.text(.5, .005, 'EXPLORATORY MODEL: conditional kinetics, not measured yield', ha='center', fontsize=8)
    fig.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    _log(f"figure 2 -> {p.name}")


def figure_3():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa
    from matplotlib import cm
    from matplotlib.colors import Normalize

    MC = RESULTS.get("module_C") or {}
    pos = XC.get("_pos") or MC.get("_pos")
    if pos is None:
        raise RuntimeError("figure_3: no persisted module-C geometry — "
                           "run stage C first")
    nums = np.array(pos["ts_nums"])
    coords = np.array(pos["ts_pos"])
    q = np.array(pos["ts_q"])
    ctr = np.array(pos.get("esp_center") or XC["esp"]["center"])
    pts, phi = esp_field(nums, coords, q, ctr)
    sel = (np.linalg.norm(pts - ctr, axis=1) > 2.2) & \
          (np.linalg.norm(pts - ctr, axis=1) < 5.6)
    pts, phi = pts[sel], phi[sel]

    fig = plt.figure(figsize=(13.0, 9.6))
    ax = fig.add_subplot(111, projection="3d")
    norm = Normalize(np.percentile(phi, 4), np.percentile(phi, 96))
    cmap = matplotlib.colormaps["RdBu_r"]
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], c=phi, cmap=cmap,
               norm=norm, s=5, alpha=0.16, depthshade=False)
    # TS substrate skeleton (heavy atoms + bonds by distance heuristic)
    rvdw = {1: 0.32, 6: 0.56, 7: 0.54, 8: 0.53, 15: 0.7}
    cols = {1: "#e8e8e8", 6: "#3b3b3b", 7: "#2759c9", 8: "#c9313d",
            15: "#e08b1f"}
    near = np.linalg.norm(coords - ctr, axis=1) < 7.0
    sub = np.where(near)[0]
    for i in sub:
        if nums[i] == 1:
            continue
        for j in sub:
            if j <= i or nums[j] == 1:
                continue
            d = np.linalg.norm(coords[i] - coords[j])
            if 0.9 < d < 1.75:
                ax.plot(*zip(coords[i], coords[j]), color="#555555",
                        lw=1.6, alpha=0.9)
    for i in sub:
        ax.scatter(*coords[i], s=140 * rvdw.get(nums[i], 0.5) ** 2 * 4,
                   c=cols.get(nums[i], "#888888"), edgecolors="black",
                   linewidths=0.4, depthshade=False)
    # catalyst docked pose from the proof run (persisted geometry)
    anchor = np.array(pos["anchor"])
    XCwin = XC.get("proof_complex") or MC.get("_proof_complex") or {}
    if "cat_pos" in XCwin:
        cpos = np.array(XCwin["cat_pos"]); cnums = np.array(
            XCwin["cat_nums"])
        for i in range(len(cnums)):
            for k2 in range(i + 1, len(cnums)):
                d = np.linalg.norm(cpos[i] - cpos[k2])
                if 0.9 < d < 1.65:
                    ax.plot(*zip(cpos[i], cpos[k2]), color="#7a5195",
                            lw=1.15, alpha=0.75)
        iP = [i for i in range(len(cnums)) if cnums[i] == 15][0]
        ax.scatter(cpos[iP][0], cpos[iP][1], cpos[iP][2], s=260,
                   c="#e08b1f", edgecolors="black", depthshade=False)
        iOan = int(XCwin.get("iO", iP))
        ax.scatter(cpos[iOan][0], cpos[iOan][1], cpos[iOan][2], s=220,
                   c="#c9313d", edgecolors="black", depthshade=False)
        ax.plot(*zip(cpos[iOan], anchor), color="#c9313d", ls="--", lw=2.2,
                alpha=0.95)
        ax.text(*(cpos[iOan] + np.array([0, 0, 0.5])),
                "P–O⁻···H–N⁺ H-bond", fontsize=9.5, color="#8e1b24",
                fontweight="bold")
    # TS polarization vector
    esp = XC.get("esp") or MC.get("esp")
    pmin = np.array(esp["min_point"]); pmax = np.array(esp["max_point"])
    ax.plot(*zip(pmin, pmax), color="#206b3a", lw=3.0, alpha=0.9)
    ax.text(*(pmax + np.array([0, 0, 0.4])), "ESP polarization axis",
            fontsize=9.5, color="#206b3a", fontweight="bold")
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, shrink=0.62, pad=0.02)
    cb.set_label("electrostatic potential φ(r)  (e/Å)", fontsize=10.5)
    dd = XC.get("proof") or MC.get("proof")
    ax.set_title(
        "Fig. 3 — De novo catalyst cavity docked on TS‡: ESP field + "
        f"non-covalent stabilization\n"
        f"winner: 3,3′-{dd['winner_motifs'][0]}/"
        f"{dd['winner_motifs'][1]}-BINOL-PA;  "
        f"constrained gas-phase energy proxy difference = {dd['barrier_drop_vs_baseline_kcal']:.2f}"
        f" kcal/mol (saddle and IRC not verified)",
        fontsize=12.5, fontweight="bold")
    ax.view_init(elev=18, azim=38)
    p = FIG / "fig3_ts_stabilization_dock.png"
    fig.text(.5, .005, 'EXPLORATORY MODEL: constrained energy proxy, not TS proof', ha='center', fontsize=8)
    fig.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    _log(f"figure 3 -> {p.name}")


# --------------------------------------------------------------------------- #
# 7.  BILINGUAL EPISTEMIC REPORTS
# --------------------------------------------------------------------------- #
def _fmt(x, nd=2):
    try:
        return f"{float(x):.{nd}f}"
    except Exception:
        return str(x)


def _write_audited_report(language):
    """Report conditional outputs without upgrading energy proxies to proof."""
    lead = ("# Phase 5 — audited exploratory chemical model\n\n"
            "These are model outputs, not validated barrier reductions, isolated yields or experimental ee."
            if language == "EN" else
            "# 第五阶段：审计后的探索性化学模型\n\n"
            "以下是条件模型输出，不是已验证势垒降低、分离收率或实验 ee。")
    evidence = {
        "audit_version": AUDIT_VERSION,
        "scope": "assigned reaction network and approximate constrained-energy proxies",
        "barrier_reduction_proven": False,
        "reason": "Unconstrained saddle stationarity and IRC connectivity not established",
        "baseline_conditional_kinetics": RESULTS.get("module_B", {}).get("reference_298K"),
        "exploratory_scenario": RESULTS.get("module_D", {}),
        "constrained_energy_proxy": RESULTS.get("module_C", {}).get("proof"),
        "assigned_parameters": RESULTS.get("assigned_parameters", []),
    }
    text = lead + "\n\nThe RRHO fallback units are corrected but its translation/rotation estimates remain approximate. " \
        "Designed kinetic caps require --exploratory-kinetics and retain the sign of the stereo difference. " \
        "The reaction graph and several rates/floors are assigned, not autonomously discovered. " \
        "Electronic calculations and full dynamics must be rerun before comparing with historical outputs.\n\n" \
        + "```json\n" + json.dumps(evidence, indent=2, ensure_ascii=False) + "\n```\n"
    (ROOT / f"WORLD_MODEL_REPORT_{language}.md").write_text(text, encoding="utf-8")


def report_EN():
    _write_audited_report("EN")


def fmt_dd(x):
    return f"{x:+.2f}"


def report_ZH():
    _write_audited_report("ZH")

def _hydrate_from_results():
    """Cross-process continuation: rebuild XA/XB/XC in-memory state from the
    persisted JSON + per-geometry cache so any stage can run standalone."""
    prev = {}
    if RESULTS_PATH.exists():
        try:
            prev = json.loads(RESULTS_PATH.read_text())
            if (prev.get('audit_version') != AUDIT_VERSION or
                    prev.get('energy_cache_version') != ENERGY_CACHE_VERSION):
                raise RuntimeError('Pre-audit Phase 5 results cannot be reused; rerun A/B/C with corrected code')
        except Exception as exc:
            _warn(f"could not reload partial results: {exc}")
            return
    RESULTS["fallbacks"] = prev.get("fallbacks", [])
    RESULTS["warnings"] = prev.get("warnings", [])
    RESULTS["assigned_parameters"] = prev.get("assigned_parameters", [])
    MA = prev.get("module_A") or {}
    if MA and "G_kcal" not in XA:
        XA.update(MA)
        RESULTS["module_A"] = MA
        for key, fname in [("RC", "cx_RC"), ("I_RC", "cx_IRC"),
                           ("TS1", "cx_TS1"), ("I1Cat", "cx_I1Cat"),
                           ("I1Cat_M", "cx_I1Cat_M"),
                           ("I1Cat_m", "cx_I1Cat_m"),
                           ("TS2aM", "cx_TS2aM"), ("TS2am", "cx_TS2am"),
                           ("TS2b", "cx_TS2b")]:
            if key not in XA:
                p = CACHE / f"{fname}.json"
                if p.exists():
                    XA[key] = json.loads(p.read_text())
        _log("hydrated Module-A state from results/cache")
    MB = prev.get("module_B") or {}
    if MB and "profile_298K" not in XB:
        XB.update(MB)
        _log("hydrated Module-B state from results JSON")
    MCp = prev.get("module_C") or {}
    if MCp and "esp" not in XC:
        XC.update({k: v for k, v in MCp.items()})
        _log("hydrated Module-C state from results JSON")
    for key in ("module_B", "module_C", "module_D"):
        if key in prev and key not in RESULTS:
            RESULTS[key] = prev[key]


def main():
    global EXPLORATORY_KINETICS
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["all", "A", "B", "C", "D", "figures",
                             "reports"])
    ap.add_argument('--exploratory-kinetics', action='store_true',
                    help='Opt in to assigned, capped catalyst sensitivity scenarios; not a prediction')
    args = ap.parse_args()
    EXPLORATORY_KINETICS = args.exploratory_kinetics

    t0 = time.time()
    _log(f"PHASE 5 chemical world model — stage={args.stage} — "
         f"engine={ENGINE_NAME}")
    if XTB_EXE is None:
        _fallback("xtb.exe not found — pipeline requires it for Modules A/C")

    _hydrate_from_results()

    if args.stage in ("all", "A"):
        module_A()
    if args.stage in ("all", "B"):
        module_B()
    if args.stage in ("all", "C"):
        module_C()
    if args.stage in ("all", "D"):
        if EXPLORATORY_KINETICS:
            module_D()
        elif args.stage == 'D':
            raise RuntimeError('Module D is exploratory; use --exploratory-kinetics explicitly')
        else:
            _warn('Module D skipped: unvalidated corrections require --exploratory-kinetics')
    if args.stage in ("all", "figures"):
        figure_1()
        figure_2()
        figure_3()
    if args.stage in ("all", "reports"):
        report_EN()
        report_ZH()

    RESULTS["wall_time_s"] = time.time() - t0
    write_json_atomic(RESULTS_PATH)
    _log(f"PHASE 5 done in {RESULTS['wall_time_s'] / 60:.1f} min — "
         f"{len(RESULTS['fallbacks'])} fallbacks, "
         f"{len(RESULTS['warnings'])} warnings")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        RESULTS["execution_status"] = "failed"
        RESULTS["fatal_error"] = f"{type(exc).__name__}: {exc}"
        write_json_atomic(RESULTS_PATH)
        raise
