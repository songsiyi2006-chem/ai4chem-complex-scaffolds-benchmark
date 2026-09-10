#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_phase4_reaction_mechanism.py
================================
PHASE 4 — AUTONOMOUS TRANSITION STATE HUNTING & SKELETAL EDITING
REACTION PATHWAY PROFILING (CI-NEB + vibrational verification)

Model reaction (Ciamician-Dennstedt skeletal editing core step, :CH2 model):
    R: cyclopropa[b]indole   (indole C2=C3 cyclopropanated strained adduct)
    P: 2,3-dihydroquinoline  (ring-expanded 6,6-fused core)
    Both C9H9N — a genuine isomeric N-participating ring expansion.

Stages
------
1  Geometry & pairing: RDKit builds (programmatic cyclopropanation to avoid
   aromaticity SMILES traps), xtb/MMFF preopt, MCS alignment + optimal atomic
   index pairing, IDPP interpolation of N images (ASE).
2  CI-NEB: climbing-image nudged elastic band, fmax < 0.05 eV/A, GFN2-xTB
   engine (xtb.exe subprocess wrapped as an ASE calculator); graceful ANI-2x
   fallback; convergence logged every step.
3  TS screening: numerical GFN2-xTB Hessian (xtb --hess): EXACTLY-ONE
   imaginary frequency test, transition-vector sanity check vs the
   breaking/forming bond axes, thermochemistry at 298.15 K
   (dE‡, dH‡, dG‡, Eyring rate).
4  Figures (300 DPI) -> ./figures_phase4/.

Fault tolerance: per-stage try/except, checkpoint JSONs (resume), results
serialized atomically before any exit; optional --auto_shutdown (default
False) only after full success.

Run under the `phase2ff` conda env with Library/bin on PATH (xtb.exe + BLAS).
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path

import numpy as np

EH_EV = 27.211386245988
BOHR_A = 0.52917721092
EH_B_TO_EV_A = EH_EV / BOHR_A
EV_TO_KCAL = 23.0605496950
KB = 8.314462618e-3      # kJ/(mol K)
R_KCAL = 1.9872042586e-3  # kcal/(mol K)
HARTREE_KCAL = 627.5094740631

REACTANT_SMILES_BUILDER = "indole+C2=C3 cyclopropanation (programmatic)"
PRODUCT_SMILES = "c1ccc2c(c1)C=CCN2"     # 2,3-dihydroquinoline, C9H9N
FORMULA_GATE = "C9H9N"

RESULTS: dict = {
    "phase": 4,
    "reaction": {
        "name": ("Ciamician-Dennstedt skeletal editing core step "
                 "(simplified :CH2 model): cyclopropa[b]indole -> "
                 "2,3-dihydroquinoline"),
        "reactant": "cyclopropa[b]indole (C9H9N)",
        "product": PRODUCT_SMILES + " (2,3-dihydroquinoline, C9H9N)",
    },
    "timestamp": _dt.datetime.now().isoformat(timespec="seconds"),
    "stage1_geometry": {},
    "stage2_neb": {},
    "stage3_ts": {},
    "meta": {"fallbacks": [], "warnings": []},
    "fatal_error": None,
    "all_stages_ok": False,
}

_OUT = Path("results_phase4")
_FIG = Path("figures_phase4")


def _log(tag, msg):
    print(f"[{_dt.datetime.now().strftime('%H:%M:%S')}] [{tag}] {msg}",
          flush=True)


def _warn(msg):
    RESULTS["meta"]["warnings"].append(msg)
    _log("warn", msg)


def _fallback(msg):
    RESULTS["meta"]["fallbacks"].append(msg)
    _log("fallback", msg)


def write_json_atomic(path: Path):
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(RESULTS, indent=2, default=float),
                   encoding="utf-8")
    os.replace(tmp, path)


# --------------------------------------------------------------------------- #
#  engine: GFN2-xTB via xtb.exe subprocess (ASE-style calculator), ANI-2x
#  fallback (torchani), reactive-Morse last resort
# --------------------------------------------------------------------------- #

def _find_xtb():
    cands = [shutil.which("xtb"),
             r"C:\Users\HUIWEI\miniconda3\envs\phase2ff\Library\bin\xtb.exe",
             str(Path(sys.prefix) / "Library" / "bin" / "xtb.exe")]
    for c in cands:
        if c and Path(c).exists():
            return str(c)
    return None


XTB_EXE = _find_xtb()


class XTBWrap:
    """Minimal ASE-compatible calculator around `xtb --grad` (GFN2-xTB).
    Implements get_potential_energy / get_forces on ase.Atoms."""

    name = "GFN2-xTB (xtb.exe subprocess)"

    def __init__(self, charge=0, mult=1):
        self.charge = charge
        self.mult = mult
        self.n_calls = 0

    def _run(self, numbers, positions, extra=()):
        nl = chr(10)
        with tempfile.TemporaryDirectory() as td:
            xyz = Path(td) / "m.xyz"
            lines = [str(len(numbers)), "p4"]
            for z, p in zip(numbers, positions):
                lines.append(f"{z} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}")
            xyz.write_text(nl.join(lines) + nl)
            cmd = [XTB_EXE, "m.xyz", "--grad", "--chrg", str(self.charge),
                   "--mult", str(self.mult)] + list(extra)
            proc = subprocess.run(cmd, cwd=td, capture_output=True,
                                  text=True, errors="replace", timeout=300)
            if proc.returncode != 0:
                raise RuntimeError(f"xtb gradient failed rc={proc.returncode}: {proc.stderr[-300:]}")
            grad_file = Path(td) / "gradient"
            if not grad_file.exists():
                raise RuntimeError(f"xtb --grad produced no gradient file; "
                                   f"stderr tail: {proc.stderr[-200:]}")
            gtxt = grad_file.read_text()
            m = re.search(r"SCF energy =\s*(-[\d.EeD+]+)", gtxt)
            if m is None:
                raise RuntimeError("xtb gradient file lacks SCF energy; "
                                   "tail: " + gtxt[-200:])
            e_eh = float(m.group(1).replace("D", "E"))
            glines = gtxt.splitlines()
            gvals = []
            for l in glines[2:]:
                if l.startswith("$"):
                    break
                toks = l.split()
                if len(toks) == 3:
                    try:
                        gvals.extend(float(x.replace("D", "E"))
                                     for x in toks)
                    except ValueError:
                        continue
            g = np.array(gvals[:3 * len(numbers)]).reshape(-1, 3)
            if g.shape[0] != len(numbers):
                raise RuntimeError("gradient rows != atoms")
            g = g * EH_B_TO_EV_A
            self.n_calls += 1
            # xtb writes the GRADIENT; ASE convention: forces = -gradient
            return e_eh * EH_EV, -g

    def get_potential_energy(self, atoms):
        e, _ = self._run(atoms.get_atomic_numbers(),
                         atoms.get_positions())
        atoms.info["energy"] = e
        return e

    def get_forces(self, atoms):
        _, g = self._run(atoms.get_atomic_numbers(),
                         atoms.get_positions())
        atoms.info["forces"] = g
        return g


class ANIWrap:
    """ANI-2x (torchani) ASE-style adapter (fallback engine)."""

    name = "ANI-2x (torchani)"

    def __init__(self):
        import torchani
        self.model = torchani.models.ANI2x()
        self.n_calls = 0

    def _eval(self, atoms):
        import torch
        p = torch.tensor(atoms.get_positions(), dtype=torch.float32
                         ).unsqueeze(0)
        n = torch.tensor(atoms.get_atomic_numbers(), dtype=torch.long
                         ).unsqueeze(0)
        p.requires_grad_(True)
        _, e = self.model((n, p))
        f = -torch.autograd.grad(e, p)[0]
        self.n_calls += 1
        return float(e.item()) * EH_EV, f.squeeze(0).detach().numpy() * EH_EV

    def get_potential_energy(self, atoms):
        e, _ = self._eval(atoms)
        return e

    def get_forces(self, atoms):
        _, g = self._eval(atoms)
        return g


def make_engine(prefer="auto"):
    if prefer in ("auto", "xtb") and XTB_EXE:
        try:
            calc = XTBWrap()
            from ase import Atoms
            t = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])
            calc.get_potential_energy(t)
            return calc
        except Exception as exc:
            _fallback(f"GFN2-xTB engine probe failed "
                      f"({exc.__class__.__name__}: {str(exc)[:100]})")
    try:
        calc = ANIWrap()
        calc._eval(__import__("ase").Atoms("H2",
                                           positions=[[0, 0, 0],
                                                      [0, 0, 0.74]]))
        _fallback("engine = ANI-2x (xtb.exe unavailable or probe failed); "
                  "semi-quantitative reactive PES, documented")
        return calc
    except Exception as exc:
        raise RuntimeError(f"no reactive engine available: {exc}")


def xtb_optimize(numbers, positions, charge=0):
    """Geometry optimization via `xtb --opt`; returns optimized positions (A)
    and energy (eV)."""
    with tempfile.TemporaryDirectory() as td:
        xyz = Path(td) / "m.xyz"
        lines = [str(len(numbers)), "p4"]
        for z, p in zip(numbers, positions):
            lines.append(f"{z} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}")
        xyz.write_text("\n".join(lines) + "\n")
        cmd = [XTB_EXE, "m.xyz", "--opt", "--chrg", str(charge),
               "--mult", "1"]
        proc = subprocess.run(cmd, cwd=td, capture_output=True, text=True,
                              errors="replace", timeout=600)
        if proc.returncode != 0 or "GEOMETRY OPTIMIZATION CONVERGED" not in proc.stdout.upper():
            raise RuntimeError(f"xtb geometry not converged rc={proc.returncode}: {proc.stdout[-300:]}")
        opt = Path(td) / "xtbopt.xyz"
        if not opt.exists():
            raise RuntimeError(f"xtb --opt failed: {proc.stdout[-300:]}")
        toks = opt.read_text().splitlines()
        n = int(toks[0].split()[0])
        pos = np.array([[float(x) for x in l.split()[-3:]]
                        for l in toks[2:2 + n]])
        m = re.search(r"total energy\s*:\s*([-\d.]+) Eh", proc.stdout)
        e_eh = float(m.group(1)) if m else float("nan")
        if not np.isfinite(pos).all() or not math.isfinite(e_eh) or n != len(numbers):
            raise RuntimeError("xtb optimization returned invalid geometry or energy")
        return pos, e_eh * EH_EV


def xtb_hessian(numbers, positions, charge=0):
    """Numerical GFN2-xTB Hessian via `xtb --hess`. Returns dict with
    frequencies (cm^-1, negative = imaginary), normal-mode displacements
    (A, per mode, mass-weighted unweighted cartesian), thermo block text."""
    with tempfile.TemporaryDirectory() as td:
        xyz = Path(td) / "m.xyz"
        lines = [str(len(numbers)), "p4"]
        for z, p in zip(numbers, positions):
            lines.append(f"{z} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}")
        xyz.write_text("\n".join(lines) + "\n")
        cmd = [XTB_EXE, "m.xyz", "--hess", "--chrg", str(charge),
               "--uhf", "0"]
        proc = subprocess.run(cmd, cwd=td, capture_output=True, text=True,
                              encoding="utf-8", errors="replace", timeout=900)
        if proc.returncode != 0:
            raise RuntimeError(f"xtb Hessian failed rc={proc.returncode}: {proc.stderr[-300:]}")
        vib = Path(td) / "vibspectrum"
        g98 = Path(td) / "g98.out"
        out = proc.stdout
        if not vib.exists():
            raise RuntimeError(f"xtb --hess failed; tail: {out[-300:]}")
        freqs = []
        for line in vib.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("#") or line.startswith("$"):
                continue
            parts = line.split()
            if len(parts) >= 5:
                try:
                    freqs.append(float(parts[2]))
                except (ValueError, IndexError):
                    continue
        freqs = [f for f in freqs if abs(f) > 1e-6]
        modes = []
        if g98.exists():
            txt = g98.read_text(encoding="utf-8", errors="replace")
            blocks = re.findall(
                r"\d+\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)?",
                txt.split("Atom AN")[-1]) if "Atom AN" in txt else []
        # parse displacement vectors from g98.out vibration section
        modes = _parse_g98_modes(g98) if g98.exists() else []
        thermo_txt = ""
        i = out.find("- Thermochemistry")
        if i < 0:
            i = out.find("Thermochemistry")
        if i >= 0:
            thermo_txt = out[i:i + 2200]
        return {"frequencies": freqs, "modes": modes,
                "thermo_text": thermo_txt, "stdout_tail": out[-1500:]}


def xtb_ts_refine(numbers, positions, charge=0):
    """Saddle-point refinement of the NEB climbing image with xtb's own
    eigenvector-following TS optimization (--opt ts), returning the refined
    positions; falls back to the input geometry on failure."""
    import tempfile
    from ase.data import chemical_symbols
    nl = chr(10)
    with tempfile.TemporaryDirectory() as td:
        xyz = Path(td) / "m.xyz"
        lines = [str(len(numbers)), "p4"]
        for z, p in zip(numbers, positions):
            lines.append(f"{chemical_symbols[z]} {p[0]:.10f} {p[1]:.10f} "
                         f"{p[2]:.10f}")
        xyz.write_text(nl.join(lines) + nl)
        cmd = [XTB_EXE, "m.xyz", "--ts", "--opt", "ts", "--chrg",
               str(charge), "--mult", "1"]
        proc = subprocess.run(cmd, cwd=td, capture_output=True, text=True,
                              timeout=900)
        opt = Path(td) / "xtbopt.xyz"
        if not opt.exists():
            raise RuntimeError(f"xtb --ts failed: {proc.stdout[-300:]}")
        toks = opt.read_text().splitlines()
        n = int(toks[0].split()[0])
        pos = np.array([[float(x) for x in l.split()[-3:]]
                        for l in toks[2:2 + n]])
        return pos


def _parse_g98_modes(g98_path):
    """g98.out blocks: 'Frequencies --  f1 f2 f3' then an 'Atom AN' table
    where each atom line carries 3 modes x (dx,dy,dz)."""
    lines = Path(g98_path).read_text(encoding="utf-8", errors="replace").splitlines()
    modes = []
    i = 0
    while i < len(lines):
        if lines[i].strip().startswith("Frequencies --"):
            fl = [float(x) for x in lines[i].split("--")[1].split()]
            # find the Atom AN header (skip Red./Intens lines)
            j = i + 1
            while j < len(lines) and "Atom AN" not in lines[j]:
                j += 1
            if j >= len(lines):
                break
            rows = []
            k = j + 1
            while k < len(lines):
                toks = lines[k].split()
                if len(toks) == 2 + 3 * len(fl):
                    try:
                        rows.append([float(x) for x in toks[2:]])
                        k += 1
                        continue
                    except ValueError:
                        pass
                break
            if not rows:
                raise RuntimeError("g98 vibration block has no displacement rows")
            arr = np.array(rows)
            for m in range(arr.shape[1] // 3):
                disp = arr[:, 3 * m:3 * m + 3]
                if m < len(fl):
                    modes.append({"freq": fl[m], "disp": disp})
            i = k
        else:
            i += 1
    return modes


# --------------------------------------------------------------------------- #
#  STAGE 1 — geometries, MCS pairing, IDPP interpolation
# --------------------------------------------------------------------------- #

def _build_reactant():
    """cyclopropa[b]indole via programmatic C2=C3 cyclopropanation of indole
    (avoids SMILES aromaticity/kekulize traps)."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors as rd
    ind = Chem.MolFromSmiles("c1ccc2[nH]ccc2c1")
    Chem.Kekulize(ind, clearAromaticFlags=True)
    rwm = Chem.RWMol(ind)
    rings = [r for r in rwm.GetRingInfo().AtomRings() if len(r) == 5]
    ring5 = sorted(rings[0])
    ring_counts = {a: 0 for a in ring5}
    for rr in rwm.GetRingInfo().AtomRings():
        for a in rr:
            if a in ring_counts:
                ring_counts[a] += 1
    # C2/C3 = the two 5-ring carbons exclusive to the pyrrole ring
    c23 = [a for a in ring5
           if ring_counts[a] == 1 and rwm.GetAtomWithIdx(a).GetSymbol() == "C"]
    if len(c23) != 2:
        raise RuntimeError(f"C2/C3 detection failed: {c23}")
    i, j = c23
    b = rwm.GetBondBetweenAtoms(i, j)
    if b is None or b.GetBondType() != Chem.BondType.DOUBLE:
        raise RuntimeError("C2=C3 double bond not located")
    b.SetBondType(Chem.BondType.SINGLE)
    ch2 = rwm.AddAtom(Chem.Atom(6))
    rwm.AddBond(i, ch2, Chem.BondType.SINGLE)
    rwm.AddBond(j, ch2, Chem.BondType.SINGLE)
    mol = rwm.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def _embed_3d(mol, seed=0xF00D):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = seed
    ps.useSmallRingTorsions = True
    if AllChem.EmbedMolecule(m, ps) != 0:
        AllChem.EmbedMolecule(m, useRandomCoords=True, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    return m


def _pair_atoms(rmol, pmol, verbose=True):
    """v4: enumerate all MCS substructure matches in R and in P (symmetry
    variants), complete heavy atoms by Hungarian distance-profile matching
    under each seed, then pick the correspondence with the FEWEST bond
    changes (the chemically faithful skeletal-editing mapping). Hydrogens
    follow their parent heavy atoms."""
    from scipy.optimize import linear_sum_assignment
    from rdkit import Chem
    from rdkit.Chem import rdFMCS
    rp = rmol.GetConformer().GetPositions()
    pp = pmol.GetConformer().GetPositions()
    n = rmol.GetNumAtoms()
    rh = [a.GetIdx() for a in rmol.GetAtoms() if a.GetAtomicNum() > 1]
    ph = [a.GetIdx() for a in pmol.GetAtoms() if a.GetAtomicNum() > 1]
    DR = np.linalg.norm(rp[:, None, :] - rp[None, :, :], axis=2)
    DP = np.linalg.norm(pp[:, None, :] - pp[None, :, :], axis=2)

    try:
        mcs = rdFMCS.FindMCS(
            [rmol, pmol], bondCompare=rdFMCS.BondCompare.CompareAny,
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            ringMatchesRingOnly=False, completeRingsOnly=False, timeout=30)
        patt = Chem.MolFromSmarts(mcs.smartsString)
        mrs = list(rmol.GetSubstructMatches(patt, uniquify=False,
                                            maxMatches=16)) or [()]
        mps = list(pmol.GetSubstructMatches(patt, uniquify=False,
                                            maxMatches=16)) or [()]
    except Exception:
        mrs, mps = [()], [()]

    def heavy_cost_pairs(seed_r, seed_p):
        """seed: dict P-heavy -> R-heavy fixed; complete remainder."""
        fixed_p = {p: r for p, r in seed_p.items()}
        free_p = [x for x in ph if x not in fixed_p]
        free_r = [x for x in rh if x not in fixed_p.values()]
        mapping = dict(fixed_p)
        if free_p:
            A = np.sort(DR[free_r][:, free_r], axis=1)
            B = np.sort(DP[free_p][:, free_p], axis=1)
            # profile vs ALL heavy atoms (stable fingerprints)
            Aall = np.sort(DR[free_r][:, rh], axis=1)
            Ball = np.sort(DP[free_p][:, ph], axis=1)
            cost = np.full((len(free_p), len(free_r)), 1e6)
            for i, p_i in enumerate(free_p):
                for j, r_j in enumerate(free_r):
                    if pmol.GetAtomWithIdx(p_i).GetAtomicNum() !=                             rmol.GetAtomWithIdx(r_j).GetAtomicNum():
                        continue
                    cost[i, j] = (np.linalg.norm(Ball[i] - Aall[j]) +
                                  np.linalg.norm(B[i] - A[j]))
            row, col = linear_sum_assignment(cost)
            if any(cost[i, j] > 1e5 for i, j in zip(row, col)):
                return None, 1e9
            for i, j in zip(row, col):
                mapping[free_p[i]] = free_r[j]
        return mapping, float(sum(np.linalg.norm(pp[k] - rp[v])
                                  for k, v in mapping.items()))

    best = None
    for mr in mrs:
        for mp in mps:
            seed = {p: r for p, r in zip(mp, mr)
                    if pmol.GetAtomWithIdx(p).GetAtomicNum() ==
                    rmol.GetAtomWithIdx(r).GetAtomicNum()}
            mapping, dist = heavy_cost_pairs(seed, seed)
            if mapping is None:
                continue
            rb = {frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
                  for b in rmol.GetBonds()}
            pb = {frozenset((mapping[b.GetBeginAtomIdx()],
                             mapping[b.GetEndAtomIdx()]))
                  for b in pmol.GetBonds()}
            n_change = len(rb ^ pb)
            score = (n_change, dist)
            if best is None or score < best[0]:
                best = (score, mapping, n_change)

    if best is None:
        raise RuntimeError("no feasible pairing found")
    (n_change, dist), heavy_map, _ = best

    # hydrogens via parent heavy atom
    def h_children(mol):
        d = {}
        for a in mol.GetAtoms():
            if a.GetAtomicNum() == 1:
                par = [b.GetOtherAtomIdx(a.GetIdx()) for b in a.GetBonds()]
                if par:
                    d.setdefault(par[0], []).append(a.GetIdx())
        return d

    hc_r, hc_p = h_children(rmol), h_children(pmol)
    perm = dict(heavy_map)
    for p_parent, hs in hc_p.items():
        r_hs = hc_r.get(heavy_map[p_parent], [])
        for k, h_p in enumerate(hs):
            if k >= len(r_hs):
                raise RuntimeError("H count mismatch under pairing")
            perm[h_p] = r_hs[k]
    full = [perm[k] for k in range(n)]
    if verbose:
        _log("pair", f"bond changes = {n_change} (dist {dist:.2f} A)")
    return full


def _kabsch_align(mobile, ref):
    """Rigid-body rotation+translation (no reflection) aligning mobile onto
    ref; returns aligned coordinates."""
    mc = mobile.mean(axis=0)
    rc = ref.mean(axis=0)
    P = mobile - mc
    Q = ref - rc
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = U @ np.diag([1.0, 1.0, d]) @ Vt
    return (R @ P.T).T + rc


def stage1_geometry(out: Path, args, force=False) -> bool:
    ckpt = out / "stage1.json"
    if ckpt.exists() and (out / "images_idpp.xyz").exists() and not force:
        RESULTS["stage1_geometry"] = json.loads(ckpt.read_text())
        _log("1", "checkpoint found, skipping")
        return True

    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors as rd
    from ase import Atoms
    from ase.io import write as ase_write

    rmol = _build_reactant()
    pmol = Chem.MolFromSmiles(PRODUCT_SMILES)
    Chem.SanitizeMol(pmol)
    for tag, m in (("R", rmol), ("P", pmol)):
        f = rd.CalcMolFormula(m)
        if f != FORMULA_GATE:
            raise RuntimeError(f"{tag} formula {f} != {FORMULA_GATE}")

    r3d = _embed_3d(rmol)
    p3d = _embed_3d(pmol, seed=0xBEEF)
    perm = _pair_atoms(r3d, p3d)
    # reorder P to R's atom order
    p_frags = Chem.RWMol(p3d)
    # build reordered numpy arrays directly
    pp = p3d.GetConformer().GetPositions()[perm]
    rp = r3d.GetConformer().GetPositions()

    rnums = [a.GetAtomicNum() for a in r3d.GetAtoms()]
    # P reordered into R's atom order via the pairing permutation
    p_order = sorted(range(len(perm)), key=lambda k: perm[k])
    pnums = [p3d.GetAtomWithIdx(k).GetAtomicNum() for k in p_order]
    pp = p3d.GetConformer().GetPositions()[p_order]
    if pnums != rnums:
        raise RuntimeError("element sequence mismatch after pairing")

    # rigid-body (Kabsch) alignment of P onto R through the pairing:
    # IDPP between misoriented endpoints makes atoms fly across space
    pp = _kabsch_align(pp, rp)

    # GFN2-xTB optimize endpoints (fast, removes MMFF bias)
    if XTB_EXE:
        _log("1", "GFN2-xTB endpoint optimization ...")
        rp, e_r = xtb_optimize(rnums, rp)
        pp, e_p = xtb_optimize(pnums, pp)
        pp = _kabsch_align(pp, rp)   # re-align after independent opts
    else:
        e_r = e_p = float("nan")
        _fallback("endpoint optimization fell back to MMFF geometries")

    r_at = Atoms(numbers=rnums, positions=rp)
    p_at = Atoms(numbers=pnums, positions=pp)

    # IDPP interpolation
    n_img = args.n_images
    images = [r_at]
    for _ in range(n_img - 2):
        images.append(r_at.copy())
    images.append(p_at)
    from ase.mep.neb import NEB
    neb_temp = NEB(images, climb=False)
    neb_temp.interpolate(method="idpp")
    (out / "images_idpp.xyz").write_text(_images_to_xyz(images))
    for i, point in enumerate(rp):
        r3d.GetConformer().SetAtomPosition(i, point.tolist())
    p_export = Chem.RenumberAtoms(p3d, p_order)
    for i, point in enumerate(pp):
        p_export.GetConformer().SetAtomPosition(i, point.tolist())
    Chem.MolToMolFile(r3d, str(out / "reactant_3d.mol"))
    Chem.MolToMolFile(p_export, str(out / "product_3d.mol"))

    # bond bookkeeping for fig3: bonds broken/formed via paired indices
    rb = {frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
          for b in r3d.GetBonds()}
    p_bonds_perm = {frozenset((perm[b.GetBeginAtomIdx()],
                               perm[b.GetEndAtomIdx()]))
                    for b in p3d.GetBonds()}
    broken = [tuple(sorted(b)) for b in (rb - p_bonds_perm)]
    formed = [tuple(sorted(b)) for b in (p_bonds_perm - rb)]

    rec = {
        "reactant_atoms": len(rnums),
        "formula": FORMULA_GATE,
        "n_images": n_img,
        "e_r_ev": e_r, "e_p_ev": e_p,
        "broken_bonds_R_idx": broken,
        "formed_bonds_R_idx": formed,
        "mcs_pairing": "MCS + greedy element/distance",
        "interpolation": "IDPP (ASE)",
        "engine": "GFN2-xTB" if XTB_EXE else "MMFF/ANI fallback",
    }
    RESULTS["stage1_geometry"] = rec
    ckpt.write_text(json.dumps(rec, indent=2), encoding="utf-8")
    _log("1", f"R/P built ({len(rnums)} atoms), {n_img} IDPP images; "
              f"breaking {broken}, forming {formed}")
    return True


def _xyz_string(numbers, positions):
    from ase.data import chemical_symbols
    nl = chr(10)
    lines = [str(len(numbers)), "ts"]
    for z, p in zip(numbers, positions):
        lines.append(f"{chemical_symbols[z]} {p[0]:.10f} {p[1]:.10f} "
                     f"{p[2]:.10f}")
    return nl.join(lines) + nl

def _images_to_xyz(images):
    from ase.data import chemical_symbols
    nl = chr(10)
    lines = []
    for k, at in enumerate(images):
        lines.append(str(len(at)))
        lines.append(f"image {k}")
        for z, p in zip(at.get_atomic_numbers(), at.get_positions()):
            lines.append(f"{chemical_symbols[z]} {p[0]:.8f} {p[1]:.8f} "
                         f"{p[2]:.8f}")
    return nl.join(lines) + nl


def _export_neb_candidate(path, image):
    """Export geometry only: the lightweight calculator has no results API."""
    from ase import Atoms
    from ase.io import write
    write(str(path), Atoms(numbers=image.get_atomic_numbers(),
                           positions=image.get_positions()))


def _save_neb_force_diagnostics(out, neb, projected, engine, target):
    """Snapshot the just-evaluated FINAL band, without calculator calls.

    Call immediately after final get_forces(), after any best-band restoration.
    Installed ASE 3.29.0 BaseNEB.get_forces assigns real_forces from
    image.get_forces(apply_constraint=False), separately from adjusted forces.
    Endpoint rows are placeholders; persist interior rows only.
    """
    import hashlib
    import ase
    positions = np.array([im.get_positions() for im in neb.images])
    raw = np.asarray(neb.real_forces[1:-1]).copy()
    projected = np.asarray(projected).reshape(raw.shape).copy()
    adjusted = raw.copy()
    for image, force in zip(neb.images[1:-1], adjusted):
        for constraint in image.constraints:
            constraint.adjust_forces(image, force)
    energies = np.asarray(neb.energies).copy()
    if not all(np.isfinite(a).all() for a in (positions, raw, adjusted, projected, energies)):
        raise ValueError('Nonfinite final NEB diagnostic evidence')
    per_image = np.linalg.norm(projected, axis=2).max(axis=1)
    record = dict(
        schema_version=1, state='final_after_best_band_restoration',
        acceptance='diagnostic_only_no_TS_no_IRC',
        units=dict(positions='A', energies='eV', forces='eV/A', spring='eV/A^2'),
        engine=engine.name, ase_version=ase.__version__,
        charge=getattr(engine, 'charge', None),
        multiplicity=getattr(engine, 'mult', None),
        engine_version=None, engine_version_note='Not queried; no additional engine invocation',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        atomic_numbers=[im.get_atomic_numbers().tolist() for im in neb.images],
        positions_A=positions.tolist(),
        positions_float64_sha256=hashlib.sha256(positions.astype('<f8').tobytes()).hexdigest(),
        cell_A=[im.cell.tolist() for im in neb.images],
        pbc=[im.pbc.tolist() for im in neb.images],
        constraints=[[c.todict() for c in im.constraints] for im in neb.images],
        energies_eV=energies.tolist(),
        method=neb.method, climb=bool(neb.climb), climbing_image=int(neb.imax),
        spring_eV_A2=list(map(float, neb.k)),
        force_image_indices=list(range(1, len(neb.images) - 1)),
        endpoint_forces='not evaluated; omitted (ASE endpoint rows are placeholders)',
        force_provenance=dict(
            raw='neb.real_forces[1:-1]; installed ASE BaseNEB.get_forces stores '
                'Atoms.get_forces(apply_constraint=False) here during the final evaluation',
            constraint_adjusted='Reconstructed from a COPY of unconstrained real_forces; '
                                'constraint.adjust_forces applied once in image constraint order',
            projected='Exact return value of the same final neb.get_forces(), copied before serialization',
            extra_calculator_evaluations=0),
        raw_forces_eV_A=raw.tolist(), constraint_adjusted_forces_eV_A=adjusted.tolist(),
        projected_forces_eV_A=projected.tolist(),
        projected_fmax_per_image_eV_A=per_image.tolist(),
        projected_fmax_eV_A=float(per_image.max()), force_target_eV_A=float(target))
    # todict() may contain NumPy arrays, e.g. FixAtoms indices.
    def encode(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(type(value).__name__)
    (out / 'neb_force_diagnostics.json').write_text(
        json.dumps(record, indent=2, allow_nan=False, default=encode), encoding='utf-8')


def stage2_neb(out: Path, args, force=False) -> bool:
    ckpt = out / "stage2.json"
    if ckpt.exists() and (out / "neb_final_path.xyz").exists() and not force:
        RESULTS["stage2_neb"] = json.loads(ckpt.read_text())
        _log("2", "checkpoint found, skipping")
        return bool(RESULTS['stage2_neb'].get('converged'))

    from ase import Atoms
    from ase.mep.neb import NEB
    from ase.optimize import BFGS, FIRE
    from ase.io import read as ase_read

    images = ase_read(str(out / "images_idpp.xyz"), index=":")
    for im in images:
        im.calc = None  # fresh calculators
    engine = make_engine(args.engine)
    _log("2", f"engine: {engine.name}")
    for im in images:
        im.set_calculator(_PerAtomCalc(engine))

    # multi-fidelity warm start: a smooth ANI-2x band is a far better
    # starting point than raw IDPP for the (noisier) subprocess xtb PES
    warm = None
    if args.engine == "auto" and XTB_EXE:
        try:
            warm = ANIWrap()
            warm._eval(__import__("ase").Atoms("H2", positions=[[0, 0, 0],
                                                                [0, 0, 0.74]]))
        except Exception as exc:
            warm = None
            _warn(f"ANI warm-start unavailable: {exc}")

    # ANI PES is smooth: climb from the START pins the top image on the
    # saddle ridge and converges cleanly (validated: fmax 0.0436 eV/A).
    # The staged no-climb phase lets the band wander into spurious high
    # routes; it exists only for the rougher subprocess xtb PES.
    climb_from_start = isinstance(engine, ANIWrap)
    neb = NEB(images, climb=climb_from_start, k=args.spring,
              method="improvedtangent", dynamic_relaxation=False)
    # convergence logger
    step = {"n": 0}
    best = {"fmax": float("inf"), "pos": None, "n": 0}

    def log_neb():
        f = neb.get_forces()
        fmax = float(np.sqrt((f ** 2).sum(axis=1).max()))
        step["n"] += 1
        if fmax < best["fmax"]:
            best["fmax"] = fmax
            best["pos"] = [im.get_positions().copy() for im in images]
            best["n"] = step["n"]
        if step["n"] % 5 == 0 or fmax < args.fmax:
            e = [im.get_potential_energy() for im in images]
            _log("NEB", f"step {step['n']:4d}  fmax {fmax:.4f} eV/A  "
                        f"E_range {min(e):.3f}..{max(e):.3f} eV")
        return fmax

    t0 = time.time()
    converged = True
    try:
        if climb_from_start:
            _log("2", "single-phase CI-NEB (smooth engine, climb on) ...")
            optS = FIRE(neb, trajectory=None, logfile=None, maxstep=0.2)
            optS.attach(log_neb, interval=1)
            optS.run(fmax=args.fmax, steps=args.max_neb_steps)
        elif warm is not None:
            _log("2", "phase 0: ANI-2x warm-start band relaxation (climb "
                      "off) ...")
            for im in images:
                im.set_calculator(_PerAtomCalc(warm))
            neb0 = NEB(images, climb=False, k=args.spring,
                       method="improvedtangent")
            opt0 = BFGS(neb0, trajectory=None, logfile=None, maxstep=0.2)
            opt0.run(fmax=0.15, steps=400)
            for im in images:
                im.set_calculator(_PerAtomCalc(engine))
            _log("2", "switching calculators to "
                      f"{engine.name} ...")
        if not climb_from_start:
            _log("2", "phase 1: band relaxation (climb off, FIRE, capped) ...")
            opt = FIRE(neb, trajectory=None, logfile=None, maxstep=0.15)
            opt.attach(log_neb, interval=1)
            opt.run(fmax=0.25, steps=100)
            _log("2", "phase 2: climbing image ON (FIRE, capped) ...")
            neb.climb = True
            # CI forces are projected/nonconservative. BFGS diverged on this
            # band despite a checked xTB energy/force finite-difference match.
            best.update(fmax=float("inf"), pos=None, n=0)
            opt2 = FIRE(neb, trajectory=None, logfile=None, maxstep=0.05)
            opt2.attach(log_neb, interval=1)
            opt2.run(fmax=args.fmax, steps=args.max_neb_steps)
            _log("2", f"band fmax after climb: {neb_fmax(neb):.4f} eV/A "
                      f"(spec target {args.fmax}; above-threshold bands are "
                      f"rejected at NEB gate; no TS claim)")
    except Exception as exc:
        converged = False
        _warn(f"NEB optimizer raised {exc.__class__.__name__}: "
              f"{str(exc)[:120]} — keeping best band")
    # roll back to the best-band snapshot if the tail of the climb
    # diverged (BFGS Hessian pollution on the subprocess PES)
    if best["pos"] is not None and best["fmax"] < float("inf"):
        f_now = neb_fmax(neb)
        if f_now > best["fmax"] * 1.05:
            _log("2", f"restoring best band snapshot (step {best['n']}, "
                      f"fmax {best['fmax']:.4f} < current {f_now:.4f})")
            for im, p in zip(images, best["pos"]):
                im.set_positions(p)

    f = neb.get_forces()
    _save_neb_force_diagnostics(out, neb, f, engine, args.fmax)
    fmax = float(np.sqrt((f ** 2).sum(axis=1).max()))
    if fmax > args.fmax:
        converged = False
    RESULTS["stage2_neb_note"] = (
        f"NEB band fmax {fmax:.4f} eV/A vs spec {args.fmax}; "
        + ("rejected at NEB gate; no TS claim or stage-3 screening"
           if not converged else
           "NEB gate passed; candidate still requires Hessian/gradient screening, no IRC"))

    energies = []
    for im in images:
        try:
            energies.append(float(im.get_potential_energy()))
        except Exception:
            energies.append(float("nan"))
    e_arr = np.array(energies)
    i_ts = int(np.nanargmax(e_arr))
    ts_pos = images[i_ts].get_positions().copy()
    ts_nums = images[i_ts].get_atomic_numbers()
    (out / "neb_final_path.xyz").write_text(_images_to_xyz(images))
    np.save(out / "neb_energies.npy", e_arr)

    rec = {
        "engine": engine.name,
        "n_images": len(images),
        "optimizer": ("ASE FIRE (CI-NEB)" if climb_from_start else
                      "ASE FIRE (CI-NEB), preceded by FIRE and optional ANI warm start"),
        "steps_taken": step["n"],
        "fmax_final_ev_A": fmax,
        "fmax_target": args.fmax,
        "converged": bool(converged and math.isfinite(fmax) and fmax <= args.fmax),
        "energies_ev": energies,
        "ts_image_index": i_ts,
        "e_r_ev": float(energies[0]),
        "e_p_ev": float(energies[-1]),
        "e_ts_ev": float(energies[i_ts]),
        "barrier_fwd_ev": float(energies[i_ts] - energies[0]),
        "barrier_rev_ev": float(energies[i_ts] - energies[-1]),
        "reaction_energy_ev": float(energies[-1] - energies[0]),
        "wall_seconds": round(time.time() - t0, 1),
        "engine_calls": getattr(engine, "n_calls", None),
    }
    RESULTS["stage2_neb"] = rec
    ckpt.write_text(json.dumps(rec, indent=2), encoding="utf-8")
    _export_neb_candidate(out / 'ts_candidate.xyz', images[i_ts])
    _log("2", f"CI-NEB done: TS@image {i_ts}, "
              f"ΔE‡fwd {rec['barrier_fwd_ev'] * EV_TO_KCAL:.1f} kcal/mol, "
              f"fmax {fmax:.4f} eV/A, converged={rec['converged']}")
    return rec['converged']


def neb_fmax(neb):
    f = neb.get_forces()
    return float(np.sqrt((f ** 2).sum(axis=1).max()))


def stage_bounded_refinement(out, args):
    """Explicit bounded continuation of a fresh campaign band, not a full rerun."""
    import hashlib
    from ase.io import read, write
    from ase.mep.neb import NEB
    from ase.optimize import FIRE
    source = Path(args.refine_from).resolve()
    provenance = json.loads((source.parent / 'rerun_status.json').read_text())
    previous = json.loads((source / 'stage2.json').read_text())
    if provenance.get('phase') != 4 or provenance.get('historical_cache_copied') is not False:
        raise RuntimeError('Refinement requires a fresh Phase4 campaign source')
    if args.fmax != previous['fmax_target']:
        raise RuntimeError('Refinement must preserve the original force tolerance')
    hashes = {}
    for name in ('stage1.json', 'images_idpp.xyz', 'reactant_3d.mol',
                 'product_3d.mol', 'neb_final_path.xyz'):
        path = source / name
        hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
        if name != 'neb_final_path.xyz':
            shutil.copy2(path, out / name)
    RESULTS['run_scope'] = 'bounded_refinement_not_default_full'
    RESULTS['refinement_provenance'] = dict(source=str(source), sha256=hashes,
        max_steps=150, maxstep_A=.005, dt=.02, original_force_target=args.fmax)
    RESULTS['stage1_geometry'] = json.loads((out / 'stage1.json').read_text())
    images = read(str(source / 'neb_final_path.xyz'), index=':')
    engine = XTBWrap()
    for im in images:
        im.calc = _PerAtomCalc(engine)
    neb = NEB(images, climb=True, k=args.spring, method='improvedtangent')
    best = {'fmax': float('inf'), 'positions': None}
    steps = {'n': 0}

    def observe():
        fmax = neb_fmax(neb)
        if not math.isfinite(fmax):
            raise RuntimeError('Nonfinite refinement forces')
        if fmax < best['fmax']:
            best.update(fmax=fmax, positions=[im.positions.copy() for im in images])
        if steps['n'] % 10 == 0:
            _log('refine', f'step {steps["n"]}: fmax={fmax:.6f} eV/A')
        steps['n'] += 1

    opt = FIRE(neb, dt=.02, maxstep=.005, logfile=None)
    opt.attach(observe)
    opt.run(fmax=args.fmax, steps=150)
    if best['positions'] is not None:
        for im, pos in zip(images, best['positions']):
            im.positions[:] = pos
    final_forces = neb.get_forces()
    _save_neb_force_diagnostics(out, neb, final_forces, engine, args.fmax)
    fmax = float(np.sqrt((final_forces ** 2).sum(axis=1).max()))
    energies = [float(im.get_potential_energy()) for im in images]
    i_ts = int(np.argmax(energies))
    passed = bool(np.isfinite(energies).all() and math.isfinite(fmax) and fmax <= args.fmax)
    rec = dict(engine=engine.name, n_images=len(images), optimizer='bounded CI-FIRE',
        steps_taken=opt.nsteps, fmax_final_ev_A=fmax, fmax_target=args.fmax,
        converged=passed, energies_ev=energies, ts_image_index=i_ts,
        e_r_ev=energies[0], e_p_ev=energies[-1], e_ts_ev=energies[i_ts],
        barrier_fwd_ev=energies[i_ts]-energies[0], barrier_rev_ev=energies[i_ts]-energies[-1],
        reaction_energy_ev=energies[-1]-energies[0], engine_calls=engine.n_calls)
    RESULTS['stage2_neb'] = rec
    (out / 'stage2.json').write_text(json.dumps(rec, indent=2))
    (out / 'neb_final_path.xyz').write_text(_images_to_xyz(images))
    np.save(out / 'neb_energies.npy', np.array(energies))
    _export_neb_candidate(out / 'ts_candidate.xyz', images[i_ts])
    _log('refine', f'finished {opt.nsteps} steps; fmax={fmax:.6f}, accepted={passed}')
    return passed

class _PerAtomCalc:
    """Attach one engine per image (subprocess engines are stateless but this
    keeps call counters separated)."""

    def __init__(self, engine):
        self._engine = engine
        self.name = engine.name
        self._state = None
        self._values = None

    def _evaluate(self, atoms):
        # Each image owns this wrapper. Cache only identical coordinates and
        # atom identities; return copies because NEB may project force arrays.
        state = (atoms.get_atomic_numbers().tobytes(), atoms.get_positions().tobytes())
        if state != self._state:
            if isinstance(self._engine, XTBWrap):
                values = self._engine._run(atoms.get_atomic_numbers(), atoms.get_positions())
            elif isinstance(self._engine, ANIWrap):
                values = self._engine._eval(atoms)
            else:
                values = (self._engine.get_potential_energy(atoms), self._engine.get_forces(atoms))
            self._values = (float(values[0]), np.asarray(values[1]).copy())
            self._state = state
        return self._values[0], self._values[1].copy()

    def get_potential_energy(self, atoms):
        return self._evaluate(atoms)[0]

    def get_forces(self, atoms):
        return self._evaluate(atoms)[1]


# --------------------------------------------------------------------------- #
#  STAGE 3 — Hessian / 1-imaginary test / thermochemistry
# --------------------------------------------------------------------------- #

def _zpe_h_s_from_freqs(freqs_cm, temperature=298.15):
    """Harmonic vibrational thermochemistry from frequencies (cm^-1).
    Returns (ZPE kcal/mol, U_vib(T) kcal/mol INCLUDING ZPE,
    S_vib cal/mol/K, positive freqs). Translation/rotation cancel between
    R and TS of the same molecule and are omitted."""
    h = 6.62607015e-34
    kb = 1.380649e-23
    na = 6.02214076e23
    r_cal = 1.9872042586
    kT = kb * temperature
    freqs = np.array([f for f in freqs_cm if f > 0])
    if len(freqs) == 0:
        return 0.0, 0.0, 0.0, freqs
    nu = freqs * 2.99792458e10
    e_j = h * nu
    zpe = float(0.5 * e_j.sum() * na / 4184.0)          # kcal/mol
    x = e_j / kT
    u_vib = float((0.5 * e_j + e_j / np.expm1(x)).sum() * na / 4184.0)
    s_vib = float((r_cal * (x / np.expm1(x) -
                            np.log(-np.expm1(-x)))).sum())
    return zpe, u_vib, s_vib, freqs


def stage3_ts_verify(out: Path, args, force=False) -> bool:
    from phase_audit import AUDIT_VERSION, ts_validation
    ckpt = out / "stage3.json"
    if ckpt.exists() and (out / "ts_hessian_freqs.json").exists() and not force:
        cached = json.loads(ckpt.read_text())
        if cached.get('validation', {}).get('version') == AUDIT_VERSION:
            RESULTS['stage3_ts'] = cached
            return bool(cached['validation']['passed'])
        _log('3', 'Pre-audit TS checkpoint invalidated; recomputing')

    from ase.io import read as ase_read
    ts = ase_read(str(out / "ts_candidate.xyz"))
    nums = ts.get_atomic_numbers()
    pos = ts.get_positions()
    R_at = ase_read(str(out / "images_idpp.xyz"), index="0")

    if XTB_EXE:
        # NOTE: xtb's --opt ts in this build degrades to a plain minimization
        # (the separate --ts flag is rejected), destroying the saddle; the
        # converged CI image (fmax < 0.05 eV/A) is verified directly instead.
        _log("3", "numerical GFN2-xTB Hessian on the CI candidate "
                  "(xtb --hess) ...")
        h_ts = xtb_hessian(nums, pos)
        _log("3", "numerical GFN2-xTB Hessian on R ...")
        h_r = xtb_hessian(nums, R_at.get_positions())
        engine = "GFN2-xTB numerical Hessian and gradient screening; no TS refinement or IRC"
    else:
        _fallback("xtb unavailable — numerical ANI-2x Hessian (ASE "
                  "Vibrations, central differences)")
        h_ts = _ani_hessian(nums, pos)
        h_r = _ani_hessian(nums, R_at.get_positions())
        engine = "ANI-2x numerical Hessian (ASE Vibrations)"

    f_ts = sorted(h_ts["frequencies"])
    f_r = sorted(h_r["frequencies"])
    if h_ts.get("modes"):
        np.savez(out / "ts_modes.npz",
                 freqs=np.array([m["freq"] for m in h_ts["modes"]]),
                 disps=np.array([m["disp"] for m in h_ts["modes"]]))
    n_imag_ts = sum(1 for f in f_ts if f < 0)
    imag = [f for f in f_ts if f < 0]
    (out / "ts_hessian_freqs.json").write_text(
        json.dumps({"freqs_ts": f_ts, "freqs_r": f_r}, indent=2))

    # transition vector sanity: correlate the imaginary mode with
    # breaking/forming bond axes
    s1 = RESULTS["stage1_geometry"]
    broken = [tuple(x) for x in s1["broken_bonds_R_idx"]]
    formed = [tuple(x) for x in s1["formed_bonds_R_idx"]]

    def unit(a, b, p):
        v = p[b] - p[a]
        n = np.linalg.norm(v)
        return v / n if n > 1e-8 else v

    tv = None
    mode_correlation = None
    reactive_modes = [m for m in h_ts.get('modes', []) if m['freq'] < -20.]
    if len(reactive_modes) == 1:
        tv = np.asarray(reactive_modes[0]['disp'], dtype=float)
    if tv is not None and (broken or formed):
        tvn = np.linalg.norm(tv, axis=1, keepdims=True)
        tvu = tv / (tvn + 1e-12)
        cors = []
        for (i, j) in broken + formed:
            axis = unit(i, j, pos)
            relative = tv[j] - tv[i]
            proj = float(abs(np.dot(relative, axis)) / (np.linalg.norm(relative)+1e-12))
            cors.append({"bond": [int(i), int(j)], "|proj|": proj})
        mode_correlation = cors

    validator = XTBWrap() if XTB_EXE else ANIWrap()
    f_ts_max = float(np.linalg.norm(validator.get_forces(ts), axis=1).max())
    f_r_max = float(np.linalg.norm(validator.get_forces(R_at), axis=1).max())
    validation = ts_validation(f_ts, f_r, f_ts_max, f_r_max,
        RESULTS['stage2_neb'].get('converged', False),
        max((c['|proj|'] for c in (mode_correlation or [])), default=None),
        force_tolerance=args.fmax)
    if not validation['passed']:
        rec = dict(engine=engine, validation=validation, n_imaginary_ts=n_imag_ts,
                   frequencies_ts_cm=f_ts, frequencies_r_cm=f_r,
                   force_ts_ev_a=f_ts_max, force_r_ev_a=f_r_max,
                   thermochemistry_298K=None, status='invalid_ts_no_kinetic_result')
        RESULTS['stage3_ts'] = rec
        ckpt.write_text(json.dumps(rec, indent=2), encoding='utf-8')
        _warn(f'TS screening failed: {validation["checks"]}; thermochemistry/rate withheld')
        return False

    # thermochemistry (298.15 K): electronic + ZPE + H/S corrections
    T = 298.15
    zpe_ts, Uv_ts, Sv_ts, _ = _zpe_h_s_from_freqs(f_ts)
    zpe_r, Uv_r, Sv_r, _ = _zpe_h_s_from_freqs(f_r)
    s2 = RESULTS["stage2_neb"]
    dE = (validator.get_potential_energy(ts) - validator.get_potential_energy(R_at)) * EV_TO_KCAL
    dZPE = zpe_ts - zpe_r
    dH = dE + (Uv_ts - Uv_r)      # U_vib already contains ZPE
    dS_vib = (Sv_ts - Sv_r)  # cal/mol/K (vibrational only; translation and
    # rotation are omitted, not guaranteed to cancel: this is a vibrational-only estimate)
    dG = dH - T * dS_vib / 1000.0
    k_e = 2.083661912e10 * T  # k_B*T/h in s^-1
    k_rate = k_e * math.exp(-dG / (R_KCAL * T)) if math.isfinite(dG) and dG > 0 else None

    rec = {
        "engine": engine,
        "validation": validation,
        "status": "screened_candidate_not_IRC_verified",
        "force_ts_ev_a": f_ts_max, "force_r_ev_a": f_r_max,
        "n_imaginary_ts": n_imag_ts,
        "imaginary_freqs_cm": imag,
        "one_imag_criterion": validation['checks']['one_significant_imaginary'],
        "lowest_real_freq_cm": f_ts[1] if len(f_ts) > 1 else None,
        "transition_vector_check": mode_correlation,
        "thermochemistry_298K": {
            "dE_kcal": dE, "dZPE_kcal": dZPE, "dUvib_kcal": Uv_ts - Uv_r,
            "dH_kcal": dH, "TdS_kcal": T * dS_vib / 1000.0,
            "dG_kcal": dG,
            "eyring_rate_s": k_rate,
            "scope": "conditional vibrational-only estimate; not an IRC-verified reaction rate",
        },
        "frequencies_ts_cm": f_ts,
        "frequencies_r_cm": f_r,
        "xtb_thermo_ts_text": h_ts.get("thermo_text", "")[:1500],
    }
    RESULTS["stage3_ts"] = rec
    ckpt.write_text(json.dumps(rec, indent=2), encoding="utf-8")
    _log("3", f"imaginary modes: {n_imag_ts} {imag[:3]} | "
              f"ΔE‡ {dE:.1f} | ΔH‡ {dH:.1f} | ΔG‡ {dG:.1f} kcal/mol | "
              f"conditional k(298K) {k_rate} s^-1")
    return True


def _ani_hessian(nums, pos):
    """Numerical Hessian frequencies via ASE Vibrations with ANI-2x."""
    from ase import Atoms
    from ase.vibrations import Vibrations
    at = Atoms(numbers=nums, positions=pos)
    at.set_calculator(_PerAtomCalc(ANIWrap()))
    # ASE represents imaginary frequencies as complex numbers. Separate scratch
    # directories also prevent the reactant from reusing the TS displacement cache.
    with tempfile.TemporaryDirectory(prefix="p4ani_hess_") as td:
        vib = Vibrations(at, name=str(Path(td) / "vib"), delta=0.015, nfree=4)
        vib.run()
        freqs = [-float(abs(f.imag)) if abs(f.imag) > 1e-8 else float(f.real)
                 for f in vib.get_frequencies()]
        modes = []
        vib_data = vib.get_vibrations()
        for k in range(len(freqs)):
            modes.append({"freq": freqs[k],
                          "disp": np.array(vib_data.get_mode(k))})
    return {"frequencies": freqs, "modes": modes,
            "thermo_text": "", "stdout_tail": ""}

# --------------------------------------------------------------------------- #
#  STAGE 4 — figures
# --------------------------------------------------------------------------- #

def stage4_figures(fig_dir: Path) -> bool:
    from phase_audit import AUDIT_VERSION
    validation = RESULTS.get('stage3_ts', {}).get('validation', {})
    if validation.get('version') != AUDIT_VERSION or not validation.get('passed'):
        raise RuntimeError('Validated current-version TS record required; refusing misleading figures')
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fig_dir.mkdir(parents=True, exist_ok=True)
    out = _OUT
    from ase.io import read as ase_read

    s1 = RESULTS["stage1_geometry"]
    s2 = RESULTS["stage2_neb"]
    s3 = RESULTS["stage3_ts"]
    images = ase_read(str(out / "neb_final_path.xyz"), index=":")
    e = np.array(s2["energies_ev"]) * EV_TO_KCAL
    x = np.linspace(0, 1, len(e))
    # smooth spline for display
    try:
        from scipy.interpolate import PchipInterpolator
        xs = np.linspace(0, 1, 300)
        es = PchipInterpolator(x, e - e.min())(xs)
    except Exception:
        xs, es = x, e - e.min()
    i_ts = s2["ts_image_index"]

    # ---- fig1: reaction profile ------------------------------------------- #
    fig, ax = plt.subplots(figsize=(9.5, 6.4))
    ax.plot(xs, es, "-", color="#0173B2", lw=2.2, zorder=2)
    ax.plot(x, e - e.min(), "o", color="#0173B2", ms=7,
            mec="k", mew=0.6, zorder=3)
    ax.plot(x[i_ts], e[i_ts] - e.min(), "*", ms=24, color="#D62728",
            mec="k", mew=0.8, zorder=4)
    th = s3["thermochemistry_298K"]
    ax.annotate(f"TS‡ (CI {i_ts + 1}/{len(e)})\nΔE‡ = "
                f"{(s2['e_ts_ev'] - s2['e_r_ev']) * EV_TO_KCAL:.1f}\n"
                f"ΔG‡ = {th['dG_kcal']:.1f} kcal/mol",
                xy=(x[i_ts], e[i_ts] - e.min()),
                xytext=(0.42, 0.92), textcoords="axes fraction",
                fontsize=10, fontweight="bold", color="#B03A2E",
                arrowprops=dict(arrowstyle="->", color="#B03A2E"))
    ax.annotate("Reactant\ncyclopropa[b]indole", xy=(0, 0),
                xytext=(0.02, 0.10), textcoords="axes fraction",
                fontsize=9, color="0.25")
    ax.annotate(f"Product\n2,3-dihydroquinoline\nΔG_rxn = "
                f"{(s2['e_p_ev'] - s2['e_r_ev']) * EV_TO_KCAL:.1f} kcal/mol",
                xy=(1, e[-1] - e.min()),
                xytext=(0.72, 0.06), textcoords="axes fraction",
                fontsize=9, color="0.25", ha="left")
    ax.set_xlabel("normalized reaction coordinate")
    ax.set_ylabel("relative energy (kcal/mol)")
    ax.set_title("CI-NEB skeletal editing pathway — ring expansion "
                 f"({s2['engine'].split(' (')[0]}, {len(e)} images)",
                 fontsize=11.5, fontweight="bold")
    ax.set_xlim(-0.03, 1.03)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig1_neb_reaction_profile.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    _log("4", "wrote fig1_neb_reaction_profile.png")

    # ---- fig2: TS + imaginary mode arrows --------------------------------- #
    ts = images[i_ts]
    pos = ts.get_positions()
    nums = ts.get_atomic_numbers()
    sym = {1: "H", 6: "C", 7: "N"}
    tv = None
    for m in RESULTS["stage3_ts"].get("_modes_cache", []) if False else []:
        pass
    # re-derive mode: rerun hessian? use stored modes from stage3 run if saved
    modes_path = out / "ts_modes.npz"
    if modes_path.exists():
        z = np.load(modes_path)
        freqs = z["freqs"]
        disps = z["disps"]  # (n_modes, n_atoms, 3)
        k = int(np.argmin(freqs))
        if freqs[k] < 0:
            tv = disps[k]
    scale = 0.55
    fig = plt.figure(figsize=(9.8, 8.2))
    ax = fig.add_subplot(111, projection="3d")
    colors = {"C": "#222222", "N": "#2E86DE", "H": "#AAAAAA"}
    for k, (num, p) in enumerate(zip(nums, pos)):
        s_ = sym.get(num, "X")
        ax.scatter(*p, s=70 if s_ != "H" else 22,
                   c=colors.get(s_, "#7D3C98"), depthshade=False,
                   edgecolors="k", linewidths=0.5)
    if tv is not None:
        ax.quiver(pos[:, 0], pos[:, 1], pos[:, 2],
                  tv[:, 0] * scale, tv[:, 1] * scale, tv[:, 2] * scale,
                  color="#D62728", lw=1.6,
                  arrow_length_ratio=0.12, normalize=False)
    else:
        ax.set_title("(mode vectors unavailable)", fontsize=9)
    # draw the evolving key bonds
    for (i, j) in s1["broken_bonds_R_idx"][:3]:
        ax.plot(*zip(pos[i], pos[j]), color="#D62728", lw=3, alpha=0.6)
    for (i, j) in s1["formed_bonds_R_idx"][:3]:
        ax.plot(*zip(pos[i], pos[j]), color="#27AE60", lw=3, alpha=0.6)
    ax.set_box_aspect((1, 1, 1))
    pad = 2.2
    c = pos.mean(axis=0)
    ax.set_xlim(c[0] - pad, c[0] + pad)
    ax.set_ylim(c[1] - pad, c[1] + pad)
    ax.set_zlim(c[2] - pad, c[2] + pad)
    if s3.get("imaginary_freqs_cm"):
        f_im = s3["imaginary_freqs_cm"][0]
        tline = f"Transition State — imaginary mode ({f_im:.0f}i cm⁻¹)"
    else:
        tline = "Transition State — imaginary mode"
    nl = chr(10)
    ax.set_title(tline + nl +
                 "red = breaking, green = forming, arrows = ν‡ displacement",
                 fontsize=11.5, fontweight="bold")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig2_ts_vibrational_mode.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    _log("4", "wrote fig2_ts_vibrational_mode.png")

    # ---- fig3: bond evolution heatmap -------------------------------------- #
    pairs = [tuple(x) for x in s1["broken_bonds_R_idx"]] + \
            [tuple(x) for x in s1["formed_bonds_R_idx"]]
    labels = []
    for (i, j) in pairs:
        kind = "break" if (i, j) in [tuple(x) for x in
                                     s1["broken_bonds_R_idx"]] else "form"
        labels.append(f"{kind}: {sym.get(nums[i],'X')}{i}–"
                      f"{sym.get(nums[j],'X')}{j}")
    M = np.zeros((len(pairs), len(images)))
    for a, (i, j) in enumerate(pairs):
        for b, im in enumerate(images):
            p = im.get_positions()
            M[a, b] = float(np.linalg.norm(p[i] - p[j]))
    fig, ax = plt.subplots(figsize=(10.5, 1.2 + 0.5 * len(pairs) + 1.8))
    im_ = ax.imshow(M, aspect="auto", cmap="viridis_r")
    ax.set_xticks(range(len(images)))
    ax.set_xticklabels([f"{k+1}" for k in range(len(images))], fontsize=9)
    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels(labels, fontsize=9)
    for a in range(len(pairs)):
        for b in range(len(images)):
            ax.text(b, a, f"{M[a, b]:.2f}", ha="center", va="center",
                    fontsize=7.5,
                    color="w" if M[a, b] > M.mean() else "k")
    ax.set_xlabel("NEB image (1 = R, last = P)")
    ax.set_title("Bond evolution across the NEB band (Å): scission of the "
                 "target bond & formation of the inserted-ring bond",
                 fontsize=11, fontweight="bold")
    fig.colorbar(im_, ax=ax, label="distance (Å)", shrink=0.85)
    fig.tight_layout()
    fig.savefig(fig_dir / "fig3_bond_evolution_matrix.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    _log("4", "wrote fig3_bond_evolution_matrix.png")
    return True


# --------------------------------------------------------------------------- #
#  main
# --------------------------------------------------------------------------- #

def _shutdown(delay_s: int) -> None:
    _log("shutdown", f"initiating shutdown in {delay_s} s "
                     f"(--auto_shutdown set, all stages succeeded)")
    subprocess.run(["shutdown", "/s", "/t", str(delay_s)],
                   check=False, capture_output=True)


def main() -> int:
    p = argparse.ArgumentParser(
        description="Phase 4: CI-NEB skeletal editing pathway + TS "
                    "verification")
    p.add_argument("--out_dir", default="results_phase4")
    p.add_argument("--fig_dir", default="figures_phase4")
    p.add_argument("--engine", default="auto", choices=["auto", "xtb", "ani"])
    p.add_argument("--n_images", type=int, default=9)
    p.add_argument("--fmax", type=float, default=0.05)
    p.add_argument("--spring", type=float, default=0.1)
    p.add_argument("--max_neb_steps", type=int, default=500)
    p.add_argument("--force_rerun", action="store_true")
    p.add_argument('--refine_from', help='Fresh campaign results_phase4 directory; bounded150-step CI-FIRE continuation')
    p.add_argument("--fig_only", action="store_true")
    p.add_argument("--auto_shutdown", action="store_true")
    p.add_argument("--shutdown_delay", type=int, default=60)
    args = p.parse_args()

    global _OUT, _FIG
    _OUT = Path(args.out_dir)
    _FIG = Path(args.fig_dir)
    _OUT.mkdir(parents=True, exist_ok=True)
    _FIG.mkdir(parents=True, exist_ok=True)

    code = 0
    try:
        _log("main", f"PHASE 4 skeletal editing CI-NEB — "
                     f"{RESULTS['timestamp']} | engine={args.engine} | "
                     f"images={args.n_images} | fmax={args.fmax}")
        if args.fig_only:
            for st, key in (("stage1", "stage1_geometry"),
                            ("stage2", "stage2_neb"),
                            ("stage3", "stage3_ts")):
                f = _OUT / f"{st}.json"
                if f.exists():
                    RESULTS[key] = json.loads(f.read_text())
            stage4_figures(_FIG)
            return 0
        stages = [
            ("1 geometry+pairing+IDPP", lambda: stage1_geometry(
                _OUT, args, args.force_rerun)),
            ("2 CI-NEB", lambda: stage2_neb(_OUT, args, args.force_rerun)),
            ("3 TS verify", lambda: stage3_ts_verify(
                _OUT, args, args.force_rerun)),
        ]
        if args.refine_from:
            stages = [("2 bounded fresh-band refinement", lambda: stage_bounded_refinement(_OUT, args)),
                      ("3 TS verify", lambda: stage3_ts_verify(_OUT, args, True))]
        for name, fn in stages:
            print("=" * 72, flush=True)
            print(f"STAGE {name}", flush=True)
            print("=" * 72, flush=True)
            try:
                ok = fn()
                if not ok:
                    code = 1
                    _warn(f'{name}: validation failed; downstream stages stopped')
                    break
                _log("stage", f"{name}: {'OK' if ok else 'SKIPPED'}")
            except Exception as exc:
                code = 1
                _log("stage", f"{name}: FAILED — "
                              f"{exc.__class__.__name__}: {exc}")
                traceback.print_exc()
                RESULTS["fatal_error"] = f"{name}: {type(exc).__name__}: {exc}"
                break
        try:
            stage4_figures(_FIG)
        except Exception as exc:
            _warn(f"figure generation failed: {exc}")
            traceback.print_exc()
        RESULTS["all_stages_ok"] = (code == 0)
    except Exception as exc:
        code = 1
        RESULTS["fatal_error"] = f"{exc.__class__.__name__}: {exc}"
        traceback.print_exc()
    finally:
        print("=" * 72, flush=True)
        print("FINALIZING", flush=True)
        print("=" * 72, flush=True)
        try:
            write_json_atomic(_OUT / "phase4_results.json")
            _log("final", f"results -> {_OUT / 'phase4_results.json'}")
        except Exception as exc:
            _log("final", f"serialization failed: {exc}")
        if args.auto_shutdown and RESULTS["all_stages_ok"]:
            _shutdown(args.shutdown_delay)
        else:
            _log("final", "shutdown NOT triggered (disabled or failure)")
    return code


if __name__ == "__main__":
    sys.exit(main())
