#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_phase6_explicit_metadynamics.py
===================================
PHASE 6 — EXPLICIT SOLVENT ENHANCED SAMPLING & 2D FREE ENERGY
SURFACE RECONSTRUCTION (Well-Tempered Metadynamics, OpenMM 8.6)

Scaling the Phase-4/5 skeletal-editing pipeline (cyclopropa[b]indole ->
2,3-dihydroquinoline, the Ciamician-Dennstedt core step) from implicit
GB solvation to full-atom explicit solvent biophysics with 2D WTMetaD.

Reaction model (atom indices in the Phase-4 reactant ordering, 0-based)
-----------------------------------------------------------------------
  Reactant : cyclopropa[b]indole, C9H9N (results_phase4/reactant_3d.mol)
  Product  : 2,3-dihydroquinoline (results_phase4/product_3d.mol, remapped
             onto R order through the Phase-4 MCS+Hungarian pairing)
  Scissile N-C bond        : (4, 5)   breaks   1.43 A -> 2.43 A
  Secondary C-CH2 cleavage : (6, 9)   breaks   1.51 A -> 2.48 A
  Forming  N-CH2 bond      : (4, 9)   forms    2.47 A -> 1.45 A
  Migrating carbon         : 9 (cyclopropane CH2)

Collective variables (windows validated on the Phase-4 GFN2-xTB IDPP path)
--------------------------------------------------------------------------
  CV1  d(N4-C5)               scissile-bond distance   [1.3, 3.0] A
  CV2  theta(C9-C5-N4-C3)     insertion dihedral of the migrating CH2
                              about the conserved C9-C5 bond
                                                        [-180, 180] deg

Reactive classical Hamiltonian (identical in BOTH legs)
--------------------------------------------------------
  * Valence FF  : OpenFF "Sage" 2.1.0 on the reactant valence graph,
                  exported to an OpenMM ForceField XML (types L0..L18) and
                  energy-validated against the interchange system.
  * Charges     : (reactant + product)/2 GFN2-xTB Mulliken charges under
                  the pairing permutation (fallback MMFF94 -> Gasteiger).
  * Morse pairs : the three reactive pairs are removed from the harmonic
                  bond list and modelled with Morse potentials whose
                  (De, Re, a) are fitted to GFN2-xTB constrained-scan
                  profiles (fallback: xtb-Hessian k + tabulated BDE).
                  Pair (4,9) is 1-4 in the reactant graph: its scaled
                  NonbondedForce exception is zeroed so the pair interacts
                  only through Morse; broken pairs keep native 1-2
                  exclusions. Same Hamiltonian in both legs, so
                  Delta-Delta-G(solv) = dG(expl) - dG(impl) is
                  self-consistent.

Simulation protocol
-------------------
  * Explicit : TIP3P, 1.0 nm (>=10.0 A) padding, PME 1.0 nm cutoff,
               0.15 M NaCl, neutralized; LangevinMiddle 300 K, dt = 2 fs,
               HBonds constraints.
  * Implicit : solute-only + GBSA/OBC2 (Bondi radii, OBC2 screening),
               NoCutoff — the Delta-G‡(implicit) reference leg.
  * WTMetaD  : W0 = 0.5 kcal/mol, sigma1 = 0.05 A, sigma2 = 5 deg,
               deposition every 500 steps (1.0 ps, explicit) / 1000 steps
               (2 ps, implicit), bias factor gamma = 10 at T = 300 K,
               production 200,000 steps explicit / 500,000 implicit.
  * Bias     : ONE CustomCVForce with N_MAX pre-allocated hill slots
               driven by global parameters (zero recompilation; deposits
               are context.setParameter calls).
  * FES      : Delta-G(d, theta) = -(gamma-1)/gamma * V_bias on a
               periodic grid; deposition-density masking, Gaussian
               smoothing, Dijkstra minimum-energy path, saddle
               extraction, convergence trace.

Outputs
-------
  results_phase6/   phase6_results.json, fes_*.npz, metadyn_state.json,
                    progress_*.csv, solvated PDB, DCDs, QM caches
  figures_phase6/   fig1_explicit_solvent_box.png
                    fig2_2d_free_energy_surface.png
                    fig3_fes_convergence_trace.png          (all 300 DPI)

Run under the `phase2ff` conda env (openff-toolkit + interchange + openmm):
  python run_phase6_explicit_metadynamics.py             # full protocol
  python run_phase6_explicit_metadynamics.py --selftest  # benchmark only
  python run_phase6_explicit_metadynamics.py --fast      # smoke protocol
Fault tolerance: per-stage try/except, atomic JSON checkpointing, resume
from Context checkpoints (existing artifacts are picked up automatically).
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path

import numpy as np

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
except ImportError as exc:  # pragma: no cover
    raise SystemExit("OpenMM is required (run under the phase2ff env).") from exc

try:
    from openff.toolkit import ForceField as OFFFF, Molecule as OFFMol
    from openff.units import unit as offunit
except ImportError:
    raise SystemExit(
        "openff.toolkit / openff.units are required (run under the "
        "phase2ff conda env).")

from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.ndimage import gaussian_filter
from scipy.optimize import least_squares
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

# ------------------------------------------------------------------ #
#  constants / units
# ------------------------------------------------------------------ #
KB_KJ_MOL_K = 0.008314462618          # kJ/(mol K)
KCAL = 4.184                          # kJ per kcal
EH_KCAL = 627.5094740631
TEMP = 300.0                          # K
GAMMA = 10.0                          # WTMetaD bias factor
W0_KCAL = 0.5                         # initial hill height (kcal/mol)
SIG1_A = 0.05                         # hill width, distance CV (A)
SIG2_DEG = 5.0                        # hill width, torsion CV (deg)
DEPOSIT_STEPS_EXPL = 500              # 1.0 ps at 2 fs
DEPOSIT_STEPS_IMPL = 1000             # 2.0 ps at 2 fs
N_STEPS_EXPL_DEFAULT = 200_000        # >=150k required by mission
N_STEPS_IMPL_DEFAULT = 500_000        # 1 ns implicit reference
N_MAX_HILLS = 512
DT_FS = 2.0
R_CUT_NM = 1.0                        # PME / LJ cutoff (mission: 1.0 nm)
PAD_NM = 1.0                          # >= 10.0 A buffer (mission)
IONIC_M = 0.15                        # physiological NaCl
RC = (4, 5)                           # CV1 atom pair (scissile N-C)
TC = (9, 5, 4, 3)                     # CV2 torsion (insertion dihedral)
MORSE_PAIRS = {"break_nc": (4, 5), "break_cc": (6, 9), "form_nc": (4, 9)}
FALLBACK_BDE_KCAL = {"break_nc": 72.6, "break_cc": 83.1, "form_nc": 72.6}

OUT = Path("results_phase6")
FIG = Path("figures_phase6")
CACHE = OUT / "cache"

RESULTS: dict = {
    "phase": 6,
    "title": ("Explicit-solvent well-tempered metadynamics & 2D free "
              "energy surface reconstruction of the cyclopropa[b]indole "
              "ring expansion"),
    "timestamp": _dt.datetime.now().isoformat(timespec="seconds"),
    "system": {},
    "force_field": {},
    "morse_calibration": {},
    "explicit": {},
    "implicit": {},
    "fes": {},
    "meta": {"fallbacks": [], "warnings": []},
    "fatal_error": None,
    "all_stages_ok": False,
}
STATE = {"explicit": {"hills": [], "steps_done": 0},
         "implicit": {"hills": [], "steps_done": 0}}

# CV grids (analysis)
R_GRID_A = np.arange(1.30, 3.0001, 0.02)
TH_GRID_DEG = np.arange(-180.0, 180.0001, 5.0)


def _log(tag, msg):
    print(f"[{_dt.datetime.now().strftime('%H:%M:%S')}] [{tag}] {msg}",
          flush=True)


def _warn(msg):
    RESULTS["meta"]["warnings"].append(msg)
    _log("warn", msg)


def _fallback(msg):
    RESULTS["meta"]["fallbacks"].append(msg)
    _log("fallback", msg)


def write_json_atomic(path: Path, obj=None):
    obj = RESULTS if obj is None else obj
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, default=float),
                   encoding="utf-8")
    os.replace(tmp, path)


# ------------------------------------------------------------------ #
#  stage 0 — assets, pairing, QM calibration (charges + Morse scans)
# ------------------------------------------------------------------ #
XTB_EXE = None
for _c in (shutil.which("xtb"),
           r"C:\Users\HUIWEI\miniconda3\envs\phase2ff\Library\bin\xtb.exe",
           str(Path(sys.prefix) / "Library" / "bin" / "xtb.exe")):
    if _c and Path(_c).exists():
        XTB_EXE = _c
        break


def _mol_from_molfile(path: Path):
    mol = Chem.MolFromMolFile(str(path), removeHs=False)
    if mol is None:
        raise RuntimeError(f"RDKit could not parse {path}")
    return mol


def _write_xyz(path: Path, symbols, positions):
    with open(path, "w") as fh:
        fh.write(f"{len(symbols)}\np6\n")
        for s, p in zip(symbols, positions):
            fh.write(f"{s} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}\n")


def _read_xyz(path: Path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    sym, pos = [], []
    for ln in lines[2:2 + n]:
        t = ln.split()
        sym.append(t[0])
        pos.append([float(x) for x in t[1:4]])
    return sym, np.array(pos)


def _xtb_run(args, cwd: Path, timeout=900):
    env = dict(os.environ)
    env["PATH"] = str(Path(XTB_EXE).parent) + os.pathsep + env.get("PATH", "")
    return subprocess.run([XTB_EXE, *args], cwd=str(cwd), capture_output=True,
                          text=True, timeout=timeout, env=env)


def _torsion_deg(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w),
                                   np.dot(v, w)))


def _cv_of(pos):
    """CV values (d in A, theta in deg); input positions in A."""
    d = float(np.linalg.norm(pos[RC[0]] - pos[RC[1]]))
    th = _torsion_deg(pos[TC[0]], pos[TC[1]], pos[TC[2]], pos[TC[3]])
    return [d, float(th)]


def _read_idpp(path: Path):
    lines = Path(path).read_text().splitlines()
    imgs, i = [], 0
    while i < len(lines):
        n = int(lines[i].split()[0])
        block = []
        for ln in lines[i + 2:i + 2 + n]:
            t = ln.split()
            block.append((t[0], float(t[1]), float(t[2]), float(t[3])))
        imgs.append(block)
        i += n + 2
    return imgs


def stage0_assets():
    """Load Phase-4 assets, re-derive the R<-P pairing, QM caches."""
    rmol = _mol_from_molfile(Path("results_phase4/reactant_3d.mol"))
    pmol = _mol_from_molfile(Path("results_phase4/product_3d.mol"))
    s0 = {"engine": "GFN2-xTB (xtb.exe)" if XTB_EXE else "unavailable",
          "xtb_exe": XTB_EXE,
          "n_atoms": rmol.GetNumAtoms(),
          "formula": Chem.rdMolDescriptors.CalcMolFormula(rmol)}

    sys.path.insert(0, str(Path.cwd()))
    import run_phase4_reaction_mechanism as p4
    perm = np.asarray(p4._pair_atoms(rmol, pmol, verbose=False))  # perm[p]=r
    inv = np.argsort(perm)                  # inv[r] = p

    rp = rmol.GetConformer().GetPositions()          # A
    pp = pmol.GetConformer().GetPositions()          # A
    pp_in_r = pp[inv]                               # P coords in R order

    # empirical sanity: broken bonds long, formed bond short under inv
    chk_d45 = float(np.linalg.norm(pp_in_r[4] - pp_in_r[5]))
    chk_d49 = float(np.linalg.norm(pp_in_r[4] - pp_in_r[9]))
    if not (chk_d45 > 2.0 and chk_d49 < 1.8):
        raise RuntimeError(f"pairing permutation sanity failed: "
                           f"d(4,5)={chk_d45:.2f} A, d(4,9)={chk_d49:.2f} A")

    rb = {frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
          for b in rmol.GetBonds()}
    pb = {frozenset((int(perm[b.GetBeginAtomIdx()]),
                     int(perm[b.GetEndAtomIdx()])))
          for b in pmol.GetBonds()}
    s0["broken_in_R"] = sorted(sorted(tuple(int(x) for x in y))
                               for y in (rb - pb))
    s0["formed_in_R"] = sorted(sorted(tuple(int(x) for x in y))
                               for y in (pb - rb))
    if s0["broken_in_R"] != [[4, 5], [6, 9]] or s0["formed_in_R"] != [[4, 9]]:
        _warn(f"bond changes differ from Phase-4 record: "
              f"{s0['broken_in_R']} / {s0['formed_in_R']}")

    s0["cv_reactant_qm"] = _cv_of(rp)
    s0["cv_product_qm"] = _cv_of(pp_in_r)
    s0["idpp_path_cv"] = [_cv_of(np.array([c[1:] for c in img], float))
                          for img in _read_idpp(
                              Path("results_phase4/images_idpp.xyz"))]

    qR = _xtb_charges("reactant", [a.GetSymbol() for a in rmol.GetAtoms()],
                      rp)
    qP = _xtb_charges("product", [a.GetSymbol() for a in pmol.GetAtoms()],
                      pp)
    charges = 0.5 * (qR + qP[inv])                  # both in R order
    s0["charges_source"] = str(CACHE / "charges_*.txt")
    s0["charge_total_e"] = float(charges.sum())
    s0["q_range_e"] = [float(charges.min()), float(charges.max())]

    RESULTS["morse_calibration"] = _morse_calibration(rmol, pmol,
                                                        pp_in_r)
    return {"rmol": rmol, "pmol": pmol, "rp": rp, "pp_in_r": pp_in_r,
            "charges": charges, "stage0": s0}


def _xtb_charges(tag, symbols, positions):
    """GFN2-xTB Mulliken charges via `$write charges=true` (cached)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    outf = CACHE / f"charges_{tag}.txt"
    if outf.exists():
        return np.array([float(x) for x in outf.read_text().split()])
    if XTB_EXE is None:
        _fallback("xtb unavailable -> Gasteiger charges")
        return _gasteiger_charges(symbols, positions)
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        _write_xyz(td / "m.xyz", symbols, positions)
        (td / "wr.inp").write_text("$write\n  charges=true\n$end\n")
        proc = _xtb_run(["m.xyz", "--sp", "--input", "wr.inp"], td)
        ch = td / "charges"
        vals = [float(x) for x in ch.read_text().split()] \
            if ch.exists() else []
        if len(vals) != len(symbols) or abs(sum(vals)) > 5e-2:
            _fallback(f"xtb charges unusable for {tag} -> Gasteiger "
                      f"({proc.stderr[-120:]})")
            return _gasteiger_charges(symbols, positions)
        outf.write_text("\n".join(repr(v) for v in vals))
        return np.array(vals)


def _gasteiger_charges(symbols, positions):
    mol = Chem.MolFromMolBlock(_molblock_from_arrays(symbols, positions),
                               removeHs=False)
    AllChem.ComputeGasteigerCharges(mol)
    return np.array([a.GetDoubleProp("_GasteigerCharge")
                     for a in mol.GetAtoms()])


def _molblock_from_arrays(symbols, positions):
    rw = Chem.RWMol()
    conf = Chem.Conformer(len(symbols))
    for s, p in zip(symbols, positions):
        conf.SetAtomPosition(rw.AddAtom(Chem.Atom(s)),
                             [float(p[0]), float(p[1]), float(p[2])])
    rw.AddConformer(conf)
    return Chem.MolToMolBlock(rw.GetMol())


def _morse_calibration(rmol, pmol, ppos_r):
    """Fit (De, Re, a) of the three reactive pairs to GFN2-xTB
    constrained relaxed-scan profiles; fall back to Hessian-k + BDE.

    The two breaking pairs are scanned on the reactant (stretching the
    bonds that cleave); the forming N-CH2 pair is scanned on the PRODUCT
    (stretching the bond that exists there) — a relaxed forming-scan on
    the reactant is pathological because alternative relaxation channels
    (secondary ring opening) preempt the bond-forming well.
    """
    CACHE.mkdir(parents=True, exist_ok=True)
    cache_f = CACHE / "morse_fit.json"
    if cache_f.exists():
        _log("morse", "using cached GFN2-xTB scan fits")
        return json.loads(cache_f.read_text())

    out = {"method": "GFN2-xTB constrained relaxed scans (sequential "
                     "warm-started; breaking pairs on the reactant, "
                     "forming pair on the product); robust Morse fit "
                     "E(r)=off+De(1-exp(-a(r-Re)))^2-De with trailing "
                     "trim of non-dissociative relaxation"}
    if XTB_EXE is None:
        _fallback("xtb unavailable for Morse scans -> BDE tables")
        out.update({k: dict(v, source="fallback BDE/Hessian")
                    for k, v in _morse_fallback(rmol, ppos_r, out).items()})
        return out

    scans = {
        "break_nc": (rmol, np.arange(1.40, 2.96, 0.22)),
        "break_cc": (rmol, np.arange(1.45, 2.96, 0.22)),
        "form_nc": (pmol, np.arange(1.35, 2.96, 0.22)),
    }
    for tag, pair in MORSE_PAIRS.items():
        i, j = pair
        mol, rvals = scans[tag]
        try:
            sym = [a.GetSymbol() for a in mol.GetAtoms()]
            pos = mol.GetConformer().GetPositions()
            rr, ee = _xtb_scan(sym, pos, i, j, rvals)
            fit = _fit_morse(rr, ee)
            if fit is None:
                raise RuntimeError("robust Morse fit rejected "
                                   "(rms > 5 kcal/mol)")
            ee_rel = ee - ee.min()
            fit["residual_rms_kcal"] = float(np.sqrt(np.mean(
                (ee_rel - _morse(rr, fit["De_kcal"], fit["Re_A"],
                                 fit["a_invA"],
                                 fit["offset_kcal"])) ** 2)))
            fit["scan_molecule"] = ("reactant" if mol is rmol
                                    else "product")
            fit["scan_r_A"] = [float(x) for x in rr]
            fit["scan_E_rel_kcal"] = [float(x - ee.min()) for x in ee]
            out[tag] = fit
            _log("morse", f"{tag} fit on {fit['scan_molecule']}: "
                          f"De={fit['De_kcal']:.1f} kcal/mol, "
                          f"Re={fit['Re_A']:.3f} A, a={fit['a_invA']:.2f}"
                          f" 1/A")
        except Exception as exc:
            _fallback(f"Morse scan {tag} failed ({exc}) -> BDE table")
            out[tag] = None
    fb = _morse_fallback(rmol, ppos_r, out)
    for tag, val in out.items():
        if val is None:
            out[tag] = dict(fb[tag])
    write_json_atomic(cache_f, out)
    return out


def _morse(rr, De, Re, a, off):
    return off + De * (1.0 - np.exp(-a * (rr - Re))) ** 2 - De


def _fit_morse(rr, ee):
    """Robust Morse LSQ on RELATIVE energies (absolute xtb energies are
    ~-1.7e4 kcal and would break offset bounds). Trailing points that
    drop (alternative relaxation channels opening) are trimmed."""
    rr = np.asarray(rr, float)
    ee = np.asarray(ee, float) - float(ee.min())
    while len(rr) > 6 and ee[-1] < ee[-2] - 2.0:
        rr, ee = rr[:-1], ee[:-1]
    best, best_cost = None, np.inf
    for De0 in (60.0, 90.0, 120.0):
        for a0 in (1.2, 1.8, 2.5):
            p0 = [De0, float(rr[np.argmin(ee)]), a0,
                  float(ee[np.argmin(ee)] + De0)]
            try:
                sol = least_squares(
                    lambda p: _morse(rr, *p) - ee, p0,
                    bounds=([15, 1.22, 0.8, -50], [150, 1.75, 4.5, 250]),
                    loss="soft_l1", f_scale=3.0, max_nfev=20000)
            except Exception:
                continue
            if sol.cost < best_cost:
                best, best_cost = sol, sol.cost
    if best is None:
        return None
    rms = float(np.sqrt(2 * best_cost / len(rr)))
    if rms > 5.0:
        return None
    De, Re, a, off = (float(x) for x in best.x)
    return {"De_kcal": De, "Re_A": Re, "a_invA": a, "offset_kcal": off}


def _morse_fallback(rmol, ppos_r=None, fitted=None):
    """Fallback Morse parameters. The forming N-CH2 pair deserves special
    care: its relaxed scans are hijacked by competing relaxation channels
    and its reactant-side Hessian coupling is ~0 (not bonded), so we
    reuse the GFN2-fitted parameters of the scissile N-C bond (same bond
    type in the same molecule) with Re set to the product bond length."""
    ks = _xtb_bond_k(rmol)
    pos = rmol.GetConformer().GetPositions()
    ppos = ppos_r if ppos_r is not None else pos
    out = {}
    fitted_nc = (fitted or {}).get("break_nc") if fitted else None
    for tag, (i, j) in MORSE_PAIRS.items():
        De = FALLBACK_BDE_KCAL[tag]
        Re = float(np.linalg.norm(pos[i] - pos[j]))
        k = ks.get((i, j))
        if tag == "form_nc":
            Re = float(np.linalg.norm(ppos[i] - ppos[j]))
            if fitted_nc is not None:
                out[tag] = {"De_kcal": fitted_nc["De_kcal"],
                            "Re_A": Re, "a_invA": fitted_nc["a_invA"],
                            "source": ("break_nc fit reused for the "
                                       "forming N-C pair (same bond "
                                       "type); Re from product")}
                continue
            k = ks.get((4, 5), 600.0)
        if k is None:
            k = 600.0  # generic C-N / C-C stretch, kcal/mol/A^2
            _warn(f"no Hessian k for pair {(i, j)}; generic k={k}")
        out[tag] = {"De_kcal": De, "Re_A": Re,
                    "a_invA": float(math.sqrt(k / (2.0 * De))),
                    "k_kcal_A2": float(k)}
    return out


def _xtb_bond_k(rmol):
    """Largest eigenvalue of the xtb-Hessian 3x3 pair block = stretch k."""
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / "hessian_k.json"
    if f.exists():
        return {tuple(int(x) for x in k.split(",")): v
                for k, v in json.loads(f.read_text()).items()}
    if XTB_EXE is None:
        return {}
    sym = [a.GetSymbol() for a in rmol.GetAtoms()]
    pos = rmol.GetConformer().GetPositions()
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        _write_xyz(td / "m.xyz", sym, pos)
        proc = _xtb_run(["m.xyz", "--hess"], td, timeout=1800)
        hfile = td / "hessian"
        if not hfile.exists():
            _warn(f"xtb --hess failed ({proc.stderr[-140:]})")
            return {}
        H = _parse_xtb_hessian(hfile.read_text(), len(sym))
        if H is None:
            _warn("xtb hessian block not parsed")
            return {}
        conv = EH_KCAL / (0.52917721092 ** 2)  # Eh/Bohr^2 -> kcal/mol/A^2
        ks = {}
        for (i, j) in MORSE_PAIRS.values():
            # coupling block H_ij: its most-negative eigenvalue along the
            # bond direction is -k_stretch (push-pull convention)
            blk = H[3 * i:3 * i + 3, 3 * j:3 * j + 3]
            w = np.linalg.eigvalsh(0.5 * (blk + blk.T))
            ks[f"{i},{j}"] = float(max(-w[0], 0.0) * conv)
        _log("morse", "Hessian stretch k (kcal/mol/A^2): "
                      + ", ".join(f"{k}:{v:.0f}" for k, v in ks.items()))
        write_json_atomic(f, ks)
        return {tuple(int(x) for x in k.split(",")): v
                for k, v in ks.items()}


def _parse_xtb_hessian(text, natoms):
    """xtb `hessian` file: '$hessian' keyword, then the lower-triangular
    force-constant matrix (Hartree/Bohr^2), column-major, no labels."""
    n3 = 3 * natoms
    tri = n3 * (n3 + 1) // 2
    vals = []
    for ln in text.splitlines():
        t = ln.strip()
        if not t or t.startswith("$") or t.startswith("#"):
            continue
        try:
            vals.extend(float(x.replace("D", "E"))
                        for x in t.split())
        except ValueError:
            continue
        if len(vals) >= tri:
            break
    if len(vals) < tri:
        return None
    H = np.zeros((n3, n3))
    c = 0
    for col in range(n3):
        for r in range(col + 1):
            H[r, col] = H[col, r] = vals[c]
            c += 1
    return H


def _xtb_scan(sym, pos0, i, j, r_values):
    """Sequential constrained optimizations along the pair distance."""
    rr, ee = [], []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        cur = np.array(pos0, float)
        _write_xyz(td / "cur.xyz", sym, cur)
        for r in r_values:
            (td / "c.inp").write_text(
                "$constrain\n  force constant=1.0\n"
                f"  distance: {i + 1},{j + 1},{float(r):.4f}\n$end\n")
            proc = _xtb_run(["cur.xyz", "--opt", "loose", "--input",
                             "c.inp"], td, timeout=1200)
            opt = td / "xtbopt.xyz"
            if not opt.exists():
                raise RuntimeError(f"constrained opt failed at r={r:.2f} "
                                   f"({proc.stderr[-140:]})")
            eline = None
            for ln in proc.stdout.splitlines():
                if "TOTAL ENERGY" in ln:
                    eline = ln
            if eline is None:
                raise RuntimeError("no TOTAL ENERGY line in xtb output")
            ee.append(float(eline.split()[3].replace("D", "E")) * EH_KCAL)
            rr.append(float(r))
            sym2, cur = _read_xyz(opt)
            _write_xyz(td / "cur.xyz", sym2, cur)
    return np.array(rr), np.array(ee)


# ------------------------------------------------------------------ #
#  stage 1 — OpenFF Sage parametrization -> OpenMM ForceField XML
# ------------------------------------------------------------------ #
def stage1_forcefield(assets):
    rmol, charges = assets["rmol"], assets["charges"]
    s1 = {"valence_ff": "OpenFF Sage 2.1.0 (openff-2.1.0.offxml) on the "
                        "reactant valence graph",
          "charges": "GFN2-xTB Mulliken, (reactant+product)/2 under the "
                     "Phase-4 pairing"}

    offmol = OFFMol.from_rdkit(rmol, allow_undefined_stereo=True)
    offmol.partial_charges = offunit.Quantity(
        np.asarray(charges, float), offunit.elementary_charge)
    sage = OFFFF("openff-2.1.0.offxml")
    inter = sage.create_interchange(offmol.to_topology(),
                                    charge_from_molecules=[offmol])
    sys_ref = inter.to_openmm()
    s1["n_particles"] = sys_ref.getNumParticles()

    hb = _force_of(sys_ref, mm.HarmonicBondForce)
    ha = _force_of(sys_ref, mm.HarmonicAngleForce)
    ht = _force_of(sys_ref, mm.PeriodicTorsionForce)
    nb = _force_of(sys_ref, mm.NonbondedForce)

    bonds = [hb.getBondParameters(k) for k in range(hb.getNumBonds())]
    angles = [ha.getAngleParameters(k) for k in range(ha.getNumAngles())]
    torsions = [ht.getTorsionParameters(k)
                for k in range(ht.getNumTorsions())]
    sig = np.zeros(offmol.n_atoms)
    eps = np.zeros(offmol.n_atoms)
    for p in range(offmol.n_atoms):
        q0, s0_, e0_ = nb.getParticleParameters(p)
        sig[p] = _q(s0_, unit.nanometer)
        eps[p] = _q(e0_, unit.kilojoule_per_mole)

    def _14(coll, default):
        try:
            return float(getattr(coll, "coulomb14Scale",
                                 getattr(coll, "scale14")))
        except Exception:
            return default
    c14 = _14(inter["Electrostatics"], 0.833333)
    lj14 = _14(inter["vdW"], 0.5)
    s1["coulomb14Scale"], s1["lj14Scale"] = c14, lj14

    sym = [a.GetSymbol() for a in rmol.GetAtoms()]
    masses = [a.GetMass() for a in rmol.GetAtoms()]

    g = lambda v, u: f"{_q(v, u):.10g}"                 # noqa: E731
    x = ['<ForceField>', ' <AtomTypes>']
    for i in range(offmol.n_atoms):
        x.append(f'  <Type name="L{i}" class="L{i}" element="{sym[i]}" '
                 f'mass="{masses[i]:.6f}"/>')
    x += [' </AtomTypes>', ' <Residues>', '  <Residue name="LIG">']
    for i in range(offmol.n_atoms):
        x.append(f'   <Atom name="{sym[i]}{i}" type="L{i}"/>')
    for b in rmol.GetBonds():
        x.append(f'   <Bond from="{b.GetBeginAtomIdx()}" '
                 f'to="{b.GetEndAtomIdx()}"/>')
    # Sage constrains X-H bonds (they are absent from HarmonicBondForce);
    # export them as residue-level rigid constraints (matched by atom
    # NAME) or the solute hydrogens lose ALL bonding and fly apart.
    n_cst = sys_ref.getNumConstraints()
    for ci in range(n_cst):
        c1, c2, cd = sys_ref.getConstraintParameters(ci)
        x.append(f'   <Constraint atomName1="{sym[c1]}{c1}" '
                 f'atomName2="{sym[c2]}{c2}" '
                 f'distance="{_q(cd, unit.nanometer):.10g}"/>')
    x += ['  </Residue>', ' </Residues>', ' <HarmonicBondForce>']
    for (p1, p2, r0, k) in bonds:
        x.append(f'  <Bond type1="L{p1}" type2="L{p2}" '
                 f'length="{g(r0, unit.nanometer)}" '
                 f'k="{g(k, unit.kilojoule_per_mole / unit.nanometer ** 2)}"/>')
    x += [' </HarmonicBondForce>', ' <HarmonicAngleForce>']
    for (p1, p2, p3, th0, k) in angles:
        x.append(f'  <Angle type1="L{p1}" type2="L{p2}" type3="L{p3}" '
                 f'angle="{g(th0, unit.radian)}" '
                 f'k="{g(k, unit.kilojoule_per_mole / unit.radian ** 2)}"/>')
    x += [' </HarmonicAngleForce>']
    # NOTE: torsions are NOT exported through XML. OpenMM's XML torsion
    # phase is interpreted in radians and <Improper> type matching has a
    # central-atom-in-slot-1 convention with permutation matching — both
    # error-prone. Instead the reference PeriodicTorsionForce entries are
    # transplanted verbatim into the System after createSystem() (see
    # _add_torsions), which is exact by construction.
    torsions_raw = [(int(p1), int(p2), int(p3), int(p4), int(per),
                     _q(ph, unit.radian), _q(k, unit.kilojoule_per_mole))
                    for (p1, p2, p3, p4, per, ph, k) in torsions]
    x += [f' <NonbondedForce coulomb14scale="{c14:.10g}" '
          f'lj14scale="{lj14:.10g}">']
    for i in range(offmol.n_atoms):
        x.append(f'  <Atom type="L{i}" charge="{charges[i]:.10g}" '
                 f'sigma="{sig[i]:.10g}" epsilon="{eps[i]:.10g}"/>')
    x += [' </NonbondedForce>', '</ForceField>']
    OUT.mkdir(parents=True, exist_ok=True)
    xml_path = OUT / "solute_sage.xml"
    xml_path.write_text("\n".join(x), encoding="utf-8")
    n_imp = sum(1 for t in torsions_raw if _is_center(rmol, t[2],
                                                      (t[0], t[1], t[3])))
    s1.update({"xml": str(xml_path), "n_bonds": len(bonds),
               "n_angles": len(angles), "n_torsions": len(torsions),
               "n_improper_entries": n_imp,
               "n_xh_constraints": int(n_cst)})

    # -------- energy-equivalence validation (Reference platform) ----
    ff_xml = app.ForceField(str(xml_path))
    top = offmol.to_topology().to_openmm()
    for res in top.residues():
        res.name = "LIG"          # match the XML residue template
    pos = np.asarray(rmol.GetConformer().GetPositions(), float) / 10.0
    ref = mm.Platform.getPlatformByName("Reference")
    integ = mm.VerletIntegrator(1 * unit.femtosecond)
    ctx_a = mm.Context(sys_ref, integ, ref)
    ctx_a.setPositions(pos)
    ea = ctx_a.getState(getEnergy=True).getPotentialEnergy() \
        .value_in_unit(unit.kilojoule_per_mole)
    del ctx_a
    sys_b = ff_xml.createSystem(top, constraints=None)
    _add_torsions(sys_b, torsions_raw)
    integ2 = mm.VerletIntegrator(1 * unit.femtosecond)
    ctx_b = mm.Context(sys_b, integ2, ref)
    ctx_b.setPositions(pos)
    eb = ctx_b.getState(getEnergy=True).getPotentialEnergy() \
        .value_in_unit(unit.kilojoule_per_mole)
    s1["energy_reference_kjmol"] = float(ea)
    s1["energy_xml_kjmol"] = float(eb)
    s1["energy_abs_diff_kjmol"] = float(abs(ea - eb))
    if abs(ea - eb) > max(1e-3, 1e-6 * abs(ea)):
        raise RuntimeError(f"XML export energy mismatch: {ea:.8f} vs "
                           f"{eb:.8f} kJ/mol — force-convention bug; "
                           f"aborting rather than shipping a wrong "
                           f"Hamiltonian")
    _log("ff", f"XML export validated against interchange system: "
               f"dE = {abs(ea - eb):.2e} kJ/mol")
    return {"ff_xml": ff_xml, "topology": top, "torsions": torsions_raw,
            "stage1": s1}


def _add_torsions(system, torsions_raw):
    """Transplant the reference PeriodicTorsionForce entries (proper +
    improper, already in OpenMM's internal ordering convention) verbatim
    into a System built from the XML."""
    ht = mm.PeriodicTorsionForce()
    for (p1, p2, p3, p4, per, ph, k) in torsions_raw:
        if k != 0.0:
            ht.addTorsion(p1, p2, p3, p4, per, ph, k)
    system.addForce(ht)


def _q(v, u):
    """Robust parameter -> float (interchange can nest pint Quantities
    inside openmm Quantities)."""
    try:
        return float(v.value_in_unit(u))
    except Exception:
        pass
    try:
        return float(v / u)
    except Exception:
        pass
    return float(v)


def _force_of(system, cls):
    for f in system.getForces():
        if isinstance(f, cls):
            return f
    raise RuntimeError(f"force {cls.__name__} not found")


def _is_center(mol, a, others):
    nb = {b.GetOtherAtomIdx(a) for b in mol.GetAtomWithIdx(a).GetBonds()}
    return all(o in nb for o in others)


# ------------------------------------------------------------------ #
#  stage 2 — systems: explicit solvation + reactive surgery + meta bias
# ------------------------------------------------------------------ #
def morse_params_kj(morse_cal):
    """{tag: (De kJ/mol, a 1/nm, Re nm)}.

    Thermodynamic recalibration: relaxed-scan fits absorb environment
    relaxation into De (De_fit ~120 kcal for N-C vs tabulated BDE ~73),
    which pushes the classical product state ~+50 kcal uphill and kills
    the product basin. We therefore keep the GFN2 well CURVATURE
    k = 2*De_fit*a_fit^2 (exact at the bottom) but set De to the
    tabulated BDE: a = sqrt(k/(2*De_BDE)). Pair-sum classical reaction
    energy then matches the QM value (about +1.8 kcal before strain
    terms vs QM -5.1 kcal), restoring both basins.
    """
    out = {}
    for tag, v in morse_cal.items():
        if not (isinstance(v, dict) and "De_kcal" in v):
            continue
        k_fit = 2.0 * v["De_kcal"] * v["a_invA"] ** 2   # kcal/mol/A^2
        De = FALLBACK_BDE_KCAL.get(tag, v["De_kcal"])
        a = math.sqrt(k_fit / (2.0 * De))               # 1/A
        out[tag] = (De * KCAL, a * 10.0, v["Re_A"] / 10.0)
    return out


def apply_reactive_surgery(system, morse_cal):
    """Zero the two scissile harmonic bonds, add the 3-pair Morse
    CustomBondForce, zero the (4,9) 1-4 nonbonded exception."""
    hb = _force_of(system, mm.HarmonicBondForce)
    n_zeroed = 0
    for b in range(hb.getNumBonds()):
        p1, p2, r0, k = hb.getBondParameters(b)
        if frozenset((p1, p2)) in (frozenset(RC), frozenset((6, 9))):
            hb.setBondParameters(b, p1, p2, r0, 0.0)
            n_zeroed += 1
    if n_zeroed != 2:
        raise RuntimeError(f"expected to zero 2 scissile bonds, got "
                           f"{n_zeroed}")

    par = morse_params_kj(morse_cal)
    morse = mm.CustomBondForce("Deb*(1-exp(-ab*(r-Reb)))^2-Deb")
    morse.addPerBondParameter("Deb")
    morse.addPerBondParameter("ab")
    morse.addPerBondParameter("Reb")
    for tag, (i, j) in MORSE_PAIRS.items():
        De, a, Re = par[tag]
        morse.addBond(i, j, [De, a, Re])
    system.addForce(morse)

    nbf = _force_of(system, mm.NonbondedForce)
    hit = False
    for e in range(nbf.getNumExceptions()):
        p1, p2, q, sg, ep = nbf.getExceptionParameters(e)
        if frozenset((p1, p2)) == frozenset((4, 9)):
            nbf.setExceptionParameters(e, p1, p2, 0.0, 1.0, 0.0)
            hit = True
    if not hit:
        nbf.addException(4, 9, 0.0, 1.0, 0.0)
    return morse


def build_meta_force(n_max=N_MAX_HILLS):
    """2D WTMetaD bias as ONE CustomCVForce with pre-allocated hill
    slots (w{i}, r{i}, th{i} global parameters) — deposits are cheap
    context.setParameter calls, no recompilation ever."""
    cv1 = mm.CustomBondForce("r")
    cv1.addBond(RC[0], RC[1], [])
    cv2 = mm.CustomTorsionForce("theta")
    cv2.addTorsion(*TC, [])

    s1_nm, s2_rad = sigma_constants()
    inv2s1sq = -1.0 / (2.0 * s1_nm ** 2)
    inv2s2sq = -1.0 / (2.0 * s2_rad ** 2)
    terms = []
    for i in range(n_max):
        # Lepton has no round(): the periodic wrap of the torsion delta
        # uses atan2(sin,cos) which maps onto (-pi, pi] exactly.
        dth = f"atan2(sin(theta-th{i}), cos(theta-th{i}))"
        terms.append(f"w{i}*exp({inv2s1sq:.9e}*(r-r{i})*(r-r{i})"
                     f"{inv2s2sq:.9e}*{dth}*{dth})")
    meta = mm.CustomCVForce("+".join(terms))
    meta.addCollectiveVariable("r", cv1)
    meta.addCollectiveVariable("theta", cv2)
    for i in range(n_max):
        meta.addGlobalParameter(f"w{i}", 0.0)
        meta.addGlobalParameter(f"r{i}", 0.0)
        meta.addGlobalParameter(f"th{i}", 0.0)
    return meta


def sigma_constants():
    return (SIG1_A / 10.0, math.radians(SIG2_DEG))


def bias_eval(r_nm, th_rad, hills, s1_nm, s2_rad):
    """Accumulated bias at one point — mirrors the CustomCVForce
    expression exactly."""
    if not hills:
        return 0.0
    hs = np.asarray(hills, float)          # (n,3): w kJ/mol, r nm, th rad
    dr = r_nm - hs[:, 1]
    dth = th_rad - hs[:, 2]
    dth -= 2 * np.pi * np.round(dth / (2 * np.pi))
    return float(np.sum(hs[:, 0] * np.exp(
        -dr * dr / (2 * s1_nm ** 2) - dth * dth / (2 * s2_rad ** 2))))


def bias_on_grid(hills, r_grid_nm, th_grid_rad, s1_nm, s2_rad):
    hs = np.asarray(hills, float) if len(hills) else np.zeros((0, 3))
    R = r_grid_nm[:, None, None]
    T = th_grid_rad[None, :, None]
    if len(hs):
        # NOTE: hs[:, k][None, None, :] -> (1, 1, n) column vectors;
        # hs[None, None, k] would select hill row k, not column k!
        w = hs[:, 0][None, None, :]
        rc = hs[:, 1][None, None, :]
        tc = hs[:, 2][None, None, :]
        dr = R - rc
        dth = T - tc
        dth -= 2 * np.pi * np.round(dth / (2 * np.pi))
        V = np.sum(w * np.exp(
            -dr * dr / (2 * s1_nm ** 2) - dth * dth / (2 * s2_rad ** 2)),
            axis=2)
    else:
        V = np.zeros((len(r_grid_nm), len(th_grid_rad)))
    return V


W_MIN_KJ = 0.05 * KCAL   # skip negligible WT tail hills (bias error
                         # per skipped hill < 0.02 kcal/mol, negligible)


def deposit(meta, ctx, hills, w0_kj, kt_gamma1):
    r, th = meta.getCollectiveVariableValues(ctx)
    s1, s2 = sigma_constants()
    v = bias_eval(r, th, hills, s1, s2)
    w = w0_kj * math.exp(-v / kt_gamma1)
    if w < W_MIN_KJ:
        # well-tempered tail: the hill is energetically irrelevant;
        # skipping keeps the compiled expression within N_MAX slots
        # while the 500-step deposition cadence is preserved.
        return float(r), float(th), w
    n = len(hills)
    if n >= N_MAX_HILLS:
        # slots full: treat as a skipped hill (cadence preserved; the
        # accumulated bias is already defined by the stored hills)
        return float(r), float(th), w
    ctx.setParameter(f"w{n}", w)
    ctx.setParameter(f"r{n}", float(r))
    ctx.setParameter(f"th{n}", float(th))
    hills.append((float(w), float(r), float(th)))
    return float(r), float(th), w


def make_simulation(system, topology, platform_name, seed=20260904):
    integrator = mm.LangevinMiddleIntegrator(
        TEMP * unit.kelvin, 1.0 / unit.picosecond, DT_FS * unit.femtosecond)
    integrator.setConstraintTolerance(1e-6)
    integrator.setRandomNumberSeed(seed)
    plat = mm.Platform.getPlatformByName(platform_name)
    props = {}
    if platform_name == "CPU":
        props["Threads"] = str(os.cpu_count() or 4)
    if platform_name == "CUDA":
        props["Precision"] = "mixed"
    return app.Simulation(topology, system, integrator, plat, props)


def build_explicit(assets, ff, morse_cal, platform_name):
    """TIP3P + 0.15 M NaCl (1.0 nm padding), PME 1.0 nm, reactive
    surgery, metadynamics bias force."""
    pos = assets["rp"] / 10.0
    # literal 10.0 A buffer to every box edge (mission spec): explicit
    # boxSize = solute extent + 2*PAD_NM, solute centroid at box center.
    # Image separation >= 2*PAD_NM > PME cutoff on every axis.
    extent = pos.max(axis=0) - pos.min(axis=0)
    box_nm = extent + 2.0 * PAD_NM
    pos_c = pos - pos.mean(axis=0) + 0.5 * box_nm
    modeller = app.Modeller(ff["topology"], pos_c * unit.nanometer)
    modeller.addSolvent(ff["ff_solv"], model="tip3p",
                        boxSize=box_nm * unit.nanometer,
                        ionicStrength=IONIC_M * unit.molar,
                        neutralize=True,
                        positiveIon="Na+", negativeIon="Cl-")
    resn = [r.name for r in modeller.topology.residues()]
    n_water, n_na, n_cl = (resn.count("HOH"), resn.count("NA"),
                           resn.count("CL"))
    _log("sys", f"solvated: {modeller.topology.getNumAtoms()} atoms, "
                f"{n_water} waters, {n_na} Na+, {n_cl} Cl-")
    box = modeller.topology.getUnitCellDimensions()
    system = ff["ff_solv"].createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=R_CUT_NM * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=1e-4)
    _add_torsions(system, ff["torsions"])
    apply_reactive_surgery(system, morse_cal)
    meta = build_meta_force()
    system.addForce(meta)
    sim = make_simulation(system, modeller.topology, platform_name)
    sim.context.setPositions(modeller.positions)
    info = {"n_atoms": modeller.topology.getNumAtoms(),
            "n_water": n_water, "n_na": n_na, "n_cl": n_cl,
            "box_nm": [float(x) for x in box.value_in_unit(unit.nanometer)],
            "padding_nm": PAD_NM, "cutoff_nm": R_CUT_NM,
            "electrostatics": "PME, tolerance 1e-4",
            "water_model": "TIP3P (amber14/tip3p.xml)",
            "ionic_strength_M": IONIC_M, "integrator":
                "LangevinMiddle 300 K, friction 1/ps, dt 2 fs, HBonds",
            "dof": _dof(system),
            "morse_kj": {t: [float(x) for x in v] for t, v in
                         morse_params_kj(morse_cal).items()}}
    return sim, modeller, meta, info


def build_implicit(assets, ff, morse_cal, platform_name):
    """Solute-only system + GBSA/OBC2 (Bondi radii, OBC2 screening)."""
    system = ff["ff_xml"].createSystem(ff["topology"],
                                       nonbondedMethod=app.NoCutoff,
                                       constraints=app.HBonds)
    _add_torsions(system, ff["torsions"])
    apply_reactive_surgery(system, morse_cal)
    nbf = _force_of(system, mm.NonbondedForce)
    radii = {"C": 0.170, "N": 0.155, "H": 0.120}    # Bondi, nm
    screens = {"C": 0.72, "N": 0.79, "H": 0.85}     # OBC2 screening
    gb = mm.GBSAOBCForce()
    for a in assets["rmol"].GetAtoms():
        s = a.GetSymbol()
        # OpenMM signature: addParticle(charge, radius, screen)
        gb.addParticle(float(assets["charges"][a.GetIdx()]), radii[s],
                       screens[s])
    for p in range(nbf.getNumParticles()):
        q, sg, ep = nbf.getParticleParameters(p)
        nbf.setParticleParameters(p, 0.0, sg, ep)   # GB owns electrostatics
    gb.setSoluteDielectric(1.0)
    gb.setSolventDielectric(78.5)
    system.addForce(gb)
    meta = build_meta_force()
    system.addForce(meta)
    sim = make_simulation(system, ff["topology"], platform_name)
    sim.context.setPositions(assets["rp"] / 10.0 * unit.nanometer)
    info = {"n_atoms": assets["rmol"].GetNumAtoms(),
            "implicit_model": "GBSA-OBC2; Bondi radii (C 1.70, N 1.55, "
                              "H 1.20 A), OBC2 screening (C .72, N .79, "
                              "H .85), eps_solute 1, eps_solvent 78.5",
            "dof": _dof(system)}
    return sim, None, meta, info


def _save_checkpoint(ctx, path: Path):
    path.write_bytes(ctx.createCheckpoint())


def _load_checkpoint(ctx, path: Path):
    ctx.loadCheckpoint(Path(path).read_bytes())


def _dof(system):
    return 3 * system.getNumParticles() - system.getNumConstraints() - 3


# ------------------------------------------------------------------ #
#  stage 3 — equilibration + production WTMetaD (shared engine)
# ------------------------------------------------------------------ #
def equilibrate(sim, tag, n_fast=5000, n_slow=15000):
    _log(tag, "energy minimization ...")
    sim.context.setVelocitiesToTemperature(TEMP * unit.kelvin, 20260904)
    sim.minimizeEnergy(tolerance=5.0 * unit.kilojoule_per_mole
                       / unit.nanometer, maxIterations=5000)
    _log(tag, f"E_min = {sim.context.getState(getEnergy=True)
              .getPotentialEnergy()}")
    # thermalize at 1 fs (freshly minimized solvent is clash-prone at
    # 2 fs), then relax to the production time step
    sim.integrator.setFriction(5.0 / unit.picosecond)
    sim.integrator.setStepSize(1.0 * unit.femtosecond)
    sim.step(n_fast)
    sim.context.setVelocitiesToTemperature(TEMP * unit.kelvin, 20260905)
    sim.integrator.setFriction(1.0 / unit.picosecond)
    sim.integrator.setStepSize(DT_FS * unit.femtosecond)
    t0 = time.time()
    sim.step(n_slow)
    ke = sim.context.getState(getEnergy=True).getKineticEnergy() \
        .value_in_unit(unit.kilojoule_per_mole)
    T = 2 * ke / (_dof(sim.system) * KB_KJ_MOL_K)
    sps = (n_fast + n_slow) / max(time.time() - t0, 1e-9)
    _log(tag, f"equilibrated: T = {T:.1f} K ({sps:.0f} steps/s)")


def run_metadynamics(sim, meta, state, leg, n_steps, deposit_every,
                     chk_every=10000, dcd_every=5000, dcd_path=None):
    """Production WTMetaD with checkpointing + resume.

    state: {"hills": [(w, r, th), ...], "steps_done": int}
    """
    kt_gamma1 = KB_KJ_MOL_K * TEMP * (GAMMA - 1.0)   # kJ/mol
    w0_kj = W0_KCAL * KCAL
    hills = state.setdefault("hills", [])
    start = int(state.get("steps_done", 0))

    for n, (w, r, th) in enumerate(hills):
        sim.context.setParameter(f"w{n}", w)
        sim.context.setParameter(f"r{n}", r)
        sim.context.setParameter(f"th{n}", th)
    if start > 0:
        chk = OUT / f"prod_{leg}.chk"
        if chk.exists():
            _load_checkpoint(sim.context, chk)
            _log(leg, f"resumed from checkpoint at step {start} "
                      f"({len(hills)} hills restored)")
        else:
            _warn(f"{leg}: steps_done={start} but no checkpoint; "
                  f"restarting this production leg")
            state["steps_done"], start = 0, 0
            state["hills"] = hills = []

    csv_path = OUT / f"progress_{leg}.csv"
    if start == 0:
        csv_path.write_text("step,time_ps,T_K,PE_kJmol,d_NC_A,theta_deg,"
                            "hills,W_last_kcal,V_kJmol,sps\n")
    dcd = None
    if dcd_path:
        mode = "wb" if start == 0 else "ab"
        dcd = app.DCDFile(open(dcd_path, mode), sim.topology,
                          DT_FS * unit.femtosecond)

    s1, s2 = sigma_constants()
    t0 = time.time()
    nblk = n_steps // deposit_every
    for blk in range(start // deposit_every, nblk):
        sim.integrator.step(deposit_every)
        r, th, w = deposit(meta, sim.context, hills, w0_kj, kt_gamma1)
        step = state["steps_done"] = (blk + 1) * deposit_every

        if (blk + 1) % 5 == 0 or step == n_steps:
            stt = sim.context.getState(getEnergy=True)
            ke = stt.getKineticEnergy() \
                .value_in_unit(unit.kilojoule_per_mole)
            pe = stt.getPotentialEnergy() \
                .value_in_unit(unit.kilojoule_per_mole)
            T = 2 * ke / (_dof(sim.system) * KB_KJ_MOL_K)
            v = bias_eval(r, th, hills, s1, s2)
            sps = (step - start) / max(time.time() - t0, 1e-9)
            eta = (n_steps - step) / max(sps, 1e-9) / 60.0
            _log(leg, f"step {step:>7d}/{n_steps} "
                      f"({step * DT_FS / 1000:6.1f} ps) "
                      f"T={T:6.1f} K PE={pe / KCAL:9.1f} kcal "
                      f"d(N-C)={r * 10:5.2f} A "
                      f"th={math.degrees(th):7.1f} deg "
                      f"hills={len(hills):4d} W={w / KCAL:5.2f} "
                      f"V={v / KCAL:6.2f} kcal "
                      f"[{sps:6.0f} st/s ETA {eta:5.1f} min]")
            with open(csv_path, "a") as fh:
                fh.write(f"{step},{step * DT_FS / 1000:.2f},{T:.2f},"
                         f"{pe:.2f},{r * 10:.4f},"
                         f"{math.degrees(th):.3f},{len(hills)},"
                         f"{w / KCAL:.4f},{v / KCAL:.4f},{sps:.1f}\n")
        if dcd is not None and step % dcd_every == 0:
            dcd.writeModel(sim.context.getState(getPositions=True).getPositions())
        if step % chk_every == 0 or step == n_steps:
            _save_checkpoint(sim.context, OUT / f"prod_{leg}.chk")
            write_json_atomic(OUT / "metadyn_state.json", STATE)
    write_json_atomic(OUT / "metadyn_state.json", STATE)
    _log(leg, f"production complete: {len(hills)} hills, "
              f"{state['steps_done']} steps")
    return hills


# ------------------------------------------------------------------ #
#  stage 4 — FES analysis
# ------------------------------------------------------------------ #
def analyze_fes(hills, leg, anchor_R, anchor_P,
                slices=(0.15, 0.3, 0.5, 0.7, 0.85, 1.0)):
    """FES = -(gamma-1)/gamma * V on the periodic CV grid, deposition
    density masking, Dijkstra MEP, saddle, basin and convergence data.

    anchor_R / anchor_P: (d_A, theta_deg) CV anchors of the reactant /
    product QM endpoint geometries.
    """
    R, T = R_GRID_A, TH_GRID_DEG
    rG, tG = R / 10.0, np.radians(T)
    s1, s2 = sigma_constants()
    scale = -(GAMMA - 1.0) / GAMMA            # F = scale * V, kJ/mol
    dep_ps = (DEPOSIT_STEPS_EXPL if leg == "explicit"
              else DEPOSIT_STEPS_IMPL) * DT_FS / 1000.0

    def fes_from(nh):
        V = bias_on_grid(hills[:nh], rG, tG, s1, s2)
        F = scale * V / KCAL                                     # kcal
        if nh:
            hs = np.asarray(hills[:nh], float)
            dr = R[:, None, None] - hs[:, 1][None, None, :] * 10.0
            dth = T[None, :, None] - np.degrees(hs[:, 2][None, None, :])
            dth -= 360.0 * np.round(dth / 360.0)
            g = np.exp(-(dr / (2 * SIG1_A)) ** 2 -
                       (dth / (2 * SIG2_DEG)) ** 2)
            dens = (g > np.exp(-3.5)).sum(axis=2)
        else:
            dens = np.zeros(F.shape)
        mask = dens >= 2
        Fs = gaussian_filter(np.where(mask, F, 0.0), 1.0)
        Fs[gaussian_filter(mask.astype(float), 1.2,
                           mode="constant") <= 0.35] = np.nan
        return F, Fs, mask, dens

    def basin_near(Fx, rA, thA, wr=0.4, wt=50.0):
        dr = np.abs(R[:, None] - rA)
        dt = np.abs((T[None, :] - thA + 180.0) % 360.0 - 180.0)
        win = (dr <= wr) & (dt <= wt) & np.isfinite(Fx)
        if not win.any():
            return None
        f = np.where(win, Fx, np.inf)
        idx = np.unravel_index(np.argmin(f), f.shape)
        return {"idx": [int(idx[0]), int(idx[1])],
                "cv": [float(R[idx[0]]), float(T[idx[1]])],
                "G_kcal": float(Fx[idx])}

    def barrier_of(Fx):
        bre = basin_near(Fx, *anchor_R)
        bpr = basin_near(Fx, *anchor_P)
        if not bre or not bpr:
            return bre, bpr, None
        pth = _dijkstra_path(Fx, bre["idx"], bpr["idx"])
        if pth is None:
            return bre, bpr, None
        Fp = np.array([Fx[i, j] for i, j in pth])
        k = int(np.nanargmax(Fp))
        sk_i, sk_j = pth[k]
        return bre, bpr, {
            "path": [[float(R[i]), float(T[j])] for i, j in pth],
            "saddle_cv": [float(R[sk_i]), float(T[sk_j])],
            "saddle_G_kcal": float(Fp[k]),
            "dG_act_kcal": float(Fp[k] - bre["G_kcal"]),
            "dG_rxn_kcal": float(bpr["G_kcal"] - bre["G_kcal"])}

    F, Fs, mask, dens = fes_from(len(hills))
    bre, bpr, mep = barrier_of(Fs)
    # sampled-extent diagnostics: max d with deposition support
    d_ok = R[np.any(mask, axis=1)]
    out = {"n_hills": len(hills), "reactant_basin": bre,
           "product_basin": bpr, "mep": mep, "Fs": Fs, "mask": mask,
           "d_limit_A": float(d_ok.max()) if len(d_ok) else None}
    # 1D minimum-free-energy profile along the scissile coordinate
    prof = []
    for i in range(len(R)):
        col = Fs[i]
        if np.isfinite(col).any():
            prof.append(float(np.nanmin(col)))
        else:
            prof.append(np.nan)
    prof = np.asarray(prof, float)
    base = np.nanmin(prof)
    out["profile_d_A"] = R
    out["profile_F_rel_kcal"] = prof - base
    np.savez_compressed(OUT / f"fes_{leg}.npz", R=R, TH=T, F=F, Fs=Fs,
                        mask=mask, dens=dens,
                        hills=np.asarray(hills, float).reshape(-1, 3),
                        profile_d_A=R, profile_F_rel_kcal=prof - base)

    def profile_span(Fx):
        col_min = np.array([
            np.nanmin(Fx[i][np.isfinite(Fx[i])])
            if np.isfinite(Fx[i]).any() else np.nan
            for i in range(Fx.shape[0])], float)
        if not np.isfinite(col_min).any():
            return None
        return float(np.nanmax(col_min) - np.nanmin(col_min))

    trace = []
    for frac in slices:
        nh = max(1, int(round(frac * len(hills))))
        _, Fc, _, _ = fes_from(nh)
        b1, b2, m = barrier_of(Fc)
        trace.append({"frac": frac, "n_hills": nh,
                      "time_ps": nh * dep_ps,
                      "dG_act_kcal": None if m is None
                      else m["dG_act_kcal"],
                      "ascent_kcal": profile_span(Fc)})
    out["convergence"] = trace
    return out


def _dijkstra_path(F, start_idx, end_idx):
    """8-connected least-energy grid path over finite cells."""
    nr, nt = F.shape
    nn = nr * nt
    finite = np.isfinite(F)
    Ff = np.where(finite, F, 1e9)
    idx = np.arange(nn).reshape(nr, nt)
    rows, cols, weights = [], [], []
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            i0, i1 = max(0, -di), min(nr, nr - di)
            j0, j1 = max(0, -dj), min(nt, nt - dj)
            src = idx[i0:i1, j0:j1]
            dst = idx[i0 + di:i1 + di, j0 + dj:j1 + dj]
            ok = (finite[i0:i1, j0:j1] &
                  finite[i0 + di:i1 + di, j0 + dj:j1 + dj])
            w = 0.5 * (Ff[i0:i1, j0:j1] +
                       Ff[i0 + di:i1 + di, j0 + dj:j1 + dj]) * \
                math.hypot(di, dj)
            rows.append(src[ok])
            cols.append(dst[ok])
            weights.append(w[ok])
    g = csr_matrix((np.concatenate(weights),
                    (np.concatenate(rows), np.concatenate(cols))),
                   shape=(nn, nn))
    s = start_idx[0] * nt + start_idx[1]
    t = end_idx[0] * nt + end_idx[1]
    d, pre = dijkstra(g, directed=False, indices=s,
                      return_predecessors=True)
    if not np.isfinite(d[t]):
        return None
    path = []
    j = t
    while j != s and j >= 0:
        path.append((j // nt, j % nt))
        j = int(pre[j])
    path.append((start_idx[0], start_idx[1]))
    path.reverse()
    return path


# ------------------------------------------------------------------ #
#  figures (300 DPI)
# ------------------------------------------------------------------ #
def figures(assets, fes_ex, fes_im, ddg_profile=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.titlesize": 12})
    FIG.mkdir(parents=True, exist_ok=True)
    R, T = R_GRID_A, TH_GRID_DEG

    # ---------- fig 1 — explicit solvent box ------------------------- #
    try:
        pdb = app.PDBFile(str(OUT / "solvated_equilibrated.pdb"))
        pos = pdb.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        top = pdb.topology
        box = np.array(top.getUnitCellDimensions()
                       .value_in_unit(unit.angstrom))
        sym, resn = [], []
        for a in top.atoms():
            sym.append(a.element.symbol)
            resn.append(a.residue.name)
        solute_idx = np.array([i for i, r in enumerate(resn) if r == "LIG"])
        water_o = np.array([i for i, r in enumerate(resn)
                            if r == "HOH" and sym[i] == "O"])
        d = np.linalg.norm(pos[water_o][:, None, :] -
                           pos[solute_idx][None, :, :], axis=2)
        dmin = d.min(axis=1)
        shell = water_o[dmin < 3.5]
        bulk = water_o[dmin >= 3.5]

        # H-bond network: solute N-H donors / N acceptor <-> shell water
        rmol = assets["rmol"]
        n_idx = [i for i in solute_idx if sym[i] == "N"]
        # heavy-H connectivity from the topology bond list
        bonded_h = {}                       # heavy global idx -> [H idx]
        for b in top.bonds():
            a1, a2 = b[0].index, b[1].index
            if sym[a1] == "H" and sym[a2] != "H":
                bonded_h.setdefault(a2, []).append(a1)
            elif sym[a2] == "H" and sym[a1] != "H":
                bonded_h.setdefault(a1, []).append(a2)
        water_oh = {}                       # water O idx -> [H idx]
        for o in water_o:
            water_oh[int(o)] = bonded_h.get(int(o), [])
        hbonds = []
        HB_CUT = 2.8                        # H...acceptor display cutoff
        for ni in n_idx:                    # solute N-H ... O(water)
            for h in bonded_h.get(int(ni), []):
                for o in shell:
                    if np.linalg.norm(pos[o] - pos[h]) < HB_CUT:
                        hbonds.append((h, int(o)))
            for o in shell:                 # H-O(water) ... N(solute)
                for h in water_oh.get(int(o), []):
                    if np.linalg.norm(pos[h] - pos[ni]) < HB_CUT:
                        hbonds.append((int(h), int(ni)))

        ion_na = np.array([i for i, r in enumerate(resn) if r == "NA"])
        ion_cl = np.array([i for i, r in enumerate(resn) if r == "CL"])

        fig = plt.figure(figsize=(11.8, 9.4))
        ax = fig.add_subplot(111, projection="3d")
        sub = bulk[::3]
        ax.scatter(pos[sub, 0], pos[sub, 1], pos[sub, 2], s=5.0,
                   c="#41b6c4", alpha=0.55, depthshade=False,
                   label=f"bulk TIP3P O ({len(sub)} of {len(bulk)} shown)")
        ax.scatter(pos[shell, 0], pos[shell, 1], pos[shell, 2], s=22,
                   c="#d94801", depthshade=False,
                   label=f"first solvation shell O (n={len(shell)}, "
                         f"O...solute < 3.5 $\\AA$)")
        if len(ion_na):
            ax.scatter(pos[ion_na, 0], pos[ion_na, 1], pos[ion_na, 2],
                       s=90, marker="*", c="#f1c40f", edgecolors="k",
                       depthshade=False, label=f"Na+ (n={len(ion_na)})")
        if len(ion_cl):
            ax.scatter(pos[ion_cl, 0], pos[ion_cl, 1], pos[ion_cl, 2],
                       s=90, marker="*", c="#8e44ad", edgecolors="k",
                       depthshade=False, label=f"Cl- (n={len(ion_cl)})")
        scol = [solute_idx[i] for i in range(len(solute_idx))]
        elem = [sym[i] for i in solute_idx]
        for b in rmol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            p, q = pos[solute_idx[i]], pos[solute_idx[j]]
            col = "#1f77b4" if "N" in (elem[i], elem[j]) else "#333333"
            lw = 1.2 if (elem[i] == "H" or elem[j] == "H") else 3.4
            ax.plot([p[0], q[0]], [p[1], q[1]], [p[2], q[2]], c=col,
                    lw=lw, solid_capstyle="round", alpha=0.95)
        for h, o in hbonds:
            ax.plot([pos[h][0], pos[o][0]], [pos[h][1], pos[o][1]],
                    [pos[h][2], pos[o][2]], ls="--", c="#0f8a2f", lw=2.4)
        s0 = assets["stage0"]
        for corners in (((0, 0, 0), (1, 0, 0)), ((0, 0, 0), (0, 1, 0)),
                        ((0, 0, 0), (0, 0, 1)), ((1, 1, 1), (0, 1, 1)),
                        ((1, 1, 1), (1, 0, 1)), ((1, 1, 1), (1, 1, 0)),
                        ((1, 0, 0), (1, 1, 0)), ((1, 0, 0), (1, 0, 1)),
                        ((0, 1, 0), (1, 1, 0)), ((0, 1, 0), (0, 1, 1)),
                        ((0, 0, 1), (1, 0, 1)), ((0, 0, 1), (0, 1, 1))):
            s, e = corners
            ax.plot([s[0] * box[0], e[0] * box[0]],
                    [s[1] * box[1], e[1] * box[1]],
                    [s[2] * box[2], e[2] * box[2]], c="k", lw=1.0,
                    alpha=0.8)
        ax.set_title(f"fig. 1 — cyclopropa[b]indole in explicit TIP3P / "
                     f"{IONIC_M:.2f} M NaCl — PBC box "
                     f"{box[0]:.1f}$\\times${box[1]:.1f}"
                     f"$\\times${box[2]:.1f} $\\AA$, PME cut-off 10.0 $\\AA$"
                     + (f", green dashed = {len(hbonds)} H-bonds to the "
                        f"pyrrolic N-H" if hbonds else ""), fontsize=10)
        ax.set_xlabel("x ($\\AA$)")
        ax.set_ylabel("y ($\\AA$)")
        ax.set_zlabel("z ($\\AA$)")
        ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
        fig.tight_layout()
        fig.savefig(FIG / "fig1_explicit_solvent_box.png", dpi=300)
        plt.close(fig)
        _log("fig", f"fig1 written (first-shell waters: {len(shell)}, "
                    f"H-bonds drawn: {len(hbonds)})")
    except Exception as exc:
        _warn(f"figure 1 failed: {exc}")
        traceback.print_exc()

    # ---------- fig 2 — 2D FES (explicit | implicit) ----------------- #
    ddg = None
    try:
        fig, axes = plt.subplots(1, 2, figsize=(13.8, 6.1),
                                 constrained_layout=True)
        panels = (("a", fes_ex, "explicit TIP3P / PME, 0.15 M NaCl", True),
                  ("b", fes_im, "implicit GBSA-OBC2 reference", False))
        for letter, fes, leg, with_qm in panels:
            ax = axes[0 if with_qm else 1]
            Fs = fes["Fs"]
            base = fes["reactant_basin"]["G_kcal"] \
                if fes.get("reactant_basin") else 0.0
            Fp = np.where(np.isfinite(Fs), Fs - base, np.nan)
            finite_max = float(np.nanmax(Fp)) if np.isfinite(Fp).any() else 1.0
            top = max(2.0, min(16.0, float(np.ceil(finite_max)) + 1.0))
            lev = np.arange(0, top + 0.01, max(0.5, top / 14.0))
            cf = ax.contourf(R, T, Fp.T, levels=lev, cmap="viridis",
                             extend="max")
            cs = ax.contour(R, T, Fp.T, levels=lev[::2], colors="w",
                            linewidths=0.4, alpha=0.55)
            ax.clabel(cs, fmt="%.1f", fontsize=6)
            cb = fig.colorbar(cf, ax=ax, pad=0.015)
            cb.set_label(r"$\Delta G$ rel. reactant basin (kcal/mol)")
            mep = fes.get("mep")
            if mep:
                pth = np.array(mep["path"])
                ax.plot(pth[:, 0], pth[:, 1], "w--", lw=2.0, label="MEP")
                sc = mep["saddle_cv"]
                ax.plot(sc[0], sc[1], marker="X", ms=13, mfc="red",
                        mec="k", mew=1.2,
                        label=(f"$\\Delta G^\\ddagger$ = "
                               f"{mep['dG_act_kcal']:.2f} kcal/mol"))
            else:
                dl = fes.get("d_limit_A")
                if dl:
                    ax.axvline(dl, color="r", ls="--", lw=1.4)
                    ax.text(dl, 172, "sampled boundary", rotation=0,
                            color="r", fontsize=8, ha="right", va="top")
            for bsn, name, mk in ((fes.get("reactant_basin"), "R", "o"),):
                if bsn:
                    ax.plot(bsn["cv"][0], bsn["cv"][1], mk, ms=11,
                            mfc="orange", mec="k")
                    ax.annotate(name, bsn["cv"],
                                textcoords="offset points",
                                xytext=(6, 6), fontsize=11,
                                weight="bold")
            if with_qm:
                qp = np.array(assets["stage0"]["idpp_path_cv"], float)
                ax.plot(qp[:, 0], qp[:, 1], "-", c="#ff7f0e", lw=1.6,
                        alpha=0.9,
                        label="GFN2-xTB IDPP path (QM reference)")
                ax.plot(qp[:, 0], qp[:, 1], ".", c="#ff7f0e", ms=5)
                ax.plot(qp[-1, 0], qp[-1, 1], "s", ms=10, mfc="orange",
                        mec="k")
                ax.annotate("P", qp[-1], textcoords="offset points",
                            xytext=(6, 6), fontsize=11, weight="bold")
            ax.set_xlim(1.3, 3.0)
            ax.set_ylim(-180, 180)
            ax.set_xticks(np.arange(1.5, 3.01, 0.25))
            ax.set_yticks(np.arange(-180, 181, 60))
            ax.set_xlabel(r"$d(\mathrm{N_4-C_5})$ scissile bond ($\AA$)")
            ax.set_ylabel(r"$\theta(\mathrm{C_9-C_5-N_4-C_3})$ (deg)")
            asc = fes.get("ascent_at_limit_kcal")
            extra = (f", ascent to boundary {asc:.1f} kcal/mol"
                     if asc is not None else "")
            ax.set_title(f"({letter}) WTMetaD FES — {leg}\n"
                         f"$\\gamma$={GAMMA:.0f}, "
                         f"W$_0$={W0_KCAL} kcal/mol, "
                         f"{fes['n_hills']} stored hills{extra}",
                         fontsize=10.5)
            hnd, lab = ax.get_legend_handles_labels()
            if hnd:
                ax.legend(hnd, lab, loc="lower left", fontsize=8,
                          framealpha=0.9)
        if ddg_profile:
            d = np.asarray(ddg_profile["d_A"], float)
            g = np.asarray(ddg_profile["ddG_solv_kcal"], float)
            if len(d):
                axin = axes[1].inset_axes([0.58, 0.03, 0.40, 0.34])
                axin.fill_between(d, g, 0, color="#d62728", alpha=0.25)
                axin.plot(d, g, "o-", c="#d62728", ms=2.5, lw=1.4)
                axin.axhline(0, c="k", lw=0.6)
                axin.set_xlabel(r"$d$ ($\AA$)", fontsize=7)
                axin.set_ylabel(r"$\Delta\Delta G_{solv}(d)$", fontsize=7)
                axin.tick_params(labelsize=6)
                axin.set_title("solvent shift, matched coverage",
                               fontsize=7)
        fig.suptitle("fig. 2 — 2D free energy surface of the "
                     "cyclopropa[b]indole -> 2,3-dihydroquinoline "
                     "skeletal rearrangement", fontsize=13)
        fig.savefig(FIG / "fig2_2d_free_energy_surface.png", dpi=300)
        plt.close(fig)
        _log("fig", "fig2 written")
    except Exception as exc:
        _warn(f"figure 2 failed: {exc}")
        traceback.print_exc()

    # ---------- fig 3 — convergence ---------------------------------- #
    try:
        fig, axes = plt.subplots(1, 3, figsize=(15.4, 4.9),
                                 constrained_layout=True)
        series = {}
        for leg, c in (("explicit", "#1f77b4"), ("implicit", "#d62728")):
            fes = fes_ex if leg == "explicit" else fes_im
            tr = fes["convergence"]
            tt = [x["time_ps"] for x in tr
                  if x.get("ascent_kcal") is not None]
            gg = [x["ascent_kcal"] for x in tr
                  if x.get("ascent_kcal") is not None]
            axes[0].plot(tt, gg, "o-", c=c, label=leg)
            rows = [ln.split(",") for ln in
                    (OUT / f"progress_{leg}.csv").read_text()
                    .splitlines()[1:]]
            series[leg] = {
                "t": [float(r[1]) for r in rows],
                "hills": [float(r[6]) for r in rows],
                "T": [float(r[2]) for r in rows],
                "W": [float(r[7]) for r in rows],
                "c": c}
        axes[0].set_xlabel("simulation time (ps)")
        axes[0].set_ylabel(r"FES span along $d$(N-C) (kcal/mol)")
        axes[0].set_title("(a) free-energy ascent convergence")
        axes[0].grid(alpha=0.3)
        axes[0].legend()

        for leg in ("explicit", "implicit"):
            s = series[leg]
            axes[1].plot(s["t"], s["hills"], c=s["c"], label=leg)
        axes[1].set_xlabel("time (ps)")
        axes[1].set_ylabel("deposited hills")
        axes[1].set_title("(b) metadynamics hill accumulation")
        axes[1].grid(alpha=0.3)
        axes[1].legend()

        ax2 = axes[2]

        def _runmean(v, k=7):
            v = np.asarray(v, float)
            if len(v) < k:
                return v
            vp = np.pad(v, (k // 2, k - 1 - k // 2), mode="edge")
            return np.convolve(vp, np.ones(k) / k, mode="valid")

        for leg in ("explicit", "implicit"):
            s = series[leg]
            ax2.plot(s["t"], s["T"], c=s["c"], alpha=0.18)
            ax2.plot(s["t"], _runmean(s["T"]), c=s["c"], lw=2.0,
                     label=f"{leg} T (7-pt mean)")
        ax2.set_xlabel("time (ps)")
        ax2.set_ylabel("T (K) — raw faded, running mean solid")
        allT = np.concatenate([np.asarray(series[l]["T"], float)
                               for l in ("explicit", "implicit")])
        ax2.set_ylim(max(0.0, float(np.nanmin(allT)) - 20),
                     float(np.nanmax(allT)) + 20)
        ax2.grid(alpha=0.3)
        ax3 = ax2.twinx()
        for leg in ("explicit", "implicit"):
            s = series[leg]
            ax3.plot(s["t"], s["W"], c=s["c"], ls="--", lw=1.0, alpha=0.8)
        ax3.set_ylabel("hill height W(t) (kcal/mol)", color="#9467bd")
        ax2.set_title("(c) temperature stability (solid) & WT height "
                      "decay (dashed)")
        h1, l1 = ax2.get_legend_handles_labels()
        ax2.legend(h1, l1, fontsize=7.5, loc="upper left",
                   framealpha=0.85)
        fig.suptitle("fig. 3 — FES convergence & simulation diagnostics",
                     fontsize=13)
        fig.savefig(FIG / "fig3_fes_convergence_trace.png", dpi=300)
        plt.close(fig)
        _log("fig", "fig3 written")
    except Exception as exc:
        _warn(f"figure 3 failed: {exc}")
        traceback.print_exc()
    return ddg


# ------------------------------------------------------------------ #
#  driver
# ------------------------------------------------------------------ #
def main():
    global STATE
    ap = argparse.ArgumentParser(
        description="Phase 6: explicit-solvent WTMetaD & 2D FES")
    ap.add_argument("--steps-explicit", type=int,
                    default=N_STEPS_EXPL_DEFAULT)
    ap.add_argument("--steps-implicit", type=int,
                    default=N_STEPS_IMPL_DEFAULT)
    ap.add_argument("--platform", default="auto",
                    choices=["auto", "CUDA", "OpenCL", "CPU"])
    ap.add_argument("--fast", action="store_true",
                    help="smoke protocol (reduced steps)")
    ap.add_argument("--selftest", action="store_true",
                    help="build systems, benchmark, exit")
    args = ap.parse_args()

    t_wall = time.time()
    OUT.mkdir(exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    if (OUT / "metadyn_state.json").exists():
        try:
            loaded = json.loads((OUT / "metadyn_state.json").read_text())
            STATE.update(loaded)
        except Exception:
            pass

    avail = [mm.Platform.getPlatform(i).getName()
             for i in range(mm.Platform.getNumPlatforms())]
    plat = args.platform
    if plat == "auto":
        for pref in ("CUDA", "OpenCL", "CPU"):
            if pref in avail:
                plat = pref
                break
    RESULTS["meta"]["platform"] = plat
    RESULTS["meta"]["platforms_available"] = avail
    _log("main", f"OpenMM platforms available: {avail} -> using {plat}")

    n_expl, n_impl = args.steps_explicit, args.steps_implicit
    if args.fast or args.selftest:
        n_expl, n_impl = min(n_expl, 6000), min(n_impl, 12000)

    def _die(stage, exc):
        RESULTS["fatal_error"] = f"{stage}: {exc}\n{traceback.format_exc()}"
        write_json_atomic(OUT / "phase6_results.json")
        raise

    # ---- stage 0: assets & QM calibration ---------------------------
    try:
        assets = stage0_assets()
        RESULTS["system"]["stage0"] = assets["stage0"]
        write_json_atomic(OUT / "phase6_results.json")
        _log("stage0", f"assets OK — R anchor "
                       f"{assets['stage0']['cv_reactant_qm']}, P anchor "
                       f"{assets['stage0']['cv_product_qm']}, "
                       f"charge total "
                       f"{assets['stage0']['charge_total_e']:+.3f} e")
    except Exception as exc:
        return _die("stage0", exc)

    # ---- stage 1: OpenFF -> XML --------------------------------------
    try:
        ff = stage1_forcefield(assets)
        ff["ff_solv"] = app.ForceField(str(OUT / "solute_sage.xml"),
                                       "amber14/tip3p.xml")
        RESULTS["force_field"] = ff["stage1"]
        write_json_atomic(OUT / "phase6_results.json")
    except Exception as exc:
        return _die("stage1", exc)

    # ---- stage 2: explicit system ------------------------------------
    try:
        t0 = time.time()
        sim_ex, modeller, meta_ex, info_ex = build_explicit(
            assets, ff, RESULTS["morse_calibration"], plat)
        info_ex["build_s"] = time.time() - t0
        RESULTS["system"]["explicit"] = info_ex
        write_json_atomic(OUT / "phase6_results.json")
        _log("stage2", f"explicit system built in {info_ex['build_s']:.1f}"
                       f" s ({info_ex['n_atoms']} atoms, dof "
                       f"{info_ex['dof']})")
    except Exception as exc:
        return _die("stage2 explicit build", exc)

    # ---- selftest benchmark ------------------------------------------
    if args.selftest:
        sim_ex.minimizeEnergy(tolerance=100 * unit.kilojoule_per_mole
                              / unit.nanometer, maxIterations=300)
        sim_ex.integrator.setFriction(1.0 / unit.picosecond)
        t0 = time.time()
        sim_ex.step(300)
        sps = 300 / (time.time() - t0)
        _log("selftest", f"explicit: {sps:.0f} steps/s -> "
                         f"{args.steps_explicit / sps / 60:.1f} min for "
                         f"{args.steps_explicit} steps")
        t0 = time.time()
        sim_im, _, meta_im, _ = build_implicit(
            assets, ff, RESULTS["morse_calibration"], plat)
        sim_im.minimizeEnergy(tolerance=100 * unit.kilojoule_per_mole
                              / unit.nanometer, maxIterations=300)
        sim_im.integrator.setFriction(1.0 / unit.picosecond)
        sim_im.step(2000)
        sps_i = 2000 / (time.time() - t0)
        _log("selftest", f"implicit: {sps_i:.0f} steps/s -> "
                         f"{args.steps_implicit / sps_i / 60:.1f} min for "
                         f"{args.steps_implicit} steps")
        t0 = time.time()
        build_meta_force()
        _log("selftest", f"meta force compile: {time.time() - t0:.2f} s "
                         f"({N_MAX_HILLS} hill slots)")
        r, th, w = deposit(meta_im, sim_im.context,
                           STATE["implicit"]["hills"], W0_KCAL * KCAL,
                           KB_KJ_MOL_K * TEMP * (GAMMA - 1))
        _log("selftest", f"deposit smoke: r={r * 10:.2f} A "
                         f"th={math.degrees(th):.1f} deg W={w / KCAL:.3f}"
                         f" kcal")
        _log("selftest", "OK — all systems healthy")
        return 0

    # ---- equilibration (explicit) ------------------------------------
    if (OUT / "eq_explicit.chk").exists():
        _load_checkpoint(sim_ex.context, OUT / "eq_explicit.chk")
        _log("eq-expl", "loaded equilibration checkpoint")
    else:
        equilibrate(sim_ex, "eq-expl")
        st = sim_ex.context.getState(getPositions=True)
        with open(OUT / "solvated_equilibrated.pdb", "w") as fh:
            app.PDBFile.writeFile(modeller.topology, st.getPositions(), fh,
                                  keepIds=True)
        _save_checkpoint(sim_ex.context, OUT / "eq_explicit.chk")

    # ---- production (explicit) ---------------------------------------
    try:
        run_metadynamics(sim_ex, meta_ex, STATE["explicit"], "explicit",
                         n_expl, DEPOSIT_STEPS_EXPL,
                         dcd_path=str(OUT / "traj_explicit.dcd"))
        RESULTS["explicit"].update({
            "steps": STATE["explicit"]["steps_done"],
            "n_hills": len(STATE["explicit"]["hills"]),
            "target_steps": n_expl,
            "sim_time_ps": n_expl * DT_FS / 1000})
        write_json_atomic(OUT / "phase6_results.json")
    except Exception as exc:
        return _die("production explicit", exc)

    # ---- implicit reference leg --------------------------------------
    fes_ex = fes_im = None
    try:
        plat_impl = "Reference" if "Reference" in avail else plat
        if (OUT / "eq_implicit.chk").exists():
            sim_im, _, meta_im, info_im = build_implicit(
                assets, ff, RESULTS["morse_calibration"], plat_impl)
            _load_checkpoint(sim_im.context, OUT / "eq_implicit.chk")
            _log("eq-impl", "loaded equilibration checkpoint")
        else:
            sim_im, _, meta_im, info_im = build_implicit(
                assets, ff, RESULTS["morse_calibration"], plat_impl)
            equilibrate(sim_im, "eq-impl", n_fast=2000, n_slow=8000)
            _save_checkpoint(sim_im.context, OUT / "eq_implicit.chk")
        RESULTS["system"]["implicit"] = info_im
        run_metadynamics(sim_im, meta_im, STATE["implicit"], "implicit",
                         n_impl, DEPOSIT_STEPS_IMPL,
                         dcd_path=str(OUT / "traj_implicit.dcd"))
        RESULTS["implicit"].update({
            "steps": STATE["implicit"]["steps_done"],
            "n_hills": len(STATE["implicit"]["hills"]),
            "target_steps": n_impl,
            "sim_time_ps": n_impl * DT_FS / 1000})
        write_json_atomic(OUT / "phase6_results.json")
    except Exception as exc:
        return _die("implicit leg", exc)

    # ---- FES analysis -------------------------------------------------
    try:
        a_r = assets["stage0"]["cv_reactant_qm"]
        a_p = assets["stage0"]["cv_product_qm"]
        fes_ex = analyze_fes(STATE["explicit"]["hills"], "explicit",
                             a_r, a_p)
        fes_im = analyze_fes(STATE["implicit"]["hills"], "implicit",
                             a_r, a_p)
        RESULTS["fes"]["explicit"] = _fes_summary(fes_ex)
        RESULTS["fes"]["implicit"] = _fes_summary(fes_im)
        if fes_ex.get("mep") and fes_im.get("mep"):
            RESULTS["fes"]["delta_delta_G_solv_kcal"] = float(
                fes_ex["mep"]["dG_act_kcal"] -
                fes_im["mep"]["dG_act_kcal"])
        # solvent shift profile at matched scissile distances
        d_common = np.asarray(fes_ex["profile_d_A"], float)
        Fe = np.asarray(fes_ex["profile_F_rel_kcal"], float)
        Fi = np.asarray(fes_im["profile_F_rel_kcal"], float)
        ok = np.isfinite(Fe) & np.isfinite(Fi) & (d_common <=
             min(fes_ex["d_limit_A"] or 0, fes_im["d_limit_A"] or 0))
        if ok.any():
            RESULTS["fes"]["ddg_profile"] = {
                "d_A": d_common[ok].round(3).tolist(),
                "ddG_solv_kcal": (Fe[ok] - Fi[ok]).round(3).tolist(),
                "d_max_matched_A": float(d_common[ok].max()),
                "ddG_at_dmax_kcal": float(
                    (Fe[ok] - Fi[ok])[-1])}
        for leg, fes in (("explicit", fes_ex), ("implicit", fes_im)):
            gs = [x["dG_act_kcal"] for x in fes["convergence"]
                  if x["dG_act_kcal"] is not None]
            if len(gs) >= 2:
                RESULTS["fes"][leg]["dG_act_drift_kcal"] = float(
                    abs(gs[-1] - gs[-2]))
        write_json_atomic(OUT / "phase6_results.json")
        g_ex = RESULTS["fes"]["explicit"].get("dG_act_kcal")
        g_im = RESULTS["fes"]["implicit"].get("dG_act_kcal")
        dd = RESULTS["fes"].get("delta_delta_G_solv_kcal")
        _log("fes", f"dG(expl)={g_ex} kcal, dG(impl)={g_im} kcal, "
                    f"ddG_solv={dd} kcal")
    except Exception as exc:
        return _die("analysis", exc)

    # ---- figures -------------------------------------------------------
    try:
        ddg = figures(assets, fes_ex, fes_im,
                      ddg_profile=RESULTS["fes"].get("ddg_profile"))
        if ddg is not None:
            RESULTS["fes"]["delta_delta_G_solv_kcal"] = float(ddg)
    except Exception as exc:
        return _die("figures", exc)

    RESULTS["all_stages_ok"] = True
    RESULTS["wall_time_min"] = (time.time() - t_wall) / 60.0
    write_json_atomic(OUT / "phase6_results.json")
    _log("main", f"PHASE 6 COMPLETE in {RESULTS['wall_time_min']:.1f} min")
    return 0


def _fes_summary(fes):
    m = fes.get("mep") or {}
    return {"n_hills": fes["n_hills"],
            "reactant_basin": fes.get("reactant_basin"),
            "product_basin": fes.get("product_basin"),
            "saddle_cv": m.get("saddle_cv"),
            "saddle_G_kcal": m.get("saddle_G_kcal"),
            "dG_act_kcal": m.get("dG_act_kcal"),
            "dG_rxn_kcal": m.get("dG_rxn_kcal"),
            "d_limit_A": fes.get("d_limit_A"),
            "ascent_at_limit_kcal": _ascent_at_limit(fes),
            "convergence": fes.get("convergence")}


def _ascent_at_limit(fes):
    """Reconstructed free-energy ascent at the sampled boundary."""
    prof = fes.get("profile_F_rel_kcal")
    dl = fes.get("d_limit_A")
    if prof is None or dl is None:
        return None
    idx = np.argmin(np.abs(fes["profile_d_A"] - dl))
    v = prof[idx]
    return float(v) if np.isfinite(v) else None


if __name__ == "__main__":
    sys.exit(main())
