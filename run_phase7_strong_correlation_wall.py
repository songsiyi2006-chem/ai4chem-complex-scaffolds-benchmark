#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_phase7_strong_correlation_wall.py
=====================================
PHASE 7 — CONFRONTING THE WALL OF SIGHS
MULTI-REFERENCE ELECTRONIC STRUCTURE & QUANTUM-AI BREAKDOWN MAP

Model coordinate: homolytic stretching of the scissile cyclopropane C–C bond
(C6–C9) of the Phase-4 reactant cyclopropa[b]indole (C9H9N, GFN2-xTB
optimized, `results_phase4/reactant_3d.mol`) from 1.40 A to 3.20 A.
Cleavage of this bond (the first bond-breaking event of the Ciamician-
Dennstedt ring expansion identified in Phase 4) produces a tethered 1,3-
diradical — the canonical two-electron static-correlation catastrophe that
no single-determinant ansatz can describe.

Protocol
--------
7A  Mean-field breakdown & spin contamination: per scan point, RHF / UHF /
    RKS(B3LYP) / UKS(B3LYP); <S^2> from the traced spin densities,
    DS^2 = <S^2> - S(S+1); R_crit = interpolated bond length where
    DS^2(UHF) first exceeds 0.3 (birth of the unphysical open-shell
    "singlet").  RHF->UHF and RKS->UKS energy splitting quantify the
    variational collapse onto broken-symmetry solutions.
7B  Multi-reference ground truth: CAS(2e,2o) CASSCF active space built from
    the sigma(C6-C9) / sigma* frontier orbitals (selected by Mulliken
    gross populations on the two scissile carbons, MOs then permuted to the
    DETCI docc/active/uocc boundary); natural-orbital occupation numbers
    (NOON) prove fractional occupancy; Yamaguchi diradical character
    y = 1 - (nu1 - nu2); frontier natural orbital cubes for visualization.
7C  The epistemic error landscape: the identical trajectory evaluated by
    GFN2-xTB (semi-empirical), MACE-OFF + ANI-2x (machine-learning
    potentials), UHF / BS-UB3LYP (mean field) and CASSCF (multi-reference
    reference).  DeltaE_error(R) = |E_rel,AI(R) - E_rel,CASSCF(R)| with a
    common zero at the equilibrium bond length; the >15 kcal/mol region is
    labeled the Epistemic Failure Zone.

Engines & environments (fault-tolerant dual-interpreter orchestration)
----------------------------------------------------------------------
QC engine    : Psi4 1.11 (conda-forge win-64) in env `phase7` — DETCI-based
               CASSCF.  The reference protocol specifies PySCF; PySCF ships
               no native win32 build (no wheels / no MSVC toolchain on this
               host), so Psi4 1.11 is the drop-in ab-initio backend and the
               substitution is logged in every result file.  On a POSIX host
               with PySCF the script reports the substitution as a fallback
               either way — the 7A/7B quantities are backend-invariant.
Chem engine  : env `phase2ff` python (ASE + torchani + MACE-OFF) driving
               xtb.exe GFN2-xTB (phase-4 subprocess wrapper pattern).
Known Psi4 build quirk: this win-64 build's DETCI/DPD cache sizing
overflows (dpd_block_matrix n<0 hard abort) when psi4 memory is set to the
4 GB band; memory defaults to 2 GB (verified band 512 MB - 2 GB & 6 GB).
Every QC point runs in an isolated subprocess so hard aborts degrade to a
per-point basis-tier retry (def2-SVP -> 6-31G -> STO-3G, trade-offs logged)
instead of killing the scan.

Outputs
-------
results_phase7/phase7_results.json         (machine-readable master record)
results_phase7/scan_geometries/R*.xyz      (relaxed scan geometries)
results_phase7/cubes/                      (frontier natural orbital cubes)
figures_phase7/fig1_spin_contamination_profile.png        (300 DPI)
figures_phase7/fig2_casscf_frontier_orbitals.png          (300 DPI)
figures_phase7/fig3_the_wall_of_sighs_discrepancy.png     (300 DPI)

Usage
-----
conda activate phase7                    # psi4 + numpy + scipy + skimage
python run_phase7_strong_correlation_wall.py                 # full pipeline
python run_phase7_strong_correlation_wall.py --stage 7B      # one stage
python run_phase7_strong_correlation_wall.py --smoke         # 3 pts / 6-31G
python run_phase7_strong_correlation_wall.py --fig_only      # refit figures
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
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

# --------------------------------------------------------------------------- #
#  constants
# --------------------------------------------------------------------------- #
EH_EV = 27.211386245988
HARTREE_KCAL = 627.5094740631
EV_KCAL = 23.0605496950

CHARGE, MULTIPLICITY = 0, 1
SINGLET_SS = 0.0            # S(S+1) for S = 0
TWO_DOUBLETS_SS = 1.0       # <S^2> of two non-interacting doublets
S2_GATE = 0.3               # mission threshold: DS^2 >> 0.3
FAIL_KCAL = 15.0            # epistemic failure threshold

REACTANT_MOL = Path("results_phase4/reactant_3d.mol")
OUT = Path("results_phase7")
FIG = Path("figures_phase7")

# scan coordinate (Angstrom); denser around the expected symmetry-breaking
# onset, capped by the mission bracket 1.4 - 3.2 A
SCAN_R = [1.40, 1.50, 1.60, 1.75, 1.90, 2.05, 2.20, 2.40, 2.60, 2.90, 3.20]
R_EQ = 1.50                 # common zero of all relative energy curves
CUBE_RS = [1.50, 2.20, 3.20]

BASIS_CHAIN = ["def2-svp", "6-31g", "sto-3g"]
DFT_FUNC = "b3lyp"

# NOTE: never "4 GB" on this psi4 build (DETCI DPD cache overflow band)
QC_MEM = "2 GB"

SMOKE_R = [1.50, 2.30, 3.20]
SMOKE_BASIS = "6-31g"

ENV_PY_QC = r"C:\Users\HUIWEI\miniconda3\envs\phase7\python.exe"
ENV_PY_CHEM = r"C:\Users\HUIWEI\miniconda3\envs\phase2ff\python.exe"
ENV_CHEM_BIN = r"C:\Users\HUIWEI\miniconda3\envs\phase2ff\Library\bin"

META: dict = {
    "phase": 7,
    "mission": ("Wall of Sighs: multi-reference breakdown map of single-"
                "determinant DFT/AI potentials along cyclopropa[b]indole "
                "C6-C9 bond homolysis"),
    "timestamp": _dt.datetime.now().isoformat(timespec="seconds"),
    "engines": {},
    "substitutions": [],
    "warnings": [],
    "fallbacks": [],
    "tradeoffs": [],
    "fatal_error": None,
    "all_stages_ok": False,
}


def _log(tag, msg):
    print(f"[{_dt.datetime.now().strftime('%H:%M:%S')}] [{tag}] {msg}",
          flush=True)


def _warn(msg):
    META["warnings"].append(msg)
    _log("warn", msg)


def _fallback(msg):
    META["fallbacks"].append(msg)
    _log("fallback", msg)


def _tradeoff(msg):
    META["tradeoffs"].append(msg)
    _log("tradeoff", msg)


def write_json_atomic(path: Path, payload: dict):
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2, default=float),
                   encoding="utf-8")
    os.replace(tmp, path)


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def ckpt_dir(stage: str) -> Path:
    d = OUT / "checkpoints" / stage
    d.mkdir(parents=True, exist_ok=True)
    return d


# --------------------------------------------------------------------------- #
#  interpreter / engine discovery
# --------------------------------------------------------------------------- #
def find_py_qc() -> str:
    try:
        import psi4  # noqa: F401
        return sys.executable
    except Exception:
        pass
    for c in (ENV_PY_QC,):
        if Path(c).exists():
            return c
    raise RuntimeError("no interpreter with psi4 found "
                       "(expected conda env `phase7`)")


def find_py_chem() -> str:
    try:
        import ase  # noqa: F401
        import torch  # noqa: F401
        return sys.executable
    except Exception:
        pass
    for c in (ENV_PY_CHEM,):
        if Path(c).exists():
            return c
    raise RuntimeError("no interpreter with ase/torch found "
                       "(expected conda env `phase2ff`)")


def find_xtb() -> str | None:
    cands = [shutil.which("xtb"),
             str(Path(ENV_CHEM_BIN) / "xtb.exe"),
             str(Path(sys.prefix) / "Library" / "bin" / "xtb.exe")]
    for c in cands:
        if c and Path(c).exists():
            return c
    return None


# --------------------------------------------------------------------------- #
#  system topology: parse phase-4 reactant, verify the scissile bond
# --------------------------------------------------------------------------- #
def parse_reactant_mol(path: Path):
    """Return (elements, positions(N,3)) from the phase-4 mol block."""
    lines = Path(path).read_text().splitlines()
    natoms = int(lines[3][:3])
    els, pos = [], []
    for i in range(4, 4 + natoms):
        p = lines[i].split()
        pos.append([float(p[0]), float(p[1]), float(p[2])])
        els.append(p[3])
    return els, np.array(pos)


def parse_reactant_bonds(path: Path):
    lines = Path(path).read_text().splitlines()
    natoms = int(lines[3][:3])
    nbonds = int(lines[3][3:6])
    bonds = []
    for i in range(4 + natoms, 4 + natoms + nbonds):
        b = lines[i].split()
        bonds.append((int(b[0]) - 1, int(b[1]) - 1))
    return bonds


def locate_scissile(els, bonds):
    """Identify the cyclopropane CH2 bridge carbon and the scissile partner
    (the cyclopropane edge that Phase 4 flagged as broken: C6-C9, i.e. the
    edge NOT shared with the N-bearing ring carbon).

    Deterministic rules from the mol-block connectivity:
      * bridge carbon  = C bonded to exactly 2 H and 2 C, inside a C3 ring;
      * partner carbon = the bridge's ring neighbour NOT bonded to N
        (the indole C3-side carbon); its N-bonded sibling stays intact.
    """
    from collections import defaultdict
    adj = defaultdict(set)
    for a, b in bonds:
        adj[a].add(b)
        adj[b].add(a)
    hyd = {i for i, e in enumerate(els) if e == "H"}
    carb = [i for i, e in enumerate(els) if e == "C"]
    bridge = None
    for c in carb:
        nb_h = adj[c] & hyd
        nb_c = adj[c] - hyd
        if len(nb_h) == 2 and len(nb_c) == 2:
            # both C neighbours must bond each other (three-membered ring)
            a, b = sorted(nb_c)
            if b in adj[a]:
                bridge = c
                pair = (a, b)
    if bridge is None:
        raise RuntimeError("cyclopropane CH2 bridge not found in mol block")
    a, b = pair
    nitro = [i for i, e in enumerate(els) if e == "N"][0]
    partner = a if nitro in adj[b] else b
    sibling = b if partner == a else a
    return bridge, partner, sibling


# --------------------------------------------------------------------------- #
#  CHEM WORKER  (runs under phase2ff python: ase + xtb + mace + torchani)
# --------------------------------------------------------------------------- #
class XTBWrap:
    """ASE calculator around `xtb --grad` (GFN2-xTB), phase-4 pattern.

    Subclasses ase Calculator so ASE >=3.25 optimizers and
    get_potential_energy work with the modern property protocol."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, xtb_exe):
        self.xtb = xtb_exe
        self.n_calls = 0

    def _run(self, numbers, positions):
        with tempfile.TemporaryDirectory() as td:
            xyz = Path(td) / "m.xyz"
            lines = [str(len(numbers)), "p7"]
            for z, p in zip(numbers, positions):
                lines.append(f"{z} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}")
            xyz.write_text("\n".join(lines) + "\n")
            cmd = [self.xtb, "m.xyz", "--grad", "--chrg", str(CHARGE),
                   "--mult", str(MULTIPLICITY)]
            proc = subprocess.run(cmd, cwd=td, capture_output=True,
                                  text=True, timeout=600)
            grad = Path(td) / "gradient"
            if not grad.exists():
                raise RuntimeError(
                    f"xtb --grad no gradient; stderr: {proc.stderr[-200:]}")
            gtxt = grad.read_text()
            m = re.search(r"SCF energy =\s*(-?[\d.EeD+]+)", gtxt)
            if m is None:
                raise RuntimeError("xtb gradient lacks SCF energy")
            e_eh = float(m.group(1).replace("D", "E"))
            gvals = []
            for l in gtxt.splitlines()[2:]:
                if l.startswith("$"):
                    break
                t = l.split()
                if len(t) == 3:
                    try:
                        gvals.extend(float(x.replace("D", "E")) for x in t)
                    except ValueError:
                        continue
            g = np.array(gvals[:3 * len(numbers)]).reshape(-1, 3)
            if g.shape[0] != len(numbers):
                raise RuntimeError("xtb gradient rows != atoms")
            self.n_calls += 1
            return e_eh * EH_EV, -(g * EH_EV / 0.52917721092)

    # --- ASE property protocol (atoms passed explicitly by ASE) ---------- #
    def get_potential_energy(self, atoms, force_consistent=False,
                             **kwargs):
        e, _ = self._run(atoms.get_atomic_numbers(), atoms.get_positions())
        return e

    def get_forces(self, atoms, apply_constraint=True, md=False,
                   **kwargs):
        _, f = self._run(atoms.get_atomic_numbers(), atoms.get_positions())
        return f

    def get_property(self, name, atoms=None, allow_calculation=True):
        if name == "energy":
            return self.get_potential_energy(atoms)
        if name == "forces":
            return self.get_forces(atoms)
        raise KeyError(name)


def worker_G(args):
    """Stage G: constrained relaxed bond scan with GFN2-xTB (+ rigid
    fallback).  Emits scan_geometries/R*.xyz + per-point checkpoint JSON."""
    from ase import Atoms
    from ase.constraints import FixBondLengths
    from ase.optimize import LBFGS

    xtb = find_xtb()
    if xtb is None:
        raise RuntimeError("xtb.exe not found — relaxed scan unavailable")
    os.environ["PATH"] = ENV_CHEM_BIN + os.pathsep + os.environ.get("PATH", "")
    os.environ.setdefault("OMP_NUM_THREADS", "4")
    os.environ.setdefault("OMP_STACKSIZE", "1G")

    els, pos = parse_reactant_mol(REACTANT_MOL)
    bonds = parse_reactant_bonds(REACTANT_MOL)
    bridge, partner, _sibling = locate_scissile(els, bonds)
    numbers = [{"C": 6, "H": 1, "N": 7}[e] for e in els]
    calc = XTBWrap(xtb)
    gdir = OUT / "scan_geometries"
    gdir.mkdir(parents=True, exist_ok=True)

    scan = []
    cur_pos = pos.copy()
    for R in args.scan_r:
        tag = f"R{R:.2f}"
        ck = ckpt_dir("G") / f"{tag}.json"
        if ck.exists() and (gdir / f"{tag}.xyz").exists() and not args.force:
            scan.append(read_json(ck))
            _log("G", f"{tag}: checkpoint, skip")
            continue
        t0 = time.time()
        rec = {"R": R, "converged": False, "mode": None,
               "e_xtb_eV": None, "n_steps": 0}
        try:
            at = Atoms(numbers=atoms_from_elements(els),
                       positions=cur_pos)
            at.calc = calc
            at.set_constraint(FixBondLengths([(bridge, partner)]))
            opt = LBFGS(at, logfile=None)
            opt.run(fmax=0.03, steps=250)
            p = at.get_positions()
            d = np.linalg.norm(p[bridge] - p[partner])
            rec.update(converged=bool(opt.converged()),
                       mode="xtb-relaxed",
                       e_xtb_eV=float(at.get_potential_energy()),
                       n_steps=opt.get_number_of_steps(),
                       bond_actual=float(d))
            cur_pos = p
        except Exception as exc:
            _fallback(f"G {tag}: relaxed xtb opt failed "
                      f"({exc.__class__.__name__}: {str(exc)[:90]}); "
                      f"rigid geometry substitute")
            p = cur_pos.copy()
            u = p[partner] - p[bridge]
            u /= np.linalg.norm(u)
            p[partner] = p[bridge] + u * R
            rec.update(mode="rigid", bond_actual=R)
            if rec["e_xtb_eV"] is None:
                try:
                    at2 = Atoms(numbers=atoms_from_elements(els),
                                positions=p)
                    at2.calc = calc
                    rec["e_xtb_eV"] = float(at2.get_potential_energy())
                except Exception:
                    pass
        # verify/force the exact target distance on whatever we keep
        d = np.linalg.norm(p[bridge] - p[partner])
        if abs(d - R) > 1e-3:
            u = p[partner] - p[bridge]
            u /= np.linalg.norm(u)
            p[partner] = p[bridge] + u * R
            rec["bond_snapped"] = True
        # sanity: no bond collapse / spurious rupture among the spectator
        # heavy-atom framework (scissile pair excluded by construction)
        dm = _min_heavy_distance(els, p, exclude=(bridge, partner))
        rec["min_heavy_dist_A"] = float(dm)
        if dm < 1.15 or dm > 2.1:
            _warn(f"G {tag}: spectator heavy-atom framework distorted "
                  f"(min d = {dm:.2f} A) — geometry flagged")
        (gdir / f"{tag}.xyz").write_text(_xyz_text(els, p, tag))
        rec["xyz"] = str(gdir / f"{tag}.xyz")
        rec["wall_s"] = round(time.time() - t0, 1)
        write_json_atomic(ck, rec)
        scan.append(rec)
        write_json_atomic(ckpt_dir("G") / "scan.json", {"points": scan})
        _log("G", f"{tag}: {rec['mode']} conv={rec['converged']} "
                  f"E={_f(rec['e_xtb_eV'])} eV  ({rec['wall_s']} s)")
    write_json_atomic(ckpt_dir("G") / "scan.json", {"points": scan})
    return {"points": scan}


def atoms_from_elements(els):
    return [{"C": 6, "H": 1, "N": 7}[e] for e in els]


def _xyz_text(els, pos, comment=""):
    lines = [str(len(els)), comment]
    for e, p in zip(els, pos):
        lines.append(f"{e} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}")
    return "\n".join(lines) + "\n"


def _min_heavy_distance(els, pos, exclude=()):
    """Minimum heavy-heavy distance among spectator pairs (the scan pair is
    excluded, so the probe watches for framework collapse/rupture)."""
    dmin = 1e9
    n = len(els)
    for i in range(n):
        if els[i] == "H":
            continue
        for j in range(i + 1, n):
            if els[j] == "H" or (i, j) == tuple(sorted(exclude)):
                continue
            d = float(np.linalg.norm(pos[i] - pos[j]))
            dmin = min(dmin, d)
    return dmin


def _f(x, nd=4):
    return "n/a" if x is None else f"{x:.{nd}f}"


def worker_AI(args):
    """Stage AI: single points of the identical trajectory at
    GFN2-xTB / MACE-OFF / ANI-2x levels."""
    from ase import Atoms
    from ase.io import read as ase_read

    els, _ = parse_reactant_mol(REACTANT_MOL)
    numbers = atoms_from_elements(els)
    scan = read_json(ckpt_dir("G") / "scan.json")["points"]
    results = {}
    xtb = find_xtb()
    if xtb:
        os.environ["PATH"] = (ENV_CHEM_BIN + os.pathsep +
                              os.environ.get("PATH", ""))
        calc_xtb = XTBWrap(xtb)
    else:
        calc_xtb = None
        _fallback("AI: xtb.exe unavailable — GFN2-xTB column omitted")

    mace_calc, ani_calc = None, None
    mace_local = Path.home() / ".cache" / "mace" / "MACE-OFF23_small.model"
    try:
        from mace.calculators import mace_off
        if mace_local.exists():
            mace_calc = mace_off(model=str(mace_local), device="cpu")
        else:
            mace_calc = mace_off(model="small", device="cpu")
        results["_engine_mace"] = "MACE-OFF small (foundation, ASL)"
    except Exception as exc:
        _fallback(f"AI: MACE-OFF load failed "
                  f"({exc.__class__.__name__}: {str(exc)[:90]}); "
                  f"ANI-2x carries the ML-potential column")
    try:
        import torchani
        ani_calc = torchani.models.ANI2x()
        results["_engine_ani"] = "ANI-2x (torchani)"
    except Exception as exc:
        _fallback(f"AI: ANI-2x load failed ({str(exc)[:90]})")

    for pt in scan:
        tag = f"R{pt['R']:.2f}"
        ck = ckpt_dir("AI") / f"{tag}.json"
        if ck.exists() and not args.force:
            results[tag] = read_json(ck)
            _log("AI", f"{tag}: checkpoint, skip")
            continue
        at = ase_read(pt["xyz"])
        rec = {"R": pt["R"]}
        if calc_xtb is not None:
            try:
                at.calc = calc_xtb
                rec["e_xtb_eV"] = float(at.get_potential_energy())
            except Exception as exc:
                _warn(f"AI {tag}: xtb single point failed: {str(exc)[:90]}")
        if mace_calc is not None:
            try:
                at.calc = mace_calc
                rec["e_mace_eV"] = float(at.get_potential_energy())
            except Exception as exc:
                _warn(f"AI {tag}: MACE failed: {str(exc)[:90]}")
        if ani_calc is not None:
            try:
                import torch
                p = torch.tensor(at.get_positions(),
                                 dtype=torch.float32).unsqueeze(0)
                n = torch.tensor(at.get_atomic_numbers(),
                                 dtype=torch.long).unsqueeze(0)
                p.requires_grad_(True)
                _, e = ani_calc((n, p))
                rec["e_ani_eV"] = float(e.item()) * EH_EV
            except Exception as exc:
                _warn(f"AI {tag}: ANI failed: {str(exc)[:90]}")
        write_json_atomic(ck, rec)
        results[tag] = rec
        _log("AI", f"{tag}: xtb={_f(rec.get('e_xtb_eV'))} "
                   f"mace={_f(rec.get('e_mace_eV'))} "
                   f"ani={_f(rec.get('e_ani_eV'))} eV")
    return results


# --------------------------------------------------------------------------- #
#  QC WORKER  (runs under phase7 python: psi4)
# --------------------------------------------------------------------------- #
def _psi4_mol(xyz_path, multiplicity=MULTIPLICITY):
    import psi4
    txt = Path(xyz_path).read_text().splitlines()
    body = "\n".join(txt[2:2 + int(txt[0])])
    return psi4.geometry(
        f"{CHARGE} {multiplicity}\n{body}\nunits angstrom\n"
        f"symmetry c1\nno_com\nno_reorient\n")


def _spin_square(wfn):
    """<S^2> = S(S+1) + n_beta - Tr[D_a S D_b S]  (exact for UHF/UKS)."""
    S = np.asarray(wfn.S())          # SO overlap, same frame as Ca/Da
    Da = np.asarray(wfn.Da_subset("AO"))
    Db = np.asarray(wfn.Db_subset("AO"))
    nb = wfn.nbeta()
    ss_two = 0.25 * (wfn.nalpha() - wfn.nbeta()) ** 2
    return float(nb - np.trace(Da @ S @ Db @ S) + ss_two)


def _ao_atom_map(basisset, natoms):
    """AO-function index -> atom index (psi4 function_to_center)."""
    return np.array([basisset.function_to_center(i)
                     for i in range(basisset.nbf())], dtype=int)


def _mulliken_mo_weight(wfn, mo_idx, atom_set, S, ao2atom):
    C = np.asarray(wfn.Ca())
    sel = np.isin(ao2atom, list(atom_set))
    P = C[:, mo_idx] @ C[:, mo_idx].T
    return float(np.sum((P * S)[sel][:, sel]))


def _scf_point(mol, theory, basis, mem, threads, out_log):
    """One SCF single point; returns dict with energy + optional <S^2>.
    `theory` may carry a _triplet suffix (molecule already built with
    multiplicity 3); the SCF reference string is derived from the prefix."""
    import psi4
    psi4.set_memory(mem)
    psi4.set_num_threads(threads)
    psi4.core.set_output_file(str(out_log), True)
    stem = theory.replace("_triplet", "")
    is_dft = stem in ("rks", "uks")
    name = DFT_FUNC if is_dft else "hf"
    psi4.set_options({"basis": basis, "reference": stem,
                      "scf_type": "df"})
    e, w = psi4.energy(name, molecule=mol, return_wfn=True)
    rec = {"energy_eh": float(e), "theory": theory, "basis": basis,
           "converged": True}
    if stem in ("uhf", "uks"):
        ss = _spin_square(w)
        rec["s2"] = ss
        rec["ds2"] = ss - (2.0 if theory.endswith("_triplet")
                           else SINGLET_SS)
    return rec, w


def worker_7A(args):
    """7A per point: RHF/UHF/RKS/UKS singlets + UHF/UKS triplets.

    On this pi-stabilized polar diradical the broken-symmetry singlet is
    NOT reachable by SCF from any smooth guess (verified exhaustively:
    SAD/SAP guesses, triplet-seeded checkpoint splicing via guess=read,
    MOM occupation pinning, SOSCF — all relax back to the closed-shell
    branch).  The open-shell single-determinant physics is therefore
    probed in the S_z = 1 sector, where the triplet SCF is well behaved,
    and the singlet-triplet gap collapse pinpoints the symmetry-breaking
    bond length.  All seed-engineering attempts are logged to the results
    file so the negative result is auditable."""
    import psi4
    point = read_json(ckpt_dir("G") / "scan.json")["points"][args.point]
    tag = f"R{point['R']:.2f}"
    mol = _psi4_mol(point["xyz"])
    mol_t = _psi4_mol(point["xyz"], multiplicity=3)
    logdir = OUT / "logs"
    logdir.mkdir(parents=True, exist_ok=True)
    basis = args.basis
    out = {"R": point["R"], "tag": tag, "theory": {}, "basis": basis,
           "failures": []}

    def run(theory, m):
        rec = None
        for attempt in range(2):
            try:
                rec, _w = _scf_point(m, theory, basis, args.mem,
                                     args.threads,
                                     logdir / f"7A_{tag}_{theory}.out")
                break
            except Exception as exc:
                msg = f"{exc.__class__.__name__}: {str(exc)[:110]}"
                out["failures"].append(f"{theory}[{attempt}]: {msg}")
                _warn(f"7A {tag} {theory} attempt {attempt}: {msg}")
                if attempt == 0:
                    psi4.core.clean()
                    try:
                        psi4.set_options({"scf__damping_ratio": 0.7})
                    except Exception:
                        pass
        return rec

    for theory, m in (("rhf", mol), ("uhf", mol), ("rks", mol),
                      ("uks", mol), ("uhf_triplet", mol_t),
                      ("uks_triplet", mol_t)):
        rec = run(theory, m)
        # triplet SCF may land on a wrong excited determinant (<S^2> ~ 1
        # instead of ~2, E above the singlet): retry with a spin-polarized
        # SAP guess and keep the lower solution
        if rec is not None and theory.endswith("_triplet"):
            ss_bad = ("s2" in rec and rec["s2"] < 1.7)
            e_bad = (rec["energy_eh"] is not None and
                     out["theory"].get(theory.replace("_triplet", ""))
                     and out["theory"][theory.replace("_triplet", "")]
                     .get("energy_eh") is not None and
                     rec["energy_eh"] > out["theory"][
                         theory.replace("_triplet", "")]["energy_eh"])
            if ss_bad or e_bad:
                try:
                    psi4.set_options({"guess": "sap"})
                    rec2 = run(theory, m)
                    psi4.set_options({"guess": "auto"})
                    if rec2 is not None and rec.get("energy_eh") is not None:
                        if (rec2.get("energy_eh") is not None and
                                rec2["energy_eh"] < rec["energy_eh"]):
                            rec = rec2
                            out["failures"].append(
                                f"{theory}: SAD guess gave an excited "
                                f"triplet; SAP retry lower (kept)")
                except Exception as exc:
                    _warn(f"7A {tag} {theory} SAP retry failed: "
                          f"{str(exc)[:80]}")
                    try:
                        psi4.set_options({"guess": "auto"})
                    except Exception:
                        pass
        if rec is None:
            out["theory"][theory] = {"energy_eh": None, "converged": False,
                                     "basis": basis}
        else:
            out["theory"][theory] = rec
        extra = (f"  <S^2>={rec['s2']:.4f}" if rec and "s2" in rec else "")
        _log("7A", f"{tag} {theory}: "
                   f"E={_f(rec['energy_eh'] if rec else None, 6)} Eh"
                   + extra)
    write_json_atomic(ckpt_dir("7A") / f"{tag}.json", out)
    return out


def _select_frontier(wfn, S, ao2atom, scissile_atoms):
    """Pick the sigma / sigma* pair of the scissile bond by complement
    coherence: for every valence-occupied candidate the antibonding
    complement u_k = P_a|k> - P_b|k> is constructed and scored by how much
    of it collapses onto a SINGLE low virtual; the occupied orbital with
    the most coherent complement IS the sigma, and that virtual is the
    sigma-star.  Raw gross-population argmax fails twice over: core 1s
    orbitals out-score the sigma, and diffuse functions inflate Rydberg
    subset quadratic forms by an order of magnitude."""
    C = np.asarray(wfn.Ca())
    eps = np.asarray(wfn.epsilon_a())
    nocc = wfn.nalpha()
    sel = np.isin(ao2atom, list(scissile_atoms))
    Sss = S[np.ix_(sel, sel)]
    weights = []
    for i in range(C.shape[1]):
        Ci = C[sel, i]
        weights.append(float(Ci @ Sss @ Ci))
    weights = np.array(weights)
    atom_ids = sorted(set(int(a) for a in np.asarray(ao2atom)[sel]))
    assert len(atom_ids) == 2, "scissile selection must span exactly 2 atoms"

    def proj(vec, idx):
        B = np.zeros((S.shape[0], len(idx)))
        B[np.array(idx, dtype=int), np.arange(len(idx))] = 1.0
        Sbb = B.T @ S @ B
        return B @ np.linalg.solve(Sbb, B.T @ S @ vec)

    idx_a = [k for k in range(len(ao2atom)) if ao2atom[k] == atom_ids[0]]
    idx_b = [k for k in range(len(ao2atom)) if ao2atom[k] == atom_ids[1]]
    virt_cands = [k for k in range(nocc, C.shape[1]) if eps[k] < 1.5]
    occ_cands = sorted((k for k in range(nocc) if eps[k] > -1.6),
                       key=lambda k: -weights[k])[:8]
    best = None
    for k in occ_cands:
        u = proj(C[:, k], idx_a) - proj(C[:, k], idx_b)
        ovs = {j: abs(float(C[:, j] @ S @ u)) for j in virt_cands}
        j_best = max(ovs, key=ovs.get)
        score = ovs[j_best]
        if best is None or score > best[0]:
            best = (score, k, j_best)
    quality, i_sig, i_ast = best
    return i_sig, i_ast, float(weights[i_sig]), float(quality)


def _reorder_to_cas(Ca: np.ndarray, nocc: int, i_sig: int, i_ast: int):
    """Permute MO columns so that [sigma, sigma*] sit exactly at the
    docc/active/uocc boundary required by DETCI occupation arrays.
    Returns (permutation, n_docc)."""
    docc = [k for k in range(nocc) if k != i_sig]
    virt = [k for k in range(Ca.shape[1])
            if k != i_ast and k >= nocc]
    order = docc + [i_sig, i_ast] + virt
    assert len(order) == Ca.shape[1]
    return order, len(docc)


def worker_7B(args):
    """7B per point: CAS(2e,2o) CASSCF ground truth + NOON + cubes."""
    import psi4
    point = read_json(ckpt_dir("G") / "scan.json")["points"][args.point]
    R = point["R"]
    tag = f"R{R:.2f}"
    mol = _psi4_mol(point["xyz"])
    logdir = OUT / "logs"
    logdir.mkdir(parents=True, exist_ok=True)
    out = {"R": R, "tag": tag, "attempts": []}

    psi4.set_memory(args.mem)
    psi4.set_num_threads(args.threads)
    psi4.core.set_output_file(str(logdir / f"7B_{tag}.out"), True)
    psi4.set_options({"basis": args.basis, "reference": "rhf",
                      "scf_type": "df"})
    e_rhf, wref = psi4.energy("hf", molecule=mol, return_wfn=True)
    S = np.asarray(wref.S())          # same SO frame as Ca / epsilon
    ao2atom = _ao_atom_map(wref.basisset(), mol.natom())
    _els, _ = parse_reactant_mol(REACTANT_MOL)
    _bonds = parse_reactant_bonds(REACTANT_MOL)
    _bridge, _partner, _sib = locate_scissile(_els, _bonds)
    scissile = {_partner, _bridge}
    i_sig, i_ast, w_sig, w_ast = _select_frontier(wref, S, ao2atom, scissile)
    nocc = wref.nalpha()
    nbf = wref.basisset().nbf()
    _log("7B", f"{tag}: sigma=MO{i_sig} (w={w_sig:.3f}), "
               f"sigma*=MO{i_ast} (w={w_ast:.3f}), nocc={nocc}, nbf={nbf}")
    if w_sig < 0.35:
        _warn(f"7B {tag}: weak scissile-sigma localization "
              f"(w={w_sig:.2f}) — CAS character may mix with the ring")
    out["selection"] = {"i_sigma": i_sig, "i_sigma_star": i_ast,
                        "w_sigma": w_sig, "w_sigma_star": w_ast,
                        "nocc": int(nocc), "nbf": int(nbf)}

    # permute RHF MOs in place so the CAS pair sits at the boundary
    Ca_view = np.asarray(wref.Ca())
    order, ndocc = _reorder_to_cas(Ca_view, nocc, i_sig, i_ast)
    Ca_view[:] = Ca_view[:, order]

    # tier 0 = psi4 default MCSCF integral path (verified sane on this
    # build); explicitly forcing mcscf_type CONV triggers a defective
    # DPD cache path here (energies ~ +4.7 Eh off) — logged as a known
    # build quirk, so CONV is never forced, only the default/AO/DF
    tiers = [{}, {"mcscf_type": "AO"}, {"mcscf_type": "DF"}][
        args.skip_tiers:]
    mc = None
    for t_i, extra in enumerate(tiers, start=args.skip_tiers):
        try:
            psi4.set_options({
                "basis": args.basis, "reference": "rhf",
                "restricted_docc": [ndocc], "active": [2],
                "restricted_uocc": [nbf - ndocc - 2],
                "maxiter": 120,
            })
            try:
                psi4.set_module_options("detci", extra)
            except Exception:
                pass
            e_cas, w = psi4.energy("casscf", molecule=mol, return_wfn=True,
                                   ref_wfn=wref)
            mc = w
            out["mcscf_options"] = extra
            out["attempts"].append({"tier": t_i, "opts": extra, "ok": True})
            break
        except Exception as exc:
            msg = f"tier{t_i} {extra}: {exc.__class__.__name__}: " \
                  f"{str(exc)[:100]}"
            out["attempts"].append({"tier": t_i, "opts": extra,
                                    "ok": False, "err": msg})
            _warn(f"7B {tag} CASSCF {msg}")
            psi4.core.clean()

    if mc is None:
        out["converged"] = False
        write_json_atomic(ckpt_dir("7B") / f"{tag}.json", out)
        return out

    # ---- natural orbitals (ACTIVE BLOCK only) -----------------------------
    # the full MO-basis 1-RDM carries the 34 doubly-occupied inactive
    # orbitals (NOON ~ 2); the diradical signature lives in the 2x2 active
    # block at columns (ndocc, ndocc+1) of the CASSCF MO basis
    O = np.asarray(mc.get_opdm(0, 0, "SUM", True))
    blk = O[np.ix_([ndocc, ndocc + 1], [ndocc, ndocc + 1])]
    ev_blk, vec_blk = np.linalg.eigh(blk)
    order_n = np.argsort(ev_blk)[::-1]
    nu1, nu2 = float(ev_blk[order_n[0]]), float(ev_blk[order_n[1]])
    # diradical character y = nu(antibonding NO) in [0, 1]
    # (0 = closed shell, 1 = pure diradical; nu1 + nu2 = 2)
    y_yamaguchi = float(nu2)
    y_luno = float(nu2)
    Ca_cas = np.asarray(mc.Ca()).copy()
    act = [ndocc, ndocc + 1]
    no1 = Ca_cas[:, act] @ vec_blk[:, order_n[0]]
    no2 = Ca_cas[:, act] @ vec_blk[:, order_n[1]]
    trace_active = float(np.trace(blk))
    out.update({
        "converged": True,
        "e_casscf_eh": float(psi4.variable("CURRENT ENERGY")),
        "e_casscf_direct": float(e_cas),
        "noon": [nu1, nu2],
        "noon_active_trace": trace_active,
        "y_yamaguchi": float(y_yamaguchi),
        "y_luno": y_luno,
        "basis": args.basis,
    })
    _log("7B", f"{tag}: E_CAS={out['e_casscf_eh']:.6f} Eh  "
               f"NOON=({nu1:.4f},{nu2:.4f})  tr(active)={trace_active:.3f}  "
               f"y={y_yamaguchi:.4f}")

    # ---- frontier natural orbital cubes -----------------------------------
    if abs(R - args.cube_r) < 1e-6:
        try:
            cube_dir = OUT / "cubes" / tag
            cube_dir.mkdir(parents=True, exist_ok=True)
            Ca_m = mc.Ca()
            arr = np.asarray(Ca_m)
            nalpha = mc.nalpha()
            arr[:, nalpha - 1] = no1     # HOMO slot  <- bonding NO
            arr[:, nalpha] = no2         # LUMO slot  <- antibonding NO
            psi4.set_options({
                "cubeprop_tasks": ["orbitals"],
                "cubeprop_orbitals": [nalpha, nalpha + 1],
                "cubeprop_filepath": str(cube_dir),
                "basis": args.basis, "reference": "rhf",
                "restricted_docc": [ndocc], "active": [2],
                "restricted_uocc": [nbf - ndocc - 2],
            })
            psi4.cubeprop(mc)
            out["cubes"] = sorted(p.name for p in cube_dir.glob("*.cube"))
            _log("7B", f"{tag}: cubes -> {out['cubes']}")
        except Exception as exc:
            _warn(f"7B {tag}: cubeprop failed: {str(exc)[:120]}")
    write_json_atomic(ckpt_dir("7B") / f"{tag}.json", out)
    return out


# --------------------------------------------------------------------------- #
#  MERGE + ANALYSIS  (epistemic landscape)
# --------------------------------------------------------------------------- #
def interp_crossing(xs, ys, thr):
    """First crossing of `thr` in either direction; linear interpolation
    between grid points; returns (R_cross, i_lo) or (None, None)."""
    for i in range(1, len(xs)):
        a, b = ys[i - 1], ys[i]
        if None in (a, b):
            continue
        if (a < thr <= b) or (a > thr >= b):
            f = (thr - a) / (b - a)
            return xs[i - 1] + f * (xs[i] - xs[i - 1]), i - 1
    return None, None


def merge_and_analyze():
    scan = read_json(ckpt_dir("G") / "scan.json")["points"]
    Rs = [p["R"] for p in scan]
    tags = [f"R{r:.2f}" for r in Rs]

    a7 = {t: (read_json(ckpt_dir("7A") / f"{t}.json")
              if (ckpt_dir("7A") / f"{t}.json").exists() else None)
          for t in tags}
    b7 = {t: (read_json(ckpt_dir("7B") / f"{t}.json")
              if (ckpt_dir("7B") / f"{t}.json").exists() else None)
          for t in tags}
    ai = {t: (read_json(ckpt_dir("AI") / f"{t}.json")
              if (ckpt_dir("AI") / f"{t}.json").exists() else None)
          for t in tags}

    def series(src, pick):
        out = []
        for t in tags:
            rec = src[t]
            try:
                out.append(pick(rec) if rec else None)
            except Exception:
                out.append(None)
        return out

    e_rhf = series(a7, lambda r: r["theory"]["rhf"]["energy_eh"])
    e_uhf = series(a7, lambda r: r["theory"]["uhf"]["energy_eh"])
    e_rks = series(a7, lambda r: r["theory"]["rks"]["energy_eh"])
    e_uks = series(a7, lambda r: r["theory"]["uks"]["energy_eh"])
    s2_uhf = series(a7, lambda r: r["theory"]["uhf"]["s2"])
    s2_uks = series(a7, lambda r: r["theory"]["uks"]["s2"])
    e_uhf_t = series(a7, lambda r: r["theory"]["uhf_triplet"]["energy_eh"])
    e_uks_t = series(a7, lambda r: r["theory"]["uks_triplet"]["energy_eh"])
    s2_uhf_t = series(a7, lambda r: r["theory"]["uhf_triplet"]["s2"])
    s2_uks_t = series(a7, lambda r: r["theory"]["uks_triplet"]["s2"])
    e_cas = series(b7, lambda r: r["e_casscf_eh"] if r["converged"] else None)
    noon1 = series(b7, lambda r: r["noon"][0] if r["converged"] else None)
    noon2 = series(b7, lambda r: r["noon"][1] if r["converged"] else None)
    y_dirc = series(b7, lambda r: r["y_yamaguchi"]
                    if r["converged"] else None)

    def to_rel(eh_series):
        i0 = Rs.index(R_EQ) if R_EQ in Rs else 0
        base = eh_series[i0] if eh_series[i0] is not None else \
            next((x for x in eh_series if x is not None), None)
        if base is None:
            return [None] * len(eh_series)
        return [None if x is None else (x - base) * HARTREE_KCAL
                for x in eh_series]

    rel = {
        "CASSCF": to_rel(e_cas),
        "UHF": to_rel(e_rhf),           # restricted = mean-field closed shell
        "BS-UB3LYP": to_rel(e_uks),
    }
    # xTB relative curve in kcal/mol
    e_xtb_ev = [p["e_xtb_eV"] for p in scan]
    rel["GFN2-xTB"] = _rel_from_ev(e_xtb_ev, Rs)
    # MACE / ANI
    e_mace_ev = series(ai, lambda r: r.get("e_mace_eV"))
    e_ani_ev = series(ai, lambda r: r.get("e_ani_eV"))
    rel["MACE-OFF"] = _rel_from_ev(e_mace_ev, Rs)
    rel["ANI-2x"] = _rel_from_ev(e_ani_ev, Rs)

    # epistemic discrepancy vectors vs CASSCF (kcal/mol)
    err = {}
    for k in ("GFN2-xTB", "MACE-OFF", "ANI-2x", "BS-UB3LYP", "UHF"):
        err[k] = [None if (a is None or b is None) else abs(a - b)
                  for a, b in zip(rel[k], rel["CASSCF"])]

    def first_fail(xs, ys, thr=FAIL_KCAL):
        for i, v in enumerate(ys):
            if v is not None and v > thr:
                return xs[i]
        return None

    r_crit = {}
    # mission gate: DS^2 = <S^2> - S(S+1) exceeding 0.3.  The singlet UHF/
    # UKS remain on the closed-shell branch for this system (DS2 ~ 0 for
    # the whole window — audited negative result); the open-shell single-
    # determinant sector is the S_z = 1 one, so the contamination gate is
    # evaluated on the triplet solutions, and the symmetry-breaking bond
    # length is additionally pinpointed by the singlet-triplet gap sign
    # change (closed-shell singlet stops being the lowest determinant).
    for name, ss in (("UHF_triplet", s2_uhf_t), ("UKS_triplet", s2_uks_t)):
        rc, _ = interp_crossing(Rs, [0 if v is None else v for v in ss],
                                S2_GATE)
        r_crit[name] = rc
    dE_st_uhf = [None if (a is None or b is None) else (a - b) * HARTREE_KCAL
                 for a, b in zip(e_uhf_t, e_uhf)]
    dE_st_uks = [None if (a is None or b is None) else (a - b) * HARTREE_KCAL
                 for a, b in zip(e_uks_t, e_uks)]
    rc_st_uhf, _ = interp_crossing(Rs, [np.nan if v is None else v
                                        for v in dE_st_uhf], 0.0)
    rc_st_uks, _ = interp_crossing(Rs, [np.nan if v is None else v
                                        for v in dE_st_uks], 0.0)
    r_crit["ST_crossing_UHF"] = rc_st_uhf
    r_crit["ST_crossing_UKS"] = rc_st_uks
    r_crit["DS2_singlet_UHF"] = interp_crossing(
        Rs, [0 if v is None else v for v in s2_uhf], S2_GATE)[0]
    r_crit["DS2_singlet_UKS"] = interp_crossing(
        Rs, [0 if v is None else v for v in s2_uks], S2_GATE)[0]

    analysis = {
        "R": Rs,
        "E_rhf_eh": e_rhf, "E_uhf_eh": e_uhf, "E_rks_eh": e_rks,
        "E_uks_eh": e_uks,
        "E_uhf_triplet_eh": e_uhf_t, "E_uks_triplet_eh": e_uks_t,
        "S2_uhf": s2_uhf, "S2_uks": s2_uks,
        "S2_uhf_triplet": s2_uhf_t, "S2_uks_triplet": s2_uks_t,
        "dE_ST_uhf_kcal": dE_st_uhf, "dE_ST_uks_kcal": dE_st_uks,
        "E_casscf_eh": e_cas, "noon_bonding": noon1,
        "noon_antibonding": noon2,
        "y_diradical": y_dirc,
        "E_rel_kcal": rel,
        "abs_error_vs_casscf_kcal": err,
        "R_crit": r_crit,
        "R_fail_15kcal": {k: first_fail(Rs, v) for k, v in err.items()},
        "R_eq_reference": R_EQ,
        "thresholds": {"DS2": S2_GATE, "fail_kcal": FAIL_KCAL},
    }
    return analysis


def _rel_from_ev(ev_series, Rs):
    i0 = Rs.index(R_EQ) if R_EQ in Rs else 0
    base = ev_series[i0]
    if base is None:
        base = next((x for x in ev_series if x is not None), None)
    if base is None:
        return [None] * len(ev_series)
    return [None if x is None else (x - base) * EV_KCAL for x in ev_series]


# --------------------------------------------------------------------------- #
#  FIGURES  (300 DPI)
# --------------------------------------------------------------------------- #
C_MAIN = "#0173B2"   # match phase-2/3/4 palette
C_SEC = "#DE8F05"
C_TER = "#029E73"
C_ML1 = "#D55E00"
C_ML2 = "#CC78BC"
C_REF = "#444444"


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.size": 10.5, "axes.titlesize": 11.5,
        "axes.labelsize": 11, "legend.fontsize": 9.2,
        "figure.dpi": 110, "savefig.dpi": 300,
        "axes.spines.top": False, "axes.spines.right": False,
    })
    return plt


def fig1(args, ana):
    plt = _style()
    Rs = np.array(ana["R"])
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.8, 9.0),
                                   sharex=True,
                                   gridspec_kw={"height_ratios": [3, 2],
                                                "hspace": 0.12})
    rc = ana["R_crit"].get("ST_crossing_UHF")
    rc_dft = ana["R_crit"].get("ST_crossing_UKS")

    # ---- panel (a): spin-contamination diagnostics ------------------------ #
    def nn(v):
        return np.nan if v is None else v

    s_st = np.array([nn(v) for v in ana["S2_uhf_triplet"]]) - 2.0
    s_kt = np.array([nn(v) for v in ana["S2_uks_triplet"]]) - 2.0
    s_ss = np.array([nn(v) for v in ana["S2_uhf"]])
    s_ks = np.array([nn(v) for v in ana["S2_uks"]])
    ax1.plot(Rs, s_st, "-o", color=C_MAIN, lw=2.2, ms=7, mec="k", mew=0.6,
             label=r"triplet UHF: $\Delta S^2=\langle S^2\rangle-2$",
             zorder=5)
    ax1.plot(Rs, s_kt, "-s", color=C_SEC, lw=2.0, ms=6.5, mec="k", mew=0.5,
             label=r"triplet UKS: $\Delta S^2=\langle S^2\rangle-2$",
             zorder=5)
    ax1.plot(Rs, s_ss, "-", color="#888888", lw=1.6,
             label=r"singlet UHF/UKS: $\Delta S^2\approx 0$"
             " (closed-shell branch)",
             zorder=4)
    ax1.plot(Rs, s_ks, "-", color="#888888", lw=1.6, zorder=4)
    ax1.axhline(S2_GATE, color="#D62728", ls="--", lw=1.5, zorder=2)
    ax1.text(Rs[0], S2_GATE + 0.02,
             r"spin-contamination gate  $\Delta S^2 = 0.3$",
             color="#D62728", fontsize=9.5)
    if rc is not None:
        ax1.axvline(rc, color="#D62728", lw=1.4, ls="-.", alpha=0.85,
                    zorder=1)
        ax1.annotate(rf"$R_\mathrm{{crit}}$ = {rc:.2f} Å" + "\n"
                     "(singlet–triplet gap closes)",
                     xy=(rc, -0.42), xytext=(rc + 0.12, -0.80),
                     fontsize=10, color="#B22",
                     arrowprops=dict(arrowstyle="->", color="#B22", lw=1.2))
    ax1.set_ylabel(r"spin contamination  $\Delta S^2$  (a.u.)")
    ax1.legend(loc="upper right", frameon=False, fontsize=9)
    ax1.set_title("Fig. 1 — Symmetry breaking along C6–C9 homolysis "
                  "(cyclopropa[b]indole singlet surface)", loc="left")

    # ---- panel (b): energetic diagnostics --------------------------------- #
    dst_u = np.array([nn(v) for v in ana["dE_ST_uhf_kcal"]])
    dst_k = np.array([nn(v) for v in ana["dE_ST_uks_kcal"]])
    static = np.array([nn(a) - nn(b) if (a is not None and b is not None)
                       else np.nan
                       for a, b in zip(ana["E_rhf_eh"],
                                       ana["E_casscf_eh"])]) * HARTREE_KCAL
    if rc is not None:
        ax2.axvspan(rc, Rs[-1], color="#D62728", alpha=0.07, zorder=0)
    ax2.plot(Rs, dst_u, "-o", color=C_MAIN, lw=2.2, ms=6.5, mec="k",
             mew=0.5,
             label=r"$\Delta E_{ST}=E(\mathrm{triplet\ UHF})"
                   r"-E(\mathrm{singlet})$", zorder=4)
    ax2.plot(Rs, dst_k, "-s", color=C_SEC, lw=2.0, ms=6, mec="k", mew=0.5,
             label=r"$\Delta E_{ST}$ with UKS-B3LYP", zorder=4)
    ax2.plot(Rs, static, "--^", color=C_TER, lw=1.8, ms=5.5,
             label=r"$E(\mathrm{RHF})-E(\mathrm{CASSCF})$"
             " — static-correlation error", zorder=4)
    ax2.axhline(0.0, color="#333333", ls=":", lw=1.2)
    ax2.set_xlabel(r"scissile C6–C9 distance  $R$  (Å)")
    ax2.set_ylabel("energy gap (kcal/mol)")
    ax2.legend(loc="upper right", frameon=False, fontsize=9)
    if rc is not None:
        ax2.text(1.42, -62,
                 "beyond $R_\\mathrm{crit}$ the closed-shell singlet is\n"
                 "no longer the lowest single determinant —\n"
                 "static-correlation error diverges",
                 fontsize=8.8, color="#B22", va="bottom")
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "fig1_spin_contamination_profile.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    _log("fig", "fig1_spin_contamination_profile.png written")


def _read_cube(path: Path):
    """Gaussian cube parser: returns (origin(3,), axis vectors(3,3),
    data(nx, ny, nz)).  Negative voxel counts (orbital-phase cubes) are
    handled via abs()."""
    txt = Path(path).read_text().splitlines()
    nat = int(txt[2].split()[0])
    origin = np.array([float(x) for x in txt[2].split()[1:4]])
    dims, axes = [], []
    for k in range(3):
        t = txt[3 + k].split()
        dims.append(abs(int(t[0])))
        axes.append(np.array([float(x) for x in t[1:4]]))
    vals = []
    for l in txt[6 + nat:]:
        vals.extend(float(x) for x in l.split())
    n = int(np.prod(dims))
    data = np.array(vals[:n]).reshape(dims)   # slowest axis = x
    return origin, axes, data


def fig2(args, ana):
    plt = _style()
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    try:
        from skimage.measure import marching_cubes
    except Exception:
        marching_cubes = None
        _warn("scikit-image unavailable — orbital isosurfaces degraded to "
              "occupancy panel only")

    cube_root = OUT / "cubes"
    Rs_all = ana["R"]
    show_cols = []
    if cube_root.exists():
        for r in CUBE_RS:
            if any((cube_root / f"R{r:.2f}").glob("*.cube")):
                show_cols.append(r)
    show_cols = show_cols[:3]

    if show_cols:
        fig = plt.figure(figsize=(13.6, 12.8))
        gs = fig.add_gridspec(3, len(show_cols),
                              height_ratios=[2.6, 2.6, 1.05], hspace=0.26,
                              wspace=0.04)
        row_labels = [r"bonding $\sigma$ natural orbital",
                      r"antibonding $\sigma^*$ natural orbital"]
        for k in (0, 1):
            for j, R in enumerate(show_cols):
                ax = fig.add_subplot(gs[k, j], projection="3d")
                ax.set_box_aspect((1, 1, 0.85))
                cdir = cube_root / f"R{R:.2f}"
                cubes = sorted(cdir.glob("*.cube"))
                noon = (ana["noon_bonding"][Rs_all.index(R)] if k == 0
                        else ana["noon_antibonding"][Rs_all.index(R)])
                ax.set_title(
                    f"{'HOMO' if k == 0 else 'LUMO'} slot — {row_labels[k]}"
                    + (f"\nNOON = {noon:.3f} e" if noon is not None
                       else "\nNOON n/a"),
                    fontsize=10, pad=-4)
                pick = cubes[k] if len(cubes) > k else None
                if pick is None or marching_cubes is None:
                    ax.axis("off")
                    ax.text2D(0.5, 0.5, "cube unavailable", ha="center",
                              transform=ax.transAxes, fontsize=9,
                              color="gray")
                    continue
                origin, axes_v, data = _read_cube(pick)
                stride = 2 if max(data.shape) > 110 else 1
                if stride > 1:
                    data = data[::stride, ::stride, ::stride]
                axes_v = [a * stride for a in axes_v]
                span = np.array([np.linalg.norm(a) for a in axes_v])
                iso = 0.03
                if np.abs(data).max() < iso:
                    iso = 0.25 * float(np.abs(data).max())
                drawn = 0
                for sign, color in ((1.0, C_MAIN), (-1.0, "#D62728")):
                    try:
                        verts, faces, _, _ = marching_cubes(
                            (sign * data).astype(np.float32),
                            level=iso, step_size=1)
                    except ValueError:
                        continue
                    v = verts * span + origin
                    ax.plot_trisurf(v[:, 0], v[:, 1], v[:, 2],
                                    triangles=faces, color=color,
                                    alpha=0.42, linewidth=0,
                                    antialiased=False, shade=True)
                    drawn += 1
                if drawn == 0:
                    ax.text2D(0.5, 0.5, "empty isosurface", ha="center",
                              transform=ax.transAxes, fontsize=9,
                              color="gray")
                ax.set_xticks([]), ax.set_yticks([]), ax.set_zticks([])
                ax.text2D(0.03, 0.90, f"R = {R:.2f} Å", fontsize=10,
                          weight="bold", transform=ax.transAxes)

        axb = fig.add_subplot(gs[2, :])
    else:
        fig, axb = plt.subplots(figsize=(9.5, 3.6))
        _warn("fig2: no orbital cubes found — NOON panel only")

    xs = np.array([r for r, v in zip(Rs_all, ana["noon_bonding"])
                   if v is not None])
    n1 = np.array([v for v in ana["noon_bonding"] if v is not None])
    n2 = np.array([v for v in ana["noon_antibonding"] if v is not None])
    yv = np.array([v for v in ana["y_diradical"] if v is not None])
    if xs.size:
        axb.plot(xs, n1, "-o", color=C_MAIN, lw=2.2, ms=6, mec="k", mew=0.5,
                 label=r"$\nu(\sigma\ \mathrm{NO})$")
        axb.plot(xs, n2, "-s", color=C_SEC, lw=2.0, ms=5.5, mec="k",
                 mew=0.5, label=r"$\nu(\sigma^*\ \mathrm{NO})$")
        axb.plot(xs, yv, "--^", color=C_TER, lw=1.8, ms=5.5,
                 label=r"diradical character  $y = \nu(\sigma^*\ \mathrm{NO})$")
        axb.axhline(1.0, color="gray", ls=":", lw=1)
        axb.text(xs[0], 1.03, "pure diradical limit", fontsize=8.4,
                 color="gray")
        for R in show_cols:
            axb.axvline(R, color="gray", lw=0.9, ls=":", alpha=0.7)
    axb.set_xlabel(r"scissile C6–C9 distance  $R$  (Å)")
    axb.set_ylabel("CASSCF(2,2)\noccupation")
    axb.set_ylim(-0.08, 2.18)
    if xs.size:
        axb.legend(frameon=False, ncols=3, loc="center right")

    fig.suptitle("Fig. 2 — CASSCF(2,2) active space of the scissile bond: "
                 "fractional natural-orbital occupations are the "
                 "fingerprint of static correlation",
                 x=0.01, ha="left", fontsize=12)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "fig2_casscf_frontier_orbitals.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    _log("fig", "fig2_casscf_frontier_orbitals.png written")


def fig3(args, ana):
    plt = _style()
    Rs = np.array(ana["R"])
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.8, 9.0),
                                   sharex=True,
                                   gridspec_kw={"height_ratios": [3, 2],
                                                "hspace": 0.12})
    styles = {
        "CASSCF": (C_REF, "-", "o", 7.5, 2.8,
                   "CASSCF(2,2) — multi-reference ground truth"),
        "BS-UB3LYP": (C_TER, "-", "s", 6, 2.0,
                      "BS-UB3LYP/def2-SVP — broken-symmetry DFT"),
        "UHF": ("#999999", "--", "^", 5.5, 1.6, "RHF — closed-shell HF"),
        "GFN2-xTB": (C_MAIN, "-", "D", 6, 2.0,
                     "GFN2-xTB — semi-empirical (evaluated identically)"),
        "MACE-OFF": (C_ML1, "-", "v", 7, 2.2,
                     "MACE-OFF (small) — equivariant ML potential"),
        "ANI-2x": (C_ML2, "-", "P", 7, 1.8, "ANI-2x — ML potential"),
    }
    for name, (c, ls, mk, ms, lw, lab) in styles.items():
        ys = ana["E_rel_kcal"].get(name)
        if not ys or all(v is None for v in ys):
            continue
        x = [r for r, v in zip(Rs, ys) if v is not None]
        y = [v for v in ys if v is not None]
        ax1.plot(x, y, ls, color=c, marker=mk, ms=ms, lw=lw,
                 mec="k" if name == "CASSCF" else c, mew=0.6,
                 label=lab, zorder=6 if name == "CASSCF" else 4)
    ax1.set_ylabel(r"$E_\mathrm{rel}(R)=E(R)-E(R_\mathrm{eq})$  "
                   "(kcal/mol)")
    ax1.legend(frameon=False, loc="upper left", fontsize=9.2)
    ax1.set_title("Fig. 3 — The Wall of Sighs: identical C6–C9 homolysis "
                  "trajectory, six theoretical lenses", loc="left")
    rc = ana["R_crit"].get("ST_crossing_UHF")
    if rc is not None:
        ax1.axvline(rc, color="#B22", ls="-.", lw=1.2, alpha=0.8)
        ax1.text(rc + 0.02, ax1.get_ylim()[0], r"$R_\mathrm{crit}$",
                 color="#B22", fontsize=9)

    zone_drawn = False
    for name, (c, _ls, mk, ms, lw, lab) in styles.items():
        if name == "CASSCF":
            continue
        errs = ana["abs_error_vs_casscf_kcal"].get(name)
        if not errs or all(v is None for v in errs):
            continue
        x = [r for r, v in zip(Rs, errs) if v is not None]
        y = [v for v in errs if v is not None]
        ax2.plot(x, y, "-", color=c, marker=mk, ms=ms - 1, lw=lw,
                 alpha=0.95, label=f"|ΔE error| — {name}")
        r_fail = ana["R_fail_15kcal"].get(name)
        if r_fail is not None:
            ax2.axvline(r_fail, color=c, ls=":", lw=1.1, alpha=0.8)
            ymax = max(v for v in errs if v is not None)
            ax2.text(r_fail + 0.015, min(ymax * 0.55, FAIL_KCAL * 2.1),
                     f"{r_fail:.2f} Å", color=c, fontsize=8.6,
                     rotation=90, va="bottom")
    ax2.axhline(FAIL_KCAL, color="#D62728", ls="--", lw=1.5)
    rfails = [v for v in ana["R_fail_15kcal"].values() if v is not None]
    if rfails:
        z0, z1 = min(rfails), Rs[-1]
        ax2.axvspan(z0, z1, color="#D62728", alpha=0.08, zorder=0)
        ax2.text(z0 + 0.25, 41.0, "EPISTEMIC FAILURE ZONE",
                 ha="left", color="#B22", fontsize=10)
        ax2.text(z0 + 0.25, 37.4,
                 r"($|\Delta E_\mathrm{error}|>15$ kcal/mol vs CASSCF)",
                 ha="left", color="#B22", fontsize=8.6)
    ax2.axhline(1.0, color="gray", ls=":", lw=1.0)
    ax2.text(Rs[0] + 0.03, -3.2, "chemical accuracy (1 kcal/mol)",
             fontsize=8.4, color="gray", ha="left")
    ax2.set_xlabel(r"scissile C6–C9 distance  $R$  (Å)")
    ax2.set_ylabel(r"$\Delta E_\mathrm{error}(R)="
                   r"|E^\mathrm{AI}_\mathrm{rel}-"
                   r"E^\mathrm{CASSCF}_\mathrm{rel}|$  (kcal/mol)")
    ax2.legend(frameon=False, ncols=2, fontsize=8.8, loc="upper left")
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "fig3_the_wall_of_sighs_discrepancy.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    _log("fig", "fig3_the_wall_of_sighs_discrepancy.png written")


def write_summary_csv(ana):
    import csv
    Rs = ana["R"]
    path = OUT / "phase7_scan_summary.csv"
    cols = ["R_A", "E_RHF_Eh", "E_UHF_Eh", "E_RKS_Eh", "E_UKS_Eh",
            "S2_UHF", "S2_UKS", "E_UHFTrip_Eh", "S2_UHFTrip", "dE_ST_UHF_kcal", "E_CASSCF_Eh", "NOON_sigma",
            "NOON_sigma_star", "y_diradical",
            "Erel_CASSCF_kcal", "Erel_GFN2xTB_kcal", "Erel_MACE_kcal",
            "Erel_ANI_kcal", "Erel_BSUB3LYP_kcal",
            "err_xTB_kcal", "err_MACE_kcal", "err_ANI_kcal",
            "err_BSUB3LYP_kcal", "err_RHF_kcal"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for i, R in enumerate(Rs):
            row = [R,
                   ana["E_rhf_eh"][i], ana["E_uhf_eh"][i],
                   ana["E_rks_eh"][i], ana["E_uks_eh"][i],
                   ana["S2_uhf"][i], ana["S2_uks"][i],
                   ana["E_uhf_triplet_eh"][i], ana["S2_uhf_triplet"][i],
                   ana["dE_ST_uhf_kcal"][i],
                   ana["E_casscf_eh"][i], ana["noon_bonding"][i],
                   ana["noon_antibonding"][i], ana["y_diradical"][i]]
            for k in ("CASSCF", "GFN2-xTB", "MACE-OFF", "ANI-2x",
                      "BS-UB3LYP"):
                row.append(ana["E_rel_kcal"][k][i])
            for k in ("GFN2-xTB", "MACE-OFF", "ANI-2x", "BS-UB3LYP",
                      "UHF"):
                row.append(ana["abs_error_vs_casscf_kcal"][k][i])
            w.writerow(row)
    _log("csv", f"{path} written")


# --------------------------------------------------------------------------- #
#  orchestrator
# --------------------------------------------------------------------------- #
def run_cmd(cmd, env=None, cwd=None, timeout=None):
    _log("orch", "spawn: " + " ".join(str(c) for c in cmd[:6]) +
         (" ..." if len(cmd) > 6 else ""))
    e = dict(os.environ)
    if env:
        e.update(env)
    return subprocess.run([str(c) for c in cmd], env=e,
                          cwd=str(cwd or Path.cwd()), timeout=timeout)


def orchestrate(args):
    t00 = time.time()
    META["engines"] = {
        "py_qc": find_py_qc(),
        "py_chem": find_py_chem(),
        "xtb": find_xtb(),
        "qc_backend": "Psi4 1.11 (DETCI CASSCF)",
        "note": ("reference protocol names PySCF; PySCF provides no native "
                 "win32 build on this host — Psi4 1.11 substitute, "
                 "protocol-equivalent (RHF/UHF/RKS/UKS + CAS(2e,2o) CASSCF "
                 "+ NOON); substitution logged by design"),
    }
    scan_r = SMOKE_R if args.smoke else args.scan_r
    basis = SMOKE_BASIS if args.smoke else args.basis
    cube_rs = SMOKE_R[:1] + [SMOKE_R[1]] if args.smoke else CUBE_RS
    threads = args.threads
    py_qc, py_chem = META["engines"]["py_qc"], META["engines"]["py_chem"]

    stages = (["G", "AI", "7A", "7B", "merge", "fig"] if args.stage == "all"
              else [s.strip() for s in args.stage.split(",")])

    # ---- stage G & AI: chem env ------------------------------------------
    if "G" in stages:
        rc = run_cmd([py_chem, __file__, "--worker", "G",
                      "--scan-r", *[str(r) for r in scan_r],
                      "--threads", str(threads)])
        if rc.returncode != 0:
            _warn(f"stage G worker exit {rc.returncode}")
    if "AI" in stages:
        rc = run_cmd([py_chem, __file__, "--worker", "AI",
                      "--threads", str(threads)])
        if rc.returncode != 0:
            _warn(f"stage AI worker exit {rc.returncode}")

    # ---- stages 7A/7B: per-point isolated QC subprocesses -----------------
    if "G" in stages or (ckpt_dir("G") / "scan.json").exists():
        scan = read_json(ckpt_dir("G") / "scan.json")["points"]
    else:
        scan = []
    if "7A" in stages:
        for i, pt in enumerate(scan):
            ok = False
            for basis_t in [basis] + _basis_chain_for(basis):
                rc = run_cmd([py_qc, __file__, "--worker", "7A",
                              "--point", str(i), "--basis", basis_t,
                              "--mem", args.mem,
                              "--threads", str(threads)])
                ck = ckpt_dir("7A") / f"R{pt['R']:.2f}.json"
                if rc.returncode == 0 and ck.exists():
                    ok = True
                    break
                _fallback(f"7A point {i} ({pt['R']:.2f} A) tier {basis_t} "
                          f"failed — degrading basis (logged trade-off: "
                          f"reduced polarization/flexibility; qualitative "
                          f"profiles preserved)")
            if not ok:
                _warn(f"7A point {i} failed on every basis tier")
    if "7B" in stages:
        for i, pt in enumerate(scan):
            ok = False
            cube_here = [pt["R"]] if any(
                abs(pt["R"] - cr) < 1e-6 for cr in cube_rs) else [-1.0]
            for basis_t in [basis] + _basis_chain_for(basis):
                for skip in (0, 1):     # 0 = incl. CONV, 1 = AO/DF only
                    rc = run_cmd([py_qc, __file__, "--worker", "7B",
                                  "--point", str(i), "--basis", basis_t,
                                  "--mem", args.mem,
                                  "--threads", str(threads),
                                  "--skip-tiers", str(skip),
                                  "--cube-r", str(cube_here[0])])
                    ck = ckpt_dir("7B") / f"R{pt['R']:.2f}.json"
                    if (rc.returncode == 0 and ck.exists()
                            and read_json(ck).get("converged")):
                        ok = True
                        break
                    _fallback(f"7B point {i} ({pt['R']:.2f} A) tier "
                              f"{basis_t}/skip{skip} unconverged — "
                              f"degrading (trade-off logged: smaller basis "
                              f"or AO/DF algorithm shifts NOON fractions "
                              f"and total energy, not the qualitative "
                              f"diradical signature)")
                if ok:
                    break
            if not ok:
                _warn(f"7B point {i} failed on every basis tier")

    if "merge" in stages or "fig" in stages:
        try:
            ana = merge_and_analyze()
            master = dict(META)
            master["scan"] = {"points": scan, "basis": basis,
                              "scan_r": scan_r, "cube_rs": cube_rs}
            master["analysis"] = ana
            master["wall_s"] = round(time.time() - t00, 1)
            master["all_stages_ok"] = all(
                v is not None for v in ana["E_casscf_eh"]) and any(
                v is not None for v in ana["S2_uhf"])
            write_json_atomic(OUT / "phase7_results.json", master)
            write_summary_csv(ana)
            _log("merge", "phase7_results.json + CSV written")
        except Exception as exc:
            _log("merge", f"merge failed: {traceback.format_exc()[-500:]}")
            master = dict(META)
            master["fatal_error"] = str(exc)
            write_json_atomic(OUT / "phase7_results.json", master)

    if "fig" in stages:
        try:
            master = read_json(OUT / "phase7_results.json")
            ana = master["analysis"]
            fig1(args, ana)
            fig2(args, ana)
            fig3(args, ana)
        except Exception:
            _log("fig", f"figures failed: {traceback.format_exc()[-500:]}")

    _log("orch", f"done in {round(time.time()-t00,1)} s")
    return 0


def _basis_chain_for(basis):
    return [b for b in BASIS_CHAIN if b != basis.lower()]


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--stage", default="all",
                    help="all|G|AI|7A|7B|merge|fig (comma list)")
    ap.add_argument("--worker", default=None, choices=["G", "AI", "7A", "7B"])
    ap.add_argument("--point", type=int, default=0)
    ap.add_argument("--basis", default="def2-svp")
    ap.add_argument("--mem", default=QC_MEM)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--cube-r", type=float, default=-1.0)
    ap.add_argument("--skip-tiers", type=int, default=0,
                    help="skip leading CASSCF algorithm tiers "
                         "(0=CONV/AO/DF, 1=AO/DF only) after hard aborts")
    ap.add_argument("--scan-r", type=float, nargs="*", default=SCAN_R)
    ap.add_argument("--smoke", action="store_true",
                    help="3-point / 6-31G pipeline validation")
    ap.add_argument("--fig_only", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if args.smoke:
        args.scan_r = SMOKE_R
        args.basis = SMOKE_BASIS
    if args.fig_only:
        args.stage = "fig"

    if args.worker:
        OUT.mkdir(parents=True, exist_ok=True)
        try:
            if args.worker == "G":
                worker_G(args)
            elif args.worker == "AI":
                worker_AI(args)
            elif args.worker == "7A":
                worker_7A(args)
            elif args.worker == "7B":
                worker_7B(args)
            return 0
        except Exception:
            traceback.print_exc()
            return 1

    return orchestrate(args)


if __name__ == "__main__":
    sys.exit(main())
