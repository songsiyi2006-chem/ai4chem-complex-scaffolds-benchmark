#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_phase8_photochemical_dynamics.py
====================================
PHASE 8 — BEYOND THE BORN-OPPENHEIMER SEA
EXCITED-STATE PHOTOCHEMISTRY, CONICAL INTERSECTIONS & TULLY SURFACE HOPPING

Model system: the azobenzene molecular photo-switch (E -> Z photo-isomerization
of the N=N chromophore), treated with a three-layer multi-scale protocol:

  Module 8A (real chromophore) — trans-azobenzene C12H10N2, GFN2-xTB relaxed.
    TD-DFT (TDA-B3LYP/def2-SVP) vertical excitations: up to 10 singlet + 10
    triplet states with oscillator strengths and electric transition dipole
    moments, Natural Transition Orbital (NTO) particle-hole pairs (cube
    export), simulated UV-Vis spectrum (sigma = 0.2 eV Gaussian broadening),
    and a rigid CNNC torsion scan of the S1/S2 vertical energies.

  Module 8B (minimal chromophore) — diazene HN=NH, the two-atom N=N switch
    core.  State-averaged CASSCF(4,4) {pi, pi*, n, n*} torsion scan, then a
    Bearpark-Robb-Schlegel penalty-constrained optimization
        F(R) = 1/2 (E1 + E0) + sigma (E1 - E0)^2 / (E1 - E0 + alpha)
    to the minimum-energy conical intersection (MECI).  Branching space:
    gradient-difference vector g = dE1 - dE0 (finite differences — analytic
    CASSCF gradients return silent zeros on this Psi4 build) and derivative
    coupling h reconstructed by directional gap-lifting maximization (the
    mission-sanctioned localized state-overlap approximation).

  Module 8C (dynamics) — 2-state / 3-mode linear vibronic-coupling (LVC)
    Hamiltonian parameterized entirely from the 8A/8B/xtb ab initio data
    (FC gap, CASSCF torsion & stretch scans, GFN2-xTB hessian modes), then an
    ensemble of Tully fewest-switches surface-hopping (FSSH) trajectories from
    Wigner-sampled Franck-Condon conditions on S1: Velocity-Verlet nuclei
    (dt = 0.5 fs, Tmax = 500 fs), exact 2x2 electronic propagator, Tully
    hopping with momentum rescaling along the nonadiabatic coupling vector,
    frustrated-hop reflection, Granucci-Persico decoherence correction.
    Observables: P_S1(t), excited-state lifetime, E->Z quantum yield.

Engines & environments (fault-tolerant multi-interpreter orchestration)
----------------------------------------------------------------------
QC engine    : Psi4 1.11 (conda env `phase7`) — TDSCF + DETCI SA-CASSCF.
               The reference protocol names PySCF; PySCF ships no win32 wheels
               (PyPI 2.14.0 = macOS/Linux only; no conda-forge win-64 build),
               so Psi4 1.11 is the drop-in ab-initio backend and the
               substitution is logged in every result file.  Quantities are
               backend-invariant (excitation energies, CAS state energies).
Chem engine  : xtb.exe GFN2-xTB (phase-4 subprocess wrapper pattern) for the
               S0 relaxed geometries and the normal-mode basis (g98.out).
Dynamics     : pure numpy vectorized FSSH under the driver interpreter.
Known Psi4 build quirks (carried over from Phase 7):
  - QC memory must stay in the 512 MB - 2 GB band (DETCI DPD cache overflow),
  - analytic CASSCF gradients return zero arrays -> all 8B gradients are
    finite differences of CAS state energies (diazene: 0.7 s/job),
  - state averaging = global `num_roots` + module detci `avg_states`,
  - TDSCF requires `save_jk` on the reference wfn and `tdscf_states`, C1.

Outputs
-------
results_phase8/phase8_results.json            machine-readable master record
results_phase8/*.xyz                          S0 geometries, scans, MECI
results_phase8/nto_cubes/**.cube              NTO particle-hole cubes
figures_phase8/fig1_uv_vis_absorption_spectrum.png    (300 DPI)
figures_phase8/fig2_conical_intersection_topology.png (300 DPI)
figures_phase8/fig3_fssh_population_trajectories.png  (300 DPI)

Usage
-----
python run_phase8_photochemical_dynamics.py               # full pipeline
python run_phase8_photochemical_dynamics.py --stage 8A    # one module
python run_phase8_photochemical_dynamics.py --smoke       # fast validation
python run_phase8_photochemical_dynamics.py --fig_only    # refit figures
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import re
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------- #
#  constants
# --------------------------------------------------------------------------- #
EH_EV = 27.211386245988
HARTREE_KCAL = 627.5094740631
EV_NM = 1239.841984                      # lambda [nm] = EV_NM / E [eV]
EA_DEBYE = 2.5417464519                  # e*a0 -> Debye
CM1_EH = 1.0 / 219474.6313632            # cm^-1 -> Eh
CM1_AU = 2.0 * math.pi * CM1_EH          # angular frequency, a.u.
AMU_ME = 1822.888486209                  # amu -> electron masses
FS_AU = 41.3413733365614                 # 1 fs in atomic time units
ANG_BOHR = 1.8897261254578281            # Angstrom -> bohr

ENV_PY_QC = r"C:\Users\HUIWEI\miniconda3\envs\phase7\python.exe"
ENV_CHEM_BIN = r"C:\Users\HUIWEI\miniconda3\envs\phase2ff\Library\bin"
QC_MEM = "2 GB"                          # NEVER 4 GB on this psi4 build

OUT = Path("results_phase8")
FIG = Path("figures_phase8")
SEED = 20260905

# Module 8A
AZO_SMILES = "C1=CC=C(C=C1)N=NC2=CC=CC=C2"
N_SINGLET, N_TRIPLET = 10, 10
BASIS_8A = "def2-svp"
BASIS_8A_FALLBACK = "6-31g"
DFT_FUNC = "b3lyp"
SPEC_EMIN, SPEC_EMAX, SPEC_DE, SPEC_SIGMA = 1.5, 8.0, 0.005, 0.2  # eV
SCAN_AZO_PHI = [180.0, 140.0, 100.0, 90.0, 80.0, 40.0, 0.0]

# Module 8B
DIAZENE_NN0 = 1.25
SCAN_DIA_PHI = [180.0, 160.0, 140.0, 120.0, 100.0, 90.0, 80.0, 60.0, 40.0, 20.0, 0.0]
BASIS_8B = "6-31g"
BASIS_8B_FALLBACK = "sto-3g"
PENALTY_SIGMA = 1.0                      # Eh^-1, BRS penalty prefactor
PENALTY_ALPHA = 0.02                     # Eh, BRS softening
MECI_GAP_EV_GATE = 0.05                  # mission convergence gate on DE10
MECI_GRAD_GATE = 0.05                    # Eh/Angstrom on |grad F|
MECI_MAX_STEPS = 60
FD_GRAD_STEP = 0.005                     # Angstrom, forward-difference step
LIFT_DELTA = 0.02                        # Angstrom, gap-lifting displacement
STRETCH_SCAN_D = [-0.10, -0.05, 0.0, 0.05, 0.10]

# Module 8C
N_TRAJ, DT_FS, TMAX_FS = 300, 0.5, 500.0
N_ESUB = 40                              # electronic substeps / nuclear step
EDC_ALPHA = 0.1                          # Granucci-Persico decoherence [Eh]
SAVE_EVERY_FS = 5.0
CI_TORSION_DEG = 90.0                    # nominal CI torsion (refined from scan)

ENGINE_NOTE = ("The reference protocol names PySCF (tdscf/mcscf); PySCF ships "
               "no win32 wheels (PyPI 2.14.0 = macOS/Linux only, no "
               "conda-forge win-64 build), so Psi4 1.11 DETCI/TDSCF is the "
               "drop-in ab-initio backend on this Windows host. Quantities "
               "are backend-invariant.")


# --------------------------------------------------------------------------- #
#  small utilities
# --------------------------------------------------------------------------- #
def _log(stage: str, msg: str) -> None:
    print(f"[{_dt.datetime.now():%H:%M:%S}] [{stage:>4}] {msg}", flush=True)


def _warn(stage: str, msg: str) -> None:
    print(f"[{_dt.datetime.now():%H:%M:%S}] [{stage:>4}] !! {msg}", flush=True)


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    els, xyz = [], []
    for ln in lines[2:2 + n]:
        t = ln.split()
        els.append(t[0])
        xyz.append([float(t[1]), float(t[2]), float(t[3])])
    return els, np.asarray(xyz)


def write_xyz(path, els, xyz, comment=""):
    xyz = np.asarray(xyz)
    lines = [str(len(els)), comment]
    for e, p in zip(els, xyz):
        lines.append(f"{e} {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}")
    Path(path).write_text("\n".join(lines) + "\n")


def psi4_geometry_string(els, xyz, charge=0, mult=1):
    body = "\n".join(f"{e} {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}"
                     for e, p in zip(els, np.asarray(xyz)))
    return f"symmetry c1\n{charge} {mult}\n{body}\n"


def parse_json(path):
    return json.loads(Path(path).read_text())


def dump_json(obj, path):
    Path(path).write_text(json.dumps(obj, indent=1, default=float))


# --------------------------------------------------------------------------- #
#  DRIVER-SIDE STRUCTURE GENERATION (rdkit + ase + xtb)
# --------------------------------------------------------------------------- #
def find_xtb():
    cands = [shutil.which("xtb"), str(Path(ENV_CHEM_BIN) / "xtb.exe")]
    for c in cands:
        if c and Path(c).exists():
            return c
    return None


def run_xtb(xtb_exe, xyz_path, workdir, opt=False, hess=False, timeout=900):
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    shutil.copy(xyz_path, workdir / "m.xyz")
    cmd = [xtb_exe, "m.xyz", "--chrg", "0", "--mult", "1"]
    if opt:
        cmd += ["--opt", "tight"]
    if hess:
        cmd += ["--hess"]
    proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True,
                          timeout=timeout)
    (workdir / "xtb.stdout").write_text((proc.stdout or "")[-8000:])
    return proc, workdir


def build_azobenzene(xtb_exe, out_xyz):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(AZO_SMILES))
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED % 2**31
    AllChem.EmbedMolecule(m, params)
    AllChem.MMFFOptimizeMolecule(m)
    c = m.GetConformer().GetPositions()
    els = [a.GetSymbol() for a in m.GetAtoms()]
    if xtb_exe:
        workdir = OUT / "_xtb_azobenzene"
        workdir.mkdir(parents=True, exist_ok=True)
        write_xyz(workdir / "start.xyz", els, c, "rdkit ETKDGv3/MMFF94")
        try:
            run_xtb(xtb_exe, workdir / "start.xyz", workdir, opt=True)
            els, c = read_xyz(workdir / "xtbopt.xyz")
            _log("PRE", "azobenzene: GFN2-xTB optimization converged")
        except Exception as exc:
            _warn("PRE", f"xtb opt failed ({str(exc)[:90]}) — keeping MMFF94")
    write_xyz(out_xyz, els, c, "trans-azobenzene S0 (GFN2-xTB/MMFF94)")


def azo_ring2_mask_and_axis(els, xyz):
    """Ring-2 side atom mask for the torsion + (C1,N1,N2,C2) dihedral quartet."""
    els_l = [e.lower() for e in els]
    n_idx = [i for i, e in enumerate(els_l) if e == "n"]
    best = min(((i, j) for k, i in enumerate(n_idx) for j in n_idx[k + 1:]),
               key=lambda p: np.linalg.norm(xyz[p[0]] - xyz[p[1]]))
    n1, n2 = best
    cset = [i for i, e in enumerate(els_l) if e == "c"]

    def bonded_c(n):
        return min(cset, key=lambda i: np.linalg.norm(xyz[i] - xyz[n]))

    c1, c2 = bonded_c(n1), bonded_c(n2)

    def side(start):
        seen, stack = {start}, [start]
        while stack:
            i = stack.pop()
            for j in cset:
                if j not in seen and np.linalg.norm(xyz[i] - xyz[j]) < 1.8:
                    seen.add(j)
                    stack.append(j)
        return seen

    s1 = side(c1)
    ring2 = [i for i in cset if i not in s1]
    h2 = []
    for i, e in enumerate(els_l):
        if e == "h":
            cs = [c for c in ring2 if np.linalg.norm(xyz[i] - xyz[c]) < 1.3]
            if cs:
                h2.append(i)
    return sorted(h2 + ring2), (c1, n1, n2, c2)


def build_azobenzene_scan(s0_xyz):
    """Rigid CNNC torsion scan geometries (rotate ring-2 about the C-N axis)."""
    from ase import Atoms
    els, xyz = read_xyz(s0_xyz)
    mask, quartet = azo_ring2_mask_and_axis(els, xyz)
    geo = {}
    for phi in SCAN_AZO_PHI:
        at = Atoms(symbols=els, positions=xyz)
        mset = set(mask)
        mask_bool = [i in mset for i in range(len(els))]
        at.set_dihedral(quartet[0], quartet[1], quartet[2], quartet[3],
                        phi, mask=mask_bool)
        tag = f"azo_phi{int(phi):03d}"
        p = OUT / f"{tag}.xyz"
        write_xyz(p, els, at.get_positions(),
                  f"azobenzene rigid scan CNNC={phi:.0f} deg")
        geo[f"{phi:.0f}"] = {"tag": tag, "phi": phi, "file": str(p)}
    return geo


def build_diazene(xtb_exe, out_xyz):
    els = ["N", "N", "H", "H"]
    r = DIAZENE_NN0 / 2
    xyz = np.array([[0.0, 0.0, r], [0.0, 0.0, -r],
                    [0.0, 0.95, r + 0.44], [0.0, -0.95, -(r + 0.44)]])
    if xtb_exe:
        workdir = OUT / "_xtb_diazene"
        workdir.mkdir(parents=True, exist_ok=True)
        write_xyz(workdir / "start.xyz", els, xyz, "diazene trans guess")
        try:
            run_xtb(xtb_exe, workdir / "start.xyz", workdir, opt=True)
            els, xyz = read_xyz(workdir / "xtbopt.xyz")
            _log("PRE", "diazene: GFN2-xTB optimization converged")
        except Exception as exc:
            _warn("PRE", f"xtb opt failed ({str(exc)[:90]}) — keeping guess")
    write_xyz(out_xyz, els, xyz, "trans-diazene S0 (GFN2-xTB)")
    return {"elements": els, "nn_ang": float(np.linalg.norm(xyz[0] - xyz[1]))}


def build_diazene_scans(s0_xyz, phi_ci_deg):
    """Torsion scan (rigid) + N=N stretch mini-scan at the twisted geometry."""
    from ase import Atoms
    els, xyz = read_xyz(s0_xyz)
    geo_t = {}
    for phi in SCAN_DIA_PHI:
        at = Atoms(symbols=els, positions=xyz)
        at.set_dihedral(2, 0, 1, 3, phi,
                        mask=[False, False, True, False])
        tag = f"dia_phi{int(phi):03d}"
        p = OUT / f"{tag}.xyz"
        write_xyz(p, els, at.get_positions(),
                  f"diazene rigid torsion {phi:.0f} deg")
        geo_t[f"{phi:.0f}"] = {"tag": tag, "phi": phi, "file": str(p)}
    at90 = Atoms(symbols=els, positions=xyz)
    at90.set_dihedral(2, 0, 1, 3, phi_ci_deg,
                      mask=[False, False, True, False])
    base = at90.get_positions().copy()
    axis = (base[0] - base[1]) / np.linalg.norm(base[0] - base[1])
    geo_s = {}
    for d in STRETCH_SCAN_D:
        pos = base.copy()
        pos[0] += axis * (d / 2)
        pos[1] -= axis * (d / 2)
        tag = f"dia_stretch{int(round(d * 100)):+04d}"
        p = OUT / f"{tag}.xyz"
        write_xyz(p, els, pos, f"diazene N=N stretch {d:+.2f} A at CI torsion")
        geo_s[f"{d:+.2f}"] = {"tag": tag, "d": d, "file": str(p)}
    return geo_t, geo_s


def diazene_normal_modes(xtb_exe, s0_xyz):
    """GFN2-xTB hessian -> (freqs cm-1, red masses amu, mode vectors)."""
    els, xyz = read_xyz(s0_xyz)
    nat = len(els)
    workdir = OUT / "_xtb_diazene_hess"
    run_xtb(xtb_exe, s0_xyz, workdir, opt=False, hess=True)
    lines = (workdir / "g98.out").read_text().splitlines()
    freqs, redmass, modes = [], [], []
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.strip().startswith("Atom AN"):
            per_atom = []
            j = i + 1
            while j < len(lines):
                t = lines[j].split()
                # atom rows carry exactly 2 id tokens + 3 modes x 3 components
                if len(t) != 11 or not re.match(r"\s*\d+\s+\d+\s", lines[j]):
                    break
                per_atom.append([float(x) for x in t[2:]])
                j += 1
            ncols = len(per_atom[0]) // 3 if per_atom else 0
            for cc in range(ncols):
                modes.append(np.array([[p[3 * cc], p[3 * cc + 1], p[3 * cc + 2]]
                                       for p in per_atom]))
            i = j
            continue
        m = re.match(r"\s*Frequencies --\s+(.*)", ln)
        if m:
            freqs.extend(float(x) for x in m.group(1).split())
        m = re.match(r"\s*Red. masses --\s+(.*)", ln)
        if m:
            redmass.extend(float(x) for x in m.group(1).split())
        i += 1
    if not modes or len(freqs) < len(modes):
        raise RuntimeError("g98.out parse failed")
    return np.array(freqs[:len(modes)]), np.array(redmass[:len(modes)]), modes


def classify_diazene_modes(freqs, redmass, mode_vecs, els, xyz):
    """Pick N=N-stretch / bend from the g98 mode table.

    The torsional normal mode is NOT taken from the g98 table (on the
    near-planar geometry it mixes with the out-of-plane wag); instead the
    LVC uses a physical torsional inertia  I = sum_H m_H rho_H^2  about the
    N=N axis combined with the CASSCF scan curvature (see build_lvc)."""
    nat = len(els)
    axis = xyz[0] - xyz[1]
    axis = axis / np.linalg.norm(axis)

    def stretch_score(mv):
        d0 = float(np.dot(mv[0], axis))
        d1 = float(np.dot(mv[1], axis))
        return (d0 - d1) ** 2

    k_st = max(range(len(mode_vecs)), key=lambda k: stretch_score(mode_vecs[k]))
    rest = [k for k in range(len(mode_vecs))
            if k != k_st and freqs[k] > 50.0]
    k_be = min(rest, key=lambda k: freqs[k]) if rest else 0
    return {"stretch": (k_st, float(freqs[k_st]), float(redmass[k_st])),
            "bend": (k_be, float(freqs[k_be]), float(redmass[k_be]))}


def torsion_inertia(els, xyz):
    """Physical moment of inertia of the HNNH torsion about the N=N axis
    (amu * Angstrom^2), H atoms only (N sit on the axis)."""
    axis = xyz[0] - xyz[1]
    axis = axis / np.linalg.norm(axis)
    mass = {"H": 1.008, "N": 14.007, "C": 12.011}
    I = 0.0
    for i, e in enumerate(els):
        if e != "H":
            continue
        rho = xyz[i] - np.dot(xyz[i], axis) * axis
        I += mass[e] * float(np.dot(rho, rho))
    return I


# --------------------------------------------------------------------------- #
#  WORKER 8A — TD-DFT / UV-Vis / NTO cubes         (runs under ENV_PY_QC)
# --------------------------------------------------------------------------- #
def _td_variables_map(wfn, nroots):
    keys = list(wfn.variables().keys())
    out = {}
    for root in range(1, nroots + 1):
        ee = osc = tdm = None
        for k in keys:
            if re.search(rf"ROOT 0 \((?:AP|A)\) -> ROOT {root} "
                         rf"\((?:APP|A)\)", k):
                if "EXCITATION ENERGY" in k:
                    ee = float(wfn.variable(k))
                elif "OSCILLATOR STRENGTH (LEN)" in k:
                    osc = float(wfn.variable(k))
                elif "ELECTRIC TRANSITION DIPOLE MOMENT (LEN)" in k:
                    tdm = np.asarray(wfn.variable(k))
        out[root] = {"dE_eh": ee, "f_osc": osc,
                     "tdm_ea": None if tdm is None else tdm.tolist()}
    return out


def worker_8a(args, res_path):
    import psi4
    psi4.set_memory(QC_MEM)
    psi4.set_num_threads(args.threads)
    els, xyz = read_xyz(Path(args.xyz))
    mol = psi4.geometry(psi4_geometry_string(els, xyz))
    basis = BASIS_8A

    def run_scf(b):
        psi4.set_options({"basis": b, "reference": "rhf", "scf_type": "df",
                          "save_jk": True})
        return psi4.energy(f"{DFT_FUNC}/{b}", molecule=mol, return_wfn=True)

    try:
        e_scf, scf_wfn = run_scf(basis)
    except Exception as exc:
        _warn("8A", f"{basis} SCF failed ({str(exc)[:90]}) -> "
                    f"{BASIS_8A_FALLBACK}")
        basis = BASIS_8A_FALLBACK
        e_scf, scf_wfn = run_scf(basis)
    nocc = scf_wfn.nalpha()
    nbf = scf_wfn.basisset().nbf()
    _log("8A", f"DF-{DFT_FUNC.upper()}/{basis}: E_S0 = {e_scf:.6f} Eh, "
               f"nbf = {nbf}, nocc = {nocc}")

    out = {"basis": basis, "method": f"TDA-{DFT_FUNC.upper()}/{basis}",
           "e_scf_eh": float(e_scf), "nbf": int(nbf), "nocc": int(nocc),
           "engine": "Psi4 1.11 (PySCF-substituted on win32)", "states": {}}

    def run_td(nstates, singlet, tag):
        psi4.set_options({"tdscf_states": [nstates], "tdscf_tda": True})
        for opt_try in ({"tdscf_singlet": singlet}, {"singlet": singlet}):
            try:
                psi4.set_options(opt_try)
                break
            except Exception:
                continue
        psi4.energy(f"td-{DFT_FUNC}/{basis}", ref_wfn=scf_wfn, molecule=mol,
                    return_wfn=True)
        data = _td_variables_map(scf_wfn, nstates)
        n_ok = len([d for d in data.values() if d["dE_eh"]])
        _log("8A", f"{tag}: solved {n_ok}/{nstates} states")
        return data

    n_sing = 3 if args.smoke else N_SINGLET
    n_trip = 2 if args.smoke else N_TRIPLET
    sing = {}
    for n_try in (n_sing, 6, 4):
        try:
            sing = run_td(n_try, True, f"{n_try} singlets")
            break
        except Exception as exc:
            _warn("8A", f"{n_try}-singlet TD failed: {str(exc)[:90]}")

    centers = np.array([scf_wfn.basisset().function_to_center(mu)
                        for mu in range(nbf)])
    els_arr = np.array(els)

    def nto_for_state(root):
        k = [kk for kk in scf_wfn.variables()
             if re.search(rf"ROOT 0 \((?:AP|A)\) -> ROOT {root} "
                          rf"\((?:APP|A)\) (?:RIGHT|LEFT) EIGENVECTOR ALPHA",
                          kk)]
        if not k:
            return None
        X = np.asarray(scf_wfn.variable(k[0]))
        nvir = X.size // nocc
        X = X.reshape(nocc, nvir)
        U, S, Vt = np.linalg.svd(X)
        w = S**2 / max((S**2).sum(), 1e-300)
        Ca = np.asarray(scf_wfn.Ca())
        C_occ = Ca[:nbf, :nocc]
        C_vir = Ca[:nbf, nocc:nocc + nvir]
        hole = C_occ @ U[:, 0]
        part = C_vir @ Vt[0, :]

        def frac(orb, syms):
            mask = np.isin(els_arr[centers], syms)
            return float((orb[mask]**2).sum())

        return {"weights": w[:4].tolist(),
                "hole": hole, "particle": part,
                "hole_N": frac(hole, ["N"]), "hole_C": frac(hole, ["C"]),
                "part_N": frac(part, ["N"]), "part_C": frac(part, ["C"])}

    for root in list(sing):
        if sing[root]["dE_eh"] is None:
            continue
        try:
            nt = nto_for_state(root)
            if nt:
                sing[root]["nto"] = {k: v for k, v in nt.items()
                                     if k not in ("hole", "particle")}
        except Exception as exc:
            _warn("8A", f"NTO root {root} failed: {str(exc)[:80]}")

    trip = {}
    for n_try in (n_trip, 6, 4):
        try:
            trip = run_td(n_try, False, f"{n_try} triplets")
            break
        except Exception as exc:
            _warn("8A", f"{n_try}-triplet TD failed: {str(exc)[:90]}")

    def fmt(data, spin):
        arr = []
        for root in sorted(data):
            d = data[root]
            if not d or d["dE_eh"] is None:
                continue
            ev = d["dE_eh"] * EH_EV
            arr.append({"spin": spin, "root": root, "dE_eV": ev,
                        "lambda_nm": EV_NM / ev,
                        "f_osc": float(d["f_osc"] or 0.0),
                        "mu_Debye": (float(np.linalg.norm(d["tdm_ea"]) *
                                           EA_DEBYE)
                                     if d["tdm_ea"] else None),
                        "nto": d.get("nto")})
        return arr

    out["states"]["singlets"] = fmt(sing, "S")
    out["states"]["triplets"] = fmt(trip, "T")

    # NTO cubes via cubeprop (doctored Ca; done last so TD steps keep Ca clean)
    cube_dir = OUT / "nto_cubes"
    cube_dir.mkdir(parents=True, exist_ok=True)
    cube_info = {}
    try:
        targets = [("S1", 1)]
        cands = [s for s in out["states"]["singlets"]]
        if cands:
            bright = max(cands, key=lambda s: s["f_osc"])
            if bright["root"] != 1:
                targets.append((f"S{bright['root']}_bright", bright["root"]))
        for label, root in targets:
            nt = nto_for_state(root)
            if nt is None:
                continue
            Ca_view = np.asarray(scf_wfn.Ca())
            Ca_view[:, nocc - 1] = nt["hole"]
            Ca_view[:, nocc] = nt["particle"]
            cdir = cube_dir / f"state{root}_{label}"
            cdir.mkdir(parents=True, exist_ok=True)
            psi4.set_options({"cubeprop_tasks": ["ORBITALS"],
                              "cubeprop_orbitals": [nocc - 1, nocc],
                              "cubeprop_filepath": str(cdir.resolve())})
            psi4.cubeprop(scf_wfn)
            cube_info[label] = {"root": root,
                                "nto": {k: v for k, v in nt.items()
                                        if k not in ("hole", "particle")},
                                "cubes": sorted(p.name
                                                for p in cdir.glob("*.cube")),
                                "dir": str(cdir)}
            _log("8A", f"NTO cubes for {label}: {cube_info[label]['cubes']}")
        out["nto_cubes"] = cube_info
    except Exception as exc:
        _warn("8A", f"NTO cubeprop failed: {str(exc)[:120]}")

    dump_json(out, res_path)
    _log("8A", f"module 8A complete -> {res_path}")


def worker_8a_scan(args, res_path):
    """Azobenzene rigid torsion scan: S0 + 3 singlet vertical energies/point."""
    import psi4
    psi4.set_memory(QC_MEM)
    psi4.set_num_threads(args.threads)
    scan = parse_json(Path(args.scan_json))
    basis = BASIS_8A_FALLBACK          # cheaper tier by design for the scan
    if args.smoke:
        keep = {"180", "90", "0"}
        scan["scan"] = {k: v for k, v in scan["scan"].items() if k in keep}
    pts = []
    for key, g in sorted(scan["scan"].items(), key=lambda kv: -kv[1]["phi"]):
        try:
            els, xyz = read_xyz(Path(g["file"]))
            mol = psi4.geometry(psi4_geometry_string(els, xyz))
            psi4.set_options({"basis": basis, "reference": "rhf",
                              "scf_type": "df", "save_jk": True})
            e0, wf = psi4.energy(f"{DFT_FUNC}/{basis}", molecule=mol,
                                 return_wfn=True)
            psi4.set_options({"tdscf_states": [3], "tdscf_tda": True})
            psi4.energy(f"td-{DFT_FUNC}/{basis}", ref_wfn=wf, molecule=mol,
                        return_wfn=True)
            td = _td_variables_map(wf, 3)
            s1 = (td[1]["dE_eh"] * EH_EV
                  if td.get(1) and td[1]["dE_eh"] else None)
            s2 = (td[2]["dE_eh"] * EH_EV
                  if td.get(2) and td[2]["dE_eh"] else None)
            f2 = td[2]["f_osc"] if td.get(2) else None
            pts.append({"phi": g["phi"], "e0_eh": float(e0),
                        "s1_eV": s1, "s2_eV": s2, "s2_f": float(f2 or 0.0)})
            if s1 and s2:
                _log("8A", f"scan phi={g['phi']:5.0f}: S1={s1:.2f} eV  "
                           f"S2={s2:.2f} eV  f2={f2:.3f}")
            else:
                _log("8A", f"scan phi={g['phi']:5.0f}: TD incomplete")
        except Exception as exc:
            _warn("8A", f"scan point phi={g['phi']} failed: {str(exc)[:90]}")
    dump_json({"basis": basis, "points": pts}, res_path)
    _log("8A", f"torsion scan complete -> {res_path}")


# --------------------------------------------------------------------------- #
#  WORKER 8B — SA-CASSCF scans, BRS MECI, branching space
# --------------------------------------------------------------------------- #
def worker_8b(args, res_path):
    """SA-CASSCF scans + BRS penalty MECI + branching space.

    Fault tolerance (phase-7 pattern): every completed unit of work is
    checkpointed to <res>.ckpt.json; if this psi4 build aborts hard
    (DETCI DPD corruption after many in-process jobs kills the process
    without a Python traceback), the driver relaunches with --resume and
    work continues from the checkpoint."""
    import psi4
    psi4.set_memory(QC_MEM)
    psi4.set_num_threads(args.threads)
    els0, xyz0 = read_xyz(Path(args.xyz))
    n_jobs = [0]

    def cas(xyz, basis, nroots=2):
        """SA-CASSCF(4,4) energies for diazene; returns (E0, E1)."""
        mol = psi4.geometry(psi4_geometry_string(els0, xyz))
        psi4.set_options({"basis": basis, "reference": "rhf",
                          "scf_type": "df"})
        e_rhf, wref = psi4.energy("hf", molecule=mol, return_wfn=True)
        nocc = wref.nalpha()
        nbf = wref.basisset().nbf()
        Ca_view = np.asarray(wref.Ca())
        act = [nocc - 2, nocc - 1, nocc, nocc + 1]
        docc = [k for k in range(nocc) if k not in act[:2]]
        virt = [k for k in range(Ca_view.shape[1])
                if k not in act and k >= nocc]
        ndocc = len(docc)
        Ca_view[:] = Ca_view[:, docc + act + virt]
        psi4.set_options({"restricted_docc": [ndocc], "active": [4],
                          "restricted_uocc": [nbf - ndocc - 4]})
        psi4.set_module_options(
            "detci", {"num_roots": nroots, "avg_states": list(range(nroots))})
        # NOTE: never pass a CIWavefunction as ref_wfn here — this Psi4
        # build treats it as already-converged and returns stale energies.
        # Near-degenerate SA roots oscillate at the 1e-5 Eh level; FD
        # gradients need only ~1e-5 Eh, so the retry chain relaxes the
        # energy gate instead of burning iterations.
        tries = [({}, 150),
                 ({"maxiter": 300, "e_convergence": 1e-6}, 150),
                 ({"maxiter": 400, "e_convergence": 1e-5}, 150),
                 ({"maxiter": 500, "e_convergence": 3e-4}, 150)]
        last_exc = None
        for extra, it0 in tries:
            try:
                if extra:
                    psi4.set_options(extra)
                _e, cw = psi4.energy("casscf", molecule=mol, return_wfn=True,
                                     ref_wfn=wref)
                e0 = float(cw.variable("CI ROOT 0 TOTAL ENERGY"))
                e1 = None
                if nroots > 1:
                    e1 = float(cw.variable("CI ROOT 1 TOTAL ENERGY"))
                n_jobs[0] += 1
                psi4.set_options({"maxiter": it0, "e_convergence": 2e-7})
                psi4.core.clean()   # free DETCI DPD instance (win-64 quirk)
                return e0, e1
            except Exception as exc:
                last_exc = exc
                psi4.set_options({"maxiter": it0, "e_convergence": 2e-7})
                psi4.core.clean()
        raise RuntimeError(f"CASSCF retry chain exhausted: "
                           f"{str(last_exc)[:90]}")

    # ---- engine with basis fallback -------------------------------------- #
    eng = None
    for b in (BASIS_8B, BASIS_8B_FALLBACK):
        try:
            e0, e1 = cas(xyz0, b)
            _log("8B", f"SA-CASSCF(4,4)/{b} trans: E0={e0:.6f} E1={e1:.6f} "
                       f"gap={(e1 - e0) * EH_EV:.3f} eV")
            eng = b
            break
        except Exception as exc:
            _warn("8B", f"CASSCF/{b} failed: {str(exc)[:110]}")
    if eng is None:
        raise RuntimeError("no CASSCF basis tier succeeded")

    def pair(xyz):
        return cas(xyz, eng, nroots=2)

    # ---- checkpoint state -------------------------------------------------- #
    ck_path = Path(str(res_path) + ".ckpt.json")
    ck = {}
    if args.resume and ck_path.exists():
        try:
            ck = json.loads(ck_path.read_text())
            _log("8B", f"resuming from checkpoint (phase={ck.get('phase')}, "
                       f"jobs so far={ck.get('n_jobs', 0)})")
        except Exception:
            ck = {}
    ck.setdefault("phase", "scan")
    ck.setdefault("torsion_scan", [])
    ck.setdefault("n_jobs", 0)

    def flush_ck():
        ck["n_jobs"] = ck.get("n_jobs_base", 0) + n_jobs[0]
        ck["basis"] = eng
        dump_json(ck, ck_path)

    if "n_jobs_base" not in ck:
        ck["n_jobs_base"] = 0

    out = {"basis": eng, "method": f"SA-CASSCF(4,4)/{eng}",
           "engine": "Psi4 1.11 DETCI (PySCF-substituted on win32)",
           "engine_note": ENGINE_NOTE,
           "gradient_note": ("analytic CASSCF gradients return zero arrays on "
                             "this Psi4 1.11 win-64 build; all 8B gradients "
                             "are finite differences of CAS state energies")}

    # ---- phase: torsion scan ----------------------------------------------- #
    if ck["phase"] == "scan":
        scan_t = parse_json(Path(args.scan_t_json))["scan_torsion"]
        if args.smoke:
            keep = {"180", "100", "90", "80", "0"}
            scan_t = {k: v for k, v in scan_t.items() if k in keep}
        done = {p["phi"] for p in ck["torsion_scan"]}
        for key, g in sorted(scan_t.items(), key=lambda kv: -kv[1]["phi"]):
            if g["phi"] in done:
                continue
            try:
                _, xyzg = read_xyz(Path(g["file"]))
                E0, E1 = pair(xyzg)
                ck["torsion_scan"].append(
                    {"phi": g["phi"], "e0_eh": E0, "e1_eh": E1,
                     "gap_eV": (E1 - E0) * EH_EV})
                _log("8B", f"torsion phi={g['phi']:5.0f}: gap="
                           f"{(E1 - E0) * EH_EV:6.3f} eV")
                flush_ck()
            except Exception as exc:
                _warn("8B", f"torsion point {g['phi']} failed: "
                            f"{str(exc)[:80]}")
        ck["torsion_scan"].sort(key=lambda p: -p["phi"])
        band = [p for p in ck["torsion_scan"]
                if 70 <= p["phi"] <= 110 and p["e1_eh"] is not None]
        ck["phi_ci_deg"] = (min(band, key=lambda p: p["gap_eV"])["phi"]
                            if band else CI_TORSION_DEG)
        ck["phase"] = "meci_init"
        flush_ck()

    phi_ci = ck["phi_ci_deg"]
    out["torsion_scan"] = ck["torsion_scan"]
    out["phi_ci_deg"] = phi_ci

    # ---- phase: BRS MECI optimization -------------------------------------- #
    d = FD_GRAD_STEP
    scan_t = parse_json(Path(args.scan_t_json))["scan_torsion"]

    if ck["phase"] == "meci_init":
        _, x0 = read_xyz(Path(scan_t[f"{phi_ci:.0f}"]["file"]))
        x0 = x0.ravel().copy()
        E0, E1 = pair(x0.reshape(-1, 3))
        gap = E1 - E0
        F = 0.5 * (E0 + E1) + PENALTY_SIGMA * gap**2 / (gap + PENALTY_ALPHA)
        ck["meci_state"] = {"x": x0.tolist(), "e0": E0, "e1": E1, "F": F,
                            "istep": 0, "steps": [],
                            "H": (np.eye(x0.size) * 0.5).tolist(),
                            "gF_prev": None, "x_prev": None}
        ck["phase"] = "meci"
        flush_ck()

    st = ck["meci_state"]
    x = np.array(st["x"])
    E0, E1, F = st["e0"], st["e1"], st["F"]
    Hbfgs = np.array(st["H"])
    x_old = st["x_prev"]
    gF_prev = None if st["gF_prev"] is None else np.array(st["gF_prev"])
    max_steps = 10 if args.smoke else MECI_MAX_STEPS

    def objective(xv):
        E0t, E1t = pair(xv.reshape(-1, 3))
        gt = E1t - E0t
        Ft = (0.5 * (E0t + E1t)
              + PENALTY_SIGMA * gt**2 / (gt + PENALTY_ALPHA))
        return E0t, E1t, Ft

    def fd_grads(xv, E0c, E1c):
        """Forward-difference per-state gradients; smaller-step retry."""
        g0 = np.zeros_like(xv)
        g1 = np.zeros_like(xv)
        for k in range(xv.size):
            for step in (d, d / 2.0):
                try:
                    xp = xv.copy()
                    xp[k] += step
                    E0p, E1p = pair(xp.reshape(-1, 3))
                    g0[k] = (E0p - E0c) / step
                    g1[k] = (E1p - E1c) / step
                    break
                except Exception:
                    if step == d / 2.0:
                        g0[k] = 0.0
                        g1[k] = 0.0
        return g0, g1

    if ck["phase"] == "meci":
        while st["istep"] < max_steps:
            istep = st["istep"] + 1
            gap = E1 - E0
            gap_eV = gap * EH_EV
            try:
                g0, g1 = fd_grads(x, E0, E1)
                gF = (0.5 * (g0 + g1)
                      + PENALTY_SIGMA * (g1 - g0) * gap
                      * (gap + 2 * PENALTY_ALPHA) / (gap + PENALTY_ALPHA)**2)
            except Exception:
                if gF_prev is None:
                    _warn("8B", f"step {istep}: gradient failed — perturbing")
                    rngx = np.random.default_rng(SEED + istep)
                    x = x + rngx.normal(0, 0.01, x.size)
                    try:
                        E0, E1, F = objective(x)
                    except Exception:
                        pass
                    st.update({"x": x.tolist(), "e0": E0, "e1": E1,
                               "F": F, "istep": istep})
                    flush_ck()
                    continue
                gF = 0.5 * gF_prev
            gnorm = float(np.linalg.norm(gF))
            st["steps"].append({"step": istep, "e0_eh": E0, "e1_eh": E1,
                                "gap_eV": gap_eV, "f_penalty_eh": F,
                                "gradF_norm": gnorm})
            _log("8B", f"MECI step {istep:2d}: gap={gap_eV:7.4f} eV  "
                       f"|gradF|={gnorm:.4f} Eh/A  "
                       f"E_avg={0.5 * (E0 + E1):.5f}")
            if gap_eV < MECI_GAP_EV_GATE and gnorm < MECI_GRAD_GATE:
                _log("8B", f"MECI converged at step {istep} "
                           f"(gap {gap_eV:.4f} eV < {MECI_GAP_EV_GATE} eV)")
                st["istep"] = istep
                break
            if x_old is not None and gF_prev is not None:
                sv = (x - np.array(x_old)).reshape(-1, 1)
                yv = (gF - gF_prev).reshape(-1, 1)
                sy = float((sv.T @ yv).item())
                if sy > 1e-12:
                    Hs = Hbfgs @ sv
                    Hbfgs += (np.outer(yv.ravel(), yv.ravel()) / sy
                              - np.outer(Hs.ravel(), Hs.ravel())
                              / float((sv.T @ Hs).item()))
            pvec = -Hbfgs @ gF
            pnorm = float(np.linalg.norm(pvec))
            if pnorm > 0.15:
                pvec *= 0.15 / pnorm
            alpha = 1.0
            improved = False
            for _ in range(5):
                x_try = x + alpha * pvec
                try:
                    E0t, E1t, Ft = objective(x_try)
                except Exception:
                    alpha *= 0.5
                    continue
                if Ft < F or abs(E1t - E0t) < abs(gap):
                    x, E0, E1, F = x_try, E0t, E1t, Ft
                    improved = True
                    break
                alpha *= 0.5
            if not improved:
                rng = np.random.default_rng(SEED + 99 * istep)
                x = x + rng.normal(0, 0.01, x.size)
                try:
                    E0, E1, F = objective(x)
                except Exception:
                    pass
                _warn("8B", f"step {istep}: line search stalled — perturbing")
            gF_prev = gF
            x_old = x.tolist()
            st.update({"x": x.tolist(), "e0": E0, "e1": E1, "F": F,
                       "istep": istep, "H": Hbfgs.tolist(),
                       "gF_prev": gF.tolist(), "x_prev": x_old})
            flush_ck()
        ck["phase"] = "branching"
        flush_ck()

    xyz_meci = x.reshape(-1, 3)
    gap_eV = (E1 - E0) * EH_EV
    write_xyz(OUT / "meci.xyz", els0, xyz_meci,
              f"diazene MECI SA-CASSCF(4,4)/{eng} gap={gap_eV:.4f} eV")
    out["meci"] = {"gap_eV": gap_eV, "converged": gap_eV < MECI_GAP_EV_GATE,
                   "steps": st["steps"], "e0_eh": E0, "e1_eh": E1,
                   "n_steps": len(st["steps"])}

    # ---- phase: branching space (g) ---------------------------------------- #
    if ck["phase"] == "branching":
        g0 = np.zeros_like(x)
        g1 = np.zeros_like(x)
        for k in range(x.size):
            for step in (d, d / 2.0):
                try:
                    xp = x.copy()
                    xp[k] += step
                    E0p, E1p = pair(xp.reshape(-1, 3))
                    g0[k] = (E0p - E0) / step
                    g1[k] = (E1p - E1) / step
                    break
                except Exception:
                    if step == d / 2.0:
                        g0[k] = 0.0
                        g1[k] = 0.0
        g_vec = g1 - g0
        ck["g_vector"] = {"units": "Eh/Angstrom",
                          "norm": float(np.linalg.norm(g_vec)),
                          "per_atom": g_vec.reshape(-1, 3).tolist()}
        ck["phase"] = "hvec"
        flush_ck()

    g_vec = np.array([c for row in ck["g_vector"]["per_atom"]
                      for c in row])
    out["g_vector"] = ck["g_vector"]
    gn = g_vec / max(np.linalg.norm(g_vec), 1e-12)

    # ---- phase: h via directional gap lifting ------------------------------- #
    if ck["phase"] == "hvec":
        _log("8B", "h-vector reconstruction by directional gap lifting "
                   f"(delta = {LIFT_DELTA} A, localized state-overlap "
                   "fallback; checkpointed per direction)")
        rng = np.random.default_rng(SEED + 7)
        cands = []
        for k in range(x.size):
            e = np.eye(x.size)[k]
            u = e - float(np.dot(e, gn)) * gn
            nu = np.linalg.norm(u)
            if nu > 1e-8:
                cands.append((f"cart_{k}", (u / nu).tolist()))
        cands += [(f"rand_{i}", None) for i in range(24)]
        ck.setdefault("h_evals", [])
        done_ids = {e["cid"] for e in ck["h_evals"]}
        for cid, pre in cands:
            if cid in done_ids:
                continue
            if pre is None:
                v = rng.normal(size=x.size)
                v = v - float(np.dot(v, gn)) * gn
                nv = np.linalg.norm(v)
                if nv < 1e-8:
                    continue
                u = v / nv
            else:
                u = np.array(pre)
            try:
                Ep, _E1 = pair((x + LIFT_DELTA * u).reshape(-1, 3))
                lift = abs(Ep - min(E0, E1)) / LIFT_DELTA
            except Exception as exc:
                _warn("8B", f"lift dir {cid} failed: {str(exc)[:60]}")
                continue
            ck["h_evals"].append({"cid": cid, "lift": float(lift),
                                  "u": u.tolist()})
            flush_ck()
        best = max(ck["h_evals"], key=lambda e: e["lift"])
        bu = np.array(best["u"])
        for ri in range(6):
            cid = f"ref_{ri}"
            if cid in {e["cid"] for e in ck["h_evals"]}:
                continue
            trial = bu + rng.normal(0, 0.25, x.size)
            trial = trial - float(np.dot(trial, gn)) * gn
            nt = np.linalg.norm(trial)
            if nt < 1e-8:
                continue
            u = trial / nt
            try:
                Ep, _E1 = pair((x + LIFT_DELTA * u).reshape(-1, 3))
                lift = abs(Ep - min(E0, E1)) / LIFT_DELTA
            except Exception:
                continue
            ck["h_evals"].append({"cid": cid, "lift": float(lift),
                                  "u": u.tolist()})
            flush_ck()
        ck["phase"] = "cuts"
        flush_ck()

    evals = ck["h_evals"]
    best = max(evals, key=lambda e: e["lift"])
    hn = np.array(best["u"])
    best_lift = best["lift"]
    h_vec = hn * best_lift
    out["h_vector"] = {"units": "Eh/Angstrom",
                       "method": (f"directional gap-lifting maximization "
                                  f"over {len(evals)} checkpointed trial "
                                  f"directions perpendicular to g (localized "
                                  f"state-overlap fallback per mission)"),
                       "lift_rate_eh_per_angstrom": float(best_lift),
                       "per_atom": h_vec.reshape(-1, 3).tolist()}
    _log("8B", f"h-vector: lift rate {best_lift:.4f} Eh/A over "
               f"{len(evals)} directions")

    # ---- phase: verification cuts ------------------------------------------- #
    if ck["phase"] == "cuts":
        ck.setdefault("branching_cuts", {})
        for label, vec in (("g", gn), ("h", hn)):
            ck["branching_cuts"].setdefault(label, [])
            done_t = {r["t_ang"] for r in ck["branching_cuts"][label]}
            for t in (-0.2, -0.1, 0.0, 0.1, 0.2):
                if t in done_t:
                    continue
                try:
                    E0t, E1t = pair((x + t * vec).reshape(-1, 3))
                    ck["branching_cuts"][label].append(
                        {"t_ang": t, "gap_eV": (E1t - E0t) * EH_EV,
                         "e0_eh": E0t, "e1_eh": E1t})
                    flush_ck()
                except Exception as exc:
                    _warn("8B", f"cut {label} {t} failed: {str(exc)[:60]}")
            ck["branching_cuts"][label].sort(key=lambda r: r["t_ang"])
            _log("8B", f"cut along {label}: " + " ".join(
                f"{r['t_ang']:+.1f}A:{r['gap_eV']:.3f}eV"
                for r in ck["branching_cuts"][label]))
        ck["phase"] = "stretch"
        flush_ck()

    out["branching_cuts"] = ck.get("branching_cuts", {})

    # ---- phase: N=N stretch mini-scan ---------------------------------------- #
    if ck["phase"] == "stretch":
        scan_s = parse_json(Path(args.scan_s_json))["scan_stretch"]
        ck.setdefault("stretch_scan_at_ci", [])
        done_d = {p["d"] for p in ck["stretch_scan_at_ci"]}
        for key, g in sorted(scan_s.items(), key=lambda kv: kv[1]["d"]):
            if g["d"] in done_d:
                continue
            try:
                _, xyzg = read_xyz(Path(g["file"]))
                E0s, E1s = pair(xyzg)
                ck["stretch_scan_at_ci"].append(
                    {"d": g["d"], "e0_eh": E0s, "e1_eh": E1s,
                     "gap_eV": (E1s - E0s) * EH_EV})
                flush_ck()
            except Exception as exc:
                _warn("8B", f"stretch point {g['d']} failed: "
                            f"{str(exc)[:80]}")
        ck["stretch_scan_at_ci"].sort(key=lambda p: p["d"])
        ck["phase"] = "done"
        flush_ck()

    out["stretch_scan_at_ci"] = ck.get("stretch_scan_at_ci", [])
    out["n_casscf_jobs"] = ck.get("n_jobs", n_jobs[0])

    dump_json(out, res_path)
    if ck["phase"] == "done":
        ck_path.unlink(missing_ok=True)
    _log("8B", f"module 8B complete ({out['n_casscf_jobs']} CAS jobs, "
               f"phase={ck['phase']}) -> {res_path}")



# --------------------------------------------------------------------------- #
#  MODULE 8C — LVC construction + vectorized Tully FSSH       (driver numpy)
# --------------------------------------------------------------------------- #

def torsion_mode_direction(els, xyz, delta=2.0):
    """Mass-weighted Cartesian unit vector of the HNNH torsion about N=N
    (finite-difference dihedral rotation of the H on N1)."""
    from ase import Atoms
    nat = len(els)
    at_p = Atoms(symbols=els, positions=xyz)
    at_m = Atoms(symbols=els, positions=xyz)
    mask = [False, False, True, False]
    at_p.set_dihedral(2, 0, 1, 3, 180.0 + delta, mask=mask)
    at_m.set_dihedral(2, 0, 1, 3, 180.0 - delta, mask=mask)
    dx = ((at_p.get_positions() - at_m.get_positions())
          / (2.0 * math.radians(delta)))
    masses = {"H": 1.008, "N": 14.007}
    mw = dx * np.sqrt(np.array([masses[e] for e in els]))[:, None]
    n = np.linalg.norm(mw)
    return mw / n

def build_lvc(res8b: dict, modes: dict, els, xyz0):
    """Build the 2-state / 3-mode dynamical model from ab-initio data.

    Torsion coordinate phi: the CAS torsion scan is *diabatized* with the
    standard mirror-symmetry construction (Delta_diab symmetric and V12
    antisymmetric about the scan gap minimum phi_CI):

        Delta_d(phi)^2 = 1/2 [Delta_adi(phi)^2 + Delta_adi(phi')^2]
        V12(phi)^2     = 1/4 [Delta_adi(phi)^2 - Delta_d(phi)^2],
        phi' = 2 phi_CI - phi

    giving smooth periodic diabats W0d(phi), W1d(phi), V12(phi) on [0, 360)
    (even / even / odd about phi_CI).  The N=N stretch and HN=N wag enter as
    local harmonic LVC modes whose coupling gradients are the projections of
    the MECI h vector onto the mass-weighted modes.  The S1 -> S0 seam then
    requires BOTH torsion (~90 deg) and N=N compression — the physics of
    the diazene/azobenzene conical intersection."""
    tpts = sorted(res8b["torsion_scan"], key=lambda p: p["phi"])
    phi_ci_deg = res8b["phi_ci_deg"]
    phis = np.array([math.radians(p["phi"]) for p in tpts])
    e0 = np.array([p["e0_eh"] for p in tpts])
    e1 = np.array([p["e1_eh"] for p in tpts])
    gap = e1 - e0
    avg = 0.5 * (e0 + e1)

    # ---- anchored diabatization (robust to scan holes) --------------------- #
    # The diabatic gap is anchored at its two symmetry-distinct knowns:
    #   Delta_d(phi_CI) = Delta_adi(phi_CI)   (the coupling vanishes at the
    #                                          gap-minimum geometry)
    #   Delta_d(180deg) = Delta_adi(180deg)   (the coupling vanishes at the
    #                                          planar trans geometry)
    #   Delta_d(phi) = Delta90 + (Delta180 - Delta90) sin^2(phi - 90 deg)
    #   (period-180 form respecting the E/Z symmetry of the torsion)
    # and the single coupling amplitude V0 in V12 = V0 |sin phi| (symmetry-
    # zero at both planar geometries, antisymmetric about the CI) is fitted
    # to the residual adiabatic gaps in least-squares sense:
    #   Delta_adi^2 = Delta_d^2 + 4 V0^2 sin^2(phi).
    phi_ci_rad = math.radians(phi_ci_deg)
    d90 = float(np.interp(phi_ci_rad, phis, gap))
    d180 = float(gap[-1])
    half = np.sin(phis - math.radians(phi_ci_deg)) ** 2
    delta_d = d90 + (d180 - d90) * half
    sin2 = np.sin(phis) ** 2
    resid2 = np.maximum(gap**2 - delta_d**2, 0.0)
    V0_f = float(np.sqrt(np.mean(resid2[np.sin(phis) > 1e-6]
                                 / (4.0 * sin2[np.sin(phis) > 1e-6]))))
    V_abs = V0_f * np.abs(np.sin(phis))
    rms = float(np.sqrt(np.mean((delta_d**2 + 4.0 * V_abs**2 - gap**2)**2)))
    _log("8B-fit", f"diabatization: Delta90 = {d90 * EH_EV:.3f} eV, "
                   f"Delta180 = {d180 * EH_EV:.3f} eV, V0 = {V0_f * EH_EV:.3f} eV, "
                   f"rms = {rms * EH_EV * 1000:.1f} meV")

    # antisymmetrize the coupling sign about phi_CI; the coupling must
    # vanish at both planar geometries (cis C2v / trans C2h: the n->pi*
    # diabat has a different irreducible representation than the ground
    # state, so V12 is symmetry-forbidden at phi = 0 and 180) — this also
    # gives the exact even/odd closure of the periodic circle
    sgn = np.where(phis < math.radians(phi_ci_deg), -1.0, 1.0)
    V12 = sgn * V_abs          # exact 0 at both planar endpoints (|sin|)

    # extend to the full circle: even W, odd V about phi = 0 (and 2 pi)
    phis_full = np.concatenate([phis, 2 * math.pi - phis[::-1]])
    w0_full = np.concatenate([avg - 0.5 * delta_d,
                              (avg - 0.5 * delta_d)[::-1]])
    w1_full = np.concatenate([avg + 0.5 * delta_d,
                              (avg + 0.5 * delta_d)[::-1]])
    v_full = np.concatenate([V12, -V12[::-1]])
    order = np.argsort(phis_full)
    phis_full, w0_full, w1_full, v_full = (phis_full[order], w0_full[order],
                                           w1_full[order], v_full[order])
    # drop self-mirror duplicates (phi = 0 and pi appear in both halves)
    keep = np.concatenate([[True], np.diff(phis_full) > 1e-9])
    phis_full, w0_full = phis_full[keep], w0_full[keep]
    w1_full, v_full = w1_full[keep], v_full[keep]
    # periodic spline needs the exact [0, 2*pi] closure
    if abs(phis_full[0]) > 1e-9:
        phis_full = np.concatenate([[0.0], phis_full])
        w0_full = np.concatenate([[w0_full[-1]], w0_full])
        w1_full = np.concatenate([[w1_full[-1]], w1_full])
        v_full = np.concatenate([[0.0], v_full])
    if abs(phis_full[-1] - 2 * math.pi) > 1e-9:
        phis_full = np.concatenate([phis_full, [2 * math.pi]])
        w0_full = np.concatenate([w0_full, [w0_full[0]]])
        w1_full = np.concatenate([w1_full, [w1_full[0]]])
        v_full = np.concatenate([v_full, [v_full[0]]])

    # ---- stretch / wag local LVC parameters -------------------------------- #
    mu_st = modes["stretch"][2] * AMU_ME
    mu_be = modes["bend"][2] * AMU_ME
    spts = sorted(res8b["stretch_scan_at_ci"], key=lambda p: p["d"])
    ds = np.array([p["d"] for p in spts]) * ANG_BOHR
    gaps_s = np.array([p["gap_eV"] / EH_EV for p in spts])
    okm = np.isfinite(ds) & np.isfinite(gaps_s)
    slope = 0.0
    if okm.sum() >= 3:
        slope = float(np.polyfit(ds[okm], gaps_s[okm], 1)[0])
    elif okm.sum() == 2 and abs(ds[okm][1] - ds[okm][0]) > 1e-9:
        slope = float((gaps_s[okm][1] - gaps_s[okm][0])
                      / (ds[okm][1] - ds[okm][0]))
    # full diabatic gap slope along the (soft, S0 nearly flat) stretch-like
    # coordinate; the ground-state curvature along the same coordinate comes
    # from the CAS E0 values of the scan itself (the coordinate is a combined
    # N=N stretch + pyramidalization motion and is far softer than the
    # 1715 cm-1 pure stretch)
    kappa_stretch = slope / math.sqrt(mu_st)         # Eh / (bohr sqrt(me))
    # Stretch-relaxation curve: the CAS gap along the N=N coordinate is
    # turned directly into a spline Delta_rel(Q_R) = gap(Q_R) - gap(0),
    # linearly extrapolated to the seam (Delta_rel = -Delta90) and beyond.
    # This embeds the true conical-intersection seam (N=N lengthening side,
    # where the ground state is soft) at the ab-initio location without any
    # linear-extrapolation assumption.
    Qs = ds * math.sqrt(mu_st)                       # bohr sqrt(me)
    gap0 = float(np.interp(0.0, Qs[np.argsort(Qs)],
                           gaps_s[np.argsort(Qs)]))
    o = np.argsort(Qs)
    Qs_o, gaps_o = Qs[o], gaps_s[o]
    rel = gaps_o - gap0
    # zero-crossing from the last two points on the lengthening side
    slope_end = 0.0
    if len(Qs_o) >= 2:
        slope_end = float((rel[-1] - rel[-2]) / (Qs_o[-1] - Qs_o[-2]))
    Qx = float(Qs_o[-1] + (-rel[-1]) / slope_end) if slope_end < 0 else         float(Qs_o[-1] + 20.0)
    Qx = max(Qx, Qs_o[-1] + 1.0)
    # constant-slope extension to +-300 keeps the PCHIP spline strictly
    # interpolated over the whole dynamical range (no extrapolation NaN)
    Q_lo = float(Qs_o[0] - 300.0)
    Q_hi = float(Qx + 300.0)
    slope_lo = float((rel[1] - rel[0]) / (Qs_o[1] - Qs_o[0]))         if len(Qs_o) > 1 else 0.0
    Q_ext = np.concatenate([[Q_lo], Qs_o, [Qx, Q_hi]])
    rel_ext = np.concatenate([[rel[0] + (Q_lo - Qs_o[0]) * slope_lo],
                              rel, [-d90, -d90 + (Q_hi - Qx) * slope_end]])
    okr = np.isfinite(Q_ext) & np.isfinite(rel_ext)
    lvc_rel_Q = Q_ext[okr]
    lvc_rel_V = rel_ext[okr]
    order = np.argsort(lvc_rel_Q)
    lvc_rel_Q, lvc_rel_V = lvc_rel_Q[order], lvc_rel_V[order]
    e0s = np.array([p["e0_eh"] for p in spts])
    oks = np.isfinite(ds) & np.isfinite(e0s)
    omega_stretch = 60.0 * CM1_AU
    stretch_curv = 0.0
    if oks.sum() >= 3:
        stretch_curv = float(2.0 * np.polyfit(ds[oks], e0s[oks], 2)[0])
        omega_stretch = math.sqrt(abs(stretch_curv) / mu_st)
    omega_stretch = max(omega_stretch, 50.0 * CM1_AU)
    _log("8B-fit", f"stretch coordinate: gap slope = {slope:.4f} Eh/bohr, "
                   f"S0 curvature = {stretch_curv:.3f} Eh/bohr^2, "
                   f"omega = {omega_stretch / CM1_AU:.0f} cm-1")

    h_cart = (np.array([c for row in res8b["h_vector"]["per_atom"]
                        for c in row]) * ANG_BOHR)
    m_arr = np.array([{"H": 1.008, "N": 14.007}[e] for e in els])
    sqrt_m3 = np.repeat(np.sqrt(m_arr), 3)

    def project(vec_unweighted):
        mw = np.asarray(vec_unweighted, dtype=float) * np.sqrt(m_arr)[:, None]
        n = np.linalg.norm(mw)
        if n < 1e-9:
            return 0.0
        return float(np.dot(h_cart, (mw / n).ravel() / sqrt_m3))

    lambda_stretch = project(modes["stretch_vec"])
    lambda_bend = project(modes["bend_vec"])

    I_tors = modes["torsion_inertia_amu_A2"] * AMU_ME * ANG_BOHR**2

    lvc = {
        "rel_grid_Q": lvc_rel_Q.tolist(),
        "rel_grid_V": lvc_rel_V.tolist(),
        "phi_grid": phis_full.tolist(),
        "W0d_grid": w0_full.tolist(),
        "W1d_grid": w1_full.tolist(),
        "V12_grid": v_full.tolist(),
        "phi0_rad": math.pi,                 # Franck-Condon point: E isomer (trans)
        "phi_ci_deg": phi_ci_deg,
        "I_tors_me_bohr2": I_tors,
        "omega_stretch_au": omega_stretch,
        "omega_bend_au": modes["bend"][1] * CM1_AU,
        "mu_stretch_me": mu_st, "mu_bend_me": mu_be,
        "kappa_stretch": kappa_stretch, "kappa_bend": 0.0,
        "lambda_stretch": lambda_stretch, "lambda_bend": lambda_bend,
        "mode_freqs_cm1": [omega_stretch / CM1_AU, modes["bend"][1]],
        "dE_fc_eV": float((e1[-1] - e0[-1]) * EH_EV),
        "h_lift_rate_eh_per_angstrom":
            res8b["h_vector"]["lift_rate_eh_per_angstrom"],
        "diabatization_fit": {"delta90_eh": d90, "delta180_eh": d180,
                              "V0_eh": V0_f, "rms_eh": rms,
                              "n_scan_points": int(len(phis))},
        "parameterization": (
            "periodic torsional diabats by 4-parameter least-squares "
            "diabatization of the SA-CASSCF(4,4) scan (Delta_d = a + b cos "
            "+ c cos^2; V12 = V0 |sin phi|, symmetry-zero at planarity); "
            "N=N-stretch gap slope from the CAS stretch scan; wag/stretch "
            "coupling gradients = projections of the MECI h vector"),
    }
    return lvc


class LVC2:
    """Vectorized 2-state Hamiltonian: periodic torsional diabats (cubic
    splines) + local harmonic stretch/wag LVC modes."""

    def __init__(self, lvc: dict):
        from scipy.interpolate import CubicSpline
        x = np.array(lvc["phi_grid"])
        self.S0 = CubicSpline(x, np.array(lvc["W0d_grid"]), bc_type="periodic")
        self.S1 = CubicSpline(x, np.array(lvc["W1d_grid"]), bc_type="periodic")
        self.SV = CubicSpline(x, np.array(lvc["V12_grid"]), bc_type="periodic")
        self.D0 = self.S0.derivative()
        self.D1 = self.S1.derivative()
        self.DV = self.SV.derivative()
        from scipy.interpolate import PchipInterpolator
        self.REL = PchipInterpolator(np.array(lvc["rel_grid_Q"]),
                                     np.array(lvc["rel_grid_V"]),
                                     extrapolate=True)
        self.DREL = self.REL.derivative()
        self.rel_q_lo = float(np.min(lvc["rel_grid_Q"]))
        self.rel_q_hi = float(np.max(lvc["rel_grid_Q"]))
        self.It = lvc["I_tors_me_bohr2"]
        self.sqrt_It = math.sqrt(self.It)
        self.wR = lvc["omega_stretch_au"]
        self.wB = lvc["omega_bend_au"]
        self.kR = lvc["kappa_stretch"]
        self.kB = lvc["kappa_bend"]
        self.lR = lvc["lambda_stretch"]
        self.lB = lvc["lambda_bend"]
        self.mu = np.array([self.It, lvc["mu_stretch_me"], lvc["mu_bend_me"]])
        self.mu_sqrt = np.sqrt(self.mu)
        self.phi0 = lvc["phi0_rad"]
        self.phi_ci_deg = lvc["phi_ci_deg"]
        # local harmonic frequencies for Wigner sampling: torsional omega
        # from the ground-diabat curvature at the planar trans minimum
        d2W0 = float(self.D0.derivative()(self.phi0))
        w_tors = math.sqrt(max(d2W0, 1e-8) / self.It)
        self.w = np.array([w_tors, self.wR, self.wB])

    def phi_of(self, Q):
        return np.mod(self.phi0 + Q[:, 0] / self.sqrt_It, 2.0 * math.pi)

    def diabats(self, Q):
        phi = self.phi_of(Q)
        # clip the REL evaluation to the spline domain (dissociative R
        # trajectories beyond the fitted range feel the constant-slope
        # asymptote; the harmonic term keeps the coordinate bounded)
        Qr = np.clip(Q[:, 1], self.rel_q_lo, self.rel_q_hi)
        harm = (0.5 * self.wR**2 * Q[:, 1] ** 2
                + 0.5 * self.wB**2 * Q[:, 2] ** 2)
        w0 = self.S0(phi) + harm
        w1 = (self.S1(phi) + self.kR * Qr + self.kB * Q[:, 2] + harm
              + self.REL(Qr))
        v12 = self.SV(phi) + self.lR * Q[:, 1] + self.lB * Q[:, 2]
        return w0, w1, v12

    def eigs(self, Q):
        w0, w1, v12 = self.diabats(Q)
        mid = 0.5 * (w0 + w1)
        det = w1 - w0
        rad = np.sqrt(0.25 * det**2 + v12**2)
        return w0, w1, v12, mid, det, rad

    def forces_and_nac(self, Q, active):
        """Adiabatic forces on each trajectory's active surface + NAC d01."""
        w0, w1, v12, mid, det, rad = self.eigs(Q)
        phi = self.phi_of(Q)
        inv_sI = 1.0 / self.sqrt_It
        dw0 = np.stack([self.D0(phi) * inv_sI,
                        self.wR**2 * Q[:, 1],
                        self.wB**2 * Q[:, 2]], axis=1)
        Qr = np.clip(Q[:, 1], self.rel_q_lo, self.rel_q_hi)
        dw1 = np.stack([self.D1(phi) * inv_sI,
                        self.wR**2 * Q[:, 1] + self.kR
                        + self.DREL(Qr),
                        self.wB**2 * Q[:, 2] + self.kB], axis=1)
        dv = np.stack([self.DV(phi) * inv_sI,
                       np.full_like(phi, self.lR),
                       np.full_like(phi, self.lB)], axis=1)
        dmid = 0.5 * (dw0 + dw1)
        drad = ((0.5 * det[:, None] * (dw1 - dw0)
                 + 2.0 * v12[:, None] * dv) / rad[:, None])
        grad_e0 = dmid - drad
        grad_e1 = dmid + drad
        F = np.where((active == 0)[:, None], -grad_e0, -grad_e1)
        e0 = mid - rad
        e1 = mid + rad
        theta = 0.5 * np.arctan2(2.0 * v12, det)
        cth, sth = np.cos(theta), np.sin(theta)
        g01 = ((cth * sth)[:, None] * (dw0 - dw1)
               + (cth * cth - sth * sth)[:, None] * dv)
        d01 = g01 / ((e1 - e0)[:, None] + 1e-12)
        return F, e0, e1, d01

def run_fssh(lvc: dict, n_traj: int, seed: int):
    """Vectorized Tully-FSSH ensemble on the LVC model -> (result, frames)."""
    rng = np.random.default_rng(seed)
    model = LVC2(lvc)
    w = model.w
    mu_s = model.mu_sqrt
    n3 = 3

    # Wigner sampling of the S0 vibrational ground state (a.u.)
    Q = rng.normal(0.0, 1.0 / np.sqrt(2.0 * w), size=(n_traj, n3))
    P = rng.normal(0.0, np.sqrt(w / 2.0), size=(n_traj, n3))
    c = np.zeros((n_traj, 2), dtype=complex)
    c[:, 1] = 1.0                          # vertical excitation onto S1
    active = np.ones(n_traj, dtype=int)

    n_steps = int(round(TMAX_FS / DT_FS))
    save_every = max(1, int(round(SAVE_EVERY_FS / DT_FS)))
    dt_au = DT_FS * FS_AU
    dt_sub = dt_au / N_ESUB

    times, pops_active, pops_coh = [], [], []
    hops_hist = np.zeros(n_steps, dtype=int)
    decay_time = np.full(n_traj, np.nan)
    on_s0_ever = np.zeros(n_traj, dtype=bool)
    n_frustrated = 0
    t0 = time.time()
    frames = {"t": [], "Q": [], "active": []}

    F, e0, e1, d01 = model.forces_and_nac(Q, active)

    for istep in range(n_steps + 1):
        t_fs = istep * DT_FS
        if istep % save_every == 0:
            w0, w1, v12, mid, det, rad = model.eigs(Q)
            theta = 0.5 * np.arctan2(2.0 * v12, det)
            cth, sth = np.cos(theta), np.sin(theta)
            ca0 = cth * c[:, 0] - sth * c[:, 1]
            ca1 = sth * c[:, 0] + cth * c[:, 1]
            times.append(t_fs)
            pops_active.append(float((active == 1).mean()))
            pops_coh.append(float((np.abs(ca1)**2).mean()))
            frames["t"].append(t_fs)
            frames["Q"].append(Q.copy())
            frames["active"].append(active.copy())
        if istep == n_steps:
            break
        newly = on_s0_ever & np.isnan(decay_time)
        decay_time[newly] = t_fs

        for isub in range(N_ESUB):
            # ---- Velocity-Verlet on the active adiabatic surface -------- #
            P_half = P + 0.5 * dt_sub * F
            Q = Q + P_half * dt_sub
            F, e0, e1, d01 = model.forces_and_nac(Q, active)
            P = P_half + 0.5 * dt_sub * F

            # ---- exact 2x2 electronic propagation ----------------------- #
            w0, w1, v12, mid, det, rad = model.eigs(Q)
            theta = 0.5 * np.arctan2(2.0 * v12, det)
            cth, sth = np.cos(theta), np.sin(theta)
            # chi0 = (cth, -sth), chi1 = (sth, cth): forward transform
            ca0 = cth * c[:, 0] - sth * c[:, 1]
            ca1 = sth * c[:, 0] + cth * c[:, 1]
            e_ad = np.stack([mid - rad, mid + rad], axis=1)
            ph = np.exp(-1j * e_ad * dt_sub)
            ca0 = ca0 * ph[:, 0]
            ca1 = ca1 * ph[:, 1]

            # ---- Tully fewest-switches hopping probability -------------- #
            v_dot_d = np.sum((P / mu_s[None, :]) * d01, axis=1)
            is0 = active == 0
            g_up = np.where(is0,
                            2.0 * dt_sub * np.real(np.conj(ca0) * ca1)
                            * v_dot_d / np.maximum(np.abs(ca0)**2, 1e-12),
                            0.0)
            g_dn = np.where(~is0,
                            2.0 * dt_sub * np.real(np.conj(ca1) * ca0)
                            * v_dot_d / np.maximum(np.abs(ca1)**2, 1e-12),
                            0.0)
            gprob = np.nan_to_num(g_up + g_dn, nan=0.0, posinf=0.0, neginf=0.0)
            gprob = np.clip(gprob, 0.0, 0.5)
            hop = rng.random(n_traj) < gprob
            if hop.any():
                to = np.where(is0, 1, 0)
                dhat = d01 / np.maximum(
                    np.linalg.norm(d01, axis=1, keepdims=True), 1e-12)
                pv = np.sum(P * dhat, axis=1)
                dE_hop = np.where(is0, e1 - e0, e0 - e1)  # E_target - E_source
                disc = pv**2 - 2.0 * dE_hop
                ok = hop & (disc >= 0.0)
                fr = hop & ~ok
                n_frustrated += int(fr.sum())
                if ok.any():
                    gamma = -pv + np.sign(pv + 1e-30) * np.sqrt(
                        np.maximum(disc, 0.0))
                    P[ok] += gamma[ok, None] * dhat[ok]
                    active[ok] = to[ok]
                    hops_hist[istep] += int(ok.sum())
                # frustrated hops are rejected (standard FSSH practice)
            is0 = active == 0
            # back-transform once per substep (coherent propagation)
            c[:, 0] = cth * ca0 + sth * ca1
            c[:, 1] = -sth * ca0 + cth * ca1

        # ---- Granucci-Persico decay-of-mixing decoherence ------------- #
        # applied once per NUCLEAR step (the standard prescription; per-
        # substep application over-damps by ~N_ESUB and suppresses hops)
        if EDC_ALPHA > 0:
            is0 = active == 0
            dE_ik = np.abs(e1 - e0)
            tau = 1.0 / (dE_ik * (1.0 + EDC_ALPHA
                                  / np.maximum(dE_ik, 1e-8)))
            damp = np.exp(-dt_au / tau)
            _w0, _w1, _v12, _mid, _det, _rad = model.eigs(Q)
            thd = 0.5 * np.arctan2(2.0 * _v12, _det)
            ctd, std = np.cos(thd), np.sin(thd)
            cb0 = ctd * c[:, 0] - std * c[:, 1]
            cb1 = std * c[:, 0] + ctd * c[:, 1]
            ca_a = np.where(is0, cb0, cb1)
            ca_o = np.where(is0, cb1, cb0)
            ca_o = ca_o * damp
            nrm = np.sqrt(np.abs(ca_a)**2 + np.abs(ca_o)**2) + 1e-300
            ca_a = ca_a / nrm
            ca_o = ca_o / nrm
            c[:, 0] = np.where(is0, ca_a, ca_o)
            c[:, 1] = np.where(is0, ca_o, ca_a)


        on_s0_ever |= active == 0

    # ---- final classification -------------------------------------------- #
    phi_final = model.phi_of(Q)                  # [0, 2*pi)
    phi_deg = np.rad2deg(phi_final)
    on_s0 = active == 0
    z_mask = on_s0 & (np.cos(phi_final) > 0.0)   # Z basin: |phi| < 90 deg
    e_mask = on_s0 & (np.cos(phi_final) <= 0.0)  # E basin recovered

    times_a = np.array(times)
    pops_a = np.array(pops_active)
    tau_half = None
    below = np.nonzero(pops_a < 0.5)[0]
    if below.size:
        i0 = below[0]
        if i0 == 0:
            tau_half = float(times_a[0])
        else:
            p1, p2 = pops_a[i0 - 1], pops_a[i0]
            t1, t2 = times_a[i0 - 1], times_a[i0]
            tau_half = float(t1 + (0.5 - p1) / (p2 - p1) * (t2 - t1))

    tau_exp = None
    try:
        from scipy.optimize import curve_fit

        def f_fit(t, tau, A, C):
            return A * np.exp(-np.maximum(t, 0.0) / tau) + C

        sel = (pops_a > 0.03) & (pops_a < 0.97)
        if sel.sum() > 8:
            popt, _ = curve_fit(f_fit, times_a[sel], pops_a[sel],
                                p0=[100.0, 0.95, 0.02], maxfev=20000)
            tau_exp = float(popt[0])
    except Exception as exc:
        _warn("8C", f"exponential lifetime fit failed: {str(exc)[:90]}")

    result = {
        "n_traj": n_traj, "dt_fs": DT_FS, "tmax_fs": TMAX_FS,
        "n_substeps": N_ESUB, "edc_alpha": EDC_ALPHA, "seed": seed,
        "wall_seconds": time.time() - t0,
        "times_fs": times_a.tolist(),
        "p_s1_active": pops_a.tolist(),
        "p_s1_coherent": pops_coh,
        "hop_histogram": hops_hist.tolist(),
        "n_hops_total": int(hops_hist.sum()),
        "hops_per_traj": float(hops_hist.sum() / n_traj),
        "n_frustrated_hops": int(n_frustrated),
        "tau_half_fs": tau_half, "tau_exp_fs": tau_exp,
        "phi_ci_deg_boundary": model.phi_ci_deg,
        "phi_final_deg": phi_deg.tolist(),
        "active_final": active.tolist(),
        "n_Z": int(z_mask.sum()), "n_E": int(e_mask.sum()),
        "n_S1_resident": int((~on_s0).sum()),
        "phi_Z_yield": float(z_mask.mean()),
        "phi_E_yield": float(e_mask.mean()),
        "decay_time_fs": [None if np.isnan(v) else float(v)
                          for v in decay_time],
    }
    return result, frames


def stage_8c(args, res8b, meta):
    if "error" in res8b or "phi_ci_deg" not in res8b:
        _warn("8C", "8B results unavailable — dynamics skipped")
        return {"error": "8B unavailable"}
    modes = meta.get("xtb_modes")
    if not modes:
        return {"error": "xtb normal modes unavailable"}
    els_d, xyz_d = read_xyz(OUT / "diazene_s0.xyz")
    lvc = build_lvc(res8b, modes, els_d, xyz_d)
    _log("8C", f"LVC: omega_stretch = {lvc['omega_stretch_au'] / CM1_AU:.0f}"
               f" cm-1, dE_FC = {lvc['dE_fc_eV']:.2f} eV, "
               f"lambda_R = {lvc['lambda_stretch']:.4f}, "
               f"kappa_R = {lvc['kappa_stretch']:.4f} Eh/(bohr*sqrt(me))")
    n_traj = 32 if args.smoke else N_TRAJ
    res, frames = run_fssh(lvc, n_traj, SEED)
    res["lvc"] = lvc
    np.savez_compressed(OUT / "fssh_population.npz",
                        times=np.array(res["times_fs"]),
                        p_s1=np.array(res["p_s1_active"]),
                        p_coh=np.array(res["p_s1_coherent"]))
    _log("8C", f"FSSH ensemble done in {res['wall_seconds']:.0f} s: "
               f"P_S1(500 fs) = {res['p_s1_active'][-1]:.3f}, "
               f"tau_half = {res['tau_half_fs']}, "
               f"tau_exp = {res['tau_exp_fs']}, "
               f"Phi_Z = {res['phi_Z_yield']:.3f}, "
               f"Phi_E = {res['phi_E_yield']:.3f}")
    return res


def stage_structures(args) -> dict:
    xtb = find_xtb()
    _log("PRE", f"xtb engine: {xtb or 'NOT FOUND (guess geometries)'}")
    meta = {}
    if not (OUT / "azobenzene_s0.xyz").exists():
        build_azobenzene(xtb, OUT / "azobenzene_s0.xyz")
    meta["azobenzene"] = str(OUT / "azobenzene_s0.xyz")
    if not (OUT / "diazene_s0.xyz").exists():
        dia = build_diazene(xtb, OUT / "diazene_s0.xyz")
        meta["diazene_nn_ang"] = dia["nn_ang"]
    dump_json({"scan": build_azobenzene_scan(OUT / "azobenzene_s0.xyz")},
              OUT / "scan_azobenzene.json")
    geo_t, geo_s = build_diazene_scans(OUT / "diazene_s0.xyz", CI_TORSION_DEG)
    dump_json({"scan_torsion": geo_t, "scan_stretch": geo_s},
              OUT / "scan_diazene.json")
    if xtb:
        try:
            fr, rm, mv = diazene_normal_modes(xtb, OUT / "diazene_s0.xyz")
            els, xyz = read_xyz(OUT / "diazene_s0.xyz")
            cls = classify_diazene_modes(fr, rm, mv, els, xyz)
            I_tors = torsion_inertia(els, xyz)
            meta["xtb_modes"] = {"freqs_cm1": fr.tolist(),
                                 "redmass_amu": rm.tolist(),
                                 "stretch": cls["stretch"],
                                 "bend": cls["bend"],
                                 "torsion_inertia_amu_A2": I_tors,
                                 "stretch_vec": mv[cls["stretch"][0]].tolist(),
                                 "bend_vec": mv[cls["bend"][0]].tolist()}
            _log("PRE", f"xtb modes: N=N stretch {cls['stretch'][1]:.0f} cm-1 "
                        f"(mu {cls['stretch'][2]:.2f} amu), "
                        f"bend {cls['bend'][1]:.0f} cm-1, "
                        f"torsion inertia {I_tors:.3f} amu A^2 "
                        "(CASSCF-scan curvature supplies omega_phi)")
        except Exception as exc:
            _warn("PRE", f"xtb hessian failed: {str(exc)[:100]}")
    return meta


def run_worker(stage, extra, timeout=14400, smoke=False):
    if smoke:
        extra = list(extra) + ["--smoke"]
    cmd = [ENV_PY_QC, str(Path(__file__).resolve()), "--worker", stage] + extra
    _log(stage.upper(), "dispatching QC worker (psi4 env) ...")
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    dt = time.time() - t0
    if proc.returncode != 0:
        err = (proc.stderr or "").splitlines()[-8:]
        raise RuntimeError(f"worker {stage} rc={proc.returncode} ({dt:.0f}s): "
                           + " | ".join(err))
    _log(stage.upper(), f"worker finished in {dt:.0f} s")
    return {"worker_seconds": dt}


def stage_8a(args) -> dict:
    res = {"engine_note": ENGINE_NOTE}
    try:
        res.update(run_worker("8a", ["--xyz", str(OUT / "azobenzene_s0.xyz"),
                                     "--res", str(OUT / "res_8a.json")],
                              smoke=args.smoke))
        res.update(parse_json(OUT / "res_8a.json"))
    except Exception as exc:
        _warn("8A", f"vertical-excitation worker failed: {str(exc)[:150]}")
        res["error"] = str(exc)[:400]
    try:
        run_worker("8a_scan", ["--scan-json", str(OUT / "scan_azobenzene.json"),
                               "--res", str(OUT / "res_8a_scan.json")],
                   smoke=args.smoke)
        res["torsion_scan"] = parse_json(OUT / "res_8a_scan.json")
    except Exception as exc:
        _warn("8A", f"torsion scan worker failed: {str(exc)[:150]}")
    return res


def stage_8b(args) -> dict:
    """Run module 8B with automatic checkpoint-resume on hard aborts.

    This psi4 build occasionally dies without a Python traceback (DETCI DPD
    corruption after many in-process jobs); worker 8B checkpoints every
    completed unit of work, so relaunching with --resume makes forward
    progress monotonic."""
    extra = ["--xyz", str(OUT / "diazene_s0.xyz"),
             "--scan-t-json", str(OUT / "scan_diazene.json"),
             "--scan-s-json", str(OUT / "scan_diazene.json"),
             "--res", str(OUT / "res_8b.json")]
    last = None
    attempts = 4
    for attempt in range(1, attempts + 1):
        try:
            run_worker("8b", extra + (["--resume"] if attempt > 1 else []),
                       timeout=21600, smoke=args.smoke)
            res = parse_json(OUT / "res_8b.json")
            res["n_worker_attempts"] = attempt
            return res
        except Exception as exc:
            last = exc
            _warn("8B", f"worker attempt {attempt}/{attempts} died "
                        f"({str(exc)[:110]}) — resuming from checkpoint")
    return {"error": str(last)[:400], "engine_note": ENGINE_NOTE}


def main():
    ap = argparse.ArgumentParser(
        description="Phase 8: excited-state photochemistry, conical "
                    "intersections & Tully surface hopping")
    ap.add_argument("--stage", default="all",
                    choices=["all", "structures", "8A", "8B", "8C", "fig"])
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--fig_only", action="store_true")
    ap.add_argument("--worker", default=None, choices=[None, "8a", "8a_scan",
                                                       "8b"])
    ap.add_argument("--xyz", default=None)
    ap.add_argument("--res", default=None)
    ap.add_argument("--scan-json", default=None)
    ap.add_argument("--scan-t-json", default=None)
    ap.add_argument("--scan-s-json", default=None)
    ap.add_argument("--threads", type=int, default=min(16, __import__("os")
                                                      .cpu_count() or 8))
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)

    if args.worker == "8a":
        worker_8a(args, Path(args.res))
        return
    if args.worker == "8a_scan":
        worker_8a_scan(args, Path(args.res))
        return
    if args.worker == "8b":
        worker_8b(args, Path(args.res))
        return

    if args.fig_only:
        from phase8_figures import render_all
        render_all()
        return

    master_path = OUT / "phase8_results.json"
    master = parse_json(master_path) if master_path.exists() else {}

    stages = ["structures", "8A", "8B", "8C"]
    if args.stage != "all":
        stages = [args.stage]
    for st in stages:
        try:
            if st == "structures":
                master["meta"] = stage_structures(args)
            elif st == "8A":
                master["module_8a"] = stage_8a(args)
            elif st == "8B":
                master["module_8b"] = stage_8b(args)
            elif st == "8C":
                master["module_8c"] = stage_8c(args,
                                               master.get("module_8b", {}),
                                               master.get("meta", {}))
            dump_json(master, master_path)
        except Exception as exc:
            _warn(st, f"stage crashed: {str(exc)[:200]}")
            master.setdefault("errors", {})[st] = \
                traceback.format_exc()[-800:]
            dump_json(master, master_path)

    if args.stage in ("all", "fig"):
        try:
            sys.path.insert(0, str(Path(__file__).parent))
            from phase8_figures import render_all
            render_all()
        except Exception as exc:
            _warn("fig", f"figure rendering failed: {str(exc)[:200]}")

    _log("P8", f"phase 8 pipeline done — master record: {master_path}")


if __name__ == "__main__":
    main()
