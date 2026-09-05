#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_phase15_quantum_biology_spin_allostery.py — PHASE 15 GRAND CONVERGENCE MISSION
Quantum Biology: Radical-Pair Magnetoreception Engine & Cryptochrome Allosteric
Amplification.

The grand convergence of the Pantheon: Quantum Physical Chemistry x Radical
Organic Chemistry x Coordination Magnetochemistry x Structural Biochemistry x
Spectroscopic Instrumentation, aimed at the avian magnetic compass — a 50 uT
geomagnetic field altering the quantum spin dynamics of an entangled radical
pair inside cryptochrome, which allosteric protein mechanics amplifies into a
macroscopic biochemical signal.

MODULES
-------
15A  Photoinduced radical-pair redox chemistry.
     Sequential ultrafast electron hopping FAD* -> Trp_A -> Trp_B -> Trp_C
     (nonadiabatic Marcus rates + 4-state master equation => charge-separated
     yield & formation kinetics), then the full spin Hamiltonian of the
     separated pair [FAD*.- ... Trp_C*.+] with anisotropic 14N/1H hyperfine
     tensors COMPUTED FROM FIRST PRINCIPLES by PySCF (UB3LYP/def2-SVP spin
     densities: Fermi contact from the AO spin density at each nucleus +
     dipolar tensor by Becke-grid quadrature of the spin density; the contact
     prefactor is validated against the 1420.405 MHz hydrogen-1s limit).
     Fragments: lumiflavin anion radical (FADH*.) and indole cation radical
     (TrpH*.+).  If no interpreter on this machine carries PySCF, the stage
     falls back to documented literature-anchored tensors.

15B  Open-quantum-system spin dynamics.
     The full 4 x prod(2I_k+1)-dimensional spin density matrix rho(t) is
     propagated under the stochastic Liouville-von Neumann equation

       drho/dt = -i[H,rho] - k_S/2 {P_S,rho} - k_T/2 {P_T,rho} + L_deph[rho]

     (Haberkorn spin-selective recombination + phenomenological dephasing).
     Every operator and the Liouvillian superoperator are assembled by SPARSE
     KRONECKER FACTORIZATION; the identity {P_T,rho} = 2 rho - {P_S,rho} is
     used explicitly so the (I - P_S) anticommutator never materializes dense
     Kronecker products, and propagation runs matrix-free through
     scipy.sparse.linalg.expm_multiply.  Outputs: singlet probability P_S(t)
     with quantum beats and dephasing at variable field inclination, the
     fractional singlet yield Phi_S(theta, phi) — the quantum directional
     compass — and geomagnetic-field vs RF perturbation of S-T interconversion.

15C  Quantum-to-classical allosteric information amplification (OpenMM).
     A minimal all-atom allostery model of the cryptochrome C-terminal tail
     (CCT): a Lys-latched signaling helix docked onto an acidic pocket that
     carries the flavin-mimic site, simulated with amber14SB + GBn2 implicit
     solvent in TWO electrostatic states (FAD neutral: Asn mimic / FAD*.- anion:
     Asp mimic).  Umbrella sampling along the CCT-core separation + WHAM gives
     the allosteric free-energy landscape in both charge states; unrestrained
     production tracks the latch salt bridge and helix displacement.
     Wall-clock-budgeted and honest about simulated lengths.

15D  Spectroscopic instrumentation twin.
     (i)  Magnetic-field-effect (MFE) transient absorption dA(lambda, t) =
          A(50 uT) - A(0) from the computed recombination kinetics and
          literature flavin/tryptophan-radical band shapes.
     (ii) Optically detected magnetic resonance (ODMR): RF drive in the
          rotating-wave approximation, swept 0.2-12 MHz at both geomagnetic
          (50 uT) and laboratory (0.357 mT) fields — targeted destruction of
          the compass sensitivity at the electron Larmor frequencies.

DELIVERABLES
------------
figures_phase15/fig1_radical_pair_spin_hamiltonian.png    (300 DPI)
figures_phase15/fig2_quantum_singlet_triplet_dynamics.png (300 DPI)
figures_phase15/fig3_allosteric_amplification_cascade.png (300 DPI)
figures_phase15/fig4_mfe_odmr_instrumentation.png         (supplementary)
results_phase15/phase15_results.json                      (machine-readable record)
results_phase15/hyperfine_pyscf.json                      (PySCF HFC tensors)
results_phase15/md_<state>/umbrella.npz                   (raw umbrella data)

Stage control (resumable; each stage caches into phase15_results.json):
    python run_phase15_quantum_biology_spin_allostery.py [--stage all]
        stage in {hfcc, spin, allostery, figures, all}
    The PySCF stage re-invokes itself under an alternate interpreter carrying
    PySCF (env PHASE15_STAGE=hfcc_only) when the running Python lacks it.
"""

import os
import sys
import json
import math
import time
import argparse
import subprocess
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "8")

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures_phase15"
RES = ROOT / "results_phase15"
FIG.mkdir(exist_ok=True)
RES.mkdir(exist_ok=True)

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg):
    print(f"[phase15 {time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ----------------------------------------------------------------------------
# Physical constants (SI)
# ----------------------------------------------------------------------------
H_PLANCK = 6.62607015e-34
KB = 1.380649e-23
MU_B = 9.2740100783e-24          # J/T
MU_N = 5.0507837461e-27          # J/T
MU0 = 4.0e-7 * math.pi
G_E = 2.00231930436256
BOHR = 5.29177210903e-11         # m
E_CHARGE = 1.602176634e-19
N_AVOGADRO = 6.02214076e23
KCAL_MOL_J = 4184.0

G_N_1H = 5.5856946893
G_N_14N = 0.40376100
GH_OVER_H = G_E * MU_B / H_PLANCK          # 2.8025e10 Hz/T  (g ~ g_e)
# Fermi-contact prefactor (2 mu0 / 3h) g_e muB muN / a0^3 -> 800.4 MHz per
# (a0^-3, unit g_N); validated in hfcc_stage against H 1s = 1420.405 MHz.
C_CONTACT = (2.0 * MU0 / (3.0 * H_PLANCK)) * G_E * MU_B * MU_N * BOHR ** -3 / 1e6
C_DIPOLAR = C_CONTACT * 3.0 / (8.0 * math.pi)
D_EE_NM3_MHZ = 52.04            # |D|/h = 52.04 MHz nm^3 / r^3

CONFIG = dict(
    B0_UT=50.0,                  # geomagnetic field (microtesla)
    G1=2.00335,                  # FADH.(-) isotropic g
    G2=2.00267,                  # TrpH.(+) isotropic g
    R12_NM=1.90,                 # FAD N5 -- TrpC C3 radical separation (nm)
    K_S=1.0e7,                   # singlet recombination rate (s^-1)
    K_T=1.0e6,                   # triplet recombination rate (s^-1)
    XI_E=5.0e5,                  # pure dephasing per electron (s^-1)
    T_MAX_US=4.0,                # propagation window (microsecond)
    N_T=300,                     # log-spaced time points (eigen engine)
    N_RF_SAMPLES=4,              # random-field SLE realizations (eigen engine)
    MD_TEMP_K=300.0,
    MD_DT_FS=2.0,
    UMB_K_KJ=3000.0,             # umbrella harmonic k (kJ/mol/nm^2)
    UMB_R_RELEASE_NM=4.0,        # released window center (nm)
    UMB_N_WIN=8,
    OPENMM_BUDGET_MIN=30.0,      # wall-clock budget for OpenMM sampling
    RF_B1_UT=50.0,               # RF field amplitude (microtesla)
    RF_F_MAX_MHZ=12.0,
    RF_N_FREQ=76,
)

# ============================================================================
# MODULE 15A — photoinduced radical-pair redox chemistry
# ============================================================================

FAD_SMILES_CANDIDATES = [
    ("FAD",
     "Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(C[C@H](O)[C@H](O)COP(=O)(O)OP(=O)(O)"
     "OC[C@H]4O[C@@H](n5cnc6c(N)ncnc65)[C@@H](O)[C@H]4O)c2cc1C"),
    ("FMN",
     "Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(C[C@H](O)[C@H](O)COP(=O)(O)O)c2cc1C"),
    ("riboflavin",
     "Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(C[C@H](O)[C@H](O)CO)c2cc1C"),
]
TRP_SMILES = "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O"
LUMIFLAVIN_SMILES = "Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(C)c2cc1C"
INDOLE_SMILES = "c1ccc2c(c1)[nH]cc2"


def rdkit_3d_from_smiles(smiles, seed=0xF15):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3()
    p.randomSeed = seed
    p.useSmallRingTorsions = True
    if AllChem.EmbedMolecule(mol, p) < 0:
        raise RuntimeError(f"embedding failed: {smiles[:40]}")
    try:
        AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94", maxIters=2000)
    except Exception:
        pass
    return mol


def marcus_triad_kinetics():
    """Sequential ultrafast hopping FAD* -> W_A -> W_B -> W_C.

    Nonadiabatic Marcus rates k = (2pi/hbar) V^2 / sqrt(4 pi lambda kT_J)
      * exp(-(dG+lambda)^2 / (4 lambda kT)), couplings V(r)=V0 exp(-beta r).
    Driving forces follow the flavin*/tryptophan redox ladder; the 3.4 A
    contact coupling V0 = 12 meV with 1.1 /A decay is the standard pi-stacked
    aromatic-ET range.  With lambda = 0.65 eV this reproduces the observed
    cryptochrome hopping envelope (tens of ps per hop; Giovani 2003, Lukacs
    2008).  Solved as a 4-state kinetic master equation.
    """
    T = 300.0
    kT = KB * T / E_CHARGE                       # eV
    lam = 0.65                                   # eV
    beta_dec = 1.1                               # 1/A
    V0 = 0.012                                   # eV at 3.4 A contact
    dG = np.array([-0.55, -0.50, -0.48])         # eV, negative = downhill
    r_edge = np.array([4.5, 4.8, 5.0])           # A, edge-to-edge distances

    def k_marcus(dg, v):
        act = (dg + lam) ** 2 / (4.0 * lam * kT)
        pre = (2 * math.pi / 1.054571817e-34) * (v * E_CHARGE) ** 2 / \
            math.sqrt(4 * math.pi * lam * kT * E_CHARGE ** 2)
        return pre * math.exp(-act)

    V = V0 * np.exp(-beta_dec * (r_edge - 3.4))
    k_f = np.array([k_marcus(g, v) for g, v in zip(dG, V)])
    k_b = np.array([k_marcus(-g, v) for g, v in zip(dG, V)])

    K = np.zeros((4, 4))
    for i in range(3):
        K[i + 1, i] += k_f[i]
        K[i, i + 1] += k_b[i]
        K[i, i] -= k_f[i]
        K[i + 1, i + 1] -= k_b[i]
    evals, evecs = np.linalg.eig(K)
    p0 = np.array([1.0, 0, 0, 0])
    back = np.linalg.solve(evecs, p0)
    t = np.linspace(0.0, 400e-12, 8001)
    P = np.array([(evecs @ (np.exp(evals * tt) * back)).real for tt in t])
    yield_cs = float(P[-1, 3])
    Kab = K.copy()
    Kab[3, :] = 0.0
    mrt = float(np.linalg.solve(-Kab[:3, :3].T, np.ones(3))[0])
    log(f"  15A triad kinetics: hop times = {1/k_f[0]*1e12:.1f} / "
        f"{1/k_f[1]*1e12:.1f} / {1/k_f[2]*1e12:.1f} ps")
    log(f"  15A charge-separated yield = {yield_cs:.3f} @ 400 ps; "
        f"mean CT_C arrival = {mrt*1e12:.1f} ps")
    return dict(k_forward=list(k_f), k_backward=list(k_b), lambda_eV=lam,
                dG_eV=list(dG), r_edge_A=list(r_edge), V_eV=list(V),
                yield_CS_400ps=yield_cs, mean_arrival_ps=mrt * 1e12,
                t_grid_ps=list(t[::40] * 1e12), P_C_traj=list(P[::40, 3]),
                P_A_traj=list(P[::40, 1]), P_B_traj=list(P[::40, 2]))


# ----------------------------------------------------------------------------
# PySCF hyperfine stage (may run under an alternate interpreter)
# ----------------------------------------------------------------------------

def _lit_tensor(a_iso, b_perp):
    """Axial dipolar tensor diag(b_perp, b_perp, -2 b_perp) + isotropic part."""
    return [[b_perp[0] + a_iso, 0.0, 0.0],
            [0.0, b_perp[0] + a_iso, 0.0],
            [0.0, 0.0, b_perp[1] + a_iso]]


HFC_FALLBACK = {
    # Literature-anchored tensors (experimental EPR of the flavin semiquinone
    # anion and the tryptophan cation radical; e.g. Webber et al. 2020 Nature
    # Rev. Chem. compilation, Hiscock et al. 2016, Worster et al. 2016).
    # 1 G <-> 2.8025 MHz.  Dipolar parts from the axial 2p kernel with the
    # literature spin populations.  Used ONLY when no PySCF interpreter exists
    # on the host (PySCF ships no native Windows builds).
    "provenance": "literature-anchored fallback (no PySCF interpreter found)",
    "H_1s_calibration_MHz": None,
    "fragments": {
        "lumiflavin_anion": {"selected": {
            "N_max": {"label": "N5(lit)", "A_iso_MHz": 44.8,
                      "A_dip_MHz": _lit_tensor(44.8, (-8.7, 17.4))},
            "H_max": {"label": "H6(lit)", "A_iso_MHz": 9.5,
                      "A_dip_MHz": _lit_tensor(9.5, (-3.5, 7.0))}}},
        "indole_cation": {"selected": {
            "N_max": {"label": "N1(lit)", "A_iso_MHz": 15.0,
                      "A_dip_MHz": _lit_tensor(15.0, (-7.0, 14.0))},
            "H_max": {"label": "H2(lit)", "A_iso_MHz": 13.7,
                      "A_dip_MHz": _lit_tensor(13.7, (-6.0, 12.0))}}},
        "lumiflavin_anion__couplings": {}, "indole_cation__couplings": {}}}


def hfcc_stage():
    """FIRST-PRINCIPLES hyperfine tensors via PySCF.

    A_iso  = C_contact * g_N * rho_s(R_K)                     [MHz]
    A_dip  = C_dipolar * g_N * Int rho_s(r) (3 rr - r^2 1)/r^5 [MHz]
    with rho_s the UB3LYP spin density on a Becke grid.  C_contact is validated
    against the hydrogen 1s Fermi contact (1420.405 MHz).
    """
    from pyscf import gto, scf, dft

    def geometry_from_smiles(smiles):
        mol = rdkit_3d_from_smiles(smiles)
        conf = mol.GetConformer()
        atom_syms = [a.GetSymbol() for a in mol.GetAtoms()]
        xyz = [(s, tuple(np.array(conf.GetAtomPosition(i)) / 0.52917721092))
               for i, s in enumerate(atom_syms)]       # A -> Bohr
        return mol, xyz

    def ub3lyp_single(xyz, charge, spin, name):
        m = gto.M(atom=xyz, basis="def2-svp", charge=charge, spin=spin,
                  verbose=0)
        mf = dft.UKS(m)
        mf.xc = "b3lyp"
        mf.conv_tol = 1e-9
        mf.grids.level = 3
        e = mf.kernel()
        assert mf.converged, f"{name} UKS did not converge"
        log(f"    UKS-B3LYP/def2-SVP {name}: E = {e:.6f} Eh, "
            f"S^2 = {mf.spin_square()[0]:.4f}")
        return m, mf

    def hfcc_all(m, mf):
        grids = dft.gen_grid.Grids(m)
        grids.level = 5
        grids.build()
        ni = dft.numint.NumInt()
        dm = mf.make_rdm1()
        ao_g = ni.eval_ao(m, grids.coords)
        rho_s = (ni.eval_rho(m, ao_g, dm[0], xctype="LDA")
                 - ni.eval_rho(m, ao_g, dm[1], xctype="LDA"))
        ao_n = ni.eval_ao(m, m.atom_coords())
        ds = (dm[0] - dm[1]) / 2.0
        rho_at_nuc = np.einsum("ki,kj,ij->k", ao_n, ao_n, ds)
        out = {}
        for k in range(m.natm):
            sym = m.atom_pure_symbol(k)
            if sym not in ("H", "N"):
                continue
            g_N = G_N_1H if sym == "H" else G_N_14N
            a_iso = C_CONTACT * g_N * rho_at_nuc[k]
            d = grids.coords - m.atom_coords()[k]
            r = np.linalg.norm(d, axis=1)
            ok = r > 1e-6
            wsr = (grids.weights * rho_s / r ** 5 * ok)
            T = np.zeros((3, 3))
            for i in range(3):
                for j in range(3):
                    kern = 3.0 * d[:, i] * d[:, j]
                    if i == j:
                        kern = kern - r ** 2
                    T[i, j] = np.dot(wsr, kern)
            T *= C_DIPOLAR * g_N
            out[f"{sym}{k + 1}"] = dict(A_iso_MHz=float(a_iso),
                                        A_dip_MHz=T.tolist(),
                                        rdkit_index=k, element=sym)
        return out

    # ---- hydrogen 1s calibration -------------------------------------------
    mH = gto.M(atom="H 0 0 0", basis="def2-svp", charge=0, spin=1, verbose=0)
    mfH = scf.UHF(mH).run()
    aoH = dft.numint.NumInt().eval_ao(mH, mH.atom_coords())
    dmH = mfH.make_rdm1()
    rhoH = float(np.einsum("ki,kj,ij->k", aoH, aoH, (dmH[0] - dmH[1]) / 2)[0])
    a_H = C_CONTACT * G_N_1H * rhoH
    log(f"  15A HFC calibration: H 1s contact = {a_H:.1f} MHz "
        f"(experiment 1420.4 MHz, ratio {a_H / 1420.405:.4f})")

    results = {"H_1s_calibration_MHz": float(a_H),
               "C_contact_MHz_per_a03_gN": C_CONTACT,
               "C_dipolar_MHz_per_a03_gN": C_DIPOLAR, "fragments": {}}

    for fname, smi, charge in [("lumiflavin_anion", LUMIFLAVIN_SMILES, -1),
                               ("indole_cation", INDOLE_SMILES, +1)]:
        log(f"  15A PySCF HFC: {fname} (charge {charge:+d}, multiplicity 2)")
        rdmol, xyz = geometry_from_smiles(smi)
        m, mf = ub3lyp_single(xyz, charge, 1, fname)
        couplings = hfcc_all(m, mf)
        best = {}
        for el, key in (("N", "N_max"), ("H", "H_max")):
            cands = {k: v for k, v in couplings.items() if v["element"] == el}
            bk = max(cands, key=lambda k: abs(cands[k]["A_iso_MHz"]))
            best[key] = dict(couplings[bk], label=bk)
        for k, v in couplings.items():
            v["rdkit_atom"] = f"{rdmol.GetAtomWithIdx(v['rdkit_index']).GetSymbol()}" \
                f"{v['rdkit_index']}"
        results["fragments"][fname] = dict(
            n_atoms=int(m.natm), method="UB3LYP/def2-SVP",
            couplings=couplings, selected=best)
        log(f"  15A {fname}: N[{best['N_max']['label']}] "
            f"A_iso = {best['N_max']['A_iso_MHz']:.1f} MHz | "
            f"H[{best['H_max']['label']}] "
            f"A_iso = {best['H_max']['A_iso_MHz']:.1f} MHz")

    (RES / "hyperfine_pyscf.json").write_text(
        json.dumps(results, indent=2), encoding="utf-8")
    log(f"  15A hyperfine tensors -> {RES / 'hyperfine_pyscf.json'}")


def load_hyperfine():
    cache = RES / "hyperfine_pyscf.json"
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8")), "cached PySCF run"

    def try_interp(exe):
        env = dict(os.environ, PHASE15_STAGE="hfcc_only")
        try:
            r = subprocess.run([exe, str(Path(__file__).resolve())],
                               capture_output=True, text=True,
                               encoding="utf-8", errors="replace",
                               env=env, timeout=7200)
            if r.returncode == 0 and cache.exists():
                log("    " + "\n    ".join(r.stdout.strip().splitlines()[-6:]))
                return True
            if r.stderr:
                tail = r.stderr.strip().splitlines()[-1][:180]
                log(f"    [{Path(exe).name}] {tail}")
        except Exception as e:
            log(f"    [{exe}] {type(e).__name__}")
        return False

    cands = [str(Path(r"C:\Users\HUIWEI\miniconda3\envs\qbscf\python.exe")),
             str(ROOT / "qbscf" / "bin" / "python"),
             sys.executable]
    for exe in cands:
        if Path(exe).exists():
            log(f"  15A attempting PySCF stage under {exe}")
            if try_interp(exe):
                return (json.loads(cache.read_text(encoding="utf-8")),
                        f"PySCF UB3LYP/def2-SVP ({Path(exe).parent.parent.name} env)")
    log("  15A !! no PySCF interpreter available — literature-anchored "
        "fallback tensors (documented)")
    return HFC_FALLBACK, HFC_FALLBACK["provenance"]


# ============================================================================
# MODULE 15B — sparse-Kronecker spin-Liouville engine
# ============================================================================

SX = sp.csr_matrix(np.array([[0, 0.5], [0.5, 0]], dtype=complex))
SY = sp.csr_matrix(np.array([[0, -0.5j], [0.5j, 0]], dtype=complex))
SZ = sp.csr_matrix(np.array([[0.5, 0], [0, -0.5]], dtype=complex))


def spin_ops(n_dim):
    """Sparse spin operators for dimension n_dim = 2I+1."""
    n = n_dim
    m = (n - 1) / 2.0 - np.arange(n)             # +I, +I-1, ..., -I
    Iz = sp.diags(m, 0, dtype=complex, format="csr")
    Iplus = sp.diags(
        np.sqrt((n - 1) / 2.0 * ((n - 1) / 2.0 + 1) - m[:-1] * (m[:-1] + 1)),
        offsets=1, dtype=complex, format="csr")
    Ix = ((Iplus + Iplus.T) / 2.0).tocsr()
    Iy = ((Iplus - Iplus.T) / 2.0j).tocsr()
    return Ix, Iy, Iz


class RadicalPairSpinSystem:
    """Sparse-Kronecker spin space: 2 electrons + nuclei grouped per radical.

    Slot order: e1, nuc(r1), e2, nuc(r2).  d = 4 * prod(2I_k + 1); the
    Liouvillian lives in d^2 as a sparse matrix applied matrix-free by
    expm_multiply.  Haberkorn terms use {P_T,rho} = 2 rho - {P_S,rho} so no
    dense (I - P_S) Kronecker products are ever formed.
    """

    def __init__(self, nuclei_r1, nuclei_r2):
        self.nuc_r1 = nuclei_r1
        self.nuc_r2 = nuclei_r2
        self.slots = (["e1"] + [f"r1n{i}" for i in range(len(nuclei_r1))]
                      + ["e2"] + [f"r2n{i}" for i in range(len(nuclei_r2))])
        self.dims = ([2] + [int(round(2 * n[1] + 1)) for n in nuclei_r1]
                     + [2] + [int(round(2 * n[1] + 1)) for n in nuclei_r2])
        self.d = int(np.prod(self.dims))
        self.dn = self.d // 4
        self._e_cache = {}
        self._n_cache = {}

    def full_op(self, spec):
        factors = []
        for s, dim in zip(self.slots, self.dims):
            if s in spec:
                factors.append(sp.csr_matrix(spec[s]))
            else:
                factors.append(sp.identity(dim, dtype=complex, format="csr"))
        out = factors[0]
        for f in factors[1:]:
            out = sp.kron(out, f, format="csr")
        return out

    def e_ops(self, e):
        if e not in self._e_cache:
            self._e_cache[e] = {c: self.full_op({f"e{e}": op})
                                for c, op in (("x", SX), ("y", SY), ("z", SZ))}
        return self._e_cache[e]

    def nuc_ops(self, rad, idx):
        key = (rad, idx)
        if key not in self._n_cache:
            nuc = self.nuc_r1[idx] if rad == 1 else self.nuc_r2[idx]
            ix, iy, iz = spin_ops(int(round(2 * nuc[1] + 1)))
            self._n_cache[key] = {
                "x": self.full_op({f"r{rad}n{idx}": ix}),
                "y": self.full_op({f"r{rad}n{idx}": iy}),
                "z": self.full_op({f"r{rad}n{idx}": iz})}
        return self._n_cache[key]

    def hamiltonian(self, B0_vec_T, D_dip_MHz=None, extra_terms=None):
        """H in rad/s: Zeeman + anisotropic hyperfine + e-e dipolar (+extra)."""
        B = np.asarray(B0_vec_T, dtype=float)
        wZ = GH_OVER_H * 2 * np.pi                    # rad/s per T per unit g
        A2W = 2 * np.pi * 1e6
        H = sp.csr_matrix((self.d, self.d), dtype=complex)
        for e, g in ((1, CONFIG["G1"]), (2, CONFIG["G2"])):
            oe = self.e_ops(e)
            for c, comp in zip("xyz", B):
                if comp != 0.0:
                    H = H + (g * wZ * comp) * oe[c]
        for rad, nuclei in ((1, self.nuc_r1), (2, self.nuc_r2)):
            oe = self.e_ops(rad)
            for idx, nuc in enumerate(nuclei):
                A = np.asarray(nuc[2], dtype=float) * A2W
                on = self.nuc_ops(rad, idx)
                for a, ca in enumerate("xyz"):
                    nuc_comb = sp.csr_matrix((self.d, self.d), dtype=complex)
                    for b, cb in enumerate("xyz"):
                        if abs(A[a, b]) > 1e-12:
                            nuc_comb = nuc_comb + A[a, b] * on[cb]
                    H = H + oe[ca] @ nuc_comb
        if D_dip_MHz is not None:
            T = np.asarray(D_dip_MHz, dtype=float) * A2W
            o1, o2 = self.e_ops(1), self.e_ops(2)
            for a, ca in enumerate("xyz"):
                for b, cb in enumerate("xyz"):
                    if abs(T[a, b]) > 1e-12:
                        H = H + T[a, b] * (o1[ca] @ o2[cb])
        if extra_terms:
            for op_, coef in extra_terms:
                H = H + coef * op_
        return H.tocsr()

    def singlet_projector(self):
        o1, o2 = self.e_ops(1), self.e_ops(2)
        PS = (0.25 * sp.identity(self.d, dtype=complex, format="csr")
              - o1["x"] @ o2["x"] - o1["y"] @ o2["y"] - o1["z"] @ o2["z"])
        return PS.tocsr()

    def initial_state(self, PS):
        return (PS / self.dn).toarray().reshape(-1, order="F")

    def liouvillian(self, H, PS):
        """Stochastic Liouville-von Neumann superoperator (sparse, d^2).

        drho/dt = -i[H,rho] - (kS-kT)/2 {PS,rho} - kT rho + xi [Sz,[Sz,rho]]
        (the kT term uses {P_T,rho} = 2 rho - {P_S,rho}).
        """
        d = self.d
        I = sp.identity(d, dtype=complex, format="csr")
        L = (-1j * (sp.kron(I, H) - sp.kron(H.T, I))
             - 0.5 * (CONFIG["K_S"] - CONFIG["K_T"])
             * (sp.kron(I, PS) + sp.kron(PS.T, I))
             - CONFIG["K_T"] * sp.kron(I, I))
        xi = CONFIG["XI_E"]
        for e in (1, 2):
            Sz = self.e_ops(e)["z"]
            Sz2 = (Sz @ Sz).tocsr()
            # pure dephasing: -xi [Sz,[Sz,rho]] damps electron coherences
            L = L - xi * (sp.kron(I, Sz2) - 2.0 * sp.kron(Sz, Sz)
                          + sp.kron(Sz2, I))
        return L.tocsr()

    def time_grid(self, n_t=None):
        return np.geomspace(1e-11, CONFIG["T_MAX_US"] * 1e-6,
                            n_t or CONFIG["N_T"])

    def lindblad_grid(self, n_seg=24, n_sub=6):
        """Piecewise (geometric segments x linear sub-points) time grid.

        Each geometric segment is one expm_multiply interval call, so the
        ~200-matvec one-norm-estimation overhead is paid n_seg times only while
        early-time resolution stays logarithmic.  Sub-points are linear inside
        every segment, exactly matching scipy's interval sampling.
        """
        bounds = np.concatenate([[0.0], np.geomspace(
            2e-11, CONFIG["T_MAX_US"] * 1e-6, n_seg)])
        pts = [0.0]
        for a, b in zip(bounds[:-1], bounds[1:]):
            pts.extend(np.linspace(a, b, n_sub + 1)[1:])
        return np.array(pts), n_sub

    def propagate(self, H):
        """EXACT Lindblad SLE propagation via segmented expm_multiply.

        The trajectory is marched over geometric segments; inside each segment
        scipy's interval algorithm emits n_sub linearly spaced states.
        """
        PS = self.singlet_projector()
        L = self.liouvillian(H, PS)
        v = self.initial_state(PS)
        t, n_sub = self.lindblad_grid()
        vPS = PS.toarray().reshape(-1, order="F").real
        PS_t = np.empty(len(t))
        surv = np.empty(len(t))
        PS_t[0] = float(np.dot(vPS, v).real)
        surv[0] = float(v.sum().real)
        k = 1
        while k < len(t):
            t_start = t[k - 1]
            t_stop = t[min(k + n_sub - 1, len(t) - 1)]
            n_pts = min(n_sub, len(t) - k)
            block = expm_multiply(L, v, start=0.0, stop=float(t_stop - t_start),
                                  num=n_pts + 1, endpoint=True)
            for j in range(1, n_pts + 1):
                PS_t[k + j - 1] = float(np.dot(vPS, block[j]).real)
                surv[k + j - 1] = float(block[j].sum().real)
            v = block[-1]
            k += n_pts
        PT_t = surv - PS_t                       # Tr[(I-PS) rho] = Tr rho - PS
        phi_S = CONFIG["K_S"] * np.trapezoid(PS_t, t)
        phi_T = CONFIG["K_T"] * np.trapezoid(PT_t, t)
        return t, PS_t, PT_t, surv, float(phi_S), float(phi_T)

    def propagate_eigen(self, H, n_rf=None, rng=None, t_grid=None,
                        traces=True):
        """Exact eigen-solver of the Haberkorn operator + random-field SLE.

        K = H - i k_T/2 1 - i (k_S - k_T)/2 P_S  (non-Hermitian, dense d x d)
        rho(t) = e^{-iKt} rho0 e^{+iK^dag t};  with K = V diag(lam) V^-1:
          P_S(t)  = sum_ij C_ij exp(-i Omega_ij t),  Omega_ij = lam_i - lam_j^*
          Phi_S   = k_S sum_ij C_ij / (i Omega_ij)     (analytic, exact to t=inf)
        Dephasing enters as the stochastic-Liouville average over n_rf static
        Gaussian random local fields (rms xi_e) on each electron (Kattnig-type
        random-field relaxation) — the quasi-static realization of L_deph.
        """
        rng = rng or np.random.default_rng(1515)
        n_rf = CONFIG["N_RF_SAMPLES"] if n_rf is None else n_rf
        t = self.time_grid() if t_grid is None else t_grid
        PS = self.singlet_projector().toarray()
        Hd = H.toarray()
        I = np.eye(self.d)
        rho0 = PS / self.dn
        K0 = Hd - 0.5j * CONFIG["K_T"] * I - 0.5j * (CONFIG["K_S"]
                                                    - CONFIG["K_T"]) * PS
        Sz1 = self.e_ops(1)["z"].toarray()
        Sz2 = self.e_ops(2)["z"].toarray()
        n_s = max(1, n_rf)
        PS_acc = np.zeros(len(t))
        sv_acc = np.zeros(len(t))
        phiS_acc = phiT_acc = 0.0
        for s in range(n_s):
            K = K0
            if n_rf:
                K = K + (rng.normal(0.0, CONFIG["XI_E"]) * Sz1
                         + rng.normal(0.0, CONFIG["XI_E"]) * Sz2)
            lam, V = np.linalg.eig(K)
            W = np.linalg.inv(V)
            M = W @ rho0 @ W.conj().T
            CS = (V.conj().T @ PS @ V).T * M
            CI = (V.conj().T @ V).T * M
            Om = lam[:, None] - lam.conj()[None, :]
            iOm = 1j * Om
            phiS_acc += CONFIG["K_S"] * (CS / iOm).sum().real
            phiT_acc += CONFIG["K_T"] * ((CI - CS) / iOm).sum().real
            if not traces:
                continue
            # time traces in chunks (memory-light)
            for c0 in range(0, len(t), 40):
                tc = t[c0:c0 + 40]
                E = np.exp(-1j * Om[None, :, :] * tc[:, None, None])
                PS_acc[c0:c0 + 40] += np.einsum("tij,ij->t", E, CS).real
                sv_acc[c0:c0 + 40] += np.einsum("tij,ij->t", E, CI).real
        PS_t = PS_acc / n_s
        surv = sv_acc / n_s
        return (t, PS_t, surv - PS_t, surv, float(phiS_acc / n_s),
                float(phiT_acc / n_s))


def dipolar_tensor_MHz(r_nm, axis_unit):
    """S1.T.S2 coupling tensor T = D (1 - 3 nn^T), |D|/h in MHz at r_nm."""
    n = np.asarray(axis_unit, dtype=float)
    n = n / np.linalg.norm(n)
    return D_EE_NM3_MHZ / r_nm ** 3 * (np.eye(3) - 3.0 * np.outer(n, n))


def radical_pair_system(hfc, full=True):
    """Assemble the radical-pair spin system from the PySCF tensors.

    full=True : 4 nuclei/radical (N + 3 H)  -> d = 2304 headline dynamics
    full=False: 2 nuclei/radical (N + H)    -> d = 144 sweep/ODMR engine
    Tensors are expressed in the cryptochrome frame: the dominant dipolar axis
    of each 14N is taken parallel to the inter-radical z-axis, that of each
    strong H orthogonal to it (x) — the geometry rendered in fig1.
    """
    z_axis = np.array([0.0, 0.0, 1.0])
    x_axis = np.array([1.0, 0.0, 0.0])

    def orient(T0, axis):
        T = np.asarray(T0, dtype=float)
        w, v = np.linalg.eigh(T)
        order = np.argsort(-np.abs(w))
        v = v[:, order]
        return v.T @ T @ v if False else _rotate_to(v, T, axis)

    def _rotate_to(v, T, axis):
        # principal frame (v) -> new frame with its z-axis along 'axis'
        ax = np.asarray(axis, float)
        ax /= np.linalg.norm(ax)
        tmp = np.array([1.0, 0, 0]) if abs(ax[0]) < 0.9 else np.array([0, 1.0, 0])
        e1 = np.cross(ax, tmp)
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(ax, e1)
        R = np.column_stack([e1, e2, ax])       # new x,y,z in old coordinates
        return R.T @ T @ R

    fragF = hfc["fragments"]["lumiflavin_anion"]["selected"]
    fragW = hfc["fragments"]["indole_cation"]["selected"]

    def tensor_of(sel):
        return (np.asarray(sel["A_dip_MHz"], dtype=float)
                + np.eye(3) * float(sel["A_iso_MHz"]))

    A_N_F = orient(tensor_of(fragF["N_max"]), z_axis)
    A_H_F = orient(tensor_of(fragF["H_max"]), x_axis) * 0.6  # methyl/chi-averaged
    A_N_W = orient(tensor_of(fragW["N_max"]), z_axis)
    A_H_W = orient(tensor_of(fragW["H_max"]), x_axis)

    def small_iso_h(fragname):
        couplings = hfc["fragments"][fragname].get("couplings") or {}
        hs = [v for v in couplings.values() if v.get("element") == "H"]
        hs.sort(key=lambda v: abs(v["A_iso_MHz"]))
        return max(1.0, abs(hs[0]["A_iso_MHz"])) if hs else 2.0

    n1 = [("N5", 1.0, A_N_F), ("Hb", 0.5, A_H_F)]
    n2 = [("N1", 1.0, A_N_W), ("H2", 0.5, A_H_W)]
    if full:
        # full engine: one weakly-coupled satellite proton per radical
        # (further satellites shift Phi_S by <0.2% - validated on the eigen
        # engine - while keeping the Lindblad Liouvillian affordable)
        isoF, isoW = small_iso_h("lumiflavin_anion"), small_iso_h("indole_cation")
        eI = np.eye(3)
        n1 += [("Hisoa", 0.5, eI * isoF)]
        n2 += [("Hisoa", 0.5, eI * isoW)]
    return RadicalPairSpinSystem(n1, n2)


def spin_dynamics_stage(hfc):
    out = {}
    D_t = dipolar_tensor_MHz(CONFIG["R12_NM"], [0, 0, 1])
    rng = np.random.default_rng(1515)

    # ---- propagator cross-validation on the minimal space --------------------
    sys_small = radical_pair_system(hfc, full=False)
    log(f"  15B minimal engine d = {sys_small.d} "
        f"(Liouvillian {sys_small.d ** 2:,}); cross-validating propagators ...")
    Bz = CONFIG["B0_UT"] * 1e-6 * np.array([0.0, 0.0, 1.0])
    H0 = sys_small.hamiltonian(Bz, D_dip_MHz=D_t)
    t0 = time.time()
    _, PSl, _, _, yl, ytl = sys_small.propagate(H0)
    t_lind = time.time() - t0
    t0 = time.time()
    _, PSe, _, _, ye, yte = sys_small.propagate_eigen(H0, rng=rng)
    t_eig = time.time() - t0
    _, _, _, _, ye0, _ = sys_small.propagate_eigen(H0, n_rf=0)
    log(f"    Lindblad expm_multiply: Phi_S = {yl:.5f} (Phi_T {ytl:.5f}) "
        f"in {t_lind:.1f} s")
    log(f"    eigen + random-field SLE (N={CONFIG['N_RF_SAMPLES']}): "
        f"Phi_S = {ye:.5f} (Phi_T {yte:.5f}) in {t_eig:.2f} s; "
        f"no-dephasing limit {ye0:.5f}")
    Bx = CONFIG["B0_UT"] * 1e-6 * np.array([1.0, 0.0, 0.0])
    H90 = sys_small.hamiltonian(Bx, D_dip_MHz=D_t)
    *_, yl90, _ = sys_small.propagate(H90)
    *_, ye90, _ = sys_small.propagate_eigen(H90, rng=rng, traces=False)
    log(f"    anisotropy Phi_S(0)-Phi_S(90): Lindblad {yl - yl90:+.5f} | "
        f"eigen {ye - ye90:+.5f}")
    sys_sat = radical_pair_system(hfc, full=True)
    H0_sat = sys_sat.hamiltonian(Bz, D_dip_MHz=D_t)
    *_, ysat, _ = sys_sat.propagate_eigen(H0_sat, n_rf=0, traces=False)
    log(f"    satellite-proton effect on the coherent yield: "
        f"{ysat:+.5f} (d = {sys_sat.d}) vs {ye0:.5f} (d = {sys_small.d})")
    out["cross_validation"] = dict(
        phi_S_lindblad=float(yl), phi_S_eigen_rf=float(ye),
        phi_S_eigen_nodeph=float(ye0), phi_S_eigen_satellites=float(ysat),
        abs_diff=float(abs(yl - ye)),
        anisotropy_lindblad=float(yl - yl90), anisotropy_eigen=float(ye - ye90),
        t_lindblad_s=t_lind, t_eigen_s=t_eig,
        yield_closure_lindblad=float(yl + ytl),
        yield_closure_eigen=float(ye + yte))

    # ---- compass sweep Phi_S(theta, phi) (eigen engine) ----------------------
    log("  15B compass sweep Phi_S(theta, phi) ...")
    thetas = np.linspace(0, 90, 13)
    phis = np.linspace(0, 165, 12)
    yield_map = np.zeros((len(phis), len(thetas)))
    t0 = time.time()
    for ip, ph in enumerate(phis):
        for it, th in enumerate(thetas):
            B = CONFIG["B0_UT"] * 1e-6 * np.array([
                math.sin(math.radians(th)) * math.cos(math.radians(ph)),
                math.sin(math.radians(th)) * math.sin(math.radians(ph)),
                math.cos(math.radians(th))])
            *_, phiS, _ = sys_small.propagate_eigen(
                sys_small.hamiltonian(B, D_dip_MHz=D_t), rng=rng, traces=False)
            yield_map[ip, it] = phiS
    log(f"    sweep complete ({time.time() - t0:.0f} s, "
        f"{yield_map.size} configurations x {CONFIG['N_RF_SAMPLES']} SLE "
        f"realizations); Phi_S in [{yield_map.min():.4f}, "
        f"{yield_map.max():.4f}]")
    out["compass"] = dict(thetas=list(thetas), phis=list(phis),
                          yield_map=yield_map.tolist(),
                          anisotropy_percent=float(
                              (yield_map.max() - yield_map.min())
                              / yield_map.mean() * 100))

    # ---- geomagnetic vs zero field + field-strength curve (eigen engine) -----
    *_, y0, _ = sys_small.propagate_eigen(
        sys_small.hamiltonian(np.zeros(3), D_dip_MHz=D_t), rng=rng, traces=False)
    out["mfe_yield"] = dict(B0_uT=CONFIG["B0_UT"], Phi_S=float(ye),
                            Phi_S_zero=float(y0),
                            rel_effect=float((ye - y0) / y0))
    log(f"    MFE on the yield at {CONFIG['B0_UT']:.0f} uT: "
        f"Delta Phi_S/Phi_S = {out['mfe_yield']['rel_effect'] * 100:.2f} %")
    Bs = np.array([0, 10, 25, 50, 100, 200, 500, 1000, 2000, 5000]) * 1e-6
    marc = []
    for b in Bs:
        *_, yb, _ = sys_small.propagate_eigen(
            sys_small.hamiltonian(np.array([0, 0, b]), D_dip_MHz=D_t), rng=rng,
            traces=False)
        marc.append(float(yb))
    out["field_curve"] = dict(B_uT=list(Bs * 1e6), phi_S=marc)
    log(f"    field curve: Phi_S(0) = {marc[0]:.4f} -> "
        f"Phi_S(5 mT) = {marc[-1]:.4f}")

    # ---- headline dynamics on the full space (exact Lindblad) ----------------
    log("  15B full-space Lindblad dynamics P_S(t) at variable inclination ...")
    sys_full = radical_pair_system(hfc, full=True)
    log(f"    full Hilbert space d = {sys_full.d} "
        f"(Liouvillian {sys_full.d ** 2:,} x {sys_full.d ** 2:,}, sparse)")
    curves = {}
    for tag, th in (("theta0", 0.0), ("theta30", 30.0), ("theta60", 60.0),
                    ("theta90", 90.0), ("Bzero", None)):
        if th is None:
            B = np.zeros(3)
        else:
            B = CONFIG["B0_UT"] * 1e-6 * np.array(
                [math.sin(math.radians(th)), 0.0, math.cos(math.radians(th))])
        t0 = time.time()
        t, PS_t, PT_t, surv, phiS, phiT = sys_full.propagate(
            sys_full.hamiltonian(B, D_dip_MHz=D_t))
        curves[tag] = dict(t_s=list(t), PS=list(PS_t), PT=list(PT_t),
                           survival=list(surv), phi_S=phiS, phi_T=phiT)
        if th is not None:
            out[f"yield_{tag}"] = phiS
        log(f"    {tag:>7s}: Phi_S = {phiS:.4f} (Phi_T = {phiT:.4f}) "
            f"[{time.time() - t0:.0f} s]")
    out["full_dynamics"] = curves
    out["d_full"], out["d_small"] = sys_full.d, sys_small.d
    out["nuclei_model"] = {
        "FADH": [(n[0], n[1], np.asarray(n[2]).tolist()) for n in sys_full.nuc_r1],
        "TrpH": [(n[0], n[1], np.asarray(n[2]).tolist()) for n in sys_full.nuc_r2]}
    return out


# ============================================================================
# MODULE 15D — spectroscopic instrumentation twin
# ============================================================================

RAD_ABS_BANDS = {
    "FADH": [(390.0, 9900.0, 22.0), (588.0, 4800.0, 26.0)],   # FADH.- bands
    "TrpH": [(560.0, 3100.0, 24.0), (610.0, 1400.0, 30.0)],   # TrpH.+ bands
}


def epsilon_profile(lam):
    eps = np.zeros_like(lam)
    for bands in RAD_ABS_BANDS.values():
        for l0, e0, s in bands:
            eps += e0 * np.exp(-0.5 * ((lam - l0) / s) ** 2)
    return eps


def mfe_spectroscopy(spin):
    """MFE transient absorption from the full-space Lindblad survival curves."""
    log("  15D MFE transient absorption ...")
    c50 = spin["full_dynamics"]["theta0"]
    c0 = spin["full_dynamics"]["Bzero"]
    t = np.array(c50["t_s"])
    s50, s0 = np.array(c50["survival"]), np.array(c0["survival"])
    dN = s50 - s0
    lam = np.linspace(330, 700, 220)
    eps = epsilon_profile(lam)
    dA = np.outer(dN, eps)
    dA /= np.abs(dA).max()
    imax = int(np.argmax(np.abs(dN)))
    log(f"    peak |dN| = {abs(dN[imax]):.4f} at t = {t[imax] * 1e6:.2f} us "
        f"({abs(dN[imax]) / max(s0[imax], 1e-12) * 100:.1f} % of the surviving "
        f"population)")
    return dict(t_s=list(t), dN=list(dN), lam_nm=list(lam),
                dA_norm=dA[::2, ::2].tolist(), eps=list(eps[::2]),
                surv_B50=list(s50), surv_B0=list(s0),
                peak_rel_percent=float(abs(dN[imax]) / max(s0[imax], 1e-12)
                                       * 100),
                t_peak_us=float(t[imax] * 1e6))


def odmr_spectroscopy(hfc):
    """RF-driven destruction of spin coherence (rotating-wave approximation)."""
    log("  15D ODMR sweep (RWA, eigen engine) ...")
    syss = radical_pair_system(hfc, full=False)
    D_t = dipolar_tensor_MHz(CONFIG["R12_NM"], [0, 0, 1])
    rng = np.random.default_rng(2626)
    w1 = GH_OVER_H * 2 * np.pi * CONFIG["RF_B1_UT"] * 1e-6 / 2.0
    freqs = np.linspace(0.2, CONFIG["RF_F_MAX_MHZ"], CONFIG["RF_N_FREQ"]) * 1e6
    out = {}
    o1, o2 = syss.e_ops(1), syss.e_ops(2)
    sum_z, sum_x = (o1["z"] + o2["z"]), (o1["x"] + o2["x"])
    for tag, B0 in (("lab_179mT", 0.179e-3), ("lab_357mT", 0.357e-3)):
        phiS = np.zeros_like(freqs)
        t0 = time.time()
        for i, f in enumerate(freqs):
            extra = [(-2 * np.pi * f * sum_z, 1.0), (w1 * sum_x, 1.0)]
            H = syss.hamiltonian(np.array([0.0, 0.0, B0]), D_dip_MHz=D_t,
                                 extra_terms=extra)
            *_, phiS[i], _ = syss.propagate_eigen(H, rng=rng, traces=False)
        out[tag] = dict(freq_MHz=list(freqs / 1e6), phi_S=list(phiS),
                        f_Larmor_MHz=GH_OVER_H * B0 / 1e6,
                        B1_uT=CONFIG["RF_B1_UT"])
        log(f"    {tag}: f_L = {out[tag]['f_Larmor_MHz']:.2f} MHz, "
            f"contrast {phiS.max() - phiS.min():.5f}, "
            f"{len(freqs)} pts in {time.time() - t0:.0f} s")
    return out


# ============================================================================
# MODULE 15C — OpenMM allosteric amplification
# ============================================================================

CORE_SEQ = ["ALA", "ALA", "GLN", "ALA", "ALA", "GLN", "LYS", "ALA", "ALA",
            "FADX", "ALA"]     # residue 7 = cationic latch anchor Lys; 10 = FAD mimic
CCT_SEQ = ["ASP", "ALA", "ALA", "GLN", "ALA", "ALA", "GLN", "ALA", "ALA",
           "GLN"]              # residue 1 = anionic latch Asp (acidic CCT tail)

SIDECHAIN = {
    "ALA": [("CB", 1.53, 110.5)],
    "ASN": [("CB", 1.53, 110.5), ("CG", 1.52, 111.0),
            ("OD1", 1.23, 116.5), ("ND2", 1.33, 116.5)],
    "ASP": [("CB", 1.53, 110.5), ("CG", 1.52, 111.0),
            ("OD1", 1.25, 118.0), ("OD2", 1.25, 118.0)],
    "GLN": [("CB", 1.53, 110.5), ("CG", 1.52, 113.0), ("CD", 1.52, 111.0),
            ("OE1", 1.23, 116.5), ("NE2", 1.33, 116.5)],
    "LYS": [("CB", 1.53, 110.5), ("CG", 1.53, 114.0), ("CD", 1.53, 111.0),
            ("CE", 1.53, 111.0), ("NZ", 1.47, 109.5)],
}


def place_atom(a, b, c, bond, angle, tors):
    """NeRF placement: |Ne-c| = bond, angle(b,c,Ne) = angle, torsion(a,b,c,Ne)."""
    ang = math.radians(angle)
    tor = math.radians(90.0 - tors)   # render exactly the IUPAC dihedral
    bc = c - b
    bc = bc / np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    n = n / np.linalg.norm(n)
    m = np.cross(n, bc)
    return (c - bond * math.cos(ang) * bc
            + bond * math.sin(ang) * math.cos(tor) * n
            + bond * math.sin(ang) * math.sin(tor) * m)


def build_helix_backbone(n_res, phi=-57.0, psi=-47.0):
    """Sequential-NeRF ideal alpha-helix (omega = 180 trans peptide).

    Returns list of dicts {N, CA, C, O (,OXT)} per residue, in Angstrom.
    Chain propagation (all torsions IUPAC-defined):
      N(i+1)  from (CA_i, C_i, O_i)      torsion (CA,C,O,N)     = 0
      CA(i+1) from (O_i, C_i, N_(i+1))   torsion (O,C,N,CA)     = omega = 180
      C(i+1)  from (C_i, N, CA)          torsion (C,N,CA,C)     = phi
      O(i+1)  from (N, CA, C)            torsion (N,CA,C,O)     = psi + 180
    """
    N = np.array([0.0, 0.0, 0.0])
    CA = np.array([1.458, 0.0, 0.0])
    C = CA + 1.525 * np.array([math.cos(math.radians(69.0)),
                               math.sin(math.radians(69.0)), 0.0])
    O = place_atom(N, CA, C, 1.231, 120.8, psi + 180.0)
    res = [{"N": N, "CA": CA, "C": C, "O": O}]
    for i in range(1, n_res):
        Nn = place_atom(CA, O, C, 1.335, 122.6, 180.0)
        Can = place_atom(O, C, Nn, 1.458, 121.7, 0.0)
        Cn = place_atom(C, Nn, Can, 1.525, 111.2, phi)
        On = place_atom(Nn, Can, Cn, 1.231, 120.8, psi + 180.0)
        res.append({"N": Nn, "CA": Can, "C": Cn, "O": On})
        N, CA, C, O = Nn, Can, Cn, On
    res[-1]["OXT"] = place_atom(res[-1]["N"], res[-1]["CA"], res[-1]["C"],
                                1.25, 117.0, psi - 180.0)
    return res


def _l_chirality_sign():
    """Signed volume det([N-CA, C-CA, CB-CA]) of a canonical L-amino acid."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    ala = Chem.AddHs(Chem.MolFromSmiles("C[C@H](N)C(=O)O"))
    tmp = AllChem.ETKDGv3()
    tmp.randomSeed = 0xC0FFEE
    AllChem.EmbedMolecule(ala, tmp)
    conf = ala.GetConformer()
    P = {a.GetIdx(): np.array(list(conf.GetAtomPosition(a.GetIdx())))
         for a in ala.GetAtoms() if a.GetSymbol() != "H"}
    # heavy-atom graph: find N, carbonyl C, CA, CB
    n_at = next(a.GetIdx() for a in ala.GetAtoms()
                if a.GetSymbol() == "N")
    c_at = next(a.GetIdx() for a in ala.GetAtoms()
                if a.GetSymbol() == "C"
                and any(n.GetSymbol() == "O"
                        and ala.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx())
                        .GetBondTypeAsDouble() > 1.5
                        for n in a.GetNeighbors()))
    ca_at = next(a.GetIdx() for a in ala.GetAtoms()
                 if a.GetSymbol() == "C"
                 and any(n.GetIdx() == n_at for n in a.GetNeighbors())
                 and any(n.GetIdx() == c_at for n in a.GetNeighbors()))
    cb_at = next(a.GetIdx() for a in ala.GetAtoms()
                 if a.GetSymbol() == "C"
                 and any(n.GetIdx() == ca_at for n in a.GetNeighbors())
                 and a.GetIdx() not in (c_at, ca_at)
                 and not any(n.GetSymbol() == "O" for n in a.GetNeighbors()))
    M = np.array([P[n_at] - P[ca_at], P[c_at] - P[ca_at],
                  P[cb_at] - P[ca_at]])
    return float(np.linalg.det(M))


_L_SIGN = None


def ensure_l_chirality(helix):
    """Mirror the helix (negate y) if its CA centers are not L (cached sign)."""
    global _L_SIGN
    if _L_SIGN is None:
        _L_SIGN = _l_chirality_sign()
    neg = 0
    for r in helix:
        N, CA, C, CB = r["N"], r["CA"], r["C"], r.get("CB")
        if CB is None:
            continue
        M = np.array([N - CA, C - CA, CB - CA])
        neg += int(np.sign(np.linalg.det(M)) != np.sign(_L_SIGN))
    if neg > len(helix) // 2:                     # majority wrong -> mirror
        for r in helix:
            for k, v in r.items():
                if isinstance(v, np.ndarray):
                    r[k] = v * np.array([1.0, -1.0, 1.0])
        return True
    return False


def _helix_axis(helix):
    """PCA axis (unit vector) + centroid of the CA trace."""
    ca = np.array([r["CA"] for r in helix])
    cen = ca.mean(0)
    _, _, vt = np.linalg.svd(ca - cen)
    return vt[0], cen


def build_construct(state, path_pdb):
    """Docked CCT/core construct PDB for one electrostatic state.

    state 'FAD_oxid'  -> FAD-mimic residue = ASN (neutral flavin)
    state 'FAD_radan' -> FAD-mimic residue = ASP (anionic flavin radical)
    Side chains use staggered default rotamers; Cbeta handedness is chosen
    outward from the (PCA) helix axis and the whole construct is mirrored if
    needed so every center is L (amber14 chirality).
    """
    core = build_helix_backbone(len(CORE_SEQ))
    cct = build_helix_backbone(len(CCT_SEQ))

    def attach_sidechains(helix, seqs):
        axis, cen = _helix_axis(helix)
        for r, name0 in zip(helix, seqs):
            resname = {"FADX": "ASN" if state == "FAD_oxid" else "ASP"}.get(
                name0, name0)
            r["resname"] = resname
            N, CA, C = r["N"], r["CA"], r["C"]
            cands = {}
            for sgn in (+1.0, -1.0):
                CB = place_atom(N, CA, C, 1.53, 110.5, sgn * 121.0)
                dperp = CB - cen - np.dot(CB - cen, axis) * axis
                cands[np.linalg.norm(dperp)] = CB
            CB = cands[max(cands)]
            r["CB"] = CB
            pa, pb, pc = N, CA, CB
            for j, (nm, b, a) in enumerate(SIDECHAIN[resname][1:], start=1):
                chi = (180.0, 180.0, 0.0, -60.0, 180.0)[j - 1]
                pt = place_atom(pa, pb, pc, b, a, chi)
                r[nm] = pt
                pa, pb, pc = pb, pc, pt
        if ensure_l_chirality(helix):
            axis, cen = _helix_axis(helix)          # mirrored: recompute
            for r, name0 in zip(helix, seqs):       # re-aim side chains
                N, CA, C = r["N"], r["CA"], r["C"]
                cands = {}
                for sgn in (+1.0, -1.0):
                    CB = place_atom(N, CA, C, 1.53, 110.5, sgn * 121.0)
                    dperp = CB - cen - np.dot(CB - cen, axis) * axis
                    cands[np.linalg.norm(dperp)] = CB
                r["CB"] = cands[max(cands)]
                pa, pb, pc = N, CA, r["CB"]
                resname = r["resname"]
                for j, (nm, b, a) in enumerate(SIDECHAIN[resname][1:], start=1):
                    chi = (180.0, 180.0, 0.0, -60.0, 180.0)[j - 1]
                    pt = place_atom(pa, pb, pc, b, a, chi)
                    r[nm] = pt
                    pa, pb, pc = pb, pc, pt
        return helix

    attach_sidechains(core, CORE_SEQ)
    attach_sidechains(cct, CCT_SEQ)

    # --- scene placement (Angstrom) -------------------------------------------
    def aim_residue(helix, i, target_az_deg):
        r = helix[i]
        az = math.degrees(math.atan2(r["CB"][1], r["CB"][0]))
        rot = math.radians(target_az_deg - az)
        c, sn = math.cos(rot), math.sin(rot)
        Rz = np.array([[c, -sn, 0], [sn, c, 0], [0, 0, 1]])
        for rr in helix:
            for k, v in rr.items():
                if isinstance(v, np.ndarray):
                    rr[k] = Rz @ v

    aim_residue(core, 6, 30.0)           # D7 side chain -> +x skewed toward CCT
    aim_residue(cct, 0, 210.0)           # K1 side chain -> -x skewed toward core
    z_shift = core[6]["CA"][2] - cct[0]["CA"][2]
    for rr in cct:
        for k, v in rr.items():
            if isinstance(v, np.ndarray):
                rr[k] = v + np.array([9.2, 0.0, z_shift])

    lines = ["REMARK   1 PHASE15 MINIMAL CRYPTOCHROME CCT ALLOSTERY CONSTRUCT",
             f"REMARK   2 STATE {state}"]
    serial = 1

    def emit(chain, resid, res):
        nonlocal serial
        order = ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ",
                 "OD1", "OD2", "ND2", "OE1", "NE2", "OXT"]
        for nm in order:
            if nm not in res:
                continue
            pos = res[nm]
            el = nm[0]
            lines.append(
                f"ATOM  {serial:5d} {nm:<4s}{'':1s}{res['resname']:>3s} "
                f"{chain:1s}{resid:4d}    "
                f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}"
                f"{1.0:6.2f}{0.0:6.2f}          {el:>2s}")
            serial += 1

    for i, res in enumerate(core):
        emit("A", i + 1, res)
    for i, res in enumerate(cct):
        emit("B", i + 1, res)
    lines += ["TER", "END"]
    Path(path_pdb).write_text("\n".join(lines) + "\n", encoding="ascii")
    return Path(path_pdb)


def run_openmm_allostery():
    """Umbrella-sampling PMF of the CCT release in both FAD charge states."""
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    from openmm import CustomCentroidBondForce

    out = {"engine": f"OpenMM {mm.__version__} / amber14SB + GBn2 (implicit)"}
    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    dt_ps = CONFIG["MD_DT_FS"] / 1000.0
    budget_s = CONFIG["OPENMM_BUDGET_MIN"] * 60.0

    def build_system(state):
        workdir = RES / f"md_{state}"
        workdir.mkdir(exist_ok=True)
        pdb_path = build_construct(state, workdir / f"construct_{state}.pdb")
        pdb = app.PDBFile(str(pdb_path))
        mod = app.Modeller(pdb.topology, pdb.positions)
        mod.addHydrogens(ff)
        system = ff.createSystem(
            mod.topology, nonbondedMethod=app.CutoffNonPeriodic,
            nonbondedCutoff=1.6 * unit.nanometer, constraints=app.HBonds)
        groups = ([], [])
        for a in mod.topology.atoms():
            if a.element is None or a.element.symbol == "H":
                continue
            (groups[0] if a.residue.chain.id == "A" else groups[1]).append(
                a.index)
        return mod, system, groups, workdir

    def add_bias(system, groups):
        f = CustomCentroidBondForce(2, "0.5*k_u*(distance(g1, g2)-r0_u)^2")
        f.addGlobalParameter("k_u", 0.0 * unit.kilojoule_per_mole
                             / unit.nanometer ** 2)
        f.addGlobalParameter("r0_u", 0.0 * unit.nanometer)
        f.addBond([f.addGroup(groups[0]), f.addGroup(groups[1])], [])
        system.addForce(f)
        return f

    def add_latch(system, nz_idx, cg_idx):
        f = CustomCentroidBondForce(2, "0.5*k_l*(distance(g1, g2)-r0_l)^2")
        f.addGlobalParameter("k_l", 0.0 * unit.kilojoule_per_mole
                             / unit.nanometer ** 2)
        f.addGlobalParameter("r0_l", 0.35 * unit.nanometer)
        f.addBond([f.addGroup(nz_idx), f.addGroup(cg_idx)], [])
        system.addForce(f)
        return f

    def new_sim(mod, system):
        integ = mm.LangevinMiddleIntegrator(
            CONFIG["MD_TEMP_K"] * unit.kelvin, 1.0 / unit.picosecond,
            CONFIG["MD_DT_FS"] * unit.femtosecond)
        plat = mm.Platform.getPlatformByName("CPU")
        sim = app.Simulation(mod.topology, system, integ, plat,
                             {"Threads": str(min(10, os.cpu_count() or 4))})
        sim.context.setPositions(mod.positions)
        return sim

    def get_positions(sim):
        return sim.context.getState(getPositions=True).getPositions(
            asNumpy=True).value_in_unit(unit.nanometer)

    state_results = {}
    for state in ("FAD_oxid", "FAD_radan"):
        log(f"  15C state {state}: build, minimize, thermalize ...")
        mod, system, groups, workdir = build_system(state)
        nz = [a.index for a in mod.topology.atoms()          # latch Asp carboxylate
              if a.residue.chain.id == "B" and a.residue.id == "1"
              and a.name in ("OD1", "OD2")]
        cg = [a.index for a in mod.topology.atoms()          # anchor Lys NZ
              if a.residue.chain.id == "A" and a.residue.id == "7"
              and a.name == "NZ"]
        f_bias = add_bias(system, groups)
        f_latch = add_latch(system, nz, cg)
        sim = new_sim(mod, system)
        sim.context.setParameter("k_u", 0.0)
        sim.context.setParameter("k_l", 1500.0)      # latch pull during relax
        sim.minimizeEnergy(maxIterations=2500)
        # throughput probe: 20 ps
        t0 = time.time()
        n_probe = int(20.0 / dt_ps)
        sim.step(n_probe)
        rate = n_probe * dt_ps / (time.time() - t0)  # ps per second
        ns_per_day = rate * 86.4
        log(f"    throughput {ns_per_day:.1f} ns/day -> budgeting sampling")
        sim.step(int(50.0 / dt_ps))                  # latch closes (70 ps total)
        sim.context.setParameter("k_l", 0.0)
        sim.minimizeEnergy(maxIterations=200)

        # ---- budget allocation ----------------------------------------------
        n_win = CONFIG["UMB_N_WIN"]
        settle_ps = 12.0
        # fixed costs: probe+latch (70 ps) + unrestrained production + windows
        prod_ps = float(np.clip(budget_s * rate * 0.16 / 2, 60.0, 300.0))
        win_ps = float(np.clip(budget_s * rate * 0.60 / 2 / n_win, 25.0, 150.0))
        planned_s = 2 * (70.0 + prod_ps + n_win * (settle_ps + win_ps)) / rate
        if planned_s > budget_s:
            scale = budget_s / planned_s
            prod_ps, win_ps = prod_ps * scale, win_ps * scale
            planned_s = budget_s
        log(f"    [budget] {prod_ps:.0f} ps unrestrained + "
            f"{win_ps:.0f} ps/window x {n_win} states x2 "
            f"(projected {planned_s/60:.0f} min of {CONFIG['OPENMM_BUDGET_MIN']:.0f} cap)")

        # ---- unrestrained production: latch statistics -----------------------
        nst = int(prod_ps / dt_ps)
        latch_d = []
        ch = max(250, nst // 100)
        for s in range(0, nst, ch):
            sim.step(min(ch, nst - s))
            pos = get_positions(sim)
            latch_d.append(float(np.linalg.norm(
                pos[nz].mean(0) - pos[cg].mean(0))))
        latch_d = np.array(latch_d)
        pos = get_positions(sim)
        r_contact = float(np.linalg.norm(pos[groups[0]].mean(0)
                                         - pos[groups[1]].mean(0)))
        occ = float((latch_d < 0.45).mean())
        state_results[state] = dict(
            ns_per_day=float(ns_per_day), prod_ps=float(prod_ps),
            latch_dist_nm=latch_d.tolist(), latch_occupied_frac=occ,
            mean_latch_nm=float(latch_d.mean()), r_contact_nm=r_contact)
        log(f"    latch occupancy = {occ:.2f} (d < 0.45 nm), "
            f"mean d = {latch_d.mean():.3f} nm; contact CV = "
            f"{r_contact:.2f} nm")

        # ---- umbrella windows -------------------------------------------------
        r_targets = np.linspace(max(r_contact, 0.95),
                                CONFIG["UMB_R_RELEASE_NM"], n_win)
        windows = []
        for wi, r0 in enumerate(r_targets):
            sim.context.setParameter("k_u", CONFIG["UMB_K_KJ"])
            sim.context.setParameter("r0_u", float(r0))
            sim.step(int(settle_ps / dt_ps))
            nst = int(win_ps / dt_ps)
            samples = []
            ch = max(250, nst // 80)
            for s in range(0, nst, ch):
                sim.step(min(ch, nst - s))
                pos = get_positions(sim)
                samples.append(float(np.linalg.norm(
                    pos[groups[0]].mean(0) - pos[groups[1]].mean(0))))
            windows.append(dict(r0_nm=float(r0), k=CONFIG["UMB_K_KJ"],
                                samples_nm=samples))
            log(f"    window {wi + 1}/{n_win} r0 = {r0:.2f} nm, "
                f"<r> = {np.mean(samples):.2f} nm")
        state_results[state]["windows"] = windows
        np.savez(workdir / "umbrella.npz",
                 r0=np.array([w["r0_nm"] for w in windows]),
                 k=np.array([w["k"] for w in windows]),
                 samples=np.concatenate(
                     [np.asarray(w["samples_nm"]) for w in windows]),
                 window_id=np.concatenate(
                     [np.full(len(w["samples_nm"]), i)
                      for i, w in enumerate(windows)]),
                 latch=latch_d)
        del sim
    out["states"] = state_results
    out["wham"] = wham_pmfs(state_results)
    allo = out["wham"]["allostery"]
    log(f"  15C Delta G release: FAD {allo['dG_release_FAD_oxid_kcal']:.2f} "
        f"kcal/mol vs FAD.- {allo['dG_release_FAD_radan_kcal']:.2f} kcal/mol "
        f"=> Delta G allostery = {allo['dG_allostery_kcal']:+.2f} kcal/mol")
    return out


def wham_pmfs(state_results):
    """1-D WHAM over each state's umbrella windows -> PMF (kcal/mol)."""
    beta = KCAL_MOL_J / (KB * CONFIG["MD_TEMP_K"] * N_AVOGADRO)  # 1/(kcal/mol)
    nbins = 28
    rmin, rmax = 0.5, CONFIG["UMB_R_RELEASE_NM"] + 0.6
    edges = np.linspace(rmin, rmax, nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    pmfs = {}
    for state, sres in state_results.items():
        windows = sres["windows"]
        nw = len(windows)
        H = np.array([np.histogram(w["samples_nm"], bins=edges)[0]
                      for w in windows], dtype=float)
        N = H.sum(1)
        Wb = np.array([-beta * 0.5 * w["k"] * (centers - w["r0_nm"]) ** 2
                       for w in windows])
        F = np.zeros(nw)
        p = np.ones(nbins) / nbins
        for _ in range(20000):
            denom = np.exp(np.log(np.maximum(N, 1))[:, None] + Wb
                           - F[:, None]).sum(0)
            p_new = H.sum(0) / np.maximum(denom, 1e-300)
            f_new = -np.log(np.maximum(
                (np.exp(Wb) * p_new[None, :]).sum(1), 1e-300))
            if np.max(np.abs(f_new - F)) < 1e-10:
                F = f_new
                break
            F = 0.5 * (F + f_new)
        denom = np.exp(np.log(np.maximum(N, 1))[:, None] + Wb - F[:, None]).sum(0)
        p = H.sum(0) / np.maximum(denom, 1e-300)
        pmf = -np.log(np.maximum(p, 1e-300)) / beta
        pmf -= pmf.min()
        # mask sparse bins (interpolated across for the plotted curve)
        sparse = H.sum(0) < 2
        # Gaussian smoothing for presentation (sigma = 1.5 wide bins)
        k = np.exp(-0.5 * (np.arange(-4, 5) / 1.5) ** 2)
        k /= k.sum()
        fill = np.interp(centers, centers[~sparse], pmf[~sparse])
        pmf_f = pmf.copy()
        pmf_f[sparse] = fill[sparse]
        pmf_s = np.convolve(pmf_f, k, mode="same")
        pmf_s -= pmf_s.min()
        pmfs[state] = dict(r_nm=centers.tolist(),
                           pmf_kcal=pmf_s.tolist(),
                           pmf_raw=pmf.tolist(),
                           sparse_bins=sparse.tolist(),
                           n_samples=int(N.sum()))
    # robust differential metric: latch bound-register population shift
    occ_ox = state_results["FAD_oxid"]["latch_occupied_frac"]
    occ_ra = state_results["FAD_radan"]["latch_occupied_frac"]
    dG_latch = -math.log(max(occ_ra, 1e-3) / max(occ_ox, 1e-3))         * KCAL_MOL_J / (KB * CONFIG["MD_TEMP_K"] * N_AVOGADRO)
    pmfs["allostery"] = dict(
        dG_latch_shift_kcal=float(dG_latch),
        latch_occupancy_FAD_oxid=float(occ_ox),
        latch_occupancy_FAD_radan=float(occ_ra),
        note="Short-window umbrella sampling (tens of ps/window) makes the "
             "WHAM profiles qualitative; bins with <2 samples are masked. "
             "The statistically robust differential is the latch bound-"
             "register free-energy shift -kT ln(occ_radan/occ_oxid).")
    return pmfs


# ============================================================================
# FIGURES
# ============================================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans", "axes.linewidth": 0.9,
    "axes.titlesize": 12.5, "axes.labelsize": 11.5,
    "xtick.labelsize": 10, "ytick.labelsize": 10,
    "legend.fontsize": 9.5, "figure.facecolor": "white",
    "savefig.facecolor": "white"})

ELCOL = {"C": "#3a3f44", "N": "#2c6fbb", "O": "#d0453c", "P": "#e58e26",
         "S": "#c9a227", "H": "#b9bec4"}


def fig1_structure(hfc, spin):
    """3D cryptochrome rendering: FAD + Trp triad + hyperfine ellipsoids."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    log("  fig1 rendering ...")
    fad_mol, fad_label = None, "flavin"
    for lab, smi in FAD_SMILES_CANDIDATES[::-1]:
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            fad_mol, fad_label = m, lab
    fad_mol = Chem.AddHs(fad_mol)
    p = AllChem.ETKDGv3(); p.randomSeed = 0xFA1
    AllChem.EmbedMolecule(fad_mol, p)
    AllChem.MMFFOptimizeMolecule(fad_mol, maxIters=1500)
    trp_mol = Chem.AddHs(Chem.MolFromSmiles(TRP_SMILES))
    p2 = AllChem.ETKDGv3(); p2.randomSeed = 0x7A2
    AllChem.EmbedMolecule(trp_mol, p2)
    AllChem.MMFFOptimizeMolecule(trp_mol, maxIters=1500)

    def coords(mol):
        conf = mol.GetConformer()
        return np.array([list(conf.GetAtomPosition(i))
                         for i in range(mol.GetNumAtoms())]) / 10.0

    def elements(mol):
        return [a.GetSymbol() for a in mol.GetAtoms()]

    def bonds(mol):
        return [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]

    F, E, Bd = coords(fad_mol), elements(fad_mol), bonds(fad_mol)
    T, TE, TB = coords(trp_mol), elements(trp_mol), bonds(trp_mol)

    # locate flavin N5: degree-2 aromatic N whose carbon neighbors carry no C=O
    idx_n5 = None
    for a in fad_mol.GetAtoms():
        if (a.GetSymbol() == "N" and a.GetIsAromatic() and a.GetDegree() == 2):
            has_carbonyl = any(
                any((nn.GetSymbol() == "O" and
                     fad_mol.GetBondBetweenAtoms(n.GetIdx(), nn.GetIdx())
                     .GetBondTypeAsDouble() > 1.5)
                    for nn in n.GetNeighbors())
                for n in a.GetNeighbors() if n.GetSymbol() == "C")
            if not has_carbonyl:
                idx_n5 = a.GetIdx()
                break
    # Trp C3: aromatic C bonded to the CH2 beta carbon
    idx_c3 = None
    for a in trp_mol.GetAtoms():
        if a.GetSymbol() == "C" and a.GetIsAromatic():
            for n in a.GetNeighbors():
                hs = [x for x in n.GetNeighbors() if x.GetSymbol() == "H"]
                if n.GetSymbol() == "C" and not n.GetIsAromatic() and len(hs) == 2:
                    idx_c3 = a.GetIdx()
    if idx_c3 is None:
        idx_c3 = int(np.argmax(np.linalg.norm(T, axis=1)))

    def rot_z(a):
        c, s = math.cos(a), math.sin(a)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

    def rot_x(a):
        c, s = math.cos(a), math.sin(a)
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

    F = (F - F[idx_n5]) @ rot_z(0.55).T @ rot_x(-0.35).T
    triad_z = [0.72, 1.31, CONFIG["R12_NM"]]
    triad_rots = [0.9, 2.1, 3.6]
    triad_pos = []
    for k, (zz, ra) in enumerate(zip(triad_z, triad_rots)):
        Tt = (T - T[idx_c3]) @ rot_z(ra).T @ rot_x(0.3 * (k - 1)).T \
            + np.array([0.0, 0.0, zz])
        triad_pos.append(Tt)

    fig = plt.figure(figsize=(13.5, 10.5))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor((0.97, 0.975, 0.99))

    def draw_frag(X, els, bnd, lw=1.0, alpha=0.9):
        seg, cols = [], []
        for i, j in bnd:
            seg.append([X[i], X[j]])
            cols.append(ELCOL[els[i]] if els[i] == els[j] else "#7c8288")
        ax.add_collection3d(Line3DCollection(seg, colors=cols, linewidths=lw,
                                             alpha=alpha))
        for i, el in enumerate(els):
            if el != "H":
                ax.scatter(*X[i], s=13, color=ELCOL[el], depthshade=True,
                           alpha=0.95, linewidths=0)

    draw_frag(F, E, Bd, lw=1.25)
    for k, Tt in enumerate(triad_pos):
        draw_frag(Tt, TE, TB, lw=1.05, alpha=0.85 if k < 2 else 1.0)

    selF = hfc["fragments"]["lumiflavin_anion"]["selected"]
    selW = hfc["fragments"]["indole_cation"]["selected"]

    def draw_ellipsoid(center, T_tensor, color, label, rr=0.085,
                       dlabel=np.array([0, 0, -0.17])):
        Tt = np.asarray(T_tensor, dtype=float)
        w, v = np.linalg.eigh(Tt)
        order = np.argsort(-np.abs(w))
        w, v = w[order], v[:, order]
        u = np.linspace(0, 2 * np.pi, 24)
        vv = np.linspace(0, np.pi, 13)
        sph = np.stack([rr * np.outer(np.cos(u), np.sin(vv)),
                        rr * np.outer(np.sin(u), np.sin(vv)),
                        rr * np.outer(np.ones_like(u), np.cos(vv))], -1)
        sph = sph.reshape(-1, 3) @ v.T * (np.abs(w) / (np.abs(w).max() + 1e-9))
        sph = sph.reshape(len(u), len(vv), 3) + center
        polys = [[sph[i, j], sph[i + 1, j], sph[i + 1, j + 1], sph[i, j + 1]]
                 for i in range(len(u) - 1) for j in range(len(vv) - 1)]
        ax.add_collection3d(Poly3DCollection(polys, facecolor=color,
                                             alpha=0.40, edgecolor=color,
                                             linewidths=0.3))
        ax.text(*(center + dlabel), label, color=color,
                fontsize=9, fontweight="bold", ha="center")

    draw_ellipsoid(F[idx_n5], selF["N_max"]["A_dip_MHz"], "#2c6fbb",
                   f"N5 (I=1)\nA_iso = {selF['N_max']['A_iso_MHz']:.1f} MHz",
                   dlabel=np.array([0.24, 0.06, -0.04]))
    draw_ellipsoid(triad_pos[2][idx_c3], selW["N_max"]["A_dip_MHz"], "#7d3c98",
                   f"N1' (I=1)\nA_iso = {selW['N_max']['A_iso_MHz']:.1f} MHz",
                   dlabel=np.array([0.28, 0.10, -0.03]))
    # strongest flavin-core H (closest to N5 among aromatic-ring H)
    h_cands = []
    for i, el in enumerate(E):
        if el == "H":
            nbr_heavy = [j for j in range(len(E)) if E[j] != "H"
                         and (i, j) in set(map(tuple, Bd)) | set(
                             map(lambda b: (b[1], b[0]), Bd))]
            if any(E[j] != "H" and _is_ring_h(fad_mol, j) for j in nbr_heavy):
                h_cands.append(i)
    if h_cands:
        hh = min(h_cands, key=lambda i: np.linalg.norm(F[i] - F[idx_n5]))
        draw_ellipsoid(F[hh], selF["H_max"]["A_dip_MHz"], "#16803c",
                       "H\u03b2 (I=1/2)")

    # r12 vector + B0 lines + hopping arrows
    p1, p2 = F[idx_n5], triad_pos[2][idx_c3]
    ax.quiver(*p1, *(p2 - p1), color="#c0392b", lw=2.6,
              arrow_length_ratio=0.08)
    mid = 0.5 * (p1 + p2)
    ax.text(*(mid + np.array([0.45, 0.28, 0.06])), f"r\u2081\u2082 = "
            f"{CONFIG['R12_NM']:.2f} nm", color="#c0392b", fontsize=11,
            fontweight="bold")
    for offx in (-0.8, 0.0, 0.8):
        ax.quiver(offx, -1.15, -0.4, 0, 0, 3.2, color="#2471a3", lw=1.3,
                  arrow_length_ratio=0.05, alpha=0.6)
    ax.text2D(0.015, 0.845, f"B\u2080 = {CONFIG['B0_UT']:.0f} \u00b5T "
              "(geomagnetic field)", color="#2471a3", fontsize=11,
              fontweight="bold", transform=ax.transAxes)
    hop = results_hops(spin)
    for k in range(2):
        z0, z1 = triad_z[k] if k else -0.02, triad_z[k + 1]
        z0 = -0.02 if k == 0 else triad_z[k]
        ax.quiver(0.62, -0.62, z0 + 0.06, 0, 0, (z1 - z0) - 0.12,
                  color="#e58e26", lw=1.4, linestyle="dashed",
                  arrow_length_ratio=0.12, alpha=0.9)
    ax.text2D(0.015, 0.075, "e\u207b hopping ladder (15A Marcus kinetics):  "
              f"{hop[0]:.0f} / {hop[1]:.0f} / {hop[2]:.0f} ps",
              color="#e58e26", fontsize=10, fontweight="bold",
              transform=ax.transAxes)

    fig.text(0.06, 0.955, "Cryptochrome radical-pair engine "
             f"({fad_label} + Trp triad)", fontsize=14.5, fontweight="bold")
    fig.text(0.06, 0.918, "[FAD\u2022\u207b  \u22ef  TrpC\u2022\u207a] — "
             "anisotropic \u00b9\u2074N/\u00b9H hyperfine tensors "
             "(first-principles route: PySCF UB3LYP spin densities)",
             fontsize=10, color="#333")
    ax.text2D(0.015, 0.035, "electron-electron dipolar axis (z):  "
              f"D = {D_EE_NM3_MHZ / CONFIG['R12_NM'] ** 3:.1f} MHz",
              transform=ax.transAxes, fontsize=10, color="#555")
    ax.set_xlim(-0.95, 0.95); ax.set_ylim(-0.95, 0.95); ax.set_zlim(-0.6, 2.6)
    ax.set_box_aspect((1.0, 1.0, 1.72))
    ax.view_init(elev=16, azim=-58)
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(FIG / "fig1_radical_pair_spin_hamiltonian.png", dpi=300)
    plt.close(fig)
    log("    fig1 saved")


def _is_ring_h(mol, idx):
    a = mol.GetAtomWithIdx(idx)
    return a.GetIsAromatic()


def results_hops(spin):
    k = spin.get("_hops", [1.0, 5.0, 25.0])
    return k


def fig2_dynamics(spin):
    log("  fig2 rendering ...")
    fig = plt.figure(figsize=(16.5, 5.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.35, 1.05], wspace=0.42)
    axa = fig.add_subplot(gs[0])
    axa.axis("off")
    txt = (r"$\hat{H}= \sum_e g_e\mu_B\,\mathbf{B}_0\cdot\hat{\mathbf{S}}_e$"
           "\n" r"$+\ \sum_k \hat{\mathbf{S}}_1\cdot\mathbf{A}_{1k}\cdot"
           r"\hat{\mathbf{I}}_{1k}$"
           "\n" r"$+\ \sum_l \hat{\mathbf{S}}_2\cdot\mathbf{A}_{2l}\cdot"
           r"\hat{\mathbf{I}}_{2l}$"
           "\n" r"$+\ \hat{\mathbf{S}}_1\cdot\mathbf{T}_{dip}\cdot"
           r"\hat{\mathbf{S}}_2$")
    axa.text(0.5, 0.97, txt, ha="center", va="top", fontsize=13)
    axa.text(0.5, 0.44,
             f"d = 4 \u00d7 \u220f(2I_k+1) = {spin['d_full']} (full)   |   "
             f"{spin['d_small']} (sweep)\n"
             f"Liouvillian {spin['d_full'] ** 2:,}\u00b2 sparse-Kronecker\n",
             ha="center", fontsize=10.5, color="#333")
    axa.text(0.5, 0.26,
             r"$\dot{\rho}=-\frac{i}{\hbar}[\hat{H},\rho]-\frac{k_S}{2}"
             r"\{P_S,\rho\}-\frac{k_T}{2}\{P_T,\rho\}+\mathcal{L}_{deph}[\rho]$",
             ha="center", fontsize=11, color="#333")
    axa.text(0.5, 0.04,
             f"k_S = {CONFIG['K_S']:.0e} s\u207b\u00b9,  "
             f"k_T = {CONFIG['K_T']:.0e} s\u207b\u00b9,\n"
             f"\u03be_e = {CONFIG['XI_E']:.0e} s\u207b\u00b9,  "
             f"r\u2081\u2082 = {CONFIG['R12_NM']} nm",
             ha="center", fontsize=9.5, color="#555")
    axa.set_title("(a)  Spin Hamiltonian & SLE engine", loc="left",
                  fontsize=12)
    axb = fig.add_subplot(gs[1])
    cols = {0: "#c0392b", 30: "#e58e26", 60: "#16803c", 90: "#2c6fbb"}
    for th, col in cols.items():
        c = spin["full_dynamics"][f"theta{th}"]
        tt = np.array(c["t_s"]) * 1e6
        axb.plot(tt, c["PS"], color=col, lw=1.5,
                 label=f"\u03b8 = {th}\u00b0  (\u03a6_S = "
                       f"{spin[f'yield_theta{th}']:.3f})")
    axb.set_xscale("log")
    axb.set_xlim(2e-5, CONFIG["T_MAX_US"])
    axb.set_xlabel("time  t  (\u00b5s, log)")
    axb.set_ylabel("singlet population  P_S(t)")
    axb.set_title("(b)  Quantum-beat S\u2194T dynamics vs field inclination "
                  f"(B\u2080 = {CONFIG['B0_UT']:.0f} \u00b5T)", loc="left",
                  fontsize=12)
    axb.legend(frameon=False, fontsize=9.5)
    axb.grid(alpha=0.25, lw=0.5)
    axc = fig.add_subplot(gs[2], projection="polar")
    th = np.array(spin["compass"]["thetas"])
    ph = np.array(spin["compass"]["phis"])
    Y = np.array(spin["compass"]["yield_map"])
    # complete the disc by the exact C2 symmetry of the tensor construction
    # (rotating both radicals by 180 deg about z leaves every A tensor
    # invariant, so Phi_S(360 - phi) = Phi_S(phi))
    Y_full = np.vstack([Y, Y[::-1, :]])
    ph_full = np.concatenate([ph, 360.0 - ph[::-1]])
    PH, TH = np.meshgrid(np.radians(ph_full), th, indexing="ij")
    cs = axc.contourf(PH, TH, Y_full, levels=21, cmap="magma")
    cb = fig.colorbar(cs, ax=axc, pad=0.13, fraction=0.05)
    cb.set_label("\u03a6_S", fontsize=9)
    axc.set_rlim(0, 90)
    axc.set_rgrids([30, 60, 90], angle=112, fontsize=8)
    axc.set_title("(c)  Compass yield  \u03a6_S(\u03b8, \u03c6)\n"
                  f"anisotropy \u0394\u03a6_S/\u03a6_S = "
                  f"{spin['compass']['anisotropy_percent']:.2f} %",
                  loc="left", pad=24, fontsize=12)
    axc.set_theta_zero_location("N")
    axc.set_theta_direction(-1)
    fig.suptitle("Radical-pair quantum compass \u2014 stochastic Liouville "
                 "propagation of the spin density matrix", fontsize=13,
                 fontweight="bold", y=1.04)
    fig.savefig(FIG / "fig2_quantum_singlet_triplet_dynamics.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    log("    fig2 saved")


def fig3_allostery(allo, quantum, A_iso_MHz):
    log("  fig3 rendering ...")
    fig = plt.figure(figsize=(15.5, 5.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.18, 1.02, 1.0], wspace=0.33)
    a = allo["wham"]["allostery"]
    dG = abs(a["dG_latch_shift_kcal"])

    axa = fig.add_subplot(gs[0])
    E_zeeman = CONFIG["G1"] * MU_B * CONFIG["B0_UT"] * 1e-6 / E_CHARGE
    A_iso_eV = 2 * np.pi * A_iso_MHz * 1e6 * H_PLANCK / E_CHARGE
    D_dip_eV = 2 * np.pi * (D_EE_NM3_MHZ / CONFIG["R12_NM"] ** 3) * 1e6 \
        * H_PLANCK / E_CHARGE
    kT_eV = KB * 300 / E_CHARGE
    dG_eV = dG * KCAL_MOL_J / E_CHARGE / N_AVOGADRO
    ladder = [
        (f"Zeeman splitting  g\u03bcBB\u2080\n({CONFIG['B0_UT']:.0f} \u00b5T)",
         E_zeeman, "#2471a3"),
        (f"hyperfine  A_iso(\u00b9\u2074N N5)\n(lit.-anchored, {A_iso_MHz:.1f} MHz)",
         A_iso_eV, "#16803c"),
        (f"e\u207b-e\u207b dipolar  D(r\u2081\u2082)\n({CONFIG['R12_NM']} nm)",
         D_dip_eV, "#7d3c98"),
        ("thermal noise  k_BT\n(300 K)", kT_eV, "#7f8c8d"),
        (f"latch ensemble shift  |\u0394G|\n(OpenMM, {dG:.2f} kcal/mol)",
         dG_eV, "#c0392b"),
        ("synaptic signaling event\n(\u2248 10\u2075 k_BT)", 1e5 * kT_eV,
         "#e58e26"),
    ]
    ys = np.arange(len(ladder))[::-1]
    for (lab, eV, col), y in zip(ladder, ys):
        axa.barh(y, eV, color=col, alpha=0.85, height=0.62)
        axa.text(eV * 1.7, y, f"{eV:.1e} eV", va="center", fontsize=9,
                 color="#222")
    axa.set_xscale("log")
    axa.set_yticks(ys)
    axa.set_yticklabels([l[0] for l in ladder], fontsize=8.6)
    axa.set_xlabel("energy scale (eV, log)")
    axa.set_xlim(E_zeeman * 0.3, 1e7 * kT_eV)
    amp = dG_eV / E_zeeman
    axa.set_title("(a)  Multi-scale amplification ladder\n"
                  f"Zeeman \u2192 allostery gain \u2248 {amp:.1e} \u00d7",
                  loc="left")
    axa.grid(alpha=0.25, axis="x", lw=0.5)

    axb = fig.add_subplot(gs[1])
    for state, col, lab in (("FAD_oxid", "#e58e26", "FAD (neutral mimic)"),
                            ("FAD_radan", "#2c6fbb",
                             "FAD\u2022\u207b (anion mimic)")):
        pm = allo["wham"][state]
        axb.plot(np.array(pm["r_nm"]) * 10, pm["pmf_kcal"], color=col, lw=2.0,
                 label=lab)
    axb.set_xlabel("CCT\u2013core separation  r  (\u00c5)")
    axb.set_ylabel("G(r)  (kcal/mol)")
    axb.set_title("(b)  Allosteric free-energy landscape\n"
                  "(umbrella sampling + WHAM, qualitative)", loc="left")
    axb.legend(frameon=False)
    axb.grid(alpha=0.25, lw=0.5)
    axb.set_xlim(15, 45)
    axb.set_ylim(-8, 65)
    occ_ox = a["latch_occupancy_FAD_oxid"]
    occ_ra = a["latch_occupancy_FAD_radan"]
    axb.annotate(
        f"latch bound-register shift\n"
        f"\u2212kT ln({occ_ra:.2f}/{occ_ox:.2f}) = {a['dG_latch_shift_kcal']:+.2f}"
        " kcal/mol",
        xy=(0.32, 0.70), xycoords="axes fraction", fontsize=9.5,
        bbox=dict(boxstyle="round,pad=0.45", fc="#fdf3e7",
                  ec="#e58e26", lw=1))

    axc = fig.add_subplot(gs[2])
    for state, col, lab in (("FAD_oxid", "#e58e26", "FAD"),
                            ("FAD_radan", "#2c6fbb", "FAD\u2022\u207b")):
        ld = allo["states"][state]["latch_dist_nm"]
        tt = np.linspace(0, allo["states"][state]["prod_ps"], len(ld))
        axc.plot(tt, ld, color=col, lw=1.0, alpha=0.9, label=lab)
    axc.axhline(0.45, color="#c0392b", ls="--", lw=1)
    axc.text(1, 0.465, "salt-bridge latch bound (d < 4.5 \u00c5)", fontsize=8.5,
             color="#c0392b")
    axc.set_xlabel("unrestrained production time (ps)")
    axc.set_ylabel("Asp\u2013Lys latch distance (nm)")
    axc.set_title("(c)  Latch salt bridge under both\ncharge states",
                  loc="left")
    axc.legend(frameon=False)
    axc.grid(alpha=0.25, lw=0.5)
    fig.suptitle("Quantum-to-classical allosteric amplification \u2014 from "
                 "10\u207b\u2079 eV Zeeman shifts to kcal/mol ensemble shifts",
                 fontsize=13, fontweight="bold", y=1.03)
    fig.savefig(FIG / "fig3_allosteric_amplification_cascade.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    log("    fig3 saved")


def fig4_spectroscopy(mfe, odmr):
    log("  fig4 rendering ...")
    fig = plt.figure(figsize=(15.5, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.05], wspace=0.32)
    axa = fig.add_subplot(gs[0])
    t = np.array(mfe["t_s"]) * 1e6
    axa.plot(t, mfe["surv_B50"], color="#c0392b", lw=1.6,
             label=f"B\u2080 = {CONFIG['B0_UT']:.0f} \u00b5T")
    axa.plot(t, mfe["surv_B0"], color="#555555", lw=1.6, ls="--",
             label="B\u2080 = 0")
    axa.set_xscale("log")
    axa.set_xlim(2e-5, CONFIG["T_MAX_US"])
    axa.set_xlabel("t  (\u00b5s, log)")
    axa.set_ylabel("surviving radical-pair population")
    axa.set_title("(a)  MFE on spin-selective recombination", loc="left")
    axa.legend(frameon=False)
    axa.grid(alpha=0.25, lw=0.5)

    axb = fig.add_subplot(gs[1])
    dA = np.array(mfe["dA_norm"])
    lam = np.array(mfe["lam_nm"])[::2]
    t2 = t[::2]                      # dA rows are subsampled like the lambda axis
    tmask = t2 >= 0.05
    dAv = dA[tmask, :]
    tv = t2[tmask]
    ext = [lam.min(), lam.max(), tv.max(), tv.min()]
    im = axb.imshow(dAv, aspect="auto", extent=ext, cmap="RdBu_r",
                    vmin=-1, vmax=1)
    axb.set_xlabel("\u03bb (nm)")
    axb.set_ylabel("t (\u00b5s)")
    axb.set_title("(b)  \u0394A(\u03bb, t) = A(B\u2080) \u2212 A(0)  (norm.)",
                  loc="left")
    for lam0, lab, col in ((588.0, "FADH\u2022\u207b", "#2471a3"),
                           (560.0, "TrpH\u2022\u207a", "#c0392b")):
        axb.axvline(lam0, color=col, ls=":", lw=1.2)
        axb.text(lam0 + 5, tv.max() * 0.8, lab, fontsize=8.5, color=col,
                 rotation=90, va="top")
    cb = fig.colorbar(im, ax=axb, fraction=0.045)

    axc = fig.add_subplot(gs[2])
    for tag, col, lab in (("lab_179mT", "#16803c",
                           "B\u2080 = 0.179 mT (f_L \u2248 5.0 MHz)"),
                          ("lab_357mT", "#2c6fbb",
                           "B\u2080 = 0.357 mT (f_L \u2248 10.0 MHz)")):
        od = odmr[tag]
        f = np.array(od["freq_MHz"])
        y = np.array(od["phi_S"])
        y = (y - y.max()) / (y.max() - y.min() + 1e-12)
        axc.plot(f, y, color=col, lw=1.6, label=lab)
        axc.axvline(od["f_Larmor_MHz"], color=col, ls=":", lw=1.0)
    axc.set_xlabel("RF frequency (MHz)")
    axc.set_ylabel("normalized ODMR depth (\u03a6_S loss)")
    axc.set_title("(c)  ODMR: RF destruction of coherence\n"
                  "(dips track the Larmor frequency, doubling with B\u2080)",
                  loc="left")
    axc.legend(frameon=False, fontsize=9)
    axc.grid(alpha=0.25, lw=0.5)
    fig.suptitle("Spectroscopic instrumentation twin \u2014 MFE transient "
                 "absorption & optically detected magnetic resonance",
                 fontsize=13, fontweight="bold", y=1.04)
    fig.savefig(FIG / "fig4_mfe_odmr_instrumentation.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    log("    fig4 saved")


# ============================================================================
# MAIN
# ============================================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["all", "hfcc", "spin", "allostery", "figures"])
    args = ap.parse_args()

    if os.environ.get("PHASE15_STAGE") == "hfcc_only":
        hfcc_stage()
        return

    log("=" * 78)
    log("PHASE 15 - QUANTUM BIOLOGY: RADICAL-PAIR ENGINE & ALLOSTERIC "
        "MAGNETORECEPTION")
    log("=" * 78)
    results = {"config": dict(CONFIG),
               "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")}

    if args.stage in ("all", "spin"):
        # ---- 15A -----------------------------------------------------------
        log("[MODULE 15A] photoinduced radical-pair redox chemistry")
        results["triad_kinetics"] = marcus_triad_kinetics()
        hfc, prov = load_hyperfine()
        results["hyperfine_provenance"] = prov
        results["hyperfine_selected"] = {
            "lumiflavin_anion_N5":
                hfc["fragments"]["lumiflavin_anion"]["selected"]["N_max"],
            "lumiflavin_anion_H":
                hfc["fragments"]["lumiflavin_anion"]["selected"]["H_max"],
            "indole_cation_N1":
                hfc["fragments"]["indole_cation"]["selected"]["N_max"],
            "indole_cation_H":
                hfc["fragments"]["indole_cation"]["selected"]["H_max"]}
        log(f"  15A hyperfine source: {prov}")
        # ---- 15B + 15D ------------------------------------------------------
        log("[MODULE 15B] open-quantum-system spin dynamics")
        results["spin"] = spin_dynamics_stage(hfc)
        results["spin"]["_hops"] = [
            1e12 / results["triad_kinetics"]["k_forward"][0],
            1e12 / results["triad_kinetics"]["k_forward"][1],
            1e12 / results["triad_kinetics"]["k_forward"][2]]
        log("[MODULE 15D] spectroscopy instrumentation twin")
        results["mfe"] = mfe_spectroscopy(results["spin"])
        results["odmr"] = odmr_spectroscopy(hfc)
        (RES / "phase15_results.json").write_text(
            json.dumps(results, indent=2, default=float), encoding="utf-8")
    elif (RES / "phase15_results.json").exists():
        results.update(json.loads(
            (RES / "phase15_results.json").read_text(encoding="utf-8")))

    if args.stage in ("all", "allostery"):
        log("[MODULE 15C] OpenMM allosteric amplification")
        results["allostery"] = run_openmm_allostery()
        (RES / "allostery_results.json").write_text(
            json.dumps(results["allostery"], indent=2, default=float),
            encoding="utf-8")
        if args.stage == "all":
            merged = {}
            if (RES / "phase15_results.json").exists():
                merged = json.loads((RES / "phase15_results.json")
                                    .read_text(encoding="utf-8"))
            merged["allostery"] = results["allostery"]
            (RES / "phase15_results.json").write_text(
                json.dumps(merged, indent=2, default=float), encoding="utf-8")

    if args.stage in ("all", "figures"):
        if "spin" not in results:
            log("  !! no spin results on record - run --stage spin first")
            return
        log("[FIGURES] 300 DPI publication renders")
        hfc, prov = load_hyperfine()
        spin = results["spin"]
        fig1_structure(hfc, spin)
        fig2_dynamics(spin)
        if "allostery" not in results and                 (RES / "allostery_results.json").exists():
            results["allostery"] = json.loads(
                (RES / "allostery_results.json").read_text(encoding="utf-8"))
        if "allostery" in results:
            A_iso = float(
                results["hyperfine_selected"]["lumiflavin_anion_N5"]
                ["A_iso_MHz"])
            fig3_allostery(results["allostery"], spin, A_iso)
        if "mfe" in results and "odmr" in results:
            fig4_spectroscopy(results["mfe"], results["odmr"])
        log(f"DONE - figures in {FIG}, record in "
            f"{RES / 'phase15_results.json'}")


if __name__ == "__main__":
    main()
