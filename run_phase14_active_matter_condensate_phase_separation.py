#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_phase14_active_matter_condensate_phase_separation.py
PHASE 14 GRAND CONVERGENCE MISSION — ACTIVE-MATTER PHASE SEPARATION,
NON-EQUILIBRIUM CONDENSATES & BIOCHEMICAL DISSIPATION.

Unites Physical Chemistry, Biochemistry, Organic Chemistry and Statistical
Mechanics into one multi-scale continuum model of membraneless biomolecular
organelles (nucleoli, stress granules, FUS condensates): liquid-liquid phase
separation (LLPS) driven by an active, ATP-consuming enzymatic reaction
network.

MODULE 14A — Molecular interaction grammar (Organic / Inorganic basis)
    A 20-letter amino-acid pairwise contact energy matrix eps_ij (kBT units)
    explicitly parameterizing
      * cation-pi interactions       Arg/Lys/His <-> Phe/Tyr/Trp
      * aromatic pi-pi stacking      Phe/Tyr/Trp x {Phe,Tyr,Trp}
      * hydrophobic patterning       aliphatic clusters
      * electrostatics under Debye-Huckel screening
        kappa_D = sqrt(8 pi N_A l_B * I), ionic strength I including the
        divalent cations Mg2+ / Zn2+ (which also act as carboxylate bridges)
    The matrix is contracted over the composition of a FUS-like low-complexity
    IDP (sticker-spacer architecture) into the Flory-Huggins parameter
        chi0(T, I) = chi_water - (z_c/2) * <eps>(T, I) * (300 K / T)
    with a composition-dependent correction chi(phi) = chi0 (1 + alpha phi).

MODULE 14B — Active Cahn-Hilliard phase-field engine (Physical Chemistry)
    Total free energy functional
        F[phi] = int [ f_FH(phi) + kappa/2 |grad phi|^2 ] d r
        f_FH   = phi/N ln phi + (1-phi) ln(1-phi) + chi(phi) phi (1-phi)
    coupled to an ATP-dependent phosphorylation cycle (Zwicker-Hyman-Juelicher
    active-droplet mechanism).  Two conserved fields:
        phi = fraction of droplet-forming (unmodified) protein  [phase A]
        psi = fraction of phosphorylated, soluble protein       [phase B]
        d phi/dt = div( M(phi) grad mu ) + Gamma(phi, psi; ATP)
        d psi/dt = div( D_B grad psi )   - Gamma(phi, psi; ATP)
        mu       = delta F / delta phi
        Gamma    = k_deph * psi * phi^2/(K_M^2 + phi^2)  -  k_ATP * phi
                   (dephosphorylation = condensation, Michaelis-Menten
                    phosphatase partitioned into the dense phase;
                    phosphorylation = ATP-driven dispersal, k_ATP ~ [ATP])
    Total protein phi + psi is conserved by construction; the k = 0 Fourier
    mode is touched by nothing but the exactly canceling reaction pair.

    Numerics: semi-implicit pseudo-spectral integrator on a periodic 2-D box.
    The linear part  M0 k^2 (f''(phi_bar) - kappa k^2)  is treated implicitly
    (unconditionally damped at the grid scale), the non-linear remainder, the
    concentration-dependent mobility excess M1*phi and the reaction are
    explicit.  The spectral Laplacian / divergence act in divergence form, so
    global protein mass is conserved to machine precision outside the active
    reaction bookkeeping (which itself cancels between the two fields).

    Non-equilibrium thermodynamics: the continuous entropy production rate
        S_dot = (1/T) * [ int M |grad mu|^2 d r  +  J_cycle * dG_ATP ]
    with J_cycle = <k_ATP phi> (= <k_deph psi h(phi)> at steady state) and
    dG_ATP ~ 19.4 kBT is integrated every frame, proving that the non-equilibrium
    steady-state droplet size is bought with continuous dissipation.

MODULE 14C — Analytical physical fingerprint simulation (Analytical Chemistry)
    1. FRAP: a Gaussian beam bleaches a circular core region of the largest
       droplet; the bleached-fraction density of both protein populations is
       propagated with the real mobility/turnover rates, I(t) is extracted,
       tau_1/2 -> D_app (Axelrod 1976) -> Stokes-Einstein droplet viscosity.
    2. SAXS: the macroscopic structure factor S(q) = <|phi_hat(q)|^2>_angle,
       Porod-law fit I(q) ~ q^-d_f and the microphase peak q* -> domain
       spacing d* = 2 pi / q*.

DELIVERABLES (300 DPI, ./figures_phase14/)
    fig1_active_droplet_spatiotemporal.png
    fig2_thermodynamic_entropy_dissipation.png
    fig3_analytical_frap_saxs_twin.png
plus  results_phase14/phase14_results.json  (machine-readable record).

Key references
    Flory (1942); Huggins (1941); Cahn & Hilliard, J. Chem. Phys. 28, 258
    (1958); Bray, Adv. Phys. 51, 481 (2002); Hyman, Weber & Julicher, Annu.
    Rev. Cell Dev. Biol. 30, 39 (2014); Zwicker, Hyman & Julicher, Phys. Rev.
    E 92, 012317 (2015); Weber & Brangwynne, Cell 149, 1188 (2012); Wang
    et al. Science 361, eaav4382 (2018) [FUS molecular grammar]; Axelrod
    et al., J. Membr. Biol. 13, 7 (1976); Glatter & Kratky, Small Angle
    X-ray Scattering (1982); Kratzer et al. eLife 14 (2025) [ATP-driven
    condensate dissolution].
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import cg, LinearOperator

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# ----------------------------------------------------------------------------
# global configuration
# ----------------------------------------------------------------------------
SEED            = 14
FIG_DIR         = Path("figures_phase14")
RES_DIR         = Path("results_phase14")

T_K             = 310.0          # physiological temperature [K]
KB_KCAL         = 0.0019872      # kcal/mol/K
DG_ATP_KBT      = 19.4           # Delta-G ATP hydrolysis [kBT] (~50 kJ/mol @ 310 K)
KB_J            = 1.380649e-23   # J/K (Stokes-Einstein)
Z_C             = 6              # lattice coordination number (contact grammar)

# --- grid / dynamics (nominal production settings) ---
N_GRID          = 160            # 2-D lattice points per side
L_BOX           = 24.0           # box size [um]   -> dx = 0.15 um
DT              = 0.25           # time step [s]
T_END           = 1000.0         # integration window [s] per the mission brief

NP_CHAIN        = 6              # effective sticker-block chain length N
KAPPA           = 0.02           # gradient penalty [kBT um^2] (interface ~0.3 um)
M0              = 0.010          # mobility M0 [um^2 s^-1 (kBT)^-1]
M1              = 0.005          # concentration-dependent mobility excess
DB_DIFF         = 2.0            # diffusivity of phosphorylated protein [um^2/s]
PHI_BAR         = 0.27           # mean fraction of droplet-forming protein
PSI_BAR         = 0.08           # mean fraction of phosphorylated protein
D_A_SLOW        = 0.06           # diffusivity of droplet-phase protein [um^2/s]
K_DEPH          = 0.09           # phosphatase turnover [s^-1]
KM_MM           = 0.15           # Michaelis constant of the dense-phase gate
K_ATP_NOMINAL   = 0.02           # nominal ATP-driven phosphorylation rate [s^-1]
CHI_WATER       = 0.35           # residual water-IDP incompatibility
ALPHA_CHI       = 0.35           # cooperative composition dependence of chi

SNAP_TIMES      = (1.0, 10.0, 60.0, 300.0, 1000.0)
FRAME_EVERY     = 25             # diagnostic frames (scalars)

# --- ionic conditions (Module 14A electrostatics) ---
C_NACL          = 0.150          # [M] monovalent background
C_MG            = 0.002          # [M] free Mg2+
C_ZN            = 5.0e-5         # [M] free Zn2+

RESULTS: dict = {"meta": {}, "module14A": {}, "module14B": {}, "module14C": {}}


def log(msg: str) -> None:
    print(f"[phase14] {msg}", flush=True)


# ============================================================================
# MODULE 14A — MOLECULAR INTERACTION GRAMMAR
# ============================================================================
AA20 = "ARNDCQEGHILKFPWYSTV"   # canonical one-letter alphabet (standard order)

# class, formal charge at pH 7, hydrophobicity [0..1], aromatic pi scale
AA_INFO = {
    "A": ("aliphatic", 0.0, 0.62, 0.0), "R": ("cation", +1.0, 0.00, 0.0),
    "N": ("polar",     0.0, 0.06, 0.0), "D": ("anion",  -1.0, 0.00, 0.0),
    "C": ("aliphatic", 0.0, 0.68, 0.0), "Q": ("polar",   0.0, 0.00, 0.0),
    "E": ("anion",    -1.0, 0.00, 0.0), "G": ("glycine", 0.0, 0.16, 0.0),
    "H": ("cation",   +0.1, 0.24, 0.5), "I": ("aliphatic", 0.0, 1.00, 0.0),
    "L": ("aliphatic", 0.0, 0.90, 0.0), "K": ("cation", +1.0, 0.00, 0.0),
    "M": ("aliphatic", 0.0, 0.74, 0.0), "F": ("aromatic", 0.0, 0.76, 0.8),
    "P": ("proline",   0.0, 0.44, 0.0), "W": ("aromatic", 0.0, 0.88, 1.0),
    "Y": ("aromatic",  0.0, 0.52, 0.9), "S": ("polar",    0.0, 0.10, 0.0),
    "T": ("polar",     0.0, 0.20, 0.0), "V": ("aliphatic", 0.0, 0.82, 0.0),
}

# cation-pi strength (|eps| at contact, kBT @ 300 K, before pi-scale)
CATION_PI = {"R": -3.6, "K": -2.6, "H": -1.4}
SPACER_CLASSES = ("polar", "glycine")   # G/S/Q/N/T-rich low-complexity spacers


def debye_kappa(T: float, c_nacl: float, c_mg: float, c_zn: float) -> tuple[float, float]:
    """Debye-Huckel screening kappa_D [nm^-1] and ionic strength [M].

    kappa_D = sqrt(8 pi N_A l_B * 1000 * I)  (I in mol/L),  l_B = Bjerrum
    length of water at T (0.714 nm at 298 K, scales as 298/T).
    """
    l_b = 0.714 * (298.0 / T)
    ionic = 0.5 * (c_nacl * 1**2 * 2 + c_mg * 4 + c_zn * 4)   # NaCl: 2 ions
    # kappa_D [m^-1] = sqrt(8 pi N_A l_B(m) * c(mol/m^3));  c = 1000 * I
    kappa = 1e-9 * np.sqrt(8.0 * np.pi * 6.022e23 * (l_b * 1e-9) * 1000.0 * ionic)
    return float(kappa), float(ionic)


def build_contact_matrix(T: float, c_nacl: float, c_mg: float,
                         c_zn: float) -> tuple[np.ndarray, dict]:
    """20x20 pairwise contact energies eps_ij in kBT units (300 K reference).

    Short-range chemistry (hydration, packing, pi-stacking, cation-pi) is
    modeled directly; electrostatics carry the Debye-Huckel screening factor
    exp(-kappa_D r_ij) at a contact distance r_ij = 0.35 nm (salt bridges,
    same-charge repulsion).  Divalent cations (Mg2+/Zn2+) add an extra
    carboxylate-bridge attraction between D/E pairs, saturated in [Mg2+].
    """
    kappa, ionic = debye_kappa(T, c_nacl, c_mg, c_zn)
    scr = np.exp(-kappa * 0.35)                      # screened contact factor
    # divalent carboxylate bridge: -B * [Mg]/([Mg]+Kd), Kd ~ 2 mM
    bridge = 2.4 * (c_mg / (c_mg + 0.002)) * scr * 2.0

    eps = np.zeros((20, 20))
    for a, i in enumerate(AA20):
        for b, j in enumerate(AA20):
            if b < a:
                continue
            ca, _, h_i, p_i = AA_INFO[i]
            cb, q_j, h_j, p_j = AA_INFO[j]
            e = 0.0
            # spacer cohesion: Q/N ladders + backbone H-bonding + poor-solvent
            # quality of the G/S/Q-rich low-complexity spacers (the second,
            # sticker-independent driver of FUS-family LLPS)
            if ca in SPACER_CLASSES and cb in SPACER_CLASSES:
                e += -0.60
            # hydrophobic patterning (aliphatic clustering)
            if {"aliphatic"} & {ca, cb}:
                e += -0.95 * h_i * h_j
            # aromatic pi-pi stacking (Trp strongest)
            if p_i > 0.0 and p_j > 0.0:
                e += -2.2 * np.sqrt(p_i * p_j)
            # cation-pi (Arg ~ 3x Lys ~ 8x His)
            pi_key = (i if ca == "cation" else j)
            pi_aa = (j if ca == "cation" else i)
            if ca == "cation" and p_j > 0.0:
                e += CATION_PI[pi_key] * p_j
            if cb == "cation" and p_i > 0.0:
                e += CATION_PI[pi_key] * p_i
            # electrostatics under Debye-Huckel screening
            if ca == "cation" and cb == "anion":       # salt bridge
                e += -3.6 * scr * 2.0
            if ca == "anion" and cb == "cation":
                e += -3.6 * scr * 2.0
            if ca == cb and ca in ("cation", "anion"): # same-charge repulsion
                e += +2.4 * scr * 2.0
            # anion-pi (weak, cation-pi competitor)
            if "aromatic" in (ca, cb) and "anion" in (ca, cb):
                e += -0.35
            if {ca, cb} == {"polar", "aliphatic"}:     # hydrophobic-polar frustration
                e += +0.45
            if {ca, cb} == {"aromatic", "polar"}:
                e += -0.25
            # Mg2+/Zn2+ carboxylate bridge between acidic side chains
            if {ca, cb} == {"anion", "anion"}:
                e += -bridge
            eps[a, b] = eps[b, a] = e

    meta = {"kappa_D_nm^-1": kappa, "ionic_strength_M": ionic,
            "screening_factor_at_0.35nm": float(scr),
            "Mg_bridge_kBT": float(bridge)}
    return eps, meta


def build_idp_sequence(seed: int = SEED) -> tuple[str, dict]:
    """FUS-like low-complexity model IDP (165 aa, sticker-spacer grammar).

    Composition and motifs follow the FUS N-terminal LC domain family
    (G/S/Q/Y-rich with SYG/GYG stickers and RG/RGG boxes).  This is a model
    sequence with FUS-like statistics, not a UniProt copy.
    """
    rng = np.random.default_rng(seed)
    stickers = ["SYG", "GYG", "SY", "YG", "GYY", "QSY"]
    spacers = ["G", "GG", "SG", "QG", "GN", "GS", "QQG", "TGG", "DG", "EG"]
    boxes = ["RGG", "RG", "RGGG"]

    seq: list[str] = []
    while len(seq) < 165:
        roll = rng.random()
        if roll < 0.12 and len(seq) < 150:
            seq.extend(list(boxes[int(rng.integers(len(boxes)))]))
        elif roll < 0.55:
            seq.extend(list(stickers[int(rng.integers(len(stickers)))]))
        else:
            seq.extend(list(spacers[int(rng.integers(len(spacers)))]))
    seq = "".join(seq[:165])

    comp = {aa: seq.count(aa) / len(seq) for aa in AA20}
    return seq, comp


_CHI0_CACHE: dict = {}


def nominal_chi0() -> float:
    """Nominal chi0 of the FUS-like model IDP at baseline ionic conditions."""
    if "chi0" not in _CHI0_CACHE:
        eps, _ = build_contact_matrix(T_K, C_NACL, C_MG, C_ZN)
        _, comp = build_idp_sequence()
        _CHI0_CACHE["chi0"] = chi_base(T_K, eps, comp)[0]
    return _CHI0_CACHE["chi0"]


def chi_base(T: float, eps: np.ndarray, comp: dict) -> tuple[float, float]:
    """Contract the contact grammar into the Flory-Huggins parameter.

    chi0(T) = chi_water - (z_c/2) * <eps> * (300 K / T);
    <eps> = sum_ij f_i f_j eps_ij (kBT@300K units).
    Favorable contacts (<eps> < 0) therefore RAISE chi (UCST behavior) and
    ionic screening / lower [Mg2+] LOWER chi, following the molecular grammar.
    """
    eps_mean = float(sum(comp[i] * comp[j] * eps[a, b]
                         for a, i in enumerate(AA20) for b, j in enumerate(AA20)))
    chi0 = CHI_WATER - 0.5 * Z_C * eps_mean * (300.0 / T)
    return float(chi0), eps_mean


# ============================================================================
# MODULE 14B — ACTIVE CAHN-HILLIARD PHASE-FIELD ENGINE
# ============================================================================
class ActiveCondensateSim:
    """Semi-implicit spectral Cahn-Hilliard solver with ATP-driven A<->B turnover."""

    def __init__(self, n: int = N_GRID, box: float = L_BOX,
                 k_atp: float = K_ATP_NOMINAL, s_chi: float = 1.0,
                 chi0_in: float | None = None, seed: int = SEED,
                 dt: float = DT, t_end: float = T_END):
        self.n, self.box, self.dx = n, box, box / n
        self.dt, self.dt_max, self.t_end = dt, dt, t_end
        self._clean_steps = 0
        self.k_atp, self.s_chi = k_atp, s_chi
        self.np_chain, self.kappa = NP_CHAIN, KAPPA
        self.m0, self.m1, self.db = M0, M1, DB_DIFF
        self.k_deph, self.km = K_DEPH, KM_MM
        self.alpha_chi = ALPHA_CHI
        self.chi0 = s_chi * (chi0_in if chi0_in is not None else nominal_chi0())

        # Fourier wavenumbers for rfft2 output shape (n, n//2+1)
        k1d = 2.0 * np.pi * np.fft.fftfreq(n, d=self.dx)
        kx = k1d[None, :]
        ky = k1d[:, None]
        self.kx, self.ky = kx[:, : n // 2 + 1], ky
        self.k2 = self.kx**2 + self.ky**2

        rng = np.random.default_rng(seed + int(k_atp * 1e4))
        self.phi = PHI_BAR + 0.02 * rng.standard_normal((n, n))
        self.psi = PSI_BAR + 0.02 * rng.standard_normal((n, n))
        self.mass0 = float((self.phi + self.psi).mean())
        self.t = 0.0

        self.snapshots: dict[float, np.ndarray] = {}
        self.frames: list[dict] = []
        self.max_mass_drift = 0.0
        self.min_dt_used = dt
        self.n_rejections = 0
        self.work = 0.0   # cumulative ATP free-energy dissipation proxy

    # ---- thermodynamics ---------------------------------------------------
    def chi(self, phi: np.ndarray) -> np.ndarray:
        return self.chi0 * (1.0 + self.alpha_chi * phi)

    def fprime(self, phi: np.ndarray) -> np.ndarray:
        pc = np.clip(phi, 1e-5, 1.0 - 1e-5)
        ch = self.chi(pc)
        dch = self.chi0 * self.alpha_chi
        return ((np.log(pc) + 1.0) / self.np_chain
                - (np.log1p(-pc) + 1.0)
                + ch * (1.0 - 2.0 * pc)
                + dch * pc * (1.0 - pc))

    def fpp(self, phi: float) -> float:
        pc = min(max(phi, 1e-3), 1.0 - 1e-3)
        ch = self.chi0 * (1.0 + self.alpha_chi * pc)
        dch = self.chi0 * self.alpha_chi
        return (1.0 / (self.np_chain * pc) + 1.0 / (1.0 - pc)
                - 2.0 * ch + 2.0 * dch * (1.0 - 2.0 * pc))

    def free_energy_density(self) -> float:
        pc = np.clip(self.phi, 1e-5, 1.0 - 1e-5)
        ch = self.chi(pc)
        f = ((pc / self.np_chain) * np.log(pc)
             + (1.0 - pc) * np.log1p(-pc)
             + ch * pc * (1.0 - pc))
        grad2 = self._grad_sq(self.phi)
        return float(np.mean(f + 0.5 * self.kappa * grad2))

    # ---- spectral helpers -------------------------------------------------
    def _grad_sq(self, f: np.ndarray) -> np.ndarray:
        fh = np.fft.rfft2(f)
        gx = np.fft.irfft2(1j * self.kx * fh, s=(self.n, self.n))
        gy = np.fft.irfft2(1j * self.ky * fh, s=(self.n, self.n))
        return gx * gx + gy * gy

    def reaction(self) -> np.ndarray:
        """Net A-production rate Gamma = k_deph psi h(phi) - k_ATP phi."""
        h = self.phi * self.phi / (self.km * self.km + self.phi * self.phi)
        return self.k_deph * h * self.psi - self.k_atp * self.phi

    def cycle_flux(self) -> float:
        """ATP-hydrolysis (phosphorylation) flux <k_ATP phi> [conc/s]."""
        return float(np.mean(self.k_atp * self.phi))

    # ---- single semi-implicit step (guarded, adaptive) ---------------------
    def step(self) -> None:
        """One adaptive semi-implicit step; dt persists across calls.

        On guard rejection dt shrinks 40% and stays small; after 40 clean
        steps it relaxes back toward the nominal value.
        """
        for _ in range(24):
            if self._attempt(self.dt):
                self.min_dt_used = min(self.min_dt_used, self.dt)
                self._clean_steps += 1
                if self._clean_steps > 40:
                    self.dt = min(self.dt_max, self.dt * 1.10)
                return
            self.n_rejections += 1
            self._clean_steps = 0
            self.dt *= 0.6
        raise RuntimeError("stability guard exhausted; dt could not be rescued")

    def _attempt(self, dt: float) -> bool:
        n, dx = self.n, self.dx
        phi, psi = self.phi, self.psi
        phihat = np.fft.rfft2(phi)
        psihat = np.fft.rfft2(psi)

        # linearization curvature around the current mean of phi
        fpp0 = self.fpp(float(phi.mean()))

        # explicit chemical-potential remainder  mu_nl = f'(phi) - fpp0*phi
        mu_nl = self.fprime(phi) - fpp0 * phi
        mu_nl_hat = np.fft.rfft2(mu_nl)

        # full mu (explicit fluxes need the complete chemical potential)
        mu_full_hat = mu_nl_hat + (fpp0 + self.kappa * self.k2) * phihat
        mu_full = np.fft.irfft2(mu_full_hat, s=(n, n))
        mux = np.fft.irfft2(1j * self.kx * mu_full_hat, s=(n, n))
        muy = np.fft.irfft2(1j * self.ky * mu_full_hat, s=(n, n))

        # concentration-dependent mobility excess: J1 = -M1 phi grad(mu)
        jx = -self.m1 * phi * mux
        jy = -self.m1 * phi * muy
        div_j1 = (np.fft.rfft2(jx) * 1j * self.kx
                  + np.fft.rfft2(jy) * 1j * self.ky)

        # explicit reaction
        gam = self.reaction()
        gam_hat = np.fft.rfft2(gam)

        # implicit denominator: 1 + dt M0 k^2 (fpp0 + kappa k^2)
        g_imp = 1.0 + dt * self.m0 * self.k2 * (fpp0 + self.kappa * self.k2)
        g_psi = 1.0 + dt * self.db * self.k2

        # stability guards
        if g_imp.min() < 0.25:
            return False
        # explicit mobility CFL: dt * M1 * kmax^2 * max|mu| < 0.8
        kmax2 = float(self.k2.max())
        if dt * self.m1 * kmax2 * float(np.abs(mu_full).max()) > 0.8:
            return False

        phi_new = np.fft.irfft2(
            (phihat + dt * (-self.m0 * self.k2 * mu_nl_hat - div_j1)) / g_imp,
            s=(n, n)) + dt * np.fft.irfft2(gam_hat, s=(n, n))
        psi_new = np.fft.irfft2((psihat - dt * gam_hat) / g_psi, s=(n, n))

        if phi_new.min() < -0.03 or phi_new.max() > 1.03:
            return False
        if psi_new.min() < -0.03 or psi_new.max() > 1.03:
            return False

        self.phi, self.psi = phi_new, psi_new
        self.t += dt

        # entropy production (this frame, mean density per um^2)
        m_loc = self.m0 * (1.0 + 0.5 * phi)     # mobility used for |grad mu|^2
        grad_mu_sq = mux * mux + muy * muy
        s_diff = float(np.mean(m_loc * grad_mu_sq))
        j_cyc = self.cycle_flux()
        s_chem = j_cyc * DG_ATP_KBT
        self.frames.append({
            "t": self.t, "phi_mean": float(phi.mean()), "psi_mean": float(psi.mean()),
            "S_diff": s_diff, "S_chem": s_chem,
            "S_total": (s_diff + s_chem) / T_K,          # kB / um^2 / s
            "cycle_flux": j_cyc, "F_density": self.free_energy_density(),
            **self.droplet_metrics(),
        })
        drift = abs(float((self.phi + self.psi).mean()) - self.mass0) / self.mass0
        self.max_mass_drift = max(self.max_mass_drift, drift)
        return True

    # ---- diagnostics ------------------------------------------------------
    def droplet_metrics(self) -> dict:
        phi = self.phi
        p5, p95 = np.percentile(phi, 5), np.percentile(phi, 95)
        if p95 - p5 < 0.12:                       # homogeneous (dissolved)
            return {"n_droplets": 0, "R_mean_um": 0.0,
                    "area_fraction": 0.0, "phase_separated": 0}
        thr = 0.5 * (p5 + p95)
        mask = phi > thr
        lab, nlab = ndimage.label(mask, structure=np.ones((3, 3)))
        if nlab == 0:
            return {"n_droplets": 0, "R_mean_um": 0.0,
                    "area_fraction": 0.0, "phase_separated": 0}
        areas = np.bincount(lab.ravel())[1:] * self.dx * self.dx
        big = areas[areas > 30 * self.dx * self.dx]      # noise filter
        r_eq = np.sqrt(big / np.pi)
        return {"n_droplets": int(len(big)),
                "R_mean_um": float(r_eq.mean()) if len(big) else 0.0,
                "area_fraction": float(mask.mean()),
                "phase_separated": 1}

    def run(self, snapshot_times=SNAP_TIMES) -> None:
        t_next_snap = 0
        pending = sorted(snapshot_times)
        while self.t < self.t_end - 1e-9:
            dt = min(self.dt, self.t_end - self.t)
            # capture snapshots that fall inside the upcoming step
            while t_next_snap < len(pending) and pending[t_next_snap] <= self.t + dt + 1e-9:
                self.snapshots[float(pending[t_next_snap])] = self.phi.copy()
                t_next_snap += 1
            self.step()
            if len(self.frames) % 400 == 1:
                fr = self.frames[-1]
                log(f"  t={fr['t']:7.1f}s  <phi>={fr['phi_mean']:.3f} "
                    f"<psi>={fr['psi_mean']:.3f}  R={fr['R_mean_um']:.2f}um "
                    f"N={fr['n_droplets']:3d}  S_dot={fr['S_total']:.3e} kB/um2/s")
        if len(self.snapshots) < len(pending):
            self.snapshots[float(pending[-1])] = self.phi.copy()

    def ness_stats(self, last_frac: float = 0.2) -> dict:
        tail = self.frames[int(len(self.frames) * (1.0 - last_frac)):]
        return {
            "S_diff_kB_um2_s": float(np.mean([f["S_diff"] for f in tail]) / T_K),
            "S_chem_kB_um2_s": float(np.mean([f["S_chem"] for f in tail]) / T_K),
            "S_total_kB_um2_s": float(np.mean([f["S_total"] for f in tail])),
            "cycle_flux_M_s": float(np.mean([f["cycle_flux"] for f in tail])),
            "R_mean_um": float(np.mean([f["R_mean_um"] for f in tail])),
            "n_droplets": float(np.mean([f["n_droplets"] for f in tail])),
            "area_fraction": float(np.mean([f["area_fraction"] for f in tail])),
        }


# ============================================================================
# MODULE 14C — ANALYTICAL FINGERPRINTS: FRAP & SAXS
# ============================================================================
def locate_largest_droplet(phi: np.ndarray, dx: float) -> tuple[float, float, float] | None:
    """Return (x, y, R_eq) of the largest droplet, or None if homogeneous."""
    p5, p95 = np.percentile(phi, 5), np.percentile(phi, 95)
    if p95 - p5 < 0.12:
        return None
    thr = 0.5 * (p5 + p95)
    lab, nlab = ndimage.label(phi > thr, structure=np.ones((3, 3)))
    if nlab == 0:
        return None
    sizes = np.bincount(lab.ravel())
    sizes[0] = 0
    k = int(sizes.argmax())
    ys, xs = np.nonzero(lab == k)
    r_eq = np.sqrt(sizes[k] * dx * dx / np.pi)
    return float(xs.mean() * dx), float(ys.mean() * dx), float(r_eq)


def simulate_frap(sim: ActiveCondensateSim, k_atp: float,
                  t_sim: float = 300.0, dt_f: float = 0.25) -> dict:
    """FRAP twin: Gaussian-beam bleaching of the largest droplet core.

    Standard reaction-diffusion FRAP model (Sprague et al. 2004) applied to
    the frozen post-run condensate field: a single bleached-fraction field
    c(r,t) in [0,1] with

        dc/dt = div( D_eff(r) grad c ) + k_turn(r) * ( f_ref - c )

    D_eff(r) = (D_A phi + D_B psi)/(phi+psi)  is the local mobility-weighted
    diffusivity (slow inside the viscous droplet, fast in the phosphorylated
    dilute sea) and k_turn(r) = (k_ATP phi + k_deph h psi)/(phi+psi) is the
    local turnover-driven exchange with the unbleached reservoir (f_ref =
    global unbleached fraction, known analytically).  Both limiting regimes
    are reproduced: passive (k_ATP=0) recovery is pure slow droplet diffusion,
    while ATP-driven turnover accelerates the apparent recovery exactly as
    measured in active condensates.  The implicit FV system is SPD (monotone,
    unconditionally stable), so c rises monotonically to f_ref.

    tau_1/2 (recovery-rate fit of ln(1-Fn)) -> D_app = 0.224 r0^2 / tau_1/2
    (Axelrod 1976) -> Stokes-Einstein droplet viscosity.
    """
    phi, psi = sim.phi, sim.psi
    n, dx = sim.n, sim.dx
    loc = locate_largest_droplet(phi, dx)
    if loc is None:   # dissolved condensate: bleach the box center instead
        cx = cy = sim.box / 2.0
        r_d = 0.8
    else:
        cx, cy, r_d = loc
    r0 = max(0.6 * r_d, 0.45)             # 1/e2 bleach radius
    w0 = 0.8 * r0

    X, Y = np.meshgrid(dx * np.arange(n), dx * np.arange(n))
    r2 = (X - cx) ** 2 + (Y - cy) ** 2

    # Gaussian beam, exp(-2 r^2/w0^2) profile, peak bleaching 78%
    bleach = np.exp(-2.0 * r2 / w0**2) * 1.5
    c = np.exp(-bleach)                    # surviving fluorophore fraction
    roi = r2 <= (1.35 * r0) ** 2

    prot = np.maximum(phi + psi, 1e-3)
    roi_prot0 = float(prot[roi].sum())
    fluo0 = float((c * prot)[roi].sum())
    fluoro_tot = float((c * prot).sum())
    f_ref = fluoro_tot / float(prot.sum())   # exact long-time asymptote

    # local mobility-weighted diffusivity and turnover rate
    phi_t = np.maximum(phi, 0.0)
    psi_t = np.maximum(psi, 0.0)
    d_eff = (D_A_SLOW * phi_t + DB_DIFF * psi_t) / np.maximum(prot, 1e-3)
    d_eff = np.clip(d_eff, D_A_SLOW, DB_DIFF)
    k_turn = (k_atp * phi_t + sim.k_deph
              * (phi_t * phi_t / (sim.km**2 + phi_t * phi_t)) * psi_t) / prot

    # implicit FV operator with harmonic-mean face weights (SPD M-matrix)
    yy, xx = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    def _faces(part):
        wx = 2.0 * part[:, :-1] * part[:, 1:] / np.maximum(part[:, :-1] + part[:, 1:], 1e-12)
        wy = 2.0 * part[:-1, :] * part[1:, :] / np.maximum(part[:-1, :] + part[1:, :], 1e-12)
        wwx = 2.0 * part[:, -1] * part[:, 0] / np.maximum(part[:, -1] + part[:, 0], 1e-12)
        wwy = 2.0 * part[-1, :] * part[0, :] / np.maximum(part[-1, :] + part[0, :], 1e-12)
        return wx, wy, wwx, wwy
    wxf, wyf, wwx, wwy = _faces(d_eff)
    axf = -wxf / dx**2
    ayf = -wyf / dx**2
    axw = -wwx / dx**2
    ayw = -wwy / dx**2
    r_list = [(yy[:, :-1] * n + xx[:, :-1]).ravel(), (yy[:, :-1] * n + xx[:, 1:]).ravel(),
              (yy[:-1, :] * n + xx[:-1, :]).ravel(), (yy[1:, :] * n + xx[1:, :]).ravel(),
              np.arange(n) * n + (n - 1), np.arange(n) * n,
              (n - 1) * n + np.arange(n), np.arange(n)]
    c_list = [(yy[:, :-1] * n + xx[:, 1:]).ravel(), (yy[:, :-1] * n + xx[:, :-1]).ravel(),
              (yy[1:, :] * n + xx[1:, :]).ravel(), (yy[:-1, :] * n + xx[:-1, :]).ravel(),
              np.arange(n) * n, np.arange(n) * n + (n - 1),
              np.arange(n), (n - 1) * n + np.arange(n)]
    d_list = [axf.ravel(), axf.ravel(), ayf.ravel(), ayf.ravel(),
              axw, axw, ayw, ayw]
    rows_f = np.concatenate(r_list)
    cols_f = np.concatenate(c_list)
    data_f = np.concatenate(d_list)
    diag_a = -np.bincount(rows_f.astype(int), weights=data_f, minlength=n * n)
    Amat = csr_matrix((np.append(data_f, diag_a),
                       (np.append(rows_f, np.arange(n * n)).astype(int),
                        np.append(cols_f, np.arange(n * n)).astype(int))),
                      shape=(n * n, n * n))
    # resolve the fastest recovery: tau_min ~ 0.224 r0^2 / D_eff_max
    dt_f = min(dt_f, 0.03 * r0**2 / float(d_eff.max()))
    inv_dt = 1.0 / dt_f
    kdiag = k_turn.ravel()
    lop = LinearOperator((n * n, n * n),
                         matvec=lambda v: inv_dt * v + kdiag * v + Amat @ v)
    precond = LinearOperator((n * n, n * n),
                             matvec=lambda v: v / (inv_dt + kdiag + diag_a))
    kbar = (k_turn * prot).sum() / prot.sum()   # area-averaged turnover

    times, fluoro = [0.0], [fluo0]
    t = 0.0
    x0 = None
    plateau_run = 0
    while t < t_sim - 1e-9:
        rhs = c.ravel() * inv_dt + kdiag * f_ref
        sol, info = cg(lop, rhs, M=precond, rtol=1e-7, atol=1e-12,
                       maxiter=250, x0=x0)
        if info != 0 or not np.all(np.isfinite(sol)):
            sol = (c.ravel() + dt_f * kdiag * f_ref) / (1.0 + dt_f * kdiag)
        x0 = sol
        c = sol.reshape(n, n)
        t += dt_f
        times.append(t)
        fluoro.append(float((c * prot)[roi].sum()))
        # early exit once the recovery has flattened (fast regimes)
        if len(fluoro) > 120:
            tail = np.asarray(fluoro[-40:])
            if tail.max() - tail.min() < 1e-5 * (tail.max() - fluo0 + 1e-12):
                plateau_run += 1
                if plateau_run >= 2:
                    break
            else:
                plateau_run = 0

    times = np.asarray(times)
    fluoro = np.asarray(fluoro)
    F = fluoro / roi_prot0
    F0 = float(F[0])
    F_end = float(F[-1])
    denom = F_end - F0
    Fn = (F - F0) / denom if denom > 1e-9 else np.full_like(F, np.nan)

    # recovery-rate fit: ln(1-Fn) is linear on the initial recovery segment
    seg = (Fn > 0.08) & (Fn < 0.60)
    if seg.sum() >= 6:
        slope, _ = np.polyfit(times[seg], np.log(np.maximum(1.0 - Fn[seg], 1e-12)), 1)
        tau_rec = -1.0 / slope if slope < -1e-6 else float("nan")
    else:
        tau_rec = float("nan")
    tau_half = float(np.log(2.0) * tau_rec) if np.isfinite(tau_rec) else float("nan")
    if np.isfinite(tau_half) and tau_half < dt_f:
        tau_half = dt_f   # under-resolved: report the resolution bound

    if np.isfinite(tau_half) and tau_half > 0:
        d_app = 0.224 * r0**2 / tau_half                   # Axelrod (2-D disk)
    else:
        d_app = float("nan")
    a_h = 2.5e-9                                          # hydrodynamic radius [m]
    eta = (KB_J * T_K / (6.0 * np.pi * d_app * 1e-12 * a_h)
           if np.isfinite(d_app) and d_app > 0 else float("nan"))
    return {"t": times, "F": F, "Fn": Fn, "tau_half_s": tau_half,
            "tau_recovery_s": float(tau_rec), "k_turn_mean": float(kbar),
            "D_app_um2_s": float(d_app), "viscosity_Pa_s": float(eta),
            "r0_um": float(r0), "roi_center": (cx, cy),
            "k_atp": k_atp, "F_ref": float(f_ref), "F0": F0}


def saxs_intensity(phi: np.ndarray, dx: float) -> dict:
    """Macroscopic structure factor S(q) = <|phi_hat(q)|^2>_azimuthal.

    q* (microphase/domain scale) is read from the maximum of the Kratky-type
    q^2 I(q); the Porod exponent d_f is fitted on I(q) ~ q^-d_f over the
    window between the peak and the interface-crossover q ~ 1/l_interface.
    """
    n = phi.shape[0]
    fl = phi - phi.mean()
    F = np.fft.fftshift(np.abs(np.fft.fft2(fl)) ** 2) / n**4
    k1d = 2.0 * np.pi * np.fft.fftfreq(n, d=dx)
    kx = np.repeat(k1d[None, :], n, 0)
    ky = np.repeat(k1d[:, None], n, 1)
    q = np.fft.fftshift(np.sqrt(kx**2 + ky**2))   # shift q identically to F
    qmax = k1d.max() * np.sqrt(2)
    edges = np.geomspace(2 * np.pi / (n * dx), qmax * 0.999, 120)
    qc, Sq = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (q >= lo) & (q < hi)
        if m.sum() >= 4:
            qc.append(0.5 * (lo + hi))
            Sq.append(F[m].mean())
    qc, Sq = np.asarray(qc), np.asarray(Sq)

    # microphase peak: argmax of q^2 I(q) over the resolved low-q window
    # (beyond ~0.5 q_nyquist the q^2 weighting amplifies the flat noise floor)
    homogenous = float(np.percentile(phi, 95) - np.percentile(phi, 5)) < 0.12
    q_lo_min = 2.0 * 2 * np.pi / (n * dx)
    q_cap = 0.5 * k1d.max()
    kr = (qc**2) * Sq
    kr[(qc < q_lo_min) | (qc > q_cap) | homogenous] = 0.0
    i_star = int(np.argmax(kr))
    if kr[i_star] <= 0.0:
        q_star, s_star = float("nan"), float("nan")
    else:
        q_star, s_star = float(qc[i_star]), float(Sq[i_star])

    # Porod-type decay exponent: search the most linear window of log I
    # vs log q inside [1.1 q*, 0.5 q_nyq] (local-slope analysis; the clean
    # 2-D Porod regime q^-3 requires droplets much larger than the interface
    # width, so the operational exponent is reported with its window)
    best = None
    if not homogenous and np.isfinite(q_star):
        q_lo_min = 1.1 * q_star
        q_hi_max = 0.5 * k1d.max()
        nb = len(qc)
        for i0 in range(nb):
            if qc[i0] < q_lo_min:
                continue
            for i1 in range(i0 + 4, nb):
                if qc[i1] > q_hi_max:
                    break
                x, y = np.log(qc[i0:i1 + 1]), np.log(Sq[i0:i1 + 1])
                if not np.all(np.isfinite(y)):
                    continue
                slope, inter = np.polyfit(x, y, 1)
                r2 = 1.0 - np.var(y - (slope * x + inter)) / max(np.var(y), 1e-30)
                if best is None or r2 > best["r2"]:
                    best = {"slope": float(slope), "d_f": float(-slope),
                            "r2": float(r2),
                            "q_range": [float(qc[i0]), float(qc[i1])]}
    porod = best if best is not None else {
        "slope": float("nan"), "d_f": float("nan"), "r2": float("nan"),
        "q_range": []}
    return {"q": qc, "S_q": Sq, "q_star": q_star,
            "d_star_um": float(2 * np.pi / q_star),
            "porod": porod}


# ============================================================================
# EXPERIMENTS
# ============================================================================
def experiment_main(s_quick: float) -> dict:
    """E1 — passive vs active LLPS, t in [0, 1000 s] (fig 1)."""
    log("E1: passive (k_ATP=0) vs active (k_ATP=0.02 s^-1) condensate dynamics")
    t_end = 300.0 if s_quick < 1 else T_END
    runs = {}
    for tag, k_atp in (("passive", 0.0), ("active", K_ATP_NOMINAL)):
        sim = ActiveCondensateSim(k_atp=k_atp, t_end=t_end,
                                  dt=max(0.25, DT * s_quick),
                                  n=N_GRID if s_quick == 1 else 112)
        log(f"  {tag}: N={sim.n}, dt={sim.dt}, T_end={t_end}")
        t0 = time.perf_counter()
        sim.run(SNAP_TIMES)
        runs[tag] = sim
        log(f"  {tag}: done in {time.perf_counter()-t0:.1f}s, "
            f"mass drift {sim.max_mass_drift:.2e}, rejections {sim.n_rejections}, "
            f"NESS {sim.ness_stats()}")
    return runs


def experiment_phase_diagram(s_quick: float) -> dict:
    """E2 — condensate stability vs ATP rate at several chi (fig 2)."""
    log("E2: phase diagram scan over (k_ATP, chi-scale)")
    if s_quick < 1:
        k_list = [0.0, 0.02, 0.08]
        s_list = [1.0, 1.15]
    else:
        k_list = [0.0, 0.004, 0.010, 0.020, 0.040, 0.080, 0.150]
        s_list = [0.85, 1.0, 1.15, 1.30]
    grid = {}
    for s_chi in s_list:
        for k_atp in k_list:
            sim = ActiveCondensateSim(k_atp=k_atp, s_chi=s_chi,
                                      n=112 if s_quick == 1 else 96,
                                      dt=0.5,
                                      t_end=T_END if s_quick == 1 else 400.0)
            t0 = time.perf_counter()
            sim.run(snapshot_times=(T_END,))
            st = sim.ness_stats()
            grid[f"s{s_chi:.2f}_k{k_atp:.3f}"] = {
                "s_chi": s_chi, "k_atp": k_atp, "chi0": sim.chi0, **st,
                "mass_drift": sim.max_mass_drift,
                "runtime_s": time.perf_counter() - t0}
            log(f"  chi-scale {s_chi:.2f} (chi0={sim.chi0:.3f}), "
                f"k_ATP={k_atp:.3f}: R={st['R_mean_um']:.2f}um, "
                f"A={st['area_fraction']:.3f}, S_dot={st['S_total_kB_um2_s']:.3e}, "
                f"({time.perf_counter()-t0:.1f}s)")
    return grid


def experiment_fingerprints(main_runs: dict, s_quick: float) -> dict:
    """E3 — FRAP + SAXS twins at three ATP levels (fig 3)."""
    log("E3: FRAP & SAXS fingerprints at k_ATP in {0, 0.02, 0.08}")
    out = {}
    for k_atp in (0.0, K_ATP_NOMINAL, 0.08):
        sim = ActiveCondensateSim(k_atp=k_atp,
                                  n=128 if s_quick == 1 else 96,
                                  dt=0.5,
                                  t_end=T_END if s_quick == 1 else 400.0)
        sim.run(snapshot_times=(T_END,))
        frap = simulate_frap(sim, k_atp,
                             t_sim=150.0 if s_quick < 1 else 300.0)
        saxs = saxs_intensity(sim.phi, sim.dx)
        out[f"k{k_atp:.3f}"] = {
            "k_atp": k_atp, "chi0": sim.chi0,
            "tau_half_s": frap["tau_half_s"],
            "D_app_um2_s": frap["D_app_um2_s"],
            "viscosity_Pa_s": frap["viscosity_Pa_s"],
            "r0_um": frap["r0_um"],
            "q_star": saxs["q_star"], "d_star_um": saxs["d_star_um"],
            "porod_slope": saxs["porod"]["slope"],
            "porod_d_f": saxs["porod"]["d_f"],
            "ness": sim.ness_stats(),
            "mass_drift": sim.max_mass_drift,
            "_frap_curve": (frap["t"], frap["Fn"]),
            "_saxs_curve": (saxs["q"], saxs["S_q"]),
            "_saxs": saxs,
        }
        log(f"  k_ATP={k_atp}: tau_1/2={frap['tau_half_s']:.2f}s, "
            f"D_app={frap['D_app_um2_s']:.4f} um^2/s, "
            f"eta={frap['viscosity_Pa_s']:.3f} Pa.s, q*={saxs['q_star']:.2f}, "
            f"d_f={saxs['porod']['d_f']:.2f}")
    return out


# ============================================================================
# FIGURES (300 DPI)
# ============================================================================
def _stylize(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def make_fig1(runs: dict, out: Path) -> None:
    tags = [("passive", "passive equilibrium (k$_{ATP}$ = 0)"),
            ("active", f"active NESS (k$_{{ATP}}$ = {K_ATP_NOMINAL} s$^{{-1}}$)")]
    times = [10.0, 60.0, 300.0, 1000.0]
    fig, axes = plt.subplots(2, 4, figsize=(13.5, 7.2), constrained_layout=True)
    ext = [0, runs["active"].box, 0, runs["active"].box]
    for r, (tag, label) in enumerate(tags):
        sim = runs[tag]
        for c, t_snap in enumerate(times):
            ax = axes[r, c]
            field = sim.snapshots.get(t_snap)
            if field is None:
                field = sim.phi
            im = ax.imshow(field, origin="lower", extent=ext, cmap="magma",
                           vmin=0.0, vmax=0.7, interpolation="bilinear")
            ax.set_title(f"t = {t_snap:g} s", fontsize=11)
            if c == 0:
                ax.set_ylabel(f"{label}\nposition [$\\mu$m]", fontsize=10)
            if r == 1:
                ax.set_xlabel("position [$\\mu$m]", fontsize=10)
            ax.set_xticks([0, 12, 24]); ax.set_yticks([0, 12, 24])
            ax.tick_params(labelsize=8)
    # annotate the physics: mean droplet radius where phase-separated
    for r, (tag, _) in enumerate(tags):
        sim = runs[tag]
        m = sim.droplet_metrics()
        note = (f"$\\langle R\\rangle$ = {m['R_mean_um']:.2f} $\\mu$m, "
                f"N = {m['n_droplets']}, area = {m['area_fraction']*100:.0f}%"
                if m["phase_separated"] else "homogeneous (dissolved)")
        axes[r, -1].text(0.985, 0.02, note, transform=axes[r, -1].transAxes,
                         ha="right", va="bottom", fontsize=9.5, color="w",
                         bbox=dict(boxstyle="round,pad=0.25", fc="k", alpha=0.55))
    cb = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.015, aspect=38)
    cb.set_label(r"droplet-forming protein fraction $\phi(\mathbf{r},t)$", fontsize=11)
    fig.suptitle(
        "Fig. 1 — Active LLPS of a FUS-like IDP: spontaneous phase separation, growth,\n"
        f"and non-equilibrium steady-state size suppression by ATP-driven turnover "
        f"(2-D Cahn-Hilliard, {N_GRID}$^2$ grid, {L_BOX:.0f} $\\mu$m box)",
        fontsize=12.5)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    log(f"saved {out}")


def make_fig2(scan: dict, out: Path) -> None:
    ks = sorted({v["k_atp"] for v in scan.values()})
    ss = sorted({v["s_chi"] for v in scan.values()})
    R = np.zeros((len(ss), len(ks)))
    A = np.zeros_like(R)
    S = np.zeros_like(R)
    for i, s_chi in enumerate(ss):
        for j, k in enumerate(ks):
            v = scan[f"s{s_chi:.2f}_k{k:.3f}"]
            R[i, j], A[i, j], S[i, j] = (v["R_mean_um"], v["area_fraction"],
                                         max(v["S_total_kB_um2_s"], 1e-12))
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.9), constrained_layout=True)
    ax = axes[0]
    for i, s_chi in enumerate(ss):
        chi0 = scan[f"s{s_chi:.2f}_k{ks[0]:.3f}"]["chi0"]
        ax.semilogx(np.array(ks) + 1e-4, R[i], "o-", lw=2, ms=5,
                    label=fr"$\chi_0$ = {chi0:.2f}")
    ax.set_xlabel(r"ATP hydrolysis rate $k_{ATP}$ [s$^{-1}$]", fontsize=11)
    ax.set_ylabel(r"steady-state mean droplet radius $\langle R\rangle$ [$\mu$m]", fontsize=11)
    ax.set_title("(a) active size regulation", fontsize=12)
    ax.legend(fontsize=9); _stylize(ax)

    ax = axes[1]
    Kgrid, Sgrid = np.meshgrid(np.array(ks) + 1e-4, np.asarray(ss))
    pc = ax.pcolormesh(Kgrid, Sgrid, A, cmap="cividis", shading="nearest",
                       vmin=0, vmax=max(0.35, A.max()))
    cs = ax.contour(Kgrid, Sgrid, np.log10(S), levels=6, colors="w", linewidths=1.2)
    ax.clabel(cs, fmt=r"$\dot S$: %.1f", fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel(r"$k_{ATP}$ [s$^{-1}$]", fontsize=11)
    ax.set_ylabel(r"$\chi$-scale $s_\chi$", fontsize=11)
    ax.set_title("(b) stability (color) vs entropy generation (contours, log$_{10}$"
                 r" $\dot S$ [k$_B\,\mu$m$^{-2}$s$^{-1}$])", fontsize=11)
    cb = fig.colorbar(pc, ax=ax, shrink=0.9)
    cb.set_label("condensate area fraction", fontsize=10)

    ax = axes[2]
    for i, s_chi in enumerate(ss):
        ax.loglog(np.array(ks)[1:] + 1e-4, S[i, 1:], "s-", lw=2, ms=5,
                  label=fr"$\chi_0$ = {scan[f's{s_chi:.2f}_k{ks[0]:.3f}']['chi0']:.2f}")
    ax.set_xlabel(r"$k_{ATP}$ [s$^{-1}$]", fontsize=11)
    ax.set_ylabel(r"$\dot S_{prod}$ [k$_B\,\mu$m$^{-2}$ s$^{-1}$]", fontsize=11)
    ax.set_title("(c) continuous entropy generation\n"
                 r"$\dot S_{chem}\approx\langle k_{ATP}\phi\rangle\,\Delta G_{ATP}/T$",
                 fontsize=12)
    ax.legend(fontsize=9); _stylize(ax)
    fig.suptitle("Fig. 2 — Condensate stability vs ATP hydrolysis: dissipation buys "
                 "the non-equilibrium steady state", fontsize=13)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    log(f"saved {out}")


def make_fig3(fp: dict, out: Path) -> None:
    keys = sorted(fp.keys())
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.9), constrained_layout=True)

    ax = axes[0]
    colors = plt.cm.viridis(np.linspace(0.15, 0.8, len(keys)))
    for col, k in zip(colors, keys):
        v = fp[k]
        t, Fn = v["_frap_curve"]
        ax.plot(t, Fn, color=col, lw=2.2,
                label=fr"$k_{{ATP}}$={v['k_atp']:g} s$^{{-1}}$, "
                      fr"$\tau_{{1/2}}$={v['tau_half_s']:.1f} s")
        ax.axhline(1.0, color="k", lw=0.6, alpha=0.3)
    ax.set_xlabel("time after bleach [s]", fontsize=11)
    ax.set_ylabel(r"normalized recovery $I(t)$", fontsize=11)
    ax.set_title("(a) FRAP recovery (droplet core)", fontsize=12)
    ax.set_xlim(0, min(60, fp[keys[-1]]["_frap_curve"][0][-1]))
    ax.legend(fontsize=8.5, loc="lower right"); _stylize(ax)

    ax = axes[1]
    kvals = [fp[k]["k_atp"] for k in keys]
    dapps = [fp[k]["D_app_um2_s"] for k in keys]
    etas = [fp[k]["viscosity_Pa_s"] for k in keys]
    ax.semilogx(np.array(kvals) + 1e-4, dapps, "o-", lw=2.2, color="#d95f02",
                label=r"$D_{app}$ (left)")
    ax.set_xlabel(r"$k_{ATP}$ [s$^{-1}$]", fontsize=11)
    ax.set_ylabel(r"$D_{app}$ [$\mu$m$^2$/s]", color="#d95f02", fontsize=11)
    ax2 = ax.twinx()
    ax2.semilogx(np.array(kvals) + 1e-4, etas, "s--", lw=2.2, color="#7570b3",
                 label=r"$\eta$ (right)")
    ax2.set_ylabel(r"droplet viscosity $\eta$ [Pa$\cdot$s]", color="#7570b3", fontsize=11)
    ax.set_title("(b) mobility & viscosity (Stokes-Einstein)", fontsize=12)
    _stylize(ax); ax2.spines["top"].set_visible(False)

    ax = axes[2]
    for col, k in zip(colors, keys):
        v = fp[k]
        q, Sq = v["_saxs_curve"]
        ax.loglog(q, Sq, color=col, lw=2,
                  label=fr"$k_{{ATP}}$={v['k_atp']:g} s$^{{-1}}$"
                        + (fr", $q^*$={v['q_star']:.1f}" if np.isfinite(v["q_star"]) else ""))
    # Porod guide slopes through the top curve
    qg = np.geomspace(6.0, 35.0, 20)
    ref = fp[keys[0]]["_saxs_curve"][1]
    scale = float(np.nanmax(ref[-6:]))
    for slope, style in ((-3.0, "--"), (-4.0, ":")):
        ax.loglog(qg, scale * (qg / qg[-1]) ** slope, "k" + style, lw=1.1)
    ax.text(0.97, 0.90, r"guide: $q^{-3}$ (Porod, 2-D)", transform=ax.transAxes,
            ha="right", fontsize=8.5)
    ax.text(0.97, 0.84, r"$q^{-4}$ (3-D Porod ref.)", transform=ax.transAxes,
            ha="right", fontsize=8.5)
    ax.set_xlabel(r"$q$ [$\mu$m$^{-1}$]", fontsize=11)
    ax.set_ylabel(r"$I(q)\sim S(q)$ [a.u.]", fontsize=11)
    ax.set_title("(c) SAXS structure factor S(q)", fontsize=12)
    ax.legend(fontsize=8.5, loc="lower left"); _stylize(ax)

    fig.suptitle("Fig. 3 — Analytical physical fingerprint twin: FRAP mobility vs "
                 "SAXS scattering of the active condensate", fontsize=13)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    log(f"saved {out}")


# ============================================================================
# MAIN
# ============================================================================
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--quick", action="store_true",
                    help="reduced grid/time smoke run for CI")
    args = ap.parse_args()
    s_quick = 0.35 if args.quick else 1.0

    FIG_DIR.mkdir(exist_ok=True)
    RES_DIR.mkdir(exist_ok=True)
    t_wall = time.perf_counter()

    # ---------------- Module 14A: molecular grammar ------------------------
    log("Module 14A: amino-acid contact grammar, Debye-Huckel screening, chi(phi,T)")
    eps, ion_meta = build_contact_matrix(T_K, C_NACL, C_MG, C_ZN)
    seq, comp = build_idp_sequence()
    chi0_nominal, eps_mean = chi_base(T_K, eps, comp)
    RESULTS["module14A"] = {
        "amino_acid_order": list(AA20),
        "contact_matrix_kBT": eps.round(4).tolist(),
        "ionic_conditions": {"NaCl_M": C_NACL, "MgCl2_M": C_MG, "ZnCl2_M": C_ZN,
                             **ion_meta},
        "idp_sequence": seq,
        "idp_composition": {k: round(v, 4) for k, v in comp.items() if v > 0},
        "chi0_310K": chi0_nominal,
        "chi_formula": "chi0(T,I) = 0.35 - (z_c/2)<eps>(T,I)(300/T); "
                       "chi(phi) = chi0 (1 + 0.35 phi)",
        "chi_crit_FH": float(0.5 * (1 + 1 / NP_CHAIN**0.5) ** 2),
        "T_scan": {f"{T:g}K": chi_base(T, eps, comp)[0]
                   for T in (290, 300, 310, 323, 340)},
        "salt_scan": {f"I={s}": chi_base(T_K, build_contact_matrix(
                          T_K, s, C_MG, C_ZN)[0], comp)[0]
                      for s in (0.01, 0.05, 0.15, 0.5, 1.0)},
        "mg_scan": {f"[Mg2+]={m}": chi_base(T_K, build_contact_matrix(
                        T_K, C_NACL, m, C_ZN)[0], comp)[0]
                    for m in (0.0, 0.001, 0.002, 0.005, 0.02)},
    }
    log(f"  sequence ({len(seq)} aa): {seq[:60]}...")
    log(f"  chi0(310K) = {chi0_nominal:.3f}  [chi_crit = "
        f"{RESULTS['module14A']['chi_crit_FH']:.3f}]  -> LLPS: "
        f"{chi0_nominal > RESULTS['module14A']['chi_crit_FH']}")
    log(f"  kappa_D = {ion_meta['kappa_D_nm^-1']:.3f} nm^-1 "
        f"(I = {ion_meta['ionic_strength_M']:.3f} M), "
        f"screening factor {ion_meta['screening_factor_at_0.35nm']:.3f}")

    # ---------------- Experiments ------------------------------------------
    main_runs = experiment_main(s_quick)
    scan = experiment_phase_diagram(s_quick)
    fp = experiment_fingerprints(main_runs, s_quick)

    # ---------------- figures ----------------------------------------------
    log("rendering 300-DPI figures")
    make_fig1(main_runs, FIG_DIR / "fig1_active_droplet_spatiotemporal.png")
    make_fig2(scan, FIG_DIR / "fig2_thermodynamic_entropy_dissipation.png")
    make_fig3(fp, FIG_DIR / "fig3_analytical_frap_saxs_twin.png")

    # ---------------- machine-readable record ------------------------------
    for tag, sim in main_runs.items():
        RESULTS["module14B"][f"E1_{tag}"] = {
            "k_atp": sim.k_atp, "dt_used": sim.dt, "n_grid": sim.n,
            "t_end": sim.t_end, "max_mass_drift": sim.max_mass_drift,
            "n_stability_rejections": sim.n_rejections,
            "min_dt_used": sim.min_dt_used,
            "ness": sim.ness_stats(),
            "frames_every": [{"t": f["t"], **{k: v for k, v in f.items() if k != "t"}}
                             for f in sim.frames[:: max(1, len(sim.frames) // 80)]],
        }
    RESULTS["module14B"]["E2_phase_diagram"] = {k: {kk: vv for kk, vv in v.items()}
                                                for k, v in scan.items()}
    for k, v in fp.items():
        RESULTS["module14C"][f"E3_{k}"] = {kk: vv for kk, vv in v.items()
                                           if not kk.startswith("_")}
        RESULTS["module14C"][f"E3_{k}"]["frap_curve"] = {
            "t": v["_frap_curve"][0].round(3).tolist()[::4],
            "Fn": v["_frap_curve"][1].round(5).tolist()[::4]}
        RESULTS["module14C"][f"E3_{k}"]["saxs_curve"] = {
            "q": v["_saxs_curve"][0].round(5).tolist(),
            "I": v["_saxs_curve"][1].round(5).tolist()}
    RESULTS["meta"] = {
        "phase": 14, "script": Path(__file__).name,
        "model": "two-field active Cahn-Hilliard (Zwicker-Hyman-Julicher) with "
                 "Flory-Huggins free energy from a 20-letter contact grammar",
        "grid": f"{N_GRID}x{N_GRID}", "box_um": L_BOX, "dt_s": DT, "t_end_s": T_END,
        "kappa_kBT_um2": KAPPA, "M0": M0, "M1": M1, "D_B": DB_DIFF,
        "phi_bar": PHI_BAR, "psi_bar": PSI_BAR,
        "k_deph": K_DEPH, "KM": KM_MM, "DG_ATP_kBT": DG_ATP_KBT, "T_K": T_K,
        "wall_time_s": round(time.perf_counter() - t_wall, 1),
        "numpy": np.__version__,
    }
    with open(RES_DIR / "phase14_results.json", "w", encoding="utf-8") as fh:
        json.dump(RESULTS, fh, indent=1)
    log(f"results -> {RES_DIR/'phase14_results.json'}  "
        f"(wall {RESULTS['meta']['wall_time_s']} s)")

    # ---------------- stability certificate --------------------------------
    worst = max([sim.max_mass_drift for sim in main_runs.values()]
                + [v["mass_drift"] for v in scan.values()]
                + [v["mass_drift"] for v in fp.values()])
    log(f"MASS-CONSERVATION CERTIFICATE: max relative drift of <phi+psi> = {worst:.3e}")
    assert worst < 1e-8, "protein mass was not conserved to machine precision"
    log("PHASE 14 COMPLETE.")


if __name__ == "__main__":
    main()
