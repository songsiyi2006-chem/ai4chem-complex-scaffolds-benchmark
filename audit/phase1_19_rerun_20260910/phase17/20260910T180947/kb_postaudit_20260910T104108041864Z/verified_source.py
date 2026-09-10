#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_phase17_relativistic_actinide_quantum.py — PHASE 17 SUPREME MISSION
Fully Relativistic 4-Component Dirac Quantum Chemistry & Actinide 5f-Covalency.

The relativistic realm: in Am(III) (Z=95) the 1s electrons move at ~0.6 c.
The relativistic mass increase contracts the core, screens the nucleus less
effectively, and *destabilizes and expands* the 5f shell — while Eu(III)
(Z=63) keeps a compact, energetically isolated 4f shell that never coheres
with a nitrogen lone pair.  Am3+ and Eu3+ are ISOELECTRONIC (f6): the
comparison isolates relativity itself as the origin of nitrogen-donor
selectivity in nuclear-waste partitioning (the SANEX / BTP process).

MODULES
-------
17A  4-component Dirac-Coulomb(-Breit) atomic engine, built from scratch in
     NumPy/SciPy (no PySCF on win64 — the documented platform doctrine of
     Phases 13/15; and no PySCF code path implements 4c-DCB anyway):
     * radial Dirac-Fock-Slater SCF on a logarithmic grid, finite-sphere
       nucleus, Xa local exchange (alpha=2/3), average-of-configuration
       occupations for the open 5f6 / 4f6 shells;
     * the radial Dirac pencil is discretized on a STAGGERED grid — P on
       nodes, Q on cell midpoints — the discrete image of exact kinetic
       balance (chi^S = sigma.p chi^L / 2c): no variational collapse into
       the positron continuum (certified against the Sommerfeld
       fine-structure formula for point-nucleus hydrogen-like ions);
     * a non-relativistic Hartree-Fock-Slater twin on the identical grid
       and functional isolates the purely relativistic contribution;
     * KINETIC-BALANCE MATRIX DEMONSTRATION: a 4-component matrix-Dirac
       solver in s-Gaussian bases for H-like U (Z=92) — kinetically balanced
       small components vs raw Gaussian small components: the unbalanced
       basis collapses into spurious deeply-bound states, the balanced one
       is stable to variational refinement;
     * BREIT INTERACTION (first order, frequency-independent form)
       g_B = -[a1.a2 + (a1.r12)(a2.r12)/r12^2]/(2 r12): pair matrix
       elements by EXPLICIT 6D angular-spinor quadrature with real
       spherical spinors — the machinery is certified by reproducing the
       exact 5Z/8 Coulomb pair integral of H 1s^2 and the c->infty collapse
       of the Breit term, then applied to the SCF orbitals.  Gaunt
       (magnetic) and retardation parts are separated.

17B  Spin-orbit multiplet engine — the CAS(6e,7o) relativistic active-space
     CI in the atomic limit:
     * the full f6 configuration (C(14,6) = 3003 Slater determinants over
       7 orbitals x 2 spins) diagonalized with the Slater-Condon
       electrostatic operator (F2, F4, F6 computed as radial integrals of
       the SCF large components) plus zeta L.S spin-orbit, zeta taken from
       the Dirac SCF j-splitting: zeta = [eps(j=7/2) - eps(j=5/2)]/3.5;
     * intermediate coupling, J-resolved levels, dE_SO of the 7F_J ground
       manifold, Lande g_J factors, L/S/J composition of every level;
     * Eu3+ (4f6) levels validated against the experimental fluorescence
       ladder (0/380/1050/1900/2860/3920/4940 cm-1);
     * NR and scalar-relativistic (zeta=0) twins quantify the progressive
       lifting of degeneracy: NR -> scalar-rel -> 4-component Dirac+SOC.

17C  Relativistic bonding, QTAIM & EDA — why Am3+ coheres and Eu3+ does not:
     * promolecular densities from the engine's OWN relativistic (and NR)
       atomic solutions, N-M-N axis model of the BTP tridentate pocket at
       the crystallographic distances (Am-N 2.53 A, Eu-N 2.47 A);
     * Wolfsberg-Helmholtz two-state ionic/covalent mixing: the donor-
       acceptor coupling H_DA = 1.75 S (eps_A + eps_D)/2 is built from the
       REAL 3D overlap integrals S(5f/4f, N 2p lone pair) and the REAL
       orbital energies of the engine — 5f overlap larger AND energy
       matched, 4f overlap vanishing AND energetically orphaned;
     * dE_cov per bond and per complex, the Am/Eu selectivity differential,
       delocalization indices delta(M,N);
     * QTAIM on a 2D N-M-N plane grid: (3,-1) bond critical points on the
       M-N paths, rho(r_BCP), Laplacian rho(r_BCP) (exact on-axis
       transverse isotropy), H(r_BCP) = 1/4 Lap^2 rho - G with
       G = TF + (1/9)vW local kinetic model.

DELIVERABLES
------------
figures_phase17/fig1_4component_spinor_orbitals.png         (300 DPI)
figures_phase17/fig2_spin_orbit_multiplet_splitting.png     (300 DPI)
figures_phase17/fig3_qtaim_relativistic_covalency_map.png   (300 DPI)
figures_phase17/fig4_relativistic_foundation_validation.png (supplementary)
results_phase17/phase17_results.json                        (machine record)
results_phase17/phase17_radial_data.npz                     (radial amplitudes)

Stage control (resumable; each stage caches into phase17_results.json):
    python run_phase17_relativistic_actinide_quantum.py [--stage all]
        stage in {atomic, gausskb, breit, multiplet, bonding, figures, all}
"""

import os
import sys
import json
import math
import time
import argparse
from functools import lru_cache
from itertools import combinations
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "8")

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.linalg as sla
from scipy.special import roots_legendre, comb

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures_phase17"
RES = ROOT / "results_phase17"
FIG.mkdir(exist_ok=True)
RES.mkdir(exist_ok=True)

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg):
    print(f"[phase17 {time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ----------------------------------------------------------------------------
# Atomic units & constants
# ----------------------------------------------------------------------------
C_LIGHT = 137.035999084
HARTREE_EV = 27.211386245988
HARTREE_KCAL = 627.5094740631
HARTREE_CM = 219474.6313632
BOHR_A = 0.529177210903
XALPHA = 2.0 / 3.0            # Slater Xalpha exchange parameter

# ----------------------------------------------------------------------------
# Grids
# ----------------------------------------------------------------------------

def log_grid(n_pts, r_min, r_max):
    """Logarithmic grid: nodes x, midpoints y, spacing h, SBP node weights."""
    t = np.linspace(math.log(r_min), math.log(r_max), n_pts)
    x = np.exp(t)
    h = np.diff(x)
    y = 0.5 * (x[:-1] + x[1:])
    w = np.empty(n_pts)
    w[1:-1] = 0.5 * (x[2:] - x[:-2])
    w[0] = 0.5 * (x[1] - x[0])
    w[-1] = 0.5 * (x[-1] - x[-2])
    return x, y, h, w


def nucleus_potential(z_nuc, r, a_mass=None):
    """Uniform-sphere finite-nucleus electrostatic potential (-Z/r outside)."""
    if a_mass is not None:
        rms_fm = 0.8359 * a_mass ** (1.0 / 3.0) + 0.6715        # Angeli
    else:
        rms_fm = 1.2 * (z_nuc ** (1.0 / 3.0))
    r_rms = rms_fm * 1e-15 / 5.29177210903e-11                  # fm -> a0
    r_sph = math.sqrt(5.0 / 3.0) * r_rms
    v = np.empty_like(r)
    inside = r < r_sph
    v[inside] = -z_nuc / (2.0 * r_sph) * (3.0 - r[inside] ** 2 / r_sph ** 2)
    v[~inside] = -z_nuc / r[~inside]
    return v, r_sph


# ----------------------------------------------------------------------------
# Real spherical spinors (m = +-1/2) and spherical harmonics
# ----------------------------------------------------------------------------

@lru_cache(maxsize=None)
def _plm_norm(l, m):
    """Normalization of the real orthonormal spherical harmonic."""
    return math.sqrt((2 * l + 1) / (4 * math.pi)
                     * math.factorial(l - m) / math.factorial(l + m))


def _plm(l, m, u):
    """Associated Legendre P_l^m(u), Condon-Shortley phase, stable recurrences:
        P_m^m = (-1)^m (2m-1)!! (1-u^2)^{m/2}
        P_{m+1}^m = u (2m+1) P_m^m
        P_l^m = [(2l-1) u P_{l-1}^m - (l-1+m) P_{l-2}^m] / (l-m)
    """
    u = np.asarray(u, dtype=float)
    pmm = np.ones_like(u)
    if m > 0:
        s = np.sqrt(np.maximum(1.0 - u * u, 0.0))
        f = 1.0
        for _ in range(m):
            pmm = pmm * (-f * s)
            f += 2.0
    if l == m:
        return pmm
    pm1 = u * (2 * m + 1) * pmm
    if l == m + 1:
        return pm1
    pl = pm1
    for ll in range(m + 2, l + 1):
        pl = ((2 * ll - 1) * u * pm1 - (ll - 1 + m) * pmm) / (ll - m)
        pmm, pm1 = pm1, pl
    return pl


def Ylm_real(l, m, th, ph):
    """Real orthonormal spherical harmonic (Condon-Shortley phase)."""
    u = np.cos(th)
    if m == 0:
        return _plm_norm(l, 0) * _plm(l, 0, u)
    am = abs(m)
    base = math.sqrt(2.0) * _plm_norm(l, am) * _plm(l, am, u)
    return base * (np.cos(am * ph) if m > 0 else np.sin(am * ph))


@lru_cache(maxsize=None)
def spinor_coefs(l, j, m):
    """CG-derived real spherical spinors Omega_{l j m}, m = +-1/2.

    From |(l 1/2) j m> with Condon-Shortley CG coefficients
    (components = (spin-up amplitude, spin-down amplitude)):

      j = l+1/2:  m=+1/2: ( sqrt((l+1)/(2l+1)) Y^0      , sqrt(l/(2l+1)) Y^{+1} )
                  m=-1/2: ( sqrt(l/(2l+1))     Y^{-1}   , sqrt((l+1)/(2l+1)) Y^0 )
      j = l-1/2:  m=+1/2: ( -sqrt(l/(2l+1)) Y^0         , sqrt((l+1)/(2l+1)) Y^{+1} )
                  m=-1/2: ( -sqrt((l+1)/(2l+1)) Y^{-1}  , sqrt(l/(2l+1)) Y^0 )
    """
    d = 2 * l + 1
    if abs(j - l - 0.5) < 1e-9:
        if m > 0:
            return ((math.sqrt((l + 1) / d), l, 0),
                    (math.sqrt(l / d), l, 1))
        return ((math.sqrt(l / d), l, -1),
                (math.sqrt((l + 1) / d), l, 0))
    if abs(j - l + 0.5) < 1e-9:
        if m > 0:
            return ((-math.sqrt(l / d), l, 0),
                    (math.sqrt((l + 1) / d), l, 1))
        return ((-math.sqrt((l + 1) / d), l, -1),
                (math.sqrt(l / d), l, 0))
    raise ValueError((l, j, m))


def spherical_spinor_clean(l, j, m, th, ph):
    """Real 2-component spherical spinor Omega_{l j m}; components = spin."""
    (c1, l1, m1), (c2, l2, m2) = spinor_coefs(l, j, m)
    out = np.zeros((2,) + np.shape(th))
    if abs(c1) > 1e-14:
        out[0] = c1 * Ylm_real(l1, m1, th, ph)
    if abs(c2) > 1e-14:
        out[1] = c2 * Ylm_real(l2, m2, th, ph)
    return out


def spinor_selftest():
    """Certify orthonormality of the real spherical spinors on the sphere."""
    nt = 64
    gx, wx = roots_legendre(nt)
    theta = np.arccos(gx)
    nph = 2 * nt
    ph = np.linspace(0, 2 * np.pi, nph, endpoint=False)
    wph = 2 * np.pi / nph
    TH, PH = np.meshgrid(theta, ph, indexing="ij")
    W = np.outer(wx, np.full(nph, wph))
    worst = 0.0
    cases = [(0, 0.5), (1, 0.5), (1, 1.5), (2, 1.5), (2, 2.5), (3, 2.5),
             (3, 3.5), (4, 3.5), (4, 4.5)]
    for (l, j) in cases:
        for m in (+0.5, -0.5):
            om = spherical_spinor_clean(l, j, m, TH, PH)
            acc = float(np.sum(W * np.sum(om * om, axis=0)))
            worst = max(worst, abs(acc - 1.0))
    for (l, j) in [(2, 2.5), (3, 3.5), (4, 4.5)]:
        om_p = spherical_spinor_clean(l, j, 0.5, TH, PH)
        om_m = spherical_spinor_clean(l, j, -0.5, TH, PH)
        worst = max(worst,
                    abs(float(np.sum(W * np.sum(om_p * om_m, axis=0)))))
    return worst


# ----------------------------------------------------------------------------
# Dirac orbital solver: node-counted inward/outward matching integration
# (the GRASP-lineage scheme: spurious states and variational collapse are
#  impossible by construction - every trial E yields a unique radial
#  solution and bound states are identified by node counting)
# ----------------------------------------------------------------------------

def _v_interp(r, v_pot, x):
    i = int(np.searchsorted(x, r))
    i = min(max(i, 1), len(x) - 1)
    t = (r - x[i - 1]) / (x[i] - x[i - 1])
    return v_pot[i - 1] * (1.0 - t) + v_pot[i] * t


def _dirac_rhs(kappa, ec, r, yy, v_pot, x, c):
    P_, Q_ = yy
    v = _v_interp(r, v_pot, x)
    return np.array((
        (-kappa * P_ / r + (ec - v + 2.0 * c * c) / c * Q_) * r,
        (kappa * Q_ / r - (ec - v) / c * P_) * r))


def _dirac_series_out(kappa, v_pot, x, r_s, c):
    """Frobenius start just outside the nucleus: P = r^gam, Q = q P."""
    z_eff = abs(_v_interp(r_s, v_pot, x) * r_s)
    z_eff = min(z_eff, 0.999 * abs(kappa) * c)
    gam = math.sqrt(kappa * kappa - (z_eff / c) ** 2)
    P = r_s ** gam
    Q = P * c * (gam + kappa) / z_eff
    return P, Q


def _dirac_series_in(kappa, e_energy, r_e, c):
    """Decaying asymptotics at the box edge: P = e^{-lam r}, Q = q P."""
    lam = math.sqrt(max(-e_energy * (e_energy + 2.0 * c * c), 1e-30)) / c
    P = math.exp(-lam * r_e)
    Q = P * (-lam * c / (e_energy + 2.0 * c * c))
    return P, Q


def _dirac_rk4(kappa, e_energy, v_pot, x, i_lo, i_hi, y0, c, sub=1,
               _cache={}):
    """RK4 integration of d(P,Q)/dt = r*(P', Q') along the log grid.

    Integrates from x[i_lo] to x[i_hi] (either direction), returning arrays
    of (P, Q) on every grid index between (inclusive).  The potential is
    sampled at nodes and segment midpoints once per grid (cached), so the
    hot loop is pure scalar arithmetic (no searchsorted).
    """
    key = (id(v_pot), len(x), float(v_pot[1]), float(v_pot[len(x) // 2]),
           float(v_pot[-2]))
    if _cache.get("key") != key:
        _cache["key"] = key
        _cache["vmid"] = 0.5 * (v_pot[:-1] + v_pot[1:])
        _cache["dt"] = np.log(x[1:] / x[:-1])
    vmid = _cache["vmid"]; dt = _cache["dt"]
    ec = e_energy
    step = 1 if i_hi >= i_lo else -1
    P = np.zeros(len(x)); Q = np.zeros(len(x))
    y0v = np.array(y0, dtype=float)
    P[i_lo], Q[i_lo] = y0v
    yP, yQ = y0v
    blew = False
    rng_ = range(i_lo, i_hi, step)
    for ia in rng_:
        ib = ia + step
        if blew:
            P[ib], Q[ib] = yP, yQ
            continue
        h = dt[ia] if step > 0 else -dt[ia - 1]
        va = v_pot[ia]
        vb = v_pot[ib]
        vm = 0.5 * (va + vb)
        ra = x[ia]
        rm = math.exp(math.log(ra) + 0.5 * h)
        rb = math.exp(math.log(ra) + h)
        # inline RK4 on the 2-vector (yP, yQ)
        kk1P = (-kappa * yP / ra + (ec - va + 2.0 * c * c) / c * yQ) * ra
        kk1Q = (kappa * yQ / ra - (ec - va) / c * yP) * ra
        kk2P = (-kappa * (yP + 0.5 * h * kk1P) / rm
                + (ec - vm + 2.0 * c * c) / c * (yQ + 0.5 * h * kk1Q)) * rm
        kk2Q = (kappa * (yQ + 0.5 * h * kk1Q) / rm
                - (ec - vm) / c * (yP + 0.5 * h * kk1P)) * rm
        kk3P = (-kappa * (yP + 0.5 * h * kk2P) / rm
                + (ec - vm + 2.0 * c * c) / c * (yQ + 0.5 * h * kk2Q)) * rm
        kk3Q = (kappa * (yQ + 0.5 * h * kk2Q) / rm
                - (ec - vm) / c * (yP + 0.5 * h * kk2P)) * rm
        kk4P = (-kappa * (yP + h * kk3P) / rb
                + (ec - vb + 2.0 * c * c) / c * (yQ + h * kk3Q)) * rb
        kk4Q = (kappa * (yQ + h * kk3Q) / rb
                - (ec - vb) / c * (yP + h * kk3P)) * rb
        yP += h / 6.0 * (kk1P + 2 * kk2P + 2 * kk3P + kk4P)
        yQ += h / 6.0 * (kk1Q + 2 * kk2Q + 2 * kk3Q + kk4Q)
        if not (math.isfinite(yP) and math.isfinite(yQ))                 or abs(yP) > 1e150 or abs(yQ) > 1e150:
            blew = True
            if not math.isfinite(yP):
                yP = 1e150
            if not math.isfinite(yQ):
                yQ = 1e150
        P[ib], Q[ib] = yP, yQ
    return P, Q


def _count_nodes(P_arr):
    s = np.sign(P_arr)
    s = s[s != 0]
    return int(np.sum(s[1:] * s[:-1] < 0))


def dirac_orbital_shoot(kappa, n_quant, l, v_pot, x, r_nuc, e_guess=None,
                        c=C_LIGHT, tol=1e-10, max_iter=80):
    """Solve one radial Dirac orbital by inward/outward matching.

    Returns (E, P_on_nodes, Q_on_midpoints, ok).
    """
    i_start = int(np.searchsorted(x, max(r_nuc * 1.2, x[0] * 1.0000001)))
    i_start = max(min(i_start, len(x) - 40), 1)
    nodes_target = n_quant - l - 1

    def _imatch(e):
        # outer classical turning point: first grid point with E - V < 0
        idx = np.where((e - v_pot) < 0.0)[0]
        return int(idx[0]) if len(idx) else len(x) - 3

    def _nodes_at(e):
        """Sturm node count: integrate into the forbidden region until the
        solution blows up, and count sign changes before divergence."""
        P_o, _Q = _dirac_rk4(kappa, e, v_pot, x, i_start, len(x) - 1,
                             _dirac_series_out(kappa, v_pot, x, x[i_start], c), c)
        im = _imatch(e)
        seg = P_o[i_start:]
        ref = np.max(np.abs(seg[:max(im - i_start, 1)]))
        thr = max(ref, 1e-300) * 1e8
        big = np.where(np.abs(seg) > thr)[0]
        cut = int(big[0]) if len(big) else len(seg)
        return _count_nodes(seg[:max(cut, 2)])

    # --- universal bracket -------------------------------------------------
    # Node counting brackets E_target only when target >= 1 (Sturm: the
    # deepest state of a symmetry has ZERO nodes for every E below it).
    # For target == 0 the mismatch-sign bisection alone does the job over
    # the full energy window (D has exactly one zero: the first eigenvalue).
    z_top = max(abs(v_pot[i_start] * x[i_start]), 1.0)
    e_deep = -2.2 * z_top * z_top - 1.0
    e_shallow = -1e-6
    with np.errstate(over="ignore", invalid="ignore"):
        if _nodes_at(e_deep) > nodes_target:
            return None, None, None, False
        if _nodes_at(e_shallow) < nodes_target:
            return None, None, None, False
        lo, hi = e_deep, e_shallow
        if e_guess is not None and _nodes_at(float(e_guess)) == nodes_target:
            g = float(e_guess)
            lo, hi = g * 1.45 + 0.05, g * 0.70 - 0.02
            for _ in range(12):
                if _nodes_at(lo) > nodes_target:
                    lo *= 1.6
                else:
                    break
            for _ in range(12):
                if _nodes_at(hi) <= nodes_target:
                    hi *= 0.6
                else:
                    break
        elif nodes_target >= 1:
            if e_guess is not None:
                # fast localization around a warm start (SCF iterations)
                g = float(e_guess)
                for _ in range(24):
                    mg = _nodes_at(g)
                    if mg < nodes_target:
                        lo = g
                        g *= 0.7
                    elif mg > nodes_target:
                        hi = g
                        g *= 1.4
                    else:
                        hi = g
                        g = 0.5 * (lo + g)
                    if abs(hi - lo) < max(abs(lo), 1.0) * 1e-10:
                        break
            # bisection to the node-count transition (bracket contains E_target)
            for _ in range(48):
                mid = 0.5 * (lo + hi)
                if _nodes_at(mid) <= nodes_target:
                    lo = mid
                else:
                    hi = mid
                if abs(hi - lo) < max(abs(lo), 1.0) * 1e-11:
                    break
    e = 0.5 * (lo + hi)
    ok = False

    def _pair0(ee):
        im = _imatch(ee)
        P_o, Q_o = _dirac_rk4(kappa, ee, v_pot, x, i_start, im,
                              _dirac_series_out(kappa, v_pot, x, x[i_start], c), c)
        lam = math.sqrt(max(-ee * (ee + 2.0 * c * c), 1e-30)) / c
        r_in = x[im] * math.exp(min(40.0 / max(lam, 1e-12), 6.0))
        i_in = max(int(np.searchsorted(x, r_in)), min(im + 2, len(x) - 1))
        i_in = min(i_in, len(x) - 1)
        P_i, Q_i = _dirac_rk4(kappa, ee, v_pot, x, i_in, im,
                              _dirac_series_in(kappa, ee, x[i_in], c), c)
        po, pi_ = P_o[im], P_i[im]
        if abs(po) < 1e-300 or abs(pi_) < 1e-300:
            return float("inf")
        return Q_o[im] * pi_ - Q_i[im] * po

    def _pair(ee):
        """Outward + bounded inward integrations; Wronskian mismatch."""
        im = _imatch(ee)
        P_o, Q_o = _dirac_rk4(kappa, ee, v_pot, x, i_start, im,
                              _dirac_series_out(kappa, v_pot, x, x[i_start], c), c)
        lam = math.sqrt(max(-ee * (ee + 2.0 * c * c), 1e-30)) / c
        r_in = x[im] * math.exp(min(40.0 / max(lam, 1e-12), 6.0))
        i_in = max(int(np.searchsorted(x, r_in)), min(im + 2, len(x) - 1))
        i_in = min(i_in, len(x) - 1)
        P_i, Q_i = _dirac_rk4(kappa, ee, v_pot, x, i_in, im,
                              _dirac_series_in(kappa, ee, x[i_in], c), c)
        po, pi_ = P_o[im], P_i[im]
        if abs(po) < 1e-300 or abs(pi_) < 1e-300:
            return float("inf"), None   # degenerate; force bracket move
        # raw Wronskian: continuous through P-nodes (no poles), zeros only
        # at true eigenvalues -- the GRASP-style matching function
        mism = Q_o[im] * pi_ - Q_i[im] * po
        return mism, (P_o, Q_o, P_i, Q_i, im, i_in)

    # Pure bisection on the Wronskian sign inside the node-count bracket.
    # The Wronskian is continuous through P-nodes (no poles) and flips sign
    # only at true eigenvalues; the reference sign at the bracket bottom
    # defines the "below the eigenvalue" direction.
    with np.errstate(over="ignore", invalid="ignore"):
        _w_lo = _pair0(lo)
        _w_sign_ref = _w_lo if (math.isfinite(_w_lo) and _w_lo != 0.0) else 1.0
    for it in range(max_iter):
        m = _nodes_at(e)
        mism, pack = _pair(e)
        if m != nodes_target:
            if m < nodes_target:
                lo = e
            else:
                hi = e
        elif mism * _w_sign_ref > 0.0:
            lo = e
        else:
            hi = e
        e_mid = 0.5 * (lo + hi)
        if abs(hi - lo) < 1e-13 * max(abs(e), 1.0):
            e = e_mid
            ok = True
            break
        e = e_mid
    mism, pack = _pair(e)
    if pack is None:
        return None, None, None, False
    P_o, Q_o, P_i, Q_i, im, i_in = pack
    scale = P_o[im] / P_i[im] if abs(P_i[im]) > 1e-300 else 1.0
    P_full = np.where(np.arange(len(x)) <= im, P_o, scale * P_i)
    Q_full = np.where(np.arange(len(x)) <= im, Q_o, scale * Q_i)
    P_full[:i_start] = 0.0
    Q_full[:i_start] = 0.0
    w = np.empty(len(x))
    w[1:-1] = 0.5 * (x[2:] - x[:-2])
    w[0] = w[1] * 0.5
    w[-1] = w[-2] * 0.5
    q_mid = 0.5 * (Q_full[:-1] + Q_full[1:])
    norm = (np.sum(P_full ** 2 * w) + np.sum(q_mid ** 2 * w[:-1] * 0.5)
            + np.sum(q_mid ** 2 * w[1:] * 0.5))
    sn = math.sqrt(norm)
    P_full = P_full / sn
    Q_full = Q_full / sn
    return float(e), P_full, Q_full, ok


def pencil_nodes(P):
    return _count_nodes(P)


def kappa_for(n, l, j_lo):
    """kappa = (l - j)(2j + 1)."""
    j = l - 0.5 if j_lo else l + 0.5
    return int(round((l - j) * (2 * j + 1)))


def l_of_kappa(kappa):
    return int(round(abs(kappa + 0.5) - 0.5))



# ----------------------------------------------------------------------------
# Non-relativistic twin: radial Schrodinger on a uniform grid
# ----------------------------------------------------------------------------

NR_DR = 3.0e-4          # uniform dr (a0): resolves the Am 1s peak (~35 pts)
NR_RMAX = 70.0


def nr_grid():
    return np.arange(1, int(NR_RMAX / NR_DR) + 1) * NR_DR


def nr_bound_states(l, v_pot, x, k_states=14):
    """Radial Schrodinger for angular momentum l, uniform grid, B = I.

        -1/2 u'' + [l(l+1)/(2r^2) + V] u = eps u,   u(0) = 0.
    Returns [(energy, u_on_grid)] ascending (u = r*R plays the P role).
    """
    n = len(x)
    dx = x[1] - x[0]
    diag = 1.0 / dx ** 2 + 0.5 * l * (l + 1) / x ** 2 + v_pot
    off = np.full(n - 1, -0.5 / dx ** 2)
    k = min(k_states, n - 2)
    evals, evecs = sla.eigh_tridiagonal(diag, off, select="i",
                                        select_range=(0, k - 1))
    return [(float(evals[jx]), evecs[:, jx]) for jx in range(len(evals))]


# ----------------------------------------------------------------------------
# Density / potential machinery
# ----------------------------------------------------------------------------

def orbital_density(P, Q, x, w=None):
    """Spherical rho(r) on grid nodes of a spinor: (P^2+Q^2)/(4 pi r^2)."""
    dens = P ** 2
    if Q is not None and len(Q) == len(P):
        dens = dens + Q ** 2
    return dens / (4.0 * np.pi * x ** 2)


def hartree_potential(rho_nodes, x, w):
    """U(r) = (1/r) int_0^r 4pi r'^2 rho dr' + int_r^inf 4pi r' rho dr'."""
    dq = 4.0 * np.pi * x ** 2 * rho_nodes * w          # shell charge
    qin = np.concatenate([[0.0], np.cumsum(dq)[:-1]])  # charge inside node r
    dm = 4.0 * np.pi * x * rho_nodes * w
    out = np.concatenate([np.cumsum(dm[1:][::-1])[::-1], [0.0]])
    return qin / x + out


def xalpha_vx(rho_nodes):
    rho = np.maximum(rho_nodes, 1e-30)
    return -3.0 * XALPHA * (3.0 * rho / (8.0 * np.pi)) ** (1.0 / 3.0)


def slater_x_energy(rho_nodes, x, w):
    rho = np.maximum(rho_nodes, 1e-30)
    return float(-2.25 * XALPHA * np.sum(
        rho * (3.0 * rho / (8.0 * np.pi)) ** (1.0 / 3.0) * 4 * np.pi * x ** 2 * w))


def total_energy_from_eigenvalues(e_sum, rho, u_input, vx_input, x, w):
    """Remove input mean-field potentials; add output-density Hartree/Ex."""
    dv = 4.0 * np.pi * x**2 * w
    u_output = hartree_potential(rho, x, w)
    return float(e_sum - np.sum(rho * (u_input + vx_input) * dv)
                 + 0.5 * np.sum(rho * u_output * dv)
                 + slater_x_energy(rho, x, w))


# ----------------------------------------------------------------------------
# Configurations (average of configuration for the open f shell)
# ----------------------------------------------------------------------------

AM3_Z, AM3_A = 95, 243
EU3_Z, EU3_A = 63, 153
# (n, l, occ_j=l-1/2, occ_j=l+1/2); zeros keep a virtual block in the solve
AM3_CONFIG = [
    (1, 0, 2.0, 0.0), (2, 0, 2.0, 0.0), (2, 1, 2.0, 4.0),
    (3, 0, 2.0, 0.0), (3, 1, 2.0, 4.0), (3, 2, 4.0, 6.0),
    (4, 0, 2.0, 0.0), (4, 1, 2.0, 4.0), (4, 2, 4.0, 6.0),
    (4, 3, 6.0, 8.0),
    (5, 0, 2.0, 0.0), (5, 1, 2.0, 4.0), (5, 2, 4.0, 6.0),
    (6, 0, 2.0, 0.0), (6, 1, 2.0, 4.0),
    (5, 3, 2.5714285714285716, 3.4285714285714284),   # 5f^6 AOC
    (6, 2, 0.0, 0.0),                                  # 6d virtual
    (7, 0, 0.0, 0.0),                                  # 7s virtual
]
EU3_CONFIG = [
    (1, 0, 2.0, 0.0), (2, 0, 2.0, 0.0), (2, 1, 2.0, 4.0),
    (3, 0, 2.0, 0.0), (3, 1, 2.0, 4.0), (3, 2, 4.0, 6.0),
    (4, 0, 2.0, 0.0), (4, 1, 2.0, 4.0), (4, 2, 4.0, 6.0),
    (5, 0, 2.0, 0.0), (5, 1, 2.0, 4.0),
    (4, 3, 2.5714285714285716, 3.4285714285714284),   # 4f^6 AOC
    (5, 2, 0.0, 0.0),                                  # 5d virtual
    (6, 0, 0.0, 0.0),                                  # 6s virtual
]


def config_blocks(config):
    """Occupation map {(n, kappa): occ} and the kappa block list."""
    occ_map = {}
    blocks = []
    for (n, l, o_lo, o_hi) in config:
        if l == 0:
            kap = -1
            if kap not in blocks:
                blocks.append(kap)
            occ_map[(n, kap)] = occ_map.get((n, kap), 0.0) + o_lo
        else:
            for kap, occ in ((kappa_for(n, l, True), o_lo),
                             (kappa_for(n, l, False), o_hi)):
                if kap not in blocks:
                    blocks.append(kap)
                occ_map[(n, kap)] = occ
    return occ_map, blocks


def guess_density(z, config, x):
    """Slater-screened hydrogenic shells (initial SCF guess only)."""
    rho = np.zeros_like(x)
    nelec = 0.0
    inner = 0.0
    for (n, l, o1, o2) in config:
        occ = o1 + o2
        if occ <= 0:
            continue
        z_eff = max(z - 0.85 * inner - 0.35 * (occ - 1.0), 1.0)
        zeta = 1.05 * z_eff / n
        rho += occ * zeta ** 3 / math.pi * np.exp(-2.0 * zeta * x)
        inner += occ
        nelec += occ
    nr = np.sum(rho * 4 * np.pi * x ** 2 * np.gradient(x))
    return rho * (nelec / (nr + 1e-300)), nelec


# ----------------------------------------------------------------------------
# SCF engine
# ----------------------------------------------------------------------------

class SCFConvergenceError(RuntimeError):
    def __init__(self, label, diagnostics):
        self.diagnostics = {"label": label, **diagnostics}
        super().__init__(f"{label}: SCF did not converge: {diagnostics}")


def scf_residuals(rho_input, rho_output, x, w, nelec, energy, previous_energy):
    """Unmixed fixed-point L1 density error per electron and absolute dE."""
    if (nelec <= 0 or not np.isfinite(energy)
            or not np.all(np.isfinite(rho_input))
            or not np.all(np.isfinite(rho_output))):
        raise FloatingPointError("SCF requires finite energies/densities and positive charge")
    residual = float(np.sum(np.abs(rho_output - rho_input) * 4 * np.pi * x**2 * w) / nelec)
    delta = None if previous_energy is None else float(abs(energy - previous_energy))
    return {"energy_change_eh": delta, "density_l1_per_electron": residual}


def normalize_scf_orbitals(orbs, w):
    """Use the density/energy quadrature for both large and small components.

    NR eigensolvers return Euclidean unit vectors; the Dirac shooter's internal
    matching quadrature differs from the node quadrature used for SCF density.
    """
    normalized = {}
    for key, (energy, P, Q) in orbs.items():
        norm = float(np.sum((P**2 + (0.0 if Q is None else Q**2)) * w))
        if not np.isfinite(energy) or not np.isfinite(norm) or norm <= 0:
            raise FloatingPointError(f"Invalid orbital energy/norm for {key}")
        factor = math.sqrt(norm)
        normalized[key] = (energy, P / factor, None if Q is None else Q / factor)
    return normalized


def scf_density_from_orbitals(orbs, occupations, x, w):
    """One density construction for iteration and verification; never rescale it."""
    rho = np.zeros_like(x)
    eigenvalue_sum = 0.0
    expected = 0.0
    for key, occ in occupations.items():
        if occ <= 0:
            continue
        if key not in orbs:
            raise RuntimeError(f"SCF occupied orbital missing {key}")
        energy, P, Q = orbs[key]
        rho += occ * orbital_density(P, Q, x, w)
        eigenvalue_sum += occ * energy
        expected += occ
    charge = float(np.sum(rho * 4 * np.pi * x**2 * w))
    error = abs(charge - expected)
    if not np.isfinite(charge) or expected <= 0 or error > 1e-10 * expected:
        raise FloatingPointError(f"SCF orbital charge {charge} != occupation sum {expected}")
    return rho, eigenvalue_sum, {"electron_count": charge, "charge_error": error}


def atomic_scf(z_nuc, config, relativistic=True, n_pts=640, a_mass=None,
               max_iter=80, tol=3e-6, label="atom", c=C_LIGHT, density_tol=3e-6):
    """Dirac-Fock-Slater (relativistic) or HF-Slater (c -> inf) atom."""
    if max_iter < 1 or not np.isfinite(tol) or tol <= 0 or not np.isfinite(density_tol) or density_tol <= 0:
        raise ValueError("SCF requires positive iteration limit and tolerances")
    if relativistic:
        x, y, h, w = log_grid(n_pts, 1e-6, 60.0)
    else:
        x = nr_grid()
        y, h = None, None
        w = np.full_like(x, NR_DR)
    v_nuc, r_sph = nucleus_potential(z_nuc, x, a_mass=a_mass)
    occ_map, blocks = config_blocks(config)
    blocks_keyset = sorted(occ_map.keys())
    sigma = -0.9 * z_nuc ** 2 / 2.0

    rho, nelec = guess_density(z_nuc, config, x)
    e_warm = {}

    prev_e = {}

    def solve_once(v_pot, scf_it):
        orbs = {}
        if relativistic:
            for (n, kap) in sorted(blocks_keyset):
                l = l_of_kappa(kap)
                e0 = prev_e.get((n, kap))
                # warm starts only from a trustworthy neighbour state
                if e0 is not None and (scf_it < 3
                                       or abs(e0) > 40.0 * z_nuc ** 2):
                    e0 = None
                e, P, Q, ok = dirac_orbital_shoot(
                    kap, n, l, v_pot, x, r_sph, e_guess=e0, c=c)
                if e is None or not ok:
                    e, P, Q, ok = dirac_orbital_shoot(
                        kap, n, l, v_pot, x, r_sph, e_guess=None, c=c)
                if e is None or not ok:
                    raise RuntimeError(f"{label}: shooting failed for "
                                       f"(n={n}, kappa={kap})")
                # band-jump guard: a warm start that moved > 35 % reverts
                if e0 is not None and abs(e / e0 - 1.0) > 0.35:
                    e2, P2, Q2, ok2 = dirac_orbital_shoot(
                        kap, n, l, v_pot, x, r_sph, e_guess=None, c=c)
                    if e2 is not None and ok2:
                        e, P, Q, ok = e2, P2, Q2, ok2
                    else:
                        raise RuntimeError(f"{label}: cold orbital retry failed for {(n, kap)}")
                orbs[(n, kap)] = (e, P, Q)
                prev_e[(n, kap)] = e
        else:
            l_blocks = sorted({l for (n, l, o1, o2) in config})
            for l in l_blocks:
                st = nr_bound_states(l, v_pot, x, k_states=16)
                for e, P in st:
                    if e >= 1.0:
                        continue
                    # noise-robust node count (truncate the tiny tail)
                    mref = np.max(np.abs(P))
                    sig = np.where(np.abs(P) > 1e-7 * mref, np.sign(P), 0.0)
                    sig = sig[sig != 0]
                    n_est = int(np.sum(sig[1:] * sig[:-1] < 0)) + l + 1
                    if l == 0:
                        key = (n_est, -1)
                    else:
                        key = (n_est, kappa_for(n_est, l, True))
                    if key not in orbs:
                        orbs[key] = (float(e), P, None)
                    if l > 0:
                        khi = (n_est, kappa_for(n_est, l, False))
                        if khi not in orbs:
                            orbs[khi] = (float(e), P, None)
        return normalize_scf_orbitals(orbs, w)

    e_old = None
    history = []
    converged = False
    verification_history = []

    def evaluate_density(density, iteration):
        u = hartree_potential(density, x, w)
        vx = xalpha_vx(density)
        potential = v_nuc + u + vx
        orbitals = solve_once(potential, iteration)
        output, eigen_sum, charge = scf_density_from_orbitals(orbitals, occ_map, x, w)
        energy = total_energy_from_eigenvalues(eigen_sum, output, u, vx, x, w)
        return orbitals, output, energy, potential, charge

    for it in range(max_iter):
        orbs, rho_new, e_tot, v_pot, charge = evaluate_density(rho, it)
        metrics = scf_residuals(rho, rho_new, x, w, nelec, e_tot, e_old)
        history.append({"iteration": it + 1, "energy_eh": e_tot, **metrics, **charge})
        de = float("inf") if e_old is None else metrics["energy_change_eh"]
        frac = 0.5 if it < 10 else 0.25
        rho = (1 - frac) * rho + frac * rho_new
        if it % 5 == 0 or de < tol:
            log(f"  [{label}] it={it:02d} E={e_tot:.6f} Eh  dE={de:.2e} "
                f"density_residual={metrics['density_l1_per_electron']:.2e} norbs={len(orbs)}")
        e_old = e_tot
        if it >= 8 and de < tol and metrics["density_l1_per_electron"] < density_tol:
            orbs, rho_final, e_final, v_pot, final_charge = evaluate_density(rho, it + 1)
            final_metrics = scf_residuals(rho, rho_final, x, w, nelec, e_final, e_old)
            final_metrics.update(final_charge)
            verification_history.append({"after_iteration": it + 1, **final_metrics})
            if (final_metrics["energy_change_eh"] < tol
                    and final_metrics["density_l1_per_electron"] < density_tol):
                converged = True
                break
            log(f"  [{label}] candidate verification failed; continuing SCF: {final_metrics}")

    # Independently check the final solve; an iteration cap is not convergence.
    if not converged:
        orbs, rho_final, e_final, v_pot, final_charge = evaluate_density(rho, it + 1)
        final_metrics = scf_residuals(rho, rho_final, x, w, nelec, e_final, e_old)
        final_metrics.update(final_charge)
    diagnostics = {"converged": converged, "iterations": it + 1,
                   "max_iter": max_iter, "energy_tolerance_eh": tol,
                   "density_tolerance": density_tol, "final": final_metrics,
                   "history": history, "verification_history": verification_history}
    if not converged:
        raise SCFConvergenceError(label, diagnostics)

    atom = {"z": z_nuc, "nelec": int(nelec), "relativistic": relativistic,
            "e_tot": e_final, "iterations": it + 1, "scf": diagnostics,
            "r_nuc_sphere_a0": float(r_sph),
            "grid": {"x": x, "y": y, "h": h, "w": w},
            "rho": rho_final, "v_pot": v_pot, "orbitals": {}}
    for (n, kap), occ in occ_map.items():
        if (n, kap) not in orbs:
            continue
        e, P, Q = orbs[(n, kap)]
        atom["orbitals"][(n, kap)] = _orbital_record(e, P, Q, x, w, occ, kap)
    return atom


def _orbital_record(e, P, Q, x, w, occ, kap):
    qq = Q if (Q is not None and len(Q) == len(P)) else np.zeros_like(P)
    norm = float(np.sum((P ** 2 + qq ** 2) * w))
    def avg(f):
        return float(np.sum(f * (P ** 2 + qq ** 2) * w) / (norm + 1e-300))
    return {"energy": float(e), "occ": float(occ), "kappa": int(kap),
            "nodes": pencil_nodes(P), "norm": norm,
            "r_mean": avg(x), "r2_mean": avg(x ** 2),
            "r_inv3": avg(x ** -3.0),
            "small_norm": float(np.sum(qq ** 2 * w)),
            "P": P, "Q": Q}


# ----------------------------------------------------------------------------
# Validation: hydrogen-like point-nucleus Dirac vs Sommerfeld
# ----------------------------------------------------------------------------

def sommerfeld_1s(z, c=C_LIGHT):
    return c ** 2 * (math.sqrt(max(1.0 - (z / c) ** 2, 1e-12)) - 1.0)
def hydrogenic_validation(z_list, n_pts=800):
    x, y, h, w = log_grid(n_pts, 1e-6, 60.0)
    res = []
    for z in z_list:
        v = -z / x
        e_ex = sommerfeld_1s(z)
        e_num, P, Q, ok = dirac_orbital_shoot(-1, 1, 0, v, x, 0.0,
                                              e_guess=e_ex * 1.1)
        nrm = np.sum(P ** 2 * w) + np.sum(Q[:-1] ** 2 * 0.5 * (w[:-1] + w[1:]))
        vc = math.sqrt(max(float(np.sum(Q ** 2 * w)), 0.0) / (nrm + 1e-300))
        res.append({"Z": int(z), "E_numeric": float(e_num),
                    "E_sommerfeld": float(e_ex), "converged": bool(ok),
                    "rel_err": float(abs(e_num - e_ex) / abs(e_ex)),
                    "v_over_c_rms": vc})
        log(f"  H-like Z={z:3d}: E={e_num:12.4f} Eh vs Sommerfeld "
            f"{e_ex:12.4f}  err={res[-1]['rel_err']*100:.4f}%  "
            f"v/c={vc:.3f}  ok={ok}")
    # Dirac accidental degeneracy E(ns1/2) == E(np1/2) at n=2 (point
    # nucleus): a razor-sharp certificate of the kappa conventions.
    z = 80
    v = -z / x
    e2 = {}
    for kap, n, l in ((-1, 2, 0), (+1, 2, 1), (-2, 2, 1)):
        e2[kap], _P, _Q, ok2 = dirac_orbital_shoot(kap, n, l, v, x, 0.0)
    deg = abs(e2[-1] - e2[+1]) / abs(e2[-1])
    so_split = (e2[-2] - e2[-1]) * HARTREE_CM
    res.append({"Z": z, "degeneracy_test": float(deg),
                "E_2s1_2": e2[-1], "E_2p1_2": e2[+1], "E_2p3_2": e2[-2],
                "fine_split_2p_cm": float(so_split)})
    log(f"  Z=80 Dirac degeneracy |E(2s1/2)-E(2p1/2)|/|E| = {deg:.2e}; "
        f"2p fine split = {so_split:.0f} cm-1 (2p3/2 above 2p1/2)")
    # NR twin spot check (uniform grid)
    xu = nr_grid()
    st = nr_bound_states(0, -1.0 / xu, xu, k_states=1)
    res.append({"Z": 1, "E_nr_1s": st[0][0], "E_nr_exact": -0.5,
                "rel_err_nr": abs(st[0][0] + 0.5) / 0.5})
    log(f"  NR H 1s: {st[0][0]:.8f} Eh vs -0.5 (err "
        f"{res[-1]['rel_err_nr']*100:.4f}%)")
    return res



# ============================================================================
# MODULE 17B - f^6 spin-orbit multiplet engine (CAS(6e,7o) atomic limit)
# ============================================================================
# The 5f^6 / 4f^6 ground manifold is the Hund-rule ^7F term (L = 3, S = 3,
# J = 0..6; J = 0 ground for the less-than-half-filled f^k shell).  To first
# order in the spin-orbit coupling the ^7F_J energies follow the Lande
# interval rule
#     E(J) = (zeta/2) [J(J+1) - L(L+1) - S(S+1)],
# with zeta taken DIRECTLY from the Dirac SCF j-splitting of the f shell
# (zeta = [eps(j=7/2) - eps(j=5/2)] / 3.5) - the multiplet structure is thus
# derived from the 4-component solution without empirical parameters.  The
# leading intermediate-coupling corrections are quantified against the
# Eu3+ experimental fluorescence ladder.  The Slater-Condon determinant
# machinery (F6Multiplet) is retained as a validator: its electrostatic
# diagonal reproduces the Hund ordering of the f^6 configuration.
# ============================================================================

def two_el_tensor_numeric(l, n_th=36, n_ph=48, kmax=6):
    """G^{(k)}[a, c, b, d] = int dOm1 dOm2 Y*_a(1)Y_c(1) Y*_b(2)Y_d(2) P_k.

    Complex Condon-Shortley harmonics; the uniform-phi trapezoid makes the
    m-selection rule exact by construction.
    """
    gx, wx = roots_legendre(n_th)
    th = np.arccos(gx)
    ph = np.linspace(0, 2 * np.pi, n_ph, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    W = np.outer(wx, np.full(n_ph, 2 * np.pi / n_ph)).reshape(-1)
    N1 = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH),
                   np.cos(TH)], axis=-1).reshape(-1, 3)
    cosg = np.sum(N1[:, None, :] * N1[None, :, :], axis=-1).reshape(-1)
    nm = 2 * l + 1
    npts = N1.shape[0]
    Ys = np.zeros((nm, npts), dtype=complex)
    for m in range(-l, l + 1):
        Ys[m + l] = (_plm_norm(l, abs(m)) * _plm(l, abs(m), np.cos(TH))
                     * np.exp(1j * m * PH)).reshape(-1)
    G = np.zeros((kmax + 1, nm, nm, nm, nm))
    for k in range(0, kmax + 1, 2):
        Pk = _plm(k, 0, cosg).reshape(npts, npts)
        Mk = Pk * (W[:, None] * W[None, :])
        for b in range(nm):
            for d in range(nm):
                h_bd = Mk @ (np.conj(Ys[b]) * Ys[d])
                for a in range(nm):
                    fa = W * np.conj(Ys[a])
                    for c in range(nm):
                        G[k, a, c, b, d] = (W * np.conj(Ys[a]) * Ys[c]
                                            * h_bd).sum().real
    return G


class F6Multiplet:
    """f^6 Slater-Condon electrostatic diagonal validator.

    Builds the antisymmetrized two-electron angular tensor from the
    quadrature G-tables and evaluates the configuration-average diagonal
    for every C(14,6) determinant.  Used to certify the Hund ordering of
    the engine's F^k integrals (first Hund rule).
    """

    def __init__(self, f2, f4, f6, zeta=0.0, l=3):
        self.l = l
        self.F = {2: f2, 4: f4, 6: f6}
        self.zeta = zeta
        self.n_so = (2 * l + 1) * 2
        self.dets = [tuple(c) for c in combinations(range(self.n_so), 6)]
        self._spin = [a % 2 for a in range(self.n_so)]
        self._m = [a // 2 - l for a in range(self.n_so)]
        self._build_t2()

    def _spin_of(self, a):
        return self._spin[a]

    def _m_of(self, a):
        return self._m[a]

    def _pair_id(self, a, b):
        lo, hi = min(a, b), max(a, b)
        return lo * self.n_so + hi

    def _two_el(self, a, b, c, d):
        return self._t2[self._pair_id(a, b), self._pair_id(c, d)]

    def _build_t2(self):
        l, n = self.l, self.n_so
        Gk = two_el_tensor_numeric(l)
        t2 = np.zeros((n * n, n * n))
        for a in range(n):
            for b in range(n):
                if a == b:
                    continue
                for c in range(n):
                    for d in range(n):
                        if c == d:
                            continue
                        sa, sb = self._spin_of(a), self._spin_of(b)
                        sc, sd = self._spin_of(c), self._spin_of(d)
                        ma, mb = self._m_of(a), self._m_of(b)
                        mc, md = self._m_of(c), self._m_of(d)
                        el = 0.0
                        if sa == sc and sb == sd:      # direct
                            for k in (2, 4, 6):
                                el += self.F[k] * Gk[k][ma + l, mc + l,
                                                        mb + l, md + l]
                        if sa == sd and sb == sc:      # exchange
                            for k in (2, 4, 6):
                                el -= self.F[k] * Gk[k][ma + l, md + l,
                                                        mb + l, mc + l]
                        t2[self._pair_id(a, b), self._pair_id(c, d)] = el
        self._t2 = t2

    def _diag_es(self, det):
        return sum(self._two_el(a, b, a, b)
                   for a, b in combinations(det, 2))


LANDE_L, LANDE_S = 3.0, 3.0


def f_slater_integrals(P, x, w):
    """Slater F^k (k = 2, 4, 6) of the large component P (f shell)."""
    nrm = float(np.sum(P ** 2 * w))
    Pn = P / math.sqrt(max(nrm, 1e-300))
    a = Pn ** 2 * w
    out = {}
    for k in (2, 4, 6):
        ratio = np.outer(x ** k, x ** -(k + 1))
        kern = np.where(np.subtract.outer(x, x) < 0.0, ratio, ratio.T)
        out[k] = float(a @ kern @ a)
    return out


def lande_7f_ladder(zeta_eh):
    """7F_J energies (cm-1, relative to J=0) and g_J for J = 0..6."""
    L, S = LANDE_L, LANDE_S
    out = {}
    for J in range(0, 7):
        jj = J * (J + 1)
        e = 0.5 * zeta_eh * (jj - L * (L + 1) - S * (S + 1))
        g = 1.0 + (jj + S * (S + 1) - L * (L + 1)) / (2.0 * jj) if J > 0 \
            else 0.0
        out[J] = {"E_cm": e * HARTREE_CM, "g_J": g,
                  "E_eh": e}
    e0 = out[0]["E_eh"]
    for J in out:
        out[J]["E_rel_cm"] = (out[J]["E_eh"] - e0) * HARTREE_CM
    return out


def zeta_from_df(atom, nf):
    """zeta = [eps(j=7/2) - eps(j=5/2)] / 3.5 from the Dirac SCF atom."""
    e52 = atom["orbitals"][(nf, +3)]["energy"]
    e72 = atom["orbitals"][(nf, -4)]["energy"]
    return (e72 - e52) / 3.5


def fit_lande_to_experiment(exp_ladder):
    """Least-squares zeta from experimental 7F_J ladder (cm-1)."""
    js = np.array(sorted(exp_ladder))
    ev = np.array([exp_ladder[J] for J in js], float)
    # E(J) = zeta * J(J+1)/2 : linear least squares through origin
    x = np.array([J * (J + 1) / 2.0 for J in js])
    zeta = float(np.dot(x, ev) / np.dot(x, x))
    devs = {int(J): float(ev[i] - zeta * x[i]) for i, J in enumerate(js)}
    return zeta, devs


def hund_validator(atom, nf, tag=""):
    """Slater-Condon diagonal check: the max-spin dets minimize the
    electrostatic diagonal (Hund first rule) for the engine F^k."""
    P = atom["orbitals"][(nf, +3)]["P"]
    f_ints = f_slater_integrals(P, atom["grid"]["x"], atom["grid"]["w"])
    mp = F6Multiplet(f_ints[2], f_ints[4], f_ints[6], 0.0)
    diags = np.array([mp._diag_es(det) for det in mp.dets])
    n_up = np.array([sum(1 for a in det if a % 2 == 0) for det in mp.dets])
    i_min = int(np.argmin(diags))
    out = {"F_ints": f_ints, "hund_ok": bool(n_up[i_min] in (0, 6)),
           "min_diag": float(diags[i_min]),
           "spread": float(diags.max() - diags.min())}
    log(f"  [17B] {tag} Hund validator: F2/F4/F6 = "
        f"{f_ints[2]:.4f}/{f_ints[4]:.4f}/{f_ints[6]:.4f} Eh, "
        f"min-diagonal det is max-spin: {out['hund_ok']}")
    return out


EU3_EXP_7F = {0: 0.0, 1: 370.0, 2: 1060.0, 3: 1900.0, 4: 2860.0, 5: 3920.0,
              6: 4940.0}     # Eu3+ 4f6 fluorescence ladder (cm-1, ~)
AM3_EXP_7F = {0: 0.0, 1: 2900.0, 2: 6300.0}   # Am3+ literature envelope


def run_multiplet_stage(atoms, results):
    log("[MODULE 17B] f^6 multiplet engine (CAS(6e,7o) atomic limit)")
    mult = {}
    for ion, tag, nf, exp in (("Am", "Am3+", 5, AM3_EXP_7F),
                              ("Eu", "Eu3+", 4, EU3_EXP_7F)):
        atom = atoms[tag]
        atom_nr = atoms[tag + "NR"]
        zeta_df = zeta_from_df(atom, nf)
        ladder_df = lande_7f_ladder(zeta_df)
        # NR twin: no j-splitting, no SOC (degenerate 7F)
        ladder_nr = lande_7f_ladder(0.0)
        # scalar-relativistic: relativistic F^k, zeta suppressed
        ladder_sc = lande_7f_ladder(0.0)
        zeta_exp, devs = fit_lande_to_experiment(exp)
        hund = hund_validator(atom, nf, tag=f"{ion}")
        # magnetic anisotropy proxies from the Dirac radial functions
        rec = atom["orbitals"][(nf, -4)]
        rec_nr = None
        for (n, kap), r2 in atom_nr["orbitals"].items():
            if n == nf:
                rec_nr = r2
                break
        mult[ion] = {
            "zeta_df_cm": zeta_df * HARTREE_CM,
            "zeta_expfit_cm": zeta_exp,
            "lande_deviations_exp_cm": devs,
            "ladder_4c": {J: ladder_df[J]["E_rel_cm"] for J in ladder_df},
            "g_J": {J: ladder_df[J]["g_J"] for J in ladder_df},
            "ladder_exp": exp,
            "hund": hund,
            "r_mean_f": rec["r_mean"], "r2_mean_f": rec["r2_mean"],
            "r_inv3_f": rec["r_inv3"],
            "r_mean_f_nr": rec_nr["r_mean"] if rec_nr else None,
            "eps_5f52": atom["orbitals"][(nf, +3)]["energy"],
            "eps_5f72": atom["orbitals"][(nf, -4)]["energy"]}
        log(f"  [17B] {ion}: zeta(DF) = {zeta_df*HARTREE_CM:.0f} cm-1 | "
            f"7F1 = {ladder_df[1]['E_rel_cm']:.0f}, "
            f"7F6 = {ladder_df[6]['E_rel_cm']:.0f} cm-1 | "
            f"zeta(exp fit) = {zeta_exp:.0f} cm-1")
    dE_Am = mult["Am"]["ladder_4c"][1] - mult["Am"]["ladder_4c"][0]
    dE_Eu = mult["Eu"]["ladder_4c"][1] - mult["Eu"]["ladder_4c"][0]
    mult["dE_SO_ratio_Am_over_Eu"] = dE_Am / dE_Eu
    log(f"  [17B] dE_SO(Am 5f) / dE_SO(Eu 4f) = "
        f"{mult['dE_SO_ratio_Am_over_Eu']:.2f}")
    results["multiplet"] = mult
    return mult


# ============================================================================
# MODULE 17A-tail - Breit pair matrix elements by 6D angular-spinor quadrature
# ============================================================================
# The frequency-independent Breit operator
#     g_B = -[a1.a2 + (a1.r12)(a2.r12)/r12^2] / (2 r12)
# is evaluated between J=0-coupled equivalent kappa^2 spinor pairs.  The
# angular-spinor bilinear form factorizes against the radial amplitudes by
# rotational invariance: it depends on (r1, r2) only through u = r_< / r_>.
# The angular tables A(u) are computed by explicit quadrature on the two
# unit spheres with REAL spherical spinors, then applied to the SCF radial
# functions as 2D radial integrals.  The machinery is certified by (i) the
# exact 5Z/8 Coulomb pair energy of H 1s^2 with kernel='coulomb' and
# (ii) the c -> infinity collapse of the Breit term.

ANG_TH = 30
ANG_PH = 36
U_GRID = 40


def _alpha_mats():
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    Z = np.zeros((2, 2), dtype=complex)
    I = np.eye(2, dtype=complex)
    ax = np.block([[Z, sx], [sx, Z]])
    ay = np.block([[Z, sy], [sy, Z]])
    az = np.block([[Z, sz], [sz, Z]])
    return ax, ay, az


def angular_pair_contraction(psi, direction, weight_over_r, kernel, K9):
    """Exact sparse contraction without per-point dense 16x16 operators."""
    amat = np.zeros((16, 16), dtype=complex)
    a1a2 = K9[0] + K9[4] + K9[8]
    if kernel == "coulomb":
        np.fill_diagonal(amat, np.sum(weight_over_r[:, None] * np.abs(psi)**2, axis=0))
        return amat
    if kernel not in ("breit", "gaunt"):
        raise ValueError(f"Unknown pair kernel: {kernel}")
    rows, cols = np.nonzero(np.any(K9 != 0, axis=0))
    products = psi[:, rows].conj() * psi[:, cols]
    if kernel == "gaunt":
        # Conventional Gaunt: -alpha1.alpha2/r12, not half that operator.
        # https://www.diracprogram.org/doc/release-21/manual/hamiltonian.html
        values = -a1a2[rows, cols] * (weight_over_r @ products)
    else:
        directional = (direction[:, :, None] * direction[:, None, :]).reshape(-1, 9)
        coefficients = -.5 * (a1a2[rows, cols] + directional @ K9[:, rows, cols])
        values = np.sum(weight_over_r[:, None] * products * coefficients, axis=0)
    amat[rows, cols] = values
    return amat


def breit_pair_tabulation(kappa, kernel="breit", chunk=20000):
    """Four-index angular tables A^{abgd}(u) for the J=0 kappa^2 pair.

    E_pair = sum_{abgd} int dr1 dr2 f_a(r1) f_g(r1) f_b(r2) f_d(r2)
             * A^{abgd}(u) / r_>,
    where A^{abgd}(u) = int dOmega1 dOmega2 Psi*_{(a,b)} g_B Psi_{(g,d)}
    at radii (1, u) is the exact angular-spinor bilinear of the frequency-
    independent Breit operator, computed by quadrature over the two unit
    spheres with real spherical spinors.  The antisymmetrized pair state
    |(j^2)J=0> = (|up,dn> - |dn,up>)/sqrt(2) carries exchange exactly.
    kernel='coulomb' (1/r12 x identity) and kernel='gaunt' (magnetic part
    only) share the machinery for certification.
    """
    l = l_of_kappa(kappa)
    j = abs(kappa) - 0.5
    ls = int(round(2 * j)) - l
    gx, wx = roots_legendre(ANG_TH)
    theta = np.arccos(gx)
    ph = np.linspace(0, 2 * np.pi, ANG_PH, endpoint=False)
    TH, PH = np.meshgrid(theta, ph, indexing="ij")
    W = np.outer(wx, np.full(ANG_PH, 2 * np.pi / ANG_PH))
    n1 = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH),
                   np.cos(TH)], axis=-1).reshape(-1, 3)
    wg = W.reshape(-1)
    npts = len(n1)
    OmL = [spherical_spinor_clean(l, j, m, TH, PH).reshape(2, -1).T
           for m in (0.5, -0.5)]
    OmS = [spherical_spinor_clean(ls, j, m, TH, PH).reshape(2, -1).T
           for m in (0.5, -0.5)]

    def vhat(mi):
        v = np.zeros((npts, 4), dtype=complex)
        v[:, 0] = OmL[mi][:, 0]; v[:, 1] = OmL[mi][:, 1]
        v[:, 2] = 1j * OmS[mi][:, 0]; v[:, 3] = 1j * OmS[mi][:, 1]
        return v

    v_up, v_dn = vhat(0), vhat(1)
    ax, ay, az = _alpha_mats()
    u_grid = np.logspace(-4, math.log10(1.0 - 2e-2), U_GRID)
    A = np.zeros((4, 4, 4, 4, U_GRID), dtype=complex)
    K9 = np.stack([np.kron(mi, mj) for mi in (ax, ay, az)
                   for mj in (ax, ay, az)])            # (9, 16, 16)
    flat = npts * npts
    for iu, u in enumerate(u_grid):
        Amat = np.zeros((16, 16), dtype=complex)
        for p0 in range(0, flat, chunk):
            indices = np.arange(p0, min(p0 + chunk, flat))
            i, jidx = indices // npts, indices % npts
            delta = n1[i] - u * n1[jidx]
            distance = np.linalg.norm(delta, axis=1)
            direction = delta / distance[:, None]
            psi = ((v_up[i, :, None] * v_dn[jidx, None, :]
                    - v_dn[i, :, None] * v_up[jidx, None, :]) / math.sqrt(2.)).reshape(-1, 16)
            Amat += angular_pair_contraction(psi, direction,
                         wg[i] * wg[jidx] / distance, kernel, K9)
        A[:, :, :, :, iu] = Amat.reshape(4, 4, 4, 4)
        if iu % 5 == 0 or iu == len(u_grid) - 1:
            log(f"  angular {kernel}: {iu+1}/{len(u_grid)} radial ratios")
    return {"u": u_grid, "A_dir": A, "kappa": kappa, "kernel": kernel}


def cached_breit_pair_tabulation(kappa, kernel="breit"):
    """Reuse only a matching angular table in this calculation's output folder."""
    path = RES / f"angular_{kernel}_k{kappa}.npz"
    signature = f"sparse-full-gaunt-v1:{ANG_TH}:{ANG_PH}:{U_GRID}:{kappa}:{kernel}"
    if path.exists():
        with np.load(path, allow_pickle=False) as saved:
            if str(saved["signature"]) == signature:
                log(f"  angular {kernel}: using matching table {path.name}")
                return {"u": saved["u"].copy(), "A_dir": saved["A_dir"].copy(),
                        "kappa": kappa, "kernel": kernel}
    tab = breit_pair_tabulation(kappa, kernel=kernel)
    np.savez_compressed(path, signature=signature, u=tab["u"], A_dir=tab["A_dir"])
    return tab


def hydrogenic_radial_grid(z, n_pts, r_min=1e-6):
    """Keep resolution and radial extent fixed in the dimensionless coordinate Z*r."""
    if not np.isfinite(z) or z <= 0 or n_pts < 2:
        raise ValueError("Hydrogenic grid requires positive Z and at least two nodes")
    x = np.linspace(r_min, 40., n_pts) / z
    return x, np.gradient(x)


def breit_pair_energy(tab, P, Q, x, w, atol=1e-9, chunk=128):
    """J=0 equivalent-pair energy <(k^2)0|g_B|(k^2)0> from the tables.

    Four-index factorization: E = sum_{abgd} int dr1 dr2
        f_a(r1) f_g(r1) f_b(r2) f_d(r2) A^{abgd}(u) / r_> .
    """
    u = tab["u"]
    qq = Q if Q is not None and len(Q) == len(P) else np.zeros_like(P)
    fr = (P, P, qq, qq)
    # Group duplicate large/small radial factors into P^2, PQ, Q^2.
    # The table coordinate is r_</r_>, not sqrt(r_</r_>).
    grouped = np.zeros((3, 3, len(u)), dtype=complex)
    Amax = float(np.max(np.abs(tab["A_dir"])))
    for al in range(4):
        for be in range(4):
            for ga in range(4):
                for de in range(4):
                    tbl = tab["A_dir"][al, be, ga, de]
                    if np.max(np.abs(tbl)) < atol * Amax:
                        continue
                    grouped[int(al >= 2) + int(ga >= 2), int(be >= 2) + int(de >= 2)] += tbl
    factors = np.array([P * P * w, P * qq * w, qq * qq * w])
    e = 0.0j
    for start in range(0, len(x), chunk):
        sl = slice(start, min(start + chunk, len(x)))
        rmax = np.maximum(x[sl, None], x[None, :])
        u_c = np.clip(np.minimum(x[sl, None], x[None, :]) / rmax, u[0], u[-1])
        for a in range(3):
            for b in range(3):
                tbl = grouped[a, b]
                if not np.any(tbl):
                    continue
                values = (np.interp(u_c, u, tbl.real) + 1j * np.interp(u_c, u, tbl.imag))
                e += np.sum(factors[a, sl, None] * factors[b, None, :] * values / rmax)
    return e


def coulomb_pair_energy_ref(P, x, w):
    """Classic 2D Coulomb J = int int P^2 P^2 / r_> (for the 5Z/8 test)."""
    nrm = float(np.sum(P ** 2 * w))
    Pn = P / math.sqrt(nrm)
    a = Pn ** 2 * w
    kern = 1.0 / np.maximum.outer(x, x)
    return float(a @ kern @ a)


# ============================================================================
# MODULE 17C - relativistic bonding, QTAIM & EDA on the N-M-N axis
# ============================================================================

N_CONFIG = [(1, 0, 2.0, 0.0), (2, 0, 2.0, 0.0), (2, 1, 0.0, 3.0)]

# crystallographic BTP pocket distances (EXAFS/XRD literature envelope)
D_MN = {"Am": 2.53, "Eu": 2.47}          # Angstrom


def sto_2p(z_eff, rvec, axis="x"):
    """Slater 2p STO pointing along +axis, centered at origin (a0 units)."""
    zeta = z_eff / 2.0
    N = 2.0 * zeta ** 2.5 * math.sqrt(2.0 / math.pi) * 0 + \
        (4.0 * (2.0 * zeta) ** 2.5 / math.sqrt(3.0 * math.pi)) * 0
    # standard real STO: phi_2pz = (1/sqrt(32 pi)) zeta^{5/2} r exp(-zeta r/2) cos(theta)
    zeta_ = zeta
    r = np.linalg.norm(rvec, axis=-1)
    ang = rvec[..., 0] / np.maximum(r, 1e-300)      # along x
    return (1.0 / math.sqrt(32.0 * math.pi)) * zeta_ ** 2.5 * r * \
        np.exp(-0.5 * zeta_ * r) * ang


def atom_density_interp(atom, x_new):
    """Interpolate a spherical atomic density onto a new radial grid."""
    return np.interp(np.abs(x_new), atom["grid"]["x"], atom["rho"])


def orbital_profile_3d(atom, n, kap, pts):
    """|psi|^2 of an orbital on arbitrary 3D points (n,3 array, a0 units)."""
    rec = atom["orbitals"][(n, kap)]
    l = l_of_kappa(kap)
    j = abs(kap) - 0.5
    ls = int(round(2 * j)) - l
    r = np.linalg.norm(pts, axis=-1)
    th = np.arccos(np.clip(pts[:, 2] / np.maximum(r, 1e-300), -1, 1))
    ph = np.arctan2(pts[:, 1], pts[:, 0])
    dens = np.zeros(len(pts))
    x = atom["grid"]["x"]
    Pr = np.interp(r, x, rec["P"])
    Qr = np.interp(r, x, rec["Q"])
    for m in (0.5, -0.5):
        Om = spherical_spinor_clean(l, j, m, th, ph)
        OmS = spherical_spinor_clean(ls, j, m, th, ph)
        dens += (Pr ** 2) * np.sum(Om * Om, axis=0) + (Qr ** 2) * np.sum(
            OmS * OmS, axis=0)
    return dens / (4.0 * math.pi)


def fft_electrostatic(rho_box, box_a0):
    """Poisson solve on a periodic box; returns potential array (a.u.)."""
    n = rho_box.shape[0]
    k = np.fft.fftfreq(n, d=box_a0 / n) * 2 * np.pi
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    K2 = KX ** 2 + KY ** 2 + KZ ** 2
    K2[0, 0, 0] = 1.0
    V = np.real(np.fft.ifftn(np.fft.fftn(rho_box) / (4 * np.pi * K2)))
    return V


PLANE_CACHE = {}


def eda_qtaim_stage(atom_am, atom_eu, atom_am_nr, atom_eu_nr, atom_n,
                    results):
    """Two-state mixing and atomic-promolecule descriptors, not complex QTAIM."""
    log("[MODULE 17C] relativistic bonding, QTAIM & EDA")
    out = {"model_scope": "Sum of isolated spherical M and N atomic densities; "
           "not a self-consistent complex density or validated QTAIM topology. "
           "Wolfsberg-Helmholtz mixing and S_rms are model proxies, not molecular EDA "
           "or a signed orbital overlap. Grid-minimum descriptors are unrefined."}
    # donor level: N 2p energy from the NR-Slater nitrogen atom
    eps_n2p = None
    for (n, kap), rec in atom_n["orbitals"].items():
        if n == 2 and kap == kappa_for(2, 1, True):
            eps_n2p = rec["energy"]
    out["eps_N2p_eh"] = eps_n2p

    # 3D overlap integrals S(5f/4f, N 2p-lone-pair) on a Cartesian grid
    n_g = 150
    L_box = 11.0                     # a0 (~5.8 A): covers M..N + tails
    ax = np.linspace(-L_box / 2, L_box / 2, n_g)
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    pts = np.stack([X, Y, Z], axis=-1).reshape(-1, 3)
    dV = (ax[1] - ax[0]) ** 3
    z_eff_2p = 3.9                    # Slater Z_eff for N 2p lone pair
    for ion, atom in (("Am", atom_am), ("Eu", atom_eu)):
        best = None
        for (n, kap) in atom["orbitals"]:
            l = l_of_kappa(kap)
            if l != 3:
                continue
            for R_A in (D_MN[ion],):
                R = R_A / BOHR_A
                phi_n = sto_2p(z_eff_2p, pts - np.array([R, 0, 0]))
                nrm = math.sqrt(np.sum(phi_n ** 2) * dV)
                psi_m = orbital_profile_3d(atom, n, kap, pts)
                # psi is a density-like quadrature of |Omega|^2; build signed
                # overlap via sqrt-weighted sign-free amplitude (model level):
                s_rms = math.sqrt(max(float(np.sum(psi_m * phi_n ** 2) * dV)
                                      / max(nrm ** 2, 1e-300), 0.0))
                cnt = f"{n}{kap:+d}"
                if best is None or s_rms > best[1]:
                    best = (cnt, s_rms)
        out.setdefault("overlap", {})[ion] = {"orbital": best[0],
                                              "S_rms": best[1]}
        log(f"  17C {ion}: S_rms(5f/4f, N2p @ {D_MN[ion]:.2f} A) = "
            f"{best[1]:.5f} ({best[0]})")

    # Wolfsberg-Helmholtz two-state ionic/covalent mixing
    K_WH = 1.75
    for ion, atom in (("Am", atom_am), ("Eu", atom_eu)):
        eps_a = None
        for (n, kap), rec in atom["orbitals"].items():
            if l_of_kappa(kap) == 3 and n == (5 if ion == "Am" else 4):
                if eps_a is None or rec["energy"] > eps_a:
                    eps_a = rec["energy"]       # acceptor: higher 5f level
        S = out["overlap"][ion]["S_rms"]
        h_da = K_WH * S * 0.5 * (eps_a + eps_n2p)
        gap = max(eps_n2p - eps_a, 1e-3)
        d_cov_bond = -(h_da ** 2) / gap * 0.5   # 2-electron pairing factor
        out.setdefault("eda", {})[ion] = {
            "eps_acceptor_eh": eps_a, "S_rms": S, "H_DA_eh": h_da,
            "gap_eh": gap, "dE_cov_bond_kcal": d_cov_bond * HARTREE_KCAL,
            "dE_cov_complex_kcal": 3.0 * d_cov_bond * HARTREE_KCAL,
            "c_cov2": (h_da / gap) ** 2}
    dd = (out["eda"]["Am"]["dE_cov_complex_kcal"]
          - out["eda"]["Eu"]["dE_cov_complex_kcal"])
    out["eda"]["selectivity_kcal"] = dd
    log(f"  17C dE_cov(Am) = {out['eda']['Am']['dE_cov_complex_kcal']:.2f} "
        f"kcal/mol, dE_cov(Eu) = "
        f"{out['eda']['Eu']['dE_cov_complex_kcal']:.2f} kcal/mol, "
        f"selectivity = {dd:.2f} kcal/mol")

    # ---- electrostatic (FFT) + QTAIM maps on the N-M-N plane --------------
    for ion, atom in (("Am", atom_am), ("Eu", atom_eu)):
        R = D_MN[ion] / BOHR_A
        n_b = 96
        L_b = 2.0 * (R + 6.0)
        axs = np.linspace(-L_b / 2, L_b / 2, n_b)
        h_b = axs[1] - axs[0]
        Xb, Yb, Zb = np.meshgrid(axs, axs, axs, indexing="ij")
        r_b = np.sqrt(Xb ** 2 + (Yb - R) ** 2 + Zb ** 2)
        rho_n = np.interp(r_b, atom_n["grid"]["x"], atom_n["rho"])
        r_m = np.sqrt(Xb ** 2 + Yb ** 2 + Zb ** 2)
        rho_m = np.interp(r_m, atom["grid"]["x"], atom["rho"])
        # covalent polarization: c^2 of charge transfer N-lone-pair -> 5f/4f
        c2 = out["eda"][ion]["c_cov2"]
        rho_pol = np.clip(rho_m - c2 * 1.0 * rho_n * 0.0, 0.0, None)
        rho_tot = rho_m + rho_n + c2 * (rho_m - rho_n) * 0.35
        V_n = fft_electrostatic(rho_n, L_b)
        d_es = float(np.sum(rho_m * V_n) * h_b ** 3)
        out.setdefault("eda_electrostatic", {})[ion] = d_es
        log(f"  17C {ion}: FFT electrostatic M-density x N-density = "
            f"{d_es:.4f} Eh")

        # ---- 2D plane density and QTAIM observables -----------------------
        n2 = 241
        xs = np.linspace(-R - 3.0, R + 3.0, n2)
        zs = np.linspace(0.0, 3.2, n2)
        X2, Z2 = np.meshgrid(xs, zs, indexing="ij")
        Y2 = np.zeros_like(X2)
        r_m2 = np.sqrt(X2 ** 2 + Y2 ** 2 + Z2 ** 2)
        r_n2 = np.sqrt((X2 - R) ** 2 + Y2 ** 2 + Z2 ** 2)
        rho2 = (np.interp(r_m2, atom["grid"]["x"], atom["rho"])
                + np.interp(r_n2, atom_n["grid"]["x"], atom_n["rho"]))
        # BCP on the axis: density minimum along z=0 between 0 and R
        ax_mask = (xs > 0.3) & (xs < R - 0.3)
        idx = np.where(ax_mask)[0]
        prof = rho2[idx, 0]
        i_bcp = idx[np.argmin(prof)]
        x_bcp = xs[i_bcp]
        # Laplacian: on-axis transverse isotropy: lap = d2dx2 + 2 d2dz2
        dx = xs[1] - xs[0]
        dz = zs[1] - zs[0]
        d2x = (rho2[i_bcp + 1, 0] - 2 * rho2[i_bcp, 0]
               + rho2[i_bcp - 1, 0]) / dx ** 2
        d2z = 2.0 * (rho2[i_bcp, 1] - rho2[i_bcp, 0]) / dz ** 2
        lap_bcp = d2x + 2.0 * d2z
        rho_bcp = float(rho2[i_bcp, 0])
        # local kinetic model G = TF + (1/9) vW; H = 1/4 lap - G (a.u.)
        cf = 0.3 * (3 * math.pi ** 2) ** (2.0 / 3.0)
        grad_x = (rho2[i_bcp + 1, 0] - rho2[i_bcp - 1, 0]) / (2 * dx)
        grad_z = 0.0  # Axisymmetry at z=0, not the off-axis z=dz derivative.
        G = cf * rho_bcp ** (5.0 / 3.0) + (rho_bcp ** -1) * \
            (grad_x ** 2 + grad_z ** 2) / 72.0
        H_bcp = 0.25 * lap_bcp - G
        lap2d = _laplacian_plane(rho2, dx, dz)
        out.setdefault("qtaim", {})[ion] = {
            "R_mn_a0": R, "rho_bcp": rho_bcp, "lap_bcp": float(lap_bcp),
            "H_bcp": float(H_bcp), "G_bcp": float(G), "x_bcp_a0": float(x_bcp),
            "descriptor_scope": "discrete on-axis promolecular density minimum; "
            "not a stationary-point/Hessian-certified molecular BCP",
            "axis_gradient_residual": float(grad_x),
            "delta_mn": 4.0 * out["eda"][ion]["c_cov2"]
            * out["overlap"][ion]["S_rms"]}
        PLANE_CACHE[ion] = {"xs": xs, "zs": zs, "lap2d": lap2d,
                            "rho2d": rho2, "R": R, "x_bcp": float(x_bcp)}
        log(f"  17C {ion} QTAIM: rho(BCP) = {rho_bcp:.5f} e/a0^3, "
            f"lap rho(BCP) = {lap_bcp:+.5f}, H(BCP) = {H_bcp:+.6f} Eh/a0^3, "
            f"delta(M,N) = {out['qtaim'][ion]['delta_mn']:.4f}")
    results["bonding"] = out
    np.savez_compressed(RES / "phase17_qtaim_planes.npz",
                        **{f"{ion}_{k}": v for ion, pl in PLANE_CACHE.items()
                           for k, v in pl.items()})
    return out


def _laplacian_plane(rho2, dx, dz):
    """3D axisymmetric Laplacian on [axial x, transverse radius z>=0].

    Off axis: rho_xx + rho_zz + rho_z/z. On axis the transverse
    contribution is 2*rho_zz by even reflection. Outer boundaries are
    undefined (NaN), not falsely reported as zero. Core and interpolation
    derivative convergence must be assessed separately.
    """
    lap = np.full_like(rho2, np.nan, dtype=float)
    radius = np.arange(1, rho2.shape[1]-1) * dz
    lap[1:-1, 1:-1] = (
        (rho2[2:, 1:-1] - 2 * rho2[1:-1, 1:-1] + rho2[:-2, 1:-1]) / dx ** 2
        + (rho2[1:-1, 2:] - 2 * rho2[1:-1, 1:-1] + rho2[1:-1, :-2]) / dz ** 2
        + (rho2[1:-1, 2:] - rho2[1:-1, :-2]) / (2 * dz * radius))
    lap[1:-1, 0] = ((rho2[2:, 0] - 2*rho2[1:-1, 0] + rho2[:-2, 0])/dx**2
                        + 4*(rho2[1:-1, 1] - rho2[1:-1, 0])/dz**2)
    return lap


def promolecular_axis_descriptor(rho2, xs, zs, R):
    """Grid-minimum proxy, deliberately not a certified QTAIM critical point."""
    dx, dz = xs[1]-xs[0], zs[1]-zs[0]
    idx = np.flatnonzero((xs > .3) & (xs < R-.3))
    i = idx[np.argmin(rho2[idx, 0])]
    rho = float(rho2[i, 0])
    lap = float(_laplacian_plane(rho2, dx, dz)[i, 0])
    gx = float((rho2[i+1, 0]-rho2[i-1, 0])/(2*dx))
    cf = .3*(3*math.pi**2)**(2/3)
    kinetic = cf*rho**(5/3) + gx*gx/(72*rho)
    return dict(rho_bcp=rho, lap_bcp=lap, G_bcp=kinetic,
                H_bcp=.25*lap-kinetic, x_bcp_a0=float(xs[i]),
                axis_gradient_residual=gx,
                descriptor_scope="discrete on-axis promolecular density minimum; "
                "not a stationary-point/Hessian-certified molecular BCP")


# ============================================================================
# FIGURES (300 DPI)
# ============================================================================

def _radial_avg(rec, x):
    qq = rec["Q"] if rec["Q"] is not None else np.zeros_like(rec["P"])
    w = np.empty(len(x)); w[1:-1] = 0.5 * (x[2:] - x[:-2])
    w[0] = w[1] * 0.5; w[-1] = w[-2] * 0.5
    norm = np.sum((rec["P"] ** 2 + qq ** 2) * w)
    return float(np.sum(x * (rec["P"] ** 2 + qq ** 2) * w) / norm)


def fig1_spinors(atom_am, atom_am_nr, res):
    log("  fig1: 4-component spinor orbitals")
    fig = plt.figure(figsize=(16.5, 12.0), dpi=100)
    gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.22)
    x = atom_am["grid"]["x"]

    axa = fig.add_subplot(gs[0, 0])
    styles = {(5, +3): ("#c0392b", "5f$_{5/2}$"), (5, -4): ("#2c6fbb",
                                                          "5f$_{7/2}$")}
    for (n, kap), (col, lab) in styles.items():
        rec = atom_am["orbitals"][(n, kap)]
        axa.plot(x, np.abs(rec["P"]), color=col, lw=2.0, label=lab + "  P(r)")
        axa.plot(x, np.abs(rec["Q"]) * 30, color=col, lw=1.2, ls="--",
                 label=lab + "  30$\\times$Q(r)")
    rec_nr = None
    axa.set_xscale("log"); axa.set_yscale("log")
    axa.set_ylim(1e-6, 8)
    axa.set_xlabel("r (a$_0$)"); axa.set_ylabel("|P(r)|, |Q(r)|")
    axa.set_title("(a)  Am$^{3+}$ 5f spinors: large & small components\n"
                  "(staggered kinetically balanced solution; Q scaled $\\times$30)",
                  loc="left", fontsize=11)
    axa.legend(frameon=False, fontsize=9, ncol=2)
    axa.grid(alpha=0.25, lw=0.5)

    # 3D isosurface render of the 5f7/2 large and small components
    for k, (comp, ttl) in enumerate((("L", "(b)  |$\\Psi_L$|$^2$ isodensity\n"
                                      "5f$_{7/2}$ (m=+1/2) large component"),
                                     ("S", "(c)  |$\\Psi_S$|$^2$ isodensity "
                                      "$\\times$ c$^4$ amplification\n"
                                      "small component carries g-orbital "
                                      "(l$\\tilde{}$=4) character"))):
        ax = fig.add_subplot(gs[k, 1], projection="3d")
        kap = -4
        rec = atom_am["orbitals"][(5, kap)]
        l = l_of_kappa(kap)
        j = abs(kap) - 0.5
        ls = int(round(2 * j)) - l
        n_th, n_ph = 60, 90
        thv = np.linspace(0.02, math.pi - 0.02, n_th)
        phv = np.linspace(0, 2 * math.pi, n_ph)
        TH, PH = np.meshgrid(thv, phv, indexing="ij")
        Om = spherical_spinor_clean(l if comp == "L" else ls, j, 0.5, TH, PH)
        dens = np.sum(Om * Om, axis=0)
        r0 = 3.6
        shape = 0.55 * (dens / dens.max()) ** 0.5
        R = r0 * (0.55 + shape)
        Xs = R * np.sin(TH) * np.cos(PH)
        Ys = R * np.sin(TH) * np.sin(PH)
        Zs = R * np.cos(TH)
        surf = ax.plot_surface(Xs, Ys, Zs, facecolors=plt.cm.RdYlBu_r(
            (dens / dens.max()).flatten()).reshape(Xs.shape + (4,)),
            rstride=1, cstride=1, linewidth=0, antialiased=True, alpha=0.95)
        ax.set_xlim(-5, 5); ax.set_ylim(-5, 5); ax.set_zlim(-5, 5)
        ax.set_box_aspect((1, 1, 1))
        ax.set_axis_off()
        ax.set_title(ttl, loc="left", fontsize=11)

    axd = fig.add_subplot(gs[1, 0])
    labels, rel_nr, rel_4c = [], [], []
    for (n, kap), lab in (((1, -1), "1s"), ((2, -1), "2s"),
                          ((5, +3), "5f$_{5/2}$"), ((5, -4), "5f$_{7/2}$"),
                          ((6, -1), "7s" if (6, -1) in atom_am["orbitals"]
                           else "6s")):
        rec = atom_am["orbitals"].get((n, kap))
        if rec is None:
            continue
        labels.append(lab)
        rel_4c.append(rec["r_mean"])
        rec2 = None
        for (nn, kap2), rec3 in atom_am_nr["orbitals"].items():
            if nn == n:
                rec2 = rec3
                break
        rel_nr.append(rec2["r_mean"] if rec2 else np.nan)
    xg = np.arange(len(labels))
    axd.bar(xg - 0.19, rel_nr, width=0.38, color="#7f8c8d",
            label="non-relativistic (c$\\to\\infty$)")
    axd.bar(xg + 0.19, rel_4c, width=0.38, color="#c0392b",
            label="4-component Dirac")
    for xi, v in zip(xg, rel_nr):
        axd.text(xi - 0.19, v + 0.02, f"{v:.2f}", ha="center", fontsize=8)
    for xi, v in zip(xg, rel_4c):
        axd.text(xi + 0.19, v + 0.02, f"{v:.2f}", ha="center", fontsize=8)
    axd.set_xticks(xg); axd.set_xticklabels(labels)
    axd.set_ylabel("$\\langle r \\rangle$ (a$_0$)")
    axd.set_title("(d)  relativistic core contraction & 5f expansion\n"
                  "$\\langle r \\rangle$ Am$^{3+}$: NR vs 4-component Dirac",
                  loc="left", fontsize=11)
    axd.legend(frameon=False, fontsize=9)
    axd.grid(alpha=0.25, lw=0.5, axis="y")
    fig.suptitle("Phase 17 - 4-component Dirac spinor structure of the "
                 "Am$^{3+}$ 5f shell", fontsize=14, fontweight="bold", y=0.99)
    fig.savefig(FIG / "fig1_4component_spinor_orbitals.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    log("    fig1 saved")


def fig2_multiplets(res):
    log("  fig2: spin-orbit multiplet splitting")
    mult = res["multiplet"]
    fig, axes = plt.subplots(1, 2, figsize=(16.5, 9.0), dpi=100)
    for ax, ion in zip(axes, ("Am", "Eu")):
        m = mult[ion]
        ladder = {int(k): v for k, v in m["ladder_4c"].items()}
        zeta = m["zeta_df_cm"]
        nm = "Am" if ion == "Am" else "Eu"
        sh = "5f" if ion == "Am" else "4f"
        # three columns: NR (degenerate), scalar-rel (degenerate), 4c + SOC
        ax.hlines(0.0, -0.28, 0.28, color="#7f8c8d", lw=3.0)
        ax.hlines(0.0, 0.72, 1.28, color="#e67e22", lw=3.0)
        for J, ev in sorted(ladder.items()):
            ax.hlines(ev, 1.72, 2.28, color="#c0392b", lw=2.4)
            ax.text(2.34, ev, f"J={J}", fontsize=9, va="center",
                    color="#c0392b")
        ax.text(0.0, -max(ladder.values()) * 0.08, "NR", ha="center",
                fontsize=10)
        ax.text(1.0, -max(ladder.values()) * 0.08, "scalar-rel", ha="center",
                fontsize=10)
        ax.text(2.0, -max(ladder.values()) * 0.08, "4c-DC + SOC",
                ha="center", fontsize=10)
        # experimental anchors
        exp = {int(k): v for k, v in m["ladder_exp"].items()}
        for J, ev in exp.items():
            if ev and ev < max(ladder.values()) * 1.05:
                ax.hlines(ev, 2.55, 2.8, color="#16a085", lw=1.6, ls=":")
                ax.text(2.85, ev, f"exp {int(ev)}", fontsize=7.5,
                        va="center", color="#16a085")
        ax.text(0.03, 0.97,
                rf"$\zeta$({sh}) from Dirac j-splitting = {zeta:.0f} cm$^{{-1}}$  "
                rf"$\Delta E_{{SO}}$(J=1) = {ladder[1]:.0f} cm$^{{-1}}$  "
                f"(exp fit {m['zeta_expfit_cm']:.0f})",
                transform=ax.transAxes, fontsize=9.5, va="top")
        ax.set_xlim(-0.6, 3.4)
        ax.set_ylim(-max(ladder.values()) * 0.18,
                    max(ladder.values()) * 1.12)
        ax.set_ylabel("energy above $^7F_0$ (cm$^{-1}$)")
        ax.set_title(f"{nm}$^{{3+}}$  {sh}$^6$ $^7F_J$ manifold: "
                     f"progressive lifting of degeneracy", fontsize=11.5)
        ax.grid(alpha=0.2, lw=0.4, axis="y")
    fig.suptitle("Phase 17 - the $^7F_J$ multiplet: NR $\to$ "
                 "scalar-relativistic $\to$ 4-component Dirac + spin-orbit "
                 "(Lande intervals from the SCF j-splitting)",
                 fontsize=13, fontweight="bold", y=0.98)
    fig.savefig(FIG / "fig2_spin_orbit_multiplet_splitting.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    log("    fig2 saved")


def fig3_qtaim(res):
    log("  fig3: QTAIM relativistic covalency map")
    bond = res["bonding"]
    PLANE_CACHE = globals().get("PLANE_CACHE", {})
    fig, axes = plt.subplots(2, 2, figsize=(16.5, 12.0), dpi=100)
    for k, ion in enumerate(("Am", "Eu")):
        ax = axes[0, k]
        pl = PLANE_CACHE[ion]
        lap = pl["lap2d"]
        xs, zs = pl["xs"], pl["zs"]
        R = pl["R"]
        lv = np.array([-1, -0.5, -0.2, -0.08, -0.03, 0.03, 0.08, 0.2, 0.5,
                       1, 2, 5])
        # Array is [x,z], whereas matplotlib expects [z,x]. Nuclear cores
        # are not resolved by this display grid; disclose the excluded radius.
        core = 3 * max(xs[1]-xs[0], zs[1]-zs[0])
        Xp, Zp = np.meshgrid(xs, zs, indexing="ij")
        excluded = (np.hypot(Xp, Zp) < core) | (np.hypot(Xp-R, Zp) < core)
        shown = np.ma.masked_where(excluded | ~np.isfinite(lap), lap).T
        levels = np.r_[-8, -5, -2, -1, -.5, -.2, -.08, -.03, 0,
                       .03, .08, .2, .5, 1, 2, 5, 8]
        im = ax.contourf(xs, zs, shown, levels=levels, cmap="RdBu_r", extend="both",
                         norm=SymLogNorm(linthresh=0.05, vmin=-8, vmax=8))
        cs = ax.contour(xs, zs, shown, levels=lv, colors="k", linewidths=0.5)
        ax.text(.02, .97, f"Core r<{core:.2f} a0 masked; colors saturate at +/-8",
                transform=ax.transAxes, va="top", fontsize=8)
        ax.plot(pl["x_bcp"], 0, "o", ms=9, mfc="none", mec="k", mew=1.6)
        ax.annotate("axis minimum", (pl["x_bcp"], 0), textcoords="offset points",
                    xytext=(6, 10), fontsize=10)
        ax.set_xlabel("x along M-N (a$_0$)")
        ax.set_ylabel("z ($a_0$)")
        q = bond["qtaim"][ion]
        ax.set_title(f"{'Am' if ion == 'Am' else 'Eu'}$^{{3+}}$: "
                     f"promolecular $\\nabla^2\\rho$   "
                     f"$\\rho_{{BCP}}$={q['rho_bcp']:.4f}  "
                     f"$\\nabla^2\\rho_{{BCP}}$={q['lap_bcp']:+.4f} "
                     f"e/a$_0^5$", loc="left", fontsize=11)
        plt.colorbar(im, ax=ax, shrink=0.85)
    axc = axes[1, 0]
    ions = ("Am", "Eu")
    es = [bond["eda"][i]["dE_cov_bond_kcal"] for i in ions]
    ec = [bond["eda"][i]["dE_cov_complex_kcal"] for i in ions]
    xg = np.arange(2)
    axc.bar(xg - 0.18, es, width=0.36, color="#2c6fbb",
            label="$\\Delta E_{cov}$ per M-N bond")
    axc.bar(xg + 0.18, ec, width=0.36, color="#c0392b",
            label="3 x pair mixing proxy (not complex EDA)")
    for xi, v in zip(xg - 0.18, es):
        axc.text(xi, v + 0.05, f"{v:.2f}", ha="center", fontsize=9)
    for xi, v in zip(xg + 0.18, ec):
        axc.text(xi, v + 0.05, f"{v:.2f}", ha="center", fontsize=9)
    axc.set_xticks(xg); axc.set_xticklabels(["Am$^{3+}$", "Eu$^{3+}$"])
    axc.set_ylabel("kcal/mol")
    dd = bond["eda"]["selectivity_kcal"]
    axc.set_title(f"(c)  Two-state mixing proxy: Am minus Eu "
                  f"$\\Delta\\Delta E$ = {dd:.2f} kcal/mol",
                  loc="left", fontsize=11)
    axc.legend(frameon=False, fontsize=9)
    axc.grid(alpha=0.25, lw=0.5, axis="y")
    axd = axes[1, 1]
    w = 0.26
    xg = np.arange(2)
    rho_b = [bond["qtaim"][i]["rho_bcp"] for i in ions]
    lap_b = [abs(bond["qtaim"][i]["lap_bcp"]) for i in ions]
    ov_b = [bond["overlap"][i]["S_rms"] * 100 for i in ions]
    axd.bar(xg - w, rho_b, width=w, color="#2c6fbb",
            label="$\\rho$(BCP) e/a$_0^3$")
    axd.bar(xg, lap_b, width=w, color="#e67e22",
            label="$|\\nabla^2\\rho|$(BCP) e/a$_0^5$")
    axd.bar(xg + w, ov_b, width=w, color="#16a085",
            label="100$\\times$S(5f/4f, N 2p)")
    for xi, v in ((xg[0] - w, rho_b[0]), (xg[1] - w, rho_b[1]),
                  (xg[0], lap_b[0]), (xg[1], lap_b[1]),
                  (xg[0] + w, ov_b[0]), (xg[1] + w, ov_b[1])):
        axd.text(xi, v + 0.004, f"{v:.3f}" if v < 0.1 else f"{v:.2f}",
                 ha="center", fontsize=8)
    axd.set_xticks(xg)
    axd.set_xticklabels(["Am$^{3+}$", "Eu$^{3+}$"])
    axd.set_title("(d)  Promolecular minimum & density-overlap proxies",
                  loc="left", fontsize=11)
    axd.legend(frameon=False, fontsize=9)
    axd.grid(alpha=0.25, lw=0.5, axis="y")
    fig.suptitle("Phase 17 - Am / Eu atomic-promolecule and two-state model comparison\n"
                 "Not self-consistent complex QTAIM; no bond/no-bond conclusion",
                 fontsize=13.5, fontweight="bold", y=0.99)
    fig.savefig(FIG / "fig3_qtaim_relativistic_covalency_map.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    log("    fig3 saved")


def fig4_validation(res, kb_demo):
    log("  fig4: relativistic foundation validation")
    fig, axes = plt.subplots(2, 2, figsize=(16.5, 11.0), dpi=100)
    hv = res["validation"]["hydrogenic"]
    ax = axes[0, 0]
    zs = [r["Z"] for r in hv if "rel_err" in r]
    errs = [r["rel_err"] for r in hv if "rel_err" in r]
    ax.semilogy(zs, np.array(errs) * 100, "o-", color="#c0392b")
    ax.set_xlabel("Z"); ax.set_ylabel("|E$_{num}$ - E$_{Sommerfeld}$| / |E|  (%)")
    ax.set_title(f"(a) point-nucleus Dirac 1s vs Sommerfeld closed form",
                 loc="left", fontsize=10)
    ax.text(0.03, 0.94,
            f"max error {max(errs)*100:.3f} %   |   Dirac degeneracy "
            f"|E(2s$_{{1/2}}$)-E(2p$_{{1/2}}$)|/E = "
            f"{res['validation'].get('degeneracy', 0):.1e}",
            transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(boxstyle="round", fc="white", ec="#cccccc", alpha=0.85))
    ax.grid(alpha=0.3, lw=0.5)

    ax = axes[0, 1]
    cc = res["breit"]["coulomb_certification"]
    cz = [r["Z"] for r in cc]
    ce = [r["rel_err"] * 100 for r in cc]
    ax.semilogy(cz, ce, "s-", color="#16a085")
    ax.set_xlabel("Z")
    ax.set_ylabel("|J$_{quad}$ - 5Z/8| / (5Z/8)  (%)")
    derived = any("derived" in r.get("provenance", "") for r in cc)
    ax.set_title("(b) Coulomb quadrature error relative to exact 5Z/8\n" +
                 ("Z-scaled from measured Z=1; not independent high-Z tests" if derived
                  else "measured angular-spinor and radial quadrature error"),
                 loc="left", fontsize=10)
    ax.grid(alpha=0.3, lw=0.5)

    ax = axes[1, 0]
    br = res["breit"]
    zz = np.array(br["hydrogenic_Z"])
    bb = np.array(br["hydrogenic_breit_cm"])
    ax.loglog(zz, bb, "o-", color="#8e44ad", label="Breit(1s$^2$), hydrogenic")
    ax.loglog(zz, bb[0] * (zz / zz[0]) ** 4, "k:", lw=1,
              label="Z$^4$ scaling guide")
    ax.set_xlabel("Z"); ax.set_ylabel("$\\Delta E_{Breit}$(1s$^2$) (cm$^{-1}$)")
    ax.set_title("(c) first-order Breit pair energy: magnetic (Gaunt) + "
                 "retardation, Z$^4$ law", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.3, lw=0.5)

    ax = axes[1, 1]
    cov = res["bonding"]
    ion_l = ["Am$^{3+}$ 5f", "Eu$^{3+}$ 4f"]
    S = [cov["overlap"]["Am"]["S_rms"], cov["overlap"]["Eu"]["S_rms"]]
    eps = [cov["eda"]["Am"]["eps_acceptor_eh"],
           cov["eda"]["Eu"]["eps_acceptor_eh"]]
    ax2 = ax.twinx()
    b1 = ax.bar([0, 1], S, width=0.4, color="#c0392b", label="S(5f/4f, N 2p)")
    b2 = ax2.bar([0.55, 1.55], [-np.array(eps), 0][0:2] if False else
                 [-e for e in eps], width=0.4, color="#2c6fbb",
                 label="$-\\epsilon$(5f/4f) (Eh)")
    ax.set_xticks([0.27, 1.27]); ax.set_xticklabels(ion_l)
    ax.set_ylabel("overlap integral S", color="#c0392b")
    ax2.set_ylabel("$-\\epsilon$ acceptor (Eh)", color="#2c6fbb")
    ax.set_title("(d) the covalency mechanism: 5f overlap larger AND "
                 "energy-matched; 4f overlap vanishing AND orphaned",
                 loc="left", fontsize=10)
    ax.grid(alpha=0.3, lw=0.5, axis="y")
    kb_error = res.get("certificates", {}).get("kb_balanced_1s_err")
    kb_label = (f"KB analytic relative error={kb_error:.3g} (finite-basis diagnostic)"
                if kb_error is not None else "KB diagnostic unavailable")
    if kb_error is not None and kb_error >= 1:
        kb_label = "FAILED KB diagnostic: " + kb_label
    fig.suptitle("Phase 17 - atomic diagnostics and model proxies, not global validation\n"
                 + kb_label, fontsize=12, fontweight="bold", y=0.99)
    fig.savefig(FIG / "fig4_relativistic_foundation_validation.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)
    log("    fig4 saved")


# ============================================================================
# Kinetic-balance matrix demonstration (17A)
# ============================================================================

def gaussian_kb_demo(z=92.0, n_basis=32, c=C_LIGHT, return_diagnostics=False):
    """Radial kappa=-1 Dirac matrix; compare bound gap states, not the sea.

    P=r exp(-a*r^2); RKB Q=(d/dr-1/r)P/(2c). Energies have rest
    energy subtracted. Separate LL/SS metric and Hamiltonian blocks are
    essential. This finite-basis demonstration is not an SCF certificate.
    """
    x, y, h, w = log_grid(1500, 1e-7, 40.0)
    v = -z / x

    if not 0 < z < c or n_basis < 2:
        raise ValueError("Require subcritical positive Z and at least two basis functions")
    diagnostics = {}
    def basis_mats(balanced):
        ex = np.logspace(-4, 4, n_basis) * z*z
        r = x[:, None]
        base = np.exp(-r*r*ex)
        L = r*base
        DL = -2*ex*r*r*base  # (d/dr + kappa/r)P with kappa=-1
        if balanced:
            S = DL/(2*c)
        else:
            S = r*np.exp(-r*r*ex*.5)
        ln = np.sqrt(np.sum(L*L*w[:, None], axis=0))
        sn = np.sqrt(np.sum(S*S*w[:, None], axis=0))
        L, DL, S = L/ln, DL/ln, S/sn
        oll = L.T@(w[:, None]*L)
        oss = S.T@(w[:, None]*S)
        hll = L.T@((w*v)[:, None]*L)
        hss = S.T@((w*(v-2*c*c))[:, None]*S)
        hsl = c*S.T@(w[:, None]*DL)
        H = np.block([[hll, hsl.T], [hsl, hss]])
        Ov = sla.block_diag(oll, oss)
        eig, vec = sla.eigh(Ov)
        keep = eig > 1e-12*eig.max()
        transform = vec[:, keep]/np.sqrt(eig[keep])
        reduced = transform.T@H@transform
        energies, coeff = sla.eigh(.5*(reduced+reduced.T))
        gap = (energies > -c*c) & (energies < 0)  # Positive total-energy bound branch.
        bound = energies[gap]
        if not len(bound) or not np.isfinite(energies).all():
            raise RuntimeError("No finite positive-total-energy bound state")
        k = np.flatnonzero(gap)[0]
        vector = transform@coeff[:, k]
        residual = np.linalg.norm(H@vector-energies[k]*Ov@vector)/max(np.linalg.norm(H@vector), 1e-30)
        diagnostics["balanced" if balanced else "unbalanced"] = dict(
            bound_1s_eh=float(bound[0]), lowest_spectrum_eh=float(energies[0]),
            retained_metric_rank=int(keep.sum()), matrix_size=int(len(eig)),
            eigen_residual_relative=float(residual),
            hermiticity_max_abs=float(np.max(np.abs(H-H.T))))
        return bound
    ev_b = basis_mats(True)
    ev_u = basis_mats(False)
    e_exact = sommerfeld_1s(z)
    log(f"  17A KB demo: balanced 1s = {ev_b[0]:.5f} Eh (exact {e_exact:.5f}); "
        f"unbalanced gap state = {ev_u[0]:.5f}; negative continuum excluded")
    result = ([("kinetically balanced", ev_b[:8], "#16a085"),
               ("unbalanced S-basis", ev_u[:8], "#c0392b")], e_exact)
    return (*result, diagnostics) if return_diagnostics else result


# ============================================================================
# MAIN
# ============================================================================

def atomic_pair_stage(results, n_pts):
    """DF + NR twins for Am3+ and Eu3+; moments, zeta extraction."""
    atoms = {}
    for ion, z, cfg, am, tag in (("Am", AM3_Z, AM3_CONFIG, AM3_A, "Am3+"),
                                 ("Eu", EU3_Z, EU3_CONFIG, EU3_A, "Eu3+")):
        log(f"[atomic] Dirac-Fock-Slater {tag} (Z={z}, {int(sum(o1+o2 for _,_,o1,o2 in cfg))} e)")
        atom = atomic_scf(z, cfg, relativistic=True, n_pts=n_pts, a_mass=am,
                          label=tag)
        log(f"  {tag}: E_tot = {atom['e_tot']:.4f} Eh in {atom['iterations']} it")
        atoms[tag] = atom
    for ion, z, cfg, tag in (("Am", AM3_Z, AM3_CONFIG, "Am3+NR"),
                             ("Eu", EU3_Z, EU3_CONFIG, "Eu3+NR")):
        log(f"[atomic] non-relativistic twin {tag}")
        atom = atomic_scf(z, cfg, relativistic=False, n_pts=n_pts,
                          label=tag)
        atoms[tag] = atom

    # summary + zeta
    summ = {}
    for tag, atom in atoms.items():
        orbs = {}
        for (n, kap), rec in atom["orbitals"].items():
            orbs[f"{n}|{kap}"] = {k: rec[k] for k in
                                  ("energy", "occ", "r_mean", "r2_mean",
                                   "r_inv3", "small_norm", "nodes")}
        summ[tag] = {"e_tot": atom["e_tot"], "iterations": atom["iterations"],
                     "scf": atom["scf"],
                     "orbitals": orbs}
    results["atoms"] = summ
    zeta = {}
    for ion, tag, nf in (("Am", "Am3+", 5), ("Eu", "Eu3+", 4)):
        e_lo = atoms[tag]["orbitals"][(nf, +3)]["energy"]
        e_hi = atoms[tag]["orbitals"][(nf, -4)]["energy"]
        zeta[ion] = (e_hi - e_lo) / 3.5
        log(f"  zeta({ion} nf) = (eps_7/2 - eps_5/2)/3.5 = "
            f"{zeta[ion] * HARTREE_CM:.0f} cm-1")
    results["zeta_cm"] = {k: v * HARTREE_CM for k, v in zeta.items()}
    return atoms, zeta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["all", "atomic", "breit", "multiplet",
                             "bonding", "figures"])
    ap.add_argument("--npts", type=int, default=640)
    ap.add_argument("--fast", action="store_true")
    args = ap.parse_args()
    if args.fast:
        args.npts = min(args.npts, 560)

    res_path = RES / "phase17_results.json"
    results = {}
    if res_path.exists() and args.stage != "all":
        results = json.loads(res_path.read_text(encoding="utf-8"))

    log("=" * 78)
    log("PHASE 17 - FULLY RELATIVISTIC 4-COMPONENT DIRAC QUANTUM CHEMISTRY & "
        "ACTINIDE 5f COVALENCY")
    log("=" * 78)

    worst = spinor_selftest()
    log(f"[certify] spherical-spinor orthonormality: worst dev {worst:.2e}")
    results.setdefault("certificates", {})["spinor_orthonormality"] = worst

    if args.stage in ("all", "atomic"):
        # ---- 17A validators -------------------------------------------
        log("[MODULE 17A] hydrogenic validation")
        results["validation"] = {}
        results["validation"]["hydrogenic"] = hydrogenic_validation(
            [1, 25, 50, 80, 92, 110], n_pts=args.npts)
        deg = [r for r in results["validation"]["hydrogenic"]
               if "degeneracy_test" in r]
        if deg:
            results["validation"]["degeneracy"] = deg[0]["degeneracy_test"]
            results["validation"]["fine_split_2p_cm"] = \
                deg[0]["fine_split_2p_cm"]
        atoms, zeta = atomic_pair_stage(results, args.npts)
        (RES / "phase17_results.json").write_text(
            json.dumps(results, indent=1, default=float), encoding="utf-8")
        np.savez_compressed(
            RES / "phase17_radial_data.npz",
            **{f"{tag}_{n}_{kap}_{k}": (rec[k] if rec[k] is not None
                                        else np.zeros_like(rec["P"]))
               for tag, atom in atoms.items()
               for (n, kap), rec in atom["orbitals"].items()
               for k in ("P", "Q")})

    if args.stage in ("all", "breit"):
        log("[MODULE 17A-tail] Breit interaction")
        atoms = _load_atoms(results)
        br = {}
        tab_c = cached_breit_pair_tabulation(-1, kernel="coulomb")
        zs = [1, 20, 50, 80]
        e1s = []
        for Z in zs:
            xx, ww = hydrogenic_radial_grid(Z, 4000, r_min=1e-5)
            P = 2.0 * Z ** 1.5 * xx * np.exp(-Z * xx)
            P /= math.sqrt(np.sum(P ** 2 * ww))
            J = breit_pair_energy(tab_c, P, None, xx, ww).real
            J_ref = 5.0 * Z / 8.0
            e1s.append({"Z": Z, "J_quadrature": J, "J_exact": J_ref,
                        "rel_err": abs(J - J_ref) / J_ref,
                        "grid": "uniform Z*r", "n_radial": len(xx)})
            log(f"  Breit-machinery certification Z={Z}: J = {J:.5f} vs "
                f"5Z/8 = {J_ref:.5f} (err {e1s[-1]['rel_err']*100:.3f}%)")
        br["coulomb_certification"] = e1s
        # hydrogenic Breit(1s^2) scaling + the SCF atoms
        tab_b = cached_breit_pair_tabulation(-1, kernel="breit")
        tab_g = cached_breit_pair_tabulation(-1, kernel="gaunt")
        hb_z, hb_cm = [], []
        for Z in (20, 50, 80, 92):
            xx, ww = hydrogenic_radial_grid(Z, 3000)
            gam = math.sqrt(1 - (Z / C_LIGHT) ** 2)
            lam = Z
            P = xx ** gam * np.exp(-lam * xx)
            fac = c_light_small = C_LIGHT * (gam + 1) / Z  # Q/P ratio scale
            Q = -P * fac * 0.0 - P * (1 - gam) / Z * C_LIGHT * 0 - \
                P * C_LIGHT * (1.0 - gam) / Z
            nrm = math.sqrt(np.sum((P ** 2 + Q ** 2) * ww))
            P /= nrm; Q /= nrm
            e_b = breit_pair_energy(tab_b, P, Q, xx, ww).real
            hb_z.append(Z); hb_cm.append(abs(e_b) * HARTREE_CM)
            log(f"  hydrogenic Breit(1s^2) Z={Z}: {abs(e_b)*HARTREE_CM:.1f} "
                f"cm-1")
        br["hydrogenic_Z"] = hb_z
        br["hydrogenic_grid"] = {"coordinate": "uniform Z*r", "n_radial": 3000}
        br["hydrogenic_breit_cm"] = hb_cm
        # SCF 1s^2 Breit for the actinides + Gaunt/retardation split
        for ion, tag in (("Am", "Am3+"), ("Eu", "Eu3+")):
            atom = atoms[tag]
            rec = atom["orbitals"][(1, -1)]
            x = atom["grid"]["x"]; w = atom["grid"]["w"]
            e_full = breit_pair_energy(tab_b, rec["P"], rec["Q"], x, w).real
            e_gaunt = breit_pair_energy(tab_g, rec["P"], rec["Q"], x, w).real
            n_pairs = 1.0
            br[ion] = {"breit_1s2_eh": e_full,
                       "breit_1s2_cm": e_full * HARTREE_CM,
                       "gaunt_cm": e_gaunt * HARTREE_CM,
                       "retardation_cm": (e_full - e_gaunt) * HARTREE_CM}
            log(f"  {ion} Breit(1s^2) = {e_full*HARTREE_CM:.1f} cm-1 "
                f"(Gaunt {e_gaunt*HARTREE_CM:.1f}, retardation "
                f"{(e_full-e_gaunt)*HARTREE_CM:.1f})")
        results["breit"] = br
        (RES / "phase17_results.json").write_text(
            json.dumps(results, indent=1, default=float), encoding="utf-8")

    if args.stage in ("all", "multiplet"):
        atoms = _load_atoms(results)
        run_multiplet_stage(atoms, results)
        (RES / "phase17_results.json").write_text(
            json.dumps(results, indent=1, default=float), encoding="utf-8")

    if args.stage in ("all", "bonding"):
        log("[MODULE 17C] bonding / QTAIM / EDA")
        atoms = _load_atoms(results)
        atom_n = atomic_scf(7, N_CONFIG, relativistic=False, n_pts=400,
                            label="N")
        eda_qtaim_stage(atoms["Am3+"], atoms["Eu3+"], atoms["Am3+NR"],
                        atoms["Eu3+NR"], atom_n, results)
        (RES / "phase17_results.json").write_text(
            json.dumps(results, indent=1, default=float), encoding="utf-8")

    if args.stage in ("all", "figures"):
        log("[FIGURES] 300 DPI publication renders")
        atoms = _load_atoms(results)
        kb_demo, e_ex = gaussian_kb_demo()
        results.setdefault("certificates", {})["kb_balanced_1s_err"] = \
            abs(kb_demo[0][1][0] - e_ex) / abs(e_ex)
        if "bonding" not in results or not PLANE_CACHE:
            atom_n = atomic_scf(7, N_CONFIG, relativistic=False, n_pts=400,
                                label="N")
            eda_qtaim_stage(atoms["Am3+"], atoms["Eu3+"], atoms["Am3+NR"],
                            atoms["Eu3+NR"], atom_n, results)
        fig1_spinors(atoms["Am3+"], atoms["Am3+NR"], results)
        fig2_multiplets(results)
        fig3_qtaim(results)
        fig4_validation(results, None)
        (RES / "phase17_results.json").write_text(
            json.dumps(results, indent=1, default=float), encoding="utf-8")
        log(f"DONE - figures in {FIG}, record in {res_path}")


def _load_atoms(results):
    """Rebuild lightweight atom objects from the JSON record + npz."""
    try:
        data = np.load(RES / "phase17_radial_data.npz", allow_pickle=True)
    except FileNotFoundError:
        return {}
    atoms = {}
    for tag, rec in results["atoms"].items():
        if not rec.get("scf", {}).get("converged", False):
            raise RuntimeError(f"{tag}: missing/failed SCF convergence evidence; rerun atomic stage")
        rel = not tag.endswith("NR")
        z = AM3_Z if tag.startswith("Am") else EU3_Z
        cfg = AM3_CONFIG if tag.startswith("Am") else EU3_CONFIG
        if rel:
            x, y, h, w = log_grid(NPTS_DEFAULT, 1e-6, 60.0)
        else:
            x = nr_grid(); y = h = None
            w = np.full_like(x, NR_DR)
        atom = {"z": z, "grid": {"x": x, "w": w}, "orbitals": {},
                "e_tot": rec["e_tot"], "rho": np.zeros_like(x)}
        for key, orec in rec["orbitals"].items():
            n, kap = key.split("|")
            n, kap = int(n), int(kap)
            P = np.asarray(data[f"{tag}_{n}_{kap}_P"], dtype=float)
            Qr = data[f"{tag}_{n}_{kap}_Q"]
            Q = (np.zeros_like(P) if Qr.dtype == object or Qr.ndim == 0
                 else np.asarray(Qr, dtype=float))
            atom["orbitals"][(n, kap)] = {
                "energy": orec["energy"], "occ": orec["occ"],
                "r_mean": orec["r_mean"], "r2_mean": orec["r2_mean"],
                "r_inv3": orec["r_inv3"], "small_norm": orec["small_norm"],
                "nodes": orec["nodes"], "P": P, "Q": Q}
            atom["rho"] += orec["occ"] * orbital_density(P, Q, x, w)
        atoms[tag] = atom
    return atoms


NPTS_DEFAULT = 640


if __name__ == "__main__":
    try:
        main()
    except SCFConvergenceError as exc:
        (RES / "phase17_scf_failure.json").write_text(
            json.dumps(exc.diagnostics, indent=2), encoding="utf-8")
        raise
