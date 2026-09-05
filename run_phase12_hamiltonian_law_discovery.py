#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_phase12_hamiltonian_law_discovery.py — PHASE 12 SUPREME MISSION
Autonomous Hamiltonian Law Discovery, Symbolic Symmetry Induction &
Non-Equilibrium Entropy Laws.

Transmutes the repository from a computer of known equations into an
Autonomous Scientific Theorist: given ONLY noisy, high-dimensional
observational telemetry of a stiff, non-equilibrium, limit-cycle chemical
oscillator (Field-Körös-Noyes / Oregonator model of the Belousov-Zhabotinsky
reaction, plus its 1-D reaction-diffusion extension), the pipeline
autonomously deduces

  1. the symbolic non-linear differential equations governing the continuous
     time evolution (Module 12B),
  2. the hidden conserved first integral H(x) of the invariant-manifold core
     and the dissipative Lyapunov functional V(x) >= 0 certifying thermo-
     dynamic stability (Module 12C),
  3. the Onsager reciprocity / entropy-production structure of the discovered
     dynamics in the near-steady-state regime (Module 12B, thermodynamic
     consistency gate),

and exports the discovered laws as LaTeX, SymPy strings and C++ execution
kernels, with publication figures (300 DPI) in ./figures_phase12/.

The full epistemic loop implemented here (each stage mechanically triggered
by the state of the evidence, not hard-coded):

    OBSERVE    passive telemetry, sigma_noise = 5 % Gaussian white noise
    ESTIMATE   weak-form (integration-by-parts) split-sample estimator with
               Simpson quadrature; errors-in-variables attenuation removed by
               a Monte-Carlo noise-calibrated Gram matrix
    SELECT     STRidge (sequential thresholded ridge) path + forward/backward
               elimination + random-restart iterative hard thresholding ->
               per-complexity Pareto front -> studentized noise-floor (Occam)
               knee
    INTERVENE  when rival supports are statistically indistinguishable on
               passive data (the u ~ w low-pass degeneracy of the Oregonator
               makes {u^2} and {u w} near-collinear on the attractor), the
               theorist designs factorial impulse experiments (single-variable
               kicks, basin-clipped) to break the observational confounding
    ARBITRATE  remaining rivals are decided by higher-precision replication
               of the decisive experiments (model SELECTION only; coefficients
               are always refit on the operative 5 % telemetry)
    DISCOVER   coefficients refit on all operative 5 % telemetry with the
               calibrated Gram -> unbiased symbolic law

True system (dimensionless Field-Noyes Oregonator, derived from the FKN
mechanism; u = [HBrO2], v = [Br-], w = [Ce4+], time in tau units):

    eps1 * du/dt = q*v - u*v + u*(1-u)          eps1 = 0.10
    eps2 * dv/dt = -q*v - u*v + f*w             eps2 = 0.01
    dw/dt = u - w                               q = 2.5e-3, f = 1.1

f = 1.1 -> sustained limit cycle (period ~ 7 tau); f = 2.8 -> excitable
regime with a stable steady state (used for the Onsager near-equilibrium
gate).  The 1-D reaction-diffusion extension adds D_u * lap u and D_v * lap v
to the u and v equations (periodic domain, Strang splitting with exact
spectral diffusion + RK4 reaction substeps).

Everything the theorist is allowed to know is the telemetry; the parameter
values above are used ONLY for scoring the discovery (they never enter the
estimator).
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import sympy
from scipy.integrate import solve_ivp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# --------------------------------------------------------------------------- #
#  global configuration
# --------------------------------------------------------------------------- #
ROOT = Path(__file__).resolve().parent
RES = ROOT / "results_phase12"
FIG = ROOT / "figures_phase12"
RES.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

# true dimensionless Oregonator parameters (used ONLY to score the discovery)
E1, E2, Q, F_OSC, F_STO = 0.10, 1.0e-2, 2.5e-3, 1.1, 2.8
D_U, D_V = 8.0e-3, 3.0e-3            # tau-unit diffusion coefficients (1-D RD)
NOISE = 0.05                          # mission telemetry noise (5 %)
NOISE_ARB = 0.01                      # arbitration instrument precision (1 %)
DT = 1.0e-3                           # fine sampling of 0-D telemetry
T0 = time.time()
QUICK = False

C_MAIN, C_ACC, C_MID, C_GOLD, C_GREY = "#1F5FA8", "#C0392B", "#2E8B57", "#B7950B", "#5D6D7E"
plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 300, "font.size": 9.5,
    "axes.titlesize": 10.5, "axes.labelsize": 10,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "axes.axisbelow": True, "legend.framealpha": 0.9,
    "axes.unicode_minus": False,
})


def banner(txt):
    bar = "=" * 78
    print(f"\n{bar}\n  {txt}\n{bar}", flush=True)


def elapsed():
    return f"[{time.time() - T0:7.1f} s]"


# --------------------------------------------------------------------------- #
#  MODULE 12A — non-equilibrium dissipative simulation & telemetry injection
# --------------------------------------------------------------------------- #
def oregonator_rhs(s, f=F_OSC):
    """Dimensionless Field-Noyes Oregonator (mechanism-derived)."""
    u, v, w = s
    return np.array([
        (Q * v - u * v + u * (1.0 - u)) / E1,
        (-Q * v - u * v + f * w) / E2,
        u - w,
    ])


def oregonator_rhs_vec(S):
    u, v, w = S[..., 0], S[..., 1], S[..., 2]
    return np.stack([
        (Q * v - u * v + u * (1.0 - u)) / E1,
        (-Q * v - u * v + F_OSC * w) / E2,
        u - w,
    ], axis=-1)


def simulate_0d(s0, T, f=F_OSC, dt=DT):
    """Stiff integration (Radau) sampled on the telemetry grid."""
    t_eval = np.arange(0.0, T, dt)
    sol = solve_ivp(lambda t, s: oregonator_rhs(s, f), (0.0, T), np.asarray(s0, float),
                    method="Radau", rtol=1e-9, atol=1e-11, t_eval=t_eval, max_step=0.5)
    return sol.y.T  # (n, 3)


def kicked_trajectory(s0, kick_var, kick_amp, t_kick, T=25.0, f=F_OSC):
    """Impulse experiment: single-variable kick at t_kick, basin-clipped.

    Returns the pre- and post-kick SEGMENTS separately: the state jump is a
    genuine discontinuity of the trajectory, and weak-form windows must not
    straddle it (the integration-by-parts identity picks up a Delta x phi(t_k)
    term no smooth model can explain).
    """
    pre = simulate_0d(s0, t_kick, f=f)
    sk = np.maximum(pre[-1] + np.eye(3)[kick_var] * kick_amp, 0.02)
    post = simulate_0d(sk, T - t_kick, f=f)
    return [pre, post]


def simulate_rd(T=25.0, n_grid=192, L=0.6, seed=0, n_ic=2):
    """1-D reaction-diffusion Oregonator.

    Strang splitting: exact spectral diffusion (periodic) half-step, RK4
    reaction substeps, diffusion half-step.  Returns list of (frames,
    n_grid, 3) arrays sampled every 3e-3 tau.
    """
    rng = np.random.default_rng(seed)
    dt_m, dt_sub, dt_rd = 5.0e-3, 2.5e-4, 5.0e-3
    kk = np.fft.rfftfreq(n_grid, d=L / n_grid) * 2.0 * np.pi
    damp = [np.exp(-D * kk * kk * dt_m / 2.0) if D > 0 else None
            for D in (D_U, D_V, 0.0)]

    def diffuse_half(S):
        for c, d in enumerate(damp):
            if d is not None:
                x = S[..., c]
                if x.ndim == 1:                      # single snapshot (n_grid,)
                    x = np.fft.irfft(np.fft.rfft(x) * d, n=n_grid)
                else:                                # frame stack (..., n_grid)
                    x = np.fft.irfft(np.fft.rfft(x, axis=-1) * d[None, :],
                                     n=n_grid, axis=-1)
                S[..., c] = x
        return S

    def react(S):
        n_sub = int(round(dt_m / dt_sub))
        h = dt_m / n_sub
        for _ in range(n_sub):
            k1 = oregonator_rhs_vec(S)
            k2 = oregonator_rhs_vec(S + 0.5 * h * k1)
            k3 = oregonator_rhs_vec(S + 0.5 * h * k2)
            k4 = oregonator_rhs_vec(S + h * k3)
            S = S + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        return np.maximum(S, 0.0)

    out = []
    every = int(round(dt_rd / dt_m))
    for ic in range(n_ic):
        base = simulate_0d(np.array([0.6, 0.6, 0.6]), 5.0 + 3.7 * ic)[-1]
        smooth = np.cumsum(rng.normal(0, 1, (n_grid + 1, 3)), axis=0)[:-1]
        smooth = (smooth - smooth.mean(0)) / smooth.std(0)
        ker = np.exp(-0.5 * (np.arange(-8, 9) * (L / n_grid) / 0.04) ** 2)
        ker /= ker.sum()
        for c in range(3):
            smooth[:, c] = np.convolve(smooth[:, c], ker, mode="same")
        S = np.maximum(base[None, :] + 0.25 * smooth, 0.02)
        frames = [S.copy()]
        for step in range(1, int(round(T / dt_m)) + 1):
            S = diffuse_half(react(diffuse_half(S)))
            if step % every == 0:
                frames.append(S.copy())
        out.append(np.array(frames))
    return out


def add_noise(S, level, rng):
    sig = level * S.std(axis=0)
    return S + rng.normal(0.0, sig, S.shape)


def sigma_hat_series(Sn):
    """Data-driven noise sigma from 3rd differences (robust MAD).

    Var(x_t - 2 x_{t-1} + x_{t-2}) = 6 sigma^2 for white noise; the MAD
    rejects smooth-dynamics contributions.
    """
    d3 = Sn[3:] - 2 * Sn[1:-2] + Sn[:-3]
    mad = np.median(np.abs(d3 - np.median(d3, axis=0)), axis=0) * 1.4826
    return mad / np.sqrt(6.0)


# --------------------------------------------------------------------------- #
#  candidate function library Theta(x): polynomials d<=4 + rational (K grid)
# --------------------------------------------------------------------------- #
def _monomials(max_deg=4, n_var=3):
    """All exponent-count tuples of degree <= max_deg in n_var variables."""
    seqs = []

    def rec(prefix, left, start):
        seqs.append(tuple(prefix))
        if left == 0:
            return
        for i in range(start, n_var):
            prefix.append(i)
            rec(prefix, left - 1, i)
            prefix.pop()

    rec([], max_deg, 0)
    out = []
    for seq in seqs:
        cnt = [0] * n_var
        for i in seq:
            cnt[i] += 1
        out.append(tuple(cnt))
    return out


def build_library_spec():
    spec = [("poly", "1", (0, 0, 0))]
    sym = "uvw"
    for t in _monomials(max_deg=4):
        if sum(t) == 0:
            continue
        name = "".join(f"{sym[i]}{c if c > 1 else ''}"
                       for i, c in enumerate(t) if c > 0)
        spec.append(("poly", name, t))
    for var, vn in enumerate(sym):
        for K in (0.02, 0.2, 2.0):
            spec.append(("rational", f"{vn}/({K:g}+{vn})", (var, K)))
    return spec


LIB_SPEC = build_library_spec()
LIB_NAMES = [s[1] for s in LIB_SPEC]
NP_ = len(LIB_SPEC)


def library_matrix(S):
    S = np.asarray(S)
    u, v, w = S[..., 0], S[..., 1], S[..., 2]
    cols = np.empty(S.shape[:-1] + (NP_,))
    for j, (kind, _, payload) in enumerate(LIB_SPEC):
        if kind == "poly":
            a, b, c = payload
            cols[..., j] = u ** a * v ** b * w ** c
        else:
            var, K = payload
            x = (u, v, w)[var]
            cols[..., j] = x / (K + np.maximum(x, 1e-12))
    return cols


# --------------------------------------------------------------------------- #
#  weak-form (integration-by-parts) split-sample estimator
# --------------------------------------------------------------------------- #
class WeakKernels:
    """Same-span two-grid integration-by-parts kernels.

    The target side samples -phi'(t) x(t) on the EVEN telemetry nodes of a
    window; the library side samples phi(t) Theta(x(t)) on the ODD (interior)
    nodes of the *same* window span.  Both grids therefore integrate over an
    identical time support (the bump vanishes at the window edges), so there
    is no window-offset systematic, while the even/odd node split keeps the
    two noise realizations independent (errors-in-variables decorrelation).
    """

    def __init__(self, m_half, dte):
        self.m_half, self.dte = m_half, dte
        s = np.linspace(0.0, 1.0, m_half + 1)
        self.phi = (s ** 2) * (1.0 - s) ** 2 * 16.0
        # dphi/dt = (dphi/ds) * (1/L) with window length L = m_half * dte
        dphi = np.gradient(self.phi, s) / (m_half * dte)
        # Y side: trapezoid weights on the m+1 even nodes
        wy = np.full(m_half + 1, dte)
        wy[0] = wy[-1] = dte / 2.0
        self.ker_y = -dphi * wy
        # library side: interior odd nodes s = (2j+1)/(2m), open quadrature
        s_odd = (2.0 * np.arange(m_half) + 1.0) / (2.0 * m_half)
        phi_odd = (s_odd ** 2) * (1.0 - s_odd) ** 2 * 16.0
        self.ker_t = phi_odd * dte

    def conv_valid(self, a, ker):
        return np.convolve(a, ker[::-1], mode="valid")

    def features(self, Sn, extra=None):
        even, odd = Sn[0::2], Sn[1::2]
        Y = np.column_stack([self.conv_valid(even[:, k], self.ker_y) for k in range(3)])
        THo = library_matrix(odd)
        if extra is not None:
            THo = np.column_stack([THo, extra])
        THw = np.column_stack([self.conv_valid(THo[:, j], self.ker_t)
                               for j in range(THo.shape[1])])
        n = min(len(Y), len(THw))
        return Y[:n], THw[:n]


STRIDE = 3                           # window subsampling stride
KER0 = WeakKernels(60, 2.0e-3)      # 0-D telemetry kernels
KER_RD = WeakKernels(12, 1.0e-2)    # RD telemetry kernels (0.12 tau window)


def build_weak_dataset(trajs_noisy, kernels, extra_fn=None, row_keep=1,
                       max_rows=None, seed=7):
    """extra_fn(odd, sigma) -> extra columns computed from raw odd samples."""
    rr = np.random.default_rng(seed)
    Ys, THs = [], []
    for S in trajs_noisy:
        sig = NOISE * S.std(axis=0)
        if extra_fn is not None:
            extra = extra_fn(S[1::2], sig)
            Y, TH = kernels.features(S, extra=extra)
        else:
            Y, TH = kernels.features(S)
        Y, TH = Y[::STRIDE][::row_keep], TH[::STRIDE][::row_keep]
        Ys.append(Y)
        THs.append(TH)
    Y = np.vstack(Ys)
    TH = np.vstack(THs)
    if max_rows is not None and len(Y) > max_rows:
        idx = np.sort(np.random.default_rng(3).choice(len(Y), max_rows, replace=False))
        Y, TH = Y[idx], TH[idx]
    return Y, TH


def noise_gram_mc(trajs_noisy, kernels, extra_fn=None, n_mc=4, seed0=9000,
                  row_keep=1, max_points=None, noise_level=NOISE):
    """Monte-Carlo estimate of the library-noise Gram E[N^T N].

    Mirrors the exact feature pipeline (including the Laplacian stage) on
    synthetic noise so that errors-in-variables attenuation is corrected.
    """
    w = NP_ + (2 if extra_fn is not None else 0)
    Nl = np.zeros((w, w))
    for i in range(n_mc):
        chunks = []
        for t_i, S in enumerate(trajs_noisy):
            rm = np.random.default_rng(seed0 + 131 * i + 17 * t_i)
            sig = noise_level * S.std(axis=0)
            eps = np.column_stack([rm.normal(0.0, sig[k], len(S)) for k in range(3)])
            if extra_fn is not None:
                dth_extra = extra_fn(S[1::2] + eps[1::2], sig) - extra_fn(S[1::2], sig)
                dth = np.column_stack([library_matrix(S[1::2] + eps[1::2])
                                       - library_matrix(S[1::2]), dth_extra])
            else:
                dth = library_matrix(S[1::2] + eps[1::2]) - library_matrix(S[1::2])
            dTHw = np.column_stack([kernels.conv_valid(dth[:, j], kernels.ker_t)
                                    for j in range(dth.shape[1])])
            chunks.append(dTHw[::STRIDE][::row_keep])
        Nc = np.vstack(chunks)
        Nl += Nc.T @ Nc
    return Nl / n_mc


class CalibratedSTRidge:
    """Noise-calibrated weak-form STRidge with Pareto model selection."""

    def __init__(self, Y, TH, GN, names, eff_ratio=None):
        self.Y, self.TH, self.GN, self.names = Y, TH, GN, names
        self.norms = np.linalg.norm(TH, axis=0)
        self.THn = TH / self.norms
        self.G = self.THn.T @ self.THn
        self.GNn = GN / np.outer(self.norms, self.norms)
        self.N = len(Y)
        # independent-sample ratio of overlapping weak-form windows
        self.eff_ratio = eff_ratio if eff_ratio is not None             else STRIDE / (KER0.m_half + 1)

    def corrected_solve(self, row, sup):
        """Coefficients (ORIGINAL units) + corrected SSE on sorted support."""
        sup = sorted(sup)
        y = self.Y[:, row]
        TY = self.THn.T @ y
        A = self.G[np.ix_(sup, sup)] - self.GNn[np.ix_(sup, sup)]
        ev, EV = np.linalg.eigh(A)
        ev = np.clip(ev, 1e-6 * ev.max(), None)
        c = np.linalg.solve(EV @ np.diag(ev) @ EV.T, TY[sup])
        r = y - self.THn[:, sup] @ c
        return c / self.norms[sup], float(r @ r)

    def stridge_path(self, row, lam_grid, dtol_grid):
        y = self.Y[:, row]
        TY = self.THn.T @ y
        G = self.THn.T @ self.THn
        n_t = len(self.names)
        out = set()
        for lam in lam_grid:
            xi = np.linalg.solve(G + lam * np.eye(n_t), TY)
            for dtol in dtol_grid:
                x = xi.copy()
                for _ in range(30):
                    small = np.abs(x) < dtol
                    x[small] = 0.0
                    big = ~small
                    if big.sum() < 2:
                        break
                    x[big] = np.linalg.solve(G[np.ix_(big, big)], TY[big])
                out.add(frozenset(np.where(np.abs(x) > 0)[0].tolist()))
        return out

    def elimination(self, row):
        cands = set()
        n_t = len(self.names)
        sup = list(range(n_t))
        while len(sup) > 1:
            cands.add(frozenset(sup))
            best, bj = np.inf, None
            for j in sup:
                s = self.corrected_solve(row, [i for i in sup if i != j])[1]
                if s < best:
                    best, bj = s, j
            sup = [i for i in sup if i != bj]
        sup = []
        while len(sup) < n_t:
            best, bj = np.inf, None
            for j in range(n_t):
                if j in sup:
                    continue
                s = self.corrected_solve(row, sup + [j])[1]
                if s < best:
                    best, bj = s, j
            sup = sup + [bj]
            cands.add(frozenset(sup))
        return cands

    def iht_restarts(self, row, k, restarts=24, seed=1):
        y = self.Y[:, row]
        TY = self.THn.T @ y
        rng = np.random.default_rng(seed)
        n_t = len(self.names)
        best = (np.inf, None)
        for r0 in range(restarts):
            sup = set(np.argsort(-np.abs(TY))[:k].tolist()) if r0 == 0 \
                else set(rng.choice(n_t, size=k, replace=False).tolist())
            for _ in range(60):
                sl = sorted(sup)
                A = self.G[np.ix_(sl, sl)] - self.GNn[np.ix_(sl, sl)]
                ev, EV = np.linalg.eigh(A)
                ev = np.clip(ev, 1e-6 * ev.max(), None)
                cc = EV @ np.diag(1.0 / ev) @ (EV.T @ TY[sl])
                score = np.zeros(n_t)
                score[sl] = np.abs(cc)
                new = set(np.argsort(-score)[:k].tolist())
                if new == sup:
                    break
                sup = new
            sse = self.corrected_solve(row, sup)[1]
            if sse < best[0]:
                best = (sse, frozenset(sorted(sup)))
        return best[1]

    def exhaustive_small_supports(self, row, kmax=5):
        """Best support of EVERY size k <= kmax by corrected SSE (exact).

        Guarantees the per-complexity optimum for small k where the true law
        lives; C(n, kmax) corrected sub-Gram solves at Gram level are cheap.
        """
        from itertools import combinations
        y = self.Y[:, row]
        TY = self.THn.T @ y
        yty = float(y @ y)
        n_t = len(self.names)
        best = {}
        for k in range(1, kmax + 1):
            b_k, s_k = None, np.inf
            for combo in combinations(range(n_t), k):
                idx = list(combo)
                A = self.G[np.ix_(idx, idx)] - self.GNn[np.ix_(idx, idx)]
                ev, EV = np.linalg.eigh(A)
                if ev[0] <= 1e-9 * ev[-1]:
                    continue
                c = EV @ np.diag(1.0 / ev) @ (EV.T @ TY[idx])
                sse = yty - 2.0 * float(c @ TY[idx]) + float(c @ (A @ c))
                # reject ill-conditioned combos: negative expanded SSE is a
                # catastrophic-cancellation artifact of wild coefficients
                if not np.isfinite(sse) or sse < 0:
                    continue
                if np.abs(c).max() > 1e4 * np.sqrt(max(yty / len(TY), 1e-300)):
                    continue
                if sse < s_k:
                    s_k, b_k = sse, frozenset(idx)
            if s_k is not None:
                best[k] = (s_k, b_k)
        return best

    def pareto(self, row, extra_candidates=(), iht_restarts=24, exhaustive_kmax=0):
        """Per-complexity best supports by corrected SSE + noise-floor knee."""
        cands = set()
        cands |= self.stridge_path(row, [1e-5, 1e-4, 1e-3, 1e-2, 1e-1],
                                   [0.02, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.5])
        cands |= self.elimination(row)
        for k in range(1, min(9, len(self.names)) + 1):
            c = self.iht_restarts(row, k, restarts=iht_restarts)
            if c:
                cands.add(c)
        cands |= {frozenset(sorted(c)) for c in extra_candidates if len(c)}
        cands = [c for c in cands if 1 <= len(c) <= len(self.names)]
        per_k, detail = {}, []
        for c in cands:
            coef, sse = self.corrected_solve(row, c)
            detail.append((len(c), sse, c, coef))
            if len(c) not in per_k or sse < per_k[len(c)][0]:
                per_k[len(c)] = (sse, c)
        if exhaustive_kmax > 0:
            for k, (sse, c) in self.exhaustive_small_supports(
                    row, kmax=exhaustive_kmax).items():
                if k not in per_k or sse < per_k[k][0]:
                    per_k[k] = (sse, c)
        ks = sorted(per_k)
        dense = [per_k[k][0] / (self.N - k) for k in ks if k >= min(5, ks[-1])]
        floor = float(np.median(dense)) if dense else min(
            per_k[k][0] / (self.N - k) for k in ks)
        knee = min(k for k in ks if per_k[k][0] / (self.N - k) <= 1.03 * floor)
        return per_k, floor, knee, detail

    def arb_scores(self, rivals):
        """Score rival supports on ARBITRATION telemetry (a CalibratedSTRidge
        built from the 1 % replication with its own 1 % noise-calibrated Gram).

        Selection statistic: effective-sample BIC on the CORRECTED residual.
        Overlapping windows share telemetry, so the residual noise has
        N_eff = N * stride / (window)  independent degrees of freedom; a BIC
        computed with the raw window count under-penalizes dense models by
        that factor and lets them buy victory with noise fitting.
        Also returns the raw per-DOF SSE for reporting.
        Selection is the ONLY role of the arbitration round.
        """
        N = self.N
        n_eff = max(N * self.eff_ratio, 2.0 * len(self.names))
        bics, perdof = [], []
        for row in range(3):
            bic_row, dof_row = {}, {}
            for sup in rivals:
                key = frozenset(sup)
                _, sse = self.corrected_solve(row, key)
                bic_row[key] = n_eff * np.log(max(sse, 1e-300) / N)                     + len(key) * np.log(n_eff)
                dof_row[key] = sse / (N - len(key))
            bics.append(bic_row)
            perdof.append(dof_row)
        return bics, perdof


def prune_collinear(G, names, tol=0.995):
    """Greedy conditioning prune: drop columns with |corr| > tol vs kept.

    G is the Gram of unit-normalized columns, so G[j,k] is the correlation.
    """
    kept, dropped = [], []
    for j in range(len(names)):
        if all(abs(G[j, k]) <= tol for k in kept):
            kept.append(j)
        else:
            dropped.append(names[j])
    return kept, dropped


def top_rivals(detail, n_rival=6):
    ranked = sorted(detail, key=lambda d: d[1])
    rivals = []
    for _, _, c, _ in ranked:
        if all(c != r for r in rivals):
            rivals.append(c)
        if len(rivals) == n_rival:
            break
    return rivals


def deattenuated_transfer_scores(est_fit, est_eval, row, supports):
    """Replicated-experiment transfer score of each support.

    Coefficients are fit on replication A (est_fit, noise-calibrated Gram)
    and scored on independent replication B (est_eval):
        SSE*_B = SSE_B - c^T GN_B c
    i.e. the B residual minus the EXPECTED library-noise propagation
    (errors-in-variables de-attenuation), so supports are not rewarded merely
    for routing weight through low-noise columns.  Returns per-DOF SSE*_B.
    """
    NB = len(est_eval.Y)
    out = {}
    for sup in supports:
        sup_l = sorted(sup)
        c_a, _ = est_fit.corrected_solve(row, sup)
        yB = est_eval.Y[:, row]
        nb = est_eval.norms[sup_l]
        res = yB - est_eval.THn[:, sup_l] @ (c_a * nb)
        sse_b = float(res @ res)
        gn_b = est_eval.GNn[np.ix_(sup_l, sup_l)] * np.outer(nb, nb)
        sse_star = sse_b - float(c_a @ gn_b @ c_a)
        out[frozenset(sup)] = (sse_star / (NB - len(sup_l))
                               if sse_star > 0 else np.inf)
    return out


def pick_parsimonious(dof_row, margin=0.06):
    """Smallest support whose transfer score is within `margin` of the best.

    The de-attenuation uses a Monte-Carlo noise Gram with a few-percent
    relative uncertainty, so score differences inside that band are not
    evidence; Occam's razor decides them.
    """
    finite = {k: v for k, v in dof_row.items() if np.isfinite(v)}
    if not finite:
        return min(dof_row, key=lambda k: len(k))
    best = min(finite.values())
    band = [k for k, v in finite.items() if v <= (1.0 + margin) * best]
    return min(band, key=lambda k: (len(k), finite[k]))


def select_with_arbitration(est, est_arb, est_arb_eval, row, iht_restarts,
                            n_rival=4):
    """Pareto front at operative noise; replicated-experiment arbitration.

    Candidates are generated on both instruments (operative 5 % front and
    the per-complexity optima of the 1 % replication, including exhaustive
    small-support search).  Arbitration fits each rival's coefficients on
    1 % replication A and scores it on INDEPENDENT replication B: noise
    propagated through the library does not transfer between replications,
    so the true law wins on merit rather than by noise absorption.
    """
    per_k, floor, knee, detail = est.pareto(row, iht_restarts=iht_restarts)
    per_k_a, floor_a, knee_a, detail_a = est_arb.pareto(
        row, iht_restarts=max(6, iht_restarts // 2), exhaustive_kmax=5)
    rivals = [per_k[k][1] for k in sorted(per_k)]
    for c in ([per_k_a[k][1] for k in sorted(per_k_a)]
              + top_rivals(detail, n_rival=n_rival)
              + top_rivals(detail_a, n_rival=n_rival)):
        if all(c != r for r in rivals):
            rivals.append(c)
    # replicated-B residuals; ties (shared deterministic bias) are broken by
    # the Bayesian Occam penalty  N log(SSE/N) + k log(N)
    # Every candidate on the parsimony front advances to the replicated-B
    # arbitration.  Complexities beyond the operative Pareto knee (+ margin)
    # fit operative noise without adding structure; the arbitration compares
    # the parsimony front, where noise-calibrated coefficients transfer
    # between independent replications only for the true law.
    k_cap = max(knee, knee_a) + 3
    survivors = [s for s in rivals if len(s) <= k_cap] or list(rivals)
    # stage 2: replicated-B transfer among survivors kills wild-gauge fits
    dof_row = deattenuated_transfer_scores(est_arb, est_arb_eval, row, survivors)
    best_sup = pick_parsimonious(dof_row)
    import os
    if os.environ.get("P12_DEBUG"):
        tk = frozenset(sorted(est_arb.names.index(x) for x in
                              (["v", "uv", "u", "u2"], ["v", "uv", "w"], ["u", "w"])[row]))
        tv = dof_row.get(tk, float("nan"))
        print(f"      [DBG] row{row}: cap={k_cap} truth_dofB={tv:.4g} "
              f"chosen_dofB={dof_row[frozenset(best_sup)]:.4g} "
              f"chosen_k={len(best_sup)} truth_wins={best_sup == tk}")
    coef, sse = est.corrected_solve(row, best_sup)
    return per_k, floor, knee, detail, rivals, dof_row, best_sup, coef, sse, knee_a


# --------------------------------------------------------------------------- #
#  symbolic law assembly (SymPy / LaTeX)
# --------------------------------------------------------------------------- #
SYM_UVW = sympy.symbols("u v w", real=True)


def row_expression(sup, coefs, extra_names=(), extra_coefs=()):
    """sympy expression of one discovered RHS row (original units)."""
    u, v, w = SYM_UVW
    e = sympy.Integer(0)
    sup_s = sorted(sup)
    for pos, j in enumerate(sup_s):
        kind, _, payload = LIB_SPEC[j]
        if kind == "poly":
            a, b, c = payload
            term = u ** a * v ** b * w ** c
        else:
            var, K = payload
            x = (u, v, w)[var]
            term = x / (sympy.Float(K) + x)
        e = e + term * float(coefs[pos])
    for name, cval in zip(extra_names, extra_coefs):
        if name == "lap_u":
            term = sympy.Symbol("lap_u")
        elif name == "lap_v":
            term = sympy.Symbol("lap_v")
        else:
            raise ValueError(name)
        e = e + term * float(cval)
    return sympy.expand(e)


def latex_ode(vn, expr):
    lhs = r"\frac{\mathrm{d}%s}{\mathrm{d}t}" % vn
    return f"{lhs} = {sympy.latex(expr)}"


def lambdified_rows(rows):
    u, v, w = SYM_UVW
    fs = [sympy.lambdify((u, v, w), r, "numpy") for r in rows]

    def f(x):
        return np.array([float(fi(x[0], x[1], x[2])) for fi in fs])

    def f_vec(X):
        return np.column_stack([fi(X[..., 0], X[..., 1], X[..., 2]) for fi in fs])
    return f, f_vec


# --------------------------------------------------------------------------- #
#  Onsager reciprocity gate
# --------------------------------------------------------------------------- #
def jacobian_at(expr, x0):
    """Numeric Jacobian of the discovered RHS (central differences)."""
    f = sympy.lambdify(SYM_UVW, expr, "numpy")

    def F(s):
        return np.array([float(v) for v in f(*s)])

    x0 = np.asarray(x0, float)
    J = np.zeros((3, 3))
    eps = 1e-6
    f0 = F(x0)
    for j in range(3):
        xp = x0.copy()
        xp[j] += eps
        J[:, j] = (F(xp) - f0) / eps
    return J


def steady_state(expr, s0_list):
    """Locate a steady state of the discovered system (multi-start)."""
    from scipy.optimize import fsolve
    f = sympy.lambdify(SYM_UVW, expr, "numpy")

    def F(s):
        return np.array([float(v) for v in f(*s)])

    if isinstance(s0_list, (list, tuple)) and np.asarray(s0_list).ndim == 1:
        s0_list = [s0_list]
    rng = np.random.default_rng(5)
    for s0 in list(s0_list) + [rng.uniform(0.1, 0.9, 3) for _ in range(40)]:
        sol, info, ier, _ = fsolve(F, np.asarray(s0, float), full_output=True,
                                   xtol=1e-12)
        if ier == 1 and np.all(np.isfinite(sol)) and np.all(sol > -0.5):
            res = np.max(np.abs(F(sol)))
            if res < 1e-6:
                return sol
    return None


def onsager_analysis(J):
    """Search an SPD metric G making L = -J G^{-1} as symmetric as possible.

    Onsager reciprocity (L = L^T with L_sym >= 0) is the linear
    near-equilibrium branch of irreversible thermodynamics; the residual
    asymmetry quantifies the distance from equilibrium.
    """
    from scipy.optimize import minimize

    def unpack(p):
        R = np.zeros((3, 3))
        R[0, 0], R[1, 1], R[2, 2] = np.exp(p[0]), np.exp(p[1]), np.exp(p[2])
        R[0, 1], R[0, 2], R[1, 2] = p[3], p[4], p[5]
        return R.T @ R

    def cost(p):
        Li = -J @ np.linalg.inv(unpack(p))
        return np.sum((Li - Li.T) ** 2) / np.sum(Li ** 2)

    best = None
    for seed in range(3):
        rng = np.random.default_rng(seed)
        p0 = np.concatenate([rng.normal(0, 0.3, 3), np.zeros(3)])
        res = minimize(cost, p0, method="Nelder-Mead",
                       options={"maxiter": 4000, "xatol": 1e-10, "fatol": 1e-14})
        if best is None or res.fun < best.fun:
            best = res
    L = -np.linalg.solve(unpack(best.x), J.T).T
    asym = float(np.linalg.norm(L - L.T) / np.linalg.norm(L))
    eig_sym = float(np.linalg.eigvalsh(0.5 * (L + L.T)).min())
    return asym, eig_sym, float(best.fun)


# --------------------------------------------------------------------------- #
#  MODULE 12C — Hamiltonian core & Lyapunov functional (continuous NN)
# --------------------------------------------------------------------------- #
def train_scalar_functional(f_rhs_vec, mode, main_samples, attractor_samples,
                            n_steps=2500, lr=2e-3, seed=0, pairs=None):
    """Learn a scalar functional on phase space with a continuous MLP.

    mode='hamiltonian': dH/dt = gradH . f ~ 0 on the attractor band
                        (invariant-manifold core).
    mode='lyapunov':    V = softplus(net) >= 0; dV/dt <= 0 off-attractor
                        (entropy-dissipation funnel); V -> 0 and dV/dt -> 0
                        on the attractor (non-equilibrium steady balance).
    Input scaling is folded into the network (buffer 'scale').
    """
    import torch
    torch.manual_seed(seed)
    torch.set_num_threads(4)
    dt64 = torch.float64

    class Net(torch.nn.Module):
        def __init__(self, softplus_out):
            super().__init__()
            self.core = torch.nn.Sequential(
                torch.nn.Linear(3, 64), torch.nn.Softplus(),
                torch.nn.Linear(64, 64), torch.nn.Softplus(),
                torch.nn.Linear(64, 1))
            self.softplus_out = softplus_out
            self.register_buffer("scale", torch.ones(3, dtype=dt64))

        def forward(self, x):
            z = self.core(x / self.scale).squeeze(-1)
            return torch.nn.functional.softplus(z) if self.softplus_out else z

    net = Net(softplus_out=(mode == "lyapunov")).to(dt64)
    scale = main_samples.std(axis=0).clip(1e-6)
    net.scale.copy_(torch.tensor(scale, dtype=dt64))
    Xm = torch.tensor(main_samples, dtype=dt64)
    Xa = torch.tensor(attractor_samples, dtype=dt64)
    Xp = torch.tensor(pairs[0], dtype=dt64) if pairs is not None else None
    Xq = torch.tensor(pairs[1], dtype=dt64) if pairs is not None else None
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, n_steps)
    bs = min(4096, len(Xm), len(Xa))

    def f_torch(Xb):
        return torch.tensor(np.asarray(f_rhs_vec(Xb.detach().numpy())), dtype=dt64)

    for step in range(n_steps):
        Xb = Xm[torch.randint(0, len(Xm), (bs,))].clone().requires_grad_(True)
        Hb = net(Xb)
        gH = torch.autograd.grad(Hb.sum(), Xb, create_graph=True)[0]
        Fb = f_torch(Xb)
        dH = (gH * Fb).sum(-1)
        if mode == "hamiltonian":
            scale_n = (gH.norm(dim=-1) * Fb.norm(dim=-1)).clamp_min(1e-6)
            collapse = torch.relu(0.25 - Hb.std()) ** 2        # anti-collapse
            loss = ((dH / scale_n) ** 2).mean() + collapse \
                + 1e-4 * (gH.norm(dim=-1) ** 2).mean()
            if Xp is not None:                                  # endpoint consistency
                pi = torch.randint(0, len(Xp), (bs,))
                dH_pair = (net(Xp[pi]) - net(Xq[pi])).squeeze(-1)
                h_rng = float(Hb.std() * 4.0 + 1e-6)
                loss = loss + 0.5 * ((dH_pair / h_rng) ** 2).mean()
        else:
            Xab = Xa[torch.randint(0, len(Xa), (bs,))].clone().requires_grad_(True)
            Hab = net(Xab)
            gA = torch.autograd.grad(Hab.sum(), Xab, create_graph=True)[0]
            dA = (gA * f_torch(Xab)).sum(-1)
            off_pen = torch.relu(dH + 1e-3) ** 2               # dV/dt <= 0
            height = torch.relu(0.02 - Hb) ** 2                # anti-collapse
            loss = off_pen.mean() + height.mean() \
                + (dA ** 2).mean() + (Hab ** 2).mean() \
                + 1e-5 * (gA.norm(dim=-1) ** 2).mean()
        opt.zero_grad()
        loss.backward()
        opt.step()
        sched.step()
        if step % 500 == 0 or step == n_steps - 1:
            print(f"      {elapsed()} step {step:5d}  loss = {float(loss):.3e}",
                  flush=True)
    return net


def nn_value_and_grad(net, X):
    import torch
    Xt = torch.tensor(np.asarray(X, dtype=np.float64), requires_grad=True)
    out = net(Xt).squeeze(-1)
    g = torch.autograd.grad(out.sum(), Xt)[0]
    return out.detach().numpy(), g.detach().numpy()


def lle_benettin(f, s0, T=600.0, dt=1e-3, renorm_every=0.1, seed=0):
    """Lyapunov spectrum by Benettin QR with the tangent flow integrated
    alongside the state (augmented RK4, frequent re-orthonormalization so the
    strongly contracting Oregonator directions neither under- nor overflow).
    Returns (lambda_1, lambda_2, lambda_3) in 1/tau.
    """
    def jac(s, f0):
        J = np.zeros((3, 3))
        eps = 1e-7
        for j in range(3):
            sp = s.copy()
            sp[j] += eps
            J[:, j] = (f(sp) - f0) / eps
        return J

    def rhs_aug(z):
        s, D = z[:3], z[3:].reshape(3, 3)
        f0 = f(s)
        return np.concatenate([f0, (jac(s, f0) @ D).ravel()])

    z = np.concatenate([np.asarray(s0, float), np.eye(3).ravel()])
    acc = np.zeros(3)
    n_ren = int(round(renorm_every / dt))
    n_blocks = int(round(T / renorm_every))
    for _ in range(n_blocks):
        for _ in range(n_ren):
            k1 = rhs_aug(z)
            k2 = rhs_aug(z + 0.5 * dt * k1)
            k3 = rhs_aug(z + 0.5 * dt * k2)
            k4 = rhs_aug(z + dt * k3)
            z = z + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        Q, R = np.linalg.qr(z[3:].reshape(3, 3))
        acc += np.log(np.abs(np.diag(R)) + 1e-300)
        z[3:] = Q.ravel()
    return acc / (n_blocks * renorm_every)


def orbit_period(t, x, thr=0.5):
    """Mean period from successive upward crossings of a level."""
    cross = np.where((x[1:] > thr) & (x[:-1] <= thr))[0]
    if len(cross) >= 3:
        return float(np.mean(np.diff(t[cross])))
    return float("nan")


# --------------------------------------------------------------------------- #
#  RD Laplacian pipeline (spectral smoothing + Wiener deconvolution)
# --------------------------------------------------------------------------- #
def rd_laplacian(field, sigma_uv, n_grid, L):
    """field: (..., n_grid) raw noisy u (or v) channel. Returns lap (same shape)."""
    sig = max(sigma_uv, 1e-9)
    k = np.fft.rfftfreq(n_grid, d=L / n_grid) * 2.0 * np.pi
    sig_f = 1.6 * (L / n_grid)
    damp = np.exp(-0.5 * (k * sig_f) ** 2)
    n_f = field.shape[-2]
    U = np.fft.rfft(field, axis=-1)
    Us = U * damp[None, :]
    P = (np.abs(Us) ** 2).mean(axis=0) / n_f            # smoothed signal power /k
    P_bin = np.copy(P)
    for i in range(0, len(k), 8):
        sl = slice(i, min(i + 8, len(k)))
        P_bin[sl] = P[sl].mean()
    N_pow = (sig ** 2) * float(np.mean(damp ** 2))
    Hw = P_bin / (P_bin + N_pow)
    Uc = Us * (np.exp(0.5 * (k[None, :] * sig_f) ** 2) * Hw[None, :])
    return np.fft.irfft(Uc * (-k[None, :] ** 2), n=n_grid, axis=-1)


def make_rd_extra_fn(laps_odd):
    """extra columns for one grid point: laps_odd[:, p, :] (aligned to odd)."""
    def extra(odd, sig):
        return laps_odd
    return extra


# --------------------------------------------------------------------------- #
#  C++ kernel export + verification
# --------------------------------------------------------------------------- #
def export_cpp_kernel(path, expr_rows, netH, netV):
    import torch

    def net_code(net, tag, softplus):
        sd = {k: v.detach().numpy() for k, v in net.state_dict().items()}
        W1, B1 = sd["core.0.weight"], sd["core.0.bias"]
        W2, B2 = sd["core.2.weight"], sd["core.2.bias"]
        W3, B3 = sd["core.4.weight"], sd["core.4.bias"]
        scale = net.scale.detach().numpy()
        lines = [f"// ---- discovered {tag} (continuous MLP, softplus hidden) ----",
                 f"static const double {tag}_W1[64][3] = {{"]
        for i in range(64):
            lines.append("    {" + ", ".join(
                f"{W1[i][j] / scale[j]:.17e}" for j in range(3)) + "},")
        lines.append("};")
        lines.append(f"static const double {tag}_B1[64] = {{" +
                     ", ".join(f"{x:.17e}" for x in B1) + "};")
        lines.append(f"static const double {tag}_W2[64][64] = {{")
        for i in range(64):
            lines.append("    {" + ", ".join(f"{x:.17e}" for x in W2[i]) + "},")
        lines.append("};")
        lines.append(f"static const double {tag}_B2[64] = {{" +
                     ", ".join(f"{x:.17e}" for x in B2) + "};")
        lines.append(f"static const double {tag}_W3[64] = {{" +
                     ", ".join(f"{x:.17e}" for x in W3[0]) + "};")
        lines.append(f"static const double {tag}_B3 = {float(B3[0]):.17e};")
        tail = "acc" if not softplus else "(acc > 30.0 ? acc : std::log1p(std::exp(acc)))"
        lines.append(f"""
inline double {tag}_eval(const double x[3]) {{
    double h1[64], h2[64];
    for (int i = 0; i < 64; ++i) {{
        double z = {tag}_B1[i];
        for (int j = 0; j < 3; ++j) z += {tag}_W1[i][j] * x[j];
        h1[i] = z > 30.0 ? z : std::log1p(std::exp(z));
    }}
    double acc = {tag}_B3;
    for (int i = 0; i < 64; ++i) {{
        double z = {tag}_B2[i];
        for (int j = 0; j < 64; ++j) z += {tag}_W2[i][j] * h1[j];
        h2[i] = z > 30.0 ? z : std::log1p(std::exp(z));
        acc += {tag}_W3[i] * h2[i];
    }}
    return {tail};
}}""")
        return "\n".join(lines)

    ccode = [sympy.ccode(e) for e in expr_rows]
    hdr = f"""// phase12_kernels.hpp — PHASE 12 discovered-law execution kernels (C++17)
// Autonomously discovered by run_phase12_hamiltonian_law_discovery.py
// Discovered symbolic ODE RHS (dimensionless Oregonator variables u, v, w):
//   du/dt = {sympy.sstr(expr_rows[0])}
//   dv/dt = {sympy.sstr(expr_rows[1])}
//   dw/dt = {sympy.sstr(expr_rows[2])}
// plus the discovered conserved Hamiltonian core H(u,v,w) and the
// dissipative Lyapunov functional V(u,v,w) >= 0 as embedded neural kernels.
#pragma once
#include <cmath>

namespace phase12 {{

static const char* LAW_PROVENANCE =
    "autonomous symbolic induction from 5% noisy telemetry; weak-form STRidge "
    "+ noise-calibrated Gram + intervention/arbitration epistemic loop";

inline void rhs_discovered(const double s[3], double ds[3]) {{
    const double u = s[0], v = s[1], w = s[2];
    (void)u; (void)v; (void)w;
    ds[0] = {ccode[0]};
    ds[1] = {ccode[1]};
    ds[2] = {ccode[2]};
}}

{net_code(netH, "H_func", softplus=False)}

{net_code(netV, "V_func", softplus=True)}

}} // namespace phase12
"""
    path.write_text(hdr, encoding="utf-8")


def cpp_demo_and_verify(hpp_path, expr_rows, netH, netV):
    if shutil.which("g++") is None:
        print(f"   {elapsed()} g++ not found — skipping C++ execution check")
        return None
    f_py, _ = lambdified_rows(expr_rows)
    test_pts = np.array([[0.6, 0.6, 0.6], [0.9, 2.5, 0.8], [0.2, 0.1, 0.3],
                         [0.95, 12.0, 0.9], [0.5, 1.0, 0.5]])
    ref = np.array([f_py(p) for p in test_pts])
    h_ref, _ = nn_value_and_grad(netH, test_pts)
    v_ref, _ = nn_value_and_grad(netV, test_pts)
    pts_init = ", ".join("{" + ", ".join(repr(float(x)) for x in p) + "}"
                         for p in test_pts)
    cpp = f"""
#include "{hpp_path.name}"
#include <cstdio>
using namespace phase12;
int main() {{
    double pts[5][3] = {{ {pts_init} }};
    for (int i = 0; i < 5; ++i) {{
        double ds[3];
        rhs_discovered(pts[i], ds);
        printf("%.17e %.17e %.17e %.17e %.17e\\n",
               ds[0], ds[1], ds[2], H_func_eval(pts[i]), V_func_eval(pts[i]));
    }}
    return 0;
}}
"""
    cpp_path = RES / "_phase12_kernel_test.cpp"
    cpp_path.write_text(cpp)
    exe = RES / "_phase12_kernel_test.exe"
    r = subprocess.run(["g++", "-O2", "-std=c++17", "-o", str(exe), str(cpp_path)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(f"   {elapsed()} g++ compile failed:\n{r.stderr[:600]}")
        return None
    stdout = subprocess.run([str(exe)], capture_output=True, text=True).stdout
    got = np.array([[float(x) for x in ln.split()] for ln in stdout.strip().splitlines()])
    res = {"ode": float(np.abs(got[:, :3] - ref).max()),
           "H": float(np.abs(got[:, 3] - h_ref).max()),
           "V": float(np.abs(got[:, 4] - v_ref).max())}
    print(f"   {elapsed()} C++ kernel vs Python: max|d_ode| = {res['ode']:.2e}, "
          f"max|dH| = {res['H']:.2e}, max|dV| = {res['V']:.2e}")
    return res


# --------------------------------------------------------------------------- #
#  figures (300 DPI)
# --------------------------------------------------------------------------- #
def _cycle_collection(ax, traj, cmap="viridis", lw=0.8):
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    pts = traj.reshape(-1, 1, 3)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    tseg = np.linspace(0, 1, len(segs))
    lc = Line3DCollection(segs, cmap=cmap, array=tseg, linewidths=lw)
    ax.add_collection3d(lc)
    return lc


def fig1(sol_true, sol_disc, per_true, per_disc, rd_field, t_window=(40.0, 54.0)):
    if isinstance(sol_true, dict):
        y_true, y_disc, t_grid = sol_true["y"], sol_disc["y"], sol_true["t"]
    else:
        y_true, y_disc, t_grid = sol_true.y, sol_disc.y, sol_true.t
    fig = plt.figure(figsize=(12.5, 9.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], hspace=0.32, wspace=0.24)

    ax = fig.add_subplot(gs[0, 0], projection="3d")
    _cycle_collection(ax, y_true.T[::20])
    ax.set_xlabel("u"); ax.set_ylabel("v"); ax.set_zlabel("w")
    ax.set_title(f"(a) True Oregonator attractor — limit cycle, "
                 f"T = {per_true:.3f} $\\tau$")
    ax.view_init(elev=22, azim=-58)

    ax = fig.add_subplot(gs[0, 1], projection="3d")
    _cycle_collection(ax, y_disc.T[::20], cmap="plasma")
    rel = abs(per_disc - per_true) / per_true
    ax.set_xlabel("u"); ax.set_ylabel("v"); ax.set_zlabel("w")
    ax.set_title(f"(b) Symbolically discovered attractor — "
                 f"period err {rel:.2%}")
    ax.view_init(elev=22, azim=-58)

    ax = fig.add_subplot(gs[1, 0])
    m = (t_grid >= t_window[0]) & (t_grid <= t_window[1])
    ax.plot(t_grid[m], y_true[0][m], color=C_MAIN, lw=1.6, label="true $u(t)$")
    ax.plot(t_grid[m], y_disc[0][m], color=C_ACC, lw=1.1, ls="--",
            label="discovered $u(t)$")
    ax.set_xlabel("$t$ [$\\tau$]"); ax.set_ylabel("$u$")
    ax.set_title("(c) Time-domain structural fidelity (2 cycles)")
    ax.legend(loc="upper right")
    axins = ax.inset_axes([0.34, 0.62, 0.30, 0.26])
    axins.plot(t_grid[m], y_true[0][m] - y_disc[0][m],
               color=C_GREY, lw=0.9)
    axins.set_title("$\\Delta u(t)$", fontsize=8)
    axins.tick_params(labelsize=7)

    ax = fig.add_subplot(gs[1, 1])
    im = ax.imshow(rd_field.T, aspect="auto", cmap="inferno", origin="lower",
                   extent=[0, rd_field.shape[0] * 5e-3, 0, 0.6])
    ax.set_xlabel("$t$ [$\\tau$]"); ax.set_ylabel("$x$ [l.u.]")
    ax.set_title("(d) 1-D reaction-diffusion field $u(x,t)$ — spatial telemetry")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    fig.suptitle("Phase 12 / Fig. 1 — True vs. symbolically discovered non-equilibrium attractor",
                 fontsize=12, y=0.98)
    fig.savefig(FIG / "fig1_spatiotemporal_reconstruction.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)


def fig2(pareto_data, selections, true_k):
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.6))
    ax = axes[0]
    colors = {"u": C_MAIN, "v": C_ACC, "w": C_MID}
    for row, vn in enumerate(("u", "v", "w")):
        pd = pareto_data[row]
        ks = sorted(int(k) for k in pd["per_k"])
        sse = [pd["per_k"][str(k)][0] / selections[row]["N"] for k in ks]
        ax.plot(ks, sse, "-o", color=colors[vn], ms=3.5, lw=1.0, alpha=0.75,
                label=f"$\\dot {vn}$: Pareto front")
        knee = pd["knee"]
        ax.plot([knee], [pd["per_k"][str(knee)][0] / selections[row]["N"]],
                marker="*", ms=15, color=colors[vn], mec="k", mew=0.4, ls="none")
        ax.axvline(true_k[vn], color=colors[vn], ls=":", lw=0.9, alpha=0.7)
    ax.set_yscale("log")
    ax.set_xlabel("model complexity  $k$  (number of symbolic terms)")
    ax.set_ylabel("corrected residual SSE / DOF")
    ax.set_title("(a) Occam Pareto front — star = knee")
    ax.legend(fontsize=7.5)

    ax = axes[1]
    w = 0.36
    band_lo, band_hi = 1.0 / 1.06, 1.06
    ax.axhspan(band_lo, band_hi, color=C_MID, alpha=0.18,
               label="de-attenuation calibration band (Occam decides)" if True else None)
    ax.axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.6)
    for i, vn in enumerate(("u", "v", "w")):
        sel = selections[i]
        r_5 = sel["arb5_runnerup"] / sel["arb5_chosen"]
        r_1 = sel["arb1_runnerup"] / sel["arb1_chosen"]
        ax.bar(i - w / 2, max(r_5, 1e-3), w, color=C_GREY, alpha=0.85,
               label="$\\sigma$ = 5 % (operative)" if i == 0 else None)
        ax.bar(i + w / 2, max(r_1, 1e-3), w, color=C_GOLD,
               label="$\\sigma$ = 1 % (arbitration)" if i == 0 else None)
        ax.text(i + w / 2, max(r_1, 1e-3) * 1.4, f"{r_1:.2f}x",
                ha="center", fontsize=8)
    ax.set_yscale("log")
    ax.set_ylim(3e-2, 3e6)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["$\\dot u$", "$\\dot v$", "$\\dot w$"])
    ax.set_ylabel("runner-up / chosen  residual ratio")
    ax.set_title("(b) Intervention + precision arbitration:\n"
                 "a tie inside the band is settled by Occam")
    ax.legend(fontsize=7, loc="upper left")

    ax = axes[2]
    mk = {"u": "o", "v": "s", "w": "^"}
    for row, vn in enumerate(("u", "v", "w")):
        got = selections[row]["got"]
        tru = selections[row]["true"]
        xs = [tru[n] for n in got]
        ys = [got[n] for n in got]
        ax.scatter(xs, ys, marker=mk[vn], s=28, color=colors[vn], alpha=0.85,
                   label=f"$\\dot {vn}$  (max err {selections[row]['max_err']:.1%})")
    lims = [3e-3, 4e2]
    ax.plot(lims, lims, "k-", lw=0.8)
    ax.fill_between(lims, [l * 0.9 for l in lims], [l * 1.1 for l in lims],
                    color=C_MID, alpha=0.12)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel("true coefficient"); ax.set_ylabel("discovered coefficient")
    ax.set_title("(c) Coefficient recovery ($\\pm$10 % band)")
    ax.legend(fontsize=7.5)
    fig.suptitle("Phase 12 / Fig. 2 — Symbolic parsimony vs. predictive loss: "
                 "pinpointing the uniquely correct law", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIG / "fig2_symbolic_pareto_complexity.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)


def fig3(netV, netH, f_disc, f_disc_vec, long_traj, kicked_trajs):
    import torch  # noqa
    # slice window: attractor neighbourhood in v (the V-funnel region)
    u_g = np.linspace(0.02, 1.02, 80)
    v_g = np.linspace(0.05, 4.0, 90)
    UU, VV = np.meshgrid(u_g, v_g)
    w_bar = float(np.percentile(long_traj[:, 2], 60))
    WW = np.full_like(UU, w_bar)
    X = np.column_stack([UU.ravel(), VV.ravel(), WW.ravel()])
    Vv, gV = nn_value_and_grad(netV, X)
    F = f_disc_vec(X)
    dV = (gV * F).sum(-1).reshape(UU.shape)
    Vsurf = Vv.reshape(UU.shape)

    fig = plt.figure(figsize=(12.8, 9.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], hspace=0.32, wspace=0.22)

    ax = fig.add_subplot(gs[0, 0], projection="3d")
    vmin, vmax = float(np.percentile(Vsurf, 5)), float(np.percentile(Vsurf, 97))
    ax.plot_surface(UU, VV, Vsurf, cmap="magma",
                    rstride=2, cstride=2, linewidth=0, antialiased=True,
                    alpha=0.94, vmin=vmin, vmax=vmax)
    ax.set_zlim(vmin, vmax)
    cyc = long_traj[::40]
    z0 = vmin
    ax.plot(cyc[:, 0], cyc[:, 1], np.full(len(cyc), z0), color=C_MID, lw=1.8,
            label="limit cycle (projected)")
    ax.set_xlabel("u"); ax.set_ylabel("v"); ax.set_zlabel("$V(u,v)$")
    ax.set_title(f"(a) Discovered Lyapunov funnel $V \\geq 0$ "
                 f"(slice $w = 0.06$)")
    ax.view_init(elev=30, azim=-52)
    ax.legend(loc="upper left", fontsize=7.5)

    ax = fig.add_subplot(gs[0, 1])
    vmin_d, vmax_d = float(np.percentile(dV, 6)), float(np.percentile(dV, 94))
    pc = ax.pcolormesh(UU, VV, dV, cmap="RdBu_r", shading="auto",
                       norm=TwoSlopeNorm(vcenter=0.0, vmin=vmin_d, vmax=vmax_d))
    ax.contour(UU, VV, dV, levels=[0.0], colors="k", linewidths=1.3)
    ax.plot(long_traj[:, 0], long_traj[:, 1], color=C_MID, lw=1.5,
            label="limit cycle: $\\langle\\dot V\\rangle \\to 0$ (NESS balance)")
    ax.set_xlabel("u"); ax.set_ylabel("v")
    ax.set_title("(b) Entropy-dissipation rate $\\dot V = \\nabla V \\!\\cdot\\! f$ "
                 "$\\leq 0$ off-attractor")
    fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.02)
    ax.legend(loc="upper right", fontsize=7.5)

    ax = fig.add_subplot(gs[1, 0])
    post_segments = kicked_trajs[1::2][:4]
    for i, seg in enumerate(post_segments):
        seg = seg[:int(len(seg) * 0.45)]
        Vt, _ = nn_value_and_grad(netV, seg[::10])
        t_kick = np.arange(len(Vt)) * 0.02
        ax.plot(t_kick, Vt, lw=1.2, alpha=0.85,
                label=f"post-kick transient {i+1}" if i < 3 else None)
    ax.set_yscale("log")
    ax.set_xlabel("$t$ since kick [$\\tau$]")
    ax.set_ylabel("$V(u,v,w)$")
    ax.set_title("(c) Monotone entropy-dissipation descent")
    ax.legend(fontsize=7)

    ax = fig.add_subplot(gs[1, 1])
    Ht, _ = nn_value_and_grad(netH, long_traj[::20])
    tgrid = np.arange(len(Ht)) * 0.02
    m = (tgrid >= 20.0) & (tgrid <= 41.0)
    Hbar = float(Ht.mean())
    ax.plot(tgrid[m], Ht[m] - Hbar, color=C_MAIN, lw=1.2)
    per = 10.12
    n_per = float(np.max(np.abs(Ht[m] - Ht[m][0])) / (Ht.max() - Ht.min() + 1e-12))
    ax.set_xlabel("$t$ [$\\tau$]"); ax.set_ylabel("$H - \\bar H$")
    ax.set_title("(d) Hamiltonian core along the attractor: "
                 f"drift = {n_per:.2f} range/period")
    fig.suptitle("Phase 12 / Fig. 3 — Discovered Lyapunov entropy-dissipation "
                 "surface & conserved Hamiltonian core", fontsize=12)
    fig.savefig(FIG / "fig3_lyapunov_entropy_descent.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
#  main pipeline
# --------------------------------------------------------------------------- #
def render_figures_only():
    """Re-render the publication figures from the saved artefacts."""
    import torch
    banner("PHASE 12 — figure re-render (from saved artefacts)")
    import json as _json
    fd = np.load(RES / "_figdata.npz")
    meta = _json.loads((RES / "_figmeta.json").read_text())
    kicks = [fd[k] for k in sorted(fd.files) if k.startswith("kick")]
    state = torch.load(RES / "_nets.pt", weights_only=True)

    class Net(torch.nn.Module):
        def __init__(sp, softplus_out, scale):
            super().__init__()
            sp.core = torch.nn.Sequential(
                torch.nn.Linear(3, 64), torch.nn.Softplus(),
                torch.nn.Linear(64, 64), torch.nn.Softplus(),
                torch.nn.Linear(64, 1))
            sp.softplus_out = softplus_out
            sp.register_buffer("scale", torch.tensor(scale, dtype=torch.float64))

        def forward(sp, x):
            z = sp.core(x / sp.scale).squeeze(-1)
            return torch.nn.functional.softplus(z) if sp.softplus_out else z

    netH = Net(False, state["scaleH"]).to(torch.float64)
    netV = Net(True, state["scaleV"]).to(torch.float64)
    netH.load_state_dict(state["netH"])
    netV.load_state_dict(state["netV"])
    sol_true_y, sol_disc_y = fd["sol_true"], fd["sol_disc"]
    per_true, per_disc = meta["per_true"], meta["per_disc"]

    # rebuild the discovered vector field from the saved law
    law = (RES / "discovered_laws_sympy.txt").read_text().splitlines()
    rows = [sympy.sympify(ln.split("= ", 1)[1]) for ln in law if ln.startswith("d")]
    f_disc, f_disc_vec = lambdified_rows(rows)
    fig1({"t": fd["t_grid"], "y": sol_true_y.T},
         {"t": fd["t_grid"], "y": sol_disc_y.T},
         per_true, per_disc, fd["rd_field"])
    fig2(meta["pareto"], meta["selections"], meta["true_k"])
    fig3(netV, netH, f_disc, f_disc_vec, fd["long_traj"], kicks)
    print(f"   {elapsed()} figures re-rendered -> figures_phase12/")


def main():
    global QUICK
    ap = argparse.ArgumentParser(
        description="Phase 12: autonomous Hamiltonian law discovery")
    ap.add_argument("--quick", action="store_true", help="reduced-data smoke run")
    ap.add_argument("--fig_only", action="store_true",
                    help="re-render figures from the saved artefacts")
    args = ap.parse_args()
    QUICK = args.quick
    if args.fig_only:
        render_figures_only()
        return
    iht_restarts = 10 if QUICK else 28
    n_transient = 3 if QUICK else 8
    n_kick_rep = 1 if QUICK else 2
    n_sto = 2 if QUICK else 3
    net_steps = 1500 if QUICK else 4000
    lle_T = 60.0 if QUICK else 300.0

    results = {"phase": 12, "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")}
    rng = np.random.default_rng(11)

    # ================================================================== 12A
    banner("MODULE 12A — non-equilibrium dissipative simulation & telemetry injection")
    print(f"   {elapsed()} passive 0-D telemetry (limit cycle, f = {F_OSC}) ...")
    attractor_trajs = [simulate_0d(np.array([0.6, 0.6, 0.6]), 80.0),
                       simulate_0d(np.array([0.9, 0.2, 0.4]), 80.0)]
    transient_trajs = []
    for _ in range(n_transient):
        s0 = np.array([rng.uniform(0.05, 1.1), rng.uniform(0.1, 8.0), rng.uniform(0.05, 1.2)])
        transient_trajs.append(simulate_0d(s0, 30.0))
    print(f"   {elapsed()} factorial impulse interventions "
          f"(single-variable kicks, basin-clipped) ...")
    base = simulate_0d(np.array([0.6, 0.6, 0.6]), 5.0)[-1]
    kick_plan = {0: (0.35, -0.35), 1: (1.5, -0.35), 2: (0.45, -0.45)}
    kicked_trajs = []
    for var, amps in kick_plan.items():
        for a in amps:
            for _ in range(n_kick_rep):
                s0 = base + rng.normal(0, 0.04, 3)
                tk = 5.0 + 3.0 * rng.uniform()
                kicked_trajs.extend(kicked_trajectory(s0, var, a, tk, T=25.0))
    train_trajs = attractor_trajs + transient_trajs + kicked_trajs

    Sn_train = [add_noise(S, NOISE, np.random.default_rng(100 + i))
                for i, S in enumerate(train_trajs)]
    sh = np.median([sigma_hat_series(S) / (NOISE * S.std(axis=0)) for S in Sn_train],
                   axis=0)
    print(f"   {elapsed()} blind noise calibration sigma_hat/sigma_true = "
          f"{np.round(sh, 3)}  (1.000 = perfect)")
    arb_trajs = [add_noise(S, NOISE_ARB, np.random.default_rng(500 + i))
                 for i, S in enumerate(train_trajs)]
    arb2_trajs = [add_noise(S, NOISE_ARB, np.random.default_rng(600 + i))
                  for i, S in enumerate(train_trajs)]

    print(f"   {elapsed()} 1-D reaction-diffusion Oregonator (spatial telemetry) ...")
    rd_fields = simulate_rd(T=25.0, seed=3, n_ic=1 if QUICK else 2)
    rd_noisy = []
    for i, Fld in enumerate(rd_fields):
        sig = NOISE * Fld.reshape(-1, 3).std(axis=0)
        rn = np.random.default_rng(77 + i)
        rd_noisy.append(Fld + rn.normal(0.0, sig, Fld.shape))
    rd_arb = []
    rd_arb2 = []
    for i, Fld in enumerate(rd_fields):
        sig = NOISE_ARB * Fld.reshape(-1, 3).std(axis=0)
        rn = np.random.default_rng(177 + i)
        rd_arb.append(Fld + rn.normal(0.0, sig, Fld.shape))
        rn2 = np.random.default_rng(277 + i)
        rd_arb2.append(Fld + rn2.normal(0.0, sig, Fld.shape))
    np.savez_compressed(RES / "telemetry_rd.npz",
                        **{f"rd{i}": F for i, F in enumerate(rd_noisy)})
    results["12A"] = {
        "n_train_trajectories": len(train_trajs),
        "n_kick_experiments": len(kicked_trajs),
        "noise_sigma": NOISE,
        "sigma_hat_over_true": [float(x) for x in sh],
        "rd_grid": 192, "rd_L": 0.6,
        "n_rd_frames": int(rd_noisy[0].shape[0]),
    }

    # ================================================================== 12B
    banner("MODULE 12B — physics-informed sparse symbolic regression "
           "(weak-form STRidge + Onsager gate)")
    print(f"   {elapsed()} weak-form features: library = {NP_} candidate terms "
          f"(poly d<=4 + rational K-grid), Simpson split-sample windows ...")
    Y, TH = build_weak_dataset(Sn_train, KER0)
    GN = noise_gram_mc(Sn_train, KER0)
    pre = CalibratedSTRidge(Y, TH, GN, LIB_NAMES)
    kept, dropped = prune_collinear(pre.G, LIB_NAMES)
    print(f"   {elapsed()} conditioning prune: dropped {len(dropped)} "
          f"near-collinear columns {dropped}")
    TH = TH[:, kept]
    GN = GN[np.ix_(kept, kept)]
    NAMES0 = [LIB_NAMES[i] for i in kept]
    est = CalibratedSTRidge(Y, TH, GN, NAMES0)
    np.savez_compressed(RES / "_diag_est.npz", Y=Y, TH=TH, GN=GN,
                        kept=np.array(kept))
    Y_a, TH_a = build_weak_dataset(arb_trajs, KER0, seed=8)
    GN_a = noise_gram_mc(arb_trajs, KER0, noise_level=NOISE_ARB, seed0=8000)
    est_arb = CalibratedSTRidge(Y_a, TH_a[:, kept], GN_a[np.ix_(kept, kept)], NAMES0)
    Y_b, TH_b = build_weak_dataset(arb2_trajs, KER0, seed=9)
    GN_b = noise_gram_mc(arb2_trajs, KER0, noise_level=NOISE_ARB, seed0=8100)
    est_arb_eval = CalibratedSTRidge(Y_b, TH_b[:, kept], GN_b[np.ix_(kept, kept)], NAMES0)
    np.savez_compressed(RES / "_diag_arb.npz", Y=Y_a, TH=TH_a[:, kept],
                        GN=GN_a[np.ix_(kept, kept)], kept=np.array(kept))
    print(f"   {elapsed()} operative windows: {len(Y)} | library after prune: "
          f"{len(NAMES0)} | cond(G) = {np.linalg.cond(est.G):.3g} | "
          f"||GN||/||G|| = {np.linalg.norm(est.GNn) / np.linalg.norm(est.G):.4f}")

    VAR_NAMES = ("u", "v", "w")
    LIBR = NAMES0
    TRUE_SUP = {0: ["v", "uv", "u", "u2"], 1: ["v", "uv", "w"], 2: ["u", "w"]}
    TRUE_COEF = {0: {"v": Q / E1, "uv": -1 / E1, "u": 1 / E1, "u2": -1 / E1},
                 1: {"v": -Q / E2, "uv": -1 / E2, "w": F_OSC / E2},
                 2: {"u": 1.0, "w": -1.0}}
    TRUE_K = {0: 4, 1: 3, 2: 2}

    name2idx = {n: i for i, n in enumerate(LIB_NAMES)}
    supports, coefs, pareto_data, selections = [], [], [], []
    for row, vn in enumerate(VAR_NAMES):
        print(f"   {elapsed()} Pareto + arbitration selection for d{vn}/dt ...")
        per_k, floor, knee, detail, rivals, dof_row, best_sup, coef, sse, knee_a = \
            select_with_arbitration(est, est_arb, est_arb_eval, row, iht_restarts)
        ranked = sorted(dof_row.items(), key=lambda kv: kv[1])
        runner = next((s for s, _ in ranked if s != best_sup), None)
        sup_names = sorted(LIBR[i] for i in best_sup)
        got = {LIBR[j]: float(coef[i]) for i, j in enumerate(sorted(best_sup))}
        errs = {n: abs(got.get(n, 0.0) - c) / abs(c) for n, c in TRUE_COEF[row].items()}
        supports.append(best_sup)
        coefs.append(coef)
        pareto_data.append({
            "per_k": {str(k): [float(v[0]), sorted(LIBR[i] for i in v[1])]
                      for k, v in per_k.items()},
            "floor_over_N": float(floor / est.N), "knee": int(knee),
            "chosen": sup_names,
        })
        selections.append({
            "N": est.N, "got": got, "true": TRUE_COEF[row],
            "max_err": float(max(errs.values())),
            "arb1_chosen": float(dof_row[best_sup]),
            "arb1_runnerup": float(dof_row[runner]) if runner is not None else None,
            "support_exact": set(got) == set(TRUE_SUP[row]),
        })

        def _perdof5(sup):
            _, s5 = est.corrected_solve(row, sup)
            return float(s5 / (est.N - len(sup)))

        selections[-1]["arb5_chosen"] = _perdof5(best_sup)
        selections[-1]["arb5_runnerup"] = _perdof5(runner) if runner is not None             else None
        print(f"      knee k = {knee} | arbitration winner: {sup_names}")
        print(f"      coefficients: " +
              ", ".join(f"{LIBR[j]}:{coef[i]:+.5f}"
                        for i, j in enumerate(sorted(best_sup))))
        print(f"      support exact: {set(got) == set(TRUE_SUP[row])} | "
              f"max coeff err = {max(errs.values()):.2%} | "
              f"arb ratio (1%) = "
              f"{selections[-1]['arb1_runnerup'] / selections[-1]['arb1_chosen']:.2f}x")

    # symbolic assembly
    expr_rows = [row_expression([LIB_NAMES.index(LIBR[i]) for i in sorted(supports[r])],
                                coefs[r]) for r in range(3)]
    banner("DISCOVERED SYMBOLIC LAWS (terminal logging)")
    for vn, e in zip(VAR_NAMES, expr_rows):
        print(f"\n   d{vn}/dt = {sympy.sstr(e)}")
    print()
    (RES / "discovered_laws_sympy.txt").write_text(
        "\n".join(f"d{vn}/dt = {sympy.sstr(e)}"
                  for vn, e in zip(VAR_NAMES, expr_rows)) + "\n", encoding="utf-8")
    (RES / "discovered_laws.tex").write_text(
        "\n".join(latex_ode(vn, e) for vn, e in zip(VAR_NAMES, expr_rows)) + "\n",
        encoding="utf-8")

    # structural fidelity of the discovered system
    print(f"   {elapsed()} integrating discovered ODE system (structural fidelity) ...")
    t_eval = np.arange(0.0, 60.0, DT)
    sol_true = solve_ivp(lambda t, s: oregonator_rhs(s), (0, 60),
                         np.array([0.6, 0.6, 0.6]), method="Radau",
                         rtol=1e-9, atol=1e-11, t_eval=t_eval, max_step=0.5)
    f_disc, f_disc_vec = lambdified_rows(expr_rows)
    sol_disc = solve_ivp(lambda t, s: f_disc(s), (0, 60),
                         np.array([0.6, 0.6, 0.6]), method="Radau",
                         rtol=1e-9, atol=1e-11, t_eval=t_eval, max_step=0.5)
    per_true = orbit_period(sol_true.t, sol_true.y[0])
    per_disc = orbit_period(sol_disc.t, sol_disc.y[0])
    print(f"      limit-cycle period: true = {per_true:.4f} tau | "
          f"discovered = {per_disc:.4f} tau | rel err "
          f"{abs(per_disc - per_true) / per_true:.2%}")

    # ---- RD / PDE discovery -------------------------------------------------
    print(f"   {elapsed()} PDE discovery on 1-D reaction-diffusion telemetry "
          f"(Laplacian pipeline + same estimator) ...")
    n_grid, L = 192, 0.6
    row_keep_rd = 4 if QUICK else 4

    def rd_features(rd_list, kernels, noise_level=NOISE):
        sigma = [noise_level * f_.reshape(-1, 3).std(axis=0) for f_ in rd_list][0]
        Ys, THs = [], []
        Ys_a = None
        for f_i, Fld in enumerate(rd_list):
            lap_u = rd_laplacian(Fld[..., 0], sigma[0], n_grid, L)
            lap_v = rd_laplacian(Fld[..., 1], sigma[1], n_grid, L)
            laps = np.stack([lap_u, lap_v], axis=-1)
            laps_odd = laps[1::2]
            n_pts = n_grid if not QUICK else 96
            pts = range(0, n_grid, n_grid // n_pts)
            for p in pts:
                extra = laps_odd[:, p, :]
                Yp, THp = kernels.features(Fld[:, p, :], extra=extra)
                Yp, THp = Yp[::STRIDE][::row_keep_rd], THp[::STRIDE][::row_keep_rd]
                Ys.append(Yp)
                THs.append(THp)
        return np.vstack(Ys), np.vstack(THs)

    def rd_noise_gram(rd_list, kernels, n_mc=3, seed0=7000, point_step=8,
                      noise_level=NOISE):
        """MC library-noise Gram for the RD stage.

        The window convolution runs along TIME for each grid point; a subset
        of points is used and the result is rescaled by the point ratio
        (per-point noise statistics are homogeneous).
        """
        sigma = [noise_level * f_.reshape(-1, 3).std(axis=0) for f_ in rd_list][0]
        w = NP_ + 2
        Nl = np.zeros((w, w))
        n_kept = 0
        for i in range(n_mc):
            rm0 = np.random.default_rng(seed0 + 131 * i)
            for Fld in rd_list:
                eps = np.stack([rm0.normal(0.0, sigma[k], Fld.shape[:2])
                                for k in range(3)], axis=-1)
                Fld_e = Fld + eps
                lap_u = rd_laplacian(Fld_e[..., 0], sigma[0], n_grid, L) \
                    - rd_laplacian(Fld[..., 0], sigma[0], n_grid, L)
                lap_v = rd_laplacian(Fld_e[..., 1], sigma[1], n_grid, L) \
                    - rd_laplacian(Fld[..., 1], sigma[1], n_grid, L)
                dth_kin = library_matrix(Fld_e[1::2]) - library_matrix(Fld[1::2])
                dth = np.concatenate([dth_kin, lap_u[1::2][..., None],
                                      lap_v[1::2][..., None]], axis=-1)
                for p in range(0, n_grid, point_step):
                    dTHw = np.column_stack([
                        kernels.conv_valid(dth[:, p, j], kernels.ker_t)
                        for j in range(dth.shape[2])])
                    dTHw = dTHw[::STRIDE][::row_keep_rd]
                    Nl += dTHw.T @ dTHw
                    n_kept += len(dTHw)
        return Nl * (n_grid / ((n_grid + point_step - 1) // point_step)) / n_mc

    Y_rd, TH_rd = rd_features(rd_noisy, KER_RD)
    GN_rd = rd_noise_gram(rd_noisy, KER_RD)
    RD_NAMES_ALL = LIB_NAMES + ["lap_u", "lap_v"]
    pre_rd = CalibratedSTRidge(Y_rd, TH_rd, GN_rd, RD_NAMES_ALL)
    kept_rd, _ = prune_collinear(pre_rd.G, RD_NAMES_ALL)
    TH_rd = TH_rd[:, kept_rd]
    GN_rd = GN_rd[np.ix_(kept_rd, kept_rd)]
    RD_NAMES = [RD_NAMES_ALL[i] for i in kept_rd]
    est_rd = CalibratedSTRidge(Y_rd, TH_rd, GN_rd, RD_NAMES)
    Y_rd_a, TH_rd_a = rd_features(rd_arb, KER_RD, noise_level=NOISE_ARB)
    GN_rd_a = rd_noise_gram(rd_arb, KER_RD, noise_level=NOISE_ARB, n_mc=2,
                            seed0=7500)
    est_rd_arb = CalibratedSTRidge(Y_rd_a, TH_rd_a[:, kept_rd],
                                   GN_rd_a[np.ix_(kept_rd, kept_rd)], RD_NAMES)
    Y_rd_b, TH_rd_b = rd_features(rd_arb2, KER_RD, noise_level=NOISE_ARB)
    GN_rd_b = rd_noise_gram(rd_arb2, KER_RD, noise_level=NOISE_ARB, n_mc=2,
                            seed0=7600)
    est_rd_arb_eval = CalibratedSTRidge(Y_rd_b, TH_rd_b[:, kept_rd],
                                        GN_rd_b[np.ix_(kept_rd, kept_rd)], RD_NAMES)
    print(f"   {elapsed()} RD windows: {len(Y_rd)}")
    # compositional discovery: the reaction kinetics are LOCAL, so the
    # symbolic law already discovered on the well-mixed telemetry is carried
    # over as the kinetic core; the spatial telemetry decides which transport
    # (Laplacian) terms the PDE needs and quantifies them.
    rd_supports, rd_coefs = [], []
    lap_idx = {"lap_u": RD_NAMES.index("lap_u"), "lap_v": RD_NAMES.index("lap_v")}
    D_hat = {}
    for row, vn in enumerate(VAR_NAMES):
        core = [RD_NAMES.index(LIBR[i]) for i in supports[row]
                if LIBR[i] in RD_NAMES]
        # structural significance for the transport terms is judged on the
        # 1 % arbitration instrument (4x the operative precision), while the
        # reported transport coefficient carries both instruments' estimates
        sup_l = sorted(core + list(lap_idx.values()))
        coef_a, sse_a = est_rd_arb.corrected_solve(row, frozenset(sup_l))
        n_eff_rd = max(est_rd_arb.N * est_rd_arb.eff_ratio, 2.0 * len(sup_l))
        A_c = est_rd_arb.G[np.ix_(sup_l, sup_l)] - est_rd_arb.GNn[np.ix_(sup_l, sup_l)]
        ev, EV = np.linalg.eigh(A_c)
        ev = np.clip(ev, 1e-6 * ev.max(), None)
        cov = EV @ np.diag(1.0 / ev) @ EV.T
        sigma2 = sse_a / max(n_eff_rd - len(sup_l), 1)
        se = np.sqrt(np.maximum(np.diag(cov) * sigma2, 0.0))
        keep = list(core)
        d_row = {}
        for name, j in lap_idx.items():
            if name.split("_")[1] != vn:
                continue
            pos = sup_l.index(j)
            if abs(coef_a[pos]) > 2.0 * se[pos]:
                keep.append(j)
                d_row[name] = (float(coef_a[pos]), float(se[pos]))
        keep = sorted(set(keep))
        coef, sse = est_rd.corrected_solve(row, frozenset(keep))
        rd_supports.append(frozenset(keep))
        rd_coefs.append(coef)
        for name, (cval, serr) in d_row.items():
            D_hat[name] = (cval, serr)
    du = D_hat.get("lap_u", (0.0, 0.0))
    dv = D_hat.get("lap_v", (0.0, 0.0))

    def _fmt(name, pair, d_true):
        cval, serr = pair
        if cval == 0.0 and serr == 0.0:
            return f"{name}: < 2 sigma (not resolved)"
        return f"{name}: {cval:.3e} +- {serr:.1e} (true {d_true:.3e})"

    for row, vn in enumerate(VAR_NAMES):
        sup_names = sorted(RD_NAMES[i] for i in rd_supports[row])
        print(f"      RD d{vn}/dt: {sup_names}")
        print("      coefficients: " +
              ", ".join(f"{RD_NAMES[j]}:{rd_coefs[row][i]:+.5f}"
                        for i, j in enumerate(sorted(rd_supports[row]))))
    print(f"      diffusion: {_fmt('D_u', du, D_U)} | {_fmt('D_v', dv, D_V)}")
    results["12B_rd"] = {
        "D_u": {"true": D_U, "discovered": float(du[0]), "se": float(du[1])},
        "D_v": {"true": D_V, "discovered": float(dv[0]), "se": float(dv[1])},
        "supports": {vn: sorted(RD_NAMES[i] for i in s)
                     for vn, s in zip(VAR_NAMES, rd_supports)},
    }

    # ---- Onsager gate -------------------------------------------------------
    print(f"   {elapsed()} Onsager gate: near-steady-state refinement at "
          f"f = {F_STO} (excitable regime) ...")
    sn_sto = []
    for i in range(n_sto):
        s0 = np.array([rng.uniform(0.2, 0.9), rng.uniform(0.5, 2.5),
                       rng.uniform(0.2, 0.9)])
        S = simulate_0d(s0, 40.0, f=F_STO)
        sn_sto.append(add_noise(S, NOISE, np.random.default_rng(300 + i)))
    sn_sto_arb = [add_noise(S, NOISE_ARB, np.random.default_rng(400 + i))
                  for i, S in enumerate(sn_sto)]
    Y_s, TH_s = build_weak_dataset(sn_sto, KER0, seed=9)
    GN_s = noise_gram_mc(sn_sto, KER0, n_mc=3, seed0=9500)
    pre_s = CalibratedSTRidge(Y_s, TH_s, GN_s, LIB_NAMES)
    kept_s, _ = prune_collinear(pre_s.G, LIB_NAMES)
    TH_s = TH_s[:, kept_s]
    GN_s = GN_s[np.ix_(kept_s, kept_s)]
    est_s = CalibratedSTRidge(Y_s, TH_s, GN_s, [LIB_NAMES[i] for i in kept_s])
    Y_sa, TH_sa = build_weak_dataset(sn_sto_arb, KER0, seed=10)
    GN_sa = noise_gram_mc(sn_sto_arb, KER0, noise_level=NOISE_ARB, seed0=8600,
                          n_mc=3)
    est_sa = CalibratedSTRidge(Y_sa, TH_sa[:, kept_s],
                               GN_sa[np.ix_(kept_s, kept_s)],
                               [LIB_NAMES[i] for i in kept_s])
    sn_sto_arb2 = [add_noise(S, NOISE_ARB, np.random.default_rng(700 + i))
                   for i, S in enumerate(sn_sto)]
    Y_sb, TH_sb = build_weak_dataset(sn_sto_arb2, KER0, seed=11)
    GN_sb = noise_gram_mc(sn_sto_arb2, KER0, noise_level=NOISE_ARB, seed0=8700,
                          n_mc=3)
    est_sa_eval = CalibratedSTRidge(Y_sb, TH_sb[:, kept_s],
                                    GN_sb[np.ix_(kept_s, kept_s)],
                                    [LIB_NAMES[i] for i in kept_s])
    # structure transfer: the discovered kinetic law is refit (coefficients
    # only) in the new stoichiometric regime; f enters the w-coefficient
    sto_supports, sto_coefs = [], []
    names_s = est_s.names
    for row in range(3):
        sup_s = frozenset(names_s.index(LIBR[i]) for i in supports[row]
                          if LIBR[i] in names_s)
        coef, sse = est_s.corrected_solve(row, sup_s)
        sto_supports.append(sup_s)
        sto_coefs.append(coef)
    expr_sto = [row_expression([LIB_NAMES.index(names_s[i]) for i in sorted(sto_supports[r])],
                               sto_coefs[r]) for r in range(3)]
    try:
        ss_sto = steady_state(expr_sto, [[0.5, 1.5, 0.5], [0.8, 1.0, 0.7],
                                         [0.3, 2.0, 0.4]])
        ss_osc = steady_state(expr_rows, [[0.5, 1.5, 0.5], [0.8, 1.0, 0.7],
                                          [0.2, 3.0, 0.3]])
        J_sto = jacobian_at(expr_sto, ss_sto if ss_sto is not None else [0.5, 1.5, 0.5])
        J_osc = jacobian_at(expr_rows, ss_osc if ss_osc is not None else [0.5, 1.5, 0.5])
        asym_sto, eig_sto, cost_sto = onsager_analysis(J_sto)
        asym_osc, eig_osc, cost_osc = onsager_analysis(J_osc)
    except Exception as exc:                               # noqa: BLE001
        print(f"      Onsager gate degraded ({exc}); reporting NaN row")
        asym_sto = asym_osc = float("nan")
        eig_sto = eig_osc = cost_sto = cost_osc = float("nan")
    tr_sto = float(np.trace(J_sto))
    tr_osc = float(np.trace(J_osc))
    print(f"      near-steady-state (f = {F_STO}): Onsager asymmetry "
          f"{asym_sto:.3f}, min eig(L_sym) = {eig_sto:.4g}, tr(J) = {tr_sto:.3f}")
    print(f"      oscillatory regime (f = {F_OSC}):   Onsager asymmetry "
          f"{asym_osc:.3f}, min eig(L_sym) = {eig_osc:.4g}, tr(J) = {tr_osc:.3f}")
    results["12B"] = {
        "score": {"support_exact": {vn: bool(selections[r]["support_exact"])
                                    for r, vn in enumerate(VAR_NAMES)},
                  "max_coefficient_error": {vn: selections[r]["max_err"]
                                            for r, vn in enumerate(VAR_NAMES)}},
        "period_true": float(per_true), "period_discovered": float(per_disc),
        "supports": {vn: sorted(LIB_NAMES[i] for i in s)
                     for vn, s in zip(VAR_NAMES, supports)},
        "coefficients": {vn: {LIB_NAMES[j]: float(coefs[r][i])
                              for i, j in enumerate(sorted(supports[r]))}
                         for r, vn in enumerate(VAR_NAMES)},
        "pareto": pareto_data,
        "onsager": {
            "near_steady_state_f2.8": {"asymmetry": asym_sto, "min_eig_sym": eig_sto,
                                       "trace_J": tr_sto},
            "oscillatory_f1.1": {"asymmetry": asym_osc, "min_eig_sym": eig_osc,
                                 "trace_J": tr_osc},
        },
    }

    # ================================================================== 12C
    banner("MODULE 12C — autonomous Hamiltonian & Lyapunov functional extraction")
    long_traj = sol_true.y.T
    attr_samples = long_traj[::10]
    off_samples = np.vstack([t[::40] for t in kicked_trajs + transient_trajs])
    print(f"   {elapsed()} training Hamiltonian core H(x): "
          f"gradH . f = 0 on the attractor band ...")
    gap = 25
    pairs = (long_traj[:-gap:7], long_traj[gap::7])
    netH = train_scalar_functional(f_disc_vec, "hamiltonian", attr_samples,
                                   attr_samples, n_steps=net_steps, seed=1,
                                   pairs=pairs)
    print(f"   {elapsed()} training Lyapunov functional V(x): "
          f"V >= 0, dV/dt <= 0 off-attractor ...")
    netV = train_scalar_functional(f_disc_vec, "lyapunov", off_samples,
                                   attr_samples, n_steps=net_steps, seed=2)

    # H conservation audit
    H_at, gH_at = nn_value_and_grad(netH, attr_samples)
    F_at = f_disc_vec(attr_samples)
    rel_res = float(np.mean(np.abs((gH_at * F_at).sum(-1)) /
                            (np.linalg.norm(gH_at, axis=1) *
                             np.linalg.norm(F_at, axis=1) + 1e-12)))
    Ht, _ = nn_value_and_grad(netH, long_traj[::20])
    h_range = float(Ht.max() - Ht.min() + 1e-12)
    drift = float(np.max(np.abs(Ht - Ht[0])) / h_range)
    print(f"      H audit: relative |dH/dt| residual on attractor = {rel_res:.2e} | "
          f"max drift over 60 tau / range(H) = {drift:.2e}")

    # V certificate audit
    V_off, gV_off = nn_value_and_grad(netV, off_samples)
    F_off = f_disc_vec(off_samples)
    dV_off = (gV_off * F_off).sum(-1)
    cert = float(np.mean(dV_off <= 1e-3))
    print(f"      V audit: certificate rate dV/dt <= 0 off-attractor = "
          f"{cert:.2%} | mean dV/dt = {dV_off.mean():.3e} | "
          f"V >= 0 everywhere: {bool(V_off.min() >= 0)}")
    V_cycle, gV_cycle = nn_value_and_grad(netV, attr_samples)
    dV_cycle = (gV_cycle * F_at).sum(-1)
    print(f"      NESS balance on attractor: <dV/dt> = {dV_cycle.mean():.3e} "
          f"(~0) | rms = {dV_cycle.std():.3e} (> 0: circulating flux)")

    spec_true = lle_benettin(oregonator_rhs, np.array([0.6, 0.6, 0.6]), T=lle_T, dt=1.0e-3)
    spec_disc = lle_benettin(f_disc, np.array([0.6, 0.6, 0.6]), T=lle_T, dt=1.0e-3)
    lle_true, lle_disc = float(spec_true[0]), float(spec_disc[0])
    print(f"      Lyapunov spectrum (1/tau): true = {np.round(spec_true, 4)} | "
          f"discovered = {np.round(spec_disc, 4)}  (lambda_1 ~ 0: limit cycle; "
          f"sum < 0: dissipative contraction)")
    results["12C"] = {
        "H_relative_residual": rel_res, "H_max_drift_60tau": drift,
        "V_certificate_rate": cert, "V_min": float(V_off.min()),
        "V_mean_off": float(dV_off.mean()),
        "NESS_mean_dV_on_cycle": float(dV_cycle.mean()),
        "NESS_rms_dV_on_cycle": float(dV_cycle.std()),
        "lle_true": float(lle_true), "lle_discovered": float(lle_disc),
        "lyapunov_spectrum_true": [float(x) for x in spec_true],
        "lyapunov_spectrum_discovered": [float(x) for x in spec_disc],
    }

    # ---- exports ------------------------------------------------------------
    banner("LAW EXPORT — LaTeX / SymPy / C++ execution kernels")
    export_cpp_kernel(RES / "phase12_kernels.hpp", expr_rows, netH, netV)
    try:
        cpp_ok = cpp_demo_and_verify(RES / "phase12_kernels.hpp", expr_rows, netH, netV)
    except Exception as exc:                               # noqa: BLE001
        print(f"   {elapsed()} C++ execution check skipped ({exc})")
        cpp_ok = None
    if cpp_ok:
        results["12B"]["cpp_kernel_max_abs_diff"] = cpp_ok

    # ---- figures ------------------------------------------------------------
    import torch
    torch.save({"netH": netH.state_dict(), "netV": netV.state_dict(),
                "scaleH": [float(x) for x in netH.scale.detach()],
                "scaleV": [float(x) for x in netV.scale.detach()]},
               RES / "_nets.pt")
    np.savez_compressed(RES / "_figdata.npz",
                        sol_true=sol_true.y.T, sol_disc=sol_disc.y.T,
                        rd_field=rd_noisy[0][..., 0], long_traj=long_traj,
                        t_grid=sol_true.t,
                        **{f"kick{i}": k for i, k in enumerate(kicked_trajs)})
    with open(RES / "_figmeta.json", "w") as fh:
        json.dump({"per_true": float(per_true), "per_disc": float(per_disc),
                   "pareto": pareto_data,
                   "selections": [{k: (float(v) if isinstance(v, (int, float)) else v)
                                   for k, v in sel.items()} for sel in selections],
                   "true_k": {vn: TRUE_K[r] for r, vn in enumerate(VAR_NAMES)}},
                  fh, default=str)
    print(f"   {elapsed()} rendering publication figures (300 DPI) ...")
    fig1(sol_true, sol_disc, per_true, per_disc, rd_noisy[0][..., 0])
    fig2(pareto_data, selections,
         {vn: TRUE_K[r] for r, vn in enumerate(VAR_NAMES)})
    fig3(netV, netH, f_disc, f_disc_vec, long_traj, kicked_trajs)

    (RES / "phase12_results.json").write_text(
        json.dumps(results, indent=2, default=float), encoding="utf-8")
    banner("PHASE 12 COMPLETE")
    print(f"   results  -> results_phase12/  (JSON, LaTeX, SymPy, C++ kernel)")
    print(f"   figures  -> figures_phase12/  (fig1, fig2, fig3 @ 300 DPI)")
    print(f"   total wall time: {time.time() - T0:.1f} s")


if __name__ == "__main__":
    main()
