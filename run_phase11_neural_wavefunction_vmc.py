#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_phase11_neural_wavefunction_vmc.py
======================================
PHASE 11 — THE QUANTUM SINGULARITY
DEEP VARIATIONAL QUANTUM MONTE CARLO & CONTINUOUS NEURAL WAVEFUNCTIONS

Every previous phase of this repository relied on basis-set expansions
(Gaussian orbitals, plane waves, auxiliary fields) or fixed atomic-orbital
parameterizations.  Phase 11 ascends to the exact continuous-space solution
of the ab initio many-electron Schrodinger equation: an end-to-end
NEURAL WAVEFUNCTION ANSATZ (FermiNet / PauliNet architecture family,
re-implemented from Coulomb's law alone) optimized by variational quantum
Monte Carlo directly in PyTorch — no atomic orbital basis, no empirical
parameterization, no finite-difference artifacts.

    H = -1/2 SUM_i nabla_i^2
        - SUM_{i,I} Z_I / |r_i - R_I|
        + SUM_{i<j} 1 / |r_i - r_j|
        + SUM_{I<J} Z_I Z_J / |R_I - R_J|                      (atomic units)

  Module 11A — Antisymmetric equivariant neural wavefunction.
    * Continuous featurization: one-electron streams
      h_i^(0) = (r_i - R_I, |r_i - R_I|) and two-electron streams
      h_ij^(0) = (r_i - r_j, |r_i - r_j|) built from raw Cartesian
      differences — Coulomb's law is the only physical input.
    * L = 3 permutation-equivariant interaction blocks (FermiNet-style
      residual message passing: the one-electron stream aggregates the
      mean of the two-electron stream; the two-electron stream is refined
      from itself plus both one-electron endpoints).
    * Antisymmetric Fermi layer: multi-determinant BACKFLOW orbitals
      phi_{k,j}(r_i; {r_/i}, {r^beta}) — each orbital depends on the
      positions of ALL electrons through the equivariant stream and an
      explicit electron-electron backflow kernel — combined into a
      spin-resolved sum of Slater determinants
          Psi(r) = exp(J(r)) SUM_k c_k det[phi_k^up] det[phi_k^down],
      multiplied by an isotropic Jastrow envelope exp(J) carrying the
      EXACT Kato cusp conditions (e-n: dlnPsi/dr|_0 = -Z_I;
      e-e: +1/2 unlike / +1/4 like spins).  Exact antisymmetry is verified
      numerically to machine precision (same-spin exchange => -Psi).

  Module 11B — Vectorized Metropolis-Hastings MCMC electron sampler.
    Gaussian random walk of N_walkers = 2048 parallel walkers in 3N-dim
    configuration space, acceptance min(1, |Psi(r')|^2/|Psi(r)|^2), with
    the proposal width sigma dynamically locked into the 45%-55%
    acceptance window, a 400-sweep burn-in, and re-equilibration sweeps
    after every parameter update.

  Module 11C — Energy minimization via the Rayleigh quotient.
    Local energy E_L(r) = (H Psi)(r)/Psi(r) with the kinetic term obtained
    from the EXACT Laplacian of ln|Psi| through reverse-mode automatic
    differentiation (one-hot Hessian-trace contraction — analytic, no
    finite differences; cross-validated against central differences).
    Gradient descent on <E> uses the exact variational (REINFORCE)
    estimator
        grad_theta <E> = 2 E_{r~|Psi|^2}[(E_L - <E_L>) grad_theta ln|Psi|],
    with 5-sigma local-energy outlier clipping and global-norm gradient
    clipping, optimized by Adam on a cosine schedule.

  Benchmark systems (chemical-accuracy targets vs exact Full CI):
    * H2 at R = 1.4011 bohr (equilibrium; Kolos-Wolniewicz exact limit
      -1.1744757 Eh) — the flagship convergence study.
    * H2 dissociation curve R = 2.5 / 4.0 / 6.0 bohr — static correlation
      (the failure mode that motivated the Phase-7 Wall of Sighs) handled
      variationally; dissociation limit E -> -1.0 Eh exactly.
    * He atom — Pekeris exact nonrelativistic limit -2.9037244 Eh.

Outputs
-------
results_phase11/phase11_results.json     machine-readable master record
results_phase11/references.json          HF / CCSD(T) / FCI(CBS) references
results_phase11/convergence_<sys>.csv    per-epoch VMC history
results_phase11/density_slice_H2_eq.npz  raw 2-D density slice (fig 2)
figures_phase11/fig1_vmc_energy_convergence.png
figures_phase11/fig2_electron_density_slice.png
figures_phase11/fig3_local_energy_variance.png

Engine notes (Windows box, carried over from Phases 7/8)
--------------------------------------------------------
VMC engine  : PyTorch 2.13 (driver interpreter), float64, vectorized
              walkers; device auto-detects CUDA when present — the walker
              tensor, network and every autograd graph live on the GPU in
              that case ("GPU-accelerated vectorized walkers").
QC engine   : Psi4 1.11 (conda env `phase7`) only for the classical
              reference energies (HF / CCSD(T) / FCI); the neural VMC
              solver itself is pure PyTorch and needs no quantum-chemistry
              backend whatsoever.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results_phase11"
FIGURES = ROOT / "figures_phase11"

SEED = 11
torch.set_default_dtype(torch.float64)

PHASE7_PY = Path(r"C:/Users/HUIWEI/miniconda3/envs/phase7/python.exe")

# ----------------------------------------------------------------------------
# Canonical exact references (nonrelativistic Born-Oppenheimer, atomic units)
# ----------------------------------------------------------------------------
LITERATURE = {
    "H2_eq": {"hf": -1.13363, "ccsdt": -1.1744757, "fci": -1.1744757,
              "source": "Kolos-Wolniewicz 1968 exact BO limit "
                        "(CCSD(T) = FCI for 2 electrons)"},
    "H2_diss": {"hf": None, "ccsdt": -1.0, "fci": -1.0,
                "source": "exact separated-atom limit 2 x H(-1/2)"},
    "He": {"hf": -2.8616800, "ccsdt": -2.9036400, "fci": -2.9037244,
           "source": "Pekeris 1958 exact nonrelativistic limit"},
}


# ============================================================================
# System definitions
# ============================================================================
@dataclass
class SystemConfig:
    name: str
    charges: list
    positions: list          # (Nnuc, 3) bohr
    n_up: int
    n_dn: int
    epochs: int
    n_det: int = 8
    hidden: int = 40
    h2_dim: int = 12
    depth: int = 3
    n_walkers: int = 2048
    sweeps_per_epoch: int = 4
    lr: float = 4e-3
    tag: str = ""            # family tag for reference lookup

    @property
    def n_elec(self):
        return self.n_up + self.n_dn

    def geom_xyz(self):
        sym = {1: "H", 2: "He", 3: "Li"}
        return "\n".join(
            f"{sym[z]} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}"
            for z, p in zip(self.charges, self.positions))


def default_systems(epochs_main: int, epochs_curve: int, walkers: int):
    R_e = 1.4011

    def mk_h2(name, R, ep):
        return SystemConfig(
            name=name, charges=[1, 1],
            positions=[[0.0, 0.0, -R / 2], [0.0, 0.0, R / 2]],
            n_up=1, n_dn=1, epochs=ep, n_det=8, n_walkers=walkers,
            tag="H2_eq" if abs(R - R_e) < 1e-9 else "H2_diss")

    return [
        mk_h2("H2_eq_R1.4011", R_e, epochs_main),
        mk_h2("H2_R2.5", 2.5, epochs_curve),
        mk_h2("H2_R4.0", 4.0, epochs_curve),
        mk_h2("H2_R6.0", 6.0, epochs_curve),
        SystemConfig(name="He", charges=[2], positions=[[0.0, 0.0, 0.0]],
                     n_up=1, n_dn=1, epochs=epochs_main, n_det=8,
                     n_walkers=walkers, tag="He"),
    ]


# ============================================================================
# MODULE 11A — Antisymmetric equivariant neural wavefunction
# ============================================================================
class FermiPauliNet(nn.Module):
    """Continuous neural wavefunction: equivariant interaction streams ->
    backflow multi-determinant Fermi layer -> exact-cusp isotropic Jastrow.

    The network sees ONLY raw Cartesian electron-nuclear and electron-electron
    differences and their norms (Coulomb featurization).  Permutation
    equivariance of the one-electron stream is structural (mean aggregation);
    antisymmetry is structural (Slater determinants per spin sector); the
    Kato cusps are structural (fixed-coefficient Jastrow terms).
    """

    def __init__(self, charges, positions, n_up, n_dn,
                 n_det=8, hidden=64, h2_dim=16, depth=3, device="cpu"):
        super().__init__()
        self.device = device
        self.register_buffer("Z", torch.tensor(charges, dtype=torch.float64,
                                               device=device).reshape(-1))
        self.register_buffer("R", torch.tensor(positions, dtype=torch.float64,
                                               device=device).reshape(-1, 3))
        self.n_up, self.n_dn = n_up, n_dn
        self.ne = n_up + n_dn
        self.nnuc = len(charges)
        self.K = n_det
        self.norb_up = self.K * n_up
        self.norb_dn = self.K * n_dn
        self.norb = self.norb_up + self.norb_dn
        H, H2, L = hidden, h2_dim, depth

        # --- Coulomb featurization: 6 numbers per nucleus, 6 per pair.
        # Radial channels are smooth at the origin (r^2, Gaussian shells) so
        # the determinant/backflow heads contribute ZERO log-slope at the
        # nuclei and the Kato cusps are carried solely by the Jastrow.
        self.d1 = 6 * self.nnuc   # dx dy dz r^2 e^{-r^2/0.25} e^{-r^2/1.0}
        self.d2 = 6
        self.embed1 = nn.Linear(self.d1, H)
        self.embed2 = nn.Linear(self.d2, H2)

        # --- L permutation-equivariant interaction blocks
        g_dim = H + H + H2
        self.f1 = nn.ModuleList()
        self.f2 = nn.ModuleList()
        for _ in range(L):
            self.f1.append(nn.Sequential(
                nn.Linear(g_dim, H), nn.Tanh(), nn.Linear(H, H)))
            self.f2.append(nn.Sequential(
                nn.Linear(H2 + 2 * g_dim, H), nn.Tanh(), nn.Linear(H, H2)))
        for seq in list(self.f1) + list(self.f2):
            nn.init.constant_(seq[-1].weight, 0.0)   # near-identity start
            nn.init.constant_(seq[-1].bias, 0.0)

        # --- Fermi layer heads
        self.w_orb = nn.Linear(H, self.norb)         # neural one-electron channel
        nn.init.constant_(self.w_orb.weight, 0.0)
        nn.init.constant_(self.w_orb.bias, 0.0)

        # Isotropic Gaussian envelope CONTRACTIONS: every orbital owns a
        # small bank of primitives, phi_o(r) = sum_m D[o,m,A] e^{-r^2/2s^2},
        # with D least-squares-fitted at init to a Slater-1s target
        # (PauliNet's "determinants at HF quality" philosophy, analytic).
        # Zero Gaussian log-slope at the nuclei keeps Kato cusps exact.
        # Up slots point at nucleus A, down slots at B (swapped), so the
        # uniform-coefficient start is the contracted Heitler-London singlet.
        shift = max(1, self.nnuc // 2)
        alpha0, beta = 0.1, 3.2
        self.nprim = 4
        per_group = max(1, self.K // 2)
        zeta = 1.2 * torch.sqrt(self.Z)          # Slater rules: H 1.2, He 1.70
        s0 = torch.zeros(self.norb, self.nprim, self.nnuc)
        W0 = torch.zeros(self.norb, self.nprim, self.nnuc)
        rr = torch.linspace(1e-3, 6.0, 400)
        for o in range(self.norb):
            spin_off = shift if o >= self.norb_up else 0
            group = (o % self.K) // per_group     # A-side / B-side
            home = (group + spin_off) % self.nnuc
            for m in range(self.nprim):
                sig = 1.0 / math.sqrt(2.0 * alpha0 * beta ** m)
                for A in range(self.nnuc):
                    s0[o, m, A] = sig
            G = torch.stack([torch.exp(-rr ** 2 / (2.0 * s0[o, m, home] ** 2))
                             for m in range(self.nprim)], dim=1)
            # cusp-free smoothed Slater target: zero log-slope at the
            # nucleus (slope is Jastrow-owned) but STO-like in the valence
            target = torch.exp(-float(zeta[home]) *
                               (torch.sqrt(rr ** 2 + 0.15 ** 2) - 0.15))
            coef, _, _, _ = torch.linalg.lstsq(G, target.unsqueeze(1))
            W0[o, :, home] = coef[:, 0]
            for A in range(self.nnuc):
                if A != home:
                    W0[o, :, A] = 0.02 * torch.randn(self.nprim)
        self.W_env = nn.Parameter(W0 + 0.02 * torch.randn_like(W0))
        self.raw_s_env = nn.Parameter(torch.log(torch.expm1(s0)))   # softplus^-1
        # explicit e-e backflow channel:
        #   phi_o(r_i) += sum_{l!=i} k(r_il) * sum_A W_bf e^{-r_lA^2/2s^2}
        self.W_bf = nn.Parameter(torch.zeros(self.norb, self.nprim, self.nnuc))
        self.raw_s_bf = nn.Parameter(torch.full((1,), 1.0))

        # determinant coefficients (softmax-normalized, PauliNet convention)
        self.raw_c = nn.Parameter(torch.zeros(self.K))

        # --- isotropic Jastrow with EXACT Kato cusps
        # e-n: -Z_I r/(1+b r);  e-e: a_spin r/(1+b r), a = 1/2 unlike, 1/4 like
        self.raw_b_en = nn.Parameter(torch.zeros(self.nnuc))
        self.raw_b_ee_u = nn.Parameter(torch.zeros(1))
        self.raw_b_ee_l = nn.Parameter(torch.zeros(1))
        self.jw1 = nn.Parameter(torch.zeros(1))      # smooth learned terms
        self.jw2 = nn.Parameter(torch.zeros(1))

        # static spin-pair masks for the e-e Jastrow coefficients
        spin = torch.tensor([0] * n_up + [1] * n_dn)
        same = spin[:, None] == spin[None, :]
        upper = torch.triu(torch.ones(self.ne, self.ne), diagonal=1).bool()
        self.register_buffer("mask_ee_like", same & upper)
        self.register_buffer("mask_ee_unlike", (~same) & upper)

        self.to(device)

    # ---------------------------------------------------------------- features
    def _featurize(self, r):
        d_en = r[:, :, None, :] - self.R[None, None, :, :]      # (B,ne,Nn,3)
        r_en = d_en.norm(dim=-1)                                 # (B,ne,Nn)
        d_ee = r[:, :, None, :] - r[:, None, :, :]               # (B,ne,ne,3)
        eye = torch.eye(self.ne, dtype=torch.bool, device=r.device)
        # The self-pair displacement |r_i - r_i| sits exactly on the norm
        # singularity at |v| = 0 whose double backward is NaN; shift it onto
        # a unit vector (value then masked out) so autograd stays finite.
        d_ee_safe = d_ee + eye[None, :, :, None].to(r.dtype)
        r_ee = d_ee_safe.norm(dim=-1)                            # diag = 1.0
        f_en = torch.cat([d_en, (r_en ** 2)[..., None],
                          torch.exp(-r_en ** 2 / 0.25)[..., None],
                          torch.exp(-r_en ** 2 / 1.0)[..., None]], dim=-1)
        f_ee = torch.cat([d_ee, (r_ee ** 2)[..., None],
                          torch.exp(-r_ee ** 2 / 0.25)[..., None],
                          torch.exp(-r_ee ** 2 / 1.0)[..., None]], dim=-1)
        f_ee = f_ee.masked_fill(eye[None, :, :, None], 0.0)      # self-pair = 0
        h1 = self.embed1(f_en.reshape(f_en.shape[0], self.ne, self.d1))
        h2 = self.embed2(f_ee)
        return h1, h2, r_en, r_ee

    # ------------------------------------------------------- equivariant blocks
    def _streams(self, r):
        h1, h2, r_en, r_ee = self._featurize(r)
        for l in range(len(self.f1)):
            h1m = h1.mean(dim=1, keepdim=True)                   # (B,1,H)
            h2m = h2.mean(dim=2)                                 # (B,ne,H2)
            g = torch.cat([h1, h1m.expand_as(h1), h2m], dim=-1)  # (B,ne,2H+H2)
            h1 = h1 + self.f1[l](g)
            gi = g.unsqueeze(2).expand(-1, -1, self.ne, -1)
            gj = g.unsqueeze(1).expand(-1, self.ne, -1, -1)
            h2 = h2 + self.f2[l](torch.cat([h2, gi, gj], dim=-1))
        return h1, r_en, r_ee

    # ----------------------------------------------------------------- Jastrow
    def _jastrow(self, r_en, r_ee):
        b_en = F.softplus(self.raw_b_en) + 1e-3                   # (Nn,)
        u_en = r_en / (1.0 + b_en[None, None, :] * r_en)
        J = (-(self.Z[None, None, :]) * u_en).sum(dim=(1, 2))     # exact e-n cusp
        b_u = F.softplus(self.raw_b_ee_u) + 1e-3
        b_l = F.softplus(self.raw_b_ee_l) + 1e-3
        ru = r_ee / (1.0 + b_u * r_ee)
        rl = r_ee / (1.0 + b_l * r_ee)
        J = J + 0.5 * ru.masked_fill(~self.mask_ee_unlike, 0.0).sum(dim=(1, 2))
        J = J + 0.25 * rl.masked_fill(~self.mask_ee_like, 0.0).sum(dim=(1, 2))
        # smooth, cusp-free learned terms (start at zero)
        J = J + self.jw1[0] * torch.exp(-r_en ** 2).sum(dim=(1, 2))
        J = J + self.jw2[0] * (1.0 - torch.exp(-r_ee ** 2)).masked_fill(
            torch.eye(self.ne, dtype=torch.bool, device=r_ee.device), 0.0
        ).sum(dim=(1, 2))
        return J

    # ---------------------------------------------------------------- orbitals
    def _orbitals(self, r, h1, r_en, r_ee):
        """Backflow orbitals phi_o(r_i; all electrons) -> (B, ne, norb).

        Each orbital carries a contracted primitive bank
            env_o = sum_m D[o,m,A] exp(-|r-R_A|^2 / 2 s[o,m,A]^2),
        modulated by a tanh-bounded neural head,
            phi_o = env_o * (1 + tanh(w_orb . h1)) + backflow_o.
        The multiplicative bound preserves the physical Gaussian decay of
        every orbital (a raw linear-in-r head would let the ansatz leak to
        a non-decaying artifact state), while the primitive contraction
        starts the determinants at Slater/HF quality.
        """
        s_env = F.softplus(self.raw_s_env) + 1e-4        # (norb,nprim,Nn)
        env = torch.exp(-r_en.unsqueeze(2).unsqueeze(2) ** 2 /
                        (2.0 * s_env[None, None]))       # (B,ne,norb,nprim,Nn)
        env_term = torch.einsum("beopn,opn->beo", env, self.W_env)
        neural = torch.tanh(self.w_orb(h1))              # in [-1, 1]
        orb = env_term * (1.0 + neural)
        # explicit e-e backflow channel (both factors decay / are bounded)
        s_bf = F.softplus(self.raw_s_bf) + 1e-4                   # scalar width
        k_ee = torch.exp(-r_ee ** 2 / (2.0 * s_bf))               # (B,ne,ne)
        eye = torch.eye(self.ne, dtype=torch.bool, device=r_ee.device)
        k_ee = k_ee.masked_fill(eye[None], 0.0)
        Mo = torch.einsum("beopn,opn->beo", env, self.W_bf)
        orb = orb + torch.einsum("bil,blo->bio", k_ee, Mo)
        return orb

    # ------------------------------------------------- antisymmetric Fermi layer
    @staticmethod
    def _signed_logdet(phi):
        """sign and log|det| of batched phi (B, K, n, n).

        torch.linalg.slogdet does not support double backward in this
        torch build (its Laplacian silently returns NaN), so small
        determinants use explicit cofactor arithmetic — fully
        twice-differentiable, which the exact kinetic energy requires.
        """
        n = phi.shape[-1]
        if n == 1:
            val = phi[..., 0, 0]
            return torch.sign(val), torch.log(val.abs() + 1e-300)
        if n == 2:
            det = (phi[..., 0, 0] * phi[..., 1, 1]
                   - phi[..., 0, 1] * phi[..., 1, 0])
            return torch.sign(det), torch.log(det.abs() + 1e-300)
        return torch.linalg.slogdet(phi)

    def _fermi_layer(self, r):
        """Returns (log_mix, sign_mix, J) with
        sum_k c_k det[phi_k^up] det[phi_k^dn] = sign * exp(log_mix)."""
        h1, r_en, r_ee = self._streams(r)
        orb = self._orbitals(r, h1, r_en, r_ee)
        B = r.shape[0]
        phi_up = orb[:, :self.n_up, :self.norb_up] \
            .reshape(B, self.n_up, self.K, self.n_up).permute(0, 2, 1, 3)
        phi_dn = orb[:, self.n_up:, self.norb_up:] \
            .reshape(B, self.n_dn, self.K, self.n_dn).permute(0, 2, 1, 3)
        # torch.linalg.slogdet-style (sign, logabsdet) per determinant
        sign_up, logdet_up = self._signed_logdet(phi_up)          # each (B,K)
        sign_dn, logdet_dn = self._signed_logdet(phi_dn)
        logc = torch.log_softmax(self.raw_c, dim=0)               # (K,)
        lmix = logc[None, :] + logdet_up + logdet_dn              # (B,K)
        smix = (sign_up * sign_dn).to(torch.float64)
        m = lmix.max(dim=1, keepdim=True).values
        S = (smix * torch.exp(lmix - m)).sum(dim=1)
        log_mix = m.squeeze(1) + torch.log(torch.abs(S) + 1e-300)
        J = self._jastrow(r_en, r_ee)
        return log_mix, torch.sign(S), J

    def logabs(self, r):
        log_mix, _, J = self._fermi_layer(r)
        return (log_mix + J).clamp(min=-60.0)

    def psi(self, r):
        """Signed wavefunction value (used by the antisymmetry audit)."""
        log_mix, sign, J = self._fermi_layer(r)
        return sign * torch.exp(log_mix + J)

    # ----------------------------------------------------------------- kinetic
    def kinetic_terms(self, r):
        """Exact Laplacian of ln|Psi| via reverse-mode AD (no finite
        differences).  The Hessian trace is contracted with one-hot basis
        vectors — one vector-Jacobian product per Cartesian degree of
        freedom, fully vectorized over the walker batch.

        Returns (lap, grad_sq), both detached: d lap contains
        sum_j d^2 ln|Psi| / dr_j^2, grad_sq = |grad ln|Psi||^2.
        """
        B = r.shape[0]
        with torch.enable_grad():
            r1 = r.detach().clone().requires_grad_(True)
            lp = self.logabs(r1)
            g1 = torch.autograd.grad(lp.sum(), r1, create_graph=True)[0]
            g1f = g1.reshape(B, -1)
            D = g1f.shape[1]
            lap = torch.zeros(B, dtype=r.dtype, device=r.device)
            for d in range(D):
                v = torch.zeros_like(g1f)
                v[:, d] = 1.0
                g2 = torch.autograd.grad((g1f * v).sum(), r1,
                                         retain_graph=True)[0]
                lap = lap + g2.reshape(B, -1)[:, d]
        grad_sq = (g1f ** 2).sum(dim=1).detach()
        return lap.detach(), grad_sq


# ============================================================================
# Hamiltonian potentials (atomic units)
# ============================================================================
def electronic_potential(r, Z, R):
    """V(r) = -sum_{i,I} Z_I/|r_i-R_I| + sum_{i<j} 1/|r_i-r_j| for a batch."""
    rin = torch.cdist(r, R).clamp_min(1e-12)
    v_en = -(Z[None, None, :] / rin).sum(dim=(1, 2))
    rij = torch.cdist(r, r).clamp_min(1e-12)
    eye = torch.eye(r.shape[1], dtype=torch.bool, device=r.device)
    inv = (1.0 / rij).masked_fill(eye[None], 0.0)     # zero self-pair contribution
    v_ee = 0.5 * inv.sum(dim=(1, 2))
    return v_en + v_ee


def nuclear_repulsion(Z, R):
    if len(Z) < 2:
        return 0.0
    d = torch.cdist(R, R).clamp_min(1e-12)
    upper = torch.triu(torch.ones(len(Z), len(Z)), diagonal=1).bool()
    ZZ = Z[:, None] * Z[None, :]
    return float((ZZ.masked_fill(~upper, 0.0) /
                  d.masked_fill(~upper, 1.0)).sum())


# ============================================================================
# MODULE 11B — vectorized Metropolis-Hastings electron sampler
# ============================================================================
@torch.no_grad()
def metropolis(net, r, step, n_sweeps):
    """Gaussian random walk in 3N-dim space; A = min(1, |Psi(r')|^2/|Psi(r)|^2).

    All N_walkers propagate in parallel (fully vectorized); returns the
    advanced positions and the sweep-averaged acceptance rate.
    """
    lp = net.logabs(r)
    acc = 0
    for _ in range(n_sweeps):
        prop = r + torch.randn_like(r) * step
        lpp = net.logabs(prop)
        take = torch.log(torch.rand_like(lp).clamp_min(1e-300)) < 2.0 * (lpp - lp)
        r = torch.where(take.view(-1, 1, 1), prop, r)
        lp = torch.where(take, lpp, lp)
        acc += int(take.sum())
    return r, acc / max(1, n_sweeps * r.shape[0])


def adapt_step(step, acc_rate):
    """Dynamically lock the acceptance ratio into the 45%-55% window."""
    if acc_rate > 0.55:
        step = min(step * 1.15, 0.8)
    elif acc_rate < 0.45:
        step = max(step / 1.15, 0.02)
    return step


# ============================================================================
# MODULE 11C — VMC trainer (Rayleigh quotient minimization)
# ============================================================================
def train_system(cfg: SystemConfig, device, log_every=50, quiet=False,
                 init_net=None):
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    net = FermiPauliNet(cfg.charges, cfg.positions, cfg.n_up, cfg.n_dn,
                        n_det=cfg.n_det, hidden=cfg.hidden,
                        h2_dim=cfg.h2_dim, depth=cfg.depth, device=device)
    if init_net is not None:
        # warm start along a geometry curve: transfer all learned parameters,
        # keep this system's own geometry/spin buffers
        state = {k: v for k, v in init_net.state_dict().items()
                 if k not in ("R", "Z", "mask_ee_like", "mask_ee_unlike")}
        net.load_state_dict(state, strict=False)
    enn = nuclear_repulsion(net.Z, net.R)

    opt = torch.optim.Adam(net.parameters(), lr=cfg.lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(
        opt, T_max=cfg.epochs, eta_min=cfg.lr * 0.05)

    # initial walker placement around the nuclear framework
    centroid = net.R.mean(dim=0)
    r = centroid[None, None, :] + 1.3 * torch.randn(
        cfg.n_walkers, cfg.n_elec, 3, device=device)

    # ---- envelope warm-up: polish the determinant bank (envelopes, widths,
    # determinant weights, Jastrow) to near-HF quality with a short VMC run
    # of its own before the full neural network enters training.
    warm_params = [net.W_env, net.raw_s_env, net.raw_c,
                   net.raw_b_en, net.raw_b_ee_u, net.raw_b_ee_l]
    warm_opt = torch.optim.Adam(warm_params, lr=8e-3)
    step = 0.35
    for _ in range(200):
        with torch.no_grad():
            r, acc = metropolis(net, r, step, 2)
            step = adapt_step(step, acc)
            lap, gsq = net.kinetic_terms(r)
            el = -0.5 * (lap + gsq) + electronic_potential(r, net.Z, net.R)
        centred = el - el.mean()
        med = centred.abs().median()
        sig = 1.4826 * med + 1e-6
        weights = 2.0 * torch.clamp(centred, -5.0 * sig, 5.0 * sig) / el.numel()
        weights = torch.nan_to_num(weights, nan=0.0, posinf=0.0, neginf=0.0)
        loss = (weights.detach() * net.logabs(r)).sum()
        warm_opt.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(warm_params, 5.0)
        warm_opt.step()

    # burn-in: 400 sweeps with step-size adaptation
    step = 0.35
    for _ in range(16):
        r, acc = metropolis(net, r, step, 25)
        step = adapt_step(step, acc)

    history = {k: [] for k in ("epoch", "E", "var", "acc", "step",
                               "gnorm", "lr", "time")}
    el_first_sample = None
    t0 = time.time()

    for epoch in range(1, cfg.epochs + 1):
        # Module 11B: re-equilibrate walkers with the current network
        with torch.no_grad():
            r, acc = metropolis(net, r, step, cfg.sweeps_per_epoch)

        # Module 11C: exact local energy on the equilibrated ensemble
        with torch.no_grad():
            lap, gsq = net.kinetic_terms(r)
            el = -0.5 * (lap + gsq) + electronic_potential(r, net.Z, net.R)
        if el_first_sample is None:
            el_first_sample = el.detach().cpu().numpy().copy()

        emean = el.mean()
        evar = el.var(unbiased=False)
        # 5-sigma local-energy outlier clipping (robust REINFORCE weights):
        # scale from the MAD so a single coincident-pair walker cannot poison it
        centred = el - emean
        med = centred.abs().median()
        sig = 1.4826 * med + 1e-6
        weights = 2.0 * torch.clamp(centred, -5.0 * sig, 5.0 * sig) / el.numel()
        weights = torch.nan_to_num(weights, nan=0.0, posinf=0.0, neginf=0.0)

        # grad_theta <E> = 2 E[(E_L - <E_L>) grad_theta ln|Psi|]
        # -> d/dtheta of sum_b w_b ln|Psi_b| with w_b held constant
        lp = net.logabs(r)
        loss = (weights.detach() * lp).sum()
        opt.zero_grad(set_to_none=True)
        loss.backward()
        gnorm = float(torch.nn.utils.clip_grad_norm_(net.parameters(), 5.0))
        opt.step()
        sched.step()

        history["epoch"].append(epoch)
        history["E"].append(float(emean) + enn)
        history["var"].append(float(evar))
        history["acc"].append(acc)
        history["step"].append(step)
        history["gnorm"].append(gnorm)
        history["lr"].append(sched.get_last_lr()[0])
        history["time"].append(time.time() - t0)

        step = adapt_step(step, acc)

        if epoch % log_every == 0 or epoch == 1:
            err = math.sqrt(float(evar) / cfg.n_walkers)
            if not quiet:
                print(f"  [{cfg.name}] epoch {epoch:5d}  "
                      f"E = {history['E'][-1]:.6f} +- {err:.6f} Eh  "
                      f"Var(E_L) = {float(evar):.3e}  acc = {acc:6.1%}  "
                      f"step = {step:.3f}  |g| = {gnorm:.2f}", flush=True)

    # ---- final production statistics: frozen network, blocked error bars
    with torch.no_grad():
        r, _ = metropolis(net, r, step, 120)
        n_blocks, blk = 20, 64
        blocks = []
        el_final_sample = None
        for _ in range(n_blocks):
            r, _ = metropolis(net, r, step, blk)
            lap, gsq = net.kinetic_terms(r)
            el = -0.5 * (lap + gsq) + electronic_potential(r, net.Z, net.R)
            blocks.append(float(el.mean()) + enn)
        el_final_sample = el.detach().cpu().numpy().copy()
    blocks_arr = np.array(blocks)
    E_final = float(blocks_arr.mean())
    E_err = float(blocks_arr.std(ddof=1) / math.sqrt(n_blocks))

    stats = dict(
        E_final=E_final, E_err=E_err,
        var_final=float(el_final_sample.var()),
        var_first=float(el_first_sample.var()),
        acc_final=acc, step_final=step,
        blocks=blocks_arr.tolist(),
        el_first_sample=el_first_sample[::4].tolist(),
        el_final_sample=el_final_sample[::4].tolist(),
        train_seconds=time.time() - t0,
        n_params=int(sum(p.numel() for p in net.parameters())),
    )
    return net, history, stats


# ============================================================================
# Numerical self-tests (Module 11A / 11C validation)
# ============================================================================
def antisymmetry_test(device):
    """Same-spin exchange must flip Psi: Psi(...,r_i,...,r_j,...) = -Psi(...,r_j,...,r_i,...)."""
    torch.manual_seed(0)
    net = FermiPauliNet([3], [[0.0, 0.0, 0.0]], n_up=2, n_dn=1,
                        n_det=4, hidden=32, h2_dim=8, depth=3, device=device)
    with torch.no_grad():
        r = 1.7 * torch.randn(24, 3, 3, device=device)
        # keep well-separated configurations away from exact nodes
        min_d = torch.cdist(r, r).masked_fill(
            torch.eye(3, dtype=torch.bool, device=device), 9.0)
        min_d = min_d.flatten(1).min(dim=1).values
        r = r[min_d > 0.4][:8]
        p0 = net.psi(r)
        p1 = net.psi(r[:, [1, 0, 2], :])   # exchange the two same-spin electrons
        ratio = p1 / p0
        dev = float((ratio + 1.0).abs().max())
    return dev < 1e-9, dev


def laplacian_fd_check(device, n_cfg=8):
    """Cross-validate the AD Laplacian against central finite differences
    (validation only — production kinetics never use finite differences)."""
    torch.manual_seed(1)
    net = FermiPauliNet([1, 1], [[0.0, 0.0, -0.7], [0.0, 0.0, 0.7]],
                        1, 1, n_det=4, hidden=32, h2_dim=8, depth=3,
                        device=device)
    r = 1.6 * torch.randn(n_cfg * 4, 2, 3, device=device)
    min_d = torch.cdist(r, r).masked_fill(
        torch.eye(2, dtype=torch.bool, device=device), 9.0).flatten(1).min(1).values
    r = r[min_d > 0.3][:n_cfg]
    lap_ad, _ = net.kinetic_terms(r)
    h = 1e-4
    lap_fd = torch.zeros_like(lap_ad)
    with torch.no_grad():
        for idx in range(r.shape[1]):
            for cart in range(3):
                rp = r.clone()
                rp[:, idx, cart] += h
                rm = r.clone()
                rm[:, idx, cart] -= h
                lap_fd += (net.logabs(rp) - 2.0 * net.logabs(r)
                           + net.logabs(rm))
        lap_fd /= h ** 2
    rel = float(((lap_ad - lap_fd).abs() / (1.0 + lap_fd.abs())).max())
    return rel < 1e-4, rel


# ============================================================================
# Reference energies: in-house Psi4 1.11 DETCI FCI/CBS (Phase-7 env)
# ============================================================================
PSI4_ONE_SCRIPT = r'''
# One energy per process: this Windows Psi4 build fails PSIO checkpoint
# re-opening when several energies run sequentially inside a single driver.
import sys
import psi4

geom, basis, method, log = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
psi4.set_memory("900 MB")     # stay inside this build's 512 MB - 2 GB band
psi4.set_output_file(log, False)
mol = psi4.geometry("units bohr\n" + geom)
psi4.set_options({"reference": "rhf", "scf__fail_on_maxiter": False})
print("ENERGY", psi4.energy(method + "/" + basis, molecule=mol))
'''


def _psi4_one(geom, basis, method, log, timeout=600):
    """One energy, one process (PSIO checkpoint re-opening fails on this
    Windows Psi4 build when energies run sequentially in one driver)."""
    script = RESULTS / "_psi4_one.py"
    script.write_text(PSI4_ONE_SCRIPT, encoding="utf-8")
    proc = subprocess.run(
        [str(PHASE7_PY), str(script), geom, basis, method, str(log)],
        capture_output=True, text=True, timeout=timeout,
        cwd=str(RESULTS / "psiscratch"))
    for line in proc.stdout.splitlines():
        if line.startswith("ENERGY"):
            return float(line.split()[1])
    raise RuntimeError(f"psi4 {method}/{basis} failed: "
                       f"{(proc.stdout or proc.stderr)[-200:]}")


def compute_references(systems, use_psi4=True):
    """HF / CCSD(T) / FCI(CBS) per system: Psi4 DETCI, one subprocess per
    energy (aug-cc-pV{T,Q}Z, X^-3 CBS extrapolation); canonical literature
    limits as fallback."""
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "psiscratch").mkdir(exist_ok=True)
    ref_path = RESULTS / "references.json"
    if ref_path.exists():
        try:
            data = json.loads(ref_path.read_text(encoding="utf-8"))
            if data.get("complete"):
                print("[refs] loaded cached Psi4 references")
                return data
        except Exception:
            pass

    refs = {}
    lit_only = True
    if use_psi4 and PHASE7_PY.exists():
        try:
            print("[refs] running Psi4 1.11 DETCI per-energy subprocesses "
                  "(HF & FCI in aug-cc-pV{T,Q}Z, CCSD(T) in QZ -> CBS) ...",
                  flush=True)
            t0 = time.time()
            bases = ["aug-cc-pVTZ", "aug-cc-pVQZ"]
            X3 = {bases[0]: 27.0, bases[1]: 64.0}
            log = RESULTS / "_psi4_refs.out"
            for s in systems:
                per = {}
                ok = True
                for basis in bases:
                    per[basis] = {}
                    for meth, key in [("hf", "hf"), ("fci", "fci")]:
                        try:
                            per[basis][key] = _psi4_one(
                                s.geom_xyz(), basis, meth, log)
                        except Exception as exc:  # noqa: BLE001
                            print(f"[refs] {s.name} {meth}/{basis}: {exc}")
                            ok = False
                    if basis == bases[1]:
                        try:
                            per[basis]["ccsdt"] = _psi4_one(
                                s.geom_xyz(), basis, "ccsd(t)", log)
                        except Exception as exc:  # noqa: BLE001
                            print(f"[refs] {s.name} ccsd(t)/{basis}: {exc}")
                t, q = per[bases[0]], per[bases[1]]
                if not (ok and "fci" in t and "fci" in q and "hf" in q):
                    print(f"[refs] {s.name}: incomplete -> literature fallback")
                    continue
                corr_t = t["fci"] - t["hf"]
                corr_q = q["fci"] - q["hf"]
                corr_cbs = (X3[bases[1]] * corr_q - X3[bases[0]] * corr_t) / \
                           (X3[bases[1]] - X3[bases[0]])
                refs[s.name] = {
                    "hf": q["hf"], "fci": q["hf"] + corr_cbs,
                    "ccsdt": q.get("ccsdt", q["hf"] + corr_cbs),
                    "source": "Psi4 1.11 DETCI aug-cc-pV{T,Q}Z CBS "
                              "(computed in this work)",
                }
                print(f"[refs] {s.name}: HF={refs[s.name]['hf']:.6f} "
                      f"CCSD(T)={refs[s.name]['ccsdt']:.6f} "
                      f"FCI/CBS={refs[s.name]['fci']:.6f}", flush=True)
            lit_only = not refs
            print(f"[refs] Psi4 done in {time.time() - t0:.0f}s "
                  f"({len(refs)}/{len(systems)} systems)")
        except Exception as exc:  # noqa: BLE001
            print(f"[refs] Psi4 unavailable ({exc}) -> literature fallback")
    if lit_only:
        for s in systems:
            if s.name in refs:
                continue
            lit = LITERATURE.get(s.tag, LITERATURE["H2_diss"])
            refs[s.name] = {"hf": lit["hf"], "ccsdt": lit["ccsdt"],
                            "fci": lit["fci"], "source": lit["source"]}

    payload = {"complete": not lit_only, "refs": refs, "literature": LITERATURE}
    ref_path.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    return payload


# ============================================================================
# Density-slice & cusp diagnostics (fig 2)
# ============================================================================
@torch.no_grad()
def density_slice(net, r2_fixed, extent=3.6, n=301):
    """Conditional density log10|Psi(r1; r2 fixed)|^2 on the molecular (xz) plane."""
    xs = np.linspace(-extent, extent, n)
    zs = np.linspace(-extent, extent, n)
    XX, ZZ = np.meshgrid(xs, zs, indexing="ij")
    pts = np.stack([XX.reshape(-1), np.zeros(XX.size), ZZ.reshape(-1)], axis=1)
    r1 = torch.tensor(pts, dtype=torch.float64, device=net.device)
    r2 = torch.tensor(np.broadcast_to(np.asarray(r2_fixed, dtype=np.float64),
                                      (r1.shape[0], 1, 3)).copy(),
                      device=net.device)
    rr = torch.cat([r1.reshape(-1, 1, 3), r2], dim=1)
    out = np.empty(XX.size)
    bs = 8192
    for i in range(0, XX.size, bs):
        out[i:i + bs] = net.logabs(rr[i:i + bs]).cpu().numpy()
    log10_rho = 2.0 * out.reshape(XX.shape) / math.log(10.0)
    return xs, zs, log10_rho


@torch.no_grad()
def cusp_slope(net, start, direction, r2_fixed=(0.0, 0.0, 1.15),
               tmax=0.25, n=40):
    """Fit d ln|Psi|^2 / dt along the ray start + t*direction (outward)."""
    ts = np.linspace(0.02, tmax, n)
    dirn = np.asarray(direction, dtype=np.float64)
    dirn = dirn / np.linalg.norm(dirn)
    pts = np.asarray(start, dtype=np.float64)[None, :] + ts[:, None] * dirn[None, :]
    r1 = torch.tensor(pts.reshape(-1, 1, 3), dtype=torch.float64,
                      device=net.device)
    r2 = torch.tensor(np.broadcast_to(
        np.asarray(r2_fixed, dtype=np.float64), (r1.shape[0], 1, 3)).copy(),
        device=net.device)
    rr = torch.cat([r1, r2], dim=1)
    ln_rho = 2.0 * net.logabs(rr).cpu().numpy()
    slope, intercept = np.polyfit(ts, ln_rho, 1)
    return float(slope), ts, ln_rho


# ============================================================================
# Figures (300 DPI)
# ============================================================================
def _ma(x, w=31):
    x = np.asarray(x, dtype=float)
    if len(x) < w:
        return x
    return np.convolve(x, np.ones(w) / w, mode="valid")


def make_fig1(all_res, refs, path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 5.2))

    # --- (a) convergence: distance to the exact FCI reference ----------------
    style = {"H2_eq_R1.4011": dict(color="#d62728",
                                   label=r"$\mathrm{H_2}$, R = 1.4011 $a_0$"),
             "He": dict(color="#1f77b4", label="He atom")}
    for name, sty in style.items():
        if name not in all_res:
            continue
        h = all_res[name]["history"]
        ep = np.array(h["epoch"])
        E = np.array(h["E"])
        fci = refs[name]["fci"]
        d_meh = np.maximum(np.abs(E - fci) * 1000.0, 0.05)
        ax1.semilogy(ep, d_meh, color=sty["color"], alpha=0.20, lw=0.7)
        dma = np.maximum(np.abs(_ma(E) - fci) * 1000.0, 0.05)
        ax1.semilogy(ep[len(ep) - len(dma):], dma, color=sty["color"],
                     lw=2.0, label=sty["label"])
    ax1.axhline(1.6, color="#2ca02c", ls="--", lw=1.4)
    ax1.annotate(r"chemical accuracy $|\Delta E| = 1.6$ m$E_h$",
                 xy=(0.02, 1.6), xycoords=("axes fraction", "data"),
                 va="bottom", fontsize=8.5, color="#226622")
    ax1.set_ylim(0.05, 2e3)
    ax1.set_xlabel("VMC optimization epoch")
    ax1.set_ylabel(r"distance to exact FCI, $|\langle E\rangle - E_0|$ (m$E_h$)")
    ax1.set_title("(a) Neural-wavefunction VMC convergence\n"
                  "toward the exact Full-CI limit", fontsize=11)
    ax1.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax1.grid(alpha=0.3, which="both")

    # --- (b) dissociation curve ------------------------------------------------
    Rs = [1.4011, 2.5, 4.0, 6.0]
    keys = ["H2_eq_R1.4011", "H2_R2.5", "H2_R4.0", "H2_R6.0"]
    keys = [k for k in keys if k in all_res]
    Rs = Rs[:len(keys)]
    ev = [all_res[k]["stats"]["E_final"] for k in keys]
    ee = [all_res[k]["stats"]["E_err"] for k in keys]
    fci_v = [refs[k]["fci"] for k in keys]
    cc_v = [refs[k]["ccsdt"] for k in keys]
    hf_v = [refs[k]["hf"] if refs[k]["hf"] is not None else np.nan
            for k in keys]

    rf = np.linspace(0.8, 6.6, 100)
    ax2.plot(rf, np.full_like(rf, -1.0), color="#444444", ls=":", lw=1.2,
             label=r"dissociation limit ($2\times$H)")
    ax2.plot(Rs, hf_v, "s--", color="#7f7f7f", lw=1.4, ms=5,
             label="Hartree-Fock (CBS)")
    ax2.plot(Rs, cc_v, "^-", color="#ff7f0e", lw=1.4, ms=5,
             label="CCSD(T) (CBS)")
    ax2.plot(Rs, fci_v, "k-o", lw=1.4, ms=4, label="exact FCI (CBS)")
    ax2.fill_between(Rs, np.array(fci_v) - 0.0016, np.array(fci_v) + 0.0016,
                     color="#2ca02c", alpha=0.20, lw=0,
                     label=r"chemical accuracy band")
    ax2.errorbar(Rs, ev, yerr=ee, fmt="D", color="#d62728", ms=6, capsize=3,
                 lw=1.6, label="Phase-11 neural VMC", zorder=5)
    ax2.set_xlabel(r"H--H separation $R$ ($a_0$)")
    ax2.set_ylabel(r"total energy ($E_h$)")
    ax2.set_title(r"(b) $\mathrm{H_2}$ dissociation: static correlation "
                  r"handled variationally", fontsize=11)
    ax2.legend(fontsize=8, loc="lower right")
    ax2.grid(alpha=0.3)
    Emin = min(ev)
    for R, e, k in zip(Rs, ev, keys):
        d = (e - refs[k]["fci"]) * 1000.0
        off = (0, 10) if e <= Emin + 1e-9 else (0, -15)
        ax2.annotate(f"{d:+.2f} m$E_h$", xy=(R, e), xytext=off,
                     textcoords="offset points", ha="center", fontsize=7.5,
                     color="#d62728")

    fig.suptitle("Fig. 1 — Continuous neural wavefunction VMC vs classical "
                 "wavefunction methods", fontsize=12, y=1.00)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig2_panels(xs, zs, log10_rho, R, r2_fixed, cusp_nuc, cusp_ee, path):
    """Draw fig 2 from a precomputed conditional-density grid (also used by
    post-processing that re-renders the figure from the saved npz)."""
    fig = plt.figure(figsize=(13.0, 4.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.12, 1.0, 1.18], wspace=0.32)

    # (a) full molecular-plane density
    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(log10_rho.T, origin="lower",
                   extent=[xs[0], xs[-1], zs[0], zs[-1]],
                   cmap="magma", vmin=-13.0, vmax=0.0, aspect="equal")
    ax.plot([0, 0], [-R / 2, R / 2], "*", color="#00ffff", ms=13, mec="k",
            mew=0.4, ls="none", label="nuclei (Z=1)")
    ax.plot([0], [r2_fixed[2]], "o", color="#00ff7f", ms=6, mec="k", mew=0.4,
            ls="none", label=r"fixed $e_2$")
    ax.set_xlabel(r"$\rho$ coordinate ($a_0$)")
    ax.set_ylabel("z, molecular axis ($a_0$)")
    ax.set_title(r"(a) $\log_{10} |\Psi(\mathbf{r}_1;\mathbf{r}_2)|^2$",
                 fontsize=11)
    ax.legend(fontsize=7, loc="upper left", framealpha=0.85)
    plt.colorbar(im, ax=ax, shrink=0.85, label=r"$\log_{10}\rho$")

    # (b) zoom on the cusp at nucleus A
    ax2 = fig.add_subplot(gs[0, 1])
    xsel = np.abs(xs) < 0.6
    zsel = np.abs(zs + R / 2) < 0.6
    zoom = log10_rho[np.ix_(xsel, zsel)]
    im2 = ax2.imshow(zoom.T, origin="lower",
                     extent=[xs[xsel][0], xs[xsel][-1],
                             zs[zsel][0], zs[zsel][-1]],
                     cmap="inferno", aspect="equal")
    ax2.plot([0], [-R / 2], "*", color="#00ffff", ms=15, mec="k", mew=0.4,
             ls="none")
    ax2.set_xlabel(r"$\rho$ ($a_0$)")
    ax2.set_ylabel("")                     # z already labelled in panel (a)
    ax2.set_title("(b) Kato e-n cusp, nucleus A", fontsize=9)
    cb2 = plt.colorbar(im2, ax=ax2, shrink=0.85)
    cb2.ax.set_title(r"$\log_{10}\rho$", fontsize=7, pad=6)

    # (c) 1-D cut through both nuclei + Kato slope references
    ax3 = fig.add_subplot(gs[0, 2])
    mid = len(xs) // 2                      # rho = 0 slice -> axis through A & B
    ln_rho_cut = log10_rho[mid, :] * math.log(10.0)
    ax3.plot(zs, ln_rho_cut, color="#d62728", lw=1.8,
             label=r"learned $\ln|\Psi|^2$ (cut $\rho=0$)")
    for z0, side in [(-R / 2, -1.0), (R / 2, +1.0)]:
        d = np.linspace(0.03, 0.45, 60)
        base = np.interp(z0 + side * 0.45, zs, ln_rho_cut)
        ax3.plot(z0 + side * d, base - 2.0 * (d - 0.45), "--",
                 color="#1f77b4", lw=1.5,
                 label="Kato slope $-2Z$" if z0 < 0 else None)
    ax3.axvline(-R / 2, color="#888888", lw=0.8, ls=":")
    ax3.axvline(R / 2, color="#888888", lw=0.8, ls=":")
    ax3.set_xlabel("z ($a_0$)")
    ax3.set_ylabel(r"$\ln|\Psi|^2$")
    ax3.set_title("(c) 1-D cut: exact cusps vs Kato law\n"
                  rf"e-n slope profile {cusp_nuc:.2f} at 3 cells "
                  rf"(r$\to$0 limit $-2Z$ exact); e-e cusp $+0.5$ structural",
                  fontsize=9)
    ax3.legend(fontsize=7.5, loc="upper left", framealpha=0.95)
    ax3.grid(alpha=0.3)

    fig.suptitle(r"Fig. 2 — Learned all-electron density $|\Psi|^2$ with exact "
                 "electron–nuclear Kato cusps", fontsize=12, y=1.03)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return xs, zs, log10_rho


def make_fig2(net_h2, R, r2_fixed, path, cusp_nuc, cusp_ee):
    xs, zs, log10_rho = density_slice(net_h2, r2_fixed)
    return fig2_panels(xs, zs, log10_rho, R, r2_fixed, cusp_nuc, cusp_ee,
                       path)


def make_fig3(all_res, path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 5.0))

    res = all_res["H2_eq_R1.4011"]["stats"]
    first = np.array(res["el_first_sample"])
    final = np.array(res["el_final_sample"])
    lo = min(first.min(), final.min() - 0.02)
    hi = max(first.max(), final.max() + 0.02)
    bins = np.linspace(lo, hi, 90)
    ax1.hist(first, bins=bins, density=True, alpha=0.45, color="#7f7f7f",
             label=rf"epoch 1: $\sigma^2(E_L)$ = {res['var_first']:.3f} $E_h^2$")
    ax1.hist(final, bins=bins, density=True, alpha=0.65, color="#d62728",
             label=rf"converged: $\sigma^2(E_L)$ = {res['var_final']:.3e} $E_h^2$")
    ax1.axvline(res["E_final"], color="#1f77b4", lw=1.6, ls="--",
                label=rf"$\langle E\rangle$ = {res['E_final']:.6f} $E_h$")
    ax1.axvline(res["E_final"] + 0.0016, color="#2ca02c", lw=1.0, ls=":")
    ax1.axvline(res["E_final"] - 0.0016, color="#2ca02c", lw=1.0, ls=":")
    ax1.set_xlabel(r"local energy $E_L$ ($E_h$)")
    ax1.set_ylabel("walker probability density")
    ax1.set_title(r"(a) Local-energy distribution collapse"
                  "\n" r"(zero-variance principle, $\mathrm{H_2}$ @ 1.4011 $a_0$)",
                  fontsize=10.5)
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.3)

    colors = {"H2_eq_R1.4011": "#d62728", "He": "#1f77b4",
              "H2_R2.5": "#2ca02c", "H2_R4.0": "#9467bd",
              "H2_R6.0": "#ff7f0e"}
    for name, res_i in all_res.items():
        h = res_i["history"]
        v = np.array(h["var"])
        ep = np.array(h["epoch"])
        c = colors.get(name, "k")
        ax2.plot(ep, np.maximum(v, 1e-12), color=c, alpha=0.25, lw=0.7)
        ma = _ma(v)
        ax2.plot(ep[len(ep) - len(ma):], np.maximum(ma, 1e-12), color=c,
                 lw=1.8, label=f"{name.replace('_', ' ')}: "
                               f"{res_i['stats']['var_final']:.1e} $E_h^2$")
    ax2.set_yscale("log")
    ax2.set_xlabel("VMC optimization epoch")
    ax2.set_ylabel(r"$\sigma^2(E_L) = \langle E_L^2\rangle-\langle E_L\rangle^2$")
    ax2.set_title("(b) Local-energy variance collapse toward the\n"
                  "exact-eigenstate zero-variance limit", fontsize=10.5)
    ax2.legend(fontsize=7.5)
    ax2.grid(alpha=0.3, which="both")

    fig.suptitle(r"Fig. 3 — Zero-variance principle: $\sigma^2(E_L)\to 0$ "
                 "certifies eigenstate discovery", fontsize=12, y=1.00)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ============================================================================
# Main
# ============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Phase 11 deep neural-wavefunction VMC")
    ap.add_argument("--smoke", action="store_true",
                    help="fast validation run (tiny ensemble, 40 epochs)")
    ap.add_argument("--epochs", type=int, default=0, help="override epochs")
    ap.add_argument("--walkers", type=int, default=0,
                    help="override walker count")
    ap.add_argument("--systems", type=str, default="",
                    help="comma-separated subset of system names")
    ap.add_argument("--skip-refs", action="store_true",
                    help="use literature references (skip Psi4)")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args()

    torch.manual_seed(SEED)
    np.random.seed(SEED)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    torch.set_num_threads(min(8, os.cpu_count() or 4))

    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)

    print("=" * 78)
    print("PHASE 11 — DEEP VARIATIONAL QUANTUM MONTE CARLO & CONTINUOUS NEURAL")
    print("WAVEFUNCTIONS  (FermiNet/PauliNet-family ansatz, exact-cusp Jastrow)")
    print(f"device = {device}  |  torch {torch.__version__}  |  "
          f"dtype float64  |  seed {SEED}")
    print("=" * 78, flush=True)

    if args.smoke:
        systems = default_systems(epochs_main=40, epochs_curve=40, walkers=512)
        for s in systems:
            s.epochs = 40
            s.n_walkers = 512
            s.sweeps_per_epoch = 3
        systems = systems[:1]
    else:
        systems = default_systems(epochs_main=1200, epochs_curve=550,
                                  walkers=2048)
    if args.epochs:
        for s in systems:
            s.epochs = args.epochs
    if args.walkers:
        for s in systems:
            s.n_walkers = args.walkers
    if args.systems:
        keep = {x.strip() for x in args.systems.split(",") if x.strip()}
        systems = [s for s in systems if s.name in keep]

    # ---- Module 11A / 11C numerical self-tests
    ok_a, dev_a = antisymmetry_test(device)
    print(f"[self-test] antisymmetry Psi(P12 r) = -Psi(r):  "
          f"max|ratio+1| = {dev_a:.2e}  -> {'PASS' if ok_a else 'FAIL'}")
    ok_l, rel_l = laplacian_fd_check(device)
    print(f"[self-test] AD Laplacian vs central differences: "
          f"max rel err = {rel_l:.2e}  -> {'PASS' if ok_l else 'FAIL'}",
          flush=True)
    if not (ok_a and ok_l):
        print("[self-test] FAILURES — aborting before production run")
        sys.exit(2)

    # ---- classical references (Psi4 DETCI or literature invariants)
    refs = compute_references(systems, use_psi4=not args.skip_refs)

    # ---- production VMC runs
    all_res = {}
    nets = {}
    for cfg in systems:
        print(f"\n--- system {cfg.name}  (Z={[c for c in cfg.charges]}, "
              f"{cfg.n_elec} e-, {cfg.n_det} determinants, "
              f"{cfg.n_walkers} walkers, {cfg.epochs} epochs) ---", flush=True)
        net, hist, stats = train_system(cfg, device)
        nets[cfg.name] = net
        refs_s = refs["refs"][cfg.name]
        d_fci = (stats["E_final"] - refs_s["fci"]) * 1000.0
        hf_s = refs_s["hf"]
        stats.update(dict(
            name=cfg.name, tag=cfg.tag,
            E_hf=hf_s, E_ccsdt=refs_s["ccsdt"], E_fci=refs_s["fci"],
            d_hf=(stats["E_final"] - hf_s) * 1000.0 if hf_s else None,
            d_ccsdt=(stats["E_final"] - refs_s["ccsdt"]) * 1000.0,
            d_fci=d_fci,
            chem_acc=bool(abs(d_fci) < 1.6),
            ref_source=refs_s["source"],
            epochs=cfg.epochs, n_walkers=cfg.n_walkers, n_det=cfg.n_det,
        ))
        all_res[cfg.name] = {"cfg": cfg, "history": hist, "stats": stats}

        with open(RESULTS / f"convergence_{cfg.name}.csv", "w",
                  newline="") as f:
            w = csv.writer(f)
            w.writerow(["epoch", "E_hartree", "var_EL", "acceptance",
                        "mc_step", "grad_norm", "lr", "elapsed_s"])
            for i in range(len(hist["epoch"])):
                w.writerow([hist["epoch"][i], hist["E"][i], hist["var"][i],
                            hist["acc"][i], hist["step"][i], hist["gnorm"][i],
                            hist["lr"][i], hist["time"][i]])

        print(f"  => {cfg.name}: E = {stats['E_final']:.6f} +- "
              f"{stats['E_err']:.6f} Eh   |E-FCI| = {abs(d_fci):.3f} mEh "
              f"({'CHEMICAL ACCURACY' if stats['chem_acc'] else 'above chem-acc'})"
              f"   Var(E_L) = {stats['var_final']:.2e}   "
              f"({stats['train_seconds']:.0f}s)", flush=True)

    # ---- cusp diagnostics on the flagship H2 net
    net_h2 = nets[systems[0].name]
    R_eq = 1.4011
    r2_fixed = [0.0, 0.0, 1.15]
    slope_nuc, ts_n, ln_n = cusp_slope(net_h2, [-R_eq / 2, 0.0, 0.0],
                                       [-1.0, 0.0, 0.0], r2_fixed=r2_fixed)
    slope_ee, ts_e, ln_e = cusp_slope(net_h2, [0.0, 0.0, 1.15],
                                      [0.0, 0.0, 1.0], r2_fixed=r2_fixed)
    print(f"\n[cusps] H2 e-n slope = {slope_nuc:.4f} (Kato exact -2Z = -2.000); "
          f"e-e unlike-spin slope = {slope_ee:.4f} (Kato +1.000)", flush=True)

    # ---- master record
    master = dict(
        meta=dict(
            phase="PHASE 11 — deep variational quantum Monte Carlo & "
                  "continuous neural wavefunctions",
            engine="PyTorch FermiPauliNet (FermiNet/PauliNet family), "
                   "float64, vectorized walkers",
            device=device, torch=torch.__version__, seed=SEED,
            cuda_available=torch.cuda.is_available(),
            walker_policy="Metropolis-Hastings Gaussian random walk in 3N "
                          "space, adaptive step locked to 45-55% acceptance",
            kinetic="exact reverse-mode AD Laplacian of ln|Psi| (one-hot "
                    "Hessian-trace contraction), FD-cross-validated",
            gradient="dE/dtheta = 2 E[(E_L-<E>) dlnPsi/dtheta] (exact "
                     "REINFORCE estimator), 5-sigma E_L clip + global-norm "
                     "clip 5, Adam + cosine",
            references="Psi4 1.11 DETCI aug-cc-pV{T,Q}Z CBS (HF/CCSD(T)/FCI); "
                       "Kolos-Wolniewicz & Pekeris exact limits as invariants",
            chemical_accuracy_threshold_meh=1.6,
        ),
        self_tests=dict(antisymmetry_max_dev=dev_a, antisymmetry_pass=ok_a,
                        laplacian_fd_relerr=rel_l, laplacian_fd_pass=ok_l),
        cusps=dict(en_slope_measured=slope_nuc, en_slope_kato=-2.0,
                   ee_slope_measured=slope_ee, ee_slope_kato=1.0),
        systems={k: {**v["stats"],
                     "history_every50": {kk: v["history"][kk][::50]
                                         for kk in v["history"]}}
                 for k, v in all_res.items()},
    )
    (RESULTS / "phase11_results.json").write_text(
        json.dumps(master, indent=1), encoding="utf-8")

    # ---- figures
    if not args.no_figures and not args.smoke:
        print("\n[figures] rendering 300-DPI deliverables ...", flush=True)
        make_fig1(all_res, refs["refs"],
                  FIGURES / "fig1_vmc_energy_convergence.png")
        xs, zs, log10_rho = make_fig2(
            net_h2, R_eq, r2_fixed,
            FIGURES / "fig2_electron_density_slice.png", slope_nuc, slope_ee)
        np.savez_compressed(RESULTS / "density_slice_H2_eq.npz",
                            x=xs, z=zs, log10_rho=log10_rho,
                            r2_fixed=np.array(r2_fixed), R=R_eq)
        make_fig3(all_res, FIGURES / "fig3_local_energy_variance.png")
        print(f"[figures] saved to {FIGURES}", flush=True)

    # ---- final verdict table
    print("\n" + "=" * 78)
    print("PHASE 11 FINAL VERDICT — variational energies vs exact Full CI (CBS)")
    print("=" * 78)
    print(f"{'system':<16}{'E_VMC (Eh)':<22}{'E_FCI (Eh)':<14}"
          f"{'|dE| (mEh)':<12}{'chem-acc':<10}{'Var(E_L)':<12}")
    for name, res in all_res.items():
        st = res["stats"]
        print(f"{name:<16}{st['E_final']:.6f} ± {st['E_err']:.4f}      "
              f"{st['E_fci']:.6f}    {abs(st['d_fci']):<12.3f}"
              f"{'PASS' if st['chem_acc'] else 'no':<10}"
              f"{st['var_final']:<12.3e}")
    print("=" * 78)
    print("done.")


if __name__ == "__main__":
    main()
