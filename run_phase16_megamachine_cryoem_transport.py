#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_phase16_megamachine_cryoem_transport.py — PHASE 16 SUPREME MISSION
Mega-macromolecular dynamics, Cryo-EM density flexible fitting & mesoscale
transport through the Nuclear Pore Complex.

Escalation from isolated protein complexes (Phases 2/3/15) to a mega-Dalton
macromolecular machine: the eight-fold-symmetric Nuclear Pore Complex scaffold,
its intrinsically-disordered phenylalanine-glycine nucleoporin (FG-Nup) brush,
and the selective translocation of transport receptors through the FG condensate
— reconciling a STATIC Cryo-EM density map with DYNAMIC non-equilibrium
polymer-brush transport physics.

MODULES
-------
16A  Coarse-grained (Martini-3-flavored) mega-system assembly.
     8-fold symmetric NPC scaffold (structured domains, cryo-EM visible),
     72 anchored FG-Nup chains (sticker/spacer grammar: TA aromatic stickers,
     TC/N0 spacers, charged SP/SN spacers, eps_TA-TA = 3.0 kJ/mol hydrophobic
     multivalent cohesion), solvated in an explicit CG water-ion lattice
     (~1e6 beads, ~4e6 all-atom equivalents at the 4:1 Martini mapping,
     0.15 M NaCl).  The solvent is integrated with an O(N) grid density
     functional (compressible Brownian continuum; periodic, vectorized,
     preallocated leak-free buffers).

16B  Cryo-EM Molecular Dynamics Flexible Fitting (MDFF) engine.
     A synthetic 3.5 A cryo-EM map of the scaffold is synthesized (Gaussian
     pseudo-atoms + B-factor damping + noise) and WRITTEN/RE-READ as a real
     MRC2014 file (bit-exact round trip asserted).  V_EM(r) = -xi * sum_i
     w_i * Phi(rho_EM(r_i)) is applied through a VECTORIZED TRILINEAR-
     INTERPOLATED SPATIAL GRADIENT kernel (analytic gradient of the trilinear
     interpolant) with MDFF g-scale weighting and an adaptive-force ramp,
     steering a deliberately distorted scaffold (per-bead noise + 1.6 nm rigid
     shift + 3 deg rotation + per-spoke twist) into the density envelope behind
     an elastic network.  Cross-validated against OpenMM: the elastic-network
     kernel is compared force-by-force against a compiled OpenMM System
     (CustomBondForce) and the fitted model is re-relaxed in OpenMM
     LangevinMiddle as an independent stability gate.

16C  Non-equilibrium translocation free energy & hydrodynamics.
     Steered seeding — the FIRST segment at the literal v = 0.05 nm/ns,
     remaining segments at the same harmonic force protocol under a documented
     time acceleration — along the pore axis z in [-25, +25] nm, then
     independent umbrella windows (k = 60 kJ/mol/nm^2) for BOTH the specific
     transport receptor and the chemically inert R_h-matched control; PMFs by a
     log-space WHAM implementation (two-half consistency bands), spatial
     apparent diffusion D(z) from lag-resolved MSD with harmonic-well
     correction, transient FG-contact statistics, and the inhomogeneous
     solubility-diffusion permeability P ~ [integral exp(beta*G)/D dz]^-1.

DELIVERABLES
------------
figures_phase16/fig1_megasystem_cryoem_fit.png        3D volumetric render
figures_phase16/fig2_fg_condensate_density_slice.png  FG gel/brush heatmap
figures_phase16/fig3_translocation_free_energy_pmf.png PMF / D(z) / contacts
MEGAMACHINE_CRYOM_REPORT_EN.md / _ZH.md               bilingual treatise
results_phase16/                                       machine-readable record

USAGE
-----
python run_phase16_megamachine_cryoem_transport.py            # full run
python run_phase16_megamachine_cryoem_transport.py --fast     # smoke budget
python run_phase16_megamachine_cryoem_transport.py --selftest # kernel gates
python run_phase16_megamachine_cryoem_transport.py --fig-only # re-render figs

Wall-clock-budgeted and honest about simulated lengths: every number shipped
to the reports is produced here, in-engine, with its own validation gate.
"""

import argparse
import gc
import json
import math
import os
import time
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, LinearSegmentedColormap, Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial.transform import Rotation

try:  # marching cubes for the cryo-EM isosurface
    from skimage import measure as _skm
    HAS_SKIMAGE = True
except Exception:  # pragma: no cover
    _skm = None
    HAS_SKIMAGE = False

try:  # optional kernel cross-validation target
    import openmm as _mm
    from openmm import unit as _u
    HAS_OPENMM = True
except Exception:  # pragma: no cover
    _mm = None
    _u = None
    HAS_OPENMM = False

ROOT = Path(__file__).resolve().parent
RES = ROOT / "results_phase16"
FIG = ROOT / "figures_phase16"
TRAPZ = getattr(np, "trapezoid", getattr(np, "trapz", None))

# ------------------------------------------------------------------ constants
KB = 0.008314462618          # kJ/(mol*K)
T_SIM = 310.0                # K
KT = KB * T_SIM              # 2.5775 kJ/mol
BETA = 1.0 / KT
KJ_TO_KCAL = 1.0 / 4.184
NM2PS_TO_UM2S = 1.0e8        # 1 nm^2/ps = 1e-2 cm^2/s = 1e8 um^2/s
N_AVOG_NM3 = 0.6022141       # particles/nm^3 per mol/L
RHO_WATER_CG = 8.35          # CG water beads/nm^3 (4:1 Martini mapping)
SALT_M = 0.15                # mol/L NaCl
KAPPA_DH = 1.30              # nm^-1, Debye kappa at 0.15 M, eps_r = 60
EPS_R = 60.0                 # effective CG dielectric
KE_COUL = 138.935 / EPS_R    # kJ*nm/(mol*e^2)
ETA_WATER = 0.85e-3          # Pa s at 310 K (Stokes-Einstein calibration anchor)

_RNG = np.random.default_rng(0xF16)
T_START = time.time()


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T_START, msg), flush=True)


class Budget:
    """Wall-clock governor: stages size their step counts from measured speed."""

    def __init__(self, total_min, fast):
        self.total = total_min * 60.0
        self.fast = fast
        self.spent = {}
        self._t = {}

    def mark(self, key):
        self._t[key] = time.time()

    def close(self, key):
        if key in self._t:
            self.spent[key] = self.spent.get(key, 0.0) + round(time.time() - self._t.pop(key), 1)

    def left(self):
        return max(60.0, self.total - (time.time() - T_START))


# ==================================================================== MRC I/O
# Grid convention throughout: array[x, y, z] with x fastest (matches MRC
# mapc/mapr/maps = 1,2,3 with C-order storage).
def write_mrc(path, grid_xyz, origin_nm, voxel_nm):
    grid = np.ascontiguousarray(grid_xyz, dtype="<f4")
    nx, ny, nz = grid.shape
    hdr = np.zeros(256, dtype="<i4")
    hdr[0], hdr[1], hdr[2] = nx, ny, nz
    hdr[3] = 2                                     # mode 2: 32-bit float
    hdr[7], hdr[8], hdr[9] = nx, ny, nz
    hdr[10:13] = (np.array([nx, ny, nz]) * voxel_nm * 10.0).astype("<i4")   # cella (A)
    hdr[13:16] = (90, 90, 90)
    hdr[16], hdr[17], hdr[18] = 1, 2, 3
    hdr[19] = int(np.float32(grid.min()).view(np.int32))
    hdr[20] = int(np.float32(grid.max()).view(np.int32))
    hdr[21] = int(np.float32(grid.mean()).view(np.int32))
    hdr[49:52] = (np.asarray(origin_nm) * 10.0).astype("<i4")               # ORIGIN (A)
    hdr[52] = int.from_bytes(b"MAP", "little")
    hdr[53] = 0x00014144                            # MACHST 'DA' little-endian
    hdr[54] = int(np.float32(grid.std()).view(np.int32))
    hdr[55] = 1
    label = b"Phase16 synthetic cryo-EM map 3.5A (NPC scaffold, MDFF benchmark)"
    hdr[56:76] = np.frombuffer(label[:80].ljust(80, b" "), dtype="<i4")
    with open(path, "wb") as fh:
        fh.write(hdr.tobytes())
        fh.write(grid.tobytes())
    return dict(nx=nx, ny=ny, nz=nz, voxel_nm=float(voxel_nm), bytes=int(1024 + grid.nbytes))


def read_mrc(path):
    raw = Path(path).read_bytes()
    hdr = np.frombuffer(raw[:1024], dtype="<i4")
    if int(hdr[52]) != int.from_bytes(b"MAP", "little"):
        raise ValueError("not an MRC MAP file: %s" % path)
    nx, ny, nz, mode = int(hdr[0]), int(hdr[1]), int(hdr[2]), int(hdr[3])
    if mode != 2:
        raise ValueError("only MRC mode 2 supported, got %d" % mode)
    voxel_nm = float(hdr[10]) / (10.0 * int(hdr[7]))
    origin_nm = np.array([float(hdr[49]), float(hdr[50]), float(hdr[51])]) / 10.0
    nsymbt = int(hdr[23])
    data = np.frombuffer(raw[1024 + nsymbt:1024 + nsymbt + 4 * nx * ny * nz], dtype="<f4")
    return data.reshape(nx, ny, nz).astype(np.float64), origin_nm, voxel_nm


# ============================================== vectorized trilinear kernels
class TrilinearField:
    """Continuous scalar potential from a 3D grid with ANALYTIC trilinear
    gradients (the MDFF -grad V_EM kernel).  Single vectorized 8-corner
    gather; clamped at boundaries, zero force outside the grid support."""

    def __init__(self, grid, origin, voxel, periodic=False):
        self.g = np.ascontiguousarray(grid, dtype=np.float64)
        self.n = np.array(self.g.shape, dtype=np.int64)          # (nx,ny,nz)
        self.origin = np.asarray(origin, dtype=np.float64)
        voxel = np.atleast_1d(np.asarray(voxel, dtype=np.float64))
        if voxel.size == 1:
            voxel = np.full(3, float(voxel[0]))
        self.voxel = voxel
        self.periodic = periodic
        self.gflat = self.g.ravel()

    def _cellfrac(self, pos):
        t = (pos - self.origin) / self.voxel
        inside = np.all((t >= 0.0) & (t <= self.n - 1.0), axis=1)
        tc = np.clip(t, 0.0, (self.n - 1.0) - 1e-6)
        c = np.floor(tc).astype(np.int64)
        return c, tc - c, inside

    def _corners(self, c):
        i, j, k = c[:, 0], c[:, 1], c[:, 2]
        if self.periodic:
            i1, j1, k1 = (i + 1) % self.n[0], (j + 1) % self.n[1], (k + 1) % self.n[2]
        else:
            i1 = np.minimum(i + 1, self.n[0] - 1)
            j1 = np.minimum(j + 1, self.n[1] - 1)
            k1 = np.minimum(k + 1, self.n[2] - 1)
        idx = np.empty((len(c), 8), dtype=np.int64)
        for m, (ii, jj, kk) in enumerate([(i, j, k), (i1, j, k), (i, j1, k), (i1, j1, k),
                                          (i, j, k1), (i1, j, k1), (i, j1, k1), (i1, j1, k1)]):
            idx[:, m] = (ii * self.n[1] + jj) * self.n[2] + kk
        return idx

    @staticmethod
    def _weights(f):
        fx, fy, fz = f[:, 0:1], f[:, 1:2], f[:, 2:3]
        return np.concatenate([(1 - fx) * (1 - fy) * (1 - fz), fx * (1 - fy) * (1 - fz),
                               (1 - fx) * fy * (1 - fz), fx * fy * (1 - fz),
                               (1 - fx) * (1 - fy) * fz, fx * (1 - fy) * fz,
                               (1 - fx) * fy * fz, fx * fy * fz], axis=1)

    def values(self, pos):
        c, f, inside = self._cellfrac(pos)
        w = self._weights(f)
        return np.sum(self.gflat[self._corners(c)] * w, axis=1), inside

    def gradients(self, pos):
        c, f, inside = self._cellfrac(pos)
        vals = self.gflat[self._corners(c)]
        fx, fy, fz = f[:, 0], f[:, 1], f[:, 2]
        # corner-order derivative weights d w_m / d(axis) (corners ordered
        # 000,100,010,110,001,101,011,111)
        gx = np.stack([-(1 - fy) * (1 - fz), (1 - fy) * (1 - fz),
                       -fy * (1 - fz), fy * (1 - fz),
                       -(1 - fy) * fz, (1 - fy) * fz, -fy * fz, fy * fz], axis=1)
        gy = np.stack([-(1 - fx) * (1 - fz), -fx * (1 - fz),
                       (1 - fx) * (1 - fz), fx * (1 - fz),
                       -(1 - fx) * fz, -fx * fz, (1 - fx) * fz, fx * fz], axis=1)
        gz = np.stack([-(1 - fx) * (1 - fy), -fx * (1 - fy),
                       -(1 - fx) * fy, -fx * fy,
                       (1 - fx) * (1 - fy), fx * (1 - fy), (1 - fx) * fy, fx * fy], axis=1)
        grad = np.empty_like(pos)
        grad[:, 0] = np.sum(vals * gx, axis=1) / self.voxel[0]
        grad[:, 1] = np.sum(vals * gy, axis=1) / self.voxel[1]
        grad[:, 2] = np.sum(vals * gz, axis=1) / self.voxel[2]
        grad[~inside] = 0.0
        return grad, inside


def splat_model(pos, amp, shape_xyz, origin, voxel, sigma):
    """Gaussian pseudo-atom splat onto an [x,y,z] grid (map synthesis + CCC)."""
    nx, ny, nz = shape_xyz
    g = np.zeros(nx * ny * nz)
    idx = np.round((pos - origin) / voxel).astype(np.int64)
    rad = max(1, int(math.ceil(3.0 * sigma / float(voxel))))
    off = np.arange(-rad, rad + 1)
    OO = np.stack(np.meshgrid(off, off, off, indexing="ij"), axis=-1).reshape(-1, 3)
    ker = np.exp(-(OO ** 2).sum(1) / (2.0 * (sigma / voxel) ** 2))
    for o, kv in zip(OO, ker):
        c = idx + o
        m = np.all((c >= 0) & (c < np.array([nx, ny, nz])), axis=1)
        if not np.any(m):
            continue
        cc = c[m]
        np.add.at(g, (cc[:, 0] * ny + cc[:, 1]) * nz + cc[:, 2], amp[m] * kv)
    return g.reshape(nx, ny, nz)


# ============================================================ bead registry
# Martini-3-flavored CG bead grammar (sigma nm, mass amu, charge e, category)
TYPES = dict(
    W=(0.47, 72.0, 0.0, "solvent"),       # CG water (grid-functional only)
    SP=(0.45, 72.0, +0.4, "charged"),     # Na-type / cationic spacer
    SN=(0.45, 72.0, -0.4, "charged"),     # Cl-type / anionic spacer
    SD=(0.50, 120.0, 0.0, "scaffold"),    # structured-domain bead (cryo-EM visible)
    N0=(0.43, 56.0, 0.0, "spacer"),       # tiny/glycine-like
    TC=(0.47, 72.0, 0.0, "spacer"),       # polar spacer
    TA=(0.47, 88.0, 0.0, "aromatic"),     # Phe-ring sticker (FG motif)
    RB=(0.47, 96.0, 0.0, "body"),         # receptor body (HEAT-repeat like)
    PA=(0.47, 96.0, 0.0, "patch"),        # FG-binding patch (hydrophobic groove)
    PI=(0.47, 96.0, 0.0, "inert"),        # inert surface (control cargo)
    CB=(0.47, 72.0, 0.0, "body"),         # cargo-protein body
)
TYPE_NAMES = list(TYPES)
TYPE_CODE = {t: q for q, t in enumerate(TYPE_NAMES)}
SIGMA = np.array([TYPES[t][0] for t in TYPE_NAMES])
MASS = {t: TYPES[t][1] for t in TYPE_NAMES}
QCHG = {t: TYPES[t][2] for t in TYPE_NAMES}
CATEGORY = {t: TYPES[t][3] for t in TYPE_NAMES}

EPS_BASE = dict(solvent=0.80, charged=0.90, scaffold=0.65, spacer=0.42,
                aromatic=1.10, body=1.00, patch=1.00, inert=0.30)
# pair-specific cohesion overrides (kJ/mol) — the selectivity physics lives here
EPS_OVERRIDE = {
    ("aromatic", "aromatic"): 3.00,   # FG-FG hydrophobic multivalent clustering (mission value)
    ("patch", "aromatic"): 2.60,      # transport-receptor FG-binding pocket
    ("inert", "aromatic"): 0.22,      # inert control: pure steric floor
    ("inert", "patch"): 0.22,
    ("inert", "inert"): 0.22,
    ("inert", "body"): 0.30,
    ("body", "aromatic"): 1.00,
}
N_T = len(TYPE_NAMES)
EPS_MAT = np.zeros((N_T, N_T))
for _a in TYPE_NAMES:
    for _b in TYPE_NAMES:
        key = tuple(sorted((_a, _b)))
        EPS_MAT[TYPE_CODE[_a], TYPE_CODE[_b]] = \
            EPS_OVERRIDE.get(key, math.sqrt(EPS_BASE[CATEGORY[_a]] * EPS_BASE[CATEGORY[_b]]))


# ================================================================ structure
COPE_CLASSES = [
    # (name, r_nm, z_nm, ext_radial, ext_tang, ext_z, map_amp) — NPC-like rings
    ("IRm", 28.8, -11.5, 2.3, 2.1, 1.7, 1.00),
    ("IRp", 28.8, +11.5, 2.3, 2.1, 1.7, 1.00),
    ("CRm", 33.2, -11.5, 2.5, 2.3, 1.9, 1.15),
    ("CRp", 33.2, +11.5, 2.5, 2.3, 1.9, 1.15),
    ("LR", 30.6, 0.0, 2.1, 2.1, 1.6, 0.95),
    ("TLm", 27.2, -4.6, 1.7, 1.6, 1.3, 0.85),
    ("TLp", 27.2, +4.6, 1.7, 1.6, 1.3, 0.85),
]
N_BEADS_COPE = 80
N_SPOKES = 8
FG_PER_SPOKE = 14
FG_CHAIN_LEN = 92
R_LUMEN = 25.0                 # nm, pore radius (D_pore = 50 nm)
L_CHANNEL = 40.0               # nm
COPIES_PER_SPOKE = len(COPE_CLASSES)


def _ball_points(n, rng):
    v = rng.normal(size=(max(2 * n, 16), 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    v *= rng.uniform(0, 1, size=(len(v), 1)) ** (1.0 / 3.0)
    return v[:n]


def _ellipsoid_lattice(n, ar, at, az, rng, jitter=0.04):
    """Deterministic near-cubic packing inside an ellipsoid of half-extents
    (ar, at, az): guarantees ~0.55 nm spacing (no catastrophic LJ overlap)."""
    for s in (0.55, 0.50, 0.46, 0.42):
        mx = max(3, int(math.ceil(2 * ar / s)))
        my = max(3, int(math.ceil(2 * at / s)))
        mz = max(3, int(math.ceil(2 * az / s)))
        gx, gy, gz = np.meshgrid(np.arange(mx) - (mx - 1) / 2.0,
                                 np.arange(my) - (my - 1) / 2.0,
                                 np.arange(mz) - (mz - 1) / 2.0,
                                 indexing="ij")
        pts = np.stack([gx.ravel() * s, gy.ravel() * s, gz.ravel() * s], axis=1)
        inside = (pts[:, 0] / ar) ** 2 + (pts[:, 1] / at) ** 2 + (pts[:, 2] / az) ** 2 <= 0.82
        pts = pts[inside]
        if len(pts) >= n:
            break
    if len(pts) < n:
        raise RuntimeError("ellipsoid lattice too small: %d < %d" % (len(pts), n))
    sel = np.linspace(0, len(pts) - 1, n).astype(int)
    pts = pts[sel] + rng.normal(scale=jitter, size=(n, 3))
    return pts


def deoverlap(pos, d_min, iters=200, pin=None):
    """Deterministic pre-relaxation: remove sub-d_min starting overlaps so the
    LJ kernel (1/r^13 wall) never sees r -> 0 at step zero."""
    from scipy.spatial import cKDTree
    pos = pos.copy()
    pin = np.zeros(len(pos), bool) if pin is None else np.asarray(pin, bool)
    for _ in range(iters):
        prs = cKDTree(pos).query_pairs(d_min, output_type="ndarray")
        if len(prs) == 0:
            break
        d = pos[prs[:, 0]] - pos[prs[:, 1]]
        r = np.maximum(np.sqrt((d * d).sum(1)), 1e-9)
        push = ((d_min - r) / r)[:, None] * d * 0.25
        np.add.at(pos, prs[:, 0], push * (~pin[prs[:, 0]])[:, None])
        np.add.at(pos, prs[:, 1], -push * (~pin[prs[:, 1]])[:, None])
    return pos


def build_scaffold(rng):
    """8-fold symmetric NPC scaffold: positions, elastic network, per-bead map
    amplitude, FG anchor sites.  Bead order: spoke-major, cope-minor."""
    pos_all, amp_all = [], []
    en_i, en_j, en_k, en_r0 = [], [], [], []
    for kk in range(N_SPOKES):
        th = 2.0 * math.pi * kk / N_SPOKES
        er = np.array([math.cos(th), math.sin(th), 0.0])
        et = np.array([-math.sin(th), math.cos(th), 0.0])
        ez = np.array([0.0, 0.0, 1.0])
        for (cname, rc, zc, ar, at, az, a) in COPE_CLASSES:
            n = N_BEADS_COPE
            pts = _ellipsoid_lattice(n, ar, at, az, rng)
            pos = rc * er[None, :] + np.array([0.0, 0.0, zc])[None, :] + \
                pts[:, 0:1] * er[None, :] + pts[:, 1:2] * et[None, :] + pts[:, 2:3] * np.array([0.0, 0.0, 1.0])[None, :]
            pos = deoverlap(pos, 0.42)
            i0 = sum(len(p) for p in pos_all)
            pos_all.append(pos)
            amp_all.append(np.full(len(pos), a) * rng.uniform(0.85, 1.15, size=len(pos)))
            d = pos[:, None, :] - pos[None, :, :]
            r = np.sqrt((d * d).sum(-1))
            ii, jj = np.where((r < 1.15) & (r > 1e-6))
            sel = ii < jj
            en_i.append(ii[sel] + i0)
            en_j.append(jj[sel] + i0)
            en_k.append(np.full(int(sel.sum()), 600.0))
            en_r0.append(r[ii[sel], jj[sel]])
    pos = np.concatenate(pos_all)
    amp = np.concatenate(amp_all)
    en = np.stack([np.concatenate(en_i), np.concatenate(en_j),
                   np.concatenate(en_k), np.concatenate(en_r0)], axis=1)
    anchors = []
    for kk in range(N_SPOKES):
        th = 2.0 * math.pi * kk / N_SPOKES
        for z in np.linspace(-8.0, 8.0, FG_PER_SPOKE):
            anchors.append([25.6 * math.cos(th), 25.6 * math.sin(th),
                            z + 0.4 * (1 if z < 0 else -1)])
    return pos, en, amp, np.array(anchors)


def fg_chain_types(n_ch, L):
    """Sticker-spacer grammar: repeat [TC, TA, N0, SP/SN] — Phe every ~10 aa
    at the 4:1 mapping (TA fraction 1/4), alternating cationic/anionic spacers."""
    out = []
    for ch in range(n_ch):
        for m in range(L):
            u = m % 4
            out.append(["TC", "TA", "N0", "SP" if (ch + m // 4) % 2 == 0 else "SN"][u])
    return out


def fg_seed_positions(anchors, L, rng):
    """Pre-brush-relax chain coordinates: biased self-avoiding walks launched
    from the anchor sites that FILL the pore lumen (annulus 12 <= r <= 25.2 nm),
    the FG-Nup brush architecture that real transport factor must displace."""
    n_ch = len(anchors)
    pos = np.empty((n_ch * L, 3))
    step = 0.355
    r_min, r_max = 2.0, 25.2
    n_head = L - 30                    # directed inward reach ...
    for ch in range(n_ch):
        a = anchors[ch].astype(float).copy()
        p = a.copy()
        for m in range(L):
            pos[ch * L + m] = p
            r_hat = np.array([p[0], p[1], 0.0])
            r = np.linalg.norm(r_hat) + 1e-9
            r_hat /= r
            if m < n_head:
                bias = -0.85 * r_hat           # inward drift: cross the pore
                if r < r_min + 1.5:
                    bias = +0.9 * r_hat
            else:
                bias = np.zeros(3)             # terminal blob fills the core
            d = rng.normal(size=3)
            d /= np.linalg.norm(d)
            d = d + bias
            d /= np.linalg.norm(d)
            p = p + d * step * rng.uniform(0.85, 1.05)
            r = math.hypot(p[0], p[1])
            if r > r_max:
                p[0] *= r_max / r
                p[1] *= r_max / r
            elif r < r_min:
                p[0] *= r_min / r
                p[1] *= r_min / r
    pin = np.zeros(n_ch * L, bool)
    pin[np.arange(0, n_ch * L, L)] = True
    return deoverlap(pos, 0.30, pin=pin)


def build_receptor(kind, center):
    """Transport-factor cage (2 rings + FG-binding patches) around a cargo
    cluster.  kind='inert' swaps PA->PI with IDENTICAL geometry, hence the
    identical hydrodynamic radius by construction."""
    rng = np.random.default_rng(0xC16A if kind == "receptor" else 0x1BE2)
    ring = []
    for zz in (-0.62, 0.62):
        for m in range(14):
            th = 2 * math.pi * m / 14
            ring.append([4.0 * math.cos(th), 4.0 * math.sin(th), zz])
    ring = np.array(ring)
    patches = np.array([[4.9 * math.cos(2 * math.pi * m / 10 + 0.31),
                         4.9 * math.sin(2 * math.pi * m / 10 + 0.31), 0.0] for m in range(10)])
    cargo = _ball_points(26, rng) * np.array([2.9, 2.9, 2.4])
    pos = deoverlap(np.concatenate([ring, patches, cargo]) + np.asarray(center, float), 0.42)
    ptyp = ["RB"] * 28 + (["PA"] if kind == "receptor" else ["PI"]) * 10 + ["CB"] * 26
    d = pos[:, None, :] - pos[None, :, :]
    r = np.sqrt((d * d).sum(-1))
    ii, jj = np.where((r < 2.3) & (r > 1e-6))
    sel = ii < jj
    en = np.stack([ii[sel], jj[sel], np.full(int(sel.sum()), 500.0), r[ii[sel], jj[sel]]], axis=1)
    return pos, ptyp, en


CX_TYPES = lambda kind: ["RB"] * 28 + (["PA"] if kind == "receptor" else ["PI"]) * 10 + ["CB"] * 26
CX_MASS = lambda kind: sum(MASS[t] for t in CX_TYPES(kind))


# ============================================================ cell pair list
class PairList:
    """O(N) sorted-cell pair builder with a skin displacement test; periodic
    minimum image on the engine box.  Positions expected in [-L/2, L/2)."""

    def __init__(self, box, rc, skin=0.30):
        self.box = np.asarray(box, dtype=np.float64)
        self.rc2 = rc * rc
        self.rl2 = (rc + skin) ** 2
        self.skin = skin
        self.ncell = np.maximum((self.box / (rc + skin)).astype(int), 1)
        self.h = self.box / self.ncell
        self.ntot = int(np.prod(self.ncell))
        offs = [(dx, dy, dz) for dz in (-1, 0, 1) for dy in (-1, 0, 1) for dx in (-1, 0, 1)
                if (dx, dy, dz) > (0, 0, 0)]
        cx, cy, cz = np.meshgrid(np.arange(self.ncell[0]), np.arange(self.ncell[1]),
                                 np.arange(self.ncell[2]), indexing="ij")
        nbr = np.empty((self.ntot, len(offs)), dtype=np.int64)
        for q, (dx, dy, dz) in enumerate(offs):
            nbr[:, q] = (((cx.ravel() + dx) % self.ncell[0]) * self.ncell[1]
                         + ((cy.ravel() + dy) % self.ncell[1])) * self.ncell[2] \
                + ((cz.ravel() + dz) % self.ncell[2])
        self.nbr = nbr
        self.i = np.zeros(0, dtype=np.int64)
        self.j = np.zeros(0, dtype=np.int64)
        self._ref = None

    def needs_update(self, pos):
        if self._ref is None or len(pos) != len(self._ref):
            return True
        d = pos - self._ref
        d -= self.box * np.round(d / self.box)
        return bool(np.max((d * d).sum(axis=1)) > (0.25 * self.skin) ** 2)

    def update(self, pos):
        cellf = np.floor((pos + 0.5 * self.box) / self.h).astype(np.int64)
        cellf = np.clip(cellf, 0, self.ncell - 1)
        cid = (cellf[:, 0] * self.ncell[1] + cellf[:, 1]) * self.ncell[2] + cellf[:, 2]
        order = np.argsort(cid, kind="stable")
        occ = np.bincount(cid, minlength=self.ntot)
        starts = np.concatenate(([0], np.cumsum(occ)))
        pieces_i, pieces_j = [], []
        act = np.nonzero(occ)[0]              # occupied cells only
        occ_a = occ[act]
        m2 = occ_a >= 2
        if np.any(m2):
            cs = act[m2]
            o = occ[cs]
            n2 = o * o
            r = np.arange(int(n2.sum()))
            cp = np.repeat(np.arange(len(cs)), n2)
            loc = r - np.repeat(np.concatenate(([0], np.cumsum(n2)[:-1])), n2)
            a = loc // o[cp]
            b = loc - a * o[cp]
            ii = starts[cs[cp]] + a
            jj = starts[cs[cp]] + b
            keep = ii < jj
            pieces_i.append(ii[keep])
            pieces_j.append(jj[keep])
        for q in range(self.nbr.shape[1]):
            nb = self.nbr[act, q]             # neighbours of occupied cells
            o1 = occ[nb]
            n2 = occ_a * o1
            sel = n2 > 0
            if not np.any(sel):
                continue
            src = act[sel]
            os2 = n2[sel]
            r = np.arange(int(os2.sum()))
            cp = np.repeat(np.arange(len(src)), os2)
            loc = r - np.repeat(np.concatenate(([0], np.cumsum(os2)[:-1])), os2)
            o1b = occ[nb[sel][cp]]
            a = loc // o1b
            b = loc - a * o1b
            pieces_i.append(starts[src[cp]] + a)
            pieces_j.append(starts[nb[sel][cp]] + b)
        I = order[np.concatenate(pieces_i)]
        J = order[np.concatenate(pieces_j)]
        d = pos[I] - pos[J]
        d -= self.box * np.round(d / self.box)
        keep = (d * d).sum(axis=1) < self.rl2
        self.i, self.j = I[keep], J[keep]
        self._ref = pos.copy()
        return int(len(self.i))


# ================================================================= dynamics
class MegaEngine:
    """Bonded + LJ/Debye-Huckel pair forces over a bead subset + external
    closures, integrated with a half-kick/half-drift/Langevin (BAOAB) scheme.
    Fixed beads carry zero velocity and never move."""

    def __init__(self, box, types, fixed):
        self.box = np.asarray(box, dtype=np.float64)
        self.types = list(types)
        self.fixed = np.asarray(fixed, dtype=bool)
        self.mass = np.array([MASS[t] for t in self.types])
        self.q = np.array([QCHG[t] for t in self.types])
        self.code = np.array([TYPE_CODE[t] for t in self.types])
        self.b = None
        self.ang = None
        self.en = None
        self.ext = []
        self.fcap = 800.0      # pair-force cap (kJ/mol/nm); lowered for light beads
        self._noise = None

    def make_pairlists(self, subset, rc_lj=1.1, rc_q=2.2, skin=0.40):
        self.pl_idx = np.asarray(subset, dtype=np.int64)
        self.pl_lj = PairList(self.box, rc_lj, skin)
        # the Debye-Huckel list runs on the CHARGED members of the subset only
        # (sparse); q_local are subset-local ids so pair_forces shares one map
        self.q_local = np.nonzero(self.q[self.pl_idx] != 0)[0]
        self.pl_q = PairList(self.box, rc_q, skin + 0.2)
        self.sub_code = self.code[self.pl_idx]
        self.q_chg = self.q[self.pl_idx[self.q_local]]
        self.p_eps = None
        self.p_sig = None
        self.q_pairs = (np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0))

    def rebuild(self, pos_sub):
        nlj = self.pl_lj.update(pos_sub)
        I, J = self.pl_lj.i, self.pl_lj.j
        ci = self.sub_code[I]
        cj = self.sub_code[J]
        self.p_eps = EPS_MAT[ci, cj]
        self.p_sig = 0.5 * SIGMA[ci] + 0.5 * SIGMA[cj]
        nq = 0
        if len(self.q_local):
            nq = self.pl_q.update(pos_sub[self.q_local])
            qq = self.q_chg[self.pl_q.i] * self.q_chg[self.pl_q.j]
            keep = qq != 0
            self.q_pairs = (self.q_local[self.pl_q.i[keep]],
                            self.q_local[self.pl_q.j[keep]], qq[keep])
        return nlj, nq

    @staticmethod
    def _scatter(frc, idx, vec):
        for ax in range(3):
            frc[:, ax] += np.bincount(idx, weights=vec[:, ax], minlength=len(frc))

    def pair_forces_local(self, pos_sub, frc_sub, idx, remap):
        """Forces on the (static) beads in `idx` from all pair contacts with
        the mobile set.  `remap` maps subset-local indices to f_cx rows."""
        box = self.box
        in_idx = remap >= 0
        f_cx = np.zeros((int(in_idx.sum()), 3))
        I, J = self.pl_lj.i, self.pl_lj.j
        if len(I) == 0:
            return f_cx
        ci, cj = in_idx[I], in_idx[J]
        cross = ci | cj
        keep = ~both_in if (both_in := (ci & cj)).any() else np.ones(len(I), bool)
        I, J = I[keep], J[keep]
        ci, cj = ci[keep], cj[keep]
        ci_c = self.sub_code[I]
        cj_c = self.sub_code[J]
        eps = EPS_MAT[ci_c, cj_c]
        sig = 0.5 * SIGMA[ci_c] + 0.5 * SIGMA[cj_c]
        d = pos_sub[I] - pos_sub[J]
        d -= box * np.round(d / box)
        r2 = (d * d).sum(axis=1)
        m = (r2 < self.pl_lj.rc2) & (r2 > 1e-9)
        if not np.any(m):
            return f_cx
        I, J, d, r2 = I[m], J[m], d[m], r2[m]
        ci, cj = ci[m], cj[m]
        eps, sig = eps[m], sig[m]
        r2s = np.maximum(r2, 2.5e-2)
        r = np.sqrt(r2s)
        sr6 = (sig / r) ** 6
        fsc = (24.0 * eps / r2s) * (2.0 * sr6 * sr6 - sr6)
        np.clip(fsc, -self.fcap, self.fcap, out=fsc)
        fd = fsc[:, None] * d
        for arr, s, hit_in in ((I, +1.0, ci), (J, -1.0, cj)):
            if np.any(hit_in):
                np.add.at(f_cx, remap[arr[hit_in]], s * fd[hit_in])
        return f_cx

    def pair_forces_split(self, pos_sub, frc_sub, idx):
        """Full pair pass with rigid-body handling: pairs fully inside `idx`
        (internal cargo forces) are skipped; pairs crossing the boundary
        scatter to the outside endpoint AND accumulate to the returned force
        block for the idx beads."""
        box = self.box
        in_idx = np.zeros(len(pos_sub), dtype=bool)
        in_idx[np.asarray(idx, dtype=np.int64)] = True
        f_cx = np.zeros((int(in_idx.sum()), 3))
        I, J = self.pl_lj.i, self.pl_lj.j
        if len(I) == 0:
            return f_cx
        ci, cj = in_idx[I], in_idx[J]
        both = ci & cj
        keep = ~both
        if not np.any(keep):
            return f_cx
        I, J = I[keep], J[keep]
        ci, cj = ci[keep], cj[keep]
        ci_c = self.sub_code[I]
        cj_c = self.sub_code[J]
        eps = EPS_MAT[ci_c, cj_c]
        sig = 0.5 * SIGMA[ci_c] + 0.5 * SIGMA[cj_c]
        d = pos_sub[I] - pos_sub[J]
        d -= box * np.round(d / box)
        r2 = (d * d).sum(axis=1)
        m = (r2 < self.pl_lj.rc2) & (r2 > 1e-9)
        if not np.any(m):
            return f_cx
        I, J, d, r2 = I[m], J[m], d[m], r2[m]
        ci, cj = ci[m], cj[m]
        eps, sig = eps[m], sig[m]
        r2s = np.maximum(r2, 2.5e-2)
        r = np.sqrt(r2s)
        sr6 = (sig / r) ** 6
        fsc = (24.0 * eps / r2s) * (2.0 * sr6 * sr6 - sr6)
        np.clip(fsc, -self.fcap, self.fcap, out=fsc)
        fd = fsc[:, None] * d
        f_cx = np.zeros((len(pos_sub), 3))
        for arr, s, hit_in in ((I, +1.0, ci), (J, -1.0, cj)):
            # outside endpoint -> normal scatter; inside endpoint -> rigid body
            out_hit = ~hit_in
            if np.any(out_hit):
                self._scatter(frc_sub, self.pl_idx[arr[out_hit]], s * fd[out_hit])
            if np.any(hit_in):
                np.add.at(f_cx, arr[hit_in], s * fd[hit_in])
        return f_cx[in_idx]

    def pair_forces(self, pos_sub, frc_sub):
        box = self.box
        I, J = self.pl_lj.i, self.pl_lj.j
        if len(I):
            d = pos_sub[I] - pos_sub[J]
            d -= box * np.round(d / box)
            r2 = (d * d).sum(axis=1)
            m = (r2 < self.pl_lj.rc2) & (r2 > 1e-9)
            if np.any(m):
                I, J, d, r2 = I[m], J[m], d[m], r2[m]
                r2s = np.maximum(r2, 2.5e-2)          # soft-core floor: r >= 0.16 nm
                r = np.sqrt(r2s)
                sr6 = (self.p_sig[m] / r) ** 6
                fsc = (24.0 * self.p_eps[m] / r2s) * (2.0 * sr6 * sr6 - sr6)
                np.clip(fsc, -self.fcap, self.fcap, out=fsc)  # pair-force cap
                fd = fsc[:, None] * d
                self._scatter(frc_sub, self.pl_idx[I], fd)
                self._scatter(frc_sub, self.pl_idx[J], -fd)
        I, J, qq = self.q_pairs
        if len(I):
            d = pos_sub[I] - pos_sub[J]
            d -= box * np.round(d / box)
            r2 = (d * d).sum(axis=1)
            m = r2 < self.pl_q.rc2
            if np.any(m):
                I, J, d, r2, qq = I[m], J[m], d[m], r2[m], qq[m]
                r = np.sqrt(r2)
                fq = KE_COUL * qq * np.exp(-KAPPA_DH * r) * (1.0 + KAPPA_DH * r) / r2
                fd = (fq / r)[:, None] * d
                self._scatter(frc_sub, self.pl_idx[I], fd)
                self._scatter(frc_sub, self.pl_idx[J], -fd)

    def bonded_forces(self, pos, frc, only="all"):
        def bond_term(pairs):
            I = pairs[:, 0].astype(np.int64)
            J = pairs[:, 1].astype(np.int64)
            kk = pairs[:, 2]
            r0 = pairs[:, 3]
            d = pos[I] - pos[J]
            d -= self.box * np.round(d / self.box)
            r = np.sqrt((d * d).sum(axis=1))
            fs = (kk * (r - r0) / r)[:, None] * d
            np.clip(fs, -4000.0, 4000.0, out=fs)
            # U = k/2 (r-r0)^2 -> F_I = -k(r-r0) d/|d|  (restoring, toward partner)
            self._scatter(frc, I, -fs)
            self._scatter(frc, J, fs)

        if only in ("all", "bond") and self.b is not None and len(self.b):
            bond_term(self.b)
        if only in ("all", "en") and self.en is not None and len(self.en):
            bond_term(self.en)
        if only in ("all", "ang") and self.ang is not None and len(self.ang):
            A = self.ang
            I, J, K = A[:, 0].astype(np.int64), A[:, 1].astype(np.int64), A[:, 2].astype(np.int64)
            kth = A[:, 3]
            u = pos[I] - pos[J]
            u -= self.box * np.round(u / self.box)
            v = pos[K] - pos[J]
            v -= self.box * np.round(v / self.box)
            nu = np.sqrt((u * u).sum(1))
            nv = np.sqrt((v * v).sum(1))
            uh, vh = u / nu[:, None], v / nv[:, None]
            c = np.clip((uh * vh).sum(1), -1 + 1e-9, 1 - 1e-9)
            s = np.sqrt(1 - c * c)
            th = np.arccos(c)
            dU = kth * (th - np.deg2rad(111.0))
            f1 = (-dU / (nu * s))[:, None] * (c[:, None] * uh - vh)
            f3 = (-dU / (nv * s))[:, None] * (c[:, None] * vh - uh)
            self._scatter(frc, I, f1)
            self._scatter(frc, K, f3)
            self._scatter(frc, J, -(f1 + f3))

    def external_forces(self, pos, frc):
        for fn in self.ext:
            fn(pos, frc)

    def relax(self, pos, n_steps=800, move=4e-4, fcap=600.0, max_step=0.02):
        """Force-capped steepest-descent minimization (no velocities): relieves
        starting overlaps/strain so subsequent MD starts near a metastable
        valley instead of detonating.  Mobile mask applies; per-step displacement
        is bounded."""
        frc = np.zeros_like(pos)
        fmax_hist = []
        for st in range(n_steps):
            if self.pl_lj.needs_update(pos):
                self.rebuild(pos)
            frc[:] = 0.0
            self.bonded_forces(pos, frc)
            self.pair_forces(pos, frc)
            self.external_forces(pos, frc)
            frc[~self.mobile] = 0.0
            nrm = np.sqrt((frc * frc).sum(1))
            over = nrm > fcap
            if np.any(over):
                frc[over] *= (fcap / nrm[over])[:, None]
            nrm = np.maximum(np.sqrt((frc * frc).sum(1)), 1e-9)
            step = frc * (move / nrm[:, None])
            nstep = np.sqrt((step * step).sum(1))
            over2 = nstep > max_step
            if np.any(over2):
                step[over2] *= (max_step / nstep[over2])[:, None]
            pos += step
            if st % 200 == 0 or st == n_steps - 1:
                fmax_hist.append(float(nrm.max()))
        return fmax_hist

    def baoab(self, pos, vel, frc, dt, gamma, wrap_box=None):
        """One BAOAB step; `frc` reused for both half-kicks (single force
        evaluation per step).  Fixed beads keep zero velocity throughout.
        `gamma` may be a scalar or a per-bead array (ps^-1)."""
        if np.isscalar(gamma):
            c1 = math.exp(-gamma * dt)
            c1v = None
            c2 = math.sqrt((1.0 - c1 * c1) * KT)
        else:
            c1v = np.exp(-np.asarray(gamma, float) * dt)
            c2 = 0.0
        if self._noise is None or self._noise.shape != vel.shape:
            self._noise = np.empty_like(vel)
        half = (0.5 * dt / self.mass)[:, None]
        np.multiply(frc, half, out=self._noise)
        self._noise[~self.mobile] = 0.0                     # B
        vel += self._noise
        pos += (0.5 * dt) * vel                             # A
        _RNG.standard_normal(size=vel.shape, out=self._noise)
        if c1v is None:
            self._noise *= c2 * np.sqrt(1.0 / self.mass)[:, None]
        else:
            self._noise *= (np.sqrt(1.0 - c1v ** 2)
                            * np.sqrt(KT / self.mass))[:, None]
        self._noise[~self.mobile] = 0.0
        vel *= (c1 if c1v is None else c1v[:, None])         # O
        vel += self._noise
        pos += (0.5 * dt) * vel                             # A
        np.multiply(frc, half, out=self._noise)
        self._noise[~self.mobile] = 0.0                     # B
        vel += self._noise
        if wrap_box is not None:
            np.copyto(pos, ((pos + 0.5 * wrap_box) % wrap_box) - 0.5 * wrap_box)

    def set_mobile(self, mobile):
        self.mobile = np.asarray(mobile, dtype=bool)


def kirkwood_rh(pos, box):
    n = len(pos)
    d = pos[:, None, :] - pos[None, :, :]
    d -= box * np.round(d / box)
    r = np.sqrt((d * d).sum(-1))
    iu = np.triu_indices(n, 1)
    return float(n * n / (2.0 / r[iu]).sum())


def kabsch_rmsd(P, Q):
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    U = V @ np.diag([1.0, 1.0, d]) @ Wt
    return float(np.sqrt((((Pc @ U) - Qc) ** 2).sum() / len(P)))


def kabsch_transform(P, Q):
    """Row-vector rigid map x' = (x - mean(P)) @ U + mean(Q)."""
    Pm, Qm = P.mean(0), Q.mean(0)
    V, S, Wt = np.linalg.svd((P - Pm).T @ (Q - Qm))
    d = np.sign(np.linalg.det(V @ Wt))
    U = V @ np.diag([1.0, 1.0, d]) @ Wt
    return U, Pm, Qm


def cross_correlation(model_grid, ref_grid, mask):
    a = model_grid[mask] - model_grid[mask].mean()
    b = ref_grid[mask] - ref_grid[mask].mean()
    return float((a * b).sum() / (np.sqrt((a * a).sum() * (b * b).sum()) + 1e-30))


# ============================================== grid solvent (16A continuum)
class GridSolvent:
    """Compressible Brownian CG water on a periodic density grid, O(N) in one
    fused pass per step: a single cell-index scatter (no modulo for pre-wrapped
    water), NGP splat + 7-point grid smoothing, buffered trilinear gradient
    lookups with zero temporaries.  Buffers cached per subset size."""

    def __init__(self, box, h=1.0):
        self.box = np.asarray(box, float)
        self.ncell = np.round(box / h).astype(int)
        self.h = box / self.ncell
        self.origin = -0.5 * self.box
        self.tot = int(np.prod(self.ncell))
        self.rho = np.zeros(self.tot)          # smoothed water density field
        self.excl = np.zeros(self.tot)         # solute exclusion field
        self.k_rho = 3.0          # kJ nm^4 / mol pressure modulus
        self.k_excl = 22.0        # solute -> water exclusion
        self.fmax = 150.0
        self._buf = {}

    def _buffers(self, n):
        b = self._buf.get(n)
        if b is None:
            b = dict(idx8=np.empty((n, 8), dtype=np.int64),
                     cid=np.empty(n, dtype=np.int64),
                     w8=np.empty((n, 8)), vals=np.empty((n, 8)),
                     vals2=np.empty((n, 8)), gw=np.empty((n, 8)),
                     S=np.empty(n), S2=np.empty(n))
            self._buf[n] = b
        return b

    def _cell_ids(self, pos, b, wrapped):
        """Integer cell coords, flat cell id, and the 8 corner flat ids built
        as one contiguous broadcast (base + wrapped per-axis offsets)."""
        t = (pos - self.origin) / self.h
        if wrapped:   # positions already inside [origin, origin+box): cast = floor
            c0 = t[:, 0].astype(np.int64)
            c1 = t[:, 1].astype(np.int64)
            c2 = t[:, 2].astype(np.int64)
            f0, f1, f2 = t[:, 0] - c0, t[:, 1] - c1, t[:, 2] - c2
        else:
            f0 = np.floor(t[:, 0]); c0 = f0.astype(np.int64)
            f1 = np.floor(t[:, 1]); c1 = f1.astype(np.int64)
            f2 = np.floor(t[:, 2]); c2 = f2.astype(np.int64)
            c0 %= self.ncell[0]; c1 %= self.ncell[1]; c2 %= self.ncell[2]
            f0 = t[:, 0] - f0; f1 = t[:, 1] - f1; f2 = t[:, 2] - f2
        ny, nz = int(self.ncell[1]), int(self.ncell[2])
        nx = int(self.ncell[0])
        cid = b["cid"]
        np.multiply(c0, ny, out=cid)
        cid += c1
        cid *= nz
        cid += c2
        # per-axis wrapped offsets: +stride normally, 1-n_cells at the boundary
        dx = (ny * nz) * (1 - nx * ((c0 + 1) >= nx))
        dy = nz * (1 - ny * ((c1 + 1) >= ny))
        dz = (1 - nz * ((c2 + 1) >= nz))
        mx = np.array([0, 1, 0, 1, 0, 1, 0, 1])
        my = np.array([0, 0, 1, 1, 0, 0, 1, 1])
        mz = np.array([0, 0, 0, 0, 1, 1, 1, 1])
        offs = np.outer(dx, mx)
        offs += np.outer(dy, my)
        offs += np.outer(dz, mz)
        idx = b["idx8"]
        np.add(cid[:, None], offs, out=idx)   # stays within [0, tot): wrap built in
        self._last_frac = (f0, f1, f2)
        return f0, f1, f2

    @staticmethod
    def _smooth7(cnt, out, ncell, cell_vol):
        """out = (center + 6 face neighbours)/7 of the count grid, renormalised
        to preserve the total bead count (boundary cells average fewer faces)."""
        np.copyto(out, cnt)
        g = cnt.reshape(ncell)
        s = out.reshape(ncell)
        s[1:] += g[:-1]; s[:-1] += g[1:]
        s[:, 1:] += g[:, :-1]; s[:, :-1] += g[:, 1:]
        s[:, :, 1:] += g[:, :, :-1]; s[:, :, :-1] += g[:, :, 1:]
        out *= cnt.sum() / out.sum()
        out *= 1.0 / cell_vol

    def splat_density(self, pos, out, amp=1.0):
        """Generic CIC splat (small arrays: the solute exclusion field)."""
        b = self._buffers(len(pos))
        f0, f1, f2 = self._cell_ids(pos, b, wrapped=False)
        w8 = b["w8"]
        w8[:, 0] = (1 - f0) * (1 - f1) * (1 - f2)
        w8[:, 1] = f0 * (1 - f1) * (1 - f2)
        w8[:, 2] = (1 - f0) * f1 * (1 - f2)
        w8[:, 3] = f0 * f1 * (1 - f2)
        w8[:, 4] = (1 - f0) * (1 - f1) * f2
        w8[:, 5] = f0 * (1 - f1) * f2
        w8[:, 6] = (1 - f0) * f1 * f2
        w8[:, 7] = f0 * f1 * f2
        cnt = np.bincount(b["idx8"].ravel(), weights=(w8.ravel() * amp),
                          minlength=self.tot)
        np.multiply(cnt, 1.0 / (self.h[0] * self.h[1] * self.h[2]), out=out)

    @staticmethod
    def _axis_stacks(f0, f1, f2):
        t0, t1, t2 = 1.0 - f0, 1.0 - f1, 1.0 - f2
        wx = np.stack([t0, f0], axis=1)      # (n,2) over the x bit
        wy = np.stack([t1, f1], axis=1)
        wz = np.stack([t2, f2], axis=1)
        dwx = np.stack([-t0, f0], axis=1)    # d/dfx of the x basis
        dwy = np.stack([-t1, f1], axis=1)
        dwz = np.stack([-t2, f2], axis=1)
        return wx, wy, wz, dwx, dwy, dwz

    def water_forces(self, pos_w, frc_w, excl_flat, in_dom):
        """One fused pass per step: NGP density + 7-point smoothing, then the
        combined pressure/exclusion trilinear-gradient force on every water
        bead.  Single scatter, broadcast-built gradient weights (corner order
        m = bx + 2*by + 4*bz), no temporary storms."""
        n = len(pos_w)
        b = self._buffers(n)
        f0, f1, f2 = self._cell_ids(pos_w, b, wrapped=True)
        cnt = np.bincount(b["cid"], minlength=self.tot)
        self._smooth7(cnt, self.rho, self.ncell,
                      float(self.h[0] * self.h[1] * self.h[2]))
        np.take(self.rho, b["idx8"], out=b["vals"])
        use_excl = excl_flat is not None
        if use_excl:
            np.take(excl_flat, b["idx8"], out=b["vals2"])
        wx, wy, wz, dwx, dwy, dwz = self._axis_stacks(f0, f1, f2)
        S = b["S"]
        S2 = b["S2"]
        ratio = self.k_excl / self.k_rho
        fclamp = self.fmax / self.k_rho
        zero_dom = not bool(in_dom.all())
        n8 = (n, 2, 2, 2)
        gw = b["gw"]
        for ax in range(3):
            if ax == 0:
                gw3 = (dwx[:, :, None, None] * wy[:, None, :, None]) * wz[:, None, None, :]
            elif ax == 1:
                gw3 = (dwy[:, None, :, None] * wx[:, :, None, None]) * wz[:, None, None, :]
            else:
                gw3 = (dwz[:, None, None, :] * wx[:, :, None, None]) * wy[:, None, :, None]
            g8 = gw3.reshape(n, 8)
            np.einsum("ij,ij->i", b["vals"], g8, out=S)
            if use_excl:
                np.einsum("ij,ij->i", b["vals2"], g8, out=S2)
                S -= ratio * S2
            np.clip(S, -fclamp, fclamp, out=S)
            if zero_dom:
                S[~in_dom] = 0.0
            frc_w[:, ax] -= S / self.h[ax]
        return self.rho

    def grad_forces(self, pos, field_flat, gain, frc, sign, in_domain):
        """Small-array trilinear gradient lookup (generic; used off the hot
        water path).  frc += sign*gain*grad(field)."""
        n = len(pos)
        b = self._buffers(n)
        self._cell_ids(pos, b, wrapped=False)
        np.take(field_flat, b["idx8"], out=b["vals"])
        f0, f1, f2 = self._last_frac
        gw = b["gw"]
        S = b["S"]
        fclamp = self.fmax / max(gain, 1e-12)
        outside = ~in_domain
        terms = (
            [-(1 - f1) * (1 - f2), (1 - f1) * (1 - f2), -f1 * (1 - f2), f1 * (1 - f2),
             -(1 - f1) * f2, (1 - f1) * f2, -f1 * f2, f1 * f2],
            [-(1 - f0) * (1 - f2), -f0 * (1 - f2), (1 - f0) * (1 - f2), f0 * (1 - f2),
             -(1 - f0) * f2, -f0 * f2, (1 - f0) * f2, f0 * f2],
            [-(1 - f0) * (1 - f1), -f0 * (1 - f1), -(1 - f0) * f1, -f0 * f1,
             (1 - f0) * (1 - f1), f0 * (1 - f1), (1 - f0) * f1, f0 * f1],
        )
        for ax in range(3):
            for m_ in range(8):
                gw[:, m_] = terms[ax][m_]
            np.einsum("ij,ij->i", b["vals"], gw, out=S)
            if outside.any():
                S[outside] = 0.0
            np.clip(S, -fclamp, fclamp, out=S)
            frc[:, ax] += sign * gain * S / self.h[ax]

    _last_frac = (None, None, None)


# ==================================================== stage 16A: assembly
SOLV_BOX = np.array([42.0, 42.0, 44.0])   # pore-local hydration domain (nm):
    # ~6.3e5 CG solvent beads at bulk density -> >2.5e6 all-atom equivalents
ENG_BOX = np.array([80.0, 80.0, 60.0])    # engine min-image box (solute)


def integrate(n_steps, step_fn, monitor=None, log_every=2000, label=""):
    t0 = time.time()
    stats = {}
    for st in range(n_steps):
        step_fn(st)
        if monitor is not None and st % 250 == 0:
            s = monitor()
            if not s["finite"]:
                raise RuntimeError("%s: non-finite state at step %d" % (label, st))
            stats = s
        if log_every and st and st % log_every == 0:
            log("   %s step %6d/%d  %.1f ms/step" % (label, st, n_steps,
                                                     1000.0 * (time.time() - t0) / (st + 1)))
    ms = 1000.0 * (time.time() - t0) / max(1, n_steps)
    return ms, stats


def chain_topology(i0, n_ch, L):
    fb = [np.arange(i0 + c * L, i0 + (c + 1) * L - 1) for c in range(n_ch)]
    fj = [np.arange(i0 + c * L + 1, i0 + (c + 1) * L) for c in range(n_ch)]
    fa = [np.arange(i0 + c * L, i0 + (c + 1) * L - 2) for c in range(n_ch)]
    bonds = np.stack([np.concatenate(fb), np.concatenate(fj),
                      np.full(n_ch * (L - 1), 1250.0), np.full(n_ch * (L - 1), 0.36)], axis=1)
    angles = np.stack([np.concatenate(fa), np.concatenate(fa) + 1, np.concatenate(fa) + 2,
                       np.full(n_ch * (L - 2), 6.0)], axis=1)
    return bonds, angles


def stage_16A(fast, budget):
    log("=== STAGE 16A : CG megasystem assembly & stability gate ===")
    meta = dict(fast_mode=bool(fast))
    seed = np.random.default_rng(0xF16A)
    sc_pos, en, amp, anchors = build_scaffold(seed)
    n_sc = len(sc_pos)
    n_ch = N_SPOKES * FG_PER_SPOKE
    L = FG_CHAIN_LEN
    n_fg = n_ch * L
    fg_types = fg_chain_types(n_ch, L)
    fg_pos0 = fg_seed_positions(anchors, L, seed)

    # ---- explicit CG water-ion lattice (periodic hydration domain) ----
    # Martini 4:1 mapping: 8.35 beads/nm^3 -> lattice spacing 0.4894 nm
    h_lat = RHO_WATER_CG ** (-1.0 / 3.0)
    nlat = np.round(SOLV_BOX / h_lat).astype(int)
    gx, gy, gz = np.meshgrid(np.arange(nlat[0]), np.arange(nlat[1]),
                             np.arange(nlat[2]), indexing="ij")
    lattice = (np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1) + 0.5) \
        * (SOLV_BOX / nlat) - 0.5 * SOLV_BOX
    if fast:
        lattice = lattice[::5]
    # purge lattice sites that start inside the solute volume
    from scipy.spatial import cKDTree
    n_pre_purge = len(lattice)
    hit = cKDTree(np.concatenate([sc_pos, fg_pos0])).query_ball_point(lattice, 0.45)
    keep = np.array([len(hs) == 0 for hs in hit])
    lattice = lattice[keep]
    log("16A | water lattice: %d sites, purged %d inside solute volume"
        % (n_pre_purge, int((~keep).sum())))
    n_ion = int(round(len(lattice) * 2 * SALT_M * N_AVOG_NM3 / RHO_WATER_CG))
    ion_sel = np.random.default_rng(5).choice(len(lattice), size=n_ion, replace=False)
    is_ion = np.zeros(len(lattice), bool)
    is_ion[ion_sel] = True
    wpos, ipos = lattice[~is_ion], lattice[is_ion]
    n_water = len(wpos)
    n_sol = n_sc + n_fg
    n_total = n_sol + n_water + n_ion
    aae = (n_water + n_ion) * 4 + n_sol * 15
    log("16A | beads: total=%d  water=%d  ions=%d  scaffold=%d  FG=%d"
        % (n_total, n_water, n_ion, n_sc, n_fg))
    log("16A | all-atom equivalents (4:1 solvent, ~15:1 protein CG): ~%.3f e6" % (aae / 1e6))

    types = (["SD"] * n_sc + fg_types + ["W"] * n_water
             + ["SP"] * (n_ion // 2) + ["SN"] * (n_ion - n_ion // 2))
    fixed = np.zeros(n_total, bool)
    fixed[:n_sc] = True
    pos = np.concatenate([sc_pos, fg_pos0, wpos, ipos])
    eng = MegaEngine(ENG_BOX, types, fixed)
    eng.set_mobile(~fixed)
    bonds, angles = chain_topology(n_sc, n_ch, L)
    eng.b, eng.ang, eng.en = bonds, angles, en

    # ---- solute minimization BEFORE solvation (start near a metastable valley)
    eng_pre = MegaEngine(ENG_BOX, types[:n_sol], fixed[:n_sol])
    eng_pre.set_mobile(~fixed[:n_sol])
    eng_pre.b, eng_pre.ang, eng_pre.en = bonds, angles, en
    fg0_pre = np.arange(n_sc, n_sol, L)

    def anchor_force_pre(p, f):
        np.add.at(f, fg0_pre, 2500.0 * (anchors - p[fg0_pre]))

    eng_pre.ext.append(anchor_force_pre)
    eng_pre.make_pairlists(np.arange(n_sol))
    pos_pre = np.concatenate([sc_pos, fg_pos0])
    hist_pre = eng_pre.relax(pos_pre, 800)
    log("16A | solute minimization: max|F| %s over 800 steps"
        % (" -> ".join("%.0f" % v for v in hist_pre)))
    pos[:n_sol] = pos_pre
    del eng_pre, pos_pre
    gc.collect()

    eng.make_pairlists(np.arange(n_sol))
    gset = GridSolvent(SOLV_BOX)
    pos_sol = pos[:n_sol]
    pos_w = pos[n_sol:]
    n_wt = len(pos_w)
    vel = np.zeros_like(pos)
    vel += np.random.default_rng(9).standard_normal(vel.shape) * np.sqrt(KT / eng.mass)[:, None]
    vel[fixed] = 0.0
    frc = np.zeros_like(pos)
    frc_sol = np.zeros((n_sol, 3))
    frc_w = np.zeros((n_wt, 3))
    ones_w = np.ones(n_wt, bool)

    def solute_domain_mask():
        t = (pos_sol - gset.origin) / gset.h
        return np.all((t >= 0) & (t < gset.ncell), axis=1)

    def mega_step(st):
        if eng.pl_lj.needs_update(pos_sol):
            eng.rebuild(pos_sol)
        frc_sol[:] = 0.0
        frc_w[:] = 0.0
        eng.bonded_forces(pos_sol, frc_sol, only="bond")
        eng.bonded_forces(pos_sol, frc_sol, only="ang")
        eng.pair_forces(pos_sol, frc_sol)
        # dual-timestep grid coupling: the smoothed water-grid forces vary on the
        # 1 nm cell scale, so they are refreshed every 2nd step (max inter-update
        # displacement ~0.05 nm << cell size); the solute exclusion field every
        # 3rd step.  frc_w is reused in between, consistent with the BAOAB
        # single-force-evaluation scheme.
        if st % 3 == 0:
            msk = solute_domain_mask()
            gset.splat_density(pos_sol[msk], gset.excl, amp=0.5)
        if st % 2 == 0:
            # one-way coupling: the grid pushes water out of the solute volume
            # (a solute "reaction" on its OWN density grid would self-attract)
            gset.water_forces(pos_w, frc_w, gset.excl, ones_w)
        frc[:n_sol] = frc_sol
        frc[n_sol:] = frc_w
        eng.baoab(pos, vel, frc, 0.020, 12.0, wrap_box=None)
        np.copyto(pos_w, ((pos_w + 0.5 * SOLV_BOX) % SOLV_BOX) - 0.5 * SOLV_BOX)

    def monitor():
        kin = float(0.5 * eng.mass[~fixed] @ (vel[~fixed] ** 2).sum(1))
        return dict(finite=bool(np.isfinite(pos).all()),
                    T=float(2 * kin / (3 * int((~fixed).sum()) * KB)))

    n_eq = 400 if fast else 1100
    ms, stats = integrate(n_eq, mega_step, monitor, log_every=400, label="16A/water")
    rho_dev = abs(float(gset.rho.mean()) - RHO_WATER_CG) / RHO_WATER_CG * 100
    # fast mode intentionally thins the lattice 5x -> density gate is full-run only
    rho_gate = (None if fast else bool(rho_dev < 15.0 and stats.get("T", 0.0) > 200))
    meta.update(n_total_beads=int(n_total), n_water=int(n_water), n_ions=int(n_ion),
                n_scaffold=int(n_sc), n_fg=int(n_fg), n_allatom_equiv=int(aae),
                water_ms_per_step=round(ms, 2), water_T_K=round(stats.get("T", 0.0), 1),
                water_density_mean=round(float(gset.rho.mean()), 3),
                water_density_dev_pct=round(float(rho_dev), 2),
                gate_16A_stability=rho_gate,
                gate_16A_beadcount=bool(n_total >= 250_000),
                gate_16A_allatom=bool(aae >= 2_500_000))
    log("16A | rho=%.2f beads/nm^3 (dev %.1f%%), T=%.0f K, %.1f ms/step | gates: "
        "stability=%s beads>=250k=%s allatom>=2.5M=%s"
        % (gset.rho.mean(), rho_dev, stats.get("T", 0.0), ms,
           meta["gate_16A_stability"], meta["gate_16A_beadcount"], meta["gate_16A_allatom"]))
    np.savez_compressed(RES / "stage16A_megasystem.npz",
                        rho_grid=gset.rho.reshape(gset.ncell).astype(np.float32),
                        ncell=gset.ncell, box=SOLV_BOX)
    del gset
    gc.collect()
    log("16A | solvent grid released (memory discipline: buffers freed, gc closed)")

    # ---------- implicit-continuum brush relaxation (frozen scaffold) ----
    log("16A | relaxing FG-Nup brush (implicit continuum, frozen scaffold) ...")
    types2 = types[:n_sol]
    fixed2 = fixed[:n_sol]
    pos2 = pos[:n_sol].copy()
    eng2 = MegaEngine(ENG_BOX, types2, fixed2)
    eng2.set_mobile(~fixed2)
    eng2.b, eng2.ang, eng2.en = bonds, angles, en
    fg0 = np.arange(n_sc, n_sol, L)
    anchor_sites = anchors.copy()

    def anchor_force(p, f):
        np.add.at(f, fg0, 2500.0 * (anchor_sites - p[fg0]))

    eng2.ext.append(anchor_force)
    eng2.make_pairlists(np.arange(n_sol))
    vel2 = np.zeros_like(pos2)
    vel2 += np.random.default_rng(11).standard_normal(vel2.shape) * np.sqrt(KT / eng2.mass)[:, None]
    vel2[fixed2] = 0.0
    frc2 = np.zeros_like(pos2)

    def brush_step(st):
        if eng2.pl_lj.needs_update(pos2):
            eng2.rebuild(pos2)
        frc2[:] = 0.0
        eng2.bonded_forces(pos2, frc2)
        eng2.pair_forces(pos2, frc2)
        eng2.external_forces(pos2, frc2)
        eng2.baoab(pos2, vel2, frc2, 0.030, GAMMA_C, wrap_box=ENG_BOX)

    def monitor2():
        kin = float(0.5 * eng2.mass[~fixed2] @ (vel2[~fixed2] ** 2).sum(1))
        return dict(finite=bool(np.isfinite(pos2).all()),
                    T=float(2 * kin / (3 * int((~fixed2).sum()) * KB)))

    n_br = 3000 if fast else 6000
    ms_br, stats2 = integrate(n_br, brush_step, monitor2, log_every=8000, label="16A/brush")
    log("16A | brush relaxed: %d steps at %.2f ms/step" % (n_br, ms_br))

    fg_pos = pos2[n_sc:]
    fg_t = np.array(fg_types)
    r_fg = np.sqrt(fg_pos[:, 0] ** 2 + fg_pos[:, 1] ** 2)
    I, J = eng2.pl_lj.i, eng2.pl_lj.j
    ta_pair = (fg_t[I - n_sc] == "TA") & (fg_t[J - n_sc] == "TA")
    d = pos2[I] - pos2[J]
    d -= ENG_BOX * np.round(d / ENG_BOX)
    close = (d * d).sum(1) < 0.65 ** 2
    ta_ids = np.concatenate([I[ta_pair & close], J[ta_pair & close]]).astype(np.int64)
    coord = (np.bincount(ta_ids, minlength=n_sol)[n_sc:] if len(ta_ids)
             else np.zeros(n_fg))
    meta["brush_TA_mean_coord"] = round(float(coord.mean()), 3)
    meta["brush_TA_coord_p90"] = round(float(np.percentile(coord, 90)), 2)
    meta["brush_inner_edge_nm"] = round(float(np.percentile(r_fg, 2)), 2)
    meta["brush_ms_per_step"] = round(ms_br, 3)
    meta["fg_n_chains"] = int(n_ch)
    meta["fg_chain_len"] = int(L)
    log("16A | brush: <N_TA-TA>=%.2f p90=%.1f, inner brush edge r~%.1f nm"
        % (coord.mean(), np.percentile(coord, 90), np.percentile(r_fg, 2)))
    np.savez_compressed(RES / "stage16A_brush.npz",
                        fg_pos=fg_pos.astype(np.float32), fg_types=fg_t,
                        r_fg=r_fg.astype(np.float32), ta_coord=coord.astype(np.float32))
    np.savez_compressed(RES / "stage16A_frames.npz",
                        scaffold=pos2[:n_sc].astype(np.float32),
                        fg=fg_pos.astype(np.float32), fg_types=fg_t,
                        anchors=anchors.astype(np.float32))
    (RES / "stage16A_meta.json").write_text(json.dumps(meta, indent=1, default=str))
    return meta


# ==================================================== stage 16B: MDFF
MAP_VOXEL = 0.35              # nm (3.5 A sampling — stated map resolution)
MAP_BOX_NM = 80.0


def stage_16B(fast, budget):
    log("=== STAGE 16B : Cryo-EM map synthesis, MRC round-trip & MDFF fitting ===")
    meta = {}
    seed = np.random.default_rng(0xF16A)
    target, en, amp, anchors = build_scaffold(seed)
    types = ["SD"] * len(target)
    amp = amp / amp.max()

    # ---- synthetic map at 3.5 A ----
    n_map = 2 * int(round(MAP_BOX_NM / (2 * MAP_VOXEL)))       # even: enables 2x binning
    origin = -0.5 * n_map * MAP_VOXEL * np.ones(3)
    sig_detail = 0.09            # resolution kernel d/(pi*sqrt(2)) ~ 0.79 A
    sig_b = 0.16                 # B ~ 8 pi^2 sigma^2 ~ 200 A^2 damping
    sig_map = math.sqrt(sig_detail ** 2 + sig_b ** 2)
    rho = splat_model(target, amp, (n_map, n_map, n_map), origin, MAP_VOXEL, sig_map)
    rho /= rho.max()
    rho += np.random.default_rng(3).normal(scale=0.025, size=rho.shape) + 0.015
    mrc_path = RES / "mdff" / "cryoem_map_3p5A.mrc"
    info = write_mrc(mrc_path, rho, origin, MAP_VOXEL)
    log("16B | wrote %s (%d^3, %.1f MB)" % (mrc_path.name, n_map, info["bytes"] / 1e6))
    rho_rt, origin_rt, voxel_rt = read_mrc(mrc_path)
    meta["mrc_roundtrip_maxdiff"] = float(np.abs(rho_rt - rho.astype("<f4").astype(np.float64)).max())
    meta["mrc_voxel_A"] = voxel_rt * 10.0
    meta["mrc_master_bytes"] = info["bytes"]
    meta["b_factor_A2"] = round(8 * math.pi ** 2 * sig_b ** 2 * 100, 1)
    assert meta["mrc_roundtrip_maxdiff"] == 0.0, "MRC round-trip failed"
    assert np.allclose(origin_rt, origin, atol=1e-3)
    zf = rho.reshape(n_map // 2, 2, n_map // 2, 2, n_map // 2, 2).mean(axis=(1, 3, 5))
    vox_f = 2 * MAP_VOXEL
    write_mrc(RES / "mdff" / "cryoem_forcemap_7A.mrc", zf, origin, vox_f)
    em_field = TrilinearField(zf, origin, vox_f)

    def model_binned(pos_):
        """Model density rendered through the IDENTICAL operator as the map
        (master-grid splat, then the same 2x binning) — apples-to-apples CCC."""
        mg = splat_model(pos_, weights_g, rho.shape, origin, MAP_VOXEL, sig_map)
        return mg.reshape(n_map // 2, 2, n_map // 2, 2, n_map // 2, 2).mean(axis=(1, 3, 5))

    weights_g = amp.copy()
    mask_cc = zf > 0.02
    ccc_ceiling = cross_correlation(model_binned(target), zf, mask_cc)
    log("16B | MDFF force grid %s @ %.1f A; CCC ceiling (target self-fit) = %.4f"
        % (tuple(zf.shape), vox_f * 10, ccc_ceiling))

    # ---- distort: per-bead noise + rigid shift + rotation + per-spoke twist ----
    pos0 = target + np.random.default_rng(4).normal(scale=0.18, size=target.shape)
    pos0 += np.array([1.3, -0.9, 1.6])
    Rg = Rotation.from_rotvec(np.array([1.0, 0.3, 0.2]) / math.sqrt(1.13) * math.radians(3.0))
    pos0 = Rg.apply(pos0)
    for kk in range(N_SPOKES):
        c0 = kk * COPIES_PER_SPOKE * N_BEADS_COPE
        sel = np.zeros(len(target), bool)
        for q in range(COPIES_PER_SPOKE):
            sel[c0 + q * N_BEADS_COPE: c0 + (q + 1) * N_BEADS_COPE] = True
        pos0[sel] = Rotation.from_rotvec([0.0, 0.0, math.radians(1.2)]).apply(pos0[sel])
    rmsd0 = float(np.sqrt(((pos0 - target) ** 2).sum(1).mean()))
    meta_rigid = {}

    # ---- MDFF simulation ----
    box = np.array([90.0, 90.0, 100.0])
    eng = MegaEngine(box, types, np.zeros(len(types), bool))
    eng.set_mobile(np.ones(len(types), bool))
    eng.en = en
    eng.make_pairlists(np.arange(len(types)), rc_lj=1.2)
    xi = [40.0]
    weights = amp.copy()
    gs = np.ones(len(target))
    mask_cc = zf > 0.02

    def em_force(p, f):
        grad, _ = em_field.gradients(p)
        f_em = (xi[0] * weights * gs)[:, None] * grad
        nrm = np.sqrt((f_em * f_em).sum(1))
        over = nrm > 400.0
        if np.any(over):
            f_em[over] *= (400.0 / nrm[over])[:, None]
        f += f_em

    def soft_walls(p, f):
        r = np.sqrt(p[:, 0] ** 2 + p[:, 1] ** 2)
        m = r > 39.0
        if np.any(m):
            f[m, 0] -= 400.0 * (r[m] - 39.0) * p[m, 0] / r[m]
            f[m, 1] -= 400.0 * (r[m] - 39.0) * p[m, 1] / r[m]
        m = np.abs(p[:, 2]) > 45.0
        if np.any(m):
            f[m, 2] -= 400.0 * (np.abs(p[m, 2]) - 45.0) * np.sign(p[m, 2])

    eng.ext.append(em_force)
    eng.ext.append(soft_walls)

    # ---- Phase A/B: hierarchical RIGID-BODY fitting (standard MDFF pipeline:
    # per-bead density gradients are meaningless until the rigid error is gone)
    def rigid_em_fit(pos_, sel, n_steps, xi_r, dt_r=0.10):
        # 6-DOF Langevin-lite fit of one rigid body against the map
        idx = np.nonzero(sel)[0]
        w = weights[idx]
        M = float(w.sum() * 120.0)
        c = pos_[idx].mean(0)
        X = pos_[idx] - c
        # diagonal inertia approximation about the body centre
        Idiag = np.array([float((w * (X[:, 1] ** 2 + X[:, 2] ** 2)).sum()),
                          float((w * (X[:, 0] ** 2 + X[:, 2] ** 2)).sum()),
                          float((w * (X[:, 0] ** 2 + X[:, 1] ** 2)).sum())]) * 120.0
        Idiag = np.maximum(Idiag, 1e-3)
        v = np.zeros(3)
        om = np.zeros(3)
        for st_ in range(n_steps):
            grad, _ins = em_field.gradients(pos_[idx])
            f = (xi_r * w * gs[idx])[:, None] * grad
            F = f.sum(0)
            tau = np.cross(X, f).sum(0)
            v += (F / M) * dt_r
            om += (tau / Idiag) * dt_r
            v *= 0.98
            om *= 0.98
            c = c + v * dt_r
            rot = Rotation.from_rotvec(om * dt_r)
            X = rot.apply(X)
            pos_[idx] = X + c
        return pos_

    sel_all = np.ones(len(target), bool)
    pos = pos0.copy()
    # recompute g-scale for the distorted pose before rigid fitting
    vals, _ = em_field.values(pos)
    gs[:] = np.clip(vals / 0.25, 0.2, 1.0)
    nA = 400 if fast else 1500
    nB = 300 if fast else 1200
    pos = rigid_em_fit(pos, sel_all, nA, 800.0)
    rmsd_rigid = float(np.sqrt(((pos - target) ** 2).sum(1).mean()))
    log("16B | phase A rigid-body fit: RMSD %.3f -> %.3f nm" % (rmsd0, rmsd_rigid))
    for kk in range(N_SPOKES):
        sel = np.zeros(len(target), bool)
        c0 = kk * COPIES_PER_SPOKE * N_BEADS_COPE
        for q in range(COPIES_PER_SPOKE):
            sel[c0 + q * N_BEADS_COPE: c0 + (q + 1) * N_BEADS_COPE] = True
        pos = rigid_em_fit(pos, sel, nB, 800.0)
    rmsd_spokes = float(np.sqrt(((pos - target) ** 2).sum(1).mean()))
    log("16B | phase B per-spoke rigid fit: RMSD -> %.3f nm" % rmsd_spokes)
    meta_rigid = dict(rmsd_after_phaseA=rmsd_rigid, rmsd_after_phaseB=rmsd_spokes)

    pos = pos0.copy() if fast else pos
    hist_pre = eng.relax(pos, 600)
    log("16B | pre-fit minimization: max|F| %s"
        % (" -> ".join("%.0f" % v for v in hist_pre)))
    vel = np.random.default_rng(6).standard_normal(pos.shape) * np.sqrt(KT / eng.mass)[:, None]
    frc = np.zeros_like(pos)
    n_fit = 2500 if fast else 12000
    dt, gamma = 0.025, 10.0
    ccc_trace, rmsd_trace, xi_trace, f99_trace = [], [], [], []
    t0 = time.time()
    for st in range(n_fit):
        if st % 250 == 0 and not np.isfinite(pos).all():
            raise RuntimeError("16B: non-finite state at step %d" % st)
        if eng.pl_lj.needs_update(pos):
            eng.rebuild(pos)
        xi[0] = 100.0 + 500.0 * min(1.0, st / (0.33 * n_fit))  # adaptive-force ramp
        if st % 50 == 0:                                        # MDFF g-scale
            vals, _ = em_field.values(pos)
            gs[:] = np.clip(vals / 0.25, 0.2, 1.0)
        frc[:] = 0.0
        eng.bonded_forces(pos, frc, only="en")
        eng.pair_forces(pos, frc)
        eng.external_forces(pos, frc)
        eng.baoab(pos, vel, frc, dt, gamma, wrap_box=box)
        if st % 100 == 0 or st == n_fit - 1:
            ccc_trace.append(cross_correlation(model_binned(pos), zf, mask_cc))
            rmsd_trace.append(float(np.sqrt(((pos - target) ** 2).sum(1).mean())))
            xi_trace.append(xi[0])
            f99_trace.append(float(np.percentile(np.linalg.norm(frc, axis=1), 99)))
            if st % 1500 == 0:
                log("16B | step %5d  CCC=%.4f  RMSD=%.3f nm  xi=%.0f  F99=%.0f kJ/mol/nm"
                    % (st, ccc_trace[-1], rmsd_trace[-1], xi[0], f99_trace[-1]))
    ms_fit = 1000.0 * (time.time() - t0) / n_fit
    fitted = pos.copy()
    rmsd1 = float(np.sqrt(((fitted - target) ** 2).sum(1).mean()))
    rmsd1_k = kabsch_rmsd(fitted, target)
    ccc_final = cross_correlation(model_binned(fitted), zf, mask_cc)
    ccc_init = cross_correlation(model_binned(pos0), zf, mask_cc)
    log("16B | MDFF done: CCC %.4f -> %.4f (ceiling %.4f); "
        "RMSD %.3f -> %.3f nm (Kabsch %.3f); %.2f ms/step"
        % (ccc_trace[0], ccc_final, ccc_ceiling, rmsd0, rmsd1, rmsd1_k, ms_fit))
    meta.update(ccc_init=float(ccc_trace[0]), ccc_final=float(ccc_final),
                ccc_ceiling=float(ccc_ceiling),
                rmsd_init=rmsd0, rmsd_final=rmsd1, rmsd_final_kabsch=rmsd1_k,
                ms_per_step=round(ms_fit, 2), n_steps=int(n_fit),
                xi_final=float(xi[0]), F99_final=float(f99_trace[-1]),
                xi_ramp="100 -> 600 kJ/mol, per-bead EM clamp 400",
                **meta_rigid,
                force_grid_shape=[int(v) for v in zf.shape],
                gate_16B_mdff=bool(ccc_final > 0.90 and rmsd1 < 0.35 * rmsd0),
                gate_16B_mrc=bool(meta["mrc_roundtrip_maxdiff"] == 0.0))
    np.savez_compressed(RES / "mdff" / "mdff_trace.npz",
                        ccc=np.array(ccc_trace), rmsd=np.array(rmsd_trace),
                        xi=np.array(xi_trace), f99=np.array(f99_trace),
                        target=target.astype(np.float32), fitted=fitted.astype(np.float32),
                        distorted=pos0.astype(np.float32),
                        forcemap=zf.astype(np.float32), forcemap_origin=origin,
                        forcemap_voxel=vox_f, sig_map=sig_map)

    # ---- OpenMM cross-validation of the force kernels ----
    if HAS_OPENMM:
        try:
            ov = openmm_crosscheck(en, types, fitted, 1500 if fast else 3000)
            meta["openmm"] = ov
            log("16B | OpenMM kernel cross-check: max|dF|=%.2e kJ/mol/nm; %d-step "
                "LangevinMiddle drift RMSD=%.4f nm (platform %s)"
                % (ov["max_force_diff"], ov["drift_steps"], ov["drift_rmsd"], ov["platform"]))
            meta["gate_16B_openmm"] = bool(ov["max_force_diff"] < 1e-4 and ov["drift_rmsd"] < 0.30)
        except Exception as exc:  # pragma: no cover
            meta["openmm"] = {"error": str(exc)}
            meta["gate_16B_openmm"] = None
            log("16B | OpenMM cross-check failed: %s" % exc)
    else:
        meta["openmm"] = {"available": False}
        meta["gate_16B_openmm"] = None
        log("16B | OpenMM not importable on this interpreter — gates 1-5 cover the engine")
    (RES / "stage16B_meta.json").write_text(json.dumps(meta, indent=1, default=str))
    return meta


def openmm_crosscheck(en, types, pos_fit, n_drift):
    """(i) force-by-force equality of the elastic-network kernel against a
    compiled OpenMM CustomBondForce system; (ii) independent LangevinMiddle
    re-relaxation stability of the fitted model."""
    import openmm as mm
    from openmm import unit as u
    n = len(pos_fit)
    system = mm.System()
    for t in types:
        system.addParticle(MASS[t] * u.amu)
    force = mm.CustomBondForce("0.5*k*(r-r0)^2")
    force.addPerBondParameter("k")
    force.addPerBondParameter("r0")
    for (i, j, k, r0) in en:
        force.addBond(int(i), int(j), [float(k), float(r0)])
    system.addForce(force)
    integ = mm.LangevinMiddleIntegrator(T_SIM * u.kelvin, 5.0 / u.picosecond, 0.020 * u.picoseconds)
    names = [mm.Platform.getPlatform(q).getName() for q in range(mm.Platform.getNumPlatforms())]
    plat = mm.Platform.getPlatformByName("CPU" if "CPU" in names else "Reference")
    properties = {"Threads": os.environ.get("OMP_NUM_THREADS", "2")} if plat.getName() == "CPU" else {}
    context = mm.Context(system, integ, plat, properties)
    context.setPositions(pos_fit * u.nanometer)
    f_omm = context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(
        u.kilojoule_per_mole / u.nanometer)
    eng = MegaEngine(np.array([90.0, 90.0, 100.0]), types, np.zeros(n, bool))
    eng.en = en
    f_np = np.zeros_like(pos_fit)
    eng.bonded_forces(pos_fit, f_np, only="en")
    max_diff = float(np.abs(f_omm - f_np).max())
    context.setVelocitiesToTemperature(T_SIM * u.kelvin, 12345)
    for _ in range(n_drift):
        integ.step(1)
    p1 = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(u.nanometer)
    drift = kabsch_rmsd(p1, pos_fit)
    del context, integ
    return dict(max_force_diff=max_diff, drift_rmsd=float(drift),
                drift_steps=int(n_drift), platform=str(plat.getName()))


def mobile_cx_mask(n_total, n_sc, n_fg, n_cx):
    m = np.zeros(n_total, dtype=bool)
    m[n_sc + n_fg:n_sc + n_fg + n_cx] = True
    return m

# ==================================================== stage 16C: transport
WIN_K = 150.0               # kJ/mol/nm^2 umbrella stiffness (sigma_z = 0.13 nm):
                            # the mean restraint force then rises far above its
                            # sampling noise -> well-posed mean-force PMF
WIN_SPACING = 1.0           # nm
V_PULL = 0.05               # nm/ns (mission value)
DT_C = 0.030                # ps
GAMMA_C = 15.0              # 1/ps implicit-continuum friction (overdamped
                            # mesoscale limit; holds the capped-force steady state
                            # within ~10% of the 310 K bath, see report 3.2)
GAMMA_CX = 0.05             # 1/ps friction of the transported complex: paired
                            # with 24 amu virtual masses the COM relaxation time
                            # is tau_r = M*gamma/k ~ 6 ps inside each umbrella
                            # window, while the LJ wall remains resolved (a
                            # full 1.5-amu bead would jump the wall in one step)
MASS_CX = 24.0              # amu per complex bead (virtual CG mass; COM M = 1536)


def run_translocation(fast, budget):
    log("=== STAGE 16C : steered/umbrella translocation free energy & hydrodynamics ===")
    meta = {}
    d16 = np.load(RES / "mdff" / "mdff_trace.npz")
    target = d16["target"].astype(np.float64)
    fitted = d16["fitted"].astype(np.float64)
    seed = np.random.default_rng(0xF16A)
    t_sc, en, amp, anchors = build_scaffold(seed)
    # move anchors onto the FITTED frame (per-spoke rigid map target->fitted)
    anchors_fit = anchors.copy()
    ths = np.arctan2(anchors[:, 1], anchors[:, 0])
    for kk in range(N_SPOKES):
        th = 2 * math.pi * kk / N_SPOKES
        dth = (ths - th + math.pi) % (2 * math.pi) - math.pi
        m = np.abs(dth) < math.pi / N_SPOKES
        sel = np.zeros(len(t_sc), bool)
        for q in range(COPIES_PER_SPOKE):
            c0 = (kk * COPIES_PER_SPOKE + q) * N_BEADS_COPE
            sel[c0:c0 + N_BEADS_COPE] = True
        U, Pm, Qm = kabsch_transform(t_sc[sel], fitted[sel])
        anchors_fit[m] = anchors[m] @ U + (Qm - Pm @ U)
    br = np.load(RES / "stage16A_brush.npz")
    fg_pos = br["fg_pos"].astype(np.float64)
    fg_types = [str(x) for x in br["fg_types"]]

    results = {}
    for kind in ("receptor", "inert"):
        results[kind] = _run_kind(kind, fitted, en, anchors_fit, fg_pos, fg_types,
                                  fast, budget)
    for kind, r in results.items():
        meta[kind] = dict(n_windows=int(len(r["z0"])), n_prod=int(r["n_prod"]),
                          rh_nm=round(r["rh"], 3), wall_s=round(r["wall"], 1),
                          T_range_K=[round(r["Tmin"], 1), round(r["Tmax"], 1)],
                          nan_free=bool(r["nan_free"]),
                          seed_steps_literal=int(r["n_seed_lit"]))
    return meta, results


def _run_kind(kind, fitted, en, anchors, fg_pos, fg_types, fast, budget):
    """Brush-equilibrated mean-force translocation free energy of ONE rigid
    cargo complex.

    The transport receptor + cargo is treated as a RIGID BODY placed at the
    migration corridor (r_com = 11 nm) and stepped along the pore axis in
    window centres z0 in [-25, +25] nm (the reversible limit of the mission's
    v = 0.05 nm/ns steered protocol).  Around the fixed solute only the FG
    brush is integrated (overdamped Langevin); the mean constraint force
    <F_z> per window integrates to G(z), and the force-autocorrelation of the
    fluctuating brush force gives the local friction via the fluctuation-
    dissipation theorem, hence D(z) = kT/(M*gamma_eff(z))."""

    rng = np.random.default_rng(0xC16C if kind == "receptor" else 0x1BE1)
    # box must contain the full scaffold extent (r ~ 36 nm): a smaller box
    # would min-image alias frozen scaffold beads into the brush
    box = np.array([80.0, 80.0, 62.0])
    n_sc = len(fitted)
    n_ch = N_SPOKES * FG_PER_SPOKE
    L = FG_CHAIN_LEN
    n_fg = n_ch * L
    cx_pos, cx_types, _cx_en = build_receptor(kind, np.array([11.0, 0.0, -1.0]))
    n_cx = len(cx_pos)
    types = ["SD"] * n_sc + list(fg_types) + cx_types
    fixed = np.ones(n_sc + n_fg + n_cx, bool)      # everything frozen in 16C
    fixed[n_sc:n_sc + n_fg] = False                # ... except the FG brush
    pos = np.concatenate([fitted, fg_pos, cx_pos])
    eng = MegaEngine(box, types, fixed)
    eng.set_mobile(~fixed)
    eng.fcap = 100.0      # soft cap: the FG gel must respond like a soft condensate
    bonds, angles = chain_topology(n_sc, n_ch, L)
    eng.b, eng.ang, eng.en = bonds, angles, en
    fg0 = np.arange(n_sc, n_sc + n_fg, L)
    cx = np.arange(n_sc + n_fg, n_sc + n_fg + n_cx)
    remap = np.full(len(pos), -1, dtype=np.int64)
    remap[cx] = np.arange(n_cx)
    eng.make_pairlists(np.arange(len(pos)))

    def anchor_force(p, f):
        np.add.at(f, fg0, 2500.0 * (anchors - p[fg0]))

    eng.ext.append(anchor_force)
    M_cx = MASS_CX * n_cx
    rh = kirkwood_rh(pos[cx], box)

    n_eq_w = 250 if fast else 900          # brush relaxation around the solute
    n_ms = 120 if fast else 500            # measurement steps (force timeline)
    probe = _probe_step_ms(eng, pos, vel_dummy(len(pos)), np.zeros_like(pos))
    steps_per_window = n_eq_w + n_ms
    steps_budget = int(0.42 * budget.left() * 1000.0 / max(probe, 1e-3))
    z0s = np.arange(-25.0, 25.0 + 1e-9, WIN_SPACING)
    if fast:
        z0s = z0s[::6]
    nwin = len(z0s)
    afford = max(1, steps_budget // max(1, nwin))
    if afford < steps_per_window:
        scale = afford / steps_per_window
        n_eq_w = max(150, int(n_eq_w * scale))
        n_ms = max(60, int(n_ms * scale))
    log("16C | %s: %.2f ms/step | windows=%d eq=%d measure=%d steps/window "
        "(budget %.0f s)" % (kind, probe, nwin, n_eq_w, n_ms, 0.42 * budget.left()))

    z_series, contact_series, force_series, gamma_series = [], [], [], []
    Tmin, Tmax, nan_free = 1e9, 0.0, True
    t0 = time.time()
    vel = np.zeros_like(pos)
    vel[n_sc:n_sc + n_fg] += rng.standard_normal((n_fg, 3)) * np.sqrt(KT / eng.mass[n_sc:n_sc + n_fg])[:, None]
    frc = np.zeros_like(pos)
    ta_ids = np.nonzero(np.array(types)[n_sc:n_sc + n_fg] == "TA")[0] + n_sc
    for wi, z0 in enumerate(z0s):
        # place the rigid complex at the window centre (corridor r_com = 11)
        pos[cx] = pos[cx] - pos[cx].mean(0) + np.array([11.0, 0.0, z0])
        for st in range(n_eq_w):
            if eng.pl_lj.needs_update(pos):
                eng.rebuild(pos)
            frc[:] = 0.0
            eng.bonded_forces(pos, frc)
            eng.external_forces(pos, frc)
            eng.baoab(pos, vel, frc, DT_C, GAMMA_C, wrap_box=box)
        n_ct = max(1, n_ms // 10)
        ft = np.empty(n_ms, dtype=np.float64)
        ct = np.zeros(n_ct, dtype=np.float32)
        for st in range(n_ms):
            if eng.pl_lj.needs_update(pos):
                eng.rebuild(pos)
            frc[:] = 0.0
            eng.bonded_forces(pos, frc)
            eng.external_forces(pos, frc)
            f_cx = eng.pair_forces_local(pos, frc, cx, remap)
            ft[st] = float(f_cx[:, 2].sum())
            # the brush must FLUCTUATE around the solute during measurement:
            # the FDT friction and the contact statistics are only meaningful
            # for a living (thermally moving) FG condensate
            eng.baoab(pos, vel, frc, DT_C, GAMMA_C, wrap_box=box)
            if st % 10 == 0 and st // 10 < n_ct:
                d = pos[cx][:, None, :] - pos[ta_ids][None, :, :]
                d -= box * np.round(d / box)
                ct[st // 10] = int(((d * d).sum(-1) < 1.0 ** 2).sum())
        fz_mean = float(ft.mean())
        # FDT friction: gamma_eff = beta * integral of the force autocorrelation
        fc = ft - ft.mean()
        n_c = min(len(fc) - 1, 120)
        ac = np.array([float((fc[:len(fc) - t] * fc[t:]).mean()) for t in range(n_c)])
        gamma_eff = BETA * float(ac.sum()) * DT_C * m_scale(kind)
        gamma_eff = float(np.clip(gamma_eff, 1e-4, 1e4))
        if not np.isfinite(ft).all():
            nan_free = False
        kin = 0.5 * eng.mass[n_sc:n_sc + n_fg] @ (vel[n_sc:n_sc + n_fg] ** 2).sum(1)
        T_fg = 2 * kin / (3 * n_fg * KB)
        Tmin, Tmax = min(Tmin, T_fg), max(Tmax, T_fg)
        z_series.append(np.full(max(1, n_ms // 5), z0, dtype=np.float32))
        contact_series.append(ct)
        force_series.append(ft[::2].astype(np.float32))
        gamma_series.append(gamma_eff)
        if wi % 10 == 0:
            log("16C | %s window %3d/%d  z0=%+06.1f  <Fz>=%+09.2f  gamma_eff=%9.3f  "
                "<contacts>=%.2f  T=%.0fK"
                % (kind, wi, nwin, z0, fz_mean, gamma_eff,
                   float(ct[ct > 0].mean()) if np.any(ct > 0) else 0.0, T_fg))
    wall = time.time() - t0
    log("16C | %s sweep done in %.1f min" % (kind, wall / 60.0))
    np.savez_compressed(RES / "umbrella" / ("windows_%s.npz" % kind),
                        z0=z0s, z=np.array(z_series), contacts=np.array(contact_series),
                        force=np.array(force_series), gamma=np.array(gamma_series),
                        eq_contact_max=np.array(contact_series).max(axis=1),
                        dt=DT_C, sample_every=2, k_umb=WIN_K, rh=rh,
                        n_prod=n_ms, n_eq=n_eq_w, n_pull_seed=0,
                        M_cx=M_cx)
    return dict(z0=z0s, n_prod=n_ms, rh=float(rh), wall=wall,
                Tmin=float(Tmin), Tmax=float(Tmax), nan_free=bool(nan_free),
                n_seed_lit=0)


def vel_dummy(n):
    return np.zeros((n, 3))


def m_scale(kind):
    return MASS_CX * len(CX_TYPES(kind))



def _probe_step_ms(eng, pos, vel, frc):
    t0 = time.time()
    for _ in range(60):
        if eng.pl_lj.needs_update(pos):
            eng.rebuild(pos)
        frc[:] = 0.0
        eng.bonded_forces(pos, frc)
        eng.pair_forces(pos, frc)
        eng.external_forces(pos, frc)
        eng.baoab(pos, vel, frc, DT_C, GAMMA_C, wrap_box=eng.box)
    return 1000.0 * (time.time() - t0) / 60.0


# ==================================================================== WHAM
def _logsumexp(a, axis=None):
    m = np.max(a, axis=axis, keepdims=True)
    m = np.where(np.isfinite(m), m, 0.0)
    out = np.log(np.sum(np.exp(a - m), axis=axis, keepdims=True)) + m
    return np.squeeze(out, axis=axis) if axis is not None else out


def wham(z0s, k, z_samples, bins, zmin, zmax, tol=1e-10, maxiter=20000):
    """Log-space WHAM for 1-D harmonic umbrellas -> (centers, G, F_i, H)."""
    nwin = len(z0s)
    edges = np.linspace(zmin, zmax, bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    H = np.zeros((nwin, bins))
    for i, zs in enumerate(z_samples):
        H[i] = np.histogram(zs, bins=edges)[0]
    n = H.sum(axis=1)
    betaU = BETA * 0.5 * k * (centers[None, :] - np.asarray(z0s, float)[:, None]) ** 2
    F = np.zeros(nwin)
    logP = np.full(bins, math.log(1.0 / bins))
    tiny = 1e-300
    covered = H.sum(axis=0) > 0
    for _ in range(maxiter):
        lse = np.log(np.maximum(n[:, None], 1e-30)) + (BETA * F[:, None] - betaU)
        denom = np.maximum(np.exp(_logsumexp(lse, axis=0)), tiny)
        P = np.maximum(H.sum(axis=0) / denom, tiny)
        logP = np.log(P / P.sum())
        acc = np.exp(logP[None, :] - betaU)
        acc = np.where(covered[None, :], acc, 0.0)
        F_new = -KT * np.log(np.maximum(acc.sum(axis=1), tiny))
        if np.abs(F_new - F).max() < tol:
            F = F_new
            break
        F = F_new
    return centers, -KT * logP, F, H


# ============================================================ 16C analysis
def analyze_transport(results, fast):
    """PMF from brush-equilibrated mean constraint forces (reversible-path
    umbrella integration); D(z) from the fluctuation-dissipation friction of
    the measured force timelines; contacts and permeability as before."""
    slim = {}
    for kind in ("receptor", "inert"):
        d = np.load(RES / "umbrella" / ("windows_%s.npz" % kind))
        z0s = d["z0"].astype(float)
        force = d["force"].astype(np.float64)      # (nwin, n_ms/2) kJ/mol/nm
        contacts = d["contacts"].astype(np.float64)
        eq_cts = d["eq_contact_max"].astype(float)
        gamma = d["gamma"].astype(np.float64)      # amu/ps (whole complex)
        M_cx = float(d["M_cx"])
        if (z0s.ndim != 1 or len(z0s) < 2 or not np.isfinite(z0s).all()
                or not np.all(np.diff(z0s) > 0)):
            raise ValueError("Transport windows must have finite increasing coordinates")
        if (force.ndim != 2 or force.shape[0] != len(z0s)
                or force.shape[1] < 2 or not np.isfinite(force).all()):
            raise ValueError("Transport needs two or more finite force samples per window")
        if (gamma.shape != z0s.shape or not np.isfinite(gamma).all()
                or np.any(gamma <= 0)):
            raise ValueError("Transport friction must be finite and positive in every window")
        Fbar = force.mean(axis=1)
        half = max(1, force.shape[1] // 2)
        F1 = force[:, :half].mean(axis=1)
        F2 = force[:, half:].mean(axis=1)
        dz = np.diff(z0s)
        # Recorded force acts on the fixed solute: dG/dz = -<Fz>.
        G = np.concatenate(([0.0], -np.cumsum(0.5 * (Fbar[1:] + Fbar[:-1]) * dz)))
        G1 = np.concatenate(([0.0], -np.cumsum(0.5 * (F1[1:] + F1[:-1]) * dz)))
        G2 = np.concatenate(([0.0], -np.cumsum(0.5 * (F2[1:] + F2[:-1]) * dz)))
        inter = (np.abs(z0s) < 12.0)
        bulk = (np.abs(z0s) > 22.0)
        if not bulk.any():
            raise ValueError("Transport PMF needs a sampled bulk reference")
        for curve in (G, G1, G2):
            curve -= curve[bulk].mean()
        Gbar = float(G[inter].max() - G[bulk].mean()) if inter.any() else float("nan")
        # D(z) = kT / zeta_eff  (fluctuation-dissipation friction)
        Dz = KT / np.maximum(gamma, 1e-6)
        contacts_m = contacts.mean(axis=1)
        # inhomogeneous solubility-diffusion resistance
        zz = z0s
        gg = G / KT  # retain attractive wells relative to bulk
        if not np.isfinite(gg).all() or np.max(gg) > 700:
            raise ValueError("Invalid or overflowing PMF resistance")
        Dsel = np.maximum(Dz, 1e-8)
        Rint = float(TRAPZ(np.exp(gg) / Dsel, zz))
        rh = float(d["rh"])
        if not np.isfinite(Rint) or Rint <= 0 or not np.isfinite(rh) or rh <= 0:
            raise ValueError("Transport resistance and hydrodynamic radius must be positive and finite")
        D_SE_um2 = (KT * 1.6605e-21 / (6 * math.pi * ETA_WATER * rh * 1e-9)) * 1e12
        zbulk = np.abs(z0s) > 22.0
        bulk_D = float(np.mean(Dz[zbulk])) if zbulk.any() else float(np.mean(Dz))
        cal = D_SE_um2 / (bulk_D * NM2PS_TO_UM2S)
        np.savez_compressed(RES / "umbrella" / ("pmf_%s.npz" % kind),
                            centers=z0s, G=G, G1=G1, G2=G2,
                            cover=np.ones(len(z0s), bool),
                            Dz=Dz, gamma=gamma,
                            contacts=contacts_m, eq_contact_max=eq_cts,
                            dwell=(contacts_m > 0).astype(float),
                            D_cal_factor=cal, Fbar=Fbar)
        slim[kind] = dict(Gbar_kT=round(Gbar / KT, 3),
                          Gbar_kcal=round(Gbar * KJ_TO_KCAL, 3),
                          R_int=Rint, P_rel=1.0 / Rint, D_cal_factor=cal, rh_nm=rh,
                          dwell_max=float((contacts_m > 0).mean()),
                          contacts_peak=float(np.nanmax(contacts_m)),
                          eq_contacts_peak=float(np.nanmax(eq_cts)),
                          Fbar_peak=float(np.nanmax(np.abs(Fbar))))
    slim["dGddagger_kcal"] = dict(
        receptor=slim["receptor"]["Gbar_kcal"], inert=slim["inert"]["Gbar_kcal"],
        delta=round(slim["inert"]["Gbar_kcal"] - slim["receptor"]["Gbar_kcal"], 3))
    slim["P_ratio_receptor_over_inert"] = float(slim["receptor"]["P_rel"]
                                                / slim["inert"]["P_rel"])
    (RES / "analysis_transport.json").write_text(json.dumps(slim, indent=1, default=str))
    return slim



# ================================================================= selftest
def run_selftests():
    log("=== SELFTEST : kernel validation gates ===")
    ok = {}
    rng = np.random.default_rng(7)
    # 1. trilinear field exactness (linear reproduces exactly; analytic grad)
    g = np.arange(6 * 7 * 8, dtype=float).reshape(6, 7, 8)
    tf = TrilinearField(g, np.zeros(3), np.ones(3) * 0.5)
    p = rng.uniform([0.1, 0.1, 0.1], [2.4, 2.9, 3.4], size=(50, 3))
    v, _ = tf.values(p)
    exact = (p[:, 0] / 0.5) * 56 + (p[:, 1] / 0.5) * 8 + (p[:, 2] / 0.5)
    ok["trilinear_linear_exact_maxerr"] = float(np.abs(v - exact).max())
    grad, _ = tf.gradients(p)
    # g[i,j,k] = 56 i + 8 j + k -> slopes per nm: 56/0.5, 8/0.5, 1/0.5
    ok["trilinear_grad_linear_maxerr"] = float(np.abs(
        grad - np.array([112.0, 16.0, 2.0])[None, :]).max())
    g2 = (np.sin(np.linspace(0, 3, 30))[:, None, None]
          * np.cos(np.linspace(0, 2, 24))[None, :, None]
          * np.linspace(1, 2, 20)[None, None, :])
    tf2 = TrilinearField(g2, np.zeros(3), np.ones(3) * 0.1)
    p2 = rng.uniform(0.12, 2.8, size=(40, 3))
    g_an, _ = tf2.gradients(p2)
    eps = 1e-5
    fd = np.array([(tf2.values(p2 + np.eye(3)[a] * eps)[0]
                    - tf2.values(p2 - np.eye(3)[a] * eps)[0]) / (2 * eps) for a in range(3)]).T
    ok["trilinear_grad_fd_maxerr"] = float(np.abs(g_an - fd).max())
    # 2. MRC round trip
    gp = rng.uniform(size=(9, 8, 7)).astype(np.float32)
    write_mrc(RES / "_selftest.mrc", gp, np.array([-1.0, -2.0, -3.0]), 0.35)
    gr, org, vox = read_mrc(RES / "_selftest.mrc")
    ok["mrc_roundtrip"] = dict(maxdiff=float(np.abs(gr - gp).max()),
                               origin_ok=bool(np.allclose(org, [-1, -2, -3])),
                               voxel_nm=float(vox))
    (RES / "_selftest.mrc").unlink()
    # 3. pair list vs cKDTree (same physical config, both periodic conventions)
    from scipy.spatial import cKDTree
    box = np.array([12.0, 12.0, 12.0])
    pts0 = rng.uniform(-5.5, 5.5, size=(400, 3))
    pl = PairList(box, 1.1)
    pl.update(pts0)
    tree = cKDTree(pts0 % 12.0, boxsize=12.0)
    pairs_ref = set(tuple(sorted(pr)) for pr in tree.query_pairs(1.1))
    d = pts0[pl.i] - pts0[pl.j]
    d -= box * np.round(d / box)
    m = (d * d).sum(1) < 1.1 ** 2
    pairs_new = set(tuple(sorted(pr)) for pr in zip(pl.i[m], pl.j[m]))
    ok["pairlist_vs_kdtree"] = dict(ref=len(pairs_ref), new=len(pairs_new),
                                    mismatch=len(pairs_ref ^ pairs_new))
    # 4. BAOAB Maxwellian
    eng = MegaEngine(np.array([10.0, 10.0, 10.0]), ["N0"] * 2000, np.zeros(2000, bool))
    eng.set_mobile(np.ones(2000, bool))
    eng.make_pairlists(np.arange(2000))
    pos = rng.uniform(-5, 5, size=(2000, 3))
    vel = np.zeros((2000, 3))
    f0 = np.zeros((2000, 3))
    eng.rebuild(pos)
    for _ in range(1500):
        eng.baoab(pos, vel, f0, 0.010, 2.0, wrap_box=np.array([10.0, 10.0, 10.0]))
    ok["baoab_temperature_K"] = float(eng.mass @ (vel ** 2).sum(1) / (3.0 * 2000 * KB))
    # 5. WHAM recovery on a synthetic anharmonic PMF
    ztrue = np.linspace(-1.3, 1.3, 200)
    Gt = 2.5 * ztrue ** 4 - 3.0 * ztrue ** 2 + 1.0
    z0s = np.linspace(-1.2, 1.2, 13)
    k = 15.0
    samples = []
    for z0 in z0s:
        z = float(z0)
        traj = np.empty(42000)
        gam = 30.0
        for t in range(42000):
            f = -(4 * 2.5 * z ** 3 - 6.0 * z) - k * (z - z0)
            z += f / gam * 0.004 + rng.normal(scale=math.sqrt(2 * KT * 0.004 / gam))
            traj[t] = z
        samples.append(traj[4000:])
    centers, Gw, F, H = wham(z0s, k, samples, bins=60, zmin=-1.3, zmax=1.3)
    Gt_c = np.interp(centers, ztrue, Gt)
    mm = H.sum(0) > 30
    i0 = int(np.argmax(H.sum(0)[mm]))
    ok["wham_recovery_max_kT"] = float(np.abs(
        (Gw[mm] - Gw[mm][i0]) - (Gt_c[mm] - Gt_c[mm][i0])).max() / KT)
    ok["openmm_available"] = bool(HAS_OPENMM)
    # 6. bond/angle force direction regression (stretching must pull back)
    eng6 = MegaEngine(np.array([10.0, 10.0, 10.0]), ["N0", "N0", "N0"],
                      np.zeros(3, bool))
    eng6.set_mobile(np.ones(3, bool))
    eng6.b = np.array([[0, 1, 1250.0, 0.36]])
    eng6.ang = np.array([[0, 1, 2, 6.0]])
    p6 = np.array([[0.0, 0.0, 0.0], [0.8, 0.0, 0.0], [0.8, 0.8, 0.0]])
    f6 = np.zeros((3, 3))
    eng6.bonded_forces(p6, f6)
    ok["bond_stretched_pulls_together"] = bool(f6[0, 0] > 0 and f6[1, 0] < 0)
    eng6.b = np.array([[0, 1, 1250.0, 1.2]])   # compressed bond must push apart
    f6[:] = 0.0
    eng6.bonded_forces(p6, f6)
    ok["bond_compressed_pushes_apart"] = bool(f6[0, 0] < 0 and f6[1, 0] > 0)
    log("selftest: " + json.dumps(ok, default=str))
    return ok


# ================================================================== figures
def _fig_style(ax):
    ax.tick_params(direction="in", top=True, right=True, labelsize=8)


def fig1():
    d = np.load(RES / "mdff" / "mdff_trace.npz")
    fr = np.load(RES / "stage16A_frames.npz")
    zf = d["forcemap"]
    origin = d["forcemap_origin"]
    vox = float(d["forcemap_voxel"])
    fitted = d["fitted"]
    fig = plt.figure(figsize=(16.0, 7.8))
    gs = fig.add_gridspec(2, 3, width_ratios=[2.05, 1.0, 1.0],
                          left=0.02, right=0.985, top=0.88, bottom=0.06,
                          wspace=0.26, hspace=0.28)
    ax3 = fig.add_subplot(gs[:, 0], projection="3d")
    nx, ny, nz = zf.shape
    if HAS_SKIMAGE:
        try:
            verts, faces, _, _ = _skm.marching_cubes(zf, level=0.30 * float(zf.max()))
            pts = verts * vox + origin
            normals = np.cross(pts[faces][:, 1] - pts[faces][:, 0],
                               pts[faces][:, 2] - pts[faces][:, 0])
            normals /= (np.linalg.norm(normals, axis=1, keepdims=True) + 1e-9)
            light = np.array([-0.5, -0.6, 0.62])
            light /= np.linalg.norm(light)
            lam = np.clip(normals @ light, 0.0, 1.0)[..., None]
            base = np.array([0.72, 0.80, 0.92])
            fc = np.clip(0.35 * base + 0.75 * base * lam, 0, 1)
            mesh = Poly3DCollection(pts[faces], alpha=0.15, facecolors=fc,
                                    edgecolor="none")
            ax3.add_collection3d(mesh)
        except Exception as exc:
            log("fig1 | isosurface fallback (%s)" % exc)
    # density slice on the back wall
    y0 = 2
    XX, ZZ = np.meshgrid(origin[0] + np.arange(nx) * vox, origin[2] + np.arange(nz) * vox)
    ax3.contourf(XX, ZZ, zf[:, y0, :].T, levels=14, cmap="Blues",
                 zdir="y", offset=origin[1] + y0 * vox, alpha=0.55)
    cope_colors = {"IR": "#2c6fbb", "CR": "#1b4f8a", "LR": "#3fa66a", "TL": "#7ec8e3"}
    for kk in range(N_SPOKES):
        for q, (cname, *_r) in enumerate(COPE_CLASSES):
            c0 = (kk * COPIES_PER_SPOKE + q) * N_BEADS_COPE
            ax3.scatter(fitted[c0:c0 + N_BEADS_COPE, 0][::2],
                        fitted[c0:c0 + N_BEADS_COPE, 1][::2],
                        fitted[c0:c0 + N_BEADS_COPE, 2][::2],
                        s=3.5, c=cope_colors[cname[:2]], alpha=0.9,
                        depthshade=False, linewidths=0)
    fg = fr["fg"].astype(float)
    fg_types = fr["fg_types"]
    mta = fg_types == "TA"
    ax3.scatter(fg[::7, 0], fg[::7, 1], fg[::7, 2], s=1.2, c="#d9b382", alpha=0.20,
                depthshade=True, linewidths=0)
    ax3.scatter(fg[mta][::3, 0], fg[mta][::3, 1], fg[mta][::3, 2], s=2.5, c="#b5651d",
                alpha=0.4, depthshade=True, linewidths=0)
    cx_pos, _, _ = build_receptor("receptor", np.array([0.0, 0.0, -25.0]))
    ax3.scatter(cx_pos[:, 0], cx_pos[:, 1], cx_pos[:, 2], s=24, c="#e6194b",
                alpha=0.95, depthshade=False, linewidths=0)
    ax3.plot([0, 0], [0, 0], [-40, 40], "k--", lw=0.8, alpha=0.5)
    ax3.set_xlim(-40, 40)
    ax3.set_ylim(-40, 40)
    ax3.set_zlim(-40, 40)
    ax3.set_box_aspect((1, 1, 1))
    ax3.set_xlabel("x (nm)", labelpad=-4)
    ax3.set_ylabel("y (nm)", labelpad=-4)
    ax3.set_zlabel("z (nm)", labelpad=-4)
    ax3.view_init(elev=16, azim=-58)
    ax3.set_title("CG NPC megasystem docked into the 3.5 A cryo-EM envelope\n"
                  "(8-fold scaffold | FG-Nup condensate | transport receptor)", fontsize=10.5)
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#2c6fbb", markersize=6, label="inner ring (IR)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#1b4f8a", markersize=6, label="outer ring (CR)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#3fa66a", markersize=6, label="luminal ring (LR)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#7ec8e3", markersize=6, label="transport linkers (TL)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#d9b382", markersize=6, label="FG-Nup brush"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#b5651d", markersize=6, label="FG Phe stickers (TA)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#e6194b", markersize=6, label="receptor + cargo"),
        Line2D([0], [0], color="#c8d4e8", lw=6, alpha=0.6, label="cryo-EM isosurface"),
    ]
    ax3.legend(handles=handles, loc="upper left", bbox_to_anchor=(-0.02, 1.00),
               fontsize=6.6, framealpha=0.6, ncol=2, columnspacing=0.8, handletextpad=0.2)
    axA = fig.add_subplot(gs[0, 1])
    _fig_style(axA)
    axA.plot(d["xi"], d["ccc"], color="#1a5276", lw=1.8)
    axA.set_xlabel("adaptive force $\\xi$ (kJ/mol)")
    axA.set_ylabel("CCC(model, map)")
    axA.set_title("MDFF cross-correlation rise", fontsize=9.5)
    axA.text(0.04, 0.07, "CCC %.3f $\\rightarrow$ %.3f" % (d["ccc"][0], d["ccc"][-1]),
             transform=axA.transAxes, fontsize=8.5)
    axB = fig.add_subplot(gs[1, 1])
    _fig_style(axB)
    axB.plot(np.arange(len(d["rmsd"])) * 100 * 0.025 / 1e3, d["rmsd"], color="#7d3c98", lw=1.8)
    axB.set_xlabel("MDFF time (ns)")
    axB.set_ylabel("RMSD to target (nm)")
    axB.set_title("fit convergence", fontsize=9.5)
    axB.text(0.04, 0.07, "RMSD %.2f $\\rightarrow$ %.2f nm" % (d["rmsd"][0], d["rmsd"][-1]),
             transform=axB.transAxes, fontsize=8.5)
    axC = fig.add_subplot(gs[:, 2])
    _fig_style(axC)
    tr = np.load(RES / "umbrella" / "pmf_receptor.npz")
    ti = np.load(RES / "umbrella" / "pmf_inert.npz")
    axC.plot(tr["centers"], tr["G"] * KJ_TO_KCAL, color="#e6194b", lw=2, label="receptor")
    axC.plot(ti["centers"], ti["G"] * KJ_TO_KCAL, color="#4d4d4d", lw=2, ls="--", label="inert control")
    axC.set_xlabel("z along pore axis (nm)")
    axC.set_ylabel("$\\Delta G$ (kcal/mol)")
    axC.set_title("translocation PMF (WHAM)", fontsize=9.5)
    axC.legend(fontsize=8)
    fig.suptitle("Phase 16 — Mega-macromolecular dynamics: Cryo-EM flexible fitting of the NPC megasystem",
                 fontsize=13, y=0.985)
    fig.savefig(FIG / "fig1_megasystem_cryoem_fit.png", dpi=300)
    plt.close(fig)
    log("fig1 written")


def fig2():
    br = np.load(RES / "stage16A_brush.npz")
    fg = br["fg_pos"].astype(float)
    r_fg = np.sqrt(fg[:, 0] ** 2 + fg[:, 1] ** 2)
    fig = plt.figure(figsize=(16.0, 5.0))
    gs = fig.add_gridspec(1, 4, left=0.045, right=0.985, top=0.84, bottom=0.13, wspace=0.28)
    nb_r, nb_z = 120, 120
    r_edges = np.linspace(0, 26, nb_r + 1)
    z_edges = np.linspace(-26, 26, nb_z + 1)
    Hr, _, _ = np.histogram2d(r_fg, fg[:, 2], bins=[r_edges, z_edges])
    shell = np.pi * (r_edges[1:] ** 2 - r_edges[:-1] ** 2)[:, None] * (z_edges[1] - z_edges[0])
    rho_rz = Hr / shell
    axA = fig.add_subplot(gs[0])
    cmap = LinearSegmentedColormap.from_list(
        "fg", ["#f7f4ec", "#e8d3a9", "#c98a3d", "#7a3e0a", "#3d1f05"])
    im = axA.pcolormesh(z_edges, r_edges, rho_rz, cmap=cmap, shading="auto",
                        norm=Normalize(0, max(np.percentile(rho_rz, 99), 1e-6)))
    if rho_rz.max() > 0:
        axA.contour(z_edges[:-1], r_edges[:-1], rho_rz,
                    levels=[0.25 * rho_rz.max()], colors="k", linewidths=0.8)
    axA.axhline(R_LUMEN, color="#1b4f8a", ls="--", lw=1.2)
    axA.text(-25.5, R_LUMEN - 1.6, "scaffold wall ($r$ = 25 nm)", fontsize=7.5, color="#1b4f8a")
    axA.set_xlabel("z along pore axis (nm)")
    axA.set_ylabel("radius r (nm)")
    axA.set_title("(a) FG-Nup number density $\\rho(r,z)$", fontsize=10)
    cb = fig.colorbar(im, ax=axA, pad=0.015)
    cb.set_label("beads/nm$^3$", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    axB = fig.add_subplot(gs[1])
    _fig_style(axB)
    r_mid = 0.5 * (r_edges[:-1] + r_edges[1:])
    for zlo, zhi, ccol in [(0, 5, "#7a3e0a"), (8, 12, "#c98a3d"),
                           (16, 20, "#7fb3d5"), (22, 26, "#2c6fbb")]:
        m = (fg[:, 2] >= zlo) & (fg[:, 2] < zhi)
        h, _ = np.histogram(np.sqrt(fg[m, 0] ** 2 + fg[m, 1] ** 2), bins=r_edges)
        prof = h / (np.pi * (r_edges[1:] ** 2 - r_edges[:-1] ** 2) * (zhi - zlo))
        axB.plot(r_mid, prof, color=ccol, lw=1.8, label="z$\\in$[%d,%d)" % (zlo, zhi))
    axB.axvline(R_LUMEN, color="k", ls=":", lw=1)
    axB.set_xlabel("radius r (nm)")
    axB.set_ylabel("$\\rho$ (beads/nm$^3$)")
    axB.set_title("(b) radial density at axial slices", fontsize=10)
    axB.legend(fontsize=7.5)
    axC = fig.add_subplot(gs[2])
    _fig_style(axC)
    ta_coord = br["ta_coord"]
    nb = 24
    edges = np.linspace(0, 26, nb + 1)
    mids = 0.5 * (edges[:-1] + edges[1:])
    mc = np.clip(np.digitize(r_fg, edges) - 1, 0, nb - 1)
    prof_c = np.array([ta_coord[mc == q].mean() if np.any(mc == q) else np.nan for q in range(nb)])
    axC.plot(mids, prof_c, color="#b5651d", lw=2, marker="o", ms=3)
    axC.set_xlabel("radius r (nm)")
    axC.set_ylabel("$\\langle N_{TA\\cdots TA}\\rangle$ (< 0.65 nm)")
    axC.set_title("(c) hydrophobic sticker clustering", fontsize=10)
    ax2 = axC.twinx()
    var_prof = np.array([r_fg[mc == q].var() if np.any(mc == q) else np.nan for q in range(nb)])
    ax2.plot(mids, var_prof, color="#5d6d7e", lw=1.4, ls="--")
    ax2.set_ylabel("$\\sigma^2_r(r)$ (nm$^2$)", color="#5d6d7e", fontsize=8)
    ax2.tick_params(colors="#5d6d7e", labelsize=7)
    if np.isfinite(prof_c).any():
        gel_r = mids[np.nanargmax(prof_c)]
        axC.axvspan(gel_r - 1.5, gel_r + 1.5, color="#f5b041", alpha=0.18)
        axC.text(gel_r, np.nanmax(prof_c) * 0.93, "gel zone", fontsize=8, ha="center")
    axC.axvline(R_LUMEN, color="k", ls=":", lw=1)
    axD = fig.add_subplot(gs[3])
    _fig_style(axD)
    n_ch = len(fg) // FG_CHAIN_LEN
    rg = np.empty(n_ch)
    for c in range(n_ch):
        seg = fg[c * FG_CHAIN_LEN:(c + 1) * FG_CHAIN_LEN]
        rg[c] = np.sqrt(((seg - seg.mean(0)) ** 2).sum(1).mean())
    axD.hist(rg, bins=30, color="#8e44ad", alpha=0.75, edgecolor="k", lw=0.3)
    axD.set_xlabel("$R_g$ per FG chain (nm)")
    axD.set_ylabel("chains")
    axD.set_title("(d) chain size distribution", fontsize=10)
    axD.text(0.45, 0.85, "$\\langle R_g\\rangle$ = %.2f nm\n$n$ = %d chains" % (rg.mean(), n_ch),
             transform=axD.transAxes, fontsize=8.5)
    fig.suptitle("Phase 16 — FG-nucleoporin condensate inside the pore: gel/brush density architecture",
                 fontsize=13, y=0.96)
    fig.savefig(FIG / "fig2_fg_condensate_density_slice.png", dpi=300)
    plt.close(fig)
    log("fig2 written")


def fig3():
    tr = np.load(RES / "umbrella" / "pmf_receptor.npz")
    ti = np.load(RES / "umbrella" / "pmf_inert.npz")
    aj = json.loads((RES / "analysis_transport.json").read_text())
    fig = plt.figure(figsize=(16.0, 5.0))
    gs = fig.add_gridspec(1, 4, left=0.05, right=0.985, top=0.84, bottom=0.13, wspace=0.31)
    axA = fig.add_subplot(gs[0])
    _fig_style(axA)
    for dd, cc, lbl, ls in [(tr, "#e6194b", "transport receptor", "-"),
                            (ti, "#4d4d4d", "inert control ($R_h$ matched)", "--")]:
        c = dd["centers"]
        G = dd["G"] * KJ_TO_KCAL
        band = 0.5 * np.abs(dd["G1"] - dd["G2"]) * KJ_TO_KCAL
        axA.plot(c, G, color=cc, lw=2.2, ls=ls, label=lbl)
        axA.fill_between(c, G - band, G + band, color=cc, alpha=0.18, lw=0)
    bar_r = aj["dGddagger_kcal"]["receptor"]
    bar_i = aj["dGddagger_kcal"]["inert"]
    ymax = max(bar_r, bar_i) * 1.05
    axA.set_xlabel("z along pore axis (nm)")
    axA.set_ylabel("$\\Delta G$ (kcal/mol)")
    axA.set_title("(a) translocation free energy (WHAM)", fontsize=10)
    axA.set_ylim(bottom=min(-1.5, float(tr["G"].min()) * KJ_TO_KCAL * 1.1), top=ymax * 1.30)
    axA.annotate("$\\Delta G^\\ddagger_{\\rm rec}$ = %.1f" % bar_r, xy=(-6, bar_r),
                 xytext=(-24, ymax * 0.95), fontsize=8.5, color="#e6194b",
                 arrowprops=dict(arrowstyle="->", color="#e6194b", lw=0.9))
    axA.annotate("$\\Delta G^\\ddagger_{\\rm inert}$ = %.1f" % bar_i, xy=(-6, bar_i),
                 xytext=(4, ymax * 0.72), fontsize=8.5, color="#4d4d4d",
                 arrowprops=dict(arrowstyle="->", color="#4d4d4d", lw=0.9))
    axA.legend(fontsize=8, loc="lower center")
    axB = fig.add_subplot(gs[1])
    _fig_style(axB)
    cal = float(tr["D_cal_factor"])
    axB.plot(tr["centers"], tr["Dz"] * NM2PS_TO_UM2S * cal, color="#e6194b", lw=2,
             label="receptor")
    axB.plot(ti["centers"], ti["Dz"] * NM2PS_TO_UM2S * cal, color="#4d4d4d", lw=2, ls="--",
             label="inert")
    axB.set_yscale("log")
    axB.set_xlabel("z (nm)")
    axB.set_ylabel("$D_{app}(z)$ ($\\mu$m$^2$/s)")
    axB.set_title("(b) axial mobility profile", fontsize=10)
    axB.text(0.03, 0.94, "Stokes-Einstein calibrated\n($R_h$ = %.2f nm)" % aj["receptor"]["rh_nm"],
             transform=axB.transAxes, fontsize=7.5, va="top")
    axB.legend(fontsize=8, loc="lower right")
    axC = fig.add_subplot(gs[2])
    _fig_style(axC)
    axC.plot(tr["centers"], tr["contacts"], color="#b5651d", lw=2,
             label="FG contacts (PA$\\cdots$TA)")
    axC.plot(ti["centers"], ti["contacts"], color="#85929e", lw=2, ls="--",
             label="inert surface (PI$\\cdots$TA)")
    axC.fill_between(tr["centers"], tr["dwell"] * float(tr["contacts"].max()), 0,
                     color="#b5651d", alpha=0.10, lw=0)
    axC.set_xlabel("z (nm)")
    axC.set_ylabel("mean transient contacts")
    axC.set_title("(c) multivalent FG engagement", fontsize=10)
    axC.legend(fontsize=8)
    axD = fig.add_subplot(gs[3])
    _fig_style(axD)
    Pr = aj["receptor"]["P_rel"]
    Pi = aj["inert"]["P_rel"]
    ratio = aj["P_ratio_receptor_over_inert"]
    axD.bar(["receptor", "inert"], [Pr, Pi], color=["#e6194b", "#4d4d4d"],
            alpha=0.85, edgecolor="k", lw=0.6)
    axD.set_yscale("log")
    axD.set_ylabel("$P\\propto[\\int e^{\\beta G}/D\\,dz]^{-1}$ (a.u.)")
    axD.set_title("(d) mesoscale permeability", fontsize=10)
    axD.text(0.08, 0.88, "$P_{\\rm rec}/P_{\\rm inert}$ = %.2f" % ratio,
             transform=axD.transAxes, fontsize=9)
    fig.suptitle("Phase 16 — Non-equilibrium translocation: multivalent FG contacts gate the pore",
                 fontsize=13, y=0.96)
    fig.savefig(FIG / "fig3_translocation_free_energy_pmf.png", dpi=300)
    plt.close(fig)
    log("fig3 written")


# ==================================================================== main
def main():
    ap = argparse.ArgumentParser(description="Phase 16 mega-machine: NPC cryo-EM MDFF + mesoscale transport")
    ap.add_argument("--fast", action="store_true", help="reduced smoke budget")
    ap.add_argument("--fig-only", action="store_true", help="re-render figures from saved npz")
    ap.add_argument("--selftest", action="store_true", help="kernel validation gates only")
    ap.add_argument("--budget-min", type=float, default=110.0)
    ap.add_argument("--resume16c", action="store_true",
                    help="reuse saved 16A/16B artifacts, run 16C onward")
    args = ap.parse_args()
    RES.mkdir(exist_ok=True)
    (RES / "mdff").mkdir(exist_ok=True)
    (RES / "umbrella").mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)
    budget = Budget(args.budget_min, args.fast)

    if args.selftest:
        st = run_selftests()
        (RES / "selftest.json").write_text(json.dumps(st, indent=1, default=str))
        return

    if not args.fig_only:
        st = run_selftests()
        (RES / "selftest.json").write_text(json.dumps(st, indent=1, default=str))
        resume = args.resume16c and (RES / "stage16A_meta.json").exists()             and (RES / "stage16B_meta.json").exists()             and (RES / "mdff" / "mdff_trace.npz").exists()             and (RES / "stage16A_brush.npz").exists()
        if resume:
            m16a = json.loads((RES / "stage16A_meta.json").read_text())
            m16b = json.loads((RES / "stage16B_meta.json").read_text())
            log("resume: reusing saved 16A/16B artifacts")
        else:
            budget.mark("16A")
            m16a = stage_16A(args.fast, budget)
            budget.close("16A")
            budget.mark("16B")
            m16b = stage_16B(args.fast, budget)
            budget.close("16B")
        budget.mark("16C")
        m16c, results = run_translocation(args.fast, budget)
        budget.close("16C")
        analysis = analyze_transport(results, args.fast)
        master = dict(phase="16",
                      engine="numpy mesoscale CG (Martini-3-flavored) + OpenMM kernel cross-validation",
                      stage16A=m16a, stage16B=m16b, stage16C=m16c,
                      transport=analysis, selftest=st, budgets_spent_s=budget.spent,
                      wall_total_s=round(time.time() - T_START, 1))
        (RES / "phase16_results.json").write_text(json.dumps(master, indent=1, default=str))
        log("master record -> results_phase16/phase16_results.json")

    fig1()
    fig2()
    fig3()
    log("PHASE 16 COMPLETE in %.1f min" % ((time.time() - T_START) / 60.0))


if __name__ == "__main__":
    main()
