#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
phase8_figures.py — publication-quality renderers for the Phase 8
photochemistry deliverables (300 DPI, figures_phase8/).

fig1  UV-Vis absorption spectrum + NTO particle-hole isosurfaces + azobenzene
      torsion scan of the vertical excitation energies.
fig2  S1/S0 conical-intersection topology: 3D double cone in the branching
      space (g, h), 2D gap map, and ab-initio verification cuts.
fig3  Tully-FSSH observables: S1 population decay, hopping distribution,
      E/Z product outcome statistics.

Reads results_phase8/phase8_results.json (+ fssh_population.npz, NTO cubes).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

OUT = Path("results_phase8")
FIG = Path("figures_phase8")

EH_EV = 27.211386245988
EV_NM = 1239.841984
SPEC_SIGMA = 0.2

C_MAIN = "#1F5FA8"     # deep blue
C_ACC = "#C0392B"      # red
C_MID = "#2E8B57"      # sea green
C_GOLD = "#B7950B"
C_GREY = "#5D6D7E"

plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 300, "font.size": 9.5,
    "axes.titlesize": 11, "axes.labelsize": 10.5,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "axes.axisbelow": True, "legend.framealpha": 0.9,
    "axes.unicode_minus": False,
})


def _load():
    res = json.loads((OUT / "phase8_results.json").read_text())
    return res


# --------------------------------------------------------------------------- #
#  cube rendering (phase-7 pattern)
# --------------------------------------------------------------------------- #
def _read_cube(path):
    txt = Path(path).read_text().splitlines()
    nat = int(txt[2].split()[0])
    origin = np.array([float(x) for x in txt[2].split()[1:4]])
    dims, axes = [], []
    for k in range(3):
        t = txt[3 + k].split()
        dims.append(abs(int(t[0])))
        axes.append(np.array([float(x) for x in t[1:4]]))
    vals = []
    for ln in txt[6 + nat:]:
        vals.extend(float(x) for x in ln.split())
    n = int(np.prod(dims))
    data = np.array(vals[:n]).reshape(dims)
    return origin, axes, data


def _draw_orbital(ax, cube_path, iso=0.03, title=""):
    from skimage.measure import marching_cubes
    origin, axes_v, data = _read_cube(cube_path)
    stride = 2 if max(data.shape) > 110 else 1
    if stride > 1:
        data = data[::stride, ::stride, ::stride]
    axes_v = [a * stride for a in axes_v]
    span = np.array([np.linalg.norm(a) for a in axes_v])
    if np.abs(data).max() < iso:
        iso = 0.25 * float(np.abs(data).max())
    drawn = 0
    for sign, color in ((1.0, C_MAIN), (-1.0, C_ACC)):
        try:
            verts, faces, _, _ = marching_cubes(
                (sign * data).astype(np.float32), level=iso, step_size=1)
        except ValueError:
            continue
        v = verts * span + origin
        ax.plot_trisurf(v[:, 0], v[:, 1], v[:, 2], triangles=faces,
                        color=color, alpha=0.42, linewidth=0,
                        antialiased=False, shade=True)
        drawn += 1
    ax.set_xticks([]), ax.set_yticks([]), ax.set_zticks([])
    ax.set_box_aspect((1, 1, 0.7))
    if title:
        ax.set_title(title, fontsize=9, pad=-2)
    return drawn


# --------------------------------------------------------------------------- #
#  figure 1 — UV-Vis spectrum + NTOs + torsion scan
# --------------------------------------------------------------------------- #
def _spectrum(states):
    grid = np.arange(1.5, 8.0 + 1e-9, 0.005)
    eps = np.zeros_like(grid)
    for s in states:
        if s["spin"] != "S":
            continue
        eps += s["f_osc"] * np.exp(-0.5 * ((grid - s["dE_eV"]) / SPEC_SIGMA) ** 2)
    return grid, eps


def fig1(res):
    m8a = res.get("module_8a", {})
    sing = [s for s in m8a.get("states", {}).get("singlets", [])]
    trip = [s for s in m8a.get("states", {}).get("triplets", [])]
    scan = m8a.get("torsion_scan", {}).get("points", [])

    fig = plt.figure(figsize=(13.6, 11.2))
    gs = gridspec.GridSpec(2, 3, height_ratios=[1.55, 1.0],
                           hspace=0.34, wspace=0.24)

    # ---- (a) spectrum ----------------------------------------------------- #
    ax = fig.add_subplot(gs[0, :])
    grid, eps = _spectrum(sing)
    if eps.max() > 0:
        eps_n = eps / eps.max() * 100.0
    else:
        eps_n = eps
    ax.plot(grid, eps_n, color=C_MAIN, lw=2.0,
            label=r"broadened $\varepsilon(E)$, $\sigma=0.2$ eV")
    for s in sing:
        ch = "n$\\to\\pi$*" if s.get("nto", {}).get("hole_N", 0) > 0.5 \
            else r"$\pi\to\pi$*"
        col = C_MID if ch.startswith("n") else C_ACC
        ax.vlines(s["dE_eV"], 0, s["f_osc"] / max(
            [x["f_osc"] for x in sing] + [1e-9]) * 55.0,
            color=col, lw=1.8)
        if s["root"] <= 4 or s["f_osc"] > 0.05:
            ax.annotate(
                f"S{s['root']} {ch}\n{s['dE_eV']:.2f} eV / "
                f"{EV_NM / s['dE_eV']:.0f} nm\nf={s['f_osc']:.3f}",
                xy=(s["dE_eV"], s["f_osc"] / max(
                    [x["f_osc"] for x in sing] + [1e-9]) * 55.0),
                xytext=(s["dE_eV"] - 0.05, 60 + 8 * (s["root"] % 3)),
                fontsize=7.6, color=col, ha="right",
                arrowprops=dict(arrowstyle="-", color=col, lw=0.7, alpha=0.6))
    ymax = max(eps_n.max(), 40) * 1.06
    for t in trip[:6]:
        ax.vlines(t["dE_eV"], -8, -2.2, color=C_GREY, lw=1.4)
        ax.annotate(f"T{t['root']}", xy=(t["dE_eV"], -8.6), fontsize=7,
                    color=C_GREY, ha="center")
    ax.axvspan(1.65, 3.26, color=C_GOLD, alpha=0.10, lw=0)
    ax.text(2.42, ymax * 0.86, "UVA (blue-LED\nisomerization band)",
            fontsize=8, color="#7D6608", ha="center")
    ax.set_xlim(1.5, 8.0)
    ax.set_ylim(-10, ymax)
    ax.set_xlabel("photon energy (eV)")
    ax.set_ylabel(r"normalized absorption $\varepsilon$ (%)")
    ax2 = ax.twiny()
    ax2.set_xlim(*[EV_NM / x for x in (8.0, 1.5)])
    ax2.set_xlabel("wavelength (nm)")
    ax2.grid(False)
    ax.set_title(f"Module 8A — simulated UV-Vis absorption of trans-azobenzene "
                 f"({m8a.get('method', 'TDA-DFT')})", pad=30)
    h1, l1 = [plt.Line2D([], [], color=c, lw=2) for c in (C_MAIN, C_MID, C_ACC, C_GREY)], \
             ["Gaussian-broadened spectrum",
              r"dark singlets (n$\to\pi$*)",
              r"bright singlets ($\pi\to\pi$*)",
              "triplet states (unscaled sticks)"]
    ax.legend(h1, l1, loc="upper right", fontsize=8)

    # ---- (b/c) NTO isosurfaces (hole | particle) --------------------------- #
    cubes = m8a.get("nto_cubes", {})
    label = "S1" if "S1" in cubes else (next(iter(cubes), None))
    ax_hole = fig.add_subplot(gs[1, 0], projection="3d")
    ax_part = fig.add_subplot(gs[1, 1], projection="3d")
    if label:
        info = cubes[label]
        cubes_f = sorted(Path(info["dir"]).glob("*.cube"),
                         key=lambda q: (len(q.name), q.name))
        nt = info.get("nto", {})
        sub = (f"{label} (root {info.get('root', '?')}) NTO pair — hole: "
               f"{nt.get('hole_N', 0):.0%} N / {nt.get('hole_C', 0):.0%} C, "
               f"particle: {nt.get('part_N', 0):.0%} N")
        if len(cubes_f) >= 2:
            _draw_orbital(ax_hole, cubes_f[0], title="hole NTO (dominant)")
            _draw_orbital(ax_part, cubes_f[1], title="particle NTO (dominant)")
            ax_hole.text2D(0.5, -0.14, sub, transform=ax_hole.transAxes,
                           ha="center", fontsize=8.4)
        else:
            ax_hole.text2D(0.5, 0.5, "cube files missing", ha="center",
                           transform=ax_hole.transAxes, color="gray")
    else:
        ax_hole.text2D(0.5, 0.5, "NTO cubes unavailable", ha="center",
                       transform=ax_hole.transAxes, color="gray")

    # ---- (d) torsion scan -------------------------------------------------- #
    axd = fig.add_subplot(gs[1, 2])
    if scan:
        ph = np.array([p["phi"] for p in scan])
        e0 = np.array([p["e0_eh"] for p in scan])
        e0 -= e0.max() if False else e0[np.argmax(ph)]
        axd.plot(ph, e0 * EH_EV, "-o", color=C_MAIN, lw=1.8, ms=4,
                 label="S0 (TD-DFT ground)")
        s1 = [p.get("s1_eV") for p in scan]
        s2 = [p.get("s2_eV") for p in scan]
        if all(v is not None for v in s1):
            axd.plot(ph, s1, "-s", color=C_MID, lw=1.8, ms=4,
                     label=r"S1 vertical (n$\to\pi$*)")
        if all(v is not None for v in s2):
            axd.plot(ph, s2, "-^", color=C_ACC, lw=1.8, ms=4,
                     label=r"S2 vertical ($\pi\to\pi$*)")
        axd.axvline(90.0, color=C_GREY, ls=":", lw=1)
        axd.text(90, axd.get_ylim()[0], " 90°", fontsize=8, color=C_GREY,
                 va="bottom")
        axd.set_xlabel("CNNC torsion $\\varphi$ (deg)")
        axd.invert_xaxis()
        axd.set_ylabel("relative energy (eV)")
        axd.set_title("E$\\to$Z isomerization coordinate", fontsize=10)
        axd.legend(fontsize=7.6)
    else:
        axd.text(0.5, 0.5, "torsion scan unavailable", ha="center",
                 transform=axd.transAxes, color="gray")
    fig.suptitle("Fig. 1 — Azobenzene photo-switch electronic structure "
                 "(Phase 8 / Module 8A)", fontsize=12.5, y=0.985)
    fig.savefig(FIG / "fig1_uv_vis_absorption_spectrum.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("[fig] fig1_uv_vis_absorption_spectrum.png written")


# --------------------------------------------------------------------------- #
#  figure 2 — conical intersection topology
# --------------------------------------------------------------------------- #
def fig2(res):
    m8b = res.get("module_8b", {})
    if "g_vector" not in m8b or "h_vector" not in m8b:
        print("[fig] fig2 skipped — no MECI branching-space data")
        return
    g_norm = m8b["g_vector"]["norm"]
    h_rate = m8b["h_vector"]["lift_rate_eh_per_angstrom"]
    gap0 = m8b["meci"]["gap_eV"] / EH_EV
    cuts = m8b.get("branching_cuts", {})

    fig = plt.figure(figsize=(13.6, 9.6))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.5, 1.0],
                           hspace=0.30, wspace=0.22)

    # ---- (a) 3D double cone ------------------------------------------------ #
    ax3 = fig.add_subplot(gs[0, 0], projection="3d")
    ext = 0.4
    uu = np.linspace(-ext, ext, 61)
    vv = np.linspace(-ext, ext, 61)
    UU, VV = np.meshgrid(uu, vv)
    sg = max(g_norm, 1e-6)
    kh = max(h_rate, 1e-6)
    lin = 0.5 * sg * UU
    sq = np.sqrt(lin**2 + (kh * VV) ** 2 + (gap0 / 2) ** 2)
    E0 = (lin - sq - gap0 / 2) * EH_EV
    E1 = (lin + sq + gap0 / 2) * EH_EV
    ax3.plot_surface(UU, VV, E1, cmap="Reds", alpha=0.55, linewidth=0,
                     antialiased=True, rstride=2, cstride=2)
    ax3.plot_surface(UU, VV, E0, cmap="Blues", alpha=0.55, linewidth=0,
                     antialiased=True, rstride=2, cstride=2)
    ax3.scatter([0], [0], [0], color="k", s=42, depthshade=False)
    ax3.text(0, 0, 0.12, "  MECI", fontsize=9, weight="bold")
    # radiationless decay path: S1 funnel -> hop -> S0
    tp = np.linspace(-0.32, 0.0, 40)
    linp = 0.5 * sg * tp
    sqp = np.sqrt(linp**2 + (kh * 0.16 * np.cos(3 * tp / ext * np.pi)) ** 2
                  + (gap0 / 2) ** 2)
    ax3.plot(tp, 0.16 * np.cos(3 * tp / ext * np.pi),
             (linp + sqp + gap0 / 2) * EH_EV, color=C_ACC, lw=2.2,
             label="S1 funnel path")
    ax3.plot(tp, 0.16 * np.cos(3 * tp / ext * np.pi),
             (linp - sqp - gap0 / 2) * EH_EV, color=C_MAIN, lw=2.2, ls="--",
             label="S0 recovery path")
    ax3.set_xlabel("g\u0302 branching coordinate (Å)")
    ax3.set_ylabel("h\u0302 coupling coordinate (Å)")
    ax3.set_zlabel("E (eV)")
    ax3.set_box_aspect((1, 1, 0.62))
    ax3.set_title("S1/S0 conical-intersection funnel\n"
                  "(penalty-MECI + branching-space reconstruction)", fontsize=10)
    ax3.legend(fontsize=8, loc="upper left")

    # ---- (b) 2D gap map ---------------------------------------------------- #
    axc = fig.add_subplot(gs[0, 1])
    dE = 2.0 * np.sqrt(lin**2 + (kh * VV) ** 2 + (gap0 / 2) ** 2) * EH_EV
    lv = [0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4]
    cs = axc.contourf(UU, VV, dE, levels=lv, cmap="magma", alpha=0.85)
    cl = axc.contour(UU, VV, dE, levels=lv, colors="w", linewidths=0.5)
    axc.clabel(cl, fontsize=6.5, fmt="%.2f")
    fig.colorbar(cs, ax=axc, shrink=0.9, label=r"$\Delta E_{10}$ (eV)")
    ar = axc.annotate("", xy=(0.32, 0), xytext=(-0.32, 0),
                      arrowprops=dict(arrowstyle="-|>", color=C_MAIN, lw=2))
    ah = axc.annotate("", xy=(0, 0.28 * kh / sg * 1.0), xytext=(0, 0),
                      arrowprops=dict(arrowstyle="-|>", color=C_ACC, lw=2))
    axc.text(0.16, 0.035, "g (slope dir.)", color=C_MAIN, fontsize=9,
             weight="bold")
    axc.text(-0.42, 0.30, "h (coupling dir.)", color=C_ACC, fontsize=9,
             weight="bold")
    axc.plot(0, 0, "wo", ms=6, mec="k")
    axc.set_xlabel("displacement along g\u0302 (Å)")
    axc.set_ylabel("displacement along h\u0302 (Å)")
    axc.set_title("degeneracy lifting in the branching plane\n"
                  "(log-spaced adiabatic gap contours)", fontsize=10)

    # ---- (c) verification cuts -------------------------------------------- #
    axv = fig.add_subplot(gs[1, 0])
    for label, color in (("g", C_MAIN), ("h", C_ACC)):
        row = cuts.get(label)
        if not row:
            continue
        t = np.array([r["t_ang"] for r in row])
        gp = np.array([r["gap_eV"] for r in row])
        axv.plot(t, gp, "o-", color=color, lw=1.8, ms=5,
                 label=f"CASSCF cut along {label}")
        ts = np.linspace(-0.22, 0.22, 80)
        if label == "g":
            model = np.sqrt((g_norm * ts) ** 2 + gap0 ** 2) * EH_EV
        else:
            model = 2 * np.sqrt((h_rate * ts) ** 2 + (gap0 / 2) ** 2) * EH_EV
        axv.plot(ts, model, "--", color=color, lw=1.2, alpha=0.8,
                 label=f"2-state model ({label})")
    axv.axhline(0.05, color=C_GREY, ls=":", lw=1)
    axv.text(0.16, 0.07, "mission MECI gate $\\Delta E_{10} < 0.05$ eV",
             fontsize=8, color=C_GREY)
    axv.set_xlabel("displacement (Å)")
    axv.set_ylabel(r"$\Delta E_{10}$ (eV)")
    axv.set_title("branching-space verification: ab-initio cuts vs local "
                  "2-state model", fontsize=10.5)
    axv.legend(fontsize=8.5)
    axv.set_yscale("symlog", linthresh=0.05)

    # ---- (d) BRS penalty-MECI convergence ---------------------------------- #
    axm = fig.add_subplot(gs[1, 1])
    steps = m8b.get("meci", {}).get("steps", [])
    if steps:
        ist = [st_["step"] for st_ in steps]
        gst = [st_["gap_eV"] for st_ in steps]
        gnorm = [st_["gradF_norm"] for st_ in steps]
        axm.semilogy(ist, gst, "-o", color=C_MAIN, lw=1.8, ms=4,
                     label=r"$\Delta E_{10}$ (eV)")
        axm.axhline(0.05, color=C_GREY, ls=":", lw=1)
        axm.text(steps[0]["step"], 0.058,
                 "mission gate 0.05 eV", fontsize=7.5, color=C_GREY)
        axm.set_xlabel("BRS penalty optimization step")
        axm.set_ylabel(r"$\Delta E_{10}$ (eV)")
        axm2 = axm.twinx()
        axm2.semilogy(ist, gnorm, "-s", color=C_GOLD, lw=1.2, ms=3,
                      alpha=0.8, label=r"$|\nabla F|$ (Eh/Å)")
        axm2.set_ylabel(r"$|\nabla F_{MECI}|$ (Eh/Å)", color=C_GOLD)
        axm2.grid(False)
        axm.set_title("Bearpark-Robb-Schlegel penalty convergence "
                      f"(gap {gst[0]:.2f} -> {gst[-1]:.4f} eV "
                      f"in {len(steps)} steps)", fontsize=9.5)
        l1m, lam1 = axm.get_legend_handles_labels()
        l2m, lam2 = axm2.get_legend_handles_labels()
        axm.legend(l1m + l2m, lam1 + lam2, fontsize=7.5, loc="upper right")
    else:
        axm.text(0.5, 0.5, "MECI steps unavailable", ha="center",
                 transform=axm.transAxes, color="gray")

    fig.suptitle("Fig. 2 — Conical intersection topology of the diazene "
                 "photo-switch (Phase 8 / Module 8B)", fontsize=12.5, y=0.99)
    fig.savefig(FIG / "fig2_conical_intersection_topology.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("[fig] fig2_conical_intersection_topology.png written")


# --------------------------------------------------------------------------- #
#  figure 3 — FSSH population dynamics
# --------------------------------------------------------------------------- #
def fig3(res):
    m8c = res.get("module_8c", {})
    if "p_s1_active" not in m8c:
        print("[fig] fig3 skipped — no FSSH data")
        return
    t = np.array(m8c["times_fs"])
    p_act = np.array(m8c["p_s1_active"])
    p_coh = np.array(m8c["p_s1_coherent"])
    hops = np.array(m8c["hop_histogram"])
    phi_f = np.array(m8c["phi_final_deg"])
    act_f = np.array(m8c["active_final"])
    dec_t = np.array([np.nan if v is None else v
                      for v in m8c["decay_time_fs"]])
    tau_exp = m8c.get("tau_exp_fs")
    tau_half = m8c.get("tau_half_fs")
    phi_b = m8c.get("phi_ci_deg_boundary", 90.0)

    fig = plt.figure(figsize=(13.8, 9.2))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.25, 1.0],
                           hspace=0.32, wspace=0.22)

    # ---- (a) population decay ---------------------------------------------- #
    ax = fig.add_subplot(gs[0, :])
    ax.plot(t, p_act, "-", color=C_ACC, lw=2.4,
            label=f"active-state population $P_{{S_1}}(t)$ "
                  f"(N = {m8c['n_traj']} trajectories)")
    ax.plot(t, p_coh, "--", color=C_MAIN, lw=1.6,
            label=r"coherent $\langle|c_1|^2\rangle$ (Tully amplitudes)")
    if tau_exp:
        C = m8c.get("tau_fit_C", 0.02)
        A = m8c.get("tau_fit_A", 0.95)
        ax.plot(t, A * np.exp(-t / tau_exp) + C, ":", color="k", lw=1.4,
                label=rf"exponential fit: $\tau = {tau_exp:.0f}$ fs")
    if tau_half:
        ax.axvline(tau_half, color=C_GREY, lw=1.1, ls="-.")
        ax.annotate(rf"$\tau_{{1/2}} = {tau_half:.0f}$ fs",
                    xy=(tau_half, 0.5), xytext=(tau_half + 18, 0.62),
                    fontsize=9.5,
                    arrowprops=dict(arrowstyle="-|>", color=C_GREY, lw=1))
    ax.set_xlabel("time (fs)")
    ax.set_ylabel("excited-state population")
    ax.set_title("S1 -> S0 non-adiabatic decay: ensemble-averaged FSSH "
                 "population (Tully fewest-switches, dt = 0.5 fs)", fontsize=10.5)
    ax.set_xlim(0, m8c["tmax_fs"])
    ax.set_ylim(-0.03, 1.03)
    ax.legend(fontsize=9, loc="upper right")

    # ---- (b) hop histogram -------------------------------------------------- #
    axh = fig.add_subplot(gs[1, 0])
    th = np.arange(len(hops)) * m8c["dt_fs"]
    axh.bar(th, hops, width=m8c["dt_fs"] * 1.02, color=C_MID, alpha=0.85,
            label="hops / 0.5 fs bin")
    axh.set_xlabel("time (fs)")
    axh.set_ylabel("non-adiabatic hop events")
    axh2 = axh.twinx()
    cum = np.cumsum(hops)
    axh2.plot(th, cum, color=C_GOLD, lw=2.0, label="cumulative")
    axh2.set_ylabel("cumulative hops", color=C_GOLD)
    axh2.grid(False)
    axh.set_title(f"S1->S0 hopping distribution "
                  f"({m8c['n_hops_total']} hops, "
                  f"{m8c['hops_per_traj']:.2f} / trajectory, "
                  f"{m8c['n_frustrated_hops']} frustrated)", fontsize=10)
    l1, la1 = axh.get_legend_handles_labels()
    l2, la2 = axh2.get_legend_handles_labels()
    axh.legend(l1 + l2, la1 + la2, fontsize=8, loc="upper right")

    # ---- (c) outcomes -------------------------------------------------------- #
    axo = fig.add_subplot(gs[1, 1])
    nZ, nE, nR = m8c["n_Z"], m8c["n_E"], m8c["n_S1_resident"]
    n = m8c["n_traj"]
    labels = ["Z product\n(isomerized)", "E recovered", "S1 resident"]
    vals = [nZ, nE, nR]
    cols = [C_ACC, C_MAIN, C_GREY]
    bars = axo.bar(labels, [v / n * 100 for v in vals], color=cols, alpha=0.88)
    for b, v in zip(bars, vals):
        axo.text(b.get_x() + b.get_width() / 2, b.get_height() + 1.2,
                 f"{v / n * 100:.1f}%\n({v})", ha="center", fontsize=8.6)
    axo.set_ylabel("trajectory fraction (%)")
    axo.set_ylim(0, 100)
    axo.set_title(rf"photochemical outcome at 500 fs "
                  rf"($\Phi_Z={nZ / n:.2f}$, $\Phi_E={nE / n:.2f}$)",
                  fontsize=10)
    twin = axo.inset_axes([0.04, 0.56, 0.44, 0.40])
    on_s0 = act_f == 0
    twin.scatter(dec_t[on_s0], phi_f[on_s0], s=7, c=C_MAIN, alpha=0.55,
                 label="decayed (S0)")
    twin.scatter(dec_t[~on_s0], phi_f[~on_s0], s=7, c=C_ACC, alpha=0.55,
                 marker="^", label="S1 resident")
    twin.axhline(phi_b, color=C_GREY, ls=":", lw=1)
    twin.text(0.03, phi_b, f"{phi_b:.0f}° CI basin boundary", fontsize=6,
              color=C_GREY, va="bottom", transform=twin.get_yaxis_transform())
    twin.set_xlabel("decay time (fs)", fontsize=7)
    twin.set_ylabel(r"$\varphi_{final}$ (°)", fontsize=7)
    twin.tick_params(labelsize=6.5)
    twin.legend(fontsize=6, loc="lower right")

    fig.suptitle("Fig. 3 — Tully surface-hopping photodynamics on the "
                 "ab-initio-parameterized LVC model (Phase 8 / Module 8C)",
                 fontsize=12.5, y=0.99)
    fig.savefig(FIG / "fig3_fssh_population_trajectories.png", dpi=300,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("[fig] fig3_fssh_population_trajectories.png written")


def render_all():
    FIG.mkdir(exist_ok=True)
    res = _load()
    fig1(res)
    fig2(res)
    fig3(res)


if __name__ == "__main__":
    render_all()
