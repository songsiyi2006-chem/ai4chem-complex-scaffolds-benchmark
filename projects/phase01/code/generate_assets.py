#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_assets.py
==================
Publication-quality figure generation for the 10-molecule hardcore benchmark.

Reads the machine-readable records produced by molecule_benchmark.py
(bench_results/benchmark_results.json) and renders three 300-DPI figures
into ./figures/:

  fig1_molecular_grid.png    2D structural grid (all targets + M09R),
                             annotated with ID / MW / MaxRing
  fig2_chemical_space.png    MW vs cLogP bubble chart (size = TPSA,
                             color = Fsp3) with the Lipinski Ro5 reference
                             zone
  fig3_radar_complexity.png  normalized complexity radar over
                             [RotB, Fsp3, Stereocenters, MaxRing, TPSA,
                             cLogP] for M01 / M02 / M03 / M05

Usage:
    python generate_assets.py                 # default paths
    python generate_assets.py --results bench_results/benchmark_results.json --out figures
"""

from __future__ import annotations

import argparse
import io
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

if sys.platform.startswith("win"):
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

sns.set_theme(style="whitegrid", context="notebook", font_scale=1.05,
              rc={"axes.unicode_minus": False, "figure.dpi": 110})

SHORT_NAMES = {
    "M01": "Macrocyclic chameleon peptide",
    "M02": "Azaspiro-cubane bioisostere",
    "M03": "PROTAC prototype (2 fragments)",
    "M04": "Atropisomeric biaryl",
    "M05": "Perfluorinated cage",
    "M06": "Oxetane polyketide mimic",
    "M07": "B-N dative macrocycle",
    "M08": "Bicyclo-acrylamide warhead",
    "M09": "Hetero-[5]-helicene (INVALID)",
    "M09R": "Aza-[5]-helicene (repaired)",
    "M10": "Tetra-ortho peptoid",
}

RADAR_IDS = ["M01", "M02", "M03", "M05"]
RADAR_AXES = ["RotB", "Fsp\u00b3", "Stereocenters", "MaxRing", "TPSA (\u00c5\u00b2)", "cLogP"]


def load_records(path: Path) -> Dict[str, Any]:
    """Load benchmark records; map id -> record."""
    if not path.exists():
        raise SystemExit(f"results file not found: {path}\n"
                         f"run molecule_benchmark.py first to produce it.")
    payload = json.loads(path.read_text(encoding="utf-8"))
    return {rec["id"]: rec for rec in payload["results"] if "id" in rec}


def _props(rec: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Return property dict if the record parsed successfully."""
    return rec.get("props") if rec.get("status") != "failed_parse" else None


# --------------------------------------------------------------------------- #
# figure 1: structural grid
# --------------------------------------------------------------------------- #
def figure_molecular_grid(records: Dict[str, Dict[str, Any]], out_dir: Path) -> Path:
    """2D depictions of every entry (M09 shown as a failure placeholder)."""
    from rdkit import Chem
    from rdkit.Chem.Draw import MolToImage

    ids = ["M01", "M02", "M03", "M04", "M05", "M06", "M07", "M08",
           "M09", "M09R", "M10"]
    ids = [i for i in ids if i in records]
    ncols = 4
    nrows = math.ceil(len(ids) / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.1 * ncols, 3.7 * nrows))
    for ax in axes.flat:
        ax.set_axis_off()

    for ax, mid in zip(axes.flat, ids):
        rec = records[mid]
        props = _props(rec)
        if props is None:
            ax.set_title(f"{mid} \u00b7 {SHORT_NAMES.get(mid, mid)}",
                         fontsize=10, fontweight="bold", color="firebrick")
            ax.text(0.5, 0.48,
                    "SMILES unkekulizable\n(kekulization failure \u2014\nsee M09R reference)",
                    ha="center", va="center", fontsize=11, color="firebrick",
                    transform=ax.transAxes)
            ax.text(0.5, -0.06, "MW \u2014 | MaxRing \u2014", ha="center",
                    fontsize=9, color="0.35", transform=ax.transAxes)
            continue

        mol = Chem.MolFromSmiles(rec["smiles"])
        if mol is None:  # defensive; cannot happen for status != failed_parse
            continue
        img = MolToImage(mol, size=(560, 420))
        ax.imshow(img)
        ax.set_title(f"{mid} \u00b7 {SHORT_NAMES.get(mid, rec.get('name', mid))}",
                     fontsize=10, fontweight="bold")
        ax.text(0.5, -0.05,
                f"MW {props['mw_full_input']:.1f} Da | "
                f"MaxRing {rec['graph']['max_ring']} | "
                f"E$_{{min}}$ {rec['conformers']['e_min']:.1f} kcal/mol",
                ha="center", fontsize=9, color="0.30", transform=ax.transAxes)

    for ax in axes.flat[len(ids):]:
        ax.set_visible(False)

    fig.suptitle("Benchmark Set: 10 Structurally Complex Molecular Entities (+1 repaired reference)",
                 fontsize=13, fontweight="bold", y=0.99)
    # manual spacing: footer texts live below each axes, so reserve generous gaps
    fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.02,
                        wspace=0.12, hspace=0.55)
    out = out_dir / "fig1_molecular_grid.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
# figure 2: chemical-space bubble plot
# --------------------------------------------------------------------------- #
def figure_chemical_space(records: Dict[str, Dict[str, Any]], out_dir: Path) -> Path:
    """MW vs cLogP bubbles; size = TPSA, color = Fsp3; Ro5 reference zone."""
    pts = []
    for mid, rec in records.items():
        p = _props(rec)
        if p is None:
            continue
        pts.append({
            "id": mid, "mw": p["mw_full_input"], "clogp": p["clogp"],
            "tpsa": p["tpsa"], "fsp3": p["fsp3"],
            "flags": rec.get("risk_flags", []),
        })

    fig, ax = plt.subplots(figsize=(9.5, 7.2))
    size_scale = 3.0

    # Lipinski Ro5 reference zone: MW <= 500 Da AND cLogP <= 5
    ax.add_patch(Rectangle((180, -2), 320, 7, facecolor="#2e7d32", alpha=0.07,
                           edgecolor="#2e7d32", linestyle="--", linewidth=1.2))
    ax.text(343, 4.6, "Lipinski Ro5 zone\n(MW \u2264 500, cLogP \u2264 5)",
            ha="center", va="top", fontsize=9, color="#2e7d32")

    sc = ax.scatter([p["mw"] for p in pts], [p["clogp"] for p in pts],
                    s=[p["tpsa"] * size_scale for p in pts],
                    c=[p["fsp3"] for p in pts], cmap="plasma",
                    vmin=0.0, vmax=0.9, alpha=0.85,
                    edgecolors="black", linewidths=0.8, zorder=3)

    for p in pts:
        dx, dy = {"M07": (-14, -16)}.get(p["id"], (6, 7))
        ax.annotate(p["id"], (p["mw"], p["clogp"]), xytext=(dx, dy),
                    textcoords="offset points", fontsize=9.5, fontweight="bold",
                    color="0.15", zorder=4)

    # bubble-size reference (TPSA) — anchored in the data-free lower-right
    size_handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=math.sqrt(t * size_scale),
               markerfacecolor="0.75", markeredgecolor="0.35", alpha=0.9,
               label=f"TPSA = {t} \u00c5\u00b2")
        for t in (50, 150, 250)
    ]
    ax.legend(handles=size_handles, title="Bubble size (TPSA)", loc="lower right",
              bbox_to_anchor=(0.985, 0.05), frameon=True, fontsize=9,
              title_fontsize=9, labelspacing=1.3, borderpad=0.8)

    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("Fraction sp\u00b3 carbons (Fsp\u00b3)", fontsize=10.5)

    ax.set_xlabel("Molecular weight (Da)", fontsize=11)
    ax.set_ylabel("cLogP (Crippen)", fontsize=11)
    ax.set_title("Chemical Space of the Benchmark Set vs. Conventional Drug-Like Space",
                 fontsize=13, fontweight="bold", pad=12)
    ax.set_xlim(180, 900)
    ax.set_ylim(-1.6, 6.4)
    ax.grid(True, alpha=0.35)

    fig.tight_layout()
    out = out_dir / "fig2_chemical_space.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
# figure 3: complexity radar
# --------------------------------------------------------------------------- #
def _radar_raw(rec: Dict[str, Any]) -> List[float]:
    p, g = rec["props"], rec["graph"]
    st = rec["stereo"]
    return [p["rotatable_bonds"], p["fsp3"], st["assigned"] + st["unassigned"],
            g["max_ring"], p["tpsa"], p["clogp"]]


def figure_radar(records: Dict[str, Dict[str, Any]], out_dir: Path) -> Path:
    """Normalized complexity radar for representative scaffolds."""
    valid = [r for r in records.values() if _props(r) is not None]
    maxima = np.max(np.array([_radar_raw(r) for r in valid], dtype=float), axis=0)
    maxima[maxima == 0] = 1.0

    n = len(RADAR_AXES)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist()
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(8.6, 7.4),
                           subplot_kw=dict(polar=True))
    palette = {"M01": "#0173B2", "M02": "#D55E00", "M03": "#009E73", "M05": "#CC78BC"}
    for mid in RADAR_IDS:
        vals = (np.array(_radar_raw(records[mid]), dtype=float) / maxima).tolist()
        vals += vals[:1]
        ax.plot(angles, vals, linewidth=2.0, label=f"{mid} \u00b7 {SHORT_NAMES[mid]}",
                color=palette[mid])
        ax.fill(angles, vals, alpha=0.12, color=palette[mid])

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(RADAR_AXES, fontsize=11)
    ax.set_rgrids([0.25, 0.50, 0.75, 1.00], angle=90,
                  labels=["0.25", "0.50", "0.75", "1.00"], fontsize=8, color="0.45")
    ax.set_ylim(0, 1.02)
    ax.grid(alpha=0.45)
    ax.spines["polar"].set_color("0.6")

    ax.set_title("Structural Complexity Fingerprint of Representative Scaffolds\n"
                 "(each axis normalized to the maximum across the benchmark set)",
                 fontsize=12.5, fontweight="bold", pad=24)
    ax.legend(loc="upper right", bbox_to_anchor=(1.32, 1.12), fontsize=9.5,
              frameon=True)

    fig.tight_layout()
    out = out_dir / "fig3_radar_complexity.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="generate benchmark figures")
    parser.add_argument("--results", default="bench_results/benchmark_results.json")
    parser.add_argument("--out", default="figures")
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    records = load_records(Path(args.results))

    for fn in (figure_molecular_grid, figure_chemical_space, figure_radar):
        path = fn(records, out_dir)
        print(f"[assets] wrote {path} ({path.stat().st_size / 1024:.0f} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
