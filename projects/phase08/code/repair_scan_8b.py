#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""repair_scan_8b.py — recompute torsion-scan points that failed with
MCSCF non-convergence during the main 8B pass (points 20/40/60/160 deg),
then patch results_phase8/res_8b.json.  Uses deeper convergence relaxation
(e_convergence down to 2e-4 Eh; FD data needs ~1e-3 Eh accuracy)."""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_phase8_photochemical_dynamics import (  # noqa: E402
    EH_EV, OUT, psi4_geometry_string, read_xyz, _log, _warn)

import psi4  # noqa: E402

BASIS = "6-31g"
TARGET = [20.0, 40.0, 60.0, 160.0]


def cas_pair(els, xyz):
    mol = psi4.geometry(psi4_geometry_string(els, xyz))
    psi4.set_options({"basis": BASIS, "reference": "rhf", "scf_type": "df",
                      "e_convergence": 2e-7})
    e_rhf, wref = psi4.energy("hf", molecule=mol, return_wfn=True)
    nocc = wref.nalpha()
    nbf = wref.basisset().nbf()
    Ca_view = np.asarray(wref.Ca())
    act = [nocc - 2, nocc - 1, nocc, nocc + 1]
    docc = [k for k in range(nocc) if k not in act[:2]]
    virt = [k for k in range(Ca_view.shape[1]) if k not in act and k >= nocc]
    Ca_view[:] = Ca_view[:, docc + act + virt]
    psi4.set_options({"restricted_docc": [len(docc)], "active": [4],
                      "restricted_uocc": [nbf - len(docc) - 4]})
    psi4.set_module_options("detci", {"num_roots": 2, "avg_states": [0, 1]})
    tries = [({}, 150),
             ({"maxiter": 300, "e_convergence": 1e-6}, 150),
             ({"maxiter": 400, "e_convergence": 1e-5}, 150),
             ({"maxiter": 600, "e_convergence": 2e-4}, 150),
             ({"maxiter": 800, "e_convergence": 1e-3}, 150)]
    last = None
    for extra, it0 in tries:
        try:
            if extra:
                psi4.set_options(extra)
            _e, cw = psi4.energy("casscf", molecule=mol, return_wfn=True,
                                 ref_wfn=wref)
            e0 = float(cw.variable("CI ROOT 0 TOTAL ENERGY"))
            e1 = float(cw.variable("CI ROOT 1 TOTAL ENERGY"))
            psi4.set_options({"maxiter": it0, "e_convergence": 2e-7})
            psi4.core.clean()
            return e0, e1
        except Exception as exc:
            last = exc
            psi4.set_options({"maxiter": it0, "e_convergence": 2e-7})
            psi4.core.clean()
    raise RuntimeError(f"all tiers failed: {str(last)[:80]}")


def main():
    psi4.set_memory("2 GB")
    psi4.set_num_threads(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
    scan = json.loads((OUT / "scan_diazene.json").read_text())["scan_torsion"]
    els = None
    for key, g in sorted(scan.items(), key=lambda kv: -kv[1]["phi"]):
        if g["phi"] not in TARGET:
            continue
        els, xyzg = read_xyz(Path(g["file"]))
        try:
            e0, e1 = cas_pair(els, xyzg)
            _log("FIX", f"phi={g['phi']:.0f}: gap={(e1 - e0) * EH_EV:.3f} eV")
            (OUT / f"repair_phi{int(g['phi'])}.json").write_text(json.dumps(
                {"phi": g["phi"], "e0_eh": e0, "e1_eh": e1,
                 "gap_eV": (e1 - e0) * EH_EV}))
        except Exception as exc:
            _warn("FIX", f"phi={g['phi']:.0f} failed: {str(exc)[:90]}")


if __name__ == "__main__":
    main()
