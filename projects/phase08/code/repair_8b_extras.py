#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""repair_8b_extras.py — recompute missing 8B extras (branching-space
verification cuts, N=N stretch points) with deep convergence relaxation and
patch results_phase8/res_8b.json in place."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_phase8_photochemical_dynamics import (  # noqa: E402
    EH_EV, OUT, psi4_geometry_string, read_xyz, _log, _warn)

import psi4  # noqa: E402

BASIS = "6-31g"
CUT_T = [-0.2, -0.1, 0.0, 0.1, 0.2]
STRETCH_D = [-0.10, -0.05, 0.0, 0.05, 0.10]


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
             ({"maxiter": 500, "e_convergence": 1e-5}, 150),
             ({"maxiter": 800, "e_convergence": 5e-4}, 150)]
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
    res_path = OUT / "res_8b.json"
    res = json.loads(res_path.read_text())
    els, x_meci = read_xyz(OUT / "meci.xyz")
    x = x_meci.ravel().copy()

    g_vec = np.array([c for row in res["g_vector"]["per_atom"]
                      for c in row])
    gn = g_vec / np.linalg.norm(g_vec)
    h_vec = np.array([c for row in res["h_vector"]["per_atom"]
                      for c in row])
    hn = h_vec / np.linalg.norm(h_vec)

    cuts = res.setdefault("branching_cuts", {})
    for label, vec in (("g", gn), ("h", hn)):
        cuts.setdefault(label, [])
        done = {r["t_ang"] for r in cuts[label]}
        for t in CUT_T:
            if t in done:
                continue
            try:
                e0, e1 = cas_pair(els, (x + t * vec).reshape(-1, 3))
                cuts[label].append({"t_ang": t, "gap_eV": (e1 - e0) * EH_EV,
                                    "e0_eh": e0, "e1_eh": e1})
                _log("FIX", f"cut {label} {t:+.1f}: "
                            f"{(e1 - e0) * EH_EV:.3f} eV")
            except Exception as exc:
                _warn("FIX", f"cut {label} {t:+.1f} failed: {str(exc)[:70]}")
        cuts[label].sort(key=lambda r: r["t_ang"])

    spts = res.setdefault("stretch_scan_at_ci", [])
    done_d = {p["d"] for p in spts}
    axis = (x.reshape(-1, 3)[0] - x.reshape(-1, 3)[1])
    axis = axis / np.linalg.norm(axis)
    for d in STRETCH_D:
        if d in done_d:
            continue
        pos = x.reshape(-1, 3).copy()
        pos[0] += axis * (d / 2)
        pos[1] -= axis * (d / 2)
        try:
            e0, e1 = cas_pair(els, pos)
            spts.append({"d": d, "e0_eh": e0, "e1_eh": e1,
                         "gap_eV": (e1 - e0) * EH_EV})
            _log("FIX", f"stretch {d:+.2f}: {(e1 - e0) * EH_EV:.3f} eV")
        except Exception as exc:
            _warn("FIX", f"stretch {d:+.2f} failed: {str(exc)[:70]}")
    spts.sort(key=lambda p: p["d"])

    res_path.write_text(json.dumps(res, indent=1, default=float))
    _log("FIX", f"patched {res_path}")


if __name__ == "__main__":
    main()
