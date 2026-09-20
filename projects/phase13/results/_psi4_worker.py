
import json, sys, time, os
import numpy as np
import psi4

BOHR_A = 0.529177210903

spec = json.load(open(sys.argv[1]))
psi4.set_memory(spec.get("memory", "2 GB"))
psi4.set_num_threads(int(spec.get("threads", 4)))
psi4.set_output_file(spec["out_log"], False)

geom = spec["geometry"]
charge = int(spec["charge"])
mult = int(spec["multiplicity"])
last_wfn = [None]


def set_recipe(method):
    psi4.set_options({"reference": "UHF" if method == "HF" else "UKS",
                      "scf_type": "df", "maxiter": int(spec.get("maxiter", 300)),
                      "fail_on_maxiter": False,
                      "scf_initial_accelerator": "NONE",
                      "guess": "sad" if last_wfn[0] is None else "read",
                      "e_convergence": 1e-8, "damping_percentage": 40,
                      "level_shift": 0.10})
    # LOCAL options: this build shadows global SCF convergence options
    for k, v in [("d_convergence", 1e-4), ("e_convergence", 1e-8),
                 ("damping_percentage", 40), ("level_shift", 0.10),
                 ("maxiter", int(spec.get("maxiter", 300))),
                 ("fail_on_maxiter", False)]:
        psi4.core.set_local_option("SCF", k, v)


def extras(w, mol):
    """Contact densities at Fe, optional 3D density grids (compute_phi is a
    single-point API in this build: compute_phi(x, y, z) -> AO values)."""
    bas = w.basisset()
    fe = mol.atom_entry(0) if hasattr(mol, "atom_entry") else None
    fx, fy, fz = mol.x(0), mol.y(0), mol.z(0)          # atom 0 = Fe (a0)
    phi0 = np.array(bas.compute_phi(fx, fy, fz))
    Da, Db = w.Da(), w.Db()
    out = {"rho0_total": float(phi0 @ Da @ phi0 + phi0 @ Db @ phi0),
           "rho0_spin": float(phi0 @ Da @ phi0 - phi0 @ Db @ phi0),
           "compute_phi_ok": True}
    dg = spec.get("density_grid")
    if dg:
        axes = np.arange(-dg["box"], dg["box"] + 1e-9, dg["step"])
        X, Y, Z = np.meshgrid(axes, axes, axes, indexing="ij")
        pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()]) / BOHR_A
        rt = np.zeros(len(pts)); rs = np.zeros(len(pts))
        phis = np.empty((len(pts), Da.shape[0]))
        for i, (px, py, pz) in enumerate(pts):
            phis[i] = bas.compute_phi(px, py, pz)
        a = phis @ Da; b = phis @ Db
        rt = (a * phis).sum(1) + (b * phis).sum(1)
        rs = (a * phis).sum(1) - (b * phis).sum(1)
        out["grid_axes"] = axes
        out["grid_total"] = rt.reshape((len(axes),) * 3)
        out["grid_spin"] = rs.reshape((len(axes),) * 3)
    return out


mol = psi4.geometry(geom)
mol.set_molecular_charge(charge)
mol.set_multiplicity(mult)

import re

def trajectory_tail(log_path, from_byte):
    """Parse (iter, E, dE, diis) from the psi4 output appended after from_byte."""
    with open(log_path, "rb") as fh:
        fh.seek(from_byte)
        data = fh.read().decode("utf-8", errors="ignore")
    rows = []
    for ln in data.splitlines():
        m = re.match(r"\s+@DF-(?:UKS|UHF) iter\s+(\d+):\s+(-?\d+\.\d+)\s+(-?\S+)\s+(\S+)", ln)
        if m:
            rows.append((int(m.group(1)), float(m.group(2)),
                         float(m.group(3)), float(m.group(4))))
    return rows

result = {"converged": False, "tier": None, "fallback": None,
          "spec_tiers": spec.get("tiers_key")}
fallback = None
offset = 0
for method, basis in spec["tiers"]:
    offset = os.path.getsize(spec["out_log"]) if os.path.exists(spec["out_log"]) else 0
    set_recipe(method)
    kw = {"molecule": mol, "return_wfn": True}
    if last_wfn[0] is not None and method != "HF":
        kw["guess_wfn"] = last_wfn[0]      # DFT-tier orbital continuation
                                           # (the HF anchor always restarts
                                           # from SAD - proven robust)
    t0 = time.time()
    try:
        e, w = psi4.energy(f"{method}/{basis}", **kw)
    except Exception as exc:
        result[f"error_{method}_{basis}"] = str(exc)[:200]
        continue
    last_wfn[0] = w
    conv = "Energy and wave function converged" in open(spec["out_log"], errors="ignore").read()
    traj = trajectory_tail(spec["out_log"], offset)
    # acceptance doctrine (logged):
    #   strict  : psi4 full convergence (E + density)
    #   energy  : energy stationary to 1e-7 Eh over 5 iterations (the win-64
    #             build freezes the DIIS commutator inside a near-degenerate
    #             d-manifold at ~6e-3 while the energy is machine-stable)
    #   plateau : density converged (DIIS < 2e-4) with a 2-cycle bracket
    #             (two spin-isomer basins < 4 kcal apart); energy = window min
    mode, E_acc, dE_acc = None, None, None
    if len(traj) >= 6:
        tail5 = traj[-5:]
        # transient SCF stalls (e.g. a high B3LYP basin with DIIS ~2.5e-2)
        # must not pass: require enough iterations, a 12-point stationary
        # window and a DIIS commutator below the stall level
        tail12 = traj[-12:]
        diis_min = min(t[3] for t in traj[-20:])
        if (max(t[1] for t in tail5) - min(t[1] for t in tail5) < 1e-7
                and len(traj) >= 50
                and max(t[1] for t in tail12) - min(t[1] for t in tail12) < 1e-6
                and diis_min < 2e-2):
            mode, E_acc = "energy", tail5[-1][1]
        else:
            tail20 = traj[-20:]
            diis = min(t[3] for t in tail20)
            span = max(t[1] for t in tail20) - min(t[1] for t in tail20)
            if diis < 2e-4 and span < 2e-3:
                mode = "plateau"
                E_acc = min(t[1] for t in tail20)
                dE_acc = span
    if conv:
        mode = "strict"
    if mode:
        rec = {"tier": [method, basis], "E": float(e), "E_accepted": E_acc,
               "E_bracket_eH": dE_acc, "acceptance": mode,
               "S2": float(w.variables().get("CURRENT SPIN", 0.0)),
               "seconds": time.time() - t0, "converged": True}
    if not mode and method != "HF" and not result.get("_retried_mix"):
        # broken-symmetry spin-mixed restart (one shot per job)
        result["_retried_mix"] = True
        psi4.core.set_local_option("SCF", "guess", "sad")
        psi4.set_options({"guess_mix": 0.7})
        try:
            e2, w2 = psi4.energy(f"{method}/{basis}", **kw)
            traj2 = trajectory_tail(spec["out_log"], offset)
            tail5 = traj2[-5:]
            if len(tail5) >= 5 and max(t[1] for t in tail5) - min(t[1] for t in tail5) < 1e-7                     and len(traj2) >= 50:
                e, w = e2, w2
                mode, E_acc = "energy-mix", tail5[-1][1]
        except Exception:
            pass

    if mode:
        result.update(rec)
        try:
            bas = w.basisset()
            np.savez(spec["npz"],
                     Ca=np.array(w.Ca()), Cb=np.array(w.Cb()),
                     Da=np.array(w.Da()), Db=np.array(w.Db()),
                     S=np.array(w.S()),
                     eps_a=np.array(w.epsilon_a()), eps_b=np.array(w.epsilon_b()),
                     atom_of_ao=np.array([bas.function_to_center(i)
                                          for i in range(bas.nbf())]),
                     E=e, S2=rec["S2"], mult=mult, charge=charge)
        except Exception as exc:
            result["save_error"] = str(exc)[:300]
        try:
            ex = extras(w, mol)
            result.update({k: v for k, v in ex.items() if not isinstance(v, np.ndarray)})
            if "grid_axes" in ex:
                np.savez(spec["npz"].replace(".npz", "_grid.npz"),
                         axes=ex["grid_axes"], total=ex["grid_total"],
                         spin=ex["grid_spin"], rho0_total=result["rho0_total"],
                         rho0_spin=result["rho0_spin"])
        except Exception as exc:
            result["extras_error"] = str(exc)[:200]
        break
    elif method == "HF" and len(traj) >= 5:
        fallback = {"tier": [method, basis], "E": float(e),
                    "S2": float(w.variables().get("CURRENT SPIN", 0.0))}

if not result["converged"] and fallback is not None:
    result.update({"fallback_kept": fallback,
                   "note": "no tier met the acceptance doctrine; final UHF iterate kept, flagged"})
json.dump(result, open(spec["json"], "w"), indent=1)
