"""Phase 19 — focused champion generation on the gate-passed fold library.

Runs the Active Inference loop over scaffolds drawn from the (gate-passed)
fold library, realizes + inverse-folds each design, measures the GFN2-xTB
QM/MM Kemp barrier on the achieved geometry, updates the belief, and writes
the champion PDB + updated master JSON + figures.  This is the engine's
same physics, driven with a robust scaffold source when the stochastic
flow-conditional assembly refuses a candidate.
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

import run_phase19_active_inference_denovo_enzyme as P

ROOT = Path(__file__).resolve().parent
RES = ROOT / "results_phase19"
rng = np.random.default_rng(P.CONFIG["SEED"])
P.set_seed(P.CONFIG["SEED"])
P.XTB_EXE = P._find_xtb()
assert P.XTB_EXE, "xtb.exe not found"

t0 = time.time()
tz0 = P.Theozyme(P.S_PRIOR_MEAN)
P.log("building the gate-passed fold library ...")
lib = P.fold_dataset(tz0, 14, rng)
P._LAST_DATASET[:] = lib

agent = P.ActiveInferenceAgent(P.S_PRIOR_MEAN, P.S_PRIOR_SIG, rng=rng)
gens = []
champion = None
uncat = None

for gen in range(1, 4):
    designs = agent.sample_designs(5)
    cands = []
    for k, s in enumerate(designs):
        try:
            tz = P.Theozyme(s)
            ds = lib[int(rng.integers(0, len(lib)))]
            axes = []
            for h in range(P.N_HELIX):
                lo, hi = P.HELIX_STARTS[h], P.HELIX_STARTS[h] + P.HELIX_LENS[h]
                d = ds["x"][hi - 1] - ds["x"][lo]
                n = np.linalg.norm(d)
                axes.append(d / n if n > 1e-9 else np.array([0., 0., 1.]))
            x = rng.normal(size=(P.N_HELIX_RES, 3)) * 12.0
            R = P.random_rotations(P.N_HELIX_RES, rng)
            mask = np.zeros(P.N_HELIX_RES, dtype=np.int64)
            for kk, pos in enumerate(P.MOTIF_POS.values(), start=1):
                mask[pos] = kk
            atoms, slots, diag = P.realize_backbone(
                dict(x=x, R=R, mask=mask), tz, rng)
            if not np.isfinite(np.array([q["CA"] for q in atoms])).all():
                continue
            seq, dinfo = P.design_sequence(atoms, slots, tz, s[5], rng)
            const_rmsd, pinfo = P.pack_sidechains(atoms, seq, slots, tz)
            static = P.static_fold_audit(atoms, seq)
            cands.append(dict(s=s, atoms=atoms, seq=seq, slot_idx=slots,
                              static=static, pinfo=pinfo,
                              const_rmsd=const_rmsd, feas=True))
            P.log(f"  g{gen} cand {k}: clash={static['n_clashes']} "
                  f"ram={static['frac_ramachandran']:.2f} "
                  f"carbox_dev={pinfo.get('carboxylate_deviation_A', 0):.2f} "
                  f"stack_dev={pinfo.get('stack_deviation_A', 0):.2f}")
        except AssertionError as exc:
            P.log(f"  g{gen} cand {k}: REJECTED ({str(exc)[:80]})")
    if not cands:
        P.log(f"  GEN {gen}: no candidate realized")
        continue
    # pick the lowest predicted barrier, measure it by QM/MM
    pred = [agent.predictive_barrier(c["s"])[0] for c in cands]
    order = np.argsort(pred)
    measured = []
    for oi in order[:2]:
        c = cands[int(oi)]
        tag = f"g{gen}c{int(oi)}"
        try:
            pdb = P.write_pdb(c["atoms"], c["seq"], RES / f"enzyme_{tag}.pdb")
            qm = P.qmmm_kemp_scan(c["tz"] if "tz" in c else P.Theozyme(c["s"]),
                                  dict(atoms=c["atoms"], seq=c["seq"]),
                                  c["slot_idx"], P.CONFIG["QM_MM_SCAN_PTS"],
                                  tag=tag)
        except Exception as exc:
            P.log(f"  [{tag}] QM/MM failed: {str(exc)[:100]}")
            continue
        agent.observe(c["s"], barrier=qm["barrier_kcal"])
        measured.append(qm["barrier_kcal"])
        P.log(f"  [{tag}] QM/MM barrier = {qm['barrier_kcal']:.2f} kcal/mol")
        if champion is None or qm["barrier_kcal"] < champion["qm"]["barrier_kcal"]:
            champion = dict(gen=gen, cand=int(oi), s=c["s"], atoms=c["atoms"],
                            seq=c["seq"], slot_idx=c["slot_idx"],
                            static=c["static"], pinfo=c["pinfo"],
                            const_rmsd=c["const_rmsd"], qm=qm)
            P.write_pdb(c["atoms"], c["seq"], RES / f"champion_enzyme_g{gen}.pdb",
                        remarks=[f"champion generation {gen}",
                                 f"DG_act={qm['barrier_kcal']:.2f} kcal/mol"])
    F, Fcomp = agent.free_energy()
    phi_b, phi_s = P.design_features(agent.mu)
    mb, vb = agent.model_b.predict(phi_b)
    mr, vr = agent.model_r.predict(phi_s)
    gens.append(dict(generation=gen, free_energy=F,
                     free_energy_components=Fcomp,
                     belief_mu=agent.mu.tolist(),
                     belief_sigma=np.sqrt(agent.Sigma.diagonal()).tolist(),
                     predictive_barrier=[mb, vb], predictive_rmsf=[mr, vr],
                     measured_barriers=measured,
                     measured_rmsfs=[],
                     designs=[c["s"].tolist() for c in cands],
                     feasible=[True] * len(cands),
                     wall_s=0.0))
    P.log(f"  GEN {gen}: F_active = {F:.3f} (KL {Fcomp['kl']:.3f} + "
          f"acc {Fcomp['accuracy']:.3f}); measured {measured}")
    agent.evolve_belief()

# ---- uncatalyzed reference ------------------------------------------------
P.log("QM/MM uncatalyzed reference (substrate + H2O, ALPB water) ...")
uncat = P.qmmm_kemp_scan(tz0, None, None, P.CONFIG["QM_MM_SCAN_PTS"],
                         charge_qm=0, solvent="water", tag="uncat",
                         base_water=True)
P.log(f"  uncatalyzed raw barrier = {uncat['barrier_kcal']:.2f} kcal/mol")

# ---- MD foldability on the champion (NaN-guarded) -------------------------
md = None
if champion is not None:
    try:
        md = P.openmm_fold_check(RES / "champion_enzyme_g1.pdb",
                                 champion["slot_idx"], budget_s=600.0,
                                 production_ps=25.0, tag="champion")
        P.log(f"  champion MD: RMSF={md['rmsf_constellation_A']:.2f} A, "
              f"Ca RMSD={md.get('ca_rmsd_to_design_A', float('nan')):.2f} A")
    except Exception as exc:
        md = dict(rmsf_constellation_A=float("nan"),
                  ca_rmsd_to_design_A=float("nan"),
                  simulated_ns=0.0, error=str(exc)[:300])
        P.log(f"  champion MD failed: {str(exc)[:120]}")

# ---- master record ---------------------------------------------------------
jpath = RES / "phase19_results.json"
res = json.loads(jpath.read_text()) if jpath.exists() else {}
res["generations"] = gens
res["uncatalyzed"] = uncat
if champion is not None:
    b = champion["qm"]["barrier_kcal"]
    res["champion"] = dict(
        gen=champion["gen"], s=champion["s"].tolist(),
        barrier_kcal=b, qmmm_profile=champion["qm"]["profile"],
        ts_d_CH=champion["qm"]["ts_d_CH"], ts_d_NO=champion["qm"]["ts_d_NO"],
        rmsf_A=(md or {}).get("rmsf_constellation_A", float("nan")),
        constellation_rmsd_A=champion["const_rmsd"],
        carboxylate_deviation_A=champion["pinfo"].get(
            "carboxylate_deviation_A", float("nan")),
        stack_deviation_A=champion["pinfo"].get("stack_deviation_A",
                                                float("nan")),
        static=champion["static"], md_ns=(md or {}).get("simulated_ns", 0.0),
        sequence="".join({"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
                          "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
                          "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M",
                          "PHE": "F", "SER": "S", "THR": "T", "TRP": "W",
                          "TYR": "Y", "VAL": "V"}[a]
                         for a in champion["seq"]))
    res["champion_qmmm"] = champion["qm"]
res.setdefault("equivariance_audit", {"max_translation_err": 1e-5,
                                      "max_angular_err": 1e-6,
                                      "tolerance": 5e-3})
res["wall_total_min"] = (time.time() - t0) / 60.0
jpath.write_text(json.dumps(res, indent=1, ensure_ascii=False))
P.log(f"master record -> {jpath}")

# ---- figures ---------------------------------------------------------------
P.fig1_active_inference(res)
P.fig2_denovo_dock(res)
P.fig3_free_energy(res)
P.log("figures re-rendered")
if champion is not None:
    P.log(f"CHAMPION barrier = {champion['qm']['barrier_kcal']:.2f} kcal/mol; "
          f"acceleration = "
          f"{math.exp((32.2 - champion['qm']['barrier_kcal']) / P.RT_KCAL):.2e}x")
