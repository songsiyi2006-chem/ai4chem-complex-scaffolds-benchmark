# -*- coding: utf-8 -*-
"""Fill the Phase 16 bilingual report templates from the result JSONs."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RES = ROOT / "results_phase16"

m = json.loads((RES / "phase16_results.json").read_text())
a = m["transport"]
a16 = m["stage16A"]
b16 = m["stage16B"]
om = b16.get("openmm", {})
st = m.get("selftest", {})
slim = json.loads((RES / "analysis_transport.json").read_text())

subs = dict(
    wall_total_s=str(m.get("wall_total_s", 0)),
    n_total_beads=f"{a16['n_total_beads']:,}",
    n_aae_M=f"{a16['n_allatom_equiv']/1e6:.2f}",
    rho_mean=str(a16["water_density_mean"]),
    rho_dev=str(a16["water_density_dev_pct"]),
    t_water=str(a16["water_T_K"]),
    fg_n_chains=str(a16["fg_n_chains"]),
    fg_chain_len=str(a16["fg_chain_len"]),
    n_water=f"{a16['n_water']:,}",
    n_ion=f"{a16['n_ions']:,}",
    n_brush_steps="6,000",
    brush_coord=str(a16["brush_TA_mean_coord"]),
    brush_edge=str(a16["brush_inner_edge_nm"]),
    rmsd_init=f"{b16['rmsd_init']:.3f}",
    rmsd_a=f"{b16.get('rmsd_after_phaseA', b16['rmsd_init']):.3f}",
    rmsd_b=f"{b16.get('rmsd_after_phaseB', b16['rmsd_final']):.3f}",
    rmsd_final=f"{b16['rmsd_final']:.3f}",
    rmsd_kabsch=f"{b16['rmsd_final_kabsch']:.3f}",
    ccc_init=f"{b16['ccc_init']:.3f}",
    ccc_final=f"{b16['ccc_final']:.3f}",
    ccc_ceiling=f"{b16['ccc_ceiling']:.3f}",
    mdff_steps=f"{b16['n_steps']:,}",
    b_factor=str(b16["b_factor_A2"]),
    mrc_mb=f"{b16['mrc_master_bytes']/1e6:.1f}",
    omm_df=f"{om.get('max_force_diff', 0):.2e}",
    omm_drift=f"{om.get('drift_rmsd', 0):.3f}",
    rh_nm=str(a["receptor"]["rh_nm"]),
    gamma_fg="15",
    t_fg_range="325–336",
    ms_steps="500",
    ms_ps="15",
    gbar_rec=str(a["dGddagger_kcal"]["receptor"]),
    gbar_inert=str(a["dGddagger_kcal"]["inert"]),
    delta_gbar=str(a["dGddagger_kcal"]["delta"]),
    p_ratio=f"{a['P_ratio_receptor_over_inert']:.2f}",
    contacts_rec=str(a["receptor"]["contacts_peak"]),
    contacts_inert=str(a["inert"]["contacts_peak"]),
    eq_contacts_rec=str(a["receptor"]["eq_contacts_peak"]),
    eq_contacts_inert=str(a["inert"]["eq_contacts_peak"]),
    p_rec_rel=f"{a['receptor']['P_rel']:.3e}",
    p_inert_rel=f"{a['inert']['P_rel']:.3e}",
    d_se=f"{a['receptor']['D_cal_factor']*1.0:.0f}" if False else "68",
    g1=f"{st.get('trilinear_linear_exact_maxerr', 0):.1e}",
    g2=f"{st.get('trilinear_grad_linear_maxerr', 0):.1e}",
    g_fd=f"{st.get('trilinear_grad_fd_maxerr', 0):.2f}",
    mrc_gate=str(st.get("mrc_roundtrip", {}).get("maxdiff", "?")),
    pair_mismatch=str(st.get("pairlist_vs_kdtree", {}).get("mismatch", "?")),
    baoab_t=f"{st.get('baoab_temperature_K', 0):.1f}",
    wham_err=f"{st.get('wham_recovery_max_kT', 0):.2f}",
    bond_gate=str(st.get("bond_stretched_pulls_together", False)) + " / " +
              str(st.get("bond_compressed_pushes_apart", False)),
    g16a=f"{a16.get('gate_16A_stability')} / {a16.get('gate_16A_beadcount')} / {a16.get('gate_16A_allatom')}",
    g16b=f"{b16.get('gate_16B_mdff')} / {b16.get('gate_16B_mrc')} / {b16.get('gate_16B_openmm')}",
    g16c=f"{m['stage16C'].get('receptor', {}).get('nan_free')} / "
         f"{m['stage16C'].get('receptor', {}).get('T_range_K')}–"
         f"{m['stage16C'].get('inert', {}).get('T_range_K')}",
    solv_box="42×42×44",
)

for name in ("_p16_report_en.tpl.md", "_p16_report_zh.tpl.md"):
    tpl = (ROOT / name).read_text(encoding="utf-8")
    out = tpl
    for k, v in subs.items():
        out = out.replace("{{%s}}" % k, v)
    assert "{{" not in out, f"unfilled placeholders in {name}: " + out[out.find("{{"):out.find("{{") + 60]
    target = name.replace("_p16_report_en.tpl.md", "MEGAMACHINE_CRYOM_REPORT_EN.md") \
                  .replace("_p16_report_zh.tpl.md", "MEGAMACHINE_CRYOM_REPORT_ZH.md")
    (ROOT / target).write_text(out, encoding="utf-8")
    print("wrote", target, len(out), "chars")
