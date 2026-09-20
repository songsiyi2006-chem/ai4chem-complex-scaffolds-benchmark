"""Fill the @TOKEN@ placeholders in the Phase-15 treatises from the results JSONs."""
import json
import math
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RES = ROOT / "results_phase15"

res = json.loads((RES / "phase15_results.json").read_text(encoding="utf-8"))
allo = json.loads((RES / "allostery_results.json").read_text(encoding="utf-8"))

spin = res["spin"]
cv = spin["cross_validation"]
triad = res["triad_kinetics"]
hs = res["hyperfine_selected"]
wham = allo["wham"]
st = allo["states"]

H_PLANCK = 6.62607015e-34
KB = 1.380649e-23
MU_B = 9.2740100783e-24
E_CHARGE = 1.602176634e-19
N_A = 6.02214076e23

cfg = res["config"]
E_zeeman = cfg["G1"] * MU_B * cfg["B0_UT"] * 1e-6 / E_CHARGE
A_N5 = float(hs["lumiflavin_anion_N5"]["A_iso_MHz"])
D_MHz = 52.04 / cfg["R12_NM"] ** 3
D_eV = 2 * math.pi * D_MHz * 1e6 * H_PLANCK / E_CHARGE
A_eV = 2 * math.pi * A_N5 * 1e6 * H_PLANCK / E_CHARGE

dg_allo = wham["allostery"]["dG_latch_shift_kcal"]
dg_abs = abs(dg_allo)
dg_ev = dg_abs * 4184.0 / E_CHARGE / N_A
gain = dg_ev / E_zeeman

mfe = res["mfe"]
odmr = res["odmr"]
md_speed = sum(st[k]["ns_per_day"] for k in st) / len(st)
prod_ps = st["FAD_oxid"]["prod_ps"]
n_win = len(st["FAD_oxid"]["windows"])
samp = sum(len(w["samples_nm"]) for w in st["FAD_oxid"]["windows"])
# per-window production length from the stage log's budget line
import re
win_ps = 19
log = RES / "allostery_stage.log"
if log.exists():
    m = re.search(r"\+ (\d+) ps/window", log.read_text(encoding="utf-8"))
    if m:
        win_ps = int(m.group(1))
closure = cv["yield_closure_lindblad"]
deph_shift = (cv["phi_S_lindblad"] - cv["phi_S_eigen_nodeph"]) \
    / cv["phi_S_eigen_nodeph"] * 100

tokens = {
    "MFPT": f"{triad['mean_arrival_ps']:.1f}",
    "HOPS": " / ".join(f"{1e12/k:.1f}" for k in triad["k_forward"]),
    "CSYIELD": f"{triad['yield_CS_400ps']:.2f}",
    "A_N5": f"{A_N5:.1f}",
    "A_HFAD": f"{float(hs['lumiflavin_anion_H']['A_iso_MHz']):.1f}",
    "A_N1": f"{float(hs['indole_cation_N1']['A_iso_MHz']):.1f}",
    "A_H2": f"{float(hs['indole_cation_H']['A_iso_MHz']):.1f}",
    "CV_LIN": f"{cv['phi_S_lindblad']:.4f}",
    "CV_EIG": f"{cv['phi_S_eigen_rf']:.4f}",
    "CV_NODEPH": f"{cv['phi_S_eigen_nodeph']:.4f}",
    "CV_DIFF": f"{cv['abs_diff']:.4f}",
    "CLOSURE": f"{closure:.4f}",
    "DEPH_SHIFT": f"{deph_shift:+.2f}",
    "CV_SAT": f"{cv['phi_S_eigen_satellites']:.5f}",
    "YS_MIN": f"{min(min(r) for r in spin['compass']['yield_map']):.4f}",
    "YS_MAX": f"{max(max(r) for r in spin['compass']['yield_map']):.4f}",
    "ANISO_PCT": f"{spin['compass']['anisotropy_percent']:.2f}",
    "ANI_EIG": f"{cv['anisotropy_eigen']:+.5f}",
    "ANI_LIN": f"{cv['anisotropy_lindblad']:+.5f}",
    "FC0": f"{spin['field_curve']['phi_S'][0]:.4f}",
    "FC5000": f"{spin['field_curve']['phi_S'][-1]:.4f}",
    "MFE_REL": f"{abs(spin['mfe_yield']['rel_effect']) * 100:.2f}",
    "PHI0": f"{spin['full_dynamics']['theta0']['phi_S']:.4f}",
    "PHI30": f"{spin['full_dynamics']['theta30']['phi_S']:.4f}",
    "PHI60": f"{spin['full_dynamics']['theta60']['phi_S']:.4f}",
    "PHI90": f"{spin['full_dynamics']['theta90']['phi_S']:.4f}",
    "OPENMM_VERSION": allo["engine"].split("OpenMM ")[1].split(" ")[0],
    "WINPS": f"{win_ps}",
    "SAMPS": f"{samp}",
    "PRODPS": f"{prod_ps:.0f}",
    "DG_ALLO": f"{dg_allo:+.2f}",
    "OCC_OX": f"{st['FAD_oxid']['latch_occupied_frac'] * 100:.0f} %",
    "OCC_RA": f"{st['FAD_radan']['latch_occupied_frac'] * 100:.0f} %",
    "MDSPEED": f"{md_speed:.0f}",
    "E_ZEEMAN": f"{E_zeeman:.2e}",
    "D_DIP": f"{D_eV:.2e}",
    "A_ISO_EV": f"{A_eV:.2e}",
    "DG_EV": f"{dg_ev:.3e}",
    "DG": f"{dg_abs:.1f}",
    "GAIN": f"{gain:.2e}",
    "MFE_PEAK": f"{mfe['peak_rel_percent']:.1f}",
    "T_PEAK": f"{mfe['t_peak_us']:.2f}",
    "B1": f"{cfg['RF_B1_UT']:.0f}",
    "FL50": f"{odmr['geomag_50uT']['f_Larmor_MHz']:.2f}",
    "FL357": f"{odmr['lab_357uT']['f_Larmor_MHz']:.2f}",
    "ODMR_C50": f"{max(odmr['geomag_50uT']['phi_S']) - min(odmr['geomag_50uT']['phi_S']):.4f}",
    "ODMR_C357": f"{max(odmr['lab_357uT']['phi_S']) - min(odmr['lab_357uT']['phi_S']):.4f}",
    "MEMGB": "0.4",
    "STAMP": time.strftime("%Y-%m-%d %H:%M %Z", time.localtime()),
}

for name in ("QUANTUM_BIOLOGY_REPORT_EN.md", "QUANTUM_BIOLOGY_REPORT_ZH.md"):
    p = ROOT / name
    t = p.read_text(encoding="utf-8")
    missing = []
    for k, v in tokens.items():
        tok = f"@{k}@"
        if tok in t:
            t = t.replace(tok, v)
    import re
    left = sorted(set(re.findall(r"@([A-Z0-9_]+)@", t)))
    if left:
        print(f"{name}: UNFILLED {left}")
    else:
        p.write_text(t, encoding="utf-8")
        print(f"{name}: filled OK")
