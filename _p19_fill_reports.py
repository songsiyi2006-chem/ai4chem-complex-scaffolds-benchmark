"""Fill the Phase 19 reports with actual run results from the master JSON."""
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
res = json.loads((ROOT / "results_phase19" / "phase19_results.json").read_text())

gens = res.get("generations", [])
ch = res.get("champion")
unc = res.get("uncatalyzed")
ea = res.get("equivariance_audit", {})


def fmt(x, n=2):
    try:
        return f"{float(x):.{n}f}"
    except (TypeError, ValueError):
        return "n/a"


lines_e = []
lines_z = []

lines_e.append("")
lines_e.append("## 6. Results")
lines_e.append("")
lines_e.append(f"- SE(3) equivariance audit: max translation-field error "
               f"**{ea.get('max_translation_err', float('nan')):.1e}**, "
               f"max body-frame angular error "
               f"**{ea.get('max_angular_err', float('nan')):.1e}** "
               f"(tolerance {ea.get('tolerance', 5e-3)}).")
if unc:
    lines_e.append(f"- Uncatalyzed reference (substrate + H2O, GFN2-xTB/ALPB "
                   f"water): raw barrier **{unc['barrier_kcal']:.2f} kcal/mol** "
                   f"(anchored to the experimental 32.2 kcal/mol for display; "
                   f"the identical offset is applied to the enzyme profile so "
                   f"DDG-act is preserved).")
if ch:
    b = ch["barrier_kcal"]
    ub = (unc or {}).get("barrier_kcal", float("nan"))
    ddg = ub - b
    acc = math.exp(ddg / (8.314462618e-3 * 298.15 / 4.184))
    lines_e.append(f"- **Champion (generation {ch['gen']})**: surrogate "
                   f"catalytic barrier **{b:.2f} kcal/mol** on the achieved "
                   f"geometry, vs the SAME-engine uncatalyzed reference "
                   f"**{ub:.2f} kcal/mol**; TS probe d(C-H) = "
                   f"{fmt(ch.get('ts_d_CH'))} A, d(N-O) = "
                   f"{fmt(ch.get('ts_d_NO'))} A.")
    lines_e.append(f"- Engine-internal barrier reduction "
                   f"**DDG-act = {ddg:.2f} kcal/mol** "
                   f"(~{acc:.0f}x rate enhancement under the identical "
                   f"surrogate). **Caveat:** the fixed-proton relaxed-scan "
                   f"surrogate overestimates absolute barriers "
                   f"(it does not include the concerted N-O cleavage in the "
                   f"scan coordinate); the 32.2 kcal/mol experimental anchor "
                   f"is therefore NOT used for the acceleration here — the "
                   f"comparison is strictly like-for-like inside the engine.")
    lines_e.append(f"- Backbone constellation RMSD (Glu CA + oxyanion N, "
                   f"exact by construction): **"
                   f"{ch['constellation_rmsd_A']:.3f} A** (gate 0.30 A); "
                   f"carboxylate placement deviation (best-effort rotamer, "
                   f"recorded, priced by the QM/MM barrier): "
                   f"**{fmt(ch.get('carboxylate_deviation_A'), 2)} A**; "
                   f"achieved pi-stack deviation from the s4 intent: "
                   f"**{fmt(ch.get('stack_deviation_A'), 2)} A**.")
    lines_e.append(f"- MD (amber14SB/GBn2, {ch['md_ns']:.2f} ns budgeted): "
                   f"catalytic-constellation RMSF **{ch['rmsf_A']:.2f} A** "
                   f"(gate 0.8 A).")
    lines_e.append(f-"" if False else f"- De novo enzyme: "
                   f"**{len(ch['sequence'])} residues**, composition "
                   f"{ch['static']['composition']}, net charge "
                   f"{ch['static'].get('composition', {}).get('+', 0) - ch['static'].get('composition', {}).get('-', 0):+d}, "
                   f"Rg = {ch['static']['radius_gyration']:.1f} A.")
if gens:
    fb = [g["free_energy"] for g in gens]
    mb = [g["predictive_barrier"][0] for g in gens]
    meas = [b for g in gens for b in g["measured_barriers"]]
    lines_e.append(f"- Variational free energy across generations: "
                   f"**{fb[0]:.2f} -> {fb[-1]:.2f}** "
                   f"(complexity + accuracy decomposition in fig1a); "
                   f"belief-predictive barrier "
                   f"**{mb[0]:.1f} -> {mb[-1]:.1f} kcal/mol**.")
    if meas:
        lines_e.append(f"- Measured QM/MM barriers: "
                       + ", ".join(f"{b:.1f}" for b in sorted(meas))
                       + " kcal/mol.")
lines_e.append("")
lines_e.append("**Interpretation.** The engine closes the loop: the flow "
               "proposes folds, the stereochemical projection realizes them "
               "with exact catalytic anchors, OpenMM certifies foldability, "
               "GFN2-xTB/QM/MM prices the achieved catalytic geometry in "
               "kcal/mol, and the Fristonian belief steers the next "
               "generation. The absolute rate claim is bounded by the "
               "embedding-MEP surrogate (Section 7 ledger); the *relative* "
               "descent across generations under the measured barrier is the "
               "engine's honest output.")

text = (ROOT / "ACTIVE_INFERENCE_ENZYME_REPORT_EN.md").read_text()
import re
new = re.sub(r"\(auto-filled from[^)]*\)[^\n]*\n(\| quantity \| value \|\n\|---\|---\|\n)?",
             "", text, count=1)
# insert results table after "## 6. Results" heading
res_block = "\n".join(lines_e)
new = new.replace("## 6. Results\n", "## 6. Results\n" + res_block + "\n", 1)
(ROOT / "ACTIVE_INFERENCE_ENZYME_REPORT_EN.md").write_text(new)

# Chinese report
lines_z.append("")
lines_z.append("## 6. 结果")
lines_z.append("")
lines_z.append(f"- SE(3) 等变性审计：平移场最大误差 "
               f"**{ea.get('max_translation_err', float('nan')):.1e}**，"
               f"体系系角速度最大误差 "
               f"**{ea.get('max_angular_err', float('nan')):.1e}**"
               f"（容差 {ea.get('tolerance', 5e-3)}）。")
if unc:
    lines_z.append(f"- 非催化参考（底物 + H₂O，GFN2-xTB/ALPB 水）：原始势垒 "
                   f"**{unc['barrier_kcal']:.2f} kcal/mol**（展示时锚定于实验值 "
                   f"32.2 kcal/mol；同一偏移用于酶曲线，保持 ΔΔG‡）。")
if ch:
    b = ch["barrier_kcal"]
    ub = (unc or {}).get("barrier_kcal", float("nan"))
    ddg = ub - b
    acc = math.exp(ddg / (8.314462618e-3 * 298.15 / 4.184))
    lines_z.append(f"- **冠军（第 {ch['gen']} 代）**：实现几何上的代理势垒 "
                   f"**{b:.2f} kcal/mol**，同引擎非催化参考 "
                   f"**{ub:.2f} kcal/mol**；过渡态探针 d(C-H) = "
                   f"{fmt(ch.get('ts_d_CH'))} Å，d(N-O) = "
                   f"{fmt(ch.get('ts_d_NO'))} Å。")
    lines_z.append(f"- 引擎内部势垒下降 **ΔΔG‡ = {ddg:.2f} kcal/mol**"
                   f"（同一代理下约 {acc:.0f} 倍速率提升）。**注意：**"
                   f"固定质子弛豫扫描代理高估绝对势垒（扫描坐标未含 N–O 协同"
                   f"断裂）；因此速率加速未使用 32.2 kcal/mol 实验锚点，"
                   f"比较严格限定在引擎内部同口径。")
    lines_z.append(f"- 骨架星座 RMSD（Glu CA + 氧阴离子 N，构造精确）："
                   f"**{ch['constellation_rmsd_A']:.3f} Å**（门 0.30 Å）；"
                   f"羧酸根放置偏差（尽力旋转异构体，已记录，由 QM/MM 势垒"
                   f"诚实计价）：**{fmt(ch.get('carboxylate_deviation_A'), 2)} Å**；"
                   f"π 堆叠相对 s4 意图的实现偏差："
                   f"**{fmt(ch.get('stack_deviation_A'), 2)} Å**。")
    lines_z.append(f"- MD（amber14SB/GBn2，{ch['md_ns']:.2f} ns 预算）：催化星座 "
                   f"RMSF **{ch['rmsf_A']:.2f} Å**（门 0.8 Å）。")
    lines_z.append(f"- 从头酶：**{len(ch['sequence'])} 残基**，组成 "
                   f"{ch['static']['composition']}，Rg = "
                   f"{ch['static']['radius_gyration']:.1f} Å。")
if gens:
    fb = [g["free_energy"] for g in gens]
    mb = [g["predictive_barrier"][0] for g in gens]
    meas = [b for g in gens for b in g["measured_barriers"]]
    lines_z.append(f"- 各代变分自由能：**{fb[0]:.2f} → {fb[-1]:.2f}**"
                   f"（复杂度+精度分解见图 1a）；信念预测势垒 "
                   f"**{mb[0]:.1f} → {mb[-1]:.1f} kcal/mol**。")
    if meas:
        lines_z.append(f"- 实测 QM/MM 势垒："
                       + "、".join(f"{b:.1f}" for b in sorted(meas))
                       + " kcal/mol。")
lines_z.append("")
lines_z.append("**解读。** 引擎闭合了循环：流模型提出折叠，立体化学投影以精确催化"
               "锚点实现之，OpenMM 认证可折叠性，GFN2-xTB/QM/MM 以 kcal/mol 为实现"
               "催化几何计价，Friston 信念引导下一代。绝对速率声明受嵌入-MEP 代理"
               "限制（第 7 节台账）；在实测势垒下跨代的*相对*下降才是引擎诚实的"
               "输出。")

text = (ROOT / "ACTIVE_INFERENCE_ENZYME_REPORT_ZH.md").read_text()
new = re.sub(r"\(由 `results_phase19/phase19_results\.json`[^\n]*\n(\| 量 \| 数值 \|\n\|---\|---\|\n)?",
             "", text, count=1)
res_block = "\n".join(lines_z)
new = new.replace("## 6. 结果\n", "## 6. 结果\n" + res_block + "\n", 1)
(ROOT / "ACTIVE_INFERENCE_ENZYME_REPORT_ZH.md").write_text(new)
print("reports filled")
