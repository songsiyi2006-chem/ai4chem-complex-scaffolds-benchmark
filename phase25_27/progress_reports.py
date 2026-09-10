"""Build bilingual reports from reviewed pilot data and a frozen execution snapshot."""
import json
from .build_inputs import ROOT,OUT


def main():
    data=json.loads((OUT/'phase25/pilot_progress.json').read_text(encoding='utf8'))
    val=json.loads((OUT/'phase25/xtb_pilot/validation.json').read_text(encoding='utf8'))
    ternary=json.loads((OUT/'phase25/ternary_diagnostics.json').read_text(encoding='utf8'))
    snap=json.loads((OUT/'phase25/hessian_checkpoint_snapshot/manifest.json').read_text(encoding='utf8'))
    dft=json.loads((OUT/'phase25/stationary_pilot/pbe_E_opt_01/result.json').read_text(encoding='utf8'))
    thermo=json.loads((OUT/'phase25/thermo_sensitivity/results.json').read_text(encoding='utf8'))
    n=data['xTB_audited_jobs'];nmin=data['xTB_accepted_preliminary_minima'];nt=len(ternary['candidates']);ng=snap['completed_native_gradient_jobs']
    eztable='| Imine | Solvent | G(Z)-G(E), cutoff 50 | G(Z)-G(E), cutoff 100 |\n|---|---|---:|---:|\n'
    for r in data['E_Z_comparisons']:
        g100=r['Z_minus_E_local_minimum_G_cutoff100_kcal_mol'];b='pending' if g100 is None else f'{g100:.3f}'
        eztable+=f"|{r['substrate']}|{r['solvation']}|{r['Z_minus_E_local_minimum_G_cutoff50_kcal_mol']:.3f}|{b}|\n"
    complex_table='| Start | Final E/Z | Proton contact state | Relative G50 (kcal/mol) | Lowest mode (cm^-1) |\n|---|---|---|---:|---:|\n'
    for r in ternary['candidates']:
        complex_table+=f"|{r['name'].replace('C01_I01_','')}|{r['final_imine_E_Z']}|{r['proton_location']}|{r['relative_local_G50_kcal_mol']:.3f}|{r['minimum_frequency_cm1']:.2f}|\n"
    failed=[r for r in val['records'] if r.get('frequencies_below_minus1_cm1')]
    failures=', '.join(f"{r['campaign']}/{r['name']}: {r['lowest_frequency_cm1']:.2f} cm^-1" for r in failed)
    elapsed=sum(r['finished_allocated_core_seconds'] for r in data['xTB_costs'])
    dfttable='| Job | Method | Stage | Status | Energy (Eh) | Final E/Z geometry |\n|---|---|---|---|---:|---|\n'
    for row in data['DFT_jobs']:
        energy='pending' if row['electronic_energy_hartree'] is None else f"{row['electronic_energy_hartree']:.12f}"
        dfttable+=f"|{row['tag']}|{row['method']}|{row['stage']}|{row['status']}|{energy}|{row.get('final_E_Z_geometry','not assessed')}|\n"
    zrow=next((r for r in data['DFT_jobs'] if r['tag']=='pbe_Z_opt_01' and r['status']=='OPTIMIZED_FREQUENCY_PENDING'),None)
    zen=zzh=''
    if zrow:
        z=json.loads((OUT/'phase25/stationary_pilot/pbe_Z_opt_01/result.json').read_text(encoding='utf8'))
        gap=(z['electronic_energy_hartree']-dft['electronic_energy_hartree'])*627.5094740631
        zen=f"Both E and Z optimizations have now converged. Z retains Z geometry; its maximum gradient is {z['max_abs_gradient_hartree_per_bohr']:.3e} Eh/bohr and wall time {z['wall_seconds']/60:.2f} min. The gas-phase PBE-D3BJ electronic E(Z)-E(E) difference is {gap:.3f} kcal/mol. This is not a Gibbs-energy difference or catalytic barrier. It must not be directly attributed to a functional effect against the solvated xTB table, since both method and environment differ."
        zzh=f"E、Z 两套优化均已收敛，Z 几何保持；Z 的最大梯度为 {z['max_abs_gradient_hartree_per_bohr']:.3e} Eh/bohr，耗时 {z['wall_seconds']/60:.2f} 分钟。气相 PBE-D3BJ 的电子能 E(Z)−E(E) 为 {gap:.3f} kcal/mol；它不是自由能差或催化势垒。与下方溶剂化 xTB 表格比较时，方法和环境同时不同，不能把差异直接归因于泛函。"
    en=f'''# Phase 25 — Local calculation progress and evidence review

**Snapshot: {snap['snapshot_utc']}. Scientific acceptance remains NOT PASSED; product selectivity remains undetermined.** This update extends the initial `93b9ab1` release with actual optimizations, Hessians, complete ternary structures and a running local DFT queue. Phase26/27 have no new production calculations in this update.

## Completed and running work

- **{n} xTB pilot jobs** are recorded, including controls and retained failed minimum candidates; **{nmin}** pass the preliminary local-minimum audit. A successful engine exit alone is insufficient.
- Eight complete **117-atom** C01/imine/dimethyl-Hantzsch-ester starts cover E/Z, both approach-face signs and two rigid-placement starts per combination. **{nt}** complete ternary candidates currently pass the xTB-level checks. The catalyst is not truncated. Face signs are construction labels, never assigned R/S product channels.
- The isolated I01-E PBE-D3BJ/def2-SVP optimization converged in **{dft['wall_seconds']/60:.2f} min**, with energy **{dft['electronic_energy_hartree']:.12f} Eh**, maximum Cartesian gradient **{dft['max_abs_gradient_hartree_per_bohr']:.3e} Eh/bohr**, and RMS gradient **{dft['rms_gradient_hartree_per_bohr']:.3e} Eh/bohr**. E geometry is retained. A converged optimization is not a frequency-validated minimum.
- The same-surface native DFT Hessian has **{ng}/163 completed atomic gradient jobs** in the frozen public checkpoint snapshot. **No complete DFT frequencies or DFT Gibbs energy are claimed.** Live output is excluded from the release manifest until separately reviewed.
- The prior four isolated-imine DFT single points remain historical evidence. The background queue adopts the current Hessian process, then runs Z optimization/frequencies and omegaB97X-D single points on the PBE-optimized E/Z geometries. Queued calculations are not counted as completed.

## DFT job ledger at the snapshot

{dfttable}
{zen}

Z optimization was additionally started when free local memory increased; the finite queue adopts or skips that same job to avoid duplicate work. A STARTED/FAILED frequency row has no accepted frequency result.

## Methods and acceptance rules

New tight-binding jobs use xTB 6.7.1, GFN2-xTB, ALPB(toluene) or ALPB(CH2Cl2), charge 0, zero unpaired electrons, accuracy 0.2, and `vtight` optimization. Imaginary-mode remediation uses `extreme`. Optimization and Hessian use the same method and solvent. Runs are isolated, commands and raw output retained, and the Hessian follows a confirmed optimization using `--hess --strict`. The [xTB manual](https://xtb-docs.readthedocs.io/en/latest/hessian.html) explains the Hessian and imaginary-mode handling.

The independent audit requires the complete 3N frequency block, agreement between repeated printouts, six projected translation/rotation modes, and all internal modes above 1 cm^-1. This last threshold is a numerical screen, not a demonstrated frequency error bound. The audit also checks atom ordering, isolated-imine E/Z identity, and, for complete ternaries, catalyst-axis sign and expected heavy-bond distances. Soft modes and all rejected candidates remain visible. The screen is not a substitute for a chemically reviewed DFT stationary point or IRC.

The xTB thermochemistry summary can report zero imaginary frequencies after its default -20 cm^-1 cutoff even when a printed mode is negative. The initial I02-E structures had negative modes. Both signs of the engine-generated displacement were reoptimized; successful replacement candidates were retained while the others stayed rejected. Negative-mode records: **{failures}**. No absolute-value conversion hides these modes. Any negative-mode ternary candidate listed here remains unresolved and is excluded from the accepted-candidate table; the imine remediation does not certify that ternary.

## Solvent and low-frequency diagnostics

The table reports **single-local-minimum E/Z differences at the xTB level**, in kcal/mol. It is not a conformer-ensemble result, activation barrier, population prediction or reaction ee. The best audited candidate by G50 is selected per E/Z species from this limited search. Standard-state constants cancel in each equal-composition 1-to-1 comparison; absolute 1 M solution free energies and binding free energies are not certified.

{eztable}
![E/Z cutoff sensitivity](results_phase25_27/phase25/figures/imine_EZ_sensitivity.png)

The 50/100 cm^-1 comparison reuses the actual saved Hessian in **native `xtb thermo`**, retaining its geometry-dependent inertia. All **{sum(r['replay50_pass'] for r in thermo)}/{len(thermo)}** reference replays agree with the original G correction within 1e-6 Eh. These are method-sensitivity contrasts, not confidence intervals. No fixed-inertia approximation is used in this reported table. The defaults are documented in [xcontrol](https://github.com/grimme-lab/xtb/blob/v6.7.1/man/xcontrol.7.adoc); the [thermodynamic implementation](https://github.com/grimme-lab/xtb/blob/v6.7.1/src/thermo.f90) uses the molecular rotational inertia. ALPB's default `gsolv` reference is retained and identified using the [reference-state documentation](https://xtb-docs.readthedocs.io/en/latest/gbsa.html#reference-states).

## Complete-system structures and controls

{complex_table}
Relative G50 refers to the lowest audited candidate in this finite set. **These are minima, not competing transition states; converting this table to ee would be invalid.** Fixed-mapping heavy-atom RMSDs are provided, but search counts are not statistical degeneracies and no ensemble weights are assigned. Absolute catalyst axial CIP and product CIP remain pending.

![Complete optimized ternary and contacts](results_phase25_27/phase25/figures/ternary_minimum.png)

For E_face-1_start0, O...H is 1.629 A and N-H is 1.048 A; the designated donor C-H remains 1.094 A. At this approximate level, the geometry supports an iminium/phosphate ion-pair candidate without completed hydride transfer. Other starts include O-bound neutral hydrogen-bond candidates. These observations motivate separate protonation-state paths; they do not establish the experimental mechanism. The exact reflected full complex has a same-geometry xTB energy difference of **{ternary['mirror_control']['delta_energy_kcal_mol']:.3e} kcal/mol at printed precision**. This parity check is not a kinetic ee-inversion test. Achiral and uncatalyzed rate controls remain uncomputed.

## Resources, failure recovery and unfinished criteria

This work uses the existing local CPU environment. Psi4 requests 500 MB and 2 threads for the completed optimization, then 3 threads for the Hessian/queue; xTB uses one thread. These are requested resources, not peak-memory guarantees. Recorded finished xTB jobs consumed **{elapsed:.1f} allocated core-seconds** (wall time times assigned threads), including unsuccessful minimum candidates; active jobs and unrelated campaigns are excluded. This is not a cost-matched active-learning comparison.

The first DFT frequency attempt failed while decoding UTF-8 native output with Windows GBK, after its first gradient had completed. The failed output is retained. `python -X utf8` repairs that interface failure without modifying the installed engine. Exact QCSchema inputs and successful AtomicResults are hash-bound and checkpointed. A real small native-gradient test verifies fresh execution, identical cache replay and rejection of corrupted bindings. The published snapshot contains only completed, checked checkpoints. The local finite queue continues while the computer/processes remain running; it does not purchase compute, use OpenAI inference GPUs, publish automatically, or guarantee eventual convergence.

There are still no accepted DFT R/S catalytic transition states or IRCs, no correlated-wavefunction validation of selectivity-determining structures, no explicit-solvent free-energy replicas, no verified Curtin-Hammett separation, and no external blind outcome labels. Therefore no rates, conversion, signed ee, confidence coverage or 0.5 kcal/mol ensemble-convergence acceptance is reported. Phase26 constant-potential sampling and Phase27 multireference dynamics remain pending.

## Reproduction and evidence

[Progress data](results_phase25_27/phase25/pilot_progress.json), [raw xTB jobs](results_phase25_27/phase25/xtb_pilot), [independent audit](results_phase25_27/phase25/xtb_pilot/validation.json), [ternary diagnostics](results_phase25_27/phase25/ternary_diagnostics.json), [native thermochemistry](results_phase25_27/phase25/thermo_sensitivity/results.json), [DFT optimization](results_phase25_27/phase25/stationary_pilot/pbe_E_opt_01/result.json), [frozen native gradients](results_phase25_27/phase25/hessian_checkpoint_snapshot/manifest.json), [execution snapshot](results_phase25_27/phase25/local_execution_snapshot.json), [publication omissions](results_phase25_27/publication_omissions.json), [commands and limits](phase25_27/README.md).
'''
    zh=f'''# Phase25 — 本地计算推进与证据复核

**快照时间：{snap['snapshot_utc']}。科学任务仍未通过整体验收，选择性未判定。** 本报告承接初版提交 `93b9ab1`，新增真实优化、频率、完整三组分结构及本地 DFT 在途计算。Phase26/27 本轮没有新增生产级计算。

## 已完成与仍在运行的工作

- 已记录 **{n} 个 xTB 作业**，包括对照和未通过的候选结构；其中 **{nmin} 个**通过该近似层级的局部极小点筛查。程序正常退出不等于结构验收通过。
- 构建 8 个 **117 原子**完整 C01 磷酸—亚胺—二甲酯型 Hantzsch 酯起点，覆盖 E/Z、两侧进攻方向及每组合两个起点；当前 **{nt} 个**完整复合物通过 xTB 检查。没有截断催化剂。face 的正负只是起点方向，不代表 R/S 产物。
- I01-E 的 **PBE-D3BJ/def2-SVP 气相优化已收敛**，耗时 {dft['wall_seconds']/60:.2f} 分钟；能量为 **{dft['electronic_energy_hartree']:.12f} Eh**，最大及均方根笛卡尔梯度分别为 **{dft['max_abs_gradient_hartree_per_bohr']:.3e}**、**{dft['rms_gradient_hartree_per_bohr']:.3e} Eh/bohr**，E 构型保持。优化收敛尚不等于通过频率验证。
- 同势能面的完整 DFT 数值 Hessian 已在公开冻结快照中保留 **{ng}/163 个已完成梯度作业**。**尚不报告完整 DFT 频率或 DFT 自由能。** 在途文件不进入发布哈希清单，完成后需单独审核。
- 历史 4 条孤立亚胺 DFT 单点保留。当前本地队列接续 E 的 Hessian，再依次运行 Z 优化、Z 频率，以及在 PBE 优化的 E/Z 几何上的 ωB97X-D 单点。排队项目不计作完成。

## 快照中的 DFT 作业状态

{dfttable}
{zzh}

空闲内存回升后，已额外启动 Z 优化；队列会接管或跳过同一作业，避免重复计算。频率行中的 STARTED/FAILED 均不代表已得到通过验收的频率结果。

## 方法与验收规则

使用 xTB 6.7.1、GFN2-xTB、ALPB（甲苯/二氯甲烷）、总电荷 0、未配对电子数 0、accuracy=0.2；常规优化使用 vtight，虚频修正使用 extreme。优化与 Hessian 的方法和溶剂保持一致；先确认优化收敛，再执行 `--hess --strict`。原始命令、日志、几何和 Hessian 均保留。[官方频率说明](https://xtb-docs.readthedocs.io/en/latest/hessian.html)

独立解析器要求 3N 个完整模式、重复打印块一致、6 个平移/转动模式已投影，以及全部内部模式高于 1 cm⁻¹；另检查原子顺序、孤立亚胺 E/Z、完整复合物的轴手性几何符号与重原子连接距离。1 cm⁻¹ 是数值筛查阈值，并非已证明的频率误差界限。该筛查也不能替代经过化学审阅的 DFT 驻点和 IRC。

初始 I02-E 在两种溶剂中均出现虚频。xTB 的热化学汇总默认使用 −20 cm⁻¹ 截断，可能在原始模式为负时仍显示“0 个虚频”。我们按原始模式拒绝这些结构，沿引擎给出的虚频方向正、负两侧重优化，保留成功替代结构和失败记录。全部负频候选：**{failures}**。没有对负频取绝对值后冒充极小点。列表中带负频的三组分候选仍未解决，已从通过筛查的候选表中排除；亚胺的修正结果不构成该三组分结构的验收。

## 溶剂与低频处理的实算结果

下表单位为 kcal/mol，是**有限搜索中单个局部极小点的 xTB 层级 G(Z)−G(E)**。按 G50 选择各 E/Z 已通过筛查的候选。它不是构象系综、活化势垒、真实 E/Z 布居或反应 ee；同组成的一对一比较中共同标准态常数抵消，尚不认证绝对 1 M 溶液自由能或结合自由能。

{eztable}
![E/Z 与低频敏感性](results_phase25_27/phase25/figures/imine_EZ_sensitivity.png)

50/100 cm⁻¹ 对照由 **xTB 原生 thermo** 重用已保存的真实 Hessian 完成，保留其依赖分子几何的转动惯量。原 50 cm⁻¹ 热化学校正重放 **{sum(r['replay50_pass'] for r in thermo)}/{len(thermo)}** 在 1e-6 Eh 内一致。图中连线是处理方法的敏感性，**不是置信区间**；表格没有使用固定转动惯量近似。依据：[参数说明](https://github.com/grimme-lab/xtb/blob/v6.7.1/man/xcontrol.7.adoc)、[热化学实现](https://github.com/grimme-lab/xtb/blob/v6.7.1/src/thermo.f90)、[ALPB 标准态定义](https://xtb-docs.readthedocs.io/en/latest/gbsa.html#reference-states)。本轮保留并注明默认 gsolv 参照。

## 完整体系与镜像对照

{complex_table}
相对 G50 以本轮有限集合中最低已通过候选为零点。**这些是复合物极小点，不是竞争过渡态，不能代入 R/S 选择性公式。** 固定原子对应的重原子 RMSD 已保存；搜索次数不当作物理简并度，没有人为赋予系综权重。催化剂绝对轴手性和产物 CIP 认证仍未完成。

![完整三组分结构与接触距离](results_phase25_27/phase25/figures/ternary_minimum.png)

以 E_face-1_start0 为例：O···H=1.629 Å、N−H=1.048 Å，指定供氢位点的 C−H 仍为 1.094 Å。在此近似层级，它是未完成氢转移的亚胺鎓—磷酸根离子对候选；其他起点也得到 O 端结合的中性氢键候选。因此后续 DFT 搜索必须区分质子化状态，但这些结果尚不能确证实验机制。对完整结构作精确镜像后，同几何 xTB 单点能差在打印精度下为 **{ternary['mirror_control']['delta_energy_kcal_mol']:.3e} kcal/mol**。这只检验能量的镜像对称性，不等于验证产物 ee 符号翻转。非手性和无催化剂的速率对照尚未完成。

## 算力、失败恢复和未完成项

使用现有 CPU 环境；已完成的 DFT 优化申请 500 MB、2 线程，Hessian/队列申请 500 MB、3 线程，xTB 使用 1 线程。申请值不是峰值内存保证。已结束 xTB 作业累计记录 **{elapsed:.1f} 核秒**（墙钟时间×分配线程数），含未通过的候选；不含仍在运行的作业，也不冒充等成本主动学习比较。

首轮 DFT 频率在一个梯度完成后，因 Windows GBK 解码 UTF-8 原始输出失败；失败日志已保留。使用 `python -X utf8` 修复接口问题，没有修改安装的量化引擎。现在逐点绑定 QCSchema 输入与成功结果的哈希，保存可复用检查点；另用真实微型梯度作业验证首次计算、精确重放和损坏检查点拒绝。公开快照只包含已完成且已核查的数据。本地有限队列需要电脑及进程保持运行；不会自动购买算力、调用 OpenAI 推理 GPU 或自动发布，也不保证后续收敛。

仍缺完整 DFT R/S 过渡态及 IRC、关键结构的相关波函数复核、显式溶剂自由能独立重复、Curtin–Hammett 时间尺度验证和真正外部盲测标签。因此不报告速率、转化率、有符号 ee、区间覆盖率或 0.5 kcal/mol 系综收敛通过。Phase26 的恒电位采样及 Phase27 的多参考非绝热动力学仍待执行。

## 数据入口

[计算进展](results_phase25_27/phase25/pilot_progress.json) · [原始 xTB 作业](results_phase25_27/phase25/xtb_pilot) · [独立频率核查](results_phase25_27/phase25/xtb_pilot/validation.json) · [三组分诊断](results_phase25_27/phase25/ternary_diagnostics.json) · [原生热化学对照](results_phase25_27/phase25/thermo_sensitivity/results.json) · [DFT 优化](results_phase25_27/phase25/stationary_pilot/pbe_E_opt_01/result.json) · [梯度冻结快照](results_phase25_27/phase25/hessian_checkpoint_snapshot/manifest.json) · [本地运行快照](results_phase25_27/phase25/local_execution_snapshot.json) · [发布排除项](results_phase25_27/publication_omissions.json) · [复现说明](phase25_27/README.md)。
'''
    for lang,text in [('EN',en),('ZH',zh)]:
        (ROOT/f'PHASE25_PROGRESS_{lang}.md').write_bytes(text.encode('utf8'))
    print('Wrote bilingual calculation-progress reports from the frozen snapshot')


if __name__=='__main__':main()
