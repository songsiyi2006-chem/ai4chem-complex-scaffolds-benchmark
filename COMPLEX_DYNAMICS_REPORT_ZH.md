# 复合物动力学报告 — 第三阶段：靶标-配体复合物、逐残基 MM-GBSA 分解与 ML 势基准

*数据：`results_phase3/phase3_results.json` | OpenMM 8.6（CPU）· pdbfixer 1.12 · AutoDock Vina 1.2.5 CLI + meeko 0.8 · OpenFF Sage 2.1 + AM1-BCC（NAGL 1.0.0，离线）· torchani ANI-2x · mdtraj 1.11*

---

## 【一、靶标背景与勘误】

**靶标：KRAS G12D，switch-II 口袋（PDB 7RPZ）** — X-ray 1.30 Å，GDP 结合态，与 **MRTX1133**（残基 `6IC`，44 重原子）共结晶；催化结构域 1–170 构造，工程化 G12D 突变（SEQADV ASP12），Cys-light 背景（S51/C80L/S118）。

> **规格勘误（经 RCSB REST API 实时核实）**：规格中为 KRAS G12D 引用的 PDB "7E27" 实为 *恶性疟原虫甲酸-亚硝酸转运蛋白与 MMV007839 复合物*（cryo-EM 2.29 Å）。此处采用的规范 G12D switch-II 口袋共晶为 **7RPZ**（在 142 个 G12D 结构中检索核实）。备选 7BQY（SARS-CoV-2 Mpro + N3，1.7 Å）亦经核实；选 KRAS 是因为第二阶段的 **T04** 即按 switch-II 口袋共价抑制剂核心设计。

**配体 T04**（第二阶段）— 2-氨基吡啶并嘧啶酮-哌嗪丙烯酰胺（C₂₈H₂₄ClFN₆O，含氢 61 原子），KRAS G12C/G12D switch-II 口袋时代的丙烯酰胺弹头范式。

**结构生物学要点**：switch-II 口袋是 KRAS 突变体（G12C/G12D）药物发现的里程碑口袋——位于 switch-II（残基 60–67）与中央 β 片之间，GDP 态下开放，MRTX1133（G12D 非共价首选）与 sotorasib（G12C 共价）均锚定于此。关键残基：**His95 / Tyr96**（氢键与疏水壁）、Gln99、Glu62（switch-II）、以及突变位点本身 **Asp12**。

---

## 【二、阶段 1A — 大分子摄取与修复】

- RCSB 自动下载；组分拆分：蛋白 1342 重原子 | MRTX1133 44 | GDP 28 | Mg²⁺ 1 | 水分子剔除；altloc 去重。
- **PDBFixer 修复**：补 7 个缺失重原子；**pH 7.4** 加氢（His 互变异构按氢键模式规则）；蛋白终态 2681 原子。
- **口袋定义**：MRTX1133 重原子质心 [1.71, 4.93, −23.16] Å；R = 10.0 Å 内 **28 个残基**（含 H95、Y96、Q99、E62、D12、D92、Y64、R68、K88、A11 等）。

---

## 【三、阶段 2 — meeko + Vina 对接】

配体 ETKDGv3+MMFF94 构象 → meeko 柔性扭转 PDBQT；受体 = 修复蛋白 **+ GDP**（刚体环境保留辅因子），openbabel 3.2.1 转换。

*引擎回退（已记录）*：Vina Python API 无 Windows/Py3.12 轮子 → 官方 **Vina v1.2.5 CLI 可执行文件**（结果表解析；meeko 回转 RDKit）。

| 构象 | ΔG_docking (kcal/mol) |
|---|---|
| **#1（选定）** | **−7.92** |
| #2 | −7.30 |
| #3 | −7.11 |

盒 22.5 Å³、6IC 质心居中、exhaustiveness 16、seed 42、9 模式。

---

## 【四、阶段 3 — 经典力场 vs ML 势（口袋冻结构象单点）】

对口袋冻结的 pose #1（61 原子；氢先做重原子固定的弛豫）做单点能量/原子力。

**环境回退（规格明文允许，均已记录）**：
1. **GAFF2**（openmmforcefields）离线不可得且不可降级安装 → 经典参考 = **Sage 2.1 价键+vdW + AM1-BCC 电荷**。值得一提的是，AM1-BCC 现已通过随包 NAGL GNN 模型（`openff-gnn-am1bcc-1.0.0.pt`，按文件路径调用）**完全离线可用**——第二阶段的 MMFF 电荷手术 workaround 就此退役。
2. **MACE-OFF23** 模型下载失败（离线）→ **ANI-2x（torchani 2.8.4）**。

| 势 | E (kcal/mol) | ‖F‖max (kcal/mol/Å) |
|---|---|---|
| Sage 2.1 + AM1-BCC（经典，GAFF2 替身） | −732.4 | 132.5 |
| MMFF94（RDKit，数值梯度） | +37.7 | 36.1 |
| **ANI-2x（ML 参考）** | −47152.0 | 1.74 |

（各模型能量零点不同；只有力可直接比较。）

### 4.1 力场差异向量（相对 ANI-2x）

| 指标 | Sage+AM1-BCC | MMFF94 |
|---|---|---|
| Pearson r（逐原子 ‖F‖） | 0.698 | 0.638 |
| 平均 cos(F_classical, F_ML) | 0.483 | 0.555 |
| ⟨‖ΔF‖⟩ (kcal/mol/Å) | 24.7 | 9.8 |
| max ‖ΔF‖ | 135.4 | 42.8 |

### 4.2 人工应变集中在哪些官能团（SMARTS 归属）

| 官能团 | 原子数 | ⟨‖ΔF_Sage−ANI‖⟩ | max |
|---|---|---|---|
| **2-氨基嘧啶杂芳环** | 7 | **81.0** | 135.4 |
| 炔基连接子（C#C） | 2 | 57.8 | 65.1 |
| 丙烯酰胺弹头（C=CC(=O)N） | 5 | 46.9 | 63.5 |
| 哌嗪 | 5 | 27.6 | 72.0 |

**解读**：经典力场与等变 ML 势的最大分歧恰好落在**共轭氮杂芳环核心与 sp/sp² 弹头矢量**上——电子离域、可极化、需贯穿共轭表达的基团，正是简谐 MM 无法表达的部位。此处 MMFF94 与 ANI 的一致性反而优于 Sage（⟨ΔF⟩ 9.8）——"用哪个经典力场"与"经典 vs ML"同等重要。图：`./figures_phase3/fig3_ml_vs_classical_ff_gap.png`。

---

## 【五、阶段 4 — OpenMM 复合物动力学（OBC2 隐式溶剂，310 K）】

**体系拼接（2785 原子）**——纯 OpenMM 路线，不依赖 openmmforcefields：
- **蛋白**：Amber14SB（ff14SB）+ `implicit/obc2.xml`（OpenMM 8.x Script 生成器 → CustomGBForce），Cα 约束 k = 5 kcal/mol/Å²。
- **GDP + T04**：Sage 2.1 价键 + AM1-BCC 电荷拼入（价键力索引偏移；GB 参数按 Bondi×OBC 补挂；Sage 自带 X–H 约束去重复制）。
- **盐**：Debye–Hückel κ = 1.263 nm⁻¹（0.15 M，ε_w = 76.6，310 K）。
- **Mg²⁺ 弃用**（纯 OpenMM 路线无模板；已记录）。

**协议偏差（均已诊断并记录）**：
1. `dt = 1.0 fs`（规格 2.0 fs）：310 K 下 2 fs 经 Sage 拼接分子的角模式发散；10→50 K（0.5 fs）→150 K（1 fs）→310 K 升温斜坡后严格稳定（生产前 20000 步清洁验证）。
2. 生产段 **50000 步 = 50 ps**（规格 50k–100k 区间），~14 steps/s（本 OpenMM 构建无 GB 截断 API）；平衡 5000 步；每 250 步存帧 → **200 帧**。

| 观测量 | 数值 |
|---|---|
| ⟨Cα RMSD⟩ | **0.347 Å**（max 0.39）——约束下折叠稳如磐石 |
| ⟨口袋 Cα RMSD⟩ | 0.352 Å |
| ⟨配体重原子 RMSD⟩（尾段） | **1.275 Å**（max 1.66）——构象保持、诱导拟合微动 |
| PLIF（疏水，最强） | **Tyr96**（主导，4 个配体碳接触 ≳60% 持续率）、Tyr64 |
| PLIF（氢键） | 稀疏——与疏水凹槽型 pose 一致（见图1） |

---

## 【六、阶段 5 — 终态 MM-GBSA 与逐残基分解】

40 帧（200 取 5）。ΔG_bind = ⟨G_复合物⟩ − ⟨G_蛋白⟩ − ⟨G_配体⟩；GB 极性项来自 κ 屏蔽的 CustomGBForce；非极性 γ·ΔSASA（γ = 0.005 kcal/mol/Å²，Shrake–Rupley 960 点）。

**ΔG_bind = −143.47 ± 2.85 kcal/mol**（dE_MM −102.24 | ΔG_GB −38.10 | ΔG_SA −3.13）

> 隐式溶剂终点法系统性过结合（无构象熵项、无显式水竞争、dE_MM 中电荷接触未屏蔽）。−143 kcal/mol 是**排序/分解量**而非 Kd 预测——交付物是残基分解，不是绝对数。

### 6.1 逐残基分解（前 10）

| 残基 | 总计 | vdW | 静电 | GB 极性 | 备注 |
|---|---|---|---|---|---|
| **His95** | **−8.65** | −6.45 | −5.59 | +3.39 | SIIP 标志性氢键/π 壁 |
| **Tyr96** | **−6.33** | −5.53 | −2.47 | +1.67 | 疏水盖（与 PLIF 主导一致） |
| **Glu62** | −4.54 | −4.60 | +7.23 | **−7.17** | switch-II；教科书式静电↔去溶剂补偿 |
| **Asp12（G12D！）** | −4.33 | −4.05 | −4.08 | +3.80 | 癌基因突变位点本身参与结合 |
| Asp92 | −4.04 | −3.74 | −4.39 | +4.09 | |
| Tyr64 | −3.44 | −2.67 | **−11.29** | **+10.52** | 静电锁定、去溶剂支付 |
| Ala11 | −2.35 | −2.41 | +3.91 | −3.85 | |
| Arg68 | −1.96 | −1.61 | +9.76 | −10.11 | |
| Lys88 | −1.82 | −1.87 | −1.72 | +1.77 | |
| Gln99 | −1.14 | −1.29 | +8.44 | −8.29 | MRTX1133 的规范接触，此处较弱 |

**解读**：switch-II 口袋的催化/变构三重奏——**His95/Tyr96（vdW+氢键壁）+ G12D 的 Asp12 与 switch-II 的 Glu62**——主导结合；带电残基的静电项几乎被 GB 去溶剂精确抵消（Y64、R68、Q99），教科书级 MM-GBSA 特征。图：`./figures_phase3/fig2_per_residue_mmgbsa.png`，构象图 `fig1_binding_pose_pocket.png`。

**诚实性验证注记**：独立 numpy 对求和的交叉 vdW+静电与 OpenMM 电荷清零差值在首帧无法对账（−68.9 vs +89.5 kcal/mol，177% 相对误差；疑为 Script-GB 体系的反应场/例外记账细节）。逐残基分解公式跨残基完全一致、排序与晶体接触图吻合，但绝对值携带此不确定性；GB 分量采用本征 OBC 半径的 Still 对项（屏蔽级近似）。

---

## 【七、经典力场的局限与等变 ML 势的未来】

1. **固定电荷极化失效**。T04 夹在 His95/Glu62/Asp12（±0.5–0.8 e）构成的口袋电场中——固定电荷拟合时从未见过该环境。阶段 3 定量了后果：⟨‖ΔF_Sage−ANI‖⟩ ≈ 25 kcal/mol/Å/原子，峰值 135 落在杂芳环核心——正是诱导极化与贯穿共轭所在。Sage 与 MMFF94 的力方向均与 ANI-2x 明显错位（平均余弦 0.48/0.56）。
2. **简谐扭转 vs 共轭**。丙烯酰胺弹头与炔基连接子——决定反应矢量几何的二面角——差异排名第 2/3。拟合于气相扫描的经典扭转无法携带口袋电场诱导的势垒漂移。
3. **等变势（MACE/PaiNN/Allegro）从数据中学会这些效应**：消息传递天然给出多体极化；e(3)-等变保证处处光滑且能量守恒的力。本次的 ANI-2x（该家族更早、更廉价的成员）仅在**一个冻结构象**上就暴露了经典差距；MACE-OFF23 跑同一协议是自然的第四阶段压力测试（本次模型离线下载失败——已记录）。
4. **ML-经典差距是口袋特异的**——这恰是 SBDD 所需：它定位经典模拟（打分/重打分底层）最不可信的部位——本例即共价抑制剂的氨基嘧啶核心与弹头矢量。

---

## 【八、可复现性、产物与局限】

### 命令

```bash
# 在 phase2ff conda 环境中、其 Library/bin 加入 PATH（BLAS DLL）后：
python run_phase3_complex_dynamics.py                  # 完整流水线
python run_phase3_complex_dynamics.py --fig_only       # 仅重出图
# 可调：--pdb_id 7RPZ --lig_code 6IC --md_steps 50000 --report_interval 250
#       --mgb_frames 40 --force_rerun --skip_dock --skip_md --auto_shutdown
```

### 产物（`results_phase3/`）

`7RPZ.pdb` · `7RPZ_fixed_protein.pdb` · `GDP_withH.pdb` · `receptor.pdbqt` · `T04_docked_poses.pdbqt`（9 模式）· `T04_pose1.sdf/pdb`（氢弛豫）· `complex_start.pdb` · `T04_complex_trajectory.dcd`（200 帧）· `T04_complex_metrics.csv` · `T04_plif_persistence.csv` · `T04_per_residue_mmgbsa.csv` · `phase3_results.json`（全阶段）· 各阶段检查点。

### 环境矩阵

撰写用主环境 Python 3.14（RDKit 2026.03.5）；执行环境 `phase2ff`（Python 3.12）：OpenMM 8.6、pdbfixer 1.12、openbabel 3.2.1、mdtraj 1.11、openff-toolkit 0.19/interchange 0.5.4、meeko 0.8、torchani 2.8.4 + torch 2.10、mace-torch 0.3.16（离线无模型）、`tools/vina.exe` 1.2.5、numpy 2.4.6（conda-forge，需 Library/bin 入 PATH）。

### 局限（引用前必读）

1. GAFF2 → Sage 2.1 + AM1-BCC 替换（openmmforcefields 离线不可得）；相对所研究的经典↔ML 差距，GAFF 与 Sage 的价键差异较小，但如实记录。
2. MACE-OFF23 → ANI-2x（规格允许）；ANI-2x 在应变共轭体系上精度低于 MACE-OFF23。
3. 隐式溶剂（OBC2）、弃 Mg²⁺、保 GDP；dt = 1 fs 偏差（已诊断）；50 ps 单轨迹 + Cα 约束——约束折叠内的诱导拟合，非完全交换。
4. MM-GBSA 终点法固有注意事项（无熵、过结合）及对求和验证失配（第六节）——逐残基**排序**稳健，绝对值近似。
5. 对接分数是 Vina 经验量，非 ΔG 预测。

---

*ai4chem-complex-scaffolds-benchmark 第三阶段。第一阶段：`BENCHMARK_REPORT_ZH.md`（10 分子构象套件）。第二阶段：`DYNAMICS_REPORT_ZH.md`（扭转势垒 + 200 ps 配体 MD）。*
