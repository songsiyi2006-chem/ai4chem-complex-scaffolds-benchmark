# 阶段导航 / Phase index

[返回首页](../README.md) · [证据说明](EVIDENCE.md) · [复现指南](REPRODUCE.md)

每个阶段已建立独立项目文件夹；点击阶段编号直接进入。按研究主题查找全部 30 个阶段。Phase 1–27 的“边界”来自历史报告、[Phase 1–19 收尾清单](../shared/audit/reports/PHASE1_19_CLOSEOUT.md)和[25–27 验收台账](../shared/phase25_27/results/acceptance.json)，本次未重跑。Phase 28–29 新增公开数据分析和本机计算；Phase30新增规格与未校准降阶原型，具体范围见各项目报告。
“数据”可能包含历史输出；较新的审计与重算说明优先于早期图表。Phase 编号是模块标识，不表示成熟度。

<a id="molecules"></a>

## 分子与反应 / Phases 1–8

| Phase | 研究问题 | 报告 | 实现 | 数据 / 图表 | 证据边界 |
|---|---|---|---|---|---|
| [01](../projects/phase01/) | 复杂骨架与构象 | [中文](../projects/phase01/reports/BENCHMARK_REPORT_ZH.md) / [EN](../projects/phase01/reports/BENCHMARK_REPORT_EN.md) | [代码](../projects/phase01/code/molecule_benchmark.py) | [数据](../projects/phase01/results) · [图表](../projects/phase01/figures) | 所列检查通过；保留输入与适用性局限 |
| [02](../projects/phase02/) | 扭转扫描与分子动力学 | [中文](../projects/phase02/reports/DYNAMICS_REPORT_ZH.md) / [EN](../projects/phase02/reports/DYNAMICS_REPORT_EN.md) | [代码](../projects/phase02/code/run_heavy_dynamics_benchmark.py) | [数据](../projects/phase02/results) · [图表](../projects/phase02/figures) | 所列检查通过；注意替代力场与电荷 |
| [03](../projects/phase03/) | KRAS 复合物与结合分析 | [中文](../projects/phase03/reports/COMPLEX_DYNAMICS_REPORT_ZH.md) / [EN](../projects/phase03/reports/COMPLEX_DYNAMICS_REPORT_EN.md) | [代码](../projects/phase03/code/run_phase3_complex_dynamics.py) | [数据](../projects/phase03/results) · [图表](../projects/phase03/figures) | 完整 MD / MM-GBSA 待验收 |
| [04](../projects/phase04/) | 骨架编辑与反应路径 | [中文](../projects/phase04/reports/SKELETAL_EDITING_REPORT_ZH.md) / [EN](../projects/phase04/reports/SKELETAL_EDITING_REPORT_EN.md) | [代码](../projects/phase04/code/run_phase4_reaction_mechanism.py) | [数据](../projects/phase04/results) · [图表](../projects/phase04/figures) | 独立精修驻点/单虚频通过；IRC 未验证 |
| [05](../projects/phase05/) | 反应网络与逆向设计 | [中文](../projects/phase05/reports/WORLD_MODEL_REPORT_ZH.md) / [EN](../projects/phase05/reports/WORLD_MODEL_REPORT_EN.md) | [代码](../projects/phase05/code/run_phase5_chemical_world_model.py) | [数据](../projects/phase05/results) · [图表](../projects/phase05/figures) | 多虚频 TS 不合格；近零产率保留 |
| [06](../projects/phase06/) | 显式溶剂与增强采样 | [中文](../projects/phase06/reports/METADYNAMICS_REPORT_ZH.md) / [EN](../projects/phase06/reports/METADYNAMICS_REPORT_EN.md) | [代码](../projects/phase06/code/run_phase6_explicit_metadynamics.py) | [数据](../projects/phase06/results) · [图表](../projects/phase06/figures) | 轨迹完成但产物采样不足 |
| [07](../projects/phase07/) | 强关联与模型适用边界 | [中文](../projects/phase07/reports/FRONTIER_EPISTEMIC_REPORT_ZH.md) / [EN](../projects/phase07/reports/FRONTIER_EPISTEMIC_REPORT_EN.md) | [代码](../projects/phase07/code/run_phase7_strong_correlation_wall.py) | [数据](../projects/phase07/results) · [图表](../projects/phase07/figures) | 全点/统一基组复核待完成 |
| [08](../projects/phase08/) | 光化学与非绝热动力学 | [中文](../projects/phase08/reports/PHOTOCHEMISTRY_REPORT_ZH.md) / [EN](../projects/phase08/reports/PHOTOCHEMISTRY_REPORT_EN.md) | [代码](../projects/phase08/code/run_phase8_photochemical_dynamics.py) | [数据](../projects/phase08/results) · [图表](../projects/phase08/figures) | 完整 QC 与动力学待验收 |

<a id="automation"></a>

## 自动化与定律发现 / Phases 9–12

| Phase | 研究问题 | 报告 | 实现 | 数据 / 图表 | 证据边界 |
|---|---|---|---|---|---|
| [09](../projects/phase09/) | 自驱动实验室模拟 | [中文](../projects/phase09/reports/SELF_DRIVING_LAB_REPORT_ZH.md) / [EN](../projects/phase09/reports/SELF_DRIVING_LAB_REPORT_EN.md) | [代码](../projects/phase09/code/run_phase9_self_driving_lab_compiler.py) | [数据](../projects/phase09/results) · [图表](../projects/phase09/figures) | 模拟与协议检查通过；无硬件验证 |
| [10](../projects/phase10/) | 连续流反应器数字孪生 | [中文](../projects/phase10/reports/FLOW_CYBERPHYSICAL_REPORT_ZH.md) / [EN](../projects/phase10/reports/FLOW_CYBERPHYSICAL_REPORT_EN.md) | [代码](../projects/phase10/code/run_phase10_cyberphysical_flow_twin.py) | [数据](../projects/phase10/results) · [图表](../projects/phase10/figures) | 所列检查通过；控制模型局限 |
| [11](../projects/phase11/) | 神经波函数与 VMC | [中文](../projects/phase11/reports/NEURAL_WAVEFUNCTION_REPORT_ZH.md) / [EN](../projects/phase11/reports/NEURAL_WAVEFUNCTION_REPORT_EN.md) | [代码](../projects/phase11/code/run_phase11_neural_wavefunction_vmc.py) | [数据](../projects/phase11/results) · [图表](../projects/phase11/figures) | 完整计算结束；He 精度目标未达 |
| [12](../projects/phase12/) | 动力学定律发现 | [中文](../projects/phase12/reports/SCIENTIFIC_AGI_MANIFESTO_ZH.md) / [EN](../projects/phase12/reports/SCIENTIFIC_AGI_MANIFESTO_EN.md) | [代码](../projects/phase12/code/run_phase12_hamiltonian_law_discovery.py) | [数据](../projects/phase12/results) · [图表](../projects/phase12/figures) | 完整训练与验收待完成 |

<a id="biological"></a>

## 生物与多尺度模型 / Phases 13–19

| Phase | 研究问题 | 报告 | 实现 | 数据 / 图表 | 证据边界 |
|---|---|---|---|---|---|
| [13](../projects/phase13/) | 金属酶与 PCET | [中文](../projects/phase13/reports/PCET_METALLOENZYME_REPORT_ZH.md) / [EN](../projects/phase13/reports/PCET_METALLOENZYME_REPORT_EN.md) | [代码](../projects/phase13/code/run_phase13_metalloenzyme_pcet_engine.py) | [数据](../projects/phase13/results) · [图表](../projects/phase13/figures) | 完整计算待验收 |
| [14](../projects/phase14/) | 活性物质与凝聚体 | [中文](../projects/phase14/reports/BIOMOLECULAR_CONDENSATE_REPORT_ZH.md) / [EN](../projects/phase14/reports/BIOMOLECULAR_CONDENSATE_REPORT_EN.md) | [代码](../projects/phase14/code/run_phase14_active_matter_condensate_phase_separation.py) | [数据](../projects/phase14/results) · [图表](../projects/phase14/figures) | 历史运行快照；待完整复核 |
| [15](../projects/phase15/) | 自由基对与变构模型 | [中文](../projects/phase15/reports/QUANTUM_BIOLOGY_REPORT_ZH.md) / [EN](../projects/phase15/reports/QUANTUM_BIOLOGY_REPORT_EN.md) | [代码](../projects/phase15/code/run_phase15_quantum_biology_spin_allostery.py) | [数据](../projects/phase15/results) · [图表](../projects/phase15/figures) | 完整生产与窗口统计待验收 |
| [16](../projects/phase16/) | 核孔复合物与输运 | [中文](../projects/phase16/reports/MEGAMACHINE_CRYOM_REPORT_ZH.md) / [EN](../projects/phase16/reports/MEGAMACHINE_CRYOM_REPORT_EN.md) | [代码](../projects/phase16/code/run_phase16_megamachine_cryoem_transport.py) | [数据](../projects/phase16/results) · [图表](../projects/phase16/figures) | 历史运行快照；待模块验收 |
| [17](../projects/phase17/) | 相对论量子化学 | [中文](../projects/phase17/reports/RELATIVISTIC_QUANTUM_REPORT_ZH.md) / [EN](../projects/phase17/reports/RELATIVISTIC_QUANTUM_REPORT_EN.md) | [代码](../projects/phase17/code/run_phase17_relativistic_actinide_quantum.py) | [数据](../projects/phase17/results) · [图表](../projects/phase17/figures) | 指定检查通过；有限基组/原子叠加 |
| [18](../projects/phase18/) | 全细胞代谢与热力学 | [中文](../projects/phase18/reports/WHOLE_CELL_METABOLISM_REPORT_ZH.md) / [EN](../projects/phase18/reports/WHOLE_CELL_METABOLISM_REPORT_EN.md) | [代码](../projects/phase18/code/run_phase18_wholecell_metabolic_thermodynamics.py) | [数据](../projects/phase18/results) · [图表](../projects/phase18/figures) | 所列数值检查通过；非实验验证 |
| [19](../projects/phase19/) | 酶设计与几何闭环 | [中文](../projects/phase19/reports/ACTIVE_INFERENCE_ENZYME_REPORT_ZH.md) / [EN](../projects/phase19/reports/ACTIVE_INFERENCE_ENZYME_REPORT_EN.md) | [代码](../projects/phase19/code/run_phase19_active_inference_denovo_enzyme.py) | [数据](../projects/phase19/results) · [图表](../projects/phase19/figures) | 几何测试不等于完整酶演化验证 |

<a id="devices"></a>

## 量子与分子器件模型 / Phases 20–23

| Phase | 研究问题 | 报告 | 实现 | 数据 / 图表 | 证据边界 |
|---|---|---|---|---|---|
| [20](../projects/phase20/) | CISS 输运与自旋电子学 | [中文](../projects/phase20/reports/CISS_QUANTUM_SPINTRONICS_REPORT_ZH.md) / [EN](../projects/phase20/reports/CISS_QUANTUM_SPINTRONICS_REPORT_EN.md) | [代码](../projects/phase20/code/run_phase20_ciss_quantum_spintronics.py) | [数据](../projects/phase20/results) · [图表](../projects/phase20/figures) | 有效输运模型；OER 与分析信号有假设 |
| [21](../projects/phase21/) | 腔 QED 与极化激元 | [中文](../projects/phase21/reports/POLARITON_CHEMISTRY_REPORT_ZH.md) / [EN](../projects/phase21/reports/POLARITON_CHEMISTRY_REPORT_EN.md) | [代码](../projects/phase21/code/run_phase21_cavity_qed_polaritonic_chemistry.py) | [数据](../projects/phase21/results) · [图表](../projects/phase21/figures) | 谐振 Pauli–Fierz；光学分裂不证明催化 |
| [22](../projects/phase22/) | 分子自旋量子比特 | [中文](../projects/phase22/reports/MOLECULAR_SPIN_QUBIT_REPORT_ZH.md) / [EN](../projects/phase22/reports/MOLECULAR_SPIN_QUBIT_REPORT_EN.md) | [代码](../projects/phase22/code/run_phase22_molecular_spin_qubits.py) | [数据](../projects/phase22/results) · [图表](../projects/phase22/figures) | 指定自旋模型与合成浴；非实测 |
| [23](../projects/phase23/) | 分子忆阻器与储备池计算 | [中文](../projects/phase23/reports/IN_MATERIO_COMPUTING_REPORT_ZH.md) / [EN](../projects/phase23/reports/IN_MATERIO_COMPUTING_REPORT_EN.md) | [代码](../projects/phase23/code/run_phase23_in_materio_neuromorphic_computing.py) | [数据](../projects/phase23/results) · [图表](../projects/phase23/figures) | 有效电路模型；原定性能目标未达 |

<a id="catalysis"></a>

## 催化研究试点 / Phases 24–27

| Phase | 研究问题 | 报告 | 实现 | 数据 / 图表 | 证据边界 |
|---|---|---|---|---|---|
| [24](../projects/phase24/) | AI 高通量催化模拟 | [中文技术报告](../projects/phase24/reports/PHASE24_GXNU_AI_HTS_REPORT_ZH.md) / [EN brief](../projects/phase24/reports/GXNU_AI_HTS_PLATFORM_PROPOSAL_EN.md) | [代码](../projects/phase24/code/run_ai_hts_platform_pilot.py) | [数据](../projects/phase24/results) · [图表](../projects/phase24/figures) | 合成反应响应；未投运实验平台 |
| [25](../projects/phase25/) | 不对称催化 | [中文](../projects/phase25/reports/PHASE25_REPORT_ZH.md) / [EN](../projects/phase25/reports/PHASE25_REPORT_EN.md) · [进展](../projects/phase25/reports/PHASE25_PROGRESS_ZH.md) | [模块指南](../shared/phase25_27/code/phase25_27/README.md) | [数据](../projects/phase25/results) | xTB/DFT 初步计算；选择性未判定 |
| [26](../projects/phase26/) | 恒电位铜界面 | [中文](../projects/phase26/reports/PHASE26_REPORT_ZH.md) / [EN](../projects/phase26/reports/PHASE26_REPORT_EN.md) | [模块指南](../shared/phase25_27/code/phase25_27/README.md) | [数据](../shared/phase25_27/results) | 设计与干表面种子；无完整界面采样 |
| [27](../projects/phase27/) | Ni 光动力学 | [中文](../projects/phase27/reports/PHASE27_REPORT_ZH.md) / [EN](../projects/phase27/reports/PHASE27_REPORT_EN.md) | [模块指南](../shared/phase25_27/code/phase25_27/README.md) | [数据](../shared/phase25_27/results) | 文献坐标与方案；无新非绝热轨迹 |

## 产业问题、公开数据与系统原型 / Phases 28–30

| Phase | 研究问题 | 任务与报告 | 实现 | 数据 / 图表 | 证据边界 |
|---|---|---|---|---|---|
| [28](../projects/phase28/) | 工业氯碱有机电极、失活识别与寿命成本 | [任务书](../projects/phase28/TASK_PROMPT.md) · [报告](../projects/phase28/reports/) | [代码](../projects/phase28/code/) | [数据](../projects/phase28/data/) · [结果](../projects/phase28/results/) · [图表](../projects/phase28/figures/) | 公开源数据复算与假设模型；未取得电极寿命、启停或工厂数据 |
| [29](../projects/phase29/) | 医药吡啶位点控制、竞争机理与可制造性 | [任务书](../projects/phase29/TASK_PROMPT.md) · [报告](../projects/phase29/reports/) | [代码](../projects/phase29/code/) | [数据](../projects/phase29/data/) · [结果](../projects/phase29/results/) · [图表](../projects/phase29/figures/) | 文献反应审计、分子计算与假设模型；尚未证明同底物同试剂 C2/C4 切换 |
| [30](../projects/phase30/) | CO₂RR 微观到工业堆数字孪生 | [完整规格](../projects/phase30/reports/phase30_technical_report.md) | [代码](../projects/phase30/code/) | [结果](../projects/phase30/results/) · [图表](../projects/phase30/figures/) | 未校准降阶原型；无DFT/AIMD、3D CFD或工业验证 |

## 跨阶段审计入口

新增 [Phase31](../projects/phase31/)：大环分子胶 QM/MM–FEP 七角色系统规格与本地工程试验。[技术报告](../projects/phase31/reports/phase31_scientific_full.md) · [真实运行边界](../projects/phase31/reports/local_validation.md)。已有 BRD4–VHL 公开结构与 MMFF 软件样例；没有执行第一性原理 QM/MM–FEP、IRC 或发现新先导。

| 范围 | 更正与来源 |
|---|---|
| 1–19 | [收尾清单](../shared/audit/reports/PHASE1_19_CLOSEOUT.md) · [有时间戳的重算表](../shared/audit/reports/PHASE1_19_RERUN_STATUS.md) · [公开证据](../shared/audit/evidence/phase1_19_rerun_20260910/README.md) |
| 1–5 | [审计更正](../shared/audit/reports/PHASE1_5_AUDIT_FIXES.md) |
| 6–10 | [审计更正](../shared/audit/reports/PHASE6_10_AUDIT_FIXES.md) |
| 11–13 | [审计更正](../shared/audit/reports/PHASE11_13_AUDIT_FIXES.md) |
| 14–17 | [审计更正](../shared/audit/reports/PHASE14_17_AUDIT_FIXES.md) |
| 18–19 | [审计更正](../shared/audit/reports/PHASE18_19_AUDIT_FIXES.md) |
| [24](../projects/phase24/) | [验证范围](../projects/phase24/reports/PHASE24_VALIDATION.md) |
| 25–27 | [来源与结果](../shared/phase25_27/results/README.md) · [逐条验收](../shared/phase25_27/results/acceptance.json) |

需要原有详细介绍和完整图表串览时，访问[历史首页](history/README_HISTORY.md)。
