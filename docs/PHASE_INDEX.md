# 阶段导航 / Phase index

[返回首页](../README.md) · [证据说明](EVIDENCE.md) · [复现指南](REPRODUCE.md)

按研究主题查找全部 27 个阶段。表内“边界”来自已提交的报告、[Phase 1–19 收尾清单](../PHASE1_19_CLOSEOUT.md)和[25–27 验收台账](../results_phase25_27/acceptance.json)，不是本次重新计算的结果。
“数据”可能包含历史输出；较新的审计与重算说明优先于早期图表。Phase 编号是模块标识，不表示成熟度。

<a id="molecules"></a>

## 分子与反应 / Phases 1–8

| Phase | 研究问题 | 报告 | 实现 | 数据 / 图表 | 证据边界 |
|---|---|---|---|---|---|
| 1 | 复杂骨架与构象 | [中文](../BENCHMARK_REPORT_ZH.md) / [EN](../BENCHMARK_REPORT_EN.md) | [代码](../molecule_benchmark.py) | [数据](../bench_results/) · [图表](../figures/) | 所列检查通过；保留输入与适用性局限 |
| 2 | 扭转扫描与分子动力学 | [中文](../DYNAMICS_REPORT_ZH.md) / [EN](../DYNAMICS_REPORT_EN.md) | [代码](../run_heavy_dynamics_benchmark.py) | [数据](../results_phase2/) · [图表](../figures_phase2/) | 所列检查通过；注意替代力场与电荷 |
| 3 | KRAS 复合物与结合分析 | [中文](../COMPLEX_DYNAMICS_REPORT_ZH.md) / [EN](../COMPLEX_DYNAMICS_REPORT_EN.md) | [代码](../run_phase3_complex_dynamics.py) | [数据](../results_phase3/) · [图表](../figures_phase3/) | 完整 MD / MM-GBSA 待验收 |
| 4 | 骨架编辑与反应路径 | [中文](../SKELETAL_EDITING_REPORT_ZH.md) / [EN](../SKELETAL_EDITING_REPORT_EN.md) | [代码](../run_phase4_reaction_mechanism.py) | [数据](../results_phase4/) · [图表](../figures_phase4/) | 独立精修驻点/单虚频通过；IRC 未验证 |
| 5 | 反应网络与逆向设计 | [中文](../WORLD_MODEL_REPORT_ZH.md) / [EN](../WORLD_MODEL_REPORT_EN.md) | [代码](../run_phase5_chemical_world_model.py) | [数据](../results_phase5/) · [图表](../figures_phase5/) | 多虚频 TS 不合格；近零产率保留 |
| 6 | 显式溶剂与增强采样 | [中文](../METADYNAMICS_REPORT_ZH.md) / [EN](../METADYNAMICS_REPORT_EN.md) | [代码](../run_phase6_explicit_metadynamics.py) | [数据](../results_phase6/) · [图表](../figures_phase6/) | 轨迹完成但产物采样不足 |
| 7 | 强关联与模型适用边界 | [中文](../FRONTIER_EPISTEMIC_REPORT_ZH.md) / [EN](../FRONTIER_EPISTEMIC_REPORT_EN.md) | [代码](../run_phase7_strong_correlation_wall.py) | [数据](../results_phase7/) · [图表](../figures_phase7/) | 全点/统一基组复核待完成 |
| 8 | 光化学与非绝热动力学 | [中文](../PHOTOCHEMISTRY_REPORT_ZH.md) / [EN](../PHOTOCHEMISTRY_REPORT_EN.md) | [代码](../run_phase8_photochemical_dynamics.py) | [数据](../results_phase8/) · [图表](../figures_phase8/) | 完整 QC 与动力学待验收 |

<a id="automation"></a>

## 自动化与定律发现 / Phases 9–12

| Phase | 研究问题 | 报告 | 实现 | 数据 / 图表 | 证据边界 |
|---|---|---|---|---|---|
| 9 | 自驱动实验室模拟 | [中文](../SELF_DRIVING_LAB_REPORT_ZH.md) / [EN](../SELF_DRIVING_LAB_REPORT_EN.md) | [代码](../run_phase9_self_driving_lab_compiler.py) | [数据](../results_phase9/) · [图表](../figures_phase9/) | 模拟与协议检查通过；无硬件验证 |
| 10 | 连续流反应器数字孪生 | [中文](../FLOW_CYBERPHYSICAL_REPORT_ZH.md) / [EN](../FLOW_CYBERPHYSICAL_REPORT_EN.md) | [代码](../run_phase10_cyberphysical_flow_twin.py) | [数据](../results_phase10/) · [图表](../figures_phase10/) | 所列检查通过；控制模型局限 |
| 11 | 神经波函数与 VMC | [中文](../NEURAL_WAVEFUNCTION_REPORT_ZH.md) / [EN](../NEURAL_WAVEFUNCTION_REPORT_EN.md) | [代码](../run_phase11_neural_wavefunction_vmc.py) | [数据](../results_phase11/) · [图表](../figures_phase11/) | 完整计算结束；He 精度目标未达 |
| 12 | 动力学定律发现 | [中文](../SCIENTIFIC_AGI_MANIFESTO_ZH.md) / [EN](../SCIENTIFIC_AGI_MANIFESTO_EN.md) | [代码](../run_phase12_hamiltonian_law_discovery.py) | [数据](../results_phase12/) · [图表](../figures_phase12/) | 完整训练与验收待完成 |

<a id="biological"></a>

## 生物与多尺度模型 / Phases 13–19

| Phase | 研究问题 | 报告 | 实现 | 数据 / 图表 | 证据边界 |
|---|---|---|---|---|---|
| 13 | 金属酶与 PCET | [中文](../PCET_METALLOENZYME_REPORT_ZH.md) / [EN](../PCET_METALLOENZYME_REPORT_EN.md) | [代码](../run_phase13_metalloenzyme_pcet_engine.py) | [数据](../results_phase13/) · [图表](../figures_phase13/) | 完整计算待验收 |
| 14 | 活性物质与凝聚体 | [中文](../BIOMOLECULAR_CONDENSATE_REPORT_ZH.md) / [EN](../BIOMOLECULAR_CONDENSATE_REPORT_EN.md) | [代码](../run_phase14_active_matter_condensate_phase_separation.py) | [数据](../results_phase14/) · [图表](../figures_phase14/) | 历史运行快照；待完整复核 |
| 15 | 自由基对与变构模型 | [中文](../QUANTUM_BIOLOGY_REPORT_ZH.md) / [EN](../QUANTUM_BIOLOGY_REPORT_EN.md) | [代码](../run_phase15_quantum_biology_spin_allostery.py) | [数据](../results_phase15/) · [图表](../figures_phase15/) | 完整生产与窗口统计待验收 |
| 16 | 核孔复合物与输运 | [中文](../MEGAMACHINE_CRYOM_REPORT_ZH.md) / [EN](../MEGAMACHINE_CRYOM_REPORT_EN.md) | [代码](../run_phase16_megamachine_cryoem_transport.py) | [数据](../results_phase16/) · [图表](../figures_phase16/) | 历史运行快照；待模块验收 |
| 17 | 相对论量子化学 | [中文](../RELATIVISTIC_QUANTUM_REPORT_ZH.md) / [EN](../RELATIVISTIC_QUANTUM_REPORT_EN.md) | [代码](../run_phase17_relativistic_actinide_quantum.py) | [数据](../results_phase17/) · [图表](../figures_phase17/) | 指定检查通过；有限基组/原子叠加 |
| 18 | 全细胞代谢与热力学 | [中文](../WHOLE_CELL_METABOLISM_REPORT_ZH.md) / [EN](../WHOLE_CELL_METABOLISM_REPORT_EN.md) | [代码](../run_phase18_wholecell_metabolic_thermodynamics.py) | [数据](../results_phase18/) · [图表](../figures_phase18/) | 所列数值检查通过；非实验验证 |
| 19 | 酶设计与几何闭环 | [中文](../ACTIVE_INFERENCE_ENZYME_REPORT_ZH.md) / [EN](../ACTIVE_INFERENCE_ENZYME_REPORT_EN.md) | [代码](../run_phase19_active_inference_denovo_enzyme.py) | [数据](../results_phase19/) · [图表](../figures_phase19/) | 几何测试不等于完整酶演化验证 |

<a id="devices"></a>

## 量子与分子器件模型 / Phases 20–23

| Phase | 研究问题 | 报告 | 实现 | 数据 / 图表 | 证据边界 |
|---|---|---|---|---|---|
| 20 | CISS 输运与自旋电子学 | [中文](../CISS_QUANTUM_SPINTRONICS_REPORT_ZH.md) / [EN](../CISS_QUANTUM_SPINTRONICS_REPORT_EN.md) | [代码](../run_phase20_ciss_quantum_spintronics.py) | [数据](../results_phase20/) · [图表](../figures_phase20/) | 有效输运模型；OER 与分析信号有假设 |
| 21 | 腔 QED 与极化激元 | [中文](../POLARITON_CHEMISTRY_REPORT_ZH.md) / [EN](../POLARITON_CHEMISTRY_REPORT_EN.md) | [代码](../run_phase21_cavity_qed_polaritonic_chemistry.py) | [数据](../results_phase21/) · [图表](../figures_phase21/) | 谐振 Pauli–Fierz；光学分裂不证明催化 |
| 22 | 分子自旋量子比特 | [中文](../MOLECULAR_SPIN_QUBIT_REPORT_ZH.md) / [EN](../MOLECULAR_SPIN_QUBIT_REPORT_EN.md) | [代码](../run_phase22_molecular_spin_qubits.py) | [数据](../results_phase22/) · [图表](../figures_phase22/) | 指定自旋模型与合成浴；非实测 |
| 23 | 分子忆阻器与储备池计算 | [中文](../IN_MATERIO_COMPUTING_REPORT_ZH.md) / [EN](../IN_MATERIO_COMPUTING_REPORT_EN.md) | [代码](../run_phase23_in_materio_neuromorphic_computing.py) | [数据](../results_phase23/) · [图表](../figures_phase23/) | 有效电路模型；原定性能目标未达 |

<a id="catalysis"></a>

## 催化研究试点 / Phases 24–27

| Phase | 研究问题 | 报告 | 实现 | 数据 / 图表 | 证据边界 |
|---|---|---|---|---|---|
| 24 | AI 高通量催化模拟 | [中文技术报告](../PHASE24_GXNU_AI_HTS_REPORT_ZH.md) / [EN brief](../GXNU_AI_HTS_PLATFORM_PROPOSAL_EN.md) | [代码](../run_ai_hts_platform_pilot.py) | [数据](../results_phase24/) · [图表](../figures_hts_pilot/) | 合成反应响应；未投运实验平台 |
| 25 | 不对称催化 | [中文](../PHASE25_REPORT_ZH.md) / [EN](../PHASE25_REPORT_EN.md) · [进展](../PHASE25_PROGRESS_ZH.md) | [模块指南](../phase25_27/README.md) | [数据](../results_phase25_27/phase25/) | xTB/DFT 初步计算；选择性未判定 |
| 26 | 恒电位铜界面 | [中文](../PHASE26_REPORT_ZH.md) / [EN](../PHASE26_REPORT_EN.md) | [模块指南](../phase25_27/README.md) | [数据](../results_phase25_27/) | 设计与干表面种子；无完整界面采样 |
| 27 | Ni 光动力学 | [中文](../PHASE27_REPORT_ZH.md) / [EN](../PHASE27_REPORT_EN.md) | [模块指南](../phase25_27/README.md) | [数据](../results_phase25_27/) | 文献坐标与方案；无新非绝热轨迹 |

## 跨阶段审计入口

| 范围 | 更正与来源 |
|---|---|
| 1–19 | [收尾清单](../PHASE1_19_CLOSEOUT.md) · [有时间戳的重算表](../PHASE1_19_RERUN_STATUS.md) · [公开证据](../audit/phase1_19_rerun_20260910/README.md) |
| 1–5 | [审计更正](../PHASE1_5_AUDIT_FIXES.md) |
| 6–10 | [审计更正](../PHASE6_10_AUDIT_FIXES.md) |
| 11–13 | [审计更正](../PHASE11_13_AUDIT_FIXES.md) |
| 14–17 | [审计更正](../PHASE14_17_AUDIT_FIXES.md) |
| 18–19 | [审计更正](../PHASE18_19_AUDIT_FIXES.md) |
| 24 | [验证范围](../PHASE24_VALIDATION.md) |
| 25–27 | [来源与结果](../results_phase25_27/README.md) · [逐条验收](../results_phase25_27/acceptance.json) |

需要原有详细介绍和完整图表串览时，访问[历史首页](../README_HISTORY.md)。
