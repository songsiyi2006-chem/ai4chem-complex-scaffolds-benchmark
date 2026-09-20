# Phase 28–29：第二轮补充计算

[首页](../README.md) · [首轮交付记录](PHASE28_29_STATUS.md) · [Phase 28](../projects/phase28/) · [Phase 29](../projects/phase29/)

2026-09-20。按用户“补充更多计算，然后继续推送”的要求，在首轮 `73b1ea0` 基础上完成以下增量。范围仍为公开数据和已有本机 CPU，未新增实验或使用 HPC。

| 阶段 | 本轮实际执行 | 报告与原始结果 |
|---|---|---|
| 28 | 新提取 50 条已发表均值/标准差记录和 777 个极化坐标；50 个小样本均值区间、48 组协方差敏感性对比、28 个共电流插值、72 组时间窗口检查 | [中文报告](../projects/phase28/reports/ADDITIONAL_PUBLIC_ANALYSIS_ZH.md) · [English](../projects/phase28/reports/ADDITIONAL_PUBLIC_ANALYSIS_EN.md) · [结果](../projects/phase28/results/additional_public_analysis/) |
| 29 | 12 个分子 × 2 个电荷/自旋状态 × 2 个 GFN 方法 × 3 个环境，实际运行 144 次 xTB 单点；其中 12 次重复首轮设置用于衔接复核 | [中文报告](../projects/phase29/reports/xtb_matrix_zh.md) · [English](../projects/phase29/reports/xtb_matrix_en.md) · [结果](../projects/phase29/results/xtb_matrix/) |
| 29 | 3 个分子 PBE0/def2-SVP 的 6 次单点，以及母体 def2-TZVP 的 2 次检查；8 次 SCF 收敛；匹配气相几何并比较中心化取代效应 | [中文报告](../projects/phase29/reports/dft_crosscheck_zh.md) · [English](../projects/phase29/reports/dft_crosscheck_en.md) · [结果](../projects/phase29/results/dft_crosscheck/) |

## 对结论有什么影响

Phase 28 的旧 +44/+17 mV 曲线端点差对早期时间窗口敏感。按分析者设定排除至少前 10 h 后，所测试窗口的范围分别缩小为 −4 至 +9.25 mV、+3 至 +6 mV。此处没有用这一规则证明物理稳态。24 h Cl₂ 均值为 0.338190 mol，条件性 95% 小样本均值区间为 0.326924–0.349456 mol；同次运行的跨时间协方差缺失，使 FE 变化的区间结论依赖假设。100/200 h 极化坐标差随电流变号，不能称作整体性能提升。

Phase 29 的 xTB 矩阵全部 SCC 收敛且正常结束，同时保留全部 144 次数值警告。GFN 方法在不同环境中有 4–6/66 对分子排序反转，同方法气相到乙腈有 7/66 对反转。原始跨方法电荷态差还含未对齐的电子参考，不能据此量化预测误差。以母体为基准后，PBE0、GFN1 和 GFN2 对甲氧基/氰基取代效应的方向一致，但氰基幅度仍有模型差异。这些结果没有给出还原电位、能垒或 C2/C4 选择性。

## 复现与证据

新增代码、配置、数值输入、原始引擎输出、双语报告及 5 张新图均位于相应阶段。首轮科学输出保留；首轮验收/哈希清单另存历史快照，当前清单覆盖新交付。扩展运行命令见各补充报告。

本轮软件验证、公开数据独立复跑比对与哈希核对记录见 [additional_compute_validation/summary.json](additional_compute_validation/summary.json)。软件测试通过、SCF 收敛、图表检查属于实现与数值证据，不能代替实验机制或工业验收。
