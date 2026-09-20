# Data card / 数据卡

版本 2026-09-20。原始实验条目数 **0**。所有科学标签携带 `evidence_layer` 或由独立文件夹规定层级，不混入训练。

| 目录/文件 | 内容与规模 | 允许用途 | 禁止推断 |
|---|---|---|---|
| `data/literature/condition_controls.csv` | JACS SI Table S1 15 条、NC SI Table S1 拆分后 25 条；40 条公开条件记录 | 来源审核、条件敏感性描述 | 不是独立重复；不作因果效应/ML 独立测试 |
| `natcomm_tableS3_energies.csv` | SI S19–S20 的 14 组 Gcorr/E/G | 单位与差分复核 | 不是本次重跑 DFT；无新的 TS/IRC 验证 |
| `reference_reaction_identity.csv` | 2 个经来源图形审核、下游手动重建重原子映射的参考反应 | 判断反应域是否一致 | 不是来源原生机器映射、也不验证反应机理 |
| `derived_electrical_metrics.csv` | JACS 标准、NC 标准及 3 mmol 放大；数值导出 | 明示电子数和分母的电量、FE 与通量 | 不跨不同产品比较工艺优越性 |
| `data/calculated/` | 12 底物及各 C2/C4/C6 提案：48 结构记录、151 MMFF 构象尝试 | 结构身份、描述符、可恢复几何 | 不是购得样品、产率、合成可行性或电极结构 |
| `data/synthetic/` | 5 假设动力学情景、120 合成家族、假设成本场景 | 流程、守恒、划分、校准算法检查 | 不作真实选择性、实验节省、ADME 或成本优势 |
| `data/experimental/reaction_template.csv` | 仅表头，无实验记录 | 实验交接模板 | 不得把 null/空表当作成功或失败反应 |
| `results/quantum_preflight/` | 12 次 GFN1/GFN2 固定几何单点，原始日志 | 分子层级方法敏感性负结果 | 不作实际界面、电位、动力学或位点选择性证据 |

## 结构和缺失值策略

环编号固定为 N1–C2–C3–C4–C5–C6，所有底物使用同一亲本方向，记录所有可用位点。P01 的 C2 与 C6 对称等价，其映射不是两个独立产品；统计总区域产率时合并同构结构。其余 3-取代底物必须区分 C2 与 C6。本文的动力学三个通道只代表一般非对称底物，不直接适用于 P01。反应映射仅守恒重原子并检查分子式；隐式氢不提供同位素转移来源证据。

新生醇位立体化学未指定；没有对映选择性结论。文件表示中性游离碱，真实酸性电解液中的质子化比例、多位点结合与盐形式未知。源条件中电解质和酸不因为未进入重原子产物映射而被认为不存在。

`Trace`、`N.R.` 保留原字符串，数值为 null；不把它们填成 0。没有原始峰面积、响应因子、LOD/LOQ、重复数或误差条，不能补造；GC-MS 相同质量数不独立证明 C2/C4 身份。NC SI Table S1 表头/正文 C1 与脚注 C1' 的冲突已单列，未静默改写原来源。JACS Table S1 各行产率口径未在表内逐项明确，保留 `reported_yield_table_basis_unresolved`。

## 划分与部署

拟议面板的 P08/P09（strong_EWG）及 P12（aryl）是冻结的取代基家族保留集，不是已完成的前瞻盲测。所有面板都基于吡啶，**不称为独立新骨架外推**。近邻家族、同底物不同条件、同组异构体、材料批次与文献来源须绑定；新增跨骨架挑战集必须另行冻结。

合成学习集仅有 120 个独立人造家族，66/24/30 用于训练/校准/测试。预处理仅在训练集拟合，固定超参数，校准残差只来自校准集。泛化误差只针对人造分布；真实部署一律 `ABSTAIN_ALL_REAL_REACTIONS_NO_EXPERIMENTAL_LABELS`。有限样本 conformal 的 90% 目标覆盖只有在交换性等假设下适用；本次未证明真实反应满足这些假设。

## 来源与使用

[来源表](../data/literature/literature.csv) 和 [下载哈希](../data/literature/source_download_provenance.json) 记录 URL、访问范围与证据。公开 SI/坐标表在仓库外 scratch 中检查，仓库仅提交少量结构化事实和导出结果，不镜像全文 PDF。JACS SI 可通过 [ACS Figshare](https://ndownloader.figshare.com/files/65594325) 获取；遵循各发布者的来源许可，不把工具许可当成数据许可。
