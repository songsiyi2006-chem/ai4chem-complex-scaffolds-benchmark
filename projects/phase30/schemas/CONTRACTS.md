# Phase30 跨尺度接口契约 v1

内部数据均为有限浮点数；未知物理量用 JSON null 加原因，不写 NaN/Infinity。
编号文件通过显式 `importlib` 加载，避免 Python 标识符不能以数字开头的问题。

## 可执行函数

| 入口 | 输入 | 输出与单位 | 失败行为 |
|---|---|---|---|
| `01_micro_surface.packet(config)` | 温度K、电位V vs bulk RHE、pH、Cu/Ag、假设配体扰动eV | 假设势垒eV、有效速率s⁻¹、计量网络 | 超出激活模型域报错 |
| `02_meso_microkinetics.solve(config,c,validate_ode)` | 溶解CO₂ mol/m³ | 覆盖度、局部pH、ψ V、Γr mol/m²/s、j A/m²、FE无量纲 | pH无根、ODE失败或位点不守恒拒绝 |
| `03_macro_flow_cfd.channel(config)` | 通道SI几何、流速m/s、气体进料mol/s | 沿程量、产品mol/s、功率W、压降Pa、碳/热残差 | 非层流、超温/过流、物料不守恒拒绝 |
| `04_tea_lca_model.assess(config,flow,stacks)` | 完整通道结果与经济情景 | 折现USD/kg、未折现kgCO₂e/kg、年度事件与费用 | 无产品或无效时间参数拒绝 |
| `05_agent_closed_loop.optimize(config)` | 冻结预算/搜索域/seed | 已查询观测、拒绝理由、提案时训练ID、观测帕累托集 | 无可行点则失败，不造冠军 |

面积基准为几何电极面积；Γ是每平方米几何面积的有效位点摩尔数，隐含粗糙度假设。
不得把按ECSA归一的速率不经转换直接传入宏观面积模型。
CDF/GP方差只供选点，不代表 DFT 或实验误差条。

## 外部微观证据包：最低解析契约

`common.validate_micro_artifact(path)` 的职责仅为以下字段/路径/数值检查：

- `schema_version = phase30.micro.v1`；`evidence_type = computed | literature`。
- `method, functional, solvation, surface_id, geometry_sha256, convergence, temperature_K`。
- `potential_reference = SHE`；`records[]` 每条包含 `reaction_id, potential_SHE_V, barrier_eV, unit="eV"`。
- 非空 `raw_files[]`：各条 `path, sha256` 指向与输入包同目录下的原始证据，不允许逃出目录；逐文件校验。
- 返回 `schema_valid/raw_hashes_verified` 与始终为 false 的 `chemistry_validated/production_adapter_implemented`，不把文件完整性升级成化学结论。

该检查器不是完整 JSON Schema 验证器，也不自动解析任何量化引擎输出。
生产适配器实施时还必须验证：slab/吸附位点/覆盖度、赝势、k网格、真空与偶极修正、
电子化学势/电荷、显式水/离子采样、温度/自由能修正、同一势能面、TS虚频与连接性、
不确定性/相关性矩阵、原始几何和波函数/日志哈希及原文定位。当前最低契约不会检查这些科学事项。

## 未来宏观/实验适配器

**CFD field packet v1（设计，未实现）**：mesh SHA256、网格单位/体积、cell/face ID、面法向、
phase fraction、p、T、u、cᵢ、φₛ/φₗ、墙面活性面积及各反应通量；时间窗口及质量/电荷/热收支。
分别记录求解器版本、离散格式、时间步/CFL、收敛残差、网格加密和边界条件。局部通量必须积分再算堆级FE，不能直接平均各网格FE。

**Operando packet v1（设计，未实现）**：instrument/calibration ID、batch/electrode ID、UTC时间、
电位参考与iR补偿策略、I/V/T/压力/流量、GC/HPLC定量与校准曲线、检出限、测量协方差、
Raman/IR/XAS原始谱及处理版本；训练/保留批次隔离。代理只提出下一条件，硬件执行需独立授权。

## 状态机与证据清单

`PREPARED → RUNNING → COMPLETED_DEMO | FAILED`。完成示例同时保留 `scientific_acceptance=false`。
独占创建 execution.lock；源/配置/已冻结副本哈希必须匹配；拒绝重入失败或完成目录。
运行过程中源码改变则失败。JSON原子替换，异常写入 failure.txt。

`artifact_manifest.json` 包含运行目录所有文件的SHA256与长度，清单本身不自哈希。
项目发布版只复制 results/reports/figures，另用 example_provenance.json 明确选中范围。
外部结果不足、解析失败或物理门槛失败，不得调用示例参数作为隐藏后备。
