# Phase31 技术规格书与科学证据审计

版本：2026-09-21；性质：架构设计、可执行骨架与本地试验，不是完成态论文。

## 1. 科学问题与系统选择

### 1.1 可证伪命题

H1：在同一 BRD4 BD2–VHL 结合机制及明确标准态下，改变大环的构象限制可改善条件结合自由能；改善不是简单的气相构象势能降低。

H2：候选的构象预组织收益能抵消结合姿态张力和去溶剂化代价。反例包括预组织不足、过度张力、溶解性/渗透性损失和新闭环造成的结合位点冲突。

H3：经独立合成和生物物理验证的三元协同性可以指导实验优先级，但不能单独预测细胞降解速率。

这些命题均未在本项目中证实。目标从“堆叠所有方法”改为“逐道证据门检验”：若基线不能复现，停止扩大候选数。

### 1.2 用户授权后选定的首个体系

选择人源 BRD4 第二溴结构域作为 Target，VHL 作为 E3 底物识别亚基，Elongin B/C 作为结构组件；macroPROTAC-1 为非共价宏环双功能基线。公开结构 6SIS 分辨率 3.5 Å，配体标识 LFE；不包含足以直接计算完整泛素化反应的 E2–泛素–全 E3 工作体系。[RCSB 6SIS](https://www.rcsb.org/structure/6SIS)

此选择利用已发表的三元复合物作为方法校验，而非宣称 BRD4 “不可成药”。宏环 PROTAC 与 molecular glue 在分子识别方式上不是同义词。“高张力”“新分子胶”“共价协同”保留为待证明的后续研究分支。

原始研究提供了宏环化、结合测量及晶体结构的先例；其成功不外推到本项目未知候选。[Testa 等原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7004083/)

### 1.3 数据与准备责任

已下载 6SIS.pdb 与 LFE_ideal.sdf，原始字节、URL、时间和 SHA256 记录于 `results/reference/provenance.json`。PDB 非对称单元含两套蛋白/配体实例；暂拟 author chains A/B/C/D，对应 BRD4/ElonginB/ElonginC/VHL。这是准备清单，不是自动确认 biological assembly。

生产前必须人工检查：生物组装、晶体接触、缺失残基/侧链、替代构象、占有率、配体键级/质子化/手性、结晶添加剂、水分子可信度、实验 pH/离子强度、端基处理。3.5 Å 不能支持对所有氢键和水桥的原子级确定性判断。

LFE 理想 SDF 是化学组分理想坐标，不是晶体结合构象；必须通过原子命名/拓扑映射转移坐标，不能把两个文件按行号拼接。

共价位点设为 null；不能为了满足标题任意添加亲电弹头或选择蛋白亲核残基。真实共价分支须另立机制、结构与实验依据。

## 2. 证据等级与验收门

| 门 | 必需证据 | 当前交付 |
|---|---|---|
| G0 来源与身份 | 原始结构哈希、化学拓扑、手性、实验来源 | 原始文件与初步身份清单已具备；准备待完成 |
| G1 工程可靠性 | 单元测试、原子领取、预算幂等、恢复、运行日志 | 本机已测试；范围见执行记录 |
| G2 构象与参数 | 参数覆盖、构象去重、异构保持、溶液多起点采样 | 仅小型 MMFF 软件样例 |
| G3 QM | 原始指定方法、收敛、稳定性、自旋、基组与活性空间敏感性 | 未执行，方法兼容性阻断 |
| G4 反应路径 | 端点极小值、收敛 NEB、鞍点 Hessian、双向 IRC | 仅接口/审计契约 |
| G5 QM/MM–FEP | 力梯度检验、周期电场、拓扑、采样和闭合热力学循环 | 未实现生产驱动，未采样 |
| G6 协同预测 | 同一结合腿、独立重复、不确定度及外部基线 | 分析函数通过解析测试，无真实新预测 |
| G7 实验 | 纯度/结构、结合、三元协同、泛素化、细胞降解/毒性 | 无新增实验 |

SUCCEEDED 指角色任务成功结束，不等于通过 G3–G7。BLOCKED 是正确保留缺失证据，不是用模拟数据填补的错误。软件测试数不能作为科学完成百分比。

证据标签固定：SPECIFICATION_ONLY、PUBLIC_STRUCTURE_NOT_PREPARED_SIMULATION、EXECUTED_MMFF_BENCHMARK、EXECUTED_QM_NOT_BINDING_FREE_ENERGY、NO_QUALIFYING_EVIDENCE。当前不存在真实采样/实验标签。

## 3. 目录与实现分层

```text
run.py -> prepare/freeze -> SQLite DAG -> isolated attempt artifacts -> Omega -> snapshots/report
                     |                 |
              source/input hashes      +-- Alpha: RDKit pilot
                     |                 +-- Beta/Gamma: QM contract + bounded healing
                     |                 +-- Delta: lambda contract + analysis, sampling blocked
                     |                 +-- Epsilon: route DAG audit
                     |                 +-- Zeta: conditional thermodynamics
                     +-- token reservations/actual provider receipts
```

编号模块经 importlib 显式加载，避免以数字开头的文件不能普通 import。CLI 的 `--out` 与仓库 `tools/run_phase.py 31 --workspace ...` 兼容。默认 plan 不进行硬计算；只有 local-pilot 运行受限构象试验。

## 4. 七角色蜂群：当前确定性实现与未来接口

| 角色 | 输入 | 输出/否决权 | 全局上限分配 |
|---|---|---|---:|
| Alpha | 原始结构哈希、研究假说、拓扑约束 | 候选身份与构象；本地只有软件样例 | 10,000,000 |
| Beta | 冻结的 Hamiltonian 与几何 | 电子结构和 TS 证据；缺引擎即阻断 | 25,000,000 |
| Gamma | 完整原始日志的索引与结构化错误 | 最多 5 次数值收敛策略；不可改化学问题 | 15,000,000 |
| Delta | 体系准备与参数清单 | QM/MM、采样计划及验收指标 | 25,000,000 |
| Epsilon | 有文献来源的路线图 | 结构规则检查；化学可行性不凭空打分 | 5,000,000 |
| Zeta | 成对自由能/解离常数及协方差 | α 和不确定度；不代替泛素化实验 | 10,000,000 |
| Omega | 全部原始证据与失败记录 | 阻断夸大结论、导出清单 | 10,000,000 |

角色提示词是未来 LLM 提供方的契约。当前 role handlers 是确定性 Python，不是七个实际联网模型。没有持久化“虚构辩论”或人为制造 Token 账单。50–100 个候选的细粒度矩阵调度仍属扩展项。

依赖图为 Alpha→Beta→Gamma，Alpha→Delta，Alpha→Epsilon，Beta/Delta→Zeta，前六者→Omega。允许父任务 BLOCKED 后审计角色执行，防止科学失败令证据日志失联；下游不得把 blocked 父结果当成数值。

### 4.1 状态机与恢复

PENDING→RUNNING→SUCCEEDED/BLOCKED/FAILED。事务领取使用 BEGIN IMMEDIATE；每次领取增加 fence。过期作业可重新领取，但旧 owner/fence 的结果被拒绝。每次 attempt 使用独立目录，失败输出不覆盖。完成任务不会在 resume 中再次运行。

工作区冻结源代码/提示词哈希和引用的结构哈希。修改代码后要求新工作区，而不是继续以原运行身份写证据。数据库 hash chain 可检测正常流程外的记录破坏，但不是第三方不可篡改公证；仍须外部签名或可信存储做生产审计。

当前 SQLite 只支持同机多进程；未实现 NFS、多主节点调度、失联 HPC job adoption、持续 heartbeat、跨机 GPU 资源控制。固定 600 s 租约只用于此轻量角色试验，不能直接拿去运行长 QM 作业。

### 4.2 集群扩展规格（未实现）

生产环境用 PostgreSQL 事务队列或受控 RPC 服务替代共享文件数据库。每个任务绑定 candidate_id、stage、immutable_input_hash、attempt_id、scheduler_job_id、resource_request、lease/fence、artifact_URI。Slurm adapter 仅允许已审计的容器/脚本模板，使用 job array；作业请求 CPU/GPU/内存/墙钟/磁盘配额，提交后保存调度器回执。

服务端 lease heartbeat 与计算作业身份分离；超时先向调度器核实状态，不能仅因控制端失联重复提交。重启后以 job_id 认领在运行作业；原子发布结果+校验值；晚到结果隔离而非覆盖。数据库日志和对象存储存 raw artifacts，不把整段 PDB 每轮塞入语言模型上下文。

### 4.3 Token 经济学

100,000,000 是可配置**上限**，不是保守最低消耗，更不等价于 CPU/GPU 小时。科学推进按证据改善衡量；上下文去重、日志索引和缓存应减少费用。不要用故意低效沟通证明“硬核”。

每个外部调用必须先按最坏输入+输出额度原子 reserve，使用提供方幂等键；真实 usage 回执 settle。cached input 已包括在 input 内，不重复累计；reasoning tokens 若提供方已计入 output，也不再次累加。未知 usage 保留预留额，不能擅自释放；超额回执全额入账，阻断后续超预算调用。唯一 receipt 防止一次调用记两次账。

当前提供方适配仅有 Protocol 和调用边界，没有实现供应商认证或真实网络调用。数据库计数 0 表示**管线未调用提供方**，不统计编写代码的本次 Codex 会话。

并发上限、模型价格、费用上限与费用实际回执还需生产部署另行配置；Token 总量不足以限制货币成本。

## 5. 模块 01：构象、手性与张力

本机执行 3 个 12–14 元手性内酰胺软件样例，各最多 6 个 ETKDGv3 构象。固定 seed、单线程、macrocycle torsions/14 configuration、enforceChirality，随后 MMFF94s 最小化。保留未收敛标记而非删掉失败；只为收敛构象计算同分子相对能量。[RDKit 官方说明](https://www.rdkit.org/docs/GettingStartedInPython.html)

这 3 个样例不是 macroPROTAC-1 衍生物，不是 50–100 个新靶向候选，未证明高张力或结合能力。MMFF 能量只作为低成本几何检查，不能跨不同分子比较总能量。

未来候选研究需保留立体异构/互变异构/质子化态唯一标识、3D 立体一致性、构象聚类、多种随机起点、溶剂环境与参数覆盖。高张力定义要明确参照反应：同构/同键反应或受控的 bound-vs-relaxed distortion；后者仍不自动等于绝对环张力。

## 6. 模块 02：QM、电子态与方法兼容性

### 6.1 指定方法不能偷偷替换

本次核对的 PySCF 官方 dispersion 实现把原始 ωB97X-D 列入不支持集合。ωB97X-D3/D4/V 不是简单别名；不得切换后仍报告原方法。本阶段直接阻断原始指定配方。[PySCF 色散接口源码](https://pyscf.org/_modules/pyscf/scf/dispersion.html)

原始方法需在有可验证实现的后端完成能量/梯度与文献基准复核，或由用户明确批准不同方法并建立新的研究配置。本机未安装 PySCF；没有为了“过测试”改用已安装 Psi4 或把 xTB 当第一性原理证据。

### 6.2 已写的适配边界

PBE0 的 RKS/UKS 单点/几何优化适配仅供**显式授权的方法替代集成测试**，并非本次执行结果。其适用范围限指定整数 charge/spin=2S、有限坐标、真空、列明基组。wavefunction.chk、stdout/stderr、数值参数与输入哈希必须保留；exit 0 不等于收敛。

生产必须检查 SCF 最终收敛、轨道稳定性、UKS 自旋污染、网格/基组/环境敏感性、梯度与几何阈值、零点及热修正适用性。实现返回 stability_verified=False；不会伪造稳定性通过。

CASSCF 使用独立 RHF/ROHF 参考，指定活性轨道与电子数，再作 SC-NEVPT2。它不是给任意 DFT 总能量叠加一个修正。[PySCF 多参考微扰文档](https://pyscf.org/user/mrpt.html)

活性空间选择、轨道跟踪、root flipping、自然占据数、state averaging、内禀 intruder/state 问题均需独立审查；当前适配未执行、未保证覆盖所有多参考体系。真实非绝热动力学还需要态间耦合与交叉缝证据，单个 CASSCF 能量不够。

### 6.3 自愈限制

5 次数值策略依次调整 max_cycle、damping、level shift、DIIS/EDIIS；只对明确 unconverged 的返回结果进行重试。缺包、未支持方法、超时、语法错误进入 BLOCKED，不自动改写任意 Python。

改变泛函、基组、溶剂、总电荷、自旋、分子组成或几何身份都要新任务/新指纹。找到较低能的另一个电子态不等于原态计算成功。收敛失败也不能推导“不可合成”或直接拉黑分子。

sandbox_executor 只运行冻结的可信量化 worker，shell=False，受限环境变量、单线程、超时与独立日志；**不是执行不可信 LLM 代码的安全沙盒**。生产需要容器/Job Objects/cgroups、只读挂载、禁网、磁盘限额和子进程树清理。

## 7. 模块 03：CI-NEB、TS 与 IRC

当前 ASE 接口要求相同原子顺序、独立 image calculators，用 IDPP 初值和 climbing-image NEB/FIRE；输出 image 势能 eV 与轨迹。没有真实 QM calculator 连接，也未执行本体系路径。

合格 TS 还必须同时满足：端点均为对应方法下的极小值；鞍点收敛；去除整体运动后只有一个显著虚频；该振动模式对应目标反应；双向质量加权 IRC 通向预期且重新优化后的两端；电子态在路径上可跟踪。程序 -20 cm⁻¹ 是契约中的显著虚频筛查阈值，不是普适物理常数。

audit_ts 验证元数据与哈希格式，不独立认证哈希所指文件的科学正确性；必须人工/独立解析器复核。IRC 驱动明确抛出阻断异常，不把几何下降轨迹重命名为 IRC。基线非共价且无指定反应，不造“共价弹头加成 TS”。

## 8. 模块 04：QM/MM 与自由能管线

### 8.1 QM/MM 体系契约

每个原子恰好属于 QM 或 MM；边界键/链接原子与电荷重分配须审查。QM/MM 应明确 electrostatic embedding、长程电势、cutoff/PME 处理、link correction 与双计数去除，而不只是把两个引擎的能量相加。先用小体系有限差分检验总力，再测试能量漂移和边界敏感性。

当前只有分区与边界输入校验，**没有**可运行的耦合能量/力服务器；不发布“真实 QM/MM 已完成”的标题。

### 8.2 热力学循环

设 T=BRD4，E=VHL，L=配体。选择 E 与 L 的结合腿，对比 T 不在场与已经形成 TL 的条件。不要混淆 VHL 的二元 Kd 与 BRD4 的二元 Kd。

FEP 变换 A→B 必须在相同原子映射定义下同时考察自由配体溶剂腿、相关二元环境和三元环境。计算相对 ΔΔG 不能无参考地获得绝对 Kd；盐浓度、质子化、标准态和 restraint 定义须一致。

示意循环：
ΔG_bind = ΔG_complex − ΔG_solvent + ΔG_restraint + ΔG_standard-state + ΔG_finite-size。
符号约定需绑定具体 coupled/decoupled 端点；以上是本项目分析函数的输入约定，不是任意 decoupling 软件可直接照抄的符号。[原始炼金自由能最佳实践](https://arxiv.org/abs/2008.03067)

### 8.3 采样计划与实际缺口

输出 17 个 λ 窗口（0 至 1）和 3 个独立种子，温度示例 298.15 K；这是计划不是运行。实际窗口数量须由重叠和端点困难调整。软核势处理 LJ/静电端点；电荷改变需有限尺寸修正；复杂约束需适当标准态/约束自由能。

REST2 是副本 Hamiltonian/有效溶质温度维度，不等同 λ 维度。Metadynamics 需指定可靠 CV、bias history、参数、重加权及独立重复；未经重加权的 biased histogram 不等于自由能。

OpenMM 核心的可导入状态不意味着已安装 PLUMED 插件。[OpenMM 官方扩展说明](https://docs.openmm.org/latest/userguide/application/05_add_on_packages.html) 当前环境缺 OpenMM–PLUMED/PyMBAR；采样器接口主动阻断，未产生假轨迹。

### 8.4 已实现的分析函数

TI：ΔG=∫₀¹〈∂H/∂λ〉dλ，以非等间距梯形权重 w 计算；均值协方差 Σ 传播为 Var(ΔG)=wᵀΣw。所需输入是时间相关性处理后的**均值协方差**，不是原始帧的方差。当前标准误未包括积分离散误差。

前向 FEP：ΔG=−RT log〈exp(−ΔU/RT)〉，用 logsumexp 防止溢出；返回权重 ESS，不声称它是时间去相关后的有效样本数，也不把一次前向估计视为收敛。MBAR/BAR、overlap matrix 和 decorrelation 驱动仍未接入。

拟议采样接受门（未达到）：每条相邻 λ 边连通、ESS/相关时间合格、至少 3 个独立重复、前后半程稳定、正反路径一致、热力学循环闭合误差在预先定义不确定度内。数值阈值须先用基线校准，不为候选结果事后调整。

## 9. 模块 05–06：合成对抗与三元协同性

路线模块检查 reaction-step DAG、重复/缺失节点、环路及最长线性步数 >12 的项目策略。总反应节点数不等于最长线性步数。UNKNOWN 不等于不可合成；图通过也不代表环化选择性、收率和纯化可行。没有伪装 ASKCOS/AiZynthFinder 的本地随机评分。

α=Kd_binary/Kd_conditional，ΔΔG=G_conditional−G_binary=−RT lnα。α>1 对应更有利的条件结合。协方差保留共享参照项，不以独立误差假设重复增加不确定度。

文献的实验值须严格保留测量类型和统计来源；不得从某个测量的 α 直接复制到新候选。本项目实际输出 alpha=null。不可逆共价捕获需 k_inact/K_I 等动力学描述，不能强制套 Kd；近接几何比例不能改名“有效泛素化碰撞概率”。

Pareto 候选集需要真实且相容的自由能、uncertainty、合成证据和特性指标。当前 `pareto_top_leads.json` 为空集合，这是正确拒绝缺失数据，不是优化出“零活性”。

## 10. 资源计划、停止条件与验收

本机 CPU/内存有限：pilot 单线程、小样例；不安装新包，不提交后台付费任务。当前不估算没有基准测量的全项目天数，也不把 1 亿 Tokens 作为已购买的算力。

扩展路线：

1. 完成基线结构准备、参数检查与小体系能量/力回归。
2. 在用户授权的服务器验证准确方法实现；做基组/网格/活性空间敏感性。
3. 小批候选验证工作流，测实际 CPUh/GPUh、RSS、I/O 和失败率。
4. 基线自由能及不确定度通过后，才分批扩展至 50–100 个候选。
5. 与独立实验组核对合成/结合/三元/细胞结果，锁定外部盲测候选。

停止条件：方法不支持、电子态不确定、关键参数缺失、结构准备失败、SCF 5 次上限、采样不重叠、闭环偏差无法解释、模型假说被反证、Token/资源配额不足。失败是保留的研究结果，不以更宽阈值或 fabricated replacement “修复”。

验收分工程验收与科学验收两张表。现阶段只对已执行的工程/低成本构象部分验收；未执行适配器不能声称无缺陷或生产级。原始代码、测试、实际产物和边界声明共同构成交付，不靠报告长度验收。

## 11. 图与产物规范

- swarm_dialogue_dynamics.png：实际本地角色事件数；标题标注非 LLM 对话，预算与实际调用分开。
- macrocycle_strain_map.png：同分子内 MMFF94s 相对构象能量与单个环二面角；不是绝对环张力。
- free_energy_landscape_3d.png：没有真实/可信重加权样本时不生成；缺失状态在 artifact_status.json 中明示。
- error_healing_log.jsonl：无实际 QM 尝试则为空，不虚构自愈案例。
- quantum_states、fep_trajectories：无科学计算成果，仅保留说明。
- 实际版本、测试与图检查覆盖：见 [local_validation.md](local_validation.md)。

所有数字按文献、解析测试、实际本地计算、待验证假说分类。不可用单元测试通过来证明药效、TS 连通或 FEP 收敛。
