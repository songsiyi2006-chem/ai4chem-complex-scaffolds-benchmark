# 新工具与数据库适用性核查

快照：2026-09-10T16:25:07.663924+00:00。本轮确认接口可用性与模型力的一致性，不改变 Phase25–27 的科学验收状态。

|资源|本轮实际核查|可用范围与限制|
|---|---|---|
|MACE-MH-1，omat_pbe|4 原子 Cu 能量 -14.969970265 eV；位移结构解析力与有限差分力误差 1.202e-09 eV/Å|可作为铜几何预处理候选；不是恒电位电子结构引擎。对称体相零力不能证明晶格常数已平衡。|
|PubChem|CID 87116 与现有映射结构去除映射后的二甲酯型 Hantzsch 酯一致|结构身份核对；不验证反应条件或 ee。名称查询 404 后，结构查询成功。|
|Materials Project OPTIMADE|mp-30 返回成功；本地识别 Fm-3m，对称容差 0.01 Å|铜参考结构；JARVIS 本次返回 HTTP500。联网可查询不等于完整数据库已下载。|
|Basis Set Exchange|取得 Ni、Cl、C、N、H 的 def2-TZVP 及引用|可用于 Phase27 基组收敛输入；它本身不执行多参考、SOC 或 NAC 计算。|
|ORD|本地 750/750 条 protobuf 记录成功解析|用于反应来源与数据管线；MCP 直接序列化 BLOB 失败，可用本地只读解析替代。没有用这些记录充当盲测标签。|
|MACE-POLAR、OrbMol-v2|模型文件大小与 SHA256 核验通过|可尝试分子预筛，须先做目标体系 DFT 验证；未建立激发态分支能力。MCP 分子接口限制 100 原子，完整 117 原子体系应从本地 Python 调用。|
|QM9、MoleculeNet、Matbench|目录合计 15 项，含 PDB 与 ORD 示例|可用于描述符或软件基线；不提供本题指定 ee、恒电位势垒或 Ni 轨迹标签。|

ORD 本地样本的反应类型：1.3.1 [N-arylation with Ar-X] Bromo Buchwald-Hartwig amination: 262; 1.3.4 [N-arylation with Ar-X] Iodo Buchwald-Hartwig amination: 76; 0.0 [Unassigned] Unrecognized: 92; 1.3.2 [N-arylation with Ar-X] Chloro Buchwald-Hartwig amination: 290; 9.7.39 [Other functional group interconversion] Chloro to amino: 1; 1.3.7 [N-arylation with Ar-X] Chloro N-arylation: 9; 1.3.9 [N-arylation with Ar-X] Iodo N-arylation: 1; 1.3.6 [N-arylation with Ar-X] Bromo N-arylation: 16; 1.3.3 [N-arylation with Ar-X] Iodo Buchwald-Hartwig amination: 3。产物测量类型计数：{'YIELD': 750}；选择性子类型：{}。这是样本范围核查，不代表全 ORD 检索。

4 个本地模型文件的哈希均一致。保留 MACE 权重的 ASL 使用条件；OrbMol 许可按安装清单记录。当前配置目录仍标记 PySCF/GPU4PySCF/周期 DFT 引擎未安装，UMA 仍需访问授权且未安装。工具增多确实改善结构准备和数据获取，但没有新增可调用的 OpenAI GPU 配额。

接入顺序：继续现有 Psi4 频率队列；先以合适 DFT 参考验证分子模型的力和相对能量，再扩大构象预筛；铜模型先做表面和吸附物独立检查，再用于几何预处理。恒电位显式溶剂采样与 Ni 多参考动力学仍需各自的实际计算，不能由上述模型检查替代。

[核查汇总](summary.json) · [ORD 范围](ord_audit.json) · [PubChem 原始响应](pubchem_structure.json) · [铜结构原始响应](copper_mp30.json) · [基组与引用](basis_Ni.json) · [力检查输入输出](force_check.json)
