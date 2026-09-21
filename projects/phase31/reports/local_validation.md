# Phase31 本机验证与交付记录

日期：2026-09-21。基础仓库版本：`db9bec08a5f92362b7b5cf393559896b749cbc69`。
本记录描述本地验证结果；发布版本及时间以仓库 Git 提交历史为准。提交或推送不改变下述科学证据边界。

## 实际执行

| 检查 | 结果 | 适用范围 |
|---|---|---|
| Phase31 单元/集成测试 | 49 项通过 | 包括两个独立 worker 进程，不是多节点集群 |
| 仓库 native phase 入口回归 | 4 项通过 | 不是全仓库科学回归 |
| Phase31 统一入口 local-pilot | 正常退出 | 7 角色，3 成功完成任务、4 科学阻断、0 软件失败 |
| 同工作区 resume | 正常退出，无角色重算 | 来源指纹一致时的幂等恢复 |
| 原始结构 | 6SIS.pdb + LFE_ideal.sdf | 公开参考；未修补/参数化/溶剂化 |
| MMFF94s 试验 | 3 个样例，18/18 构象最小化收敛 | 12、13、14 元手性内酰胺，不是靶向配体 |
| 产物哈希 | 15 个文件全部匹配 | 清单见下方 |
| 相对能量独立复算 | 18 个值符合 E−E_min | PowerShell 重新从原始能量计算，不依赖绘图函数 |
| 图形检查 | 2/2 实际 PNG 目视检查 | 修复事件热图标签重叠；黑白编码，不依赖颜色 |

测试命令：

```powershell
python projects/phase31/run.py --self-test
python -m unittest discover -s tools -p test_native_phases.py -v
python tools/run_phase.py 31 --workspace projects/phase31/results/local_pilot -- --mode local-pilot
python projects/phase31/run.py --out projects/phase31/results/local_pilot --resume
```

复现时请替换为新工作区，不能覆盖已保存例子。当前解释器为已安装的 chem-ai4s 环境；依赖版本见 [preflight.json](../results/local_pilot/preflight.json)。本机 RDKit 2026.3.5、NumPy 2.4.6、SciPy 1.18.0、Matplotlib 3.11.1、OpenMM 8.6.0、ASE 3.29.0；PySCF、OpenMM–PLUMED、PyMBAR 未发现。

## 原始结果与复核入口

- [冻结的源码、提示词和参考输入哈希](../results/local_pilot/manifest.json)
- [角色结果快照](../results/local_pilot/results/status.json)
- [导出的哈希链接事件](../results/local_pilot/reports/events.jsonl)
- [原始 SDF 与构象能量](../results/local_pilot/attempts/)
- [15 文件哈希清单](../results/local_pilot/reports/artifact_hashes.json)
- [各缺失产物的明示状态](../results/local_pilot/results/artifact_status.json)

SQLite 运行库被 .gitignore 排除，保留于本地供恢复；共享的是 JSON/JSONL 快照和原始构象。其他机器应新建运行目录，不能仅靠快照直接恢复数据库。

发布时保持 PDB/SDF 原始字节及其哈希不变；阶段内 .gitattributes 仅对这两种保留格式化空格的数据文件关闭空白风格告警，不豁免 Python/Markdown/JSON 的检查。

三个样例的收敛构象相对能量范围分别为 0–9.0702、0–6.7894、0–10.6840 kcal/mol。每个分子各自选零点，**不得跨分子比较总能量或由此宣称环张力排序**。6 构象不是充分采样，也没有统计自由能。

## 本次确实修复的问题

1. 冻结预算不一致时，SQLite 初始化异常未关闭连接；在 Windows 测试中导致文件被占用。已加入异常路径关闭并通过重测。
2. Windows 路径分隔符进入可移植结果引用；改用 POSIX 相对路径并新建运行，保留旧本地输出用于追踪。
3. 事件热图横轴标签重叠；改为自适应宽高与旋转标签，重新生成并目视确认。
4. 恢复参数检查补充识别 `--mode=...` / `--token-cap=...`，避免静默忽略显式覆盖请求。

图检查未进行专门的色觉模拟；本次图只用灰阶、文字和空心圆。没有三维自由能图，因为没有合格采样数据。

## 不在验收范围

没有实际调用外部 LLM、提交 HPC、执行 DFT/CASSCF、NEB/IRC、QM/MM、FEP/REST2 或 PLUMED。PySCF/ASE 适配器是未用真实目标体系验证的代码边界，不是经过生产负载验证的驱动。没有新分子胶、合成路线或药效结论。

提供方 Token 实际回执为 0，上限 100,000,000；此账本不统计编写项目的 Codex 对话。没有消费承诺，也未为了凑预算造假。
