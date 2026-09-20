# 资源、执行范围与后续接口 / Resources and execution boundary

用户本次条件为：没有新增实验数据和 HPC，只完成公开数据与本地计算。已经复用既有 chem-ai4s 环境，无安装或升级，无高等级界面电子结构任务。

| 层级 | 本次状态 | 内存/线程/时间预算 | 恢复及停止规则 |
|---|---|---|---|
| local 公共数据提取 | EXECUTED | ≤0.5 GiB，2 线程，预计 <2 min | 来源下载完成后校验 SHA256；单位冲突不作绝对归一化 |
| local 模型与控制 | EXECUTED | ≤0.5 GiB，2 线程，预计 <2 min | 固定 seed/config；新输出目录；可从提取 CSV 独立重跑 |
| local 测试与绘图 | EXECUTED | ≤0.5 GiB，2 线程，预计 <1 min | 任何守恒、因果性、批次泄漏检查失败即阻止结论升级 |
| HPC 溶剂化恒电位界面 | NOT_RUN / INPUTS_NOT_READY | 探索性估计 32–128 核、64–256 GiB/任务，需服务器实测预算 | 先验证结构、锚定/覆盖度、表面/电解液和电位标尺，再做收敛 pilot |
| experiment 独立电极与产物测量 | NOT_RUN | 本次不申请、不执行 | 仅提供需观测字段；外部有资质实验室独立定义实施方案 |

## 为什么不做伪界面计算

Nature Synthesis 2026 的公开 CIF 证明特定晶体结构可访问，但并不提供完整工作电极锚定、溶剂、离子、膜、覆盖度和恒电位边界。一个分子晶体或真空气相 xTB 能量不能替代界面活化自由能。本次 `configs/hpc.json` 是缺项契约，不是可以宣称 `INPUTS_READY` 的 DFT 输入。

若将来取得适用输入，先定义恒电荷与恒电位的映射、参考 NHE/RHE、pH、离子强度、电子计量和溶剂；用小体系核查泛函/溶剂/构象/覆盖度敏感性，差异超过候选优势就停止排序。关键路径需连接性和采样证据，不能把单虚频等同完整机理。

## 本地再运行

```powershell
$env:PATH = 'C:\Users\HUIWEI\miniconda3\envs\phase2ff\Library\bin;' + $env:PATH
$env:OMP_NUM_THREADS='2'
$env:MKL_NUM_THREADS='2'
& 'C:\Users\HUIWEI\.codex\tools\chem-ai4s\venv\Scripts\python.exe' projects/phase28/run.py --self-test
& 'C:\Users\HUIWEI\.codex\tools\chem-ai4s\venv\Scripts\python.exe' projects/phase28/run.py --out work/runs/phase28_new
```

跨平台：使用已有 Python 3.12、numpy、scipy、matplotlib，直接运行同一命令；仅从原 XLSX 重新提取时额外需要 openpyxl。默认不联网，输出目录已存在时拒绝覆盖。数值结果采用固定种子；运行时/平台信息不属于化学证据。
