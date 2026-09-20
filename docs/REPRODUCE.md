# 复现入口 / Reproduction guide

[返回首页](../README.md) · [阶段导航](PHASE_INDEX.md) · [证据说明](EVIDENCE.md)

从一个阶段开始，先核对其报告、依赖和输入。仓库包含不同计算引擎与原型，没有经验证的“单环境运行全部 27 阶段”流程。复用已有环境；不要为了阅读结果安装整套软件。

## 1. 只查看已保存结果

通过[阶段导航](PHASE_INDEX.md)打开报告、数据和图表即可。核对审计记录，区分历史输出与修复后重算；无需执行任何计算。

## 2. 一个独立模块的软件自检示例

Phase 23 使用 NumPy、SciPy 和 Matplotlib，提供独立的 `--self-test` 入口。
在仓库根目录、已满足对应依赖的环境中：

```sh
python run_phase23_in_materio_neuromorphic_computing.py --self-test
```

若需要新环境，可在独立检出副本中创建虚拟环境，激活后再安装[该阶段依赖](../requirements_phase23.txt)：

```sh
python -m venv .venv-phase23
# Windows PowerShell:
.venv-phase23/Scripts/Activate.ps1
# macOS/Linux instead: source .venv-phase23/bin/activate
python -m pip install -r requirements_phase23.txt
python run_phase23_in_materio_neuromorphic_computing.py --self-test
```

自检是软件检查，不是完整基准运行，也不证明实验性能。本次文档整理未重新执行这一自检；命令入口来自[脚本](../run_phase23_in_materio_neuromorphic_computing.py)。

## 3. 按任务选择运行指南

| 任务 | 依赖 / 操作入口 | 注意事项 |
|---|---|---|
| Phase 1 构象基准 | [代码](../molecule_benchmark.py) · [报告](../BENCHMARK_REPORT_EN.md) | 根目录 requirements 不覆盖所有扩展模块 |
| Phase 2–8 分子模拟 / 量子计算 | [阶段导航](PHASE_INDEX.md)中的报告和脚本 | 按阶段核对 OpenMM、OpenFF、xTB、Psi4 等环境 |
| Phase 1–19 审计重算 | [公开证据](../audit/phase1_19_rerun_20260910/README.md) · [运行器](../rerun_phase1_19_campaign.py) | 运行器含原机器路径；先适配，避免直接启动整套长任务 |
| Phase 20–23 独立模型 | 各自 `requirements_phaseNN.txt` 与脚本 | 依赖和输出目录按阶段隔离 |
| Phase 24 高通量模拟 | [验证说明](../PHASE24_VALIDATION.md) · [依赖](../requirements_phase24.txt) | 用新的输出目录；保留已有 measurements |
| Phase 25–27 | [模块详细指南](../phase25_27/README.md) · [依赖](../requirements_phase25_27.txt) | 初始化可能重置未测标签；运行前保护现有 campaign |

## 4. 输出与验收

- 部分脚本使用固定相对输出路径；完整运行应在独立副本或明确的新输出目录进行。
- 保存命令、源码提交、随机种子、依赖版本、日志和结果哈希。
- 先检查完整性、收敛、采样与对照，再解释数值。失败或负结果也属于有效记录。
- 文献数据与恢复坐标遵守原始许可；仓库代码的 MIT 许可不覆盖所有第三方来源。

完整命令和历史机器配置可在[旧首页](../README_HISTORY.md)查到，但应先与现有脚本、环境和后续更正核对。
