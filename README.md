# AI4Chem · 分阶段项目集合

**每个阶段一个文件夹，代码、技术报告、结果与图表集中存放。**

[阶段与证据总表](docs/PHASE_INDEX.md) · [复现指南](docs/REPRODUCE.md) · [证据说明](docs/EVIDENCE.md) · [共享工具](shared/) · [旧路径迁移表](docs/layout_manifest.json)

新增：[Phase 28–29 公开数据与本机计算交付](docs/PHASE28_29_STATUS.md)。

## 目录结构

```text
projects/
  phase01/ … phase29/
    README.md       阶段入口与证据链接
    run.py          统一运行入口
    code/           本阶段专属源码（共同实现见 shared/）
    reports/        中英文技术报告、阶段更正
    results/        已保存结果（若有）
    figures/        图表（若有）
    tests/          阶段专属测试（若有）
  legacy_solubility/ 历史辅助项目
shared/             跨阶段代码、测试、审计证据、文献来源
tools/              目录检查与兼容运行工具
docs/               导航、复现、证据与历史文档
```

## 打开一个项目

| 阶段 | 项目文件夹 | 代码与技术报告 |
|---|---|---|
| 01 | [复杂骨架与构象](projects/phase01/) | [代码](projects/phase01/code/) · [报告](projects/phase01/reports/) |
| 02 | [扭转扫描与分子动力学](projects/phase02/) | [代码](projects/phase02/code/) · [报告](projects/phase02/reports/) |
| 03 | [KRAS 复合物与结合分析](projects/phase03/) | [代码](projects/phase03/code/) · [报告](projects/phase03/reports/) |
| 04 | [骨架编辑与反应路径](projects/phase04/) | [代码](projects/phase04/code/) · [报告](projects/phase04/reports/) |
| 05 | [反应网络与逆向设计](projects/phase05/) | [代码](projects/phase05/code/) · [报告](projects/phase05/reports/) |
| 06 | [显式溶剂与增强采样](projects/phase06/) | [代码](projects/phase06/code/) · [报告](projects/phase06/reports/) |
| 07 | [强关联与模型适用边界](projects/phase07/) | [代码](projects/phase07/code/) · [报告](projects/phase07/reports/) |
| 08 | [光化学与非绝热动力学](projects/phase08/) | [代码](projects/phase08/code/) · [报告](projects/phase08/reports/) |
| 09 | [自驱动实验室模拟](projects/phase09/) | [代码](projects/phase09/code/) · [报告](projects/phase09/reports/) |
| 10 | [连续流反应器数字孪生](projects/phase10/) | [代码](projects/phase10/code/) · [报告](projects/phase10/reports/) |
| 11 | [神经波函数与 VMC](projects/phase11/) | [代码](projects/phase11/code/) · [报告](projects/phase11/reports/) |
| 12 | [动力学定律发现](projects/phase12/) | [代码](projects/phase12/code/) · [报告](projects/phase12/reports/) |
| 13 | [金属酶与 PCET](projects/phase13/) | [代码](projects/phase13/code/) · [报告](projects/phase13/reports/) |
| 14 | [活性物质与凝聚体](projects/phase14/) | [代码](projects/phase14/code/) · [报告](projects/phase14/reports/) |
| 15 | [自由基对与变构模型](projects/phase15/) | [代码](projects/phase15/code/) · [报告](projects/phase15/reports/) |
| 16 | [核孔复合物与输运](projects/phase16/) | [代码](projects/phase16/code/) · [报告](projects/phase16/reports/) |
| 17 | [相对论量子化学](projects/phase17/) | [代码](projects/phase17/code/) · [报告](projects/phase17/reports/) |
| 18 | [全细胞代谢与热力学](projects/phase18/) | [代码](projects/phase18/code/) · [报告](projects/phase18/reports/) |
| 19 | [酶设计与几何闭环](projects/phase19/) | [代码](projects/phase19/code/) · [报告](projects/phase19/reports/) |
| 20 | [CISS 输运与自旋电子学](projects/phase20/) | [代码](projects/phase20/code/) · [报告](projects/phase20/reports/) |
| 21 | [腔 QED 与极化激元](projects/phase21/) | [代码](projects/phase21/code/) · [报告](projects/phase21/reports/) |
| 22 | [分子自旋量子比特](projects/phase22/) | [代码](projects/phase22/code/) · [报告](projects/phase22/reports/) |
| 23 | [分子忆阻器与储备池计算](projects/phase23/) | [代码](projects/phase23/code/) · [报告](projects/phase23/reports/) |
| 24 | [AI 高通量催化模拟](projects/phase24/) | [代码](projects/phase24/code/) · [报告](projects/phase24/reports/) |
| 25 | [不对称催化](projects/phase25/) | [代码](projects/phase25/code/) · [报告](projects/phase25/reports/) |
| 26 | [恒电位铜界面](projects/phase26/) | [代码](projects/phase26/code/) · [报告](projects/phase26/reports/) |
| 27 | [Ni 光动力学](projects/phase27/) | [代码](projects/phase27/code/) · [报告](projects/phase27/reports/) |
| 28 | [工业氯碱有机电极与寿命成本](projects/phase28/) | [任务书](projects/phase28/TASK_PROMPT.md) · [代码](projects/phase28/code/) · [报告](projects/phase28/reports/) |
| 29 | [医药吡啶位点控制与可制造性](projects/phase29/) | [任务书](projects/phase29/TASK_PROMPT.md) · [代码](projects/phase29/code/) · [报告](projects/phase29/reports/) |

## 运行示例

在已具备对应科学依赖的 Python 环境中，从仓库根目录执行：

```sh
# 查看全部阶段，无需导入科学依赖
python tools/run_phase.py --list

# 在新目录准备原有运行布局，不启动计算
python projects/phase23/run.py --workspace work/runs/phase23 --prepare-only

# 在同一工作目录执行 Phase 23 的软件自检
python projects/phase23/run.py --workspace work/runs/phase23 --reuse -- --self-test

# 新阶段直接在项目布局运行，将新结果写入独立目录
python tools/run_phase.py 28 --workspace work/runs/phase28
python tools/run_phase.py 29 --workspace work/runs/phase29
```

旧脚本依赖的相对路径与跨阶段导入，通过独立工作目录兼容；源码不需要复制回仓库根目录。运行输出留在工作目录，不覆盖已提交结果。完整说明见[复现指南](docs/REPRODUCE.md)。

## 研究状态

本项目包含真实计算、有效模型、合成数据和待验证方案，各阶段成熟度不同。
Phase 24 是模拟高通量原型；Phase 25–27 尚未完成科学验收。历史“运行中”只表示报告当时的状态。
Phase 28–29 已开展公开文献数据分析与本机计算；文献测量、分子计算、假设模型和待验证方案分别标记。暂无自有实验原始数据与 HPC，不能据此声称工业电极寿命或可切换位点选择性已验证。
目录整理不等于重新计算或科学结论通过；请结合[证据总表](docs/EVIDENCE.md)阅读结果。

代码许可：[MIT](LICENSE)。第三方坐标与数据保留各自来源和许可。
