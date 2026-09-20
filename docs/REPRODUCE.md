# 复现入口 / Reproduction guide

[返回首页](../README.md) · [项目目录](../projects/) · [证据说明](EVIDENCE.md)

代码已按阶段实际归档。旧脚本内部的相对路径、包内导入和跨阶段依赖保持原样，由统一入口在独立工作目录重建原有运行布局。请不要直接在 `code/` 内启动原脚本。

## 准备与运行

复用已有、满足对应阶段依赖的 Python 环境。在仓库根目录运行：

```sh
# 查看阶段
python tools/run_phase.py --list

# 仅准备目录，不启动计算
python projects/phase23/run.py --workspace work/runs/phase23 --prepare-only

# 运行软件自检；-- 后的参数原样传给阶段程序
python projects/phase23/run.py --workspace work/runs/phase23 --reuse -- --self-test
```

也可一次完成准备与自检：

```sh
python tools/run_phase.py 23 --workspace work/runs/phase23-new -- --self-test
```

`--python PATH_TO_PYTHON` 可以指定已安装的科学计算解释器。统一入口本身只依赖 Python 标准库，不安装软件、不下载数据、不启动未指定的任务。

Windows 的 Conda 环境须先激活，让 `Library/bin` 中的数值库 DLL 可见；只指定 `python.exe` 路径不等于完整环境激活。本次检查使用已有环境，未安装或升级依赖。

## 依赖

各项目的 `requirements*.txt` 已随项目移动；例如 Phase 23 的依赖位于 [projects/phase23/requirements_phase23.txt](../projects/phase23/requirements_phase23.txt)。仅在需要的独立环境中安装：

```sh
python -m pip install -r projects/phase23/requirements_phase23.txt
```

没有经验证的单环境覆盖全部 27 阶段。Phase 2–19 的 OpenMM、OpenFF、xTB、Psi4 等要求以各阶段技术报告为准；机器专用解释器路径仍需自行适配。

## Phase 25–27 与辅助工具

三个阶段共享原来的 `phase25_27` 包。专属源文件已分别归档，运行入口会还原完整包，保持相对导入：

```sh
# 数学/结构单元测试；不提交量子化学计算
python projects/phase25/run.py --workspace work/runs/phase25-tests --module phase25_27.test_analysis

# 只准备 Phase 26 的工作目录，不运行会改写输入的初始化程序
python projects/phase26/run.py --workspace work/runs/phase26 --prepare-only
```

三项任务尚无完整生产入口，因此 25–27 必须显式指定模块，或只准备目录。可用模块与原命令见[共用指南](../shared/phase25_27/README.md)。运行 `build_inputs` 等初始化工具前仍须核对它们的覆盖行为。

使用原辅助脚本时，可显式指定迁移表中的旧脚本名：

```sh
python tools/run_phase.py 4 --workspace work/runs/phase04-help --script diagnose_phase4_neb.py -- --help
```

也可以进入已准备的工作目录，执行原报告中的命令；那里的 `run_phase*.py`、`results_phase*/` 与 `phase25_27/` 均已还原。历史命令中的机器路径和依赖条件并未自动修复。

## 工作目录与证据

- 每次首次准备都会复制完整映射文件集，包含跨阶段依赖和已提交结果，占用约一个仓库数据副本的空间；没有使用可写硬链接。
- 必须指定新目录，或对已有兼容目录明确传入 `--reuse`。来源变动或工作目录源码被改写时，复用会被拒绝。
- 默认输出位于工作目录；用户显式传入的绝对输出路径仍由原程序处理。已提交的项目数据不会被准备流程覆盖。
- 改善导航而调整过链接的 Markdown，在来源未被进一步编辑时还原原始字节，便于旧哈希审计；未来对代码的编辑会随下一次准备进入新工作目录。
- [迁移清单](layout_manifest.json)记录旧路径、新路径与原始/整理后 SHA-256。科学代码和非 Markdown 结果文件保持原字节；文档正文中的旧路径是历史运行布局。
- 本次是目录与复现兼容性验证，不是全阶段生产重算。科学验收状态见[证据说明](EVIDENCE.md)。

## 后续维护

修改已有文件后，新工作目录会使用修改后的代码；不要手动修改已准备工作目录后再把它当作未变更来源复用。新增或改名源文件时，也应更新迁移清单的 `old` / `new` 映射，以便运行入口包含它。历史哈希保留为迁移时的基线，新的科学运行另记自己的来源哈希。
