# Phase 19 · 酶设计与几何闭环

[全部项目](../../README.md) · [阶段总表](../../docs/PHASE_INDEX.md) · [复现指南](../../docs/REPRODUCE.md)

本目录集中保存本阶段的专属代码、报告与结果。跨阶段共用文件在 `shared/`，下方提供直接入口。

## 技术报告

- [ACTIVE_INFERENCE_ENZYME_REPORT_EN.md](reports/ACTIVE_INFERENCE_ENZYME_REPORT_EN.md)
- [ACTIVE_INFERENCE_ENZYME_REPORT_ZH.md](reports/ACTIVE_INFERENCE_ENZYME_REPORT_ZH.md)
- [PHASE19_CLASH_NOTES.md](reports/PHASE19_CLASH_NOTES.md)
- [PHASE19_LOOP_CLOSURE_NOTES.md](reports/PHASE19_LOOP_CLOSURE_NOTES.md)
- [PHASE19_ORIENTATION_NOTES.md](reports/PHASE19_ORIENTATION_NOTES.md)

## 代码与数据

- [代码](code/)
- [结果](results/)
- [图表](figures/)
- [测试](tests/)
- [跨阶段收尾清单](../../shared/audit/reports/PHASE1_19_CLOSEOUT.md)
- [公开重算证据](../../shared/audit/evidence/phase1_19_rerun_20260910/README.md)

## 运行

在仓库根目录执行。使用已有、满足本阶段依赖的 Python；工作目录必须是新目录。

```sh
python projects/phase19/run.py --workspace work/runs/phase19 --prepare-only
```

上面只复制运行所需的完整布局，不启动计算。代码文件保留原始内容；它们的旧相对路径和包导入在工作目录中保持兼容。

准备完成后，按原技术报告在该工作目录内执行命令；也可使用统一入口透传参数。已有目录仅在明确指定 `--reuse` 且源码匹配时复用。

**证据边界：** 文件归档与软件检查不是新的科学验收。详见[证据说明](../../docs/EVIDENCE.md)及本阶段报告。
