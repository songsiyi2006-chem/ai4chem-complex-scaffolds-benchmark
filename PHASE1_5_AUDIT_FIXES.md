# Phase 1-5 code audit and corrections / 代码审计与修正

## Status / 状态

The corrections below have been implemented and regression-tested. Full historical
conformer ensembles, production MD, NEB/xTB Hessians and catalyst searches have NOT
been rerun. Historical numerical files and figures remain for provenance, not as
post-fix validated results. This note supersedes stronger claims in the archived
Phase 1-5 reports and README result tables.

以下修复已经实施并完成针对性回归测试；未重跑全部分子构象、生产级动力学、
NEB/xTB Hessian或催化剂搜索。原始数据和图像保留用于追溯，不代表修复后的
已验证结果。本说明优先于历史报告和README中相关的较强结论。

## Implemented changes / 已实施的修改

1. **Phase 1:** only finite energies with optimization flag zero participate in
   conformer selection, ensemble spread and SDF export. All-rejected ensembles
   become a 3D failure rather than `ok=True`. Accepted counts are recorded.
2. **Phase 2:** optimization must converge and achieve the circular torsion target
   within 5.1 degrees (a +/-5 degree flat-bottom constraint is used). Energies are
   evaluated with a fresh unrestrained force field. The actual minimum's coordinate
   snapshot, energy and achieved angle are exported, not the final scan frame.
3. **Phase 3:** an independent OpenMM interaction-group kernel evaluates only
   unscreened receptor-ligand Coulomb and Lennard-Jones terms. The incorrect
   ligand-charge-zeroing difference is removed. Cross-fragment exceptions and
   alchemical parameter offsets fail closed. This validates the explicit pair sum,
   NOT GB, cutoff/PME treatment, decomposition at other cutoffs, or binding entropy.
   Pre-audit MM-GBSA checkpoints are invalidated.
4. **Phase 4:** failed NEB/TS checks stop downstream stages and set failure status.
   Numerical Hessian, independent gradient stationarity, a significant imaginary
   mode (< -20 cm^-1), a reactant-minimum screen and a relative bond-displacement
   mode check are required. Failure records contain no thermochemistry/rate.
   Passing candidates remain explicitly **not IRC-verified**; the thermal estimate
   is vibrational-only, not a complete validated reaction free energy. Electronic
   differences now use the same calculator as the Hessian screening. Fabricated
   fallback rates and nonexistent TS-refinement labels are removed. Old validation
   caches and unvalidated figure-only output are rejected.
5. **Phase 5:** the fallback entropy term converts cal to kcal once, removing the
   extra gas-constant factor. Its approximate translation/rotation entropy remains
   a limitation. Geometry caches use a new versioned directory and old result
   hydration is rejected. A constrained complex-energy difference alone can no
   longer set `claim_proven=True`: gradient, saddle and IRC verification of the
   compared systems are absent. Module D is opt-in with `--exploratory-kinetics`;
   caps preserve both negative barrier shifts and the stereo sign. Reports identify
   the network, assigned rates/floors and kinetic corrections as exploratory.

## Verification / 验证

```bash
python -m pip install numpy rdkit openmm
python -m unittest -v test_phase1_5_audit
```

Tests include Python syntax, convergence filtering through the actual Phase 1
function, a real RDKit butane torsion scan/SDF-energy consistency check, a real
OpenMM cross-interaction calculation, reaction-screen rejection cases, figure
blocking, the entropy-unit regression, signed scenario caps, Module-D opt-in and
old-cache rejection. Heavy legacy scripts are not imported wholesale: isolated
AST-extracted function bodies avoid filesystem, external-engine and shutdown
side effects. No GPU, xTB or production MD is executed by these tests.

测试使用真实RDKit和OpenMM，但不等于全部科学结果已经复现。
修复前的约束能污染、TS门槛缺失和未经验证的催化势垒结论，应先撤回，
在明确运行条件下重新计算后再恢复。

## Rerun boundaries / 重算边界

- Phase 1: rerun the benchmark into a new output directory.
- Phase 2: rerun torsion scans; do not reuse the old `*_scan_min.sdf` as minima.
- Phase 3: rerun MM-GBSA analysis on the intended trajectory; use independent
  sampling/entropy work before interpreting it as a binding free energy.
- Phase 4: use `--force_rerun` in a dedicated output directory. A failed screen
  requires better optimization, not relaxation of acceptance thresholds. Follow
  the imaginary mode/IRC in both directions before claiming reaction connectivity.
- Phase 5: rerun A/B/C with the corrected caches. D requires deliberate scenario
  opt-in; it is not a replacement for solvation calibration or verified TS searches.

No historical result file was silently overwritten with invented replacement data.
