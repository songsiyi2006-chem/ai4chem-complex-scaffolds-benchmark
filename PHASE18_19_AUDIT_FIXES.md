# Phase 18–19 audit and recalculation — 2026-09-10

## Status / 本次验证范围

Phase 18: the complete calculation was rerun (curation, LP/MILP, glucose-uptake sweep, dynamics and figures). After identifying additional dynamic stoichiometry errors, dynamics were integrated again and figures regenerated. Phase 19: two **isolated-fragment QM scans**, nine points each, were recalculated using actual GFN2-xTB energies and gradients. **The complete five-generation enzyme evolution, flow training, protein MD and enzyme-environment scan have not been rerun successfully.** Historical Phase 19 design results and figures 1–3 are not post-fix validated outputs.

Phase 18 已完成全流程重算，并在补修动力学计量后重新积分与绘图。Phase 19 完成水/乙酸根两个孤立反应片段的各 9 点 QM 重算；不是完整酶设计流程重跑。早期酶构建尝试中断于昂贵的几何/折叠分布构建阶段，未得到修复后完整流程结果。不得将本次孤立片段结果用作历史五代酶进化或催化速率提升的证据。

## Corrections

### Phase 18

- LP relaxation now relaxes only binary activation bounds; it no longer replaces formation-energy uncertainty bounds with [0, 1].
- MILP retry preserves the requested glucose uptake and checks thermodynamic signs against the effective retry epsilon. Boundary/biomass reactions are excluded from the enzymatic/transport sign test.
- Corrected three dynamic stoichiometries: one fructose bisphosphate consumed per aldolase event; ADP produced by ATP-consuming glycogen synthesis; no net CoA loss for the lumped alpha-ketoglutarate-to-succinate step.
- Stage-specific reruns preserve other result sections instead of silently discarding them.
- Driving-force figure now uses solved, uncertainty-corrected energies for the TCA panel and only thermodynamically gated reactions for its sign panel. Dynamic figures use the computed exhaustion time; time markers are no longer drawn on a concentration-axis phase portrait.

### Phase 19

- Corrected two orthogonal-projection formulas in the theozyme coordinate frame.
- Acetate QM fragment includes its three methyl hydrogens. The isolated substrate/water and substrate/acetate systems contain 19 and 23 atoms, respectively, with total charges 0 and -1.
- Corrected the xTB closed-shell argument to `--uhf 0`; subprocess decoding is explicit. Exit status, finite energy/gradient/charge arrays and optimizer convergence are checked.
- Native fixed-atom optimization failed the distance guard in this environment. Production scans now use ASE BFGS with both endpoint atoms fixed exactly, xTB energy/gradients, a 0.05 eV/angstrom free-atom force tolerance and at most 500 steps per scan point. Constraint failures are rejected, not accepted as scan data. This fixes Cartesian endpoints, which is more restrictive than a pure distance constraint.
- Replaced all-zero fragment charges with geometry-specific xTB charges for the post-hoc Coulomb calculation. This is **not self-consistent electrostatic embedding**: MM charges do not polarize the QM calculation or drive geometry optimization. The existing 3-angstrom charge-deletion shell is an additional approximation; no full enzyme-environment result was revalidated here.
- Peak C–H/N–O distances now come from the highest-energy sampled geometry, not the last scan point. Result records retain each geometry, charge vector and absolute/relative energy, and count successful gradient evaluations.
- Removed experimental free-energy anchoring and unsupported exponential rate-enhancement estimates from figure generation and the acceptance summary. Legacy fields `barrier_kcal` and `ts_*` mean sampled potential-energy peak and its geometry, **not a validated transition state**.
- `--stage figures` now loads saved results without accidentally launching evolution; corrected summary sequence lookup. Added `--stage qm-audit` and a separate figure 4 so new isolated scans cannot be confused with historical enzyme results.

## Recomputed results

| Quantity | Result | Interpretation |
|---|---:|---|
| Phase 18 LP relaxed growth bound | 0.3246662971 /h | Corrected LP bounds |
| Phase 18 tFBA growth | 0.3234185314 /h | MILP status 0 |
| Gated active sign violations | 0 | Within numerical tolerance |
| Active-network loop certificate | pass, loop flux 0 | Model-specific stoichiometric test |
| Glucose uptake sweep | 16 levels, 1–45 | Completed |
| Dynamic glucose exhaustion | 1944.849998 s | Computed, not imposed figure text |
| Relative adenylate drift | 2.4353e-15 | Previously approximately 0.003586 |
| Relative NAD-pool drift | 3.7470e-15 | Numerical conservation check |
| Water fragment scan | 9 points; peak 36.57 kcal/mol | Peak at final 1.00 angstrom endpoint |
| Acetate fragment scan | 9 points; peak 23.33 kcal/mol | Peak at 1.19375 angstrom |
| xTB evaluations | water 311; acetate 427 | Includes successful pre-relaxation |

All 18 saved scan points completed optimization. Independent checks of saved coordinates found maximum O–H constraint error 4.44e-16 angstrom and maximum atomic-charge sum error 3.00e-8 electrons; all stored energies were finite. A preceding five-point acetate scan gave 10.81 kcal/mol, demonstrating substantial grid/path sensitivity; nine points do **not** establish grid convergence. Each curve is referenced to its own first point rather than a validated reactant minimum. The water endpoint maximum is not a located saddle point. These data cannot establish activation free energies, relative catalytic rates or experimental enzyme performance. TS optimization, Hessian/IRC checks, solvent/protein sampling and free-energy calculations remain necessary.

Phase 18 remains a constructed coarse-grained model, not an experimentally calibrated whole-cell predictor. The very broad formation-energy uncertainty bounds (up to ±2000 kJ/mol), disconnected/orphan metabolites and reduced kinetic model limit biological interpretation. Numerical feasibility and pool conservation do not establish physical accuracy or conservation of every chemical element.

## Reproduction

```sh
python run_phase18_wholecell_metabolic_thermodynamics.py --stage all --milp-time 120
python run_phase18_wholecell_metabolic_thermodynamics.py --stage dynamics
python run_phase18_wholecell_metabolic_thermodynamics.py --stage figures
python run_phase19_active_inference_denovo_enzyme.py --stage qm-audit
python -m unittest -q test_phase1_5_audit test_phase6_10_audit test_phase11_13_audit test_phase14_17_audit test_phase18_19_audit
```

The first three commands used the primary NumPy/SciPy/Matplotlib Python environment with NetworkX available. The QM command used the existing `phase2ff` environment (ASE, RDKit, Torch, OpenMM, NumPy/SciPy and an accessible xTB executable), with `OMP_NUM_THREADS=2` and `MKL_NUM_THREADS=2`. An initial Phase 18 attempt in the molecular environment exited natively during least-squares curation; the primary environment completed it. Numerical library versions can affect optimization paths.

Tests: 41 total passed, including six new targeted tests for source compilation, LP bounds, thermodynamic gating, complete fragments/peak geometry/charge updates, Coulomb gradient finite differences and dynamic conserved pools. Mock-based tests do not certify the full protein pipeline; real QM evidence is separately recorded below.

## Output provenance

- [Phase 18 results](results_phase18/phase18_results.json), [curation](results_phase18/phase18_curation.json), [dynamics arrays](results_phase18/phase18_dynamics.npz).
- [Phase 18 driving forces](figures_phase18/fig2_thermodynamic_driving_forces.png), [dynamic trajectories](figures_phase18/fig3_dynamic_metabolic_rewiring.png).
- [New Phase 19 scan geometries, charges and energies](results_phase19/phase19_qm_audit.json), [new isolated-fragment plot](figures_phase19/fig4_qm_audit_recomputed.png).
- Existing `phase19_results.json`, champion PDBs, flow checkpoints and Phase 19 figures 1–3 remain historical, not newly recomputed results. The accompanying older scholarly reports are superseded by this note wherever they conflict.
