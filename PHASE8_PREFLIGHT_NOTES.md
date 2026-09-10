# Phase 8 queued-run preflight — 2026-09-10

Scope: `run_phase8_photochemical_dynamics.py`, `test_phase8_preflight.py`,
and this note only. Main owns integration, scheduling, and runtime. No campaign
snapshot, result, environment, or other project was changed; no QC/trajectory,
commit, or push was performed. No callable parallel subagent facility was exposed
to this assigned task, so no additional user-owned tasks were created.

Read chem-ai4s SKILL.md and its quantum/molecules references,
PHASE6_8_RERUN_NOTES.md and PHASE6_10_AUDIT_FIXES.md first. The memory registry
search had no relevant Phase8 entry. Current scheduling instruction supersedes
the older notes: Phase8 remains queued behind 16/12/13/15; this task did not launch it.

## Concrete corrections

- CLI honors inherited OMP_NUM_THREADS, defaults to two threads when unset,
  and rejects nonpositive threads/missing worker arguments before dispatch.
  Child interpreter remains the installed phase7 Python, with its own DLL PATH;
  both PYTHONPATH and PYTHONHOME are removed. Thread budgets propagate to workers.
- The existing pre-spawn 2.5 GiB available-memory gate and QC_MEM="2 GB" remain
  unchanged. Timeout stdout/stderr are now appended to the worker diagnostic log.
- Installed xTB `--help` exposes `--uhf`, not `--mult`: use `--uhf 0` for the
  existing neutral closed-shell molecules. Nonzero exits preserve stderr and fail.
  Structure/Hessian failure no longer silently substitutes guessed/MMFF geometry.
- Installed Psi4 `driver/procrouting/proc.py` reads TDSCF_TRIPLETS; the old
  tdscf_singlet/singlet attempts silently ignored unsupported options. Use
  NONE for singlets and ONLY for triplets. Singlet NTOs are captured before triplet
  results overwrite wavefunction variables. Cube orbitals use one-based indices
  for the intended occupied/virtual columns; the coefficient matrix is restored.
  Missing requested NTOs/cube pairs and cubeprop exceptions are failures.
- Validate complete finite 8A state/grid and 8B scan/cut/vector records, including
  both original MECI gap/gradient gates, before accepting results. Incomplete cached
  8B records cannot feed dynamics. Lower-basis or reduced-state fallback results
  are no longer accepted as successful full-default protocol results. Existing
  computational retry tiers were not expanded or used by this audit.
- Checkpoints bind exact input geometry/scan files, driver source, basis, and smoke
  mode via SHA256. Mismatched/legacy checkpoints fail explicitly instead of silently
  restarting/mixing protocols. JSON replacement is atomic within the output folder.
  Partial torsion/cut/stretch scans retain their checkpoint phase and raise before
  advancing. Resumed initial random h-vector directions consume the same random
  sequence as uninterrupted work even when previously evaluated IDs are skipped.
- Full-run master starts clean; partial reruns discard old full-run certification.
  A failed rerun cannot retain the old module result. Structure/8B reruns invalidate
  dependent module records. Successfully retried stage errors are cleared.
- --fig_only uses the recorded failure path, just like --stage fig. Missing, empty,
  or unrefreshed expected PNGs fail explicitly, including renderer silent skips.

Scientific constants, scan grids, active space, tolerances, MECI gates,
300 default trajectories, 0.5 fs timestep, 500 fs duration, and sampling seed
were preserved. No substitute calculation was launched or declared equivalent.

## Exact verification

Runtime:
`C:/Users/HUIWEI/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe`.
All test processes used OMP_NUM_THREADS=MKL_NUM_THREADS=OPENBLAS_NUM_THREADS=1,
PYTHONDONTWRITEBYTECODE=1 and `subprocess.run(..., timeout=25)`.

Final unittest arguments (passed to that Python after `-B -m unittest -v`):

```text
test_phase8_preflight
test_phase6_8_rerun.RerunTests.test_memory_gate_precedes_psi4_process_creation
test_phase6_8_rerun.RerunTests.test_inverse_bfgs_satisfies_secant_and_positive_curvature
test_phase6_8_rerun.RerunTests.test_meci_requires_finite_gap_and_gradient_with_original_gates
test_phase6_8_rerun.RerunTests.test_phase8_workers_keep_thread_budget_and_logs
test_phase6_8_rerun.RerunTests.test_phase8_worker_failure_preserves_diagnostics
test_phase6_8_rerun.RerunTests.test_phase8_dispatch_passes_threads_to_every_worker
test_phase6_10_audit
```

Result: **27 passed in 1.589 seconds**. New suite has 14 stdlib AST/mock tests;
existing tests include the LVC finite-difference force regression. Tests never
import Psi4 or execute molecular engines, FSSH ensembles, or figure rendering.
Initial new-suite failure was an AST test-loader omission of tuple-assigned
constants, corrected before the successful runs.

Additional read-only probe: installed phase2ff `Library/bin/xtb.exe --help`,
one thread, 15-second subprocess timeout, exit 0, confirming --chrg and --uhf.
Scoped `git diff --check` passed.

## Still unvalidated / integration requirements

- Use a fresh campaign output/snapshot with these edits; do not transplant them
  into already-running snapshots. No historical structure/result cache provenance
  was certified. Legacy checkpoints without the new identity intentionally reject.
- No real Psi4 SCF/TDA/CASSCF, xTB optimization/Hessian, DLL import, or numerical
  convergence validation occurred. Installed source API inspection and mocks do
  not establish successful calculations or scientific validity.
- No actual cube or PNG was rendered/visually inspected. NTO/cube consistency and
  figure output require the full run. Plot libraries/runtime availability remain
  untested. The freshness check assumes ordinary local filesystem mtimes.
- No process RSS/peak-private measurement or memory-budget optimization was done;
  2 GB is engine allocation, not total-process RAM. Main must honor existing gates.
- Native abort/restart numerical reproducibility (including h-vector refinement),
  convergence of the MECI, checkpoint job-count accuracy across multiple restarts,
  downstream lifetime/yield/FSSH physics and model adequacy remain unvalidated.
- QA assertions expose incomplete/changed-protocol results; they do not make
  fallback bases equivalent or guarantee the planned QC can converge. Main must
  inspect result gates and diagnostics before reporting full scientific completion.
