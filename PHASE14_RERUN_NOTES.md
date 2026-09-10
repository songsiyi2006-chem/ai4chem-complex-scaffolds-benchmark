# Phase14 fresh rerun follow-up

## 2026-09-10: diagnostics timing and adaptive-step averages

The ongoing full attempt `20260910T135701` must be preserved until it exits.
It has completed both 160x160 main trajectories and is progressing through the
28-point 112x112 scan, before the three FRAP/SAXS calculations. Its saved source
is not identical to the subsequently corrected root script.

Two diagnostic defects were identified in the root source and corrected:

- A frame's timestamp was advanced, but its mean concentrations still came from
  the preceding state. Means now describe the accepted endpoint; the diffusion
  dissipation proxy is explicitly labeled as a left-endpoint evaluation.
- NESS statistics averaged the last 20% of adaptive-step *records* equally.
  They now average the final 20% of elapsed physical time, weighting accepted
  intervals and clipping the first interval at the averaging-window boundary.

Seven deterministic tests in `test_phase14_rerun.py` pass, including an unequal
timestep example giving 3.6 for the full-window mean and 2.0 for the last half.
These are unit tests, not full grid/time convergence evidence.

The running attempt cannot acquire these source edits retroactively. Its output
will therefore retain the old averaging limitation. It saves only decimated E1
frames and aggregate E2/E3 statistics; these are insufficient to reconstruct an
exact corrected time-weighted average after exit. Do not silently relabel those
aggregates as corrected. A subsequent full fresh run of corrected code remains
required for current-source end-to-end validation, after existing work finishes
and resource availability permits it. Do not shorten the grids or 1000 s windows.

Earlier corrections to snapshot timing, endpoint integration and finite-state
checks are covered by the same tests. The running snapshot predates some of them;
an exit code of zero would not alone certify the latest root source.

## Checkpoint recovery (2026-09-10)

The checkpoint store commits each completed E1 state, E2 parameter point,
and E3 equilibration/FRAP unit using an atomic directory rename. Source and
scientific configuration hashes must match on explicit `--resume`; changed or
incomplete stores fail closed. Arrays are loaded without pickle. Time histories,
fields and snapshots are retained, including corrected time-weighted means.
An interrupted subexperiment must restart: this is not an integrator-step
checkpoint. Resume to a fresh output directory with the byte-identical script
and an absolute `--checkpoint-dir`. Existing final outputs are not overwritten.

Thirteen Phase14 tests passed in the molecular environment. This does not
replace the pending default-grid, full-duration scientific rerun. The old
complete run retains its statistical limitations; the later 20:29 attempt was
interrupted and had no checkpoints from this new implementation.
