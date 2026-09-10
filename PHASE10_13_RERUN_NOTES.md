# Phase 10–13 fresh rerun preparation — 2026-09-10

## Scope and current status

Only the four assigned phase runners, `test_phase10_13_rerun.py`, and this note
were edited. No stage/commit/push, hardware protocol execution, or shutdown.
Read `PHASE6_10_AUDIT_FIXES.md` and `PHASE11_13_AUDIT_FIXES.md` before changes.

Main initially held heavy jobs with available RAM below 1 GB, then explicitly
allocated one serial slot at approximately 3.3 GB free. The requested command
`rerun_phase1_19_campaign.py --phases 10,11,12,13 --profile molecular` is now
started. Phase 10 snapshot `phase10/20260910T135356` (PID 11692) completed FULL
at 14:20:07; independent retained-process-handle watcher returned exit code 0.
Queue parent PID 44852 was stopped before handoff at low memory, without stopping
the child. Main was asked to reconcile its stale runner status. Phase 11–13 now
launch one at a time, each after free RAM reaches at least 1.5 GiB.

## Corrections

- All four runners now respect `OMP_NUM_THREADS` (default two) for their
  explicit Torch, OpenMM CPU, or Psi4 worker thread settings.
- Phase 11 and 13 Psi4 child processes prepend their own environment root and
  `Library/bin`, remove inherited `PYTHONPATH`, and check nonzero exit codes.
  Phase 11 no longer disables the SCF maximum-iteration failure gate.
- Phase 11 raises on nonfinite local energies or parameter gradient norms,
  instead of zeroing nonfinite training weights. Nuclear cusp diagnostic now
  evaluates on the nuclear z axis, matching the actual H2 geometry.
- During Phase 10 execution, a Phase 11 partial-reference bug was corrected
  before its snapshot: failed systems now receive explicitly listed fallbacks
  even if another system's Psi4 calculations succeed. The result no longer
  falsely marks partial success complete or crashes on missing system keys.
  One synthetic partial-success regression passed in 0.058 seconds.
- Phase 12 loss logging detaches the scalar tensor. C++ execution checks now
  record measured/not-completed status in JSON, including exceptions, rather
  than omitting an unavailable check from results.
- Phase 13 figure-only rendering no longer discovers/initializes Psi4.
  Full/default remains production/all and still executes every module.
- Phase 13 bath charge order is O,H,H per water (`tile`, not grouped `repeat`).
  Potential evaluation uses a consistent absolute coordinate frame. The new
  `bath_observables` helper is tested against direct sums and translation.
- Phase 13 Coulomb field conversion is 1.43996455 GV/m per e/nm^2; the old
  factor 10.36427 had incompatible units. The applied field remains 10.36
  kJ/(mol nm e); its comment now correctly identifies approximately 0.1074 GV/m.
  The gap conversion remains numerically 138.935 kJ nm/(mol e^2).
- Phase 13 Langevin random seed is explicitly seven.

## Cheap verification and prerequisites

Molecular Python imports verified: Torch 2.10.0, NumPy 2.4.6, SciPy 1.18.0,
OpenMM 8.6, SymPy 1.14.0; CUDA unavailable. Psi4 1.11 imports in `phase7` with
its DLL path. `qbscf/python.exe` exists but `import pyscf` raises
`ModuleNotFoundError`. Phase 13 retains its explicit Psi4 backend substitution.
`g++` is available at `E:/mingw64/bin/g++.exe` for Phase 12's C++ check.

Focused command (PowerShell, from this source checkout):

```powershell
$env:PATH='C:/Users/HUIWEI/miniconda3/envs/phase2ff/Library/bin;C:/Users/HUIWEI/miniconda3/envs/phase2ff;'+$env:PATH
$env:OMP_NUM_THREADS='2'
$env:MKL_NUM_THREADS='2'
$env:OPENBLAS_NUM_THREADS='2'
& C:/Users/HUIWEI/miniconda3/envs/phase2ff/python.exe -X faulthandler -m unittest -v test_phase10_13_rerun test_phase11_13_audit
```

Initial result: 12 tests passed in 6.216 seconds. Tests execute selected actual
definitions via AST to avoid output-directory/import side effects. They cover
two SAC optimizer updates, VMC antisymmetry, AD vs finite-difference Laplacian,
finite parameter gradients, two training steps for each scalar mode, scalar
audit semantics, PCET overlaps/rate grid refinement and electrostatic invariance.
Final verification after subsequent edits added `test_phase6_10_audit` to that
command: **19 tests passed in 5.922 seconds**. Scoped `git diff --check` passed.

## Defaults to run only after allocation

The existing `rerun_phase1_19_campaign.py` copies top-level Python source into
`work/phase1-19-rerun-20260910/phaseXX/TIMESTAMP`, hashes sources, records the
command/PID/status, and captures `run.log`. It copies no historical results,
checkpoints, cache, or figures. Each runner constructs its results under its
snapshot directory; defaults require no pre-existing local data artifacts.
Figure-only and checkpoint reuse options must not be passed for these reruns.

Run each allocated phase separately (replace `10` with the allocated phase):

```powershell
& C:/Users/HUIWEI/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe rerun_phase1_19_campaign.py --phases 10 --profile molecular
```

The campaign runner supplies the molecular DLL path and two-thread environment.
Do not queue all four phases without main's allocation.

- Phase 10: default FULL plant/control simulation, 110 SAC training episodes and validation,
  PAT audit, technoeconomic grid and figures. No `--selftest` or `--skip_train`.
- Phase 11: five systems; H2 equilibrium and He 1200 epochs, three other H2
  geometries 550 epochs, 2048 walkers each; 200 warm-up updates and frozen-network
  production sampling remain intact. Psi4 references attempted with 600-second
  per-job timeouts; existing literature fallback remains explicitly reported.
- Phase 12: default full telemetry/discovery, 4000 steps per scalar network, Lyapunov spectra,
  C++ export/execution comparison, and figures. No `--quick`.
- Phase 13: default `--tier production --stage all`: 13C runs 1250 equilibration
  plus 12500 production MD steps before 13A electronic jobs, then 13B PCET and
  13D spectroscopy. Production QC timeout limits are 2700 s per ladder job,
  1500 s imidazole, and 3600 s per scan job; multiple fallback tiers can run.

## Limits requiring full output review

Passing these tests does not establish training convergence, chemical accuracy,
global scalar conservation/stability, full-run RAM fit, or converged PCET QC.
Review Phase 13 per-job `converged`, `acceptance`, selected method tiers and
fallback fields, including spectroscopy contact-density fallback. Its existing
acceptance doctrine and reorganization-energy clipping were not relaxed.
Process exit zero alone is not scientific acceptance. Retain logs of failures
and make another fresh source snapshot after any further fix.

## Phase 10 completed full attempt and timestep audit

110 SAC episodes, all 24 held-out controller/scenario combinations, both
3600-second demonstrations, PAT audit, technoeconomic grid and figures completed.
Held-out breach40 rates: open loop 1/8, SAC and shielded SAC 0/8 each.
Demo peak rises: open loop 71.5532 K, shielded 21.2962 K; arrest time 10 s.
PAT mean/max relative error 0.00528155/0.0468663.

At main's request, audited the slow hot open-loop segment. It advanced through
the full horizon and is not an endless timestep loop. Both recorded histories
have 721 finite frames (temperature, coolant and all concentrations).
Read-only evaluation of snapshot `_rates` across recorded states, conservatively
setting maximum fouling, gives max effective rates 42.08965/s (open) and
1.36093/s (controlled). Required chemical timesteps are at least 7.12764 ms and
220.438 ms, respectively: zero recorded frames require below the old 1 ms floor.
Substeps were not saved, so this cannot certify every intermediate microstep.
Final times are 3600.01041463 and 3600.00208529 s: small endpoint overshoot exists.

Source now removes the timestep floor, preserving the CFL/reaction/end-interval
minimum, and raises on nonfinite input/state, nonadvancing timestep or exhausted
200000-substep budget. One real-source synthetic timestep test passed in 0.075 s.
These guard changes were made after this snapshot: its outputs are not claimed
to validate the new guards. No observed critical underresolution currently
justifies prioritizing another complete 25-minute training run over Phase11–13;
this limitation and tradeoff were reported to main. No old policy was reused.

## Phase 11 running full attempt

Fresh snapshot `phase11/20260910T142310`, child PID 11196, runner session 64397.
Launch memory gate passed with 2147560 KiB free. All five Psi4 reference sets
completed in 167 seconds with no fallback. Full numerical self-tests passed:
antisymmetry deviation zero; AD Laplacian relative error 3.52e-7.

First equilibrium-H2 epoch appeared at 14:31:46 after the default 200 warm-up
updates and 400 burn-in sweeps. CPU time advanced and process threads were
running; this was a long unlogged stage, not a demonstrated deadlock. Epoch50
appeared at 14:33:08, giving 82/49 = 1.673 s per additional training epoch.
At epoch50 energy was -1.174043 +/-0.002095 Eh and variance .008993 Eh^2;
these are training diagnostics, not final production statistics.

Future source snapshots now log warm-up first/every25 updates with timing,
burn-in every100 sweeps, and mean epoch timing. Running snapshot and all default
walker/update counts are unchanged. Phase12 and13 remain pending independent
memory checks after each preceding phase completes.

At the 18:08 resume, original child11196/parent16468 were still running. The
equilibrium-H2 CSV contains all1200 epochs; production summary E=-1.173979
+/-0.000321 Eh, absolute computed-FCI deviation0.463mEh, final variance.00315.
The chemical-accuracy criterion passes for this system. H2_R2.5 has started;
no Phase12/13 snapshot or process exists. Reported12900s spans the interruption
and is wall time, not an appropriate uninterrupted compute-speed estimate.

## Durable handoff at18:42,2026-09-10

Original Phase11 runner16468/child11196 remain alive (retained session64397).
Snapshotphase11/20260910T142310; no duplicate run has been launched.
H2_R2.5 finished550epochs and production: E=-1.093561+/-0.000203Eh,
computed-FCI error0.382mEh, variance.00244, chemical-accuracy gate passed.
H2_R4.0 now at epoch300of550, E=-1.014054+/-0.001082Eh training diagnostic.
H2_R6 and He still pending after R4; wholePhase11 incomplete, no final exitcode.
Phase12/13 remain unlaunched, no fresh completion claims. Before launching each,
check preceding child real exit/status/output and freeRAM>=1.5GiB, molecular
environment PATH prefix Library/bin+root andOMP/MKL/OPENBLAS threads2. Current
user has not authorized new background automation; none was created. Main is
preparing handoff, so reconcile its current allocation before a new heavy child.
