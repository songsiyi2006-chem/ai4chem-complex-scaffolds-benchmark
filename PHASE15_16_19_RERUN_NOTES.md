# Phase 15 / 16 / 19 preparation and fresh-rerun status

2026-09-10. Coordinator initially allocated sequential **15 -> 19 -> 16**,
then explicitly allocated a second slot for concurrent Phase19.
Phase15 started at 13:54:07 in `../phase1-19-rerun-20260910/phase15/20260910T135407`.
Phase19 started at 14:12:04 in `../phase1-19-rerun-20260910/phase19/20260910T141204`
after a measured 3,023,508 KiB free RAM (above its 2 GiB launch gate).
Latest scheduling: HOLD Phase16; prioritize the fresh Phase19 retry, but
require BOTH >=2GiB free RAM and coordinator/Epicurus confirmation that the
Phase6 implicit Context transition (around14:39) is complete/resumed. Until
that confirmation, do not launch even if RAM briefly exceeds2GiB. Keep15
running. No new phase starts below500 MB free.
Preparation checks below are focused regressions, dependency imports and
bounded microbenchmarks, NOT full production results. Completed full-run
results will be recorded below; no full enzyme evolution is yet claimed.

## Owned changes

- `run_phase15_quantum_biology_spin_allostery.py`: declare Bohr units for
  already-converted PySCF coordinates; make `--stage hfcc` actually execute;
  isolate alternate interpreter PATH/PYTHONPATH; respect requested OpenMM
  CPU thread count instead of explicitly taking up to ten threads.
- `run_phase16_megamachine_cryoem_transport.py`: reject unordered/nonfinite
  window coordinates, fewer than two force samples, nonfinite forces,
  nonpositive/nonfinite friction, invalid resistance or hydrodynamic radius.
  These checks do not certify adequate sampling or PMF convergence.
- `run_phase19_active_inference_denovo_enzyme.py`: cache the 3375 Glu chi
  triples and rigidly rotate them for the same 36-roll exhaustive grid;
  evaluate collision clearance after translating stems into the pocket.
  Remove 30 identical deterministic N-cap searches (the old loop variable
  was unused). Correct perpendicular projections in rod construction,
  indole alignment and dihedrals. Attach CB at CA instead of carbonyl C;
  align the acetate's actual O-to-C vector away from the substrate with
  an orthonormal rigid transform. Read reference CA coordinates with the
  original PDB topology, not hydrogen-renumbered atom indices. Compute RMSF
  as sqrt(mean(squared displacement)), and honor the CPU thread limit.
  Log fold proposal/acceptance progress every 25 proposals.
- `test_phase15_16_19_rerun.py`: isolated real-source regressions for these
  corrections, without running production entrypoints.
- This note. No stage/commit/push, unrelated source edits or hardware actions.

The speed optimizations preserve the candidate grid and remove duplicate
work; geometry corrections necessarily change selected structures. Existing
convergence and numerical acceptance thresholds have not been relaxed.

## Dependency probes and focused verification

Molecular environment, with its own Library/bin and root first on PATH:
NumPy 2.4.6, SciPy 1.18.0, RDKit, Torch 2.10.0, ASE and OpenMM import
successfully. xTB found at
`C:/Users/HUIWEI/miniconda3/envs/phase2ff/Library/bin/xtb.EXE`.

`qbscf/python.exe` exists but `import numpy` fails with
`ModuleNotFoundError: No module named 'numpy'`. The coordinator independently
confirmed absent NumPy/PySCF. No unsupported Windows native PySCF build was
attempted. Phase15 default fallback must be reported as literature tensors,
NOT a successful first-principles hyperfine calculation.

Run focused checks from this directory in PowerShell:

```powershell
$env:PYTHONPATH='../phase24-deps;../phasefix-test-deps'
$env:OMP_NUM_THREADS='2'
$env:MKL_NUM_THREADS='2'
$env:OPENBLAS_NUM_THREADS='2'
& C:/Users/HUIWEI/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe -m unittest -q test_phase15_16_19_rerun test_phase14_17_audit test_phase18_19_audit
```

Final suite: **23 tests passed in 3.490 seconds**, including the RMSF
regression. `git diff --check` passed for all three owned production files.

## Phase19 runtime estimate: measured preparation, conditional total

Fresh process, molecular Python, CPU, two threads, 30-second hard process
guard; no checkpoint or historical result loaded. Final geometry fixes were
present for these measurements:

| Small measurement | Elapsed |
| --- | ---: |
| First Theozyme including rotamer-grid construction | 3.172 s |
| Repeated Theozyme using in-process grid | 0.265 s |
| One accepted fold, one proposal | 0.328 s |
| Three batch-1 training steps, including optimizer startup | 3.266 s |
| Separate warmed 10-step batch-1 synthetic training probe | 1.250 s |

These probes are not full training or evolution. The full default still uses
70 training structures, 1400 steps with batch size 6, five generations with
six candidates each, and nine points per QM scan.

- If the one measured proposal is representative, 70 accepted structures
  take approximately 23 seconds; the 2800-proposal cap corresponds to about
  15 minutes at that per-proposal cost. Acceptance rate and cost may vary.
- Linear batch scaling of 0.125 s/step gives about **17.5 minutes** for
  full flow training. Allow roughly **15-35 minutes** for planning; batch-6
  memory and throughput were NOT measured while RAM/slots are constrained.
- Cached theozyme construction for 30 candidates is only about eight seconds;
  this excludes ODE generation, backbone realization, annealing and rotamer
  packing. The latter still searches 10,800 Glu and 2,700 Trp chi candidates
  per design, plus continuous polishing and other sidechains.
- MD can verify up to ten candidates. Its nominal 22-minute shared budget
  is NOT a total runtime bound: minimization and 15,000 equilibration steps
  per candidate precede the budget-tested production loop, and per-candidate
  budgets have a 90-second floor. No full protein MD throughput was measured.
- Up to five enzyme scans plus one water reference require nine relaxed
  points each, up to 500 BFGS steps per point, plus pre-relaxation and line
  searches. No fresh enzyme QM throughput was measured; convergence failure
  remains possible and must be reported, not bypassed.

Practical allocation: **reserve at least a two-hour initial exclusive slot,
with possible extension to several hours**. This is a scheduling estimate,
not a demonstrated full-run ETA or upper bound. A credible narrower total
needs the first real candidate construction/MD/QM timings after allocation.

## Full-attempt commands, only after coordinator allocation

Run separately so the coordinator can allocate each phase independently:

```powershell
& C:/Users/HUIWEI/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe rerun_phase1_19_campaign.py --phases 19
& C:/Users/HUIWEI/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe rerun_phase1_19_campaign.py --phases 15
& C:/Users/HUIWEI/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe rerun_phase1_19_campaign.py --phases 16
```

The shared runner selects molecular/primary environments and sets PATH and
two-thread limits; it snapshots sources into a new
`../phase1-19-rerun-20260910/phaseXX/TIMESTAMP` with no historical cache.
Review `rerun_status.json`, `run.log` and fresh phase outputs. Process exit
zero is not scientific acceptance. Phase16 default budget is 110 minutes;
Phase15 OpenMM sampling budget is 30 minutes, excluding other stages.

## Outstanding scientific limitations

Before the full Phase19 start, the coordinator explicitly required fixing
the invalid 25 A RMSF threshold: default is now **0.8 A**. Candidate acceptance
also requires successful MD with >=2 frames and positive duration, finite
CA RMSD <=4 A, actual all-anchor constellation RMSD <=0.30 A, and a successful
finite potential scan <=12.5 kcal/mol. The previous constellation metric
included only two backbone anchors, omitting the sidechain targets; this is
corrected. Every candidate records missing/failed gates and error details.
Failed MD placeholders are excluded from model observations. Per-generation
progress is persisted. No accepted design leaves champion null and sets
`acceptance_status=no_accepted_design`; no champion PDB is exported for it.
Additional champion MD is explicitly recorded as not performed. Tests after
these acceptance changes: **24 passed in 4.272 seconds**.

Static gates remain loose; sidechain grafting and spline loops can introduce
strain. No gate was relaxed by this patch; model-gate acceptance is not proof
of biochemical activity.
QM remains an isolated substrate/acetate optimization with post-hoc protein
charges, not self-consistent embedding or a validated transition state/free
energy. The named final champion MD budget is not actually invoked by the
existing final-champion block. These limitations preclude certified enzyme
performance even if an eventual default execution completes.

## Actual full reruns

### Phase15 controlled stage-boundary continuation

Coordinator authorized stopping ONLY the tracked Phase15 PID after its
complete fresh spin/MFE/ODMR checkpoint is saved, to avoid MD with known
incorrect input geometry. A100ms lightweight watcher validates that JSON has
all five full-dynamics curves and MFE/ODMR, verifies PID9344's command line,
then stops it and records the SHA256. At14:46 the result directory was still
empty; the last Bzero curve was still in memory. No valid spin calculation
is intentionally discarded or rerun.

Corrected source additionally fixes Phase15 CA-CB attachment, per-residue
L-stereochemistry, ASN/ASP/GLN branch terminals (previously chained to the
oxygen rather than the carbonyl carbon), and C-terminal O/OXT overlap
(torsions previously differed by360degrees). Umbrella bias centroid weights
are now uniformly1, matching the arithmetic-mean CV saved in samples.
34tests passed6.770s after these changes and stage-assembly tests.

Existing `--stage allostery` runs corrected geometry/sampling in a separate
fresh snapshot without recomputing spin. Added `--stage assemble --audit-input
FRESH_SPIN_JSON --allostery-input FRESH_ALLOSTERY_JSON` to combine those fresh
stages, including separate input hashes/source hashes and explicit
`single_default_run=false`; it refuses to overwrite an existing target.
Original population arrays and inputs remain unchanged; yields are separately
reintegrated with Simpson and diagnostics. This is a staged corrected result,
not a claim that the original default attempt completed unmodified.

### Phase16 static memory estimate (no heavy probe)

Default solvent lattice85x85x89=643,025 sites before exclusion;14,784solute
beads (4,480scaffold+10,304FG). GridSolvent retains344bytes/solvent bead:
five Nx8 arrays and three N scalars, about211MiB before exclusion. Each
solvent Nx3 double array is14.7MiB; multiple positions/velocities/forces/noise
and retained lattice copies add roughly100-180MiB. Corner-index and gradient
temporaries add39.2MiB each; several can coexist. Python neighbor-query lists,
solute cell/pair tables, imports and allocator retention add further overhead.
Estimated16A process peak roughly0.6-0.9GiB, not a hard upper bound.

Phase16B uses228^3 voxels: eachfloat64 master map90.4MiB. The round-trip
comparison expression can hold roughlyfour full-size arrays simultaneously
(~362MiB), plus binned maps, model maps, imports and prior allocator retention.
Estimated16B process peak roughly0.5-0.8GiB; data-dependent pairing and memory
allocator behavior can increase peaks. Stage16C is mostly14,784solute beads
plus cargo and pairs, without the large solvent buffers.

Conclusion: Phase16 is NOT demonstrably materially smaller than Phase19's
measured515MiB peak; plan roughly0.6-1.0GiB process headroom, with uncertainty.
Do not treat the110-minute budget as a wall-clock upper bound. No full Phase16
probe/launch was used to produce this estimate.

Native concurrency issue found statically: Phase16's OpenMM Context omitted
CPU thread properties, so OMP_NUM_THREADS alone did not enforce the requested
limit. Fixed its CPU Context to explicitly use that value (default2), while
Reference uses no CPU-specific properties. Phase16 remains unlaunched pending
coordinator allocation after this estimate.

Phase15 default/full attempt running in `phase15/20260910T135407`; recorded
PySCF failures (qbscf: no NumPy; phase2ff: no PySCF) and explicit literature
tensor fallback. First d576 theta0 curve took552s and yielded PhiS=.8445,
PhiT=.1600; the quadrature closure deviation is not a scientific pass.
Phase19 first default/full attempt `phase19/20260910T141204`, PID41208, was
stopped after the coordinator flagged angular loss explosion at step875:
loss=2.437878e21, translation loss1.4471, angular loss9.7515e21. Verified
tracked command before stopping; Phase15 was left running. Runner recorded
exit4294967295 and561.6s. No candidate verification or MD was reached.
Peak working set observed539,885,568bytes. Training reached875/1400, not
completed, and this is an actual failed FULL-default attempt, not a full pass.

Root cause: the batch SO(3) logarithm divided nonzero skew components by
sin(arccos(-1)) when near-pi float32 matrices rounded their trace to-1.
Replaced scalar and batch logarithms with validated-rotation SciPy quaternion
conversion/as_rotvec, without target clipping. Reject nonfinite, reflection,
or materially nonorthonormal inputs; allowed1e-5 input roundoff is explicitly
checked. Principal log length is checked <=pi. Added catastrophic training
loss fail-fast, nonfinite-gradient rejection and max-target-angle logging.
30focusedtests passed6.339s, including float32 near/exactpi,5000random
rotation-pair roundtrips and invalid-rotation rejection.
Fresh default retrain pending the existing >=2GiB startup headroom gate
(latest availableRAM~1.3-1.7GiB); old training state will not be reused.

Current source potential-scan gate is12.0kcal/mol, matching the preference.
The failed141204 snapshot had12.5; its outputs are not silently relabeled.
Added an independent `--stage acceptance-audit --audit-input PATH` command
that reevaluates saved candidate observations against current gates with
input SHA256 and no measurement changes. A12.4 regression explicitly fails.

Phase15 source now uses composite Simpson quadrature. Equal-rate analytic
tests on its145-point segmented grid show closure error<2e-4 and >10x
improvement over trapezoid; malformed populations remain malformed, not
normalized. Current running snapshot still uses the original trapezoid
method; added `--stage yield-audit --audit-input PATH` to recompute yields
from its fresh saved traces without rerunning density dynamics or overwriting
original outputs. The audit reports endpoint survival, initial trace,
expected total from trace loss, Simpson/trapezoid closure errors and input
hash. Actual defaults are unequal kS=1e7/kT=1e6; the equal-rate exponential
is a regression oracle only. Neither individual-yield grid convergence nor
exact numerical trace propagation is claimed.
Initial observed peak working sets: Phase15 580,751,360bytes, Phase19 282,144,768bytes.
Phase16 queued pending a free slot.
No isolated QM audit substitutes for full enzyme evolution.

## 18:13 continuation update

P15 five full-space spin curves, MFE and ODMR saved in135407 checkpoint,
SHA256924cb245409e866a460b884684d73bb3ed1557cd269b48ebc505ab95e6cc4df7.
Own9344 was stopped at15C boundary, original failed runner status retained.
Saved-population Simpson audit does NOT establish convergence: closure errors
theta0 +0.0037793, theta30 +0.0042972, theta60 -0.0014225,
theta90 -0.0091530, Bzero -0.0094238. Remaining survival2.77e-5..7.32e-5.
Some Simpson errors exceed trapezoid errors; no renormalization or exact yield claim.

Corrected-geometry allostery continuation145816 was also stopped with authority
at18:10 after its throughput probe crossed a long wall-clock interruption,
shrinking production/window sampling below0.5ps. PeakWS254,988,288bytes,
private425,836,544bytes. This is invalid undersampling, not full acceptance.
Fixed sampling now declares300ps production +8*(12ps settle+150ps sample)
+70ps equilibration =1666ps or833,000steps per state at2fs, two states.
Throughput/wall interruptions cannot reduce these counts. Progress records
requested_steps, executed_steps, simulated_ps and separate wall/process CPU time.
Wall elapsed includes interruptions and is not presented as compute time.
Fixed late-stage KeyError caused by stale release-energy result fields.

Fresh allostery-only181243 launched18:12:43, PID41388, RAM1,811,800KiB.
Source453cb8ef86596889de09a29fae1a567d54e404411f02d5b51b3ba2fcd6f34bfd.
Machine-readable continuation_status.json at135407 points to this latest stage,
retains both stopped attempts and fresh checkpoint hash. Spin is not rerun.

P19 fresh FULL default181030 launched18:10:30 PID3804 at2,158,152KiB available.
70folds completed28s;1400-step training active, step0loss3.3068/maxangle3.139973.
Strict0.8A RMSF/12.0kcal/mol scan gates apply. P16 no attempt yet; remains queued.
36 focused/prior audit tests passed6.712s. No full scientific pass claimed yet.

At18:17 root-source acceptance additionally requires the full declared40ps
MD production, rejecting NaN or shorter durations. Future runs use fixed
requested production steps, never an early wall-budget break. Running181030
snapshot is unchanged: its candidate measurements will be post-audited with
input hash and explicit40ps provenance, without restarting training. Removed
remaining DG_act/experimental-anchor wording from future log/PDB labels.
36 tests passed6.776s including short-duration acceptance rejection.
P15 measured11.0ns/day wall throughput: full3.332ns would take~7.3h at this
contended rate, not a30min guarantee. No further launches without coordinator
acknowledgment plus RAM gate; existing15/19 explicitly kept running.

## CURRENT HANDOFF — 2026-09-10 18:38 (supersedes earlier running estimates)

All paths below are relative to this repository worktree unless absolute.
Campaign root: `../phase1-19-rerun-20260910`.
Only owned phase source/test/note files were edited, plus the explicitly
requested campaign `continuation_status.json`. No staging, commit or push.

### Phase15 — KEEP RUNNING; incomplete

- Active corrected fixed-sampling allostery PID **41388**, runner session58843,
  snapshot `phase15/20260910T181243`, command
  `PRIMARY_PYTHON -u rerun_phase1_19_campaign.py --phases 15 --stage allostery`.
  Child uses phase2ff Python, CPU two threads. Do not duplicate or kill it.
- Latest recorded FAD_oxid progress:92,000/833,000steps,184ps; FAD_radan
  not yet started. Progress is in
  `phase15/20260910T181243/results_phase15/md_FAD_oxid/sampling_progress.json`.
- Fixed per-state sampling:70ps equilibration+300ps production+
  eight*(12ps settle+150ps production),833,000steps at2fs=1666ps.
  Two states require3.332ns total. Current measured wall throughput11ns/day
  implies about **7.3 hours**, not a promise; contention/suspension can extend it.
  Wall time is not compute time. No throughput-dependent sampling reduction.
- Observed peak working set256,491,520bytes, private427,876,352bytes.
- Original fresh spin/MFE/ODMR checkpoint remains immutable in135407,
  hash924cb245409e866a460b884684d73bb3ed1557cd269b48ebc505ab95e6cc4df7.
  Both controlled-stop attempts retain original failed statuses/logs.
  `phase15/20260910T135407/continuation_status.json` points to active181243,
  with PID, command, snapshot source hash, checkpoint hash, log and progress paths.
- After actual allostery completion: review real sampling/WHAM convergence,
  assemble this fresh stage with135407 via `--stage assemble --audit-input
  FRESH_SPIN_JSON --allostery-input FRESH_ALLOSTERY_JSON`, then render figures.
  Do not rerun spin or call this a single corrected full-default run.
  Saved spin quadrature remains unconverged (up to0.94% closure error);
  explicit literature hyperfine fallback remains, not native PySCF success.

### Phase19 — full default attempt completed, ZERO accepted designs

- Snapshot `phase19/20260910T181030`, PID3804 **exited0** at18:36:27;
  runner elapsed1557.5s (~26min). No active19 process remains.
- Full70fold distribution (~28s),1400trainingsteps (~14min),fivegenerations,
  **30 candidates** completed. Equivariance errors7.6e-6/1.1e-6<.005.
  Step875loss1.8712, no prior near-pi explosion. PeakWS observed578,416,640bytes.
- All30 were rejected by the unchanged0.30A all-anchor gate at~3.10A.
  No candidate MD or enzyme QM/MM executed. One fresh uncatalyzed reference
  scan executed, sampled potential peak23.16kcal/mol; this is not an enzyme
  result or validated activation free energy. No champion/extra champion MD.
- Master: `phase19/20260910T181030/results_phase19/phase19_results.json`.
  Incremental/final candidate records: same directory `evolution_progress.json`.
  Independent `phase19_acceptance_audit.json` hashes the raw master and
  applies0.8A RMSF,12.0kcal/mol sampled peak AND full40ps production duration:
  **30 records,0 accepted**. No changes to original observations or snapshot.
- Actual completed runtime26min is for this all-rejected path. It is NOT an
  estimate for an accepted path with up to10candidateMD/fiveenzymeQM scans.
- The current snapshot has **NO trained-model or70fold-library checkpoint**:
  no torch.save/np.savez persistence existed, only live process memory.
  Therefore generation-only resume from this run is unavailable. Never reuse
  historical outputs or claim that the completed1400step model was saved.

### Remaining Phase19 geometry work / next attempt

- Root source has additional repairs NOT applied to completed181030: planar
  sp2 branch placement and correct angles (old helper used the plane normal
  and supplementary angle), separate tetrahedral branch construction,
  ASP branch reference usesCB rather than its own centerCG, and connectivity-
  based Trp ring naming/closure replacing an incorrect pyrrole/six-ring walk.
- Analytic bond length/angle/planarity and named Trp ring closure/unique-atom
  tests pass:25owned focused tests4.373s. Latest combined suite
  `-m unittest -q test_phase15_16_19_rerun test_phase14_17_audit test_phase18_19_audit`
  passed **38tests in5.076s**; owned-file `git diff --check` clean.
- **Systematic3.10A anchor-offset cause is still under investigation**.
  These definite geometry repairs do not yet prove it fixed. Do not loosen
  the gate, graft invalid geometry merely to match targets, or claim acceptance.
- Planar GLU repair changes theozyme geometry and thus fold/training targets.
  After target consistency and geometry tests, a new fresh training attempt
  is required (also no saved current checkpoint exists). Add durable model+
  fold-library checkpoint and source/target hashes before that next attempt;
  checkpoint persistence/resume is **not yet implemented**.
- Root future MD now uses fixed declared steps and reports requested/executed
  steps/simulated_ps. Root acceptance rejects<40ps. Completed181030 used its
  older wall-budget loop but never reached candidateMD; strict post-audit is
  the explicit provenance layer, not a retroactive source change.

### Phase16 / allocation — UNSTARTED, MANUAL ACK REQUIRED

- **No Phase16 snapshot, process, or full attempt exists.** No smoke/full claim.
- Static estimate0.6–1.0GiB peak, uncertain, not a measured bound;
  ~643025solventbeads,344bytes/bead base buffers plus temporaries and maps.
  Existing110min budget is not a runtime guarantee.
- **Do not launch Phase16 or any replacement/new15/19 independently.** Main
  must explicitly acknowledge allocation first, followed by a live RAM check
  immediately before launch. P7 needs2.5GiB and main17 may release RAM;
  avoid independent-check launch races. New19 requires>=2GiB. No new job
  below500MiB headroom. P6 Context transition hold was cleared previously.
- Once allocated, Phase16 command:
  `PRIMARY_PYTHON -u rerun_phase1_19_campaign.py --phases 16` (two threads,
  primary environment chosen by runner). Capture actual peak and real gates.
- Main's background-work approval question is pending user response. Existing
 15 remains running by explicit instruction; this note is a handoff, not new
  unattended-launch authorization. Do not create an automation or hardware action.

## Post-main-unfreeze updates (after main batch1dd9768; new edits uncommitted here)

Only owned19source/test/notes changed after the explicit unfreeze. Main owns
all staging/commits/pushes. Existing15 remains running;16 remains unallocated.

- Fixed Trp indole attachment axis: use the inward bisector of both CG ring
  bonds, not CG->NE1. The old mapping made NE1 collinear with the CB-CG
  torsion axis (cosine0.99999999999); chi2 sweep0..180deg moved NE1 only2e-5A,
  incorrectly eliminating a placement degree of freedom. Regression checks
  chi2 moves NE1>0.5A and inverse fitting to a known attainable target reaches
  <1e-6A. Named ring connectivity/planarity/bond checks remain passing.
- Fresh single geometry probe `phase19/20260910T191530` completed35.2s,
  scope **probe only**, no training/MD/QM. RMSD2.8434739A still FAILS0.30A.
  Per-anchor attribution: GLU CA/OE1/OE2=0,TRP CA=0,ASN N~8e-17A,
  **TRP NE1=6.9650602A**. Thus the systematic offset is confined to Trp
  placement, but axis repair alone is insufficient. No gate relaxed.
- Then fixed stale canonical CB/CG aiming and chirality-reference builders
  still attaching CB to C while actual packing attaches CB to CA. New
  regression equates canonical world-space CB/CG directions to actual atoms.
  Removed fictitious terminal O-O/O-N branch extension steps before replacing
  branch positions; this avoids invalid temporary NeRF geometry.
  These latter fixes are unit-tested but the next fresh geometry probe is
  **pending RAM**: blocked twice at1,004,012/1,014,740KiB available against a
  conservative1.2GiB probe launch floor. No automatic heavy-training request.
- Added `save_training_checkpoint` immediately after completed training and
  equivariance audit, before generation1. Files: `trained_flow.pt` (weights
  and Torch RNG), `training_folds.npz` (x,R,mask,fullrodN/CA/C/O arrays), and
  completion marker `training_checkpoint.json` with model/data SHA256,
  source SHA256/config/architecture/completedsteps/foldcount/versions,
  generator RNG, globalNumPy RNG and Python RNG. Marker written last;
  existing files cannot be overwritten. This is for FUTURE runs, not a
  recovered181030 model. A generation-only load/resume CLI is not implemented.
- Real molecular Torch checkpoint roundtrip test passed1test3.055s: weights,
  Torch RNG, dataset shape/data, file hashes and no-overwrite checked.
  Latest primary combined command
  `-m unittest -q test_phase15_16_19_rerun test_phase14_17_audit test_phase18_19_audit`
  ran40tests5.551s:39passed,1Torch-dependent test skipped in primary (passed
  separately in molecular). Owned-file diff check clean.
- Next: after safe RAM headroom, run the corrected single geometry probe,
  investigate residuals if any; do not start training until actual geometry
  meets unchanged criteria. Then request MAIN ACK and >2GiB live available
  RAM for a new full attempt with changed targets and durable checkpointing.

## SECOND-BATCH SOURCE FREEZE READY — 2026-09-10 19:24

- Latest corrected-source single geometry probe:
  `phase19/20260910T192253`, command
  `PRIMARY_PYTHON -u rerun_phase1_19_campaign.py --phases 19 --stage geometry-probe`.
  LaunchRAM1,600,148KiB; exited0 in28.3s. SourceSHA256
  `423e30eb7ff6c09c714021f77452ecdc449598d41020d2773d8a806bf773b847`.
- **Geometry gate FAILED**: RMSD1.013322739A>0.30A. NE1error2.482123656A;
  GLU CA4.44e-16A, GLU OE1/OE2=0, TRP CA=0, ASN N8.39e-17A.
  No division warning in this probe log. Staticclashes569,
  Ramachandranfraction0.8563. These static diagnostics do not establish MD
  stability, acceptable chemistry, or a validated enzyme design.
- Per-anchor attribution is exact for the two NEW probes: all nonzero
  deviation resides in Trp NE1. Axis-only probe191530 gave2.84347A overall /
  6.96506A NE1; consistent canonicalCB/CG directions reduce this to1.01332A /
  2.48212A. This demonstrates definite contributions from the axis/reference
  bugs, **not complete resolution of the old~3.10A systematic offset**.
  The original181030 full attempt did not save per-anchor errors/candidate
  coordinates or model/library checkpoints; its exact per-anchor attribution
  cannot be retroactively claimed from these changed-source probes.
- Remaining investigation: actual Trp CA/backbone orientation and target
  reachability under physical chi1/chi2 geometry; do not translate/graft a
  strained ring merely to satisfy the target, and do not relax0.30A.
- Latest combined regression40tests6.691s:39pass/1primaryTorchskip.
  Latest REAL molecular Torch checkpoint roundtrip separately passed1test
  in3.845s at19:24; verifies weights,RNG,dataset,hashes,no-overwrite.
  TestfileSHA256
  `d712ae127eb159e3a4dd9bd6f8046da0f85fa47627222865c6d5ebae2dbf927d`.
- New edits are unit/probe-tested, **NO new full training/evolution validates
  this source**. No request for a costly training allocation while geometry
  fails. Main owns secondcommit/push; source/test freeze-ready now.
- P15PID41388 remains running: latest272,000steps/544ps in firststate,
  requested833,000steps/state. P16 still unstarted and unallocated.

### Final acknowledged probe / updated freeze — 19:26

Main explicitly authorized one corrected geometry probe with live>=1.2GiB.
Fresh snapshot `phase19/20260910T192525` launched at1,544,276KiB available,
completed30.8s with **exit2**. Geometry JSON is preserved despite scientific
failure: RMSD1.013322739A, Trp NE1error2.482123656A, otheranchors~0.
Root CLI now exits2 after saving any false geometry gate; prior probe exit0
statuses remain unchanged historical process records, never relabeled passes.
SourceSHA256 `7ca9987f8f2085582d47b3535bbd1a096998176b8c6e55a3e2bb4b0b0e73286d`.
Header now describes fixed declared MD production steps, not wall budgeting.
Latest combined41tests9.018s:40passed/1primaryTorchskip; real molecular
Torch roundtrip independently passed3.845s before this CLI/header-only edit.
No new full training or evolution validated these fixes; geometry still fails,
so no training allocation requested. SOURCE freeze-ready again after this note.
