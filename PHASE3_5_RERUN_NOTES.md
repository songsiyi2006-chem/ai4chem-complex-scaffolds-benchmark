# Phase 3–5 campaign status — 2026-09-10, TS2b causal diagnosis update

## Current outcome

| Phase / attempt | Scope | Outcome |
| --- | --- | --- |
| 3 / 20260910T144409 | Default full attempt | **Running**, original PID14832; production began18:57:05;26 complete frames /13000 of100000 production steps at20:34:26 |
| 4 / 20260910T141256 | Default full attempt, 500 CI steps | Exit1 after774.5s; best force0.152173 >0.05 eV/A; **not accepted** |
| 4 / 20260910T181601 | Bounded continuation, 150 steps | Exit1 after219.1s; best force0.150663 >0.05 eV/A; **not accepted** |
| 5 / 20260910T142809 | Default full attempt before energy correction | Exit1 after917.6s; TS2b only2 converged frames; ModuleA incomplete |
| 5 / 20260910T182103 | Fresh default full attempt with physical-SP correction | Exit1 after846.4s; TS2b still only2 converged frames; **not accepted** |
| 5 / 20260910T201431 | Fresh default full attempt with physical-SP and TS2b pose correction | Exit1 after805.375s; TS2b4/4 converged; ModuleB250K ODE failure; stereo invalidated by Windows cache collision |

No further Phase4/5 retries unless new causal evidence appears. Existing failures
remain failures; no convergence thresholds were relaxed, no rates invented, no
historical outputs counted as fresh. No staging, commit, push or shutdown.

Campaign root (all attempt paths below are relative to this absolute directory):
`C:/Users/HUIWEI/Documents/Codex/2026-09-08/https-github-com-songsiyi2006-chem-ai4chem-2/work/phase1-19-rerun-20260910/`

Each attempt contains `run.log`, `rerun_status.json`, Python source snapshots and
source SHA256 records. Process exit and scientific acceptance are distinct.

## Phase3: preserve the existing full run

Live attempt: `phase03/20260910T144409/`; original molecular PID14832, runnerPID15352.
No duplicate or restart on18:08 resume. The approximately3-hour sleep gap is not
production runtime. Do not silently shorten5000 equilibration or100000 production
steps. Actual production timestep is1fs, so default production is100ps, despite
the legacy header's200ps wording. Target analysis is40 MMGBSA frames.

Heating completed18:25:52; full5000-step equilibration completed before production
began18:57:05. At19:15:20.659271, DCD contains6 complete frames (500steps/frame),
3000 of100000 production steps. Frame5 completed19:12:14.126007; the next500steps
took186.533264s, or2.680487steps/s. Average from production start is2.738078steps/s.
Estimated total production duration10.145h, remaining9.841h as of19:15:20;
projected production finish approximately2026-09-11 05:06 China time, **assuming
continued awake operation and similar load**. MMGBSA/figure analysis is additional.
Only the production interval is used; earlier heating and approximately3h sleep
are excluded. No timestep or step-count reduction was made.

Measurement method: count complete DCD Fortran records/header frames and use
the last completed-frame file timestamp, checking a stable file size/read.
Direct MDTraj access to the live writer's file raised a read error; an untouched
temporary snapshot read successfully (5frames/2785atoms, all coordinates finite),
and byte-level records were complete. No evidence of trajectory corruption was
found by that check. Do not confuse a live-file sharing/read issue with MD failure.

Fresh7RPZ download/curation, Vina docking and MACE force comparison completed.
Top docking scores: -7.95, -7.47, -7.34 kcal/mol. MD system:2785atoms/170residues,
Amber14SB protein plus **Sage2.1** with NAGL AM1-BCC charges for ligand/GDP.
This is not GAFF2. Binding-energy interpretations still require the MMGBSA checks;
no binding free energy or successful full MD result is claimed while running.

Original launch explicitly set `$env:OPENMM_CPU_THREADS='2'` before the runner;
the runner inherits it. OMP alone was not assumed to control OpenMM. At18:18,
PowerShell observed13 total threads,2Running/11Wait, CPU535.140625s, RSS706.6MiB.
No oversized OpenMM pool was demonstrated. The same original process is retained.

Follow-up monitor `finish-existing-phase-3-md` is **PAUSED pending explicit user
confirmation**. Main confirmed the pause through the app update and TOML readback.
Its creation was not authorized; the earlier ACTIVE statement is superseded.
Do not reactivate, duplicate or create another automation without new user approval.
Pausing this monitor does not stop the existing Phase3 process, which is untouched.
Subsequent user confirmation authorizes background continuation through the MAIN
task's existing20-minute heartbeat, reported ACTIVE by main. This owner has not
created, updated or reactivated any automation. The separate Phase3 monitor above
remains paused; main owns coordination. Phase3 completion is not claimed.

Source fixes after this snapshot: OpenMM H-relax API corrected to
`LocalEnergyMinimizer.minimize`; success/plot labels corrected from GAFF2 to Sage;
bond tuple explicitly length/spring and fallback X-H constraint uses length.
The running snapshot is unchanged: it logged H-relax failure, so its force
comparison retains that caveat; legacy figure/provenance labels need this note.
These source fixes do not retroactively validate the running snapshot.

## Phase4: full failure and bounded refinement

Full `phase04/20260910T141256/results_phase4/`:
535 logged steps include500 CI steps. Best band restored from step41;
fmax0.15217299888877053 eV/A versus unchanged0.05. `converged=false`,
`all_stages_ok=false`, `stage3_ts={}`. No Hessian screening, thermochemistry or rate.
Figure traceback is the expected validation guard, not a new independent bug.
Source prose now explicitly says failedNEB gate stops stage3; the completed
snapshot is preserved unchanged.

Bounded `phase04/20260910T181601/results_phase4/` starts from that campaign's
fresh best band, records input hashes in `refinement_provenance`, and is labelled
`bounded_refinement_not_default_full`. FIRE maxstep0.005A, dt0.02, hard150steps;
unchanged0.05gate. Best force0.15066262106165354 eV/A still fails. No TS/rate.
No more refinement is planned without new causal evidence.

Earlier attempts retained for provenance:134734 stopped after57.3s upon finding
warm-start skipped xTB CI and MOL exports held pre-xTB coordinates;134845 stopped
after1384.7s when CI-BFGS diverged (.1823 atstep35 to4.3996 at75). Corrected exports
and CI execution; changed CI optimizer to FIRE and reset best history at climb
switch. Per-image exact-state caching avoids duplicate subprocess calls and
returns force copies to protect cached values from NEB projection.

Isolated xTB finite difference at IDPP midpoint atom5/x, h1e-4A: analytic force
4.9957827814 versus central difference4.9958894243 eV/A; relative error2.1346e-5,
net force5.5667e-14 eV/A. No force-unit/sign inconsistency found in this check.

Fresh endpoint assets for Phase6/7 (sent to main for relay to Epicurus):
`phase04/20260910T141256/results_phase4/`. Hashes match corrected134845 exactly.

| File | SHA256 |
| --- | --- |
| reactant_3d.mol | E52CD7FBC8028BFBDF97841B04C1C87EC83B56C75097A894EC82E046E0DB8C78 |
| product_3d.mol | 7F3952E8F11112FC4CF3537FA66B04891A9D74E5B974C0B657EC74C5D830B0CD |
| images_idpp.xyz | 326423CE9E9159B487DA22863B2321AAF7FB79F3CB371CCF47FDCBE7CDAD9C48 |
| neb_final_path.xyz | 48DFF5D033981EC75278AA5BADE2D151376C39060B59A507E8F56ED07F473943 |
| ../run_phase4_reaction_mechanism.py | F1B21C24C8861C4FDD2587E5A6F3166DCB3C9B305337BE9D0B486ED9FCCE1BED |

Nine IDPP images; MOL vs optimized endpoint coordinate differences4.961e-5A(R),
4.742e-5A(P), consistent with MOL rounding. Starting structures are not validatedTSs.

## Phase5: restraint bias, correction and fresh rerun failure

Durable probe JSON with exact command, stdout, geometry/source hashes and limits:
`phase05/20260910T142809/restraint_probe_review.json`.
Geometry: `results_phase5/cache_phase1-5-validation-v1/cx_TS1.json` in that attempt,
SHA256 `eb282294bfc89354cf03d0b446d6c4840b8973061f1c0b2807b210a4e952d9d1`.
Reproducible probe source snapshot: `phase04/20260910T181601/test_phase3_5_rerun.py`,
SHA256 `0845E2C229A12F92F5700EC61B929AC103F63821C62FAE0AC6FA76F67AC1B489`.

Two single points on exactly the same saved refinedTS1 geometry, no optimization,
GFN2/charge0/no solvent, atoms10/2 target2A, constraintforceconstant0.8:

| Evaluation | Energy(Eh) |
| --- | --- |
| Unconstrained | -98.326353289187 |
| Constrained | -98.325879630601 |

Difference **+0.29722525018822155 kcal/mol**. Both displayed `total energy` and
`TOTAL ENERGY` include restraint potential. The test does not assume a parser
separates physical/restraint components. Original individual scan geometries
were not persisted; this saved refined geometry cannot retrocorrect every old frame.

Fix: constrained optimization now returns a separate unrestrained same-geometry
single-point energy, with identical charge/GFN2 and no guessed subtraction.
Cache version `phase1-5-validation-v1-physical-sp-v2` invalidates old hydration.
The constrained geometry and resulting energy remain a TS proxy, not verifiedTS.

Corrected fresh full `phase05/20260910T182103/` used that cache version and failed
at the unchanged minimum3 converged frames: TS2b supplied2. ModuleA incomplete;
B/C notrun; Ddisabled/notrun. No rates or acceptedTSs. TS1 scan4/6frames, one
imaginary mode; I1CatM one imaginary mode (not verified minimum); TS2a three
imaginary modes (not first-order saddle). No unrestrained gradient stationarity
or IRC connectivity demonstrated. One imaginary frequency alone is insufficient.

Fixed scriptSHA256 `DBB02293C06A14746D5479CB7ED4D842D154EFFA0A9E58014E966CB565ADE206`;
failure JSONSHA256 `5704BE044CB5F7A1670210D5D76DC53C316EB4589DC4CCECDFF713E32D8B788D`.
This rerun demonstrates the physical-energy fix but does not cure geometry
nonconvergence. No further Phase5 retry without new evidence.

### New TS2b causal evidence and algorithm correction (freeze-ready)

Read-only AST reconstruction of ModuleA using ONLY the fresh182103 cached
intermediates identified a definite assembly defect. `assemble_face` returned
catalyst-local O/P indices after merging substrate+catalyst. TS2b used local O=0
as `pos[0]`, which is a substrate **carbon**, not the actual phosphate oxygen at
global index23 (23 substrate atoms;60 atoms total). Migrating H is index15 and
C3 is index2, all zero-based. The resulting midpoint H had a0.768139171A contact
to another atom. The erroneous C3-to-selected-carbon distance was2.599352889A;
actual C3-to-phosphate-O distance was6.313371135A. Merely adding the index offset
would still put the proton halfway across an implausibly distant donor/acceptor.

Implemented solution: return global catalyst indices from both ion-pair
assemblers; for TS2b retain the substrate's actual C3-H bond and rigidly orient
and translate the phosphate to the C-H donor, with initial H...O separation1.65A
and P-O pointing toward H. This is an INITIAL pose choice, not a distance
acceptance tolerance or stationary TS. The helper validates atom elements and
fragment membership and preserves both fragment internal geometries. Read-only
reconstruction of the corrected source yields C-H=1.091738656A, H...O=1.65A,
minimum cross-fragment separation=1.65A (the intended H/O pair). No xTB, Hessian,
optimization or ML import was used for these geometry checks.

Both previous full attempts failed the1.60A and1.40A TS2b optimizations and kept
only1.25A and1.10A, below the unchanged3-frame requirement. The index/pose defect
is confirmed, but its quantitative contribution to optimizer nonconvergence is
NOT yet isolated by an engine comparison. Earlier errors retained only a timing
tail and deleted temporary directories, so the precise original stopping reason
(iteration limit versus other numerical failure) cannot be reconstructed from
those logs. New `XTBFailure` preserves full stdout/stderr and requested engine
files; every scan frame now saves constraints, input coordinates, converged
coordinates/physical-SP energy or failure evidence under
`results_phase5/scan_diagnostics/`. This enables the next bounded attempt to
diagnose residual failures instead of repeating an opaque error.

New cache version `phase1-5-validation-v1-physical-sp-ts2b-pose-v3` prevents
reusing old incorrectly assembled results. Physical unrestrained SP energy
evaluation remains mandatory. The3-frame guard and all convergence/TS validation
limitations are unchanged; no kinetic or saddle success is claimed.

Evidence input hashes (SHA256), in
`phase05/20260910T182103/results_phase5/cache_phase1-5-validation-v1-physical-sp-v2/`:
`sp_Rprot.json`: `CA1407564B9E3D00244BC6FBFC3C0766643D211F1C0354CE498C428C278AAC7C`;
`sp_CPAneg.json`: `46E41D15CDDD571B19B869BA01BADA1167AB616230C65C83A093AFAD81702262`.
Corrected sourceSHA256:
`79F67387E93D25604647B109EBA09202822B8AC91FB8DC293902BA462130EBED`.

Current batch: **27 regression tests passed in2.194s**, including wrong-fragment
index rejection, rigid-pose geometry preservation and full failure-evidence
retention with the unchanged minimum3-frame gate. No new xTB was launched during
this update (reported freeRAM below0.8GiB). New fullPhase5 remains pending main's
memory/compute allocation; test success is NOT a molecular rerun. Phase4 still
fails0.05eV/A and has no newly established cause or authorized retry. Original
Phase3 PID14832 remains alive and untouched (last observed CPU2140.94s,
working set745934848bytes). Main owns commit/push; this owner performs neither.

### Allocated fresh full Phase5 started20:14:31

Main explicitly allocated the low-memory Phase5 attempt after its128+10 regression
executions passed. Before launching, process inspection found only existing
Phase3 PID14832 and no Phase5 duplicate; free physical RAM was1022564KiB.
Command: `PRIMARY rerun_phase1_19_campaign.py --phases 5`, no extra arguments.
Fresh attempt `phase05/20260910T201431/`; unified execution session83162;
runner PID41292, molecular PID28064. Runtime metadata records all five thread
limits as2. No historical cache copied, no module/sample reduction, ModuleD
remains default opt-in/not requested. Started ModuleA species optimizations and
Hessians; reached A.4 at20:15:16. Full/scientific outcome still pending.
P3 is untouched, and no other new heavy task was started. Source/tests remain
frozen; this launch record is the only post-freeze notes edit.

### Source-only v4 cache repair discovered during201431 (running snapshot unchanged)

At20:22:11 the fresh run logged `I1CatM` followed immediately by
`cache hit: cx_I1Cat_m`. Windows case-insensitive filenames alias
`cx_I1Cat_M.json` and `cx_I1Cat_m.json`. Only the uppercase file exists in the
directory listing; reading the lowercase path returns the same `face_deg=35.0`
record. Thus **201431 stereoselectivity is not acceptable**, regardless of later
exit code. Main explicitly approved retaining the current attempt to inspect
TS2b and subsequent gates, not restarting or modifying the snapshot.

Source-only repair uses `cx_I1Cat_face_plus35` / `cx_I1Cat_face_minus35` via
shared `FACE_CACHE_KEYS`, including result hydration. Separate exact-key reuse
also skipped the intended unlocked optimization: `comp_I1Cat_unloc` previously
requested already-populated `cx_I1Cat`. It now uses `cx_I1Cat_unlocked` in both
calculation and hydration. Hydration of the deliberately shared TS2a proxy now
loads its actual `cx_TS2a` file, rather than nonexistent M/m-suffixed files.
Cache version is now `phase1-5-validation-v1-physical-sp-ts2b-pose-cache-v4`.

All cache entry points reviewed: literal complex keys, `sp_{name}` for the
eight species R/P_target/P_elim/Cat/H2/P_poly/Rprot/CPAneg, direct species reads,
and opt-in ModuleD `dface_WIN_+1/-1` and `dface_BASE_+1/-1`. No additional
case-only filename pairs were found. In-memory dictionary keys M/m remain
distinct in Python and need not be renamed. Scan diagnostics use distinct
TS1/TS2a/TS2b labels and numeric frame files.

28 focused/existing regressions passed in2.818s; the new test compares casefolded
keys, executes independent cache computations/readbacks and forbids obsolete
case-only literals/hydration keys. Future sourceSHA256:
`F41697E98B2CA9996A581ACE2E329C8E9C6F5E079674958D1EDF1860B8021B90`.
This differs from the immutable201431 v3 snapshot. No v4 full run has started;
one will require main coordination. No gates, modules or sampling were reduced.

###201431 final outcome and source-only ODE repair

The retained process exited1 after805.375s (session83162 completed). TS2b now
converged at all four original targets1.60/1.40/1.25/1.10A, directly resolving
the previous2-frame blocker with unchanged convergence gates. Selected maximum
is the endpoint1.10A, n_imag=2: NOT a stationary first-order TS or an established
barrier maximum. TS1 has4/6 scan frames and n_imag=1; TS2a3/5 and n_imag=3.
ModuleA completed its calculations but remains exploratory/stereo-invalid.
At298.15K BDF failed and Radau completed; at250K all three solvers failed.
ModuleC was not run, ModuleD remained disabled. No full Phase5 acceptance.

Preserved artifacts under `phase05/20260910T201431/` (SHA256):
- `results_phase5/phase5_results.json`:
  `95aa3cc1eaec309695e3cb4db15fafce3ed7fd0b11d371710d4105c1005017d8`
- `run.log`: `3b2b9e5eab00cb48a9a2bc1ecd40f3475f15d6973e774bc5efdb03d4e93a614c`
- `results_phase5/cache_phase1-5-validation-v1-physical-sp-ts2b-pose-v3/cx_TS2b.json`:
  `dba7869e5a678dc4e1567ece2ac58b15574780236616df118aa0eff8dd08b5c5`
- All four TS2b input/output geometries, constraints and physical-SP energies
  are in `results_phase5/scan_diagnostics/TS2b-scan_dC3-H_/000.json` through003.

ODE diagnosis uses ONLY this saved G/dG; no molecular rerun. The physical RC
electronic binding difference is+94.144024886kcal/mol and binding free energy
is+101.669457836kcal/mol. The existing expression gives finite
k_off=7.544031540045538e97s^-1 at250K and fastest time1.325551192e-98s.
At initial time1e-9s the floating-point spacing is2.067951531e-25s: BDF fails
after6 RHS calls because it cannot represent the initial step. Starting the
same autonomous equations at an internal zero time permits those tiny initial
steps and subsequent adaptive growth. A diagnostic with no RHS-call bound was
stopped after becoming unresponsive; only its verified probePID22684 was stopped,
never P3. Bounded replays then established the time-origin cause reproducibly.

Source fix translates internal time by1e-9s, preserving the original elapsed
duration,160 output times, dense-output coordinates, rate constants and
rtol1e-6/atol1e-14. No rate cap, timescale deletion or looser tolerance. BDF250K
now succeeds with382 RHS calls; independent Radau agrees. The analytic Jacobian
matched finite differences at positive concentrations (maxerror2.48e-10 in a
scaled-rate diagnostic) but not negative Newton trial RC (maxerror5), because
RHS clipped atzero while Jacobian did not. Jacobian now differentiates that same
clip, using the physical right derivative atzero.

Finite/nonnegative-rate validation is explicit. The bimolecular association,
reassociation and dimerization constants have units M^-1s^-1 (r1/r9/r10); other
constants are s^-1. The k_off expression explicitly multiplies the1M standard
concentration, numerically unchanged. JSON now uses `rate_constants` plus
`rate_constant_units`, replacing misleading `rate_constants_s`; no internal
consumer referenced the old field. Eigenvalue diagnostics no longer claim all
modes negative or prove stability: conservation imposes zero modes, and extreme
scale separation compromises numerical eigenvalues.

**30 tests passed in4.633s**, including actual SciPy replay at250..350K (11points)
plus298.15K, unchanged250K k_off,160times, dense output, finite/nonnegative
concentrations, catalyst/substrate conservation, and independent250K Radau/BDF
agreement. Scoped `git diff --check` passed. This is a targeted numerical replay,
not a new full Phase5 run or validation of gas-phase binding/TS proxies. Current
sourceSHA256 `CE58E232E29AF8324A6289890346C0E2F2587933157CCE448C71F2FF59144DA9`.
Main reported cachev4 committed/pushed as e4c3031; the subsequent ODE patch/tests
are source-only and await main review. No second full run was launched.

P3 latest read-only observation: original PID14832 alive, DCD26 complete frames,
869820bytes, last frame20:34:26.840047;13000/100000 production steps. The previous
19:15 ETA is load-dependent, not a deadline; prior sleep remains excluded.

## Commands, dependencies and verification

### 2026-09-10: default Phase5 A/B/C completed; drawing repaired separately

Fresh attempt 215504 finished Module A at 22:16:06, B at 22:16:10 and C at
22:42:25. The default correctly skipped optional exploratory Module D. It then
failed in figure 2: a label at `min(yield_percent)-2` was outside the tiny
near-zero-yield axis and tight bounding-box rendering requested a canvas
599,229,264,149 pixels tall. The label now uses data-x/axes-y coordinates.
The optional designed-panel baseline ee was already percent and its erroneous
extra factor of 100 was also removed. No calculated concentrations were changed.

`recover_phase5_figures.py` copied saved A/B/C evidence into the separate
`figure_recovery_20260910T225557` directory, hashed original inputs, rendered
three PNGs (largest dimension 4011 pixels) and generated both reports with
exit codes zero. Original result files and failed process record remain intact.
Figure 2 was visually inspected. Recovery does not establish scientific
acceptance: TS2a/TS2b retain three/two imaginary modes; calculated 298 K target
yield is approximately 4.65e-16 (fraction); the catalyst complex-energy change
is not a verified activation-barrier reduction. Do not relabel this as a
successful catalyst or a validated transition-state network.

A real near-zero-yield PNG regression prevents the oversized-canvas failure.
Full regressions at 225623 passed 140+23 executions, no skipped tests; one
Torch test is executed in both groups, so this is not 163 unique tests.

Source directory: `work/phase1-5-fixes`. PrimaryPython:
`C:/Users/HUIWEI/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe`.
MolecularPython: `C:/Users/HUIWEI/miniconda3/envs/phase2ff/python.exe`.
Molecular PATH must prefix its `Library/bin` and environment root.
Runner sets OMP/MKL/OPENBLAS/NUMEXPR to2; parent explicitly sets OPENMM_CPU_THREADS2.
No new phase below500MBfreeRAM; later user authorization allowed small native
probes alongside P3. Phase5 representative imports measured82.832MiBRSS, native
xTB sampled33.3MiB; no ML imports in Phase5 or boundedPhase4.

Commands actually used (each phase separately, with the environment above):
`PRIMARY rerun_phase1_19_campaign.py --phases 4`, `--phases 5`, `--phases 3`.
BoundedPhase4 additionally used `--engine xtb --refine_from` pointing to the
absolute `phase04/20260910T141256/results_phase4` directory. It is not a full run.

Vina CLI1.2.7 installed locally at `work/phase3-deps/vina-1.2.7/vina.exe` and
passed `--version`; parent passed absolute path through VINA_EXE.
Official asset: https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_win.exe
Size1233920bytes; SHA256 `E0C4B2715E0C1A74F6E92D0F3BE0328AC97542EAFBC111E6B1EFAD897A73CCE5`.
Official release API supplied no digest/checksum asset: this is local provenance,
not independent publisher-checksum verification. MACE loaded successfully in P3.

Latest regression command: `PRIMARY -m unittest -q test_phase3_5_rerun test_phase1_5_audit`,
with PYTHONPATH `work/phase24-deps;work/phasefix-test-deps`, OpenMMthreads2.
**25 tests passed in2.635s**. Includes real OpenMM bond/constraint units and
duplicate prevention, physical-energy selection, failed-engine rejection,
imaginary-mode handling/cache isolation and existing scientific guards.
These tests are not a substitute for full molecular calculations. Scoped
`git diff --check` also passed. No unrelated source files were edited by this owner.
