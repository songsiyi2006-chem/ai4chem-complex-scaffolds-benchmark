# Phase 3–5 campaign status — 2026-09-10, updated 18:38

## Current outcome

| Phase / attempt | Scope | Outcome |
| --- | --- | --- |
| 3 / 20260910T144409 | Default full attempt | **Running**, original PID14832; heating completed18:25:52; 5000-step equilibration active; production not yet observed |
| 4 / 20260910T141256 | Default full attempt, 500 CI steps | Exit1 after774.5s; best force0.152173 >0.05 eV/A; **not accepted** |
| 4 / 20260910T181601 | Bounded continuation, 150 steps | Exit1 after219.1s; best force0.150663 >0.05 eV/A; **not accepted** |
| 5 / 20260910T142809 | Default full attempt before energy correction | Exit1 after917.6s; TS2b only2 converged frames; ModuleA incomplete |
| 5 / 20260910T182103 | Fresh default full attempt with physical-SP correction | Exit1 after846.4s; TS2b still only2 converged frames; **not accepted** |

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

Latest known stage: heating completed18:25:52, equilibration5000steps underway.
ProductionDCD had no frames at the last review. Measure steps/sec from production
progress once it begins and exclude sleep from elapsed-time estimates.

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
No automatic follow-up is currently authorized or promised, and Phase3 completion
is not claimed.

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

## Commands, dependencies and verification

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
