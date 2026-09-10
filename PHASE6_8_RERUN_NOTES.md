# Phase 6-8 fresh rerun preparation (2026-09-10)

## Resource and sampling follow-up (20:08)

P7 PID 38388 is still the original gated orchestrator (~196 MiB RSS,
170 MiB private); no Psi4 child has launched. Available physical memory
was 1.20 GiB on this check. P8 has not started and remains sequential after
P7. No gate was lowered, process duplicated, or automation created.

The 2.5 GiB launch threshold is a conservative campaign allowance, **not
a measured minimum**. `QC_MEM="2 GB"` means 2,000,000,000 bytes (~1.863 GiB)
for major Psi4 data structures, not a total-process cap. The remaining
~0.637 GiB is overhead/headroom. Psi4 documents that this memory is per job,
not per thread: reducing two threads to one does not halve the allocation.
Sources: [memory specification](https://psicode.org/psi4manual/master/psithoninput.html),
[SCF algorithms](https://github.com/psi4/psi4/blob/master/doc/sphinxman/source/scf.rst).
Installed `psi4/driver/procrouting/scf_proc/scf_iterator.py` separately
budgets JK and DFT collocation arrays from the configured memory.

Candidate for a coordinated validation slot: lower **algorithm allocation**
to 768 MiB or 1 GiB and disk-backed DF, keeping basis, functional, active
space, scan points and convergence tolerances unchanged. This is not yet
validated or applied. Disk-backed SCF does not bound CASSCF integral/CI
memory, and lowering allocation is not proof of safe RSS. Measure whole
worker peak private/RSS and compare converged energies, spin and NOON
against the default configuration before proposing any revised safety gate.
Changing to a smaller basis or replacing CASSCF is not an equivalent
memory optimization. Main was asked to coordinate this validation; the
existing parent remains on its original 2.5 GiB gate in the meantime.

P6 cheap refit of **this campaign's** stored calibration scan explains the
apparently large N-C RMS: its final descending point was intentionally
trimmed by the unchanged fitter, while the reported 15.565 kcal/mol RMS
included that point. Retained-point RMSE is 0.577 kcal/mol (7 retained,
1 trimmed); C-C and forming N-C RMSE are 1.144 and 2.182 respectively.
Future source now labels full-scan RMS and separately reports retained
RMSE, soft-L1 score and trimmed count. The robust acceptance gate and
optimizer are unchanged. No completed trajectory, calibration cache or
historical result was rewritten. Regression suite: **22 tests passed**.

This correction does not resolve sampling: the reactive model retains
reactant valence angles/torsions and fixed nonbonded exclusions while only
two bonds become dissociative Morse terms and a third Morse pair forms.
Its two biased CVs do not directly bias all three changing distances.
These are model/CV adequacy questions, not evidence of an implementation
failure or a reason to remove forces to induce a reaction. Validating a
different reactive Hamiltonian/CV set needs a separately identified protocol
and fresh runs; the completed default still has no sampled product basin
and no defensible reaction barrier.

## Latest status (18:18)

Phase 6 **full default computation completed**, via the authorized same-campaign
continuation: 200,000 explicit steps (0.4 ns), 500,000 implicit steps (1 ns),
all 400/500 hills preserved. The continuation exited **0** at 18:05:08;
wall time was 12,270.98266029358 seconds, including a roughly three-hour log
gap after 14:59:49. It is not an active-compute benchmark. Parent termination
exit 4294967295 remains recorded without alteration. The attempt-root
`continuation_status.json` now records completion and its scientific limits.

**Scientific result remains incomplete:** neither leg sampled a product basin;
both activation barriers and the solvent barrier shift are null. The sampled
FES limits are 1.84 A explicit and 1.80 A implicit. Execution completion is
not barrier convergence; no gate was relaxed and no 10 ns claim is made.
No rerun is needed merely to recover trajectory coordinates, which are intact.

Both original DCDs were preserved. Metadata-only copies
`traj_explicit_timing_corrected.dcd` and `traj_implicit_timing_corrected.dcd`
have sidecar `.provenance.json` files in the same results directory. Only
three timing integers changed: first step 5000, stride 5000, final step
200000/500000. Times are **production-relative**, excluding equilibration:
10..400 ps (40 explicit frames) and 10..1000 ps (100 implicit frames).
The DCD header/record reader agrees with the CSV times within float32
timestamp precision. Independent MDTraj reads confirm exactly identical
coordinates and box arrays; MDTraj's low-level DCD API does not expose
timestamps, which were checked independently against the header and CSV.
Final distance CV errors vs the last deposited hill were 4.39e-8 A explicit
and 9.54e-9 A implicit. All coordinates and hills are finite. Mean recorded
temperatures were 299.71 K explicit and 290.95 K implicit (19-atom solute).

| Artifact | SHA256 |
| --- | --- |
| original explicit DCD | dcb66ad507850b12e0d7d4b5b3177e6cd0b6b9df327ba42a4cd27829c081b567 |
| corrected explicit DCD | fd5059dc8b8acef5d62eeaf23c414fd86908f9f3beb141fd2b375c5f91e799f6 |
| original implicit DCD | 608ee06e2e201e7e4d3841d9c7f21df475fe53a76ddc08d661aa4cdb1fa834d7 |
| corrected implicit DCD | f428ddeec8a6ee9b215ae4fcf9fbea61da205e43ec2af0f2528f6840207643c3 |

Latest regression suite: **21 tests passed**, including metadata-only DCD
repair and an actual OpenMM DCD-writer round trip. Future writer now emits
the correct first step/stride; this was not retroactively applied to the
already-running continuation, whose output was repaired as a separate copy.

Phase 7 default attempt `phase07/20260910T181009` started after checking
no duplicate Phase 7/8 processes or attempt directories existed. It uses
corrected Phase 4 inputs and two threads; each Psi4 worker remains gated
on at least 2.5 GiB available RAM. Phase 8 is queued sequentially.
Phase 7 completed all eleven constrained geometries with `converged=true`
and all eleven xTB/MACE-OFF/ANI-2x energy comparisons. It is held before
7A point 0 because available RAM is approximately 1 GiB. No Psi4 worker has
launched. Orchestrator PID 38388 remains alive and checks every 30 seconds;
the earlier `orch spawn` log is intent preceding the gate, not evidence
that a QC process started. QC results and the full Phase 7/8 runs remain
incomplete. Main received the priority/queued-launch notice before this gate.

Scope: only the three phase drivers, optional Phase 8 renderer, this note,
and `test_phase6_8_rerun.py`. No staging, commits, pushes, hardware operations,
or historical output reuse. Read `PHASE6_10_AUDIT_FIXES.md` before this audit.

## Fixes

- Phase 6: identify all missing Phase 4 inputs before QM work; map reactive
  reactant indices into product ordering for the product-side Morse scan;
  respect OPENMM_CPU_THREADS/OMP_NUM_THREADS instead of using every CPU;
  correct the introductory free-energy formula to match the audited code.
  Preserve small WT hills, reject saturated bias capacity and reject step
  counts that would exceed capacity or truncate the requested duration.
- Phase 7: impose the requested constrained bond length during relaxation;
  give child interpreters their own Library/bin and root PATH and remove
  inherited PYTHONPATH; forward force to geometry worker; reject missing
  prerequisites, incomplete/nonfinite QC series and failed workers rather
  than returning unconditional success. Figure failure returns nonzero.
  Subtract S(S+1)=2 before applying triplet spin-contamination thresholds;
  missing scan points are not treated as zero contamination. Full completion
  also requires converged geometries and complete MACE/ANI comparisons.
- Phase 8: forward the thread budget into every Psi4 worker, configure its
  DLL path, remove inherited PYTHONPATH, and preserve stdout/stderr logs.
  Require the existing gap AND gradient gates for MECI convergence. Failed
  finite differences cannot become zero gradients or certify convergence.
  Preserve matched coordinates/energies when perturbation evaluations fail.
  Correct inverse-BFGS update and retain the pre-step coordinate for secants.
  Keep the 2e-7 energy tolerance on all CASSCF retries (more iterations only).
  Record module/figure failures and return nonzero on failed default runs.
  Remove unsupported claims of numerical equivalence between Psi4 and PySCF.

## Cheap verification and dependency findings

`PRIMARY_PYTHON -m unittest -v test_phase6_8_rerun test_phase6_10_audit`
Final result: **16 tests passed** (9 rerun regressions plus 7 prior audit
regressions), in 1.772 seconds. Scoped `git diff --check` passed. These are
focused regressions, not full reruns; they include inverse-BFGS secants,
MECI gates, triplet contamination, small-hill preservation and saturation.

With PATH prefixed by `phase7/Library/bin;phase7`, imports report Psi4 1.11,
NumPy 2.5.2 and SciPy 1.18.0. The `qbscf` Python exists but has neither NumPy
nor PySCF; it is not a usable substitute. No packages were installed.

At the compute-slot request, free physical memory was approximately 450 MB.
No heavy computations have been launched without the main agent's allocation.

## Commands and prerequisites

Run from this checkout using
`C:/Users/HUIWEI/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe`:

```text
python rerun_phase1_19_campaign.py --phases 8
python rerun_phase1_19_campaign.py --phases 6 --phase4-attempt ../phase1-19-rerun-20260910/phase04/20260910T134845
python rerun_phase1_19_campaign.py --phases 7 --phase4-attempt ../phase1-19-rerun-20260910/phase04/20260910T134845
```

The campaign runner chooses molecular Python for these phases, sets
OMP/MKL/OPENBLAS threads to 2 and configures the molecular DLL path. Phase 8
forwards two threads to QC workers. Default steps, trajectories, scan grids
and convergence gates are not reduced. Use sequential slots because Psi4
requests 2 GB and the machine has 16 GB total RAM.

Source-only snapshots cannot run Phase 6/7 computations without fresh
Phase 4 artifacts. Phase 6 needs `results_phase4/reactant_3d.mol`,
`product_3d.mol`, and `images_idpp.xyz`; Phase 7 needs `reactant_3d.mol`.
The runner now accepts `--phase4-attempt` and copies only those inputs to
the same `results_phase4/` relative paths inside the fresh snapshot, recording
source paths and SHA256. Phase 8 needs no prior-phase artifact. Do not copy
historical output directories or QM caches.

Corrected source: `../phase1-19-rerun-20260910/phase04/20260910T134845`.
All three hashes were independently verified:

| File | SHA256 |
| --- | --- |
| reactant_3d.mol | E52CD7FBC8028BFBDF97841B04C1C87EC83B56C75097A894EC82E046E0DB8C78 |
| product_3d.mol | 7F3952E8F11112FC4CF3537FA66B04891A9D74E5B974C0B657EC74C5D830B0CD |
| images_idpp.xyz | 326423CE9E9159B487DA22863B2321AAF7FB79F3CB371CCF47FDCBE7CDAD9C48 |

These are stage-1 geometries, not TS validation. The superseded
`20260910T134734` attempt is excluded because MOL and IDPP coordinates were
inconsistent. It was read for a preliminary size estimate, subsequently
discarded, but never copied into or used for a Phase 6-8 rerun.

## Prospective Phase 6 CPU-only 10 ns assessment

Code defaults are 200,000 explicit steps (0.4 ns) and 500,000 implicit steps
(1 ns), with dt=2 fs. Default full completion is not a 10 ns run.
Ten ns for one leg requires 5,000,000 production steps. At the unchanged
cadence this requires 10,000 hills for explicit or 5,000 for implicit.
The present 512-slot force supports at most 0.512 ns explicit or 1.024 ns
implicit deposition. Previously it silently stopped adding bias beyond the
capacity and skipped small hills. It now rejects this unsupported protocol;
neither changing step counts alone nor treating frozen-bias MD as WTMetaD
is valid. A validated scalable force representation is required for 10 ns.

There is no freshly measured throughput. Conditional arithmetic for the
explicit 5-million-step leg only is:

| Sustained steps/s | Production hours |
| ---: | ---: |
| 50 | 27.8 |
| 100 | 13.9 |
| 250 | 5.6 |
| 500 | 2.8 |

These are scenarios, not a forecast or performance measurement. Add QM
charge/Morse calibration, force-field generation/compilation, equilibration,
the implicit leg and analysis. The implicit leg currently uses OpenMM's
Reference platform and needs its own benchmark. Ten ns on both legs totals
10 million steps and different throughputs for the two systems.

The FES broadcasts a 86x73 grid across hills. One full float64 intermediate
is 24.5 MiB at 512 hills, 239.5 MiB at 5,000 hills or 479 MiB at 10,000 hills.
Several intermediates may coexist; the 10,000-hill expression alone can
therefore need more than 1 GiB during analysis, before driver/QM/OpenMM/FF
allocations. Compiling thousands of analytic Gaussian terms has unmeasured
time and peak memory. Current free RAM below 1 GB is insufficient headroom
to assume safety. Exact RSS/throughput must be measured after slot allocation;
no whole-process upper memory bound can be justified from source inspection.

## Full rerun status

One compute slot allocated for sequential default Phase 6, then 7, then 8.
Before every Psi4 worker, a Windows available-physical-memory check requires
2.5 GiB; otherwise the launcher logs HOLD and rechecks after 30 seconds.
No shutdown or hardware action is part of this campaign.

- Phase 6 `20260910T135745`: full-default attempt failed after 22.3 seconds,
  before MD, because an installed CUDA plugin was selected on a CPU-only
  host (CUDA Context error 801). Also found Windows GBK decoding of xTB
  UTF-8 output causing scan fallbacks. Both issues fixed in source: probe
  actual Context usability for auto-selection; decode subprocess output
  explicitly. The latter fix also applies to Phase 7/8 subprocesses.
- Phase 6 `20260910T135824`: fresh retry using the command above plus
  `--platform CPU`; scientific defaults unchanged, 200k explicit / 500k
  implicit steps. All three Morse scans fitted; FF export delta energy
  2.69e-6 kJ/mol; 1497 atoms, 492 waters, one Na+ and one Cl-. CPU Context
  creation took 206.2 seconds, transient peak RSS about 1.73 GiB, then
  private memory fell to about 0.70 GiB. Equilibration ended at 294.9 K.
  At 14:06:31 explicit production reached 2500/200000 steps and all five
  scheduled hills, measured 112 steps/s (early explicit ETA 29.3 minutes).
  This is progress, not full-run completion or scientific acceptance.
- Phases 7/8: not yet launched; awaiting Phase 6 completion and RAM gate.

Focused suite rerun after launcher fixes: 16 tests passed, scoped diff check
passed. Actual final outputs/gates remain to be reviewed after completion.
Any future 10 ns Phase 6 run still requires scalable-bias work; this campaign
does not claim 10 ns.

### Lifecycle and authorized continuation preparation

The running `135824` snapshot retains its explicit Context while building
implicit. At the user's request, **future source only** now releases the
explicit Simulation/Context/Integrator references before constructing the
implicit system. It also adds `--resume-implicit`, rejecting incomplete
explicit step/hill counts and fingerprinting required completed artifacts
(results, state, checkpoint, trajectory, PDB, calibration files, original
snapshot source, continuation source, campaign status). It writes
`implicit_resume_provenance.json` and does not construct or recompute the
explicit leg. A partial implicit DCD continuation now uses OpenMM's append
mode with a readable/writable existing file, avoiding an appended second header.

If the running job becomes unsafe at transition, stop only tracked PID 32736
after verifying its campaign path. Preserve its completed artifacts. Launch
the updated source by absolute path under molecular Python with the **same
campaign snapshot as cwd**, `--platform CPU --resume-implicit`, and the
original two-thread/DLL environment. This is a documented same-campaign
continuation, not an independently fresh run and not historical cache reuse.
Do not stop/restart merely to apply this change while the current job is safe.

Additional Phase 7 correction: `<S^2>` needs `M_s*(M_s+1)` rather than `M_s^2`;
the old expression returned 1 instead of 2 for a pure triplet. Pure singlet,
doublet and triplet regressions pass. Latest combined focused suite: **19
tests passed** (including the memory gate and continuation provenance tests).
