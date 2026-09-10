# Phase 17 SCF convergence audit — 2026-09-10

## Active first run is untouched

Read-only audit of `../phase1-19-rerun-20260910/phase17/20260910T134422`.
No heavy computation, process termination, or snapshot edits. Preserve main's
two raw-fstring warning fixes in the shared Phase 17 source.

At the last inspection Am3+ iteration 15 had dE = 1.16e-3 Eh, versus the
3e-6 Eh criterion (iteration 10: 1.85e-2 Eh). This is encouraging decrease,
not proof that iteration 26 will pass. The running old source stops on energy
alone or exhausts 26 iterations; it never marks SCF failure and does not check
the returned Dirac orbital `ok` flag. Its final solve is not tested for
fixed-point consistency. Logs print only selected energies, not density errors.

## Corrected source

- Increase the default budget from 26 to 80 iterations. Preserve 3e-6 Eh
  energy tolerance and existing 0.5/0.25 density mixing schedule. The observed
  decreasing energy error justifies more iterations; no evidence yet justifies
  replacing mixing with DIIS or loosening a tolerance.
- Require unmixed integrated absolute density residual per electron below
  3e-6 as well as energy convergence. Testing the unmixed update prevents a
  small mixing factor from falsely making a residual appear converged.
- Check the final solve against its actual input density and previous energy.
  A failed candidate final check continues within the remaining budget;
  exhaustion raises `SCFConvergenceError`;
  dependent Breit/multiplet/bonding calculations do not continue.
- Save per-iteration energy and density residual history, tolerances, iteration
  count, final residuals and `converged` in each successful atom's JSON summary.
  An SCF convergence exception writes `results_phase17/phase17_scf_failure.json`
  before exiting nonzero. An orbital shooting failure stops with an explicit
  exception; no nonconverged orbital is accepted.
- Stage reload rejects atoms lacking positive SCF convergence evidence; legacy
  snapshots must not silently become certified inputs to corrected stages.
- Initial dE is unknown (`null` in history), not a reported zero.

## Cheap verification

```powershell
$env:PYTHONPATH='../phase24-deps;../phasefix-test-deps'
$env:OMP_NUM_THREADS='2'
$env:MKL_NUM_THREADS='2'
$env:OPENBLAS_NUM_THREADS='2'
& C:/Users/HUIWEI/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe -m unittest -v test_phase17_rerun test_phase14_17_audit
```

12 tests passed in 1.183 seconds. New tests use synthetic orbital solves to
exercise the actual SCF loop: reject a density cycle with constant energy,
accept a fixed point, reject a failing final solve, reject a false orbital
convergence flag, and reject an invalid budget. This is control-flow and
numerical-accounting verification, not a physical heavy-atom rerun.

## What can be salvaged from the first run?

Keep its original log, source hashes and artifacts. Independent hydrogenic and
spherical-spinor validation are usable as evidence of those particular checks:
the observed log reports spinor deviation 9.77e-15 and hydrogenic tests passing.
Hydrogenic Coulomb/Breit machinery tests, if subsequently produced, are likewise
independent of atomic SCF; retain their individual measured errors.

Atomic energies, relativistic/NR comparisons, atomic Breit corrections,
j-splittings/multiplets and bonding/QTAIM/EDA inherit the missing SCF evidence.
Their files may be retained as provisional numerical outputs, but neither
process success nor reaching the cap establishes convergence. Even a small
last logged dE cannot establish density consistency or successful orbital solves.

If the first run ultimately writes JSON and radial NPZ, a separate *allocated*
post-hoc fixed-point solve for each atom could establish residuals from those
orbitals. That is additional computation and has not been done. The old files
do not contain enough input-potential/residual history to certify it read-only.
If that audit fails, the SCF and every dependent result must be recomputed.
For the requested fresh campaign, use a corrected source-only rerun after
allocation; do not relabel the first run as converged or post-fix validation.

## Allocated rerun command

```powershell
& C:/Users/HUIWEI/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe rerun_phase1_19_campaign.py --phases 17 --profile primary
```

Still on compute hold. The higher budget and stricter checks are not yet
validated on Am/Eu/N; failures may expose additional numerical issues. No
staging, commits or pushes performed. Original Phase 10–13 work remains ready
and held as documented in `PHASE10_13_RERUN_NOTES.md`.

## Corrected attempt 20260910T140133: failure and quadrature repair

Main launched this full attempt; this worker only inspected it. At Am iteration
25, dE=4.7097e-5 and density residual=2.8430e-5 confirmed the old 26-iteration cap
would be insufficient on the reproduced trajectory. At iteration 35, the loop
passed both gates, but final density residual 7.0484e-6 failed (energy change
2.3658e-6 passed). Failure JSON records this honestly.

Source inspection identified inconsistent normalization: the iteration globally
rescaled density, the final solve did not; Dirac shooting uses midpoint small-
component quadrature while density uses node quadrature. NR eigensolver vectors
are Euclidean-normalized and likewise require conversion to radial quadrature.
The warm-start condition only resets for `scf_it < 3` or implausible energy;
the previous final value 99 did not explain this discrepancy.

Revision normalizes every P/Q pair with the actual density node weights and
uses `scf_density_from_orbitals` for both loop and final construction, without
global density rescaling. Charge count and absolute charge error are recorded;
an unnormalized density fails a separate relative charge check of 1e-10.
Failed final candidates continue within the 80-iteration budget, with a saved
verification history; all convergence tolerances remain unchanged.

Final cheap check: `test_phase17_rerun test_phase14_17_audit`, **16 passed in
1.514 seconds**, scoped diff check clean. Added NR/Dirac quadrature unit norms,
input-amplitude invariance, rejection of unnormalized charge, and candidate
continuation tests. Main is responsible for launching the next fresh full17.
No extra Phase17 computation was launched by this worker. The old attempt has
diagnostic JSON only, not density arrays; it cannot be restarted from that JSON.
Phase10–13 have since been allocated and are running serially (see their note).

## Third attempt 20260910T141424: saved atomic convergence verified

Main launched this fresh full attempt. Read-only JSON inspection after the
atomic stage confirmed all four `scf.converged` flags true, with unchanged
3e-6 energy/density thresholds:

| Atom | Iterations | Final dE (Eh) | Final density residual | Charge error (electrons) |
| --- | ---: | ---: | ---: | ---: |
| Am3+ | 35 | 1.923e-6 | 1.738e-6 | 1.42e-14 |
| Eu3+ | 34 | 3.090e-7 | 2.159e-6 | 7.11e-15 |
| Am3+ NR | 34 | 2.646e-8 | 2.009e-6 | 1.42e-14 |
| Eu3+ NR | 33 | 7.363e-8 | 2.201e-6 | 0 |

The loop/final normalization discrepancy is resolved in this actual calculation.
Breit started at 14:43:19; full Phase17 completion and downstream acceptance
remain main's responsibility and are not asserted by this atomic-stage check.

## Breit optimization and convention audit

Source now builds angular pairs in chunks and contracts only nonzero Dirac
matrix connections. It no longer allocates per-point dense16x16 T/K/bilinear
arrays or the entire Kronecker-product pair wavefunction. Full-grid defaults
remain30x36 angular points and40 radius ratios, with progress every5 ratios.

Same-grid benchmark against the running third-attempt source used8x10 angular
points,3 radius ratios,chunk2000 (cheap equivalence test, not full computation):

| Kernel | Old seconds | New seconds | Speedup | Old/new peak MiB | Max difference |
| --- | ---: | ---: | ---: | ---: | ---: |
| Coulomb | .175595 | .030233 | 5.81x | 33.42 / 2.19 | 1.11e-16 |
| Breit | .737325 | .117920 | 6.25x | 41.24 / 8.82 | 2.78e-17 |
| Gaunt | .168245 | .060709 | 2.77x | 33.41 / 6.60 | 9.44e-16 against twice old |

Conventional Gaunt is -alpha1.alpha2/r12, as explicitly defined in the primary
[DIRAC Hamiltonian documentation](https://www.diracprogram.org/doc/release-21/manual/hamiltonian.html).
The old branch returned half Gaunt while labeling it Gaunt. The full Breit
operator retains its one-half factor; only the reported Gaunt/gauge split changes.

Also corrected radial table interpolation from sqrt(r_less/r_greater) to
r_less/r_greater, matching the actual table coordinate. Radial integration
groups identical P²/PQ/Q² factors and chunks128 rows, avoiding full NxN arrays.
Complex sums are retained instead of discarding their imaginary component.
Dense synthetic radial contraction and exact-ratio tests verify the change.
Final test command `test_phase17_rerun test_phase14_17_audit`:19 passed2.173s.

Main was informed these changes justify stopping its old Breit calculation
and resuming Breit/multiplet/bonding/figures from only the fresh validated atomic
JSON/radialNPZ of20260910T141424, with provenance. This worker neither stopped
nor launched Phase17. The small-grid speedup does not guarantee full-grid time.

## Hydrogenic radial resolution follow-up

Fixed40a0 grids underresolve high-Z1s orbitals. Source now uses uniform Z*r
grids (4000 Coulomb nodes,3000 hydrogenic Breit nodes), preserving domain and
resolution relative to the orbital extent. This changes only the hydrogenic
diagnostics; the four actual atomic wavefunctions remain unchanged.

A cheap exact discrete Coulomb-kernel probe at4000nodes gives relative radial
errors on old fixed grids: Z1=6.669e-6,Z20=.0025494,Z50=.0115279,Z80=.0084350.
Scaled grids give6.669e-6 for every Z, verifying scale consistency. This isolates
radial error; it does not remove the angular table's remaining quadrature error.
The new tests show grid-doubling improvement and identical J/Z for Z1and80.
Combined Phase17/audit test command:20 passed2.737s.

Future runs save matching angular tables in the same result directory, with
signature including grid sizes, kappa,kernel and full-Gaunt convention version.
This permits further radial checks without recomputing angular tables. Existing
running snapshot180947 lacks this cache output. Its scaled Coulomb diagnostic
can be derived exactly as Z times the existing Z1 value, with the same.449%
table error; hydrogenic Breit gamma varies with Z and needs its angular table
for corrected radial integration. Main was advised to complete current stages,
keep those hydrogenic rows provisional, then allocate a targeted check if needed.

## Targeted post-audit helper (prepared, not launched by this worker)

`audit_phase17_scaled_hydrogenic.py` requires completed zero-exit Breit,
multiplet and figures stages, and an exact SHA256 for the NEW source. It checks
the old source against the attempt manifest and verifies that the anchor's
angular and radial contraction implementations are AST-identical to new source.
This comparison passed against actual attempt20260910T180947.

The helper creates a unique hydrogenic_postaudit_TIMESTAMP directory, preserves
read-only content-addressed original JSON and the exact new source, and writes
audit.json with source/helper/input/cache/output hashes. The canonical original
phase17_results.json is never overwritten. Corrected data are written separately
as phase17_results.corrected.json. Only the Breit angular table is recomputed
(full30x36x40 defaults) into a fresh cache; no atomic or Gaunt work is invoked.
Coulomb rows derive from measured Z1 by discrete J(Z)=Z*J(1), with explicit
linearity provenance and the unchanged measured relative error. They are not
independent high-Z certificates.

Hydrogenic Breit uses3000 scaled radial nodes, plus6000-node refinement checks.
The declared new radial-only gate is relative change<=0.001; this does not
establish angular convergence or global scientific acceptance. Nonfinite values,
significant imaginary energies or failed radial gates produce a failed audit
and no published figure. Candidate corrected JSON is retained if available.
On success only fig4 is regenerated and published, preserving the old figure
read-only with hash. Figure4 now labels derived Coulomb rows honestly.

Run from this checkout AFTER main completes (main launches, not this worker):

```powershell
$env:PYTHONPATH='../phase24-deps;../phasefix-test-deps'
$env:OMP_NUM_THREADS='2'; $env:MKL_NUM_THREADS='2'; $env:OPENBLAS_NUM_THREADS='2'
& 'C:/Users/HUIWEI/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -u audit_phase17_scaled_hydrogenic.py --attempt ../phase1-19-rerun-20260910/phase17/20260910T180947 --source run_phase17_relativistic_actinide_quantum.py --expected-source-sha256 4092d03e56f4302d5adc8745920ec5c09833a317369ab361854f3049a61a6324
```

Prepared helper SHA256:51ea7cb3bce3356a13080a4c876d765ea87203f720a2aec4cd44b405af39d537.
Cheap verification: `python -m unittest -v test_phase17_rerun test_phase14_17_audit`
passed24tests in2.453s; scoped diff whitespace check passed. No full helper or
full-grid post-audit numerical result is claimed yet.

## Figure3 / promolecular derivative post-audit

Saved180947 planes confirm the extreme negative values occur at index(120,1):
actual x=0,z=.013333a0, not the displayed z~1.6. Matplotlib received[x,z]
instead of[z,x], transposing the physical positions. Original minima were
-6.4241e8(Am) and-3.9761e8(Eu). The core is unresolved on this display grid;
these large numbers are not accepted as converged nuclear derivatives.

Corrected source transposes plotting arrays and computes the full axisymmetric
3D Laplacian rho_xx+rho_zz+rho_z/z, with axis limit rho_xx+2rho_zz.
Axis descriptors now consistently evaluate z=0, with even-reflection transverse
derivatives. Undefined outer derivative boundaries are NaN, explicitly masked
and counted rather than falsely zero. Figure uses explicit diverging levels,
discloses saturation and a3-cell nuclear-core exclusion. This does not certify
linear radial interpolation derivatives or remove their refinement uncertainty.
Labels are neutral: atomic-promolecule axis minima and WH mixing proxies,
not self-consistent complex QTAIM, signed overlap, molecular EDA or a binary
Am-bonds/Eu-does-not conclusion. The modeled plane has one M and one N, not N-M-N.

Cheap cached-density probe (no SCF): Am lap=.1616950157,H=.0321659817;
Eu lap=.1667846789,H=.0336393847. Coarsened-grid laps=.1653098831 and
.1586523071; these differences are diagnostics, not convergence passes.
Axis gradients=-.0019573142 and-.0002279955; minima are not certified stationary
points. Nitrogen SCF metadata was not persisted with these plane densities.

`audit_phase17_qtaim_figure3.py` reuses these exact cached fresh M+N planes,
so no N/Am/Eu/NR SCF, FFT, overlap, Breit or Gaunt recomputation is needed.
It requires a completed, hash-verified Breit-corrected candidate, rejects changes
outside the explicitly allowed Breit fields, preserves read-only original JSON,
candidate, source, planes and priorfig3, then writes a separate merged corrected
JSON with audit provenance. Only bonding.qtaim.Am/Eu and bonding.model_scope
change relative to that Breit candidate. Canonical JSON remains untouched.
Only figure3 is published; original planes remain untouched too.

Main launches after Breit post-audit completion, with primaryPython/deps as above:

```powershell
python -u audit_phase17_qtaim_figure3.py --attempt ../phase1-19-rerun-20260910/phase17/20260910T180947 --source run_phase17_relativistic_actinide_quantum.py --expected-source-sha256 f5cf4aed9bfadde27508dd6490c5ebc48f09a66b995a366dc342e6d26164603d --breit-candidate PATH_TO_COMPLETED_BREIT_AUDIT/phase17_results.corrected.json
```

Helper SHA256:6b45d7c0b4fd374cc781939f5f22bed3c5390c685157b33234bb4619d9c7d85f.
26 focused tests passed2.436s, including exact polynomial3D/axis Laplacian,
on-axis minimum consistency and rejection of unrelated merge changes.
This worker has not launched the figure3 helper; rendering/output review remains.

## Completed post-audits and formal acceptance at18:42

Main completed hydrogenic102402808366Z, QTAIM103025628566Z and KB104108041864Z
post-audits under180947, all exit0. Formal review is180947/acceptance_review.json:
status `specified_checks_passed_with_limitations`, whole_model_validated=false.
Authoritative reviewed candidate is kb_postaudit_20260910T104108041864Z/
phase17_results.merged_corrected.json, SHA256
2f4888159cc8cbaf69373e517d593867557f634cb6ca8910828901450547511b.
Canonical original JSON remains unchanged, SHA256
d0f731c64f66ca4f9d6858453a0ef5e97e4069b29f6845fb2ec375895f53f963.

All4 persisted atomic SCFs meet unchanged energy/density3e-6 tolerances.
All1844 numeric finalJSON values and204 atomicNPZ arrays finite; corrected planes
have721 deliberately undefined outer derivative cells perion, not interiorNaNs.
Ladder formula residual<=1.46e-11cm^-1, centered weighted trace<=2.70e-10cm^-1.
These are internalLS algebra checks; J1 experimental relative errors19.3%Am and
387%Eu disallow quantitative experimental validation. N log reaches iteration32
density2.69e-6,dE2.78e-9, then strict source returns to bonding, but finalN
verification/charge/orbitals were not persisted: process evidence only.

Original GaussianKB certificate403.599 relative error failed. Corrected separate
LL/SS blocks, P=r*exp(-a*r²), analytic kappa=-1 derivative, RKBsmallbasis and
positive-total-energy branch(-c²,0). Defaultbasis increased10to32 on measured
10/16/24/32 comparison. New32basis E=-4860.633745613Eh vs-4861.197904370,
relativeerror1.16053443569e-4;24to32change1.16711962428e-5; eigenresidual5.39e-13.
Explicit finite-basis gates1e-3/1e-3/1e-8 pass, not a global stability proof.
Newhelper audit_phase17_kb.py preserves upstreamJSON/hash/source/priorfig4,
changes only certificates KBfields, and regeneratesfig4 only. SourceSHA
d828beb1875239adc5f15e25dd9c92d535fe35bcb20ae9e92fd49b6abe5d43d1.
27 focused tests passed1.930s. Finalfig4SHA
e203a3a873e14d8d26b88f287afe6f5b64103542f1dee053daada4f0dae488da.

Remaining limits: angular convergence not certified; scaled Coulomb retains
.449268% measured angular/radialtable error; promolecular derivatives not
refinement/stationary-point certified; WH/density-overlap proxies not complex
QTAIM/EDA. Main performs final visual review. No atom redo or canonical overwrite.
