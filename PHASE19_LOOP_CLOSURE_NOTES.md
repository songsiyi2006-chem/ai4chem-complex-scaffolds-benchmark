# Phase19 bounded peptide closure sidecar

## Delivered scope

Replaced the tangent spline in `realize_backbone` with the production
`close_peptide_loop` solver. Changed only
`run_phase19_active_inference_denovo_enzyme.py`, this note, and the new
`test_phase19_loop_closure.py`. Read chem-ai4s SKILL.md, molecules.md,
quantum.md and PHASE19_CLASH_NOTES.md. Used the existing phase2ff runtime;
no installation, training, MD, full geometry probe, snapshot modification,
commit, or work on other phases. Tests compile selected production AST
functions instead of importing the entrypoint; no torch import.

## Construction and rejection contract

- Optimize only new-residue phi/psi in degrees with scipy least_squares,
  bounded to [-180, 180]. Coordinates are angstroms. Reuse production
  `step_forward`/`place_atom` with ideal lengths/angles and trans omega.
- The previous residue's entire N/CA/C/O frame is fixed, including its
  carbonyl orientation. Predict both next N and next CA after the last new
  residue and fit them to the actual fixed next rod. Neither endpoint nor
  any motif/rod coordinate is changed or snapped to a numerical solution.
- Accept only when each endpoint error is at most 1e-6 A, every checked
  nonbonded separation is at least 2.65 A, and the existing strict
  `require_constructed_backbone` accepts the actual junction coordinates.
  The complete assembled backbone still passes that unchanged guard before
  downstream processing. All pre-existing scientific thresholds remain intact.
- Collision residuals include every new N/CA/C/O against local nonbonded
  backbone pairs and every atom of all already accepted and future fixed
  rods. The local graph excludes only distances one and two in the bond
  graph. The two endpoint residues are supplied separately, not duplicated
  as external obstacles. A 2.67 A optimization buffer does not change the
  2.65 A acceptance cutoff. Exact coincident nonbonded atoms are penalized.
- Try lengths 2 through 8, with two fixed torsion starts per length.
  Each optimizer invocation has max_nfev=80. A shared default budget of
  2400 residual evaluations per junction also counts finite-difference calls;
  each start has a bounded share. Exhaustion produces BACKBONE_REJECT.
  The API permits a budget up to 10000, but production uses 2400.
- Reject gaps beyond their covalent contour upper bound without optimization.
  Reserve two residues for every remaining junction when applying the total
  200-residue cap. There is no spline fallback or Cartesian geometry repair.

This is a deterministic local feasibility search, not an exhaustive closure
algorithm or an energy calculation. A budget failure does not prove that no
solution exists. Starts do not cover all torsional conformations; torsions are
not constrained to Ramachandran regions during optimization. Existing
downstream torsion checks retain responsibility for their original criteria.
Fixed-rod/fixed-rod collisions cannot be repaired by this loop solver and are
still rejected by the full constructor guard. Greedy successive closures do
not backtrack earlier loops. Full-candidate yield/runtime and physical fold
quality have not been established.

## Verification

Every command used one thread for OMP, OpenBLAS, MKL and NumExpr, disabled
bytecode writes, and wrapped unittest in subprocess.run(timeout=30).

```powershell
$env:OMP_NUM_THREADS='1'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:NUMEXPR_NUM_THREADS='1'
$env:PYTHONDONTWRITEBYTECODE='1'
$env:PATH='C:/Users/HUIWEI/miniconda3/envs/phase2ff/Library/bin;'+$env:PATH
& 'C:/Users/HUIWEI/miniconda3/envs/phase2ff/python.exe' -c 'import subprocess,sys; subprocess.run([sys.executable,"-m","unittest","-v","test_phase19_loop_closure"],timeout=30,check=True)'
```

Final new suite: **8/8 passed, 1.261 seconds** (2.115 s shell wall time).
Coverage includes exact deterministic recovery without fixed-atom movement,
noncanonical endpoint optimization, sensitivity to both endpoint atoms,
future-rod exact-overlap rejection, impossible distance, evaluation exhaustion,
invalid inputs, and length bounds. An actual-source eight-rod realization
fixture closes seven two-residue gaps into 62 residues with zero backbone
clashes and exact motif pins. This is a synthetic helix with omitted sections,
not a full candidate/probe. A separate actual assembly-block fixture verifies
the precise 200-residue limit and reservation of later gaps; that count-only
fixture intentionally stubs the geometry solver.

Initial regression command with `test_phase19_orientation test_phase19_clashes`
(historical result, superseded by the migration below):
**12/14 passed, two obsolete fixture errors, 2.384 seconds**. All four
orientation tests and eight packing/clash/strict-constructor tests passed.
That initial aggregate command was not green. The two errors were:

1. `test_actual_loop_spline_ignores_next_n_and_can_be_collinear`: StopIteration,
   because it searches the source for the removed `build_loop_spline`.
2. `test_realization_rejects_actual_spline_before_downstream_gates`: NameError,
   because its AST namespace does not include the new `close_peptide_loop`.
   Its intended spline-specific fixture must also be replaced, not merely
   taught to import the new solver.

Those tests initially lived outside this sidecar's explicit ownership and
were left unchanged until the user authorized the migration below.
The solver is wired into production and tested
at the assembly boundary; complete enzyme geometry repair is still unproven
without an independently authorized full-candidate evaluation. The historical
661-clash snapshot and all archived results remain unchanged.

## Authorized test migration and freeze

The user subsequently authorized changing only the two obsolete tests in
`test_phase19_clashes.py`. Both were migrated without skips or weakened gates:

- `test_actual_loop_closure_depends_on_next_n_and_is_noncollinear` runs the
  actual solver, verifies both predicted endpoint coordinates and 111.2-degree
  N-CA-C angles, and passes the unchanged strict constructor guard. Observing
  real optimizer residuals proves that changing next N changes its endpoint
  residual while next CA remains fixed; the incompatible endpoint is rejected.
- `test_realization_and_probe_propagate_closure_failure_before_downstream`
  injects a specific BACKBONE_REJECT exception at closure, verifies that actual
  realization and the AST-only probe caller propagate that same exception,
  records zero downstream calls, and preserves fixed rods and motif targets.

Combined command, using the single-thread environment above:

```powershell
& 'C:/Users/HUIWEI/miniconda3/envs/phase2ff/python.exe' -c 'import subprocess,sys; subprocess.run([sys.executable,"-m","unittest","-v","test_phase19_loop_closure","test_phase19_clashes","test_phase19_orientation"],timeout=30,check=True)'
```

**22/22 passed: 8 closure + 10 clashes + 4 orientation, 3.049 seconds**
(4.301 s shell wall time). No tests skipped. The 30-second subprocess timeout
was enforced. `git diff --check -- test_phase19_clashes.py` also passed.
This was small AST-fixture verification, not a full geometry probe or enzyme
validation. The migration changed no production code or other files and made
no commit.

**Frozen for the user's final unified regression.** The test run is complete.
This documentation update records its result; no Python source is changed in
the freeze turn. This worker will make no further Python edits during the
unified regression. No commit, extra test run, training, MD, or probe is started
as part of recording the freeze.

## Main-agent integration correction

The first unified run exposed test-order dependence: the geometry test asserted
that Torch was absent from the entire interpreter, although preceding training
tests legitimately loaded it. Replaced that global-state assertion with an
import guard around the actual geometry fixture and loop solve. A Torch import
by this operation now fails regardless of earlier tests; no chemistry acceptance
check was removed. The failed run remains in local regression evidence.
The final unified run is recorded by the public regression manifest.
