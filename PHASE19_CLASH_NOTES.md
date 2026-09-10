# Phase19 bounded clash correction and diagnostic

## Result and scope

Changed only `run_phase19_active_inference_denovo_enzyme.py`,
`test_phase19_clashes.py`, and this file. No full probe, training, MD, QM,
new environment, commit, push, running snapshot edit or other project edit.
Used chem-ai4s runtime guidance and the existing phase2ff interpreter. Each
test command runs with one numerical-library thread and a 30-second subprocess
timeout. No subagent API was available in this worker; tests ran sequentially.

Two bounded packing corrections are implemented and tested: exact overlaps
with other residues now receive the existing 2.9 A penalty, and later rotamer
choices see earlier accepted sidechains, including alanine CB atoms. Candidate
rotamer sets, iteration order, seed, targets, catalytic fitting, backbone coordinates,
chirality enforcement, collision cutoffs and existing scientific gate thresholds are
unchanged. A follow-up now adds a hard constructor rejection before those gates,
as detailed below. Accepted rotamers may change, as intended by the corrected score.
No full-candidate collision reduction is claimed.

## Immutable probe evidence

Read-only source: `audit/phase1_19_rerun_20260910/phase19/20260910T233318/`
`results_phase19/geometry_probe.json` (and its sibling source snapshot).
Recorded source SHA256:
`15aaffc5af5876a4492ad42fda54afc417fb180ad87665bd07819cd8ddf5d91c`.

- Six-anchor RMSD: 0.1910461949933137 A; TRP NE1 error: 0.46796569503387686 A.
- 1477 heavy atoms, 172 residues, 661 legacy clashes below 2.65 A.
- Ramachandran fraction: 0.8245614035087719; chirality violations: zero.
- The snapshot contains a summary JSON and log, but no atom coordinate artifact.
  It cannot establish how many of these 661 pairs are backbone/backbone,
  backbone/sidechain or sidechain/sidechain without another allocated computation.

The historical audit masks **all** atom pairs within the same or adjacent
residues. Consequently, on the standard single peptide chain used here, none
of its counted 661 pairs can be direct bonds or 1-3 pairs. They are at least
661 overlaps under the existing cutoff, not a bonded-pair overcount. The
blanket residue mask instead hides local nonbonded overlaps. The saved
`geometry_gate_passed=true` reflects the existing permissive static gate
(`n_clashes <= 5000`), not a physically collision-free scaffold.

## Additional diagnostic, with existing gates preserved

`static_fold_audit` now also returns `clash_diagnostic`. Its standard-residue
bond graph includes peptide bonds, branches, aromatic ring closure and terminal
OXT. It excludes graph distance one and two only; bonds are never guessed from
short Cartesian distances. Missing coordinate intermediates retain their
chemical topology. The new report includes:

- Total nonbonded overlaps, nonlocal legacy-equivalent count and previously
  hidden local count, using the unchanged strict distance < 2.65 A rule.
- Backbone/backbone, backbone/sidechain and sidechain/sidechain counts.
- Number of close bonded/1-3 pairs excluded, and up to 12 worst atom pairs
  with distances in angstroms and zero-based residue indices.

The existing `n_clashes` field and its use in gates are deliberately preserved;
the topology-aware count is additional evidence, not substituted acceptance.
Hydrogens are excluded from the new diagnostic. Unsupported residue names and
nonfinite coordinates raise errors. This graph assumes one standard peptide
chain without disulfides, covalent ligands or modified residues. It does not
audit ligand clashes or validate bond lengths. Row-wise distance evaluation
uses O(N) temporary memory, O(N^2) time and bounded pair output. The pre-existing
legacy audit still allocates its original dense distance matrix.

## Backbone limitations found in actual source

`place_free_bundle.rod_cost` uses CA-only distances, and its fallback accepts
the minimum-cost motif trio when none is clash-clean. `realize_backbone` checks
a gross nonlocal CA separation of only 0.5 A. The spline loop selector considers
already-placed CA coordinates rather than all backbone atoms or future rods.
Its `n_next` argument is unused; loop N and C follow spline tangents rather
than enforcing peptide closure. These are plausible collision/strain mechanisms,
not a measured attribution of the saved 661 pairs.

A cheap execution of the actual nested spline function, endpoints
CA=[10,0,0] and [27.2,0,0] A with center=[0,0,0], produces three residues
with N-CA-C angles [180,180,180] degrees. Changing the supplied next N leaves
all coordinates identical. Fixing this requires constrained backbone closure
and subsequent geometric evaluation, which is beyond this worker's allocation.
The follow-up below fixes silent acceptance of invalid construction. No repaired
backbone or physical feasibility is claimed.

## Follow-up: hard actionable construction rejection

The requested bounded physical correction is not established: changing the local
N/CA/C frame alone cannot simultaneously close both junctions at fixed motif
anchors. The current spline has no peptide-closure constraint, ignores `n_next`,
and cannot be treated as a validated peptide constructor. Consequently,
`realize_backbone` now calls `require_constructed_backbone` after assembly and
before torsion evaluation, sequence design, packing or simulation. Invalid
construction raises `AssertionError` starting with:

`BACKBONE_REJECT: implement peptide closure/collision resolution; ...`

This uses the existing candidate-construction rejection path. A real probe
propagates the failure instead of writing a new success result. Existing archived
results are not reclassified, deleted or overwritten. The message contains total
internal-geometry violations, total backbone overlaps, up to three geometric
examples and the closest clash pair. Candidate logs may truncate details using
their existing limits, but the actionable rejection prefix and counts come first.

The check explicitly enforces the **ideal generated-coordinate contract**, not
experimental or minimized protein geometry: N-CA=1.458 A, CA-C=1.525 A,
C-O=1.231 A and inter-residue C-N=1.335 A; N-CA-C must match the existing
builder's 111.0-degree initial residue or 111.2-degree subsequent residue.
The 1e-4 A/degree comparison tolerance permits floating-point construction noise,
not structural relaxation. All reference values come from the existing NeRF
builder. Missing/nonfinite backbone coordinates also fail.

Backbone nonbonded overlaps below the unchanged 2.65 A cutoff are an additional
hard construction rejection, excluding topology distances 1 and 2 properly.
**Eligibility is intentionally tightened by this new prerequisite**; the old
`n_clashes <= 5000`, CA-distance, RMSD, chirality and Ramachandran gates are
not edited or weakened. No coordinates or random choices are changed. Passing
this limited contract still does not establish every bond angle, peptide
planarity, loop torsion quality, foldability or catalytic activity.

This can reject the existing spline's candidates broadly. Do not interpret it as
a successful physical repair or schedule an expensive run expecting closure to
be fixed. The actionable remaining work is a constrained internal-coordinate
loop solver that preserves the motif anchors, satisfies both peptide junctions
and resolves backbone intersections. Main owns that scientific redesign and
any larger allocation.

The small integration fixture executes actual `realize_backbone` and the probe
caller with deterministic prebuilt rods replacing expensive upstream generation.
It is **not a full geometry probe**. It measures 28 internal-geometry violations
and 160 backbone clashes; examples are C-O=1.324713 A, N-CA-C=143.943710 degrees,
and junction C-N=0.688527 A. The constructor rejects before downstream gates;
the probe caller cannot reach design, packing, static success reporting or file
output. Original rods and target coordinates remain exact. A clean actual-source
12-residue canonical helix passes unchanged. The 180-degree spline fixture,
broken-junction fixture and 124-overlap fixture all reject explicitly.

Packing remains a greedy, order-dependent approximation. The existing score
still excludes all own-residue pairs and omits CB from ordinary candidate
scoring; no claim of globally optimal packing or removal of all internal
sidechain strain is made. The new audit exposes these potential local overlaps.

## Exact verification

PowerShell environment (process-local only):

```powershell
$env:OMP_NUM_THREADS='1'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:NUMEXPR_NUM_THREADS='1'
$env:PYTHONDONTWRITEBYTECODE='1'
$env:PATH='C:\Users\HUIWEI\miniconda3\envs\phase2ff\Library\bin;'+$env:PATH
& 'C:/Users/HUIWEI/miniconda3/envs/phase2ff/python.exe' -c 'import subprocess,sys; subprocess.run([sys.executable,"-m","unittest","-v","test_phase19_clashes"],timeout=30,check=True)'
& 'C:/Users/HUIWEI/miniconda3/envs/phase2ff/python.exe' -c 'import subprocess,sys; subprocess.run([sys.executable,"-m","unittest","-v","test_phase19_orientation"],timeout=30,check=True)'
```

Final follow-up clash run: **10 tests passed in 1.083 s** (shell wall 2.314 s).
Orientation regression: **4 tests passed in 1.696 s** (shell wall 2.719 s).
Tests AST-extract actual production functions/blocks without importing the
Phase19 entrypoint, Torch or OpenMM. The packing test reuses the existing
orientation test's small RDKit/MMFF fixture and executes the production
scaffold packing block, skipping catalytic inverse fitting.

Coverage includes an independent peptide adjacency-matrix 1-3 oracle, missing
ring intermediates, retained 1-4 clashes, exact cutoff, hydrogens, invalid
inputs, bounded output, static-audit integration and unchanged coordinates.
Two actual-source six-residue helices offset by [0.1,0.1,0.1] A yield **124
backbone/backbone clashes**, identical to independently recomputed historical
count; **100 close bonded/1-3 pairs** are excluded. The clean six-residue helix
has zero clashes. These are controlled fixtures, not the saved full probe.

The packing fixture verifies every accepted sidechain is present in the
reference with its correct owner, all backbone coordinates remain exact, and
a coincident other-residue atom contributes the full existing 2.9 penalty.
The orientation regression retains the 0.467965684 A analytical NE1 bound and
0.191046191 A conditional six-anchor RMSD. These tests demonstrate software
corrections and preserve known geometry; they do not validate a complete enzyme.

Main/coordinator must allocate any new full probe and assess its new diagnostic,
and owns integration/runtime. The immutable 661-clash result remains unchanged.
