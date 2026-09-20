# Bounded Phase19 orientation correction

## Coordinator full geometry probe (not full evolution)

Fresh campaign probe `20260910T233318` completed in 18.9 s with exit 0.
The immutable six-anchor RMSD is **0.191046195 A**, below the unchanged 0.30 A
geometry gate; TRP NE1 error remains **0.467965695 A**. Target and reachability
coordinates were not modified and chirality violations are zero. However,
the realized static model has **661 heavy-atom clashes** and Ramachandran
fraction 0.82456. This is only a geometry-probe pass, not an accepted candidate,
full training/evolution, stability proof or experimental enzyme claim.
Original failed 1.013323 A probe remains preserved. New source SHA256:
`15aaffc5af5876a4492ad42fda54afc417fb180ad87665bd07819cd8ddf5d91c`.

Scope: only Phase19 source, `test_phase19_orientation.py`, and this note.
No training, complete geometry-probe, MD, QM, commit, or unrelated-file edit.

## Evidence and construction

The active Trp call in `place_free_bundle` used `htype="CB"`, hence generic
CG-direction aiming. The special 0.62/0.79 Trp offset was dead code; it must
not be described as the active cause. The active generic aiming also fails
to optimize the actual reachable NE1 band.

The existing analytical helper gives a sphere of radius r centered on CB,
with axial projection z in [zmin,zmax] along CA->CB. With b=|CA-CB|, maximum
CA reach is R=sqrt(r^2+b^2+2*b*zmax). For a target L>R from immutable CA,
the reverse triangle inequality gives NE1 error >= L-R. Equality requires
the maximal-reach NE1 vector to point toward the target, which requires
CB_hat dot target_hat = (b+zmax)/R.

A rigid roll about the supplied rod axis changes that dot product as
A*cos(theta)+B*sin(theta)+C. New `trp_optimal_anchor_roll` solves this
equation analytically and chooses the smaller absolute rotation of the two
solutions. If no solution exists, or L<=R, it keeps the original rod. This
is a global NE1-distance optimum when the equation is attainable, not an
optimizer of packing, ring stacking orientation, or full enzyme feasibility.

The Trp call now selects TRP mode; it starts with the existing generic aim
and applies this rigid roll. Existing substrate-clearance filtering runs
afterward. All atoms in the rod undergo one proper rigid rotation about its
fixed CA. Target generation, bonds, angles, stereochemistry, chi packing,
other motif rods, and all acceptance thresholds are unchanged. The unused
arbitrary Trp offset branch was removed.

## Cheap actual-source fixture

AST-extracted production builders and `_trp_ca_target` with an in-memory
RDKit/MMFF indole, seed 485; no top-level Phase19 imports. Translated stack
centroid [0,0,4] and the constructor's 25-degree in-plane direction. Rod
length 18 and anchor index 1 match production. A distant placeholder
substrate isolates orientation; a separate test exercises collision rejection.

- Original fixed-backbone analytical minimum NE1 error: 2.464216272 A.
- Optimized analytical minimum: 0.467965684 A.
- Independent three-start chi-only least-squares fit: 0.467965684 A.
- Conditional six-anchor RMSD, if the other five errors are zero: 0.191046191 A.
- This fixture's original minimum is not the saved full probe's 2.482123656 A
  observed packing error; no saved full result is overwritten or reclassified.

NE1 alone still misses 0.30 A. The aggregate lower bound is below 0.30 A,
but this does not establish a complete candidate pass. Current packing,
downstream realization, full scaffold clashes, and other scientific gates
have not been evaluated by a new full probe. Existing best full-probe RMSD
remains 1.013322739 A until a fresh authorized full probe measures otherwise.

Tests cover the analytical bound across twelve perpendicular axes, actual
chi fitting, immutable targets/CA, all pairwise backbone/sidechain distances,
chirality, conservative no-op cases, active call wiring, and collision rejection.
Use molecular Python with its Library/bin on PATH and one BLAS/OpenMP thread:
`python -m unittest -v test_phase19_orientation`. Run under a subprocess
timeout of 30 seconds. No Torch/OpenMM imports or full production entrypoints.
Final run: four tests passed in 1.309 seconds (1.944 seconds shell wall time).
Owned production diff whitespace check passed.
The first primary-Python attempt lacked SciPy; mixing in Conda site-packages
then failed NumPy DLL loading. The native molecular interpreter resolves both
without installation or environment changes.
