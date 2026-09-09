# Phase 11–13 targeted corrections — 2026-09-10

- Phase 11 removes the global log-amplitude floor. Sampling, kinetic
  differentiation and signed wavefunction evaluation now use the same expression.
  This is not a certification of exact-node handling or training convergence.
- Phase 12 moves the anti-collapse penalty from variation along an orbit to
  variation on off-orbit training states. This is a regularizer, not proof of a
  nontrivial Hamiltonian. Independent initial conditions (seed 120023) now supply
  the held-out scalar evaluation. Report strict nonincrease fractions and positive
  derivative violations rather than a global certificate. Drift is reported in
  absolute units and relative to a fixed 0.25 target scale, not its own range.
  The figure uses absolute drift too. Legacy result keys `V_certificate_rate`
  and `H_max_drift_60tau` are replaced with explicitly named diagnostics.
- Phase 13 uses the Euclidean dot product of normalized finite-difference
  eigenvectors, without a second grid-spacing factor. The potential tilt now
  satisfies V(acceptor)-V(donor)=dG at the specified reference well positions.
  The grid convergence table includes H/D rates as well as energies.

## Tests and limits

Run `python -m unittest -v test_phase11_13_audit` (NumPy/SciPy).
Seven regressions pass: source compilation, unclipped log-amplitude expression,
scalar audit signs/scales/invalid inputs, transverse constraint wiring, overlap
normalization, potential direction, and H/D rate grid refinement (10% tolerance).
The combined Phase 1–13 regression command is:

`python -m unittest -v test_phase1_5_audit test_phase6_10_audit test_phase11_13_audit`

The first suite also needs RDKit/OpenMM. The audit executes isolated actual source
functions without running legacy entrypoints. PyTorch is unavailable in the tested
runtimes, so neural training and autograd integration were not rerun. No full VMC,
law-discovery or electronic-structure calculation was performed. Held-out samples
and empirical inequalities cannot prove global stability or conservation.

Historical reports, figures, checkpoints and cached results have NOT been
regenerated. They must not be treated as post-fix evidence. Use fresh output/cache
locations for full reruns and regenerate reports before making quantitative claims.
These changes address the six reported findings, not every modeling limitation.
