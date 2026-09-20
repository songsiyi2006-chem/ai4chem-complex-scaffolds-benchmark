# Phase 6–10 targeted audit fixes

Date: 2026-09-09. Scope: the six findings from the source audit, not certification
of the complete scientific workflows or laboratory hardware.

## Changes

- Phase 6: reconstruct free energy from deposited WTMetaD bias using
  `-gamma/(gamma-1)`, rejecting gamma <= 1. The previous factor produced 81%
  of the correct energy differences at gamma=10.
- Phase 7: use the UHF energy series, rather than RHF, for the UHF comparison.
- Phase 8: correct both terms in the derivative of
  `sqrt(det**2/4 + coupling**2)`. Forces now pass the isolated smooth two-state
  finite-difference check, including coordinate-dependent coupling.
- Phase 9: change tips between substrate/catalyst and methanol/internal-standard
  stocks. Split transfers also receive a fresh tip before revisiting stock.
  This increases tip consumption; installed tip inventory, geometry, solvent
  compatibility and actual robot operation still require preflight validation.
- Phase 10: convert electricity from W to kWh/s before combining with catalyst
  emissions in kg/s and dividing by production in kg/s.
- Phase 10: encode the already-corrected observation from `step`, instead of
  drawing another noisy sensor sample with a stale temperature derivative.

## Verification

`python -m unittest -v test_phase6_10_audit test_phase1_5_audit`

Result: 21 tests passed (7 new, 14 existing). New tests require NumPy; the existing
suite additionally requires RDKit and OpenMM. Tests extract actual function bodies
with AST to avoid legacy module import side effects and launching expensive jobs.
Generated robot code is compiled/inspected and its split-transfer helper is tested
with a contamination-tracking mock. No real robot was operated.

## Historical outputs are not revalidated

No full MD, quantum-chemistry, trajectory, RL training or hardware run was made.
Existing cached results, reports, figures and trained policies are historical,
not outputs validated against these fixes. Use a fresh output/cache location for
reruns; do not treat old output files as evidence that the corrected code ran.
Phase 6 free-energy analyses, Phase 7 method comparisons, Phase 8 trajectories,
Phase 9 exported protocols and Phase 10 training/evaluation must be regenerated
before making updated quantitative claims. This patch does not certify physical
model adequacy, convergence, exact-degeneracy handling or safe laboratory use.
