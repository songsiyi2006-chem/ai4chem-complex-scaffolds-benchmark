# Phase 14–17 targeted audit corrections — 2026-09-10

## Corrected calculations

- Phase 14: remove the repeated temperature division from reduced dissipation;
  use the actual mobility M0+M1*phi. Since volume fractions have no molecular
  surface-density calibration, outputs and plot labels now explicitly describe a
  reduced model proxy, not absolute entropy/area. Summary keys change from
  `S_*_kB_um2_s` to `S_*_reduced`; regenerate old output before plotting.
- Phase 15: use the density-matrix trace for survival, not the sum of all matrix
  elements. WHAM converts kJ/mol umbrella energies to kcal/mol, uses the correct
  offset sign, stable log-space normalization and fixed gauge, and rejects empty
  windows/nonconvergence. Latch population metric uses kT, not beta. Population
  ratios remain a limited diagnostic, not a certified binding free energy.
- Phase 16: remove the extra inverse distance from all three LJ force paths.
  Integrate minus the measured force, set the PMF reference to the sampled bulk
  mean, and retain attractive wells in the resistance integral. Reject missing
  bulk and invalid/overflowing PMFs. Convert the hydrodynamic radius from nm to m
  for the SI Stokes–Einstein calibration.
- Phase 17: Ex now differentiates to the existing X-alpha potential. Total
  energy removes the input Hartree/exchange potentials from the eigenvalue sum,
  then adds output-density Hartree energy and Ex. The final energy and density
  are recomputed from the final orbitals rather than copied from the prior iterate.
  This does not independently establish SCF convergence or model adequacy.

## Verification and remaining limitations

`python -m unittest -v test_phase14_17_audit`

Seven tests exercise syntax, reduced entropy, analytical singlet decay, synthetic
harmonic WHAM and latch units, all three LJ paths against energy finite differences,
transport PMF/bulk/SI conversion, and exchange derivatives/double-counting.
Tests use NumPy/SciPy and isolated real-source function bodies, without legacy
entrypoint side effects. The earlier regression suites remain applicable.

No full active-matter simulation, umbrella sampling, hydrodynamic simulation,
radical-pair production run or heavy-element SCF was rerun. Historical figures,
reports and cached numerical claims have NOT been revalidated. Use fresh output
locations and regenerate results before making quantitative claims. Existing
soft-core force caps, finite sampling/overlap, friction estimation, kinetic
integration and scientific model approximations are not certified by this patch.
