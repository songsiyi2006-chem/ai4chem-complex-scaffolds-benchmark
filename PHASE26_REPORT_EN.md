# Phase 26 — Constant-potential C–C coupling at Cu interfaces

**Status: interface design and analysis checks prepared; no constant-potential production calculation performed.** There are no new activation free energies, partial currents, Tafel slopes, reaction orders or product Faradaic efficiencies. Neither CO–CO nor CO–CHO dominance is established.

## Completed work and boundary conditions

The [48-condition matrix](results_phase25_27/phase26/interface_matrix.csv) covers Cu(100), Cu(111), Cu(211), four potentials (−0.4, −0.7, −1.0, −1.3 V vs SHE), K⁺/Cs⁺, and 0.125/0.25 ML CO. These are proposed conditions: 298.15 K, bulk pH 6.8, 0.10 M bicarbonate electrolyte. Three independent interface seeds per condition would require **144 solvated initial states**, before path windows or convergence tests. None has been sampled.

Three ASE [dry slab seeds](results_phase25_27/phase26/dry_slab_seeds/manifest.json) were generated (Cu100: 64 atoms; Cu111: 64; Cu211: 96). They have no water, ions, adsorbed CO or defined electrode potential and are not relaxed interfaces. The provisional 3.61 Å lattice parameter must be reoptimized at the production functional. A dry slab file does not meet the interface requirement.

Define CO coverage as CO per exposed Cu site and separately report step/terrace occupations. The stepped surface's denominator must come from a declared exposed-site registry, not the number of atoms in the slab. For cation comparisons, match chemical potential, anion, ionic strength, potential, interfacial water volume, CO occupancy and equivalent adsorption configuration. A single ion placed at an arbitrary distance is not a cation-ensemble calculation. Bulk pH is a boundary condition; local pH follows from transport and reaction.

Potential conversion is `U_RHE = U_SHE + (kBT/e) ln(10) pH`. Store the reference, calibration and applied iR treatment explicitly. Do not mix SHE slopes with RHE slopes across changing pH.

## What the supplied reference does and does not establish

[Cheng et al., Nature Communications 16, 4064 (2025)](https://doi.org/10.1038/s41467-025-59267-3) performs potential-dependent transition-state searches within a grand-canonical treatment. Its Methods uses an implicit electrolyte, while its kinetic discussion adopts CO dimerization as the rate-determining step instead of constructing the complete microkinetic network requested here. It also distinguishes changes in adsorption/coverage from intrinsic coupling barriers. These are useful methodological benchmarks; they do not provide new data for the present K/Cs × facet × potential × coverage matrix. The [original supplementary information](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-59267-3/MediaObjects/41467_2025_59267_MOESM1_ESM.pdf) is linked for its geometries and electrostatic sensitivity analyses.

## Grand-canonical electronic structure

Choose a verified grand-canonical engine or a validated constant-potential equivalence at **each geometry and path region**. Store nuclear coordinates, cell, pseudopotentials, basis/cutoffs, k-grid, exchange-correlation/dispersion definition, electrolyte model and potential-calibration parameters. PBE-D3(BJ)/RPBE-D3(BJ) are proposed sensitivity levels; their difference must be reported separately from numerical sampling error.

For fixed nuclei/ions and a common electronic reference, write `Ω(U)=G(N_e)-μ_e(U)(N_e-N_e,0)`. Reservoir terms for different proton/ion/water compositions must be included explicitly. Constant `N_e,0` terms can be dropped only for comparisons of the same composition and gauge. Total cell charge in units of +e is `N_e,0-N_e`. A partitioned electrode charge is a separate observable and need not equal total cell charge because electrolyte and adsorbate charges also contribute.

The [analysis implementation](phase25_27/analysis.py) rejects missing electron numbers, inconsistent charge signs and potential deviations beyond 0.05 V. It returns a pointwise Legendre audit, **not a solvent activation PMF**. The electronic chemical potential must be aligned with the actual electrostatic calibration at every point, not inferred by a universal vacuum-reference constant without checking the solvent model.

If an equivalent fixed-charge scheme is used, sample a sufficient charge range separately for reactant, product and transition-region geometries, solve for the target potential, quantify fit/capacitance uncertainty and validate selected points against direct constant-potential calculations. One shared linear correction to all barriers is prohibited. Image-specific potential and charge should be converged during the transition-state search and solvent sampling.

## Explicit interfacial sampling and competing pathways

The [protocol](results_phase25_27/phase26/protocol.json) specifies CO–CO and CO–CHO coupling, CO→CHO/COH proton-electron transfer, and Volmer/Heyrovsky/Tafel hydrogen-evolution competition. Add CO₂ adsorption/conversion to CO and exchange with solution to predict CO coverage. Include reverse reactions and site occupancy for every elementary event.

Equilibrate at least three independent water/cation arrangements per condition. Test water thickness (initial proposal: 15/20/25 Å), doubled lateral area at fixed physical density, slab thickness, fixed-layer choice and ionic reservoir treatment. Distinguish a restrained occupancy calculation from an equilibrated coverage calculation; comparing their apparent rates requires the microkinetic state model.

Use enhanced sampling or constrained dynamics for the leading activated processes. Candidate collective variables are C–C separation, C–H/O–H coordination, water hydrogen-bond connectivity and cation coordination. Add proton donor identity or water-wire coordinates if committor tests reveal hidden barriers. A C–C distance alone can miss the solvent/proton bottleneck. Validate the dividing region with unbiased shooting trajectories, recrossing analysis and additional coordinates.

For each PMF, retain bias definitions, window centers and spring constants, trajectory timestamps, burn-in choices, autocorrelation estimates, effective samples, neighboring-window overlap, forward/reverse checks, replicate estimates and bootstrap interval. Reference constrained free energies to the same adsorption basin with consistent translational/volume and Jacobian terms. Separate window resampling from independent-interface variability; correlated windows are not independent replicates.

If a fitted potential is introduced, independently test energy, forces, electron-number response, differential capacitance and work-function/potential response over adsorbed reactant, TS and product environments. A neutral fixed-charge force model is not validated for variable-potential chemistry merely because its force MAE is small. Hold out crossing/transition-region structures and ion/water arrangements from fitting.

## Microkinetics and mechanism attribution

Construct a surface-specific model with a declared site registry and occupancy matrix. Enforce site conservation and reservoir-adjusted local detailed balance; calculate reverse rates from the same transition states. Include pair/multisite availability rather than replacing every two-site event by an unqualified `θ²` law. Model lateral interactions consistently in state energies and chemical potentials, and propagate their uncertainty.

Solve coverages and rates jointly with transport for CO₂/CO, proton donors, bicarbonate/OH⁻ and products, with stated diffusion-layer thickness or hydrodynamic conditions. Report site density and real/geometric area conversion. Cathodic-current sign conventions must be declared; electron-consuming net reaction fluxes determine current. C–C bond formation itself can transfer zero net electrons, so its flux cannot be multiplied by an arbitrary two-electron number to make a partial current.

Final-product FE requires the downstream network to each named product, its electron count and competing losses. A truncated CO–CO/CO–CHO network can report conditional intermediate formation fluxes but not ethylene/ethanol FE. Transport-limited and intrinsic rates must be distinguished.

Calculate reaction orders by perturbing relevant activities at matched potential and re-solving the coupled model. Derive apparent Tafel slopes from `dU/dlog10|j|` over an identified monotonic regime, documenting curvature, coverage changes and transport. A single barrier derivative is not automatically a macroscopic Tafel slope.

For degree of rate control, perturb each TS energy at fixed state free energies, thereby changing both forward and reverse rates consistently, and re-solve the steady state. Also assess thermodynamic control of intermediates by separate state-energy perturbations. Compare constant-charge/constant-potential, implicit/explicit-solvent, coverage and cation-position controls. Do not identify a rate-determining step from the largest uphill reaction free energy or attribute increased coverage to an intrinsic barrier decrease.

Propagate correlated electronic-structure, potential-calibration, lateral-interaction and sampling errors into fluxes and slopes. At 298.15 K, a **hypothetical** 0.05 eV barrier shift changes a transition-state-theory rate by a factor of ≈7.00; this analytical sensitivity is not a simulated current and shows why apparently small barrier uncertainty can prevent mechanism discrimination.

## Acceptance and discriminating observations

| Criterion | Current state | Required evidence |
|---|---|---|
| 26.1 comparable interfaces | Matrix/dry seeds only | Fully solvated, relaxed, calibrated structures and occupancy registry |
| 26.2 constant potential | Audit code only | Electron number, charge, actual U and Legendre transform per path point |
| 26.3 solvent rearrangement | Not run | Water/cation sampling; area/thickness/initial-state convergence |
| 26.4 competition | Specified only | Activated PMFs and reaction-coordinate validation |
| 26.5 macroscopic kinetics | Not run | Site-conserving coupled microkinetics/transport and downstream product network |
| 26.6 attribution | Protocol only | Matched controls and uncertainty-aware rate-control analysis |
| 26.7 accuracy | Not passed | Numerical ≤0.05 eV; sampling 95% halfwidth ≤0.05 eV; actual-U deviation ≤0.05 V; ≥3 starts |

No mechanism is selected. Once predicted intervals exist, choose discriminating potentials where the two pathways' observable intervals separate, rather than choosing them after inspecting experiments. Candidate tests are potential-dependent H/D kinetic isotope effects at matched activities, CO activity reaction orders at measured coverage, and time-resolved IR/Raman constraints on CO/CHO-related intermediates. Solvent isotopic changes and electrochemical Stark shifts need their own calibration. Sensitivity and detection thresholds must be computed; an assigned spectral peak alone will not establish a mechanism.

Production remains blocked by an unconnected periodic grand-canonical/explicit-solvent backend and the absence of sampled data, not by a completed computation with an inconclusive answer.
