# Phase 27 — Multireference Ni photochemistry and spin-crossing branches

**Status: primary structures and methods audited; production multireference dynamics NOT run.** No new excited-state lifetimes, state populations, Ni–C/Ni–X branching quantum yields or transient spectra are established. The trajectory count is zero; no statistical convergence claim is made.

## Molecular series and recovered evidence

The main system is the untethered **Ni(bpy)(o-tolyl)Cl**, complex 1 in [Bím et al., Inorganic Chemistry 63, 4120–4131 (2024)](https://doi.org/10.1021/acs.inorgchem.3c03822). The [six-member series](results_phase25_27/phase27/series.json) includes its Br and I analogues, 4,4′-dimethyl-bpy and 4,4′-bis(trifluoromethyl)-bpy chloride derivatives, plus the structurally constrained cyclometalated Ni(PhBpy)Cl (published complex 2). Four derivatives are prospective designs, not a claim that all six were reported experimentally.

Five full literature XYZ files were extracted from [original SI](https://ndownloader.figshare.com/files/44610159): complex 1 singlet/triplet (36 atoms, S45/S46), complex 2 singlet/triplet (31 atoms, S47/S48), and the named tethered Ni(I) intermediate (30 atoms, S49). Atom counts were checked and selected original coordinate pages visually inspected. These are **published geometries**, not this run's optimized structures. Source SHA-256 and attribution are in the [manifest](results_phase25_27/literature_structures/manifest.json). Derivative geometries, coordination-isomer registries, solvent adducts and full atom/bond maps remain incomplete.

Primary spectroscopic anchors were located in SI S7: the untethered and tethered compounds have different absorption maxima in THF and toluene. The SI S29 method uses a CAS(10e,9o) reference, dynamic correlation, and extensive singlet/triplet averaging. These facts define reference checks, not an imported set of independent validation labels. The [2020 multireference study](https://doi.org/10.1021/acs.jpca.0c08646) and its [original SI](https://ndownloader.figshare.com/files/25563604) supply active-space/state-averaging sensitivity examples and show why a TDDFT-only bond-breaking interpretation is insufficient.

## Electronic-structure benchmark to be executed

Begin with the full neutral Ni(II) complexes in THF. Enumerate coordination arrangements, zero/one-THF adducts, solvent exchange and low-lying singlet/triplet states; test quintets if the state manifold warrants them. State labels should track wavefunction character, not assume that the nth energy-sorted root represents one fixed chemical state throughout a scan.

The published CAS(10e,9o) is a **starting benchmark**, not a certified space for this expanded network. An active space for competing Ni–C and Ni–X cleavage must cover relevant Ni 3d, occupied ligand π and π*, both scissile σ/σ* pairs and important halide/ligand donor orbitals. Determine electron/orbital counts from the complete molecule. Validate them using natural occupations, dominant configuration weights and, when feasible, orbital entanglement; inspect changes from equilibrium through dissociation/crossing regions.

Use a state-averaged multireference method with dynamic correlation (e.g. multistate/quasidegenerate NEVPT2 or a suitable multistate CASPT2 variant), checking basis, relativistic treatment, active-space expansion, root number, state weights and intruder/denominator behavior as applicable. Spin–orbit terms for Br/I require an appropriate relativistic Hamiltonian. If a larger space needs DMRG or another solver, converge its independent numerical controls. Report full-vs-fragment energy surfaces, state characters and couplings before deploying any truncated dynamics model.

Production surfaces need **energies and gradients** for all relevant spin states, nonadiabatic coupling within same-spin manifolds and spin–orbit coupling across spins. Use consistent electronic levels or quantify the approximation when dynamically corrected energies are combined with lower-level gradients/couplings. Such a mixture is a model error, not hidden behind a method label.

Track roots by overlaps and orbital/configuration character; align wavefunction phases and handle rotations in near-degenerate subspaces. Locate same-spin conical-intersection regions and spin-crossing seams. A minimal-energy crossing is a structural reference, not a population-transfer probability. Do not tune SOC/NAC values to force an experimental lifetime.

## Solvated nonadiabatic dynamics

The [preregistered protocol](results_phase25_27/phase27/protocol.json) proposes excitation windows 380–400, 440–460 and 510–530 nm. These are design choices. Initial states must be selected from calculated oscillator strengths and the actual pump spectrum; a uniform choice among electronic roots is not an excitation model.

Sample independent equilibrium solvent configurations and nuclear thermal/Wigner distributions. Define the QM region, solvent polarization, boundary treatment and force-field calibration. QM/MM or explicit-solvent direct dynamics must preserve the solvent cage relevant to radical separation/recombination; gas-phase trajectories cannot validate solvent escape. Use at least three independent batches of initial conditions and distinguish repeated hopping seeds from independent nuclear/solvent samples.

Execute a spin-mixed nonadiabatic scheme with documented time steps, electronic substeps, momentum rescaling, decoherence and frustrated-hop treatment. Proposed sensitivity grids are 0.5/0.25 fs nuclear and 0.02/0.01 fs electronic steps. Compare an appropriate alternative decoherence prescription, turn SOC off as a physical-model sensitivity control, expand the active space and extend observation times. Report energy drift, norm conservation, population positivity, failed trajectory counts and reasons for exclusions.

Direct ab-initio trajectories and trajectories on fitted surfaces must have separate identifiers, sample sizes and outcome intervals. Fitted energies alone are insufficient: validate gradients, SOC/NAC or the chosen diabatic representation, derivative consistency, crossing topology and phase behavior on independently held-out crossing/transition data. Never use validation trajectory outcomes to tune the couplings.

## Bond rupture, cages and competing events

Keep separate event histories for first Ni–C elongation/rupture, first Ni–X rupture, radical-pair formation, electronic relaxation, cage recombination, solvent capture and permanent products. Bond excursion requires state/bond-order evidence and a persistence/recrossing analysis. A long Ni–C distance at one frame is not a separated radical yield.

For the **tethered** complex, Ni–C cleavage leaves the aryl attached to the ligand framework. Free aryl escape cannot be assumed without an additional detachment/bond-breaking route. Intramolecular radical formation and recombination must be distinguished from true molecular separation; the permanent-separation channel needs a connectivity-compatible definition for each member.

Classify a final absorbing outcome only after appropriate persistence and connectivity checks. Trajectories that remain in excited/radical/cage states at the end are right-censored, not discarded or labelled failures. Record censoring reasons, horizon and importance weight. The [Aalen–Johansen implementation](phase25_27/analysis.py) retains censored samples and estimates competing-risk cumulative incidence under independent/non-informative censoring. Its software tests do not supply actual trajectories.

Quantum yield is per absorbed photon. Equal-weight trajectory counting requires correctly sampled initial absorption probabilities; importance-weighted sampling needs a corresponding weighted estimator and uncertainty analysis. Avoid treating every stored time frame as an independent photon. A finite-horizon incidence curve is not an asymptotic quantum yield, and unresolved censoring may require trajectory extension or bounded conclusions.

## Spectra and comparison with experiment

Use the same electronic ensemble and dynamics to calculate absorption, populations, survival/lifetime distributions and branch yields. Generate transient absorption from applicable ground-state bleach, stimulated emission and excited-state absorption, with state-dependent transition moments, excited-to-excited transitions, line broadening and instrument-response convolution. Populations alone cannot determine an ESA spectrum. Include radical/photoproduct absorption if its predicted abundance and spectral window matter.

Keep an explicit data-role registry: equilibrium spectral calibration, pump/response calibration, model selection and independently held-out transient/branch validation. A lifetime used to tune a surface or dephasing model cannot then be cited as successful external validation. Compare observables at matched solvent, temperature, excitation fluence/window and observation horizon, and distinguish fit uncertainty from electronic/dynamical/model uncertainty.

For convergence, start with ≥3 independent initial-condition batches, retain all outcomes/censoring and increase independent samples at preregistered checkpoints. Require 95% interval halfwidth ≤0.05 for main branches or explicitly report budget-limited nonconvergence. With independent, equal-horizon, fully observed binomial outcomes near probability 0.5, 385 trials give a Wilson halfwidth ≈0.0497; **this is a planning example, not the required trajectory count for a correlated/censored multibranch ensemble**. Use hierarchical/batch-aware bootstrap or an appropriate competing-risk uncertainty estimator. Optional stopping requires fixed checkpoints with multiplicity control or time-uniform intervals, not repeated unadjusted 95% stopping tests.

## Acceptance ledger and execution boundary

| Criterion | Status | Outstanding evidence |
|---|---|---|
| 27.1 real molecular series | Partial | Derivative/adduct/isomer structures and complete atom maps |
| 27.2 multireference benchmark | Literature audit only | New active-space/basis/state-average/dynamic-correlation convergence |
| 27.3 coupling quantities | Not run | Gradients, SOC/NAC, seams and state/phase continuity |
| 27.4 solvated dynamics | Protocol only | Three-window trajectories, convergence, fitted-vs-direct tests |
| 27.5 effective photochemistry | Estimator only | Persistent/cage-resolved outcomes with censoring |
| 27.6 observables | Not run | Consistent spectra, lifetimes, branch yields and independent validation |
| 27.7 precision/sensitivity | Not passed | ≥3 batches, branch CI target, SOC/decoherence/space/time checks |

The existing Phase7 small-active-space and Phase8 reduced-model workflows can supply implementation patterns, not quantitative truth for Ni photochemistry. No SOC-aware multireference trajectory backend is connected here. Final conclusions must ultimately separate electronic-structure error, dynamical approximations, initial sampling and finite-time censoring. None of those uncertainties has been quantified for this series in the present run.
