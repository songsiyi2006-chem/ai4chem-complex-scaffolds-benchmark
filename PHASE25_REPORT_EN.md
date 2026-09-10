# Phase 25 — Conformer ensembles and blind selectivity prediction

**Later calculations:** [current progress report](PHASE25_PROGRESS_EN.md). The body below preserves the initial `93b9ab1` release snapshot; it is not the latest completed-job count. Overall scientific acceptance remains incomplete.

**Status: preliminary implementation and local pilot; scientific acceptance NOT achieved.** No catalytic rate, conversion, absolute product configuration or signed ee is currently established. All 96 proposed conditions remain **selectivity undetermined**. The local single-point ledger records actual successes, failures and interruptions separately.

## Evidence produced in this run

- Six complete BINOL-phosphate connectivities, each with a coordinate mirror; 12 intended axial entries. Substituents at both 3 and 3′ are phenyl, 4-methylphenyl, 3,5-bis(trifluoromethyl)phenyl, 2,4,6-triisopropylphenyl, 9-anthracenyl and 1-naphthyl. The fluorinated member's connectivity passes an exact canonical-graph comparison to the [supplier's original structure](https://www.sigmaaldrich.com/US/en/product/aldrich/674605).
- Four ketimines: N-phenyl derivatives of acetophenone, 4-methylacetophenone, 1-acetylnaphthalene and 2-acetylnaphthalene, each explicitly E and Z. Toluene and dichloromethane give 96 catalyst–substrate–solvent conditions; E/Z are network species, not independent experimental labels.
- 144 ETKDG/MMFF starting conformers requested: 72 free catalyst, 64 imine and 8 dimethyl Hantzsch ester starts. Coordinates and individual convergence flags are retained in [search audit](results_phase25_27/phase25/search_audit.json). Search frequency is not physical degeneracy. No thermodynamic weight is assigned to these force-field structures.
- A [19-species, 26-reversible-edge topology](results_phase25_27/phase25/network.json), algebraically checked for catalyst-, imine- and donor-derived moiety conservation. All production free energies remain null.
- Tested kernels for partition sums, reversible concentration-dependent kinetics, correlated uncertainty, conditional fast-exchange diagnostics, and same-surface/IRC metadata rejection. Analytic fixtures verify software; they are not chemical measurements.

## Identity, stoichiometry and reference provenance

### New local electronic single points

The completed calculations use the same full, isolated 28-atom N-phenyl acetophenone ketimine E/Z starting geometries, def2-SVP, gas phase and a 75×302 quadrature grid. Energies below are **electronic single points on MMFF geometries**, with no geometry/frequency/solvent validation:

| Method | E isomer, Hartree | Z isomer, Hartree | E(Z)−E(E), kcal/mol |
|---|---:|---:|---:|
| PBE-D3(BJ) | −594.9027949712 | −594.8965536426 | +3.91649 |
| ωB97X-D | −595.4352432708 | −595.4301707678 | +3.18304 |

Both methods rank the selected E starting geometry lower electronically. Their difference of about 0.73345 kcal/mol in this comparison is not a statistical confidence interval. Neither these E/Z energy differences nor their method spread establishes E/Z solution populations, catalytic barriers, R/S selectivity or the required ≤0.5 kcal/mol ensemble convergence. Successful jobs used 647.872 allocated core-seconds in total; interrupted/failed attempts are separate and total campaign cost is incomplete. No cost-matched active-learning claim follows.

The donor is **dimethyl 2,6-dimethyl-1,4-dihydropyridine-3,5-dicarboxylate**, explicitly distinguished from its diethyl analogue. Net stoichiometry is ketimine + reduced Hantzsch ester → secondary amine + aromatic pyridine ester; phosphoric acid is regenerated. Initial design conditions are 298.15 K, imine 0.10 M, donor 1.20 equivalents and catalyst 5 mol%. These are proposed computational conditions, **not recovered experimental conditions**.

Heavy-atom maps are in the component SMILES, with disjoint hydrogen-map ranges in SDFs. Catalyst axis atoms are maps 14/15; O-bearing reference carbons are maps 5/16; substitution sites are maps 6/18. Geometric torsion sign and coordinate reflection are recorded. **R_a/S_a certification, actual atom-resolved reactive H transfer, and product CIP assignment remain pending**; torsion sign is not silently equated with absolute axial configuration. The neutral acid, neutral ketimine and neutral donor are starting components. The proposed ion-pair states are net-neutral iminium/phosphate aggregates; solvent-separated ion populations must be tested.

The [Simón–Goodman primary paper](https://doi.org/10.1021/ja800793t) motivates acid/imine/donor cooperation. Its [original SI](https://ndownloader.figshare.com/files/4627135) was retrieved and hashed. The 28-atom E-imine at SI S2–S3 was extracted and visually checked against the coordinate pages. The SI also contains a reduced phosphoric-acid model: its energies cannot be treated as full-catalyst barriers. Complete catalyst-specific experimental conditions and independently assigned ee labels have not been recovered. Source-derived geometries are labelled `literature_geometry`, never new calculations. See [source manifest](results_phase25_27/literature_structures/manifest.json).

## Required production calculation

1. Generate multistart acid–imine, acid–donor and ternary complexes, neutral and ionic microstates, E/Z isomers, proton-relay states, product-bound complexes, and aggregates. Search both product faces independently, then mirror the complete system as a control. Record failed and duplicate searches, atom maps and connectivity changes.
2. Cluster conformers by mapped geometry and interaction topology. Use physical symmetry degeneracies once per basin; do not count repeated optimization hits as statistical multiplicities. Broaden energy windows until omitted basin weight and the R/S free-energy difference are stable.
3. Use two distinct dispersion-inclusive approximations, such as PBE-D3(BJ) and a range-separated hybrid, with consistent basis/solvation and optimization, frequencies and IRC on **each method's own surface**. The local pilot uses PBE-D3(BJ) and ωB97X-D/def2-SVP on isolated MMFF imine geometries, in gas phase; it is not this production protocol. [Actual pilot ledger](results_phase25_27/phase25/local_qm_pilot/results.json).
4. For selectivity-critical full-system structures, use a suitable correlated method (e.g. tight local CCSD(T), with localization/basis checks and single-reference diagnostics). Verify reduced/full corrections on both competing channels. If diagnostics reject a single-reference treatment, change the reference hierarchy; do not quote CCSD(T) by default.
5. Accept a TS only with exactly one reaction-relevant imaginary mode and forward/reverse IRC endpoints matching the intended atom-mapped states. A single point, guessed bond distance or endpoint-energy scan is insufficient. Audit the raw output in addition to the metadata prefilter.
6. Assemble consistent 1 M solution Gibbs energies, solvent treatment, standard-state changes, quasi-harmonic/hindered-rotor corrections and symmetry terms. Report 50/100 cm⁻¹ low-frequency-treatment sensitivity. Do not combine translational binding entropy twice with concentration-dependent activities. At 298.15 K, the ideal 1 atm→1 M correction is +1.89433 kcal/mol **per species**; reaction corrections depend on stoichiometry.
7. Select the two closest R/S competitions from preliminary DFT, then perform explicit-solvent umbrella/constrained sampling on both channels, with at least three independent solvent starts. Include proton-transfer and solvent-reorganization coordinates, overlap diagnostics, equilibration cuts, autocorrelation-aware effective samples, and bootstrap confidence intervals. Repeat forward/backward or alternative-CV sampling and inspect recrossings.

## From free energies to observables

For basin energies `G_i`, physical degeneracies `d_i` and temperature `T`:

`G_ensemble = -RT log Σ_i d_i exp(-G_i/RT)`.

Use stable log-sum-exp and a unified reference composition. For every reversible edge, the **same TS free energy** determines both directional rates, enforcing local detailed balance. The mass-action implementation uses dimensionless activities `c/c°`, and flux units mol L⁻¹ s⁻¹. Association may need capture/diffusion models; submerged barriers are rejected rather than clipped to zero. Fit neither an arbitrary R/S offset nor a target ee into the network.

Propagate correlated basin/TS errors through the full transient network, preserving shared reactant and method errors. Report rates, initial E/Z ratio, conversion and accumulated product amounts at explicit times and concentrations. `ee=(R-S)/(R+S)`; zero product gives undefined ee. A 95% ee interval spanning zero gives **selectivity undetermined**. An interval based only on conformer sampling omits functional/model error and must be labelled accordingly.

The conformational generator's spectral mixing rate must exceed all escape hazards over the relevant condition range. A factor of 100 is a proposed conservative screening threshold, not a proof. Check irreducibility, reversible populations, concentration dependence, initial-condition loss and full-network/reduced-model agreement before invoking Curtin–Hammett. Product inhibition, slow E/Z exchange or aggregation can defeat a single ΔΔG‡ description.

Only when validated may `ee=tanh[(G_S‡-G_R‡)/(2RT)]` be used. Otherwise `RT ln(R/S)` is an **effective product-ratio coordinate**, not a mechanistic barrier. As an analytical sensitivity illustration only, 0.5 kcal/mol corresponds to |ee|≈0.3986 at 298.15 K under the two-channel identity; a 0.5 kcal/mol convergence target therefore does not by itself secure the major enantiomer.

Controls still requiring physical calculations: full coordinate mirror (opposite ee, equal rate), achiral diphenyl phosphate (zero ee in an achiral environment), catalyst-free chemistry, catalyst-concentration series for aggregation and product-spiking series for inhibition. Constructing a mirror or proving an algebraic symmetry does not replace these controls.

## External holdout and Phase24 integration

The [fixed split](results_phase25_27/phase25/conditions.csv) has 32 training, 16 catalyst-family holdout, 32 substrate-family holdout and 16 double-holdout conditions. Fused-polyaryl 3,3′ groups and naphthyl ketimines are held out. This is a small, declared domain shift; it does not establish transfer to all catalyst backbones or N-substituent classes. Both enantiomers of every pair stay together. The [blind protocol](results_phase25_27/phase25/blind_protocol.json) fixes ten seeds before labels.

The repository's Phase24 measurements are explicitly synthetic. Its geometry, charge and orbital-gap features include proxies, and its tetrahedral `CIP_R_minus_S` feature does not encode BINOL axial chirality. Reuse its screening architecture only; add a validated axial parity representation and fit preprocessing on training data only. Enforce the expected sign reversal of a descriptor predictor under catalyst reflection.

Compare ensemble-network, lowest-conformer-network and Phase24-style descriptor models against the same external labels. Evaluate major-configuration accuracy **with abstention rate**, signed effective-ΔΔG MAE, 95% coverage and interval width by holdout family. Experimental ee at ±100% should be treated with detection limits, not transformed to infinite energies. Account for measurement uncertainty and condition matching.

For random and diversity acquisitions, count failed jobs, conformer searches, gradients, frequencies, IRC and reference calculations in cumulative electronic-structure cost. Use a common calibrated hardware basis; report GPU/CPU time separately and at identical budget checkpoints. Do not compare equal numbers of molecules when their costs differ. Ten frozen seeds alone do not constitute an executed active-learning experiment. External label custody, locked prediction hashes, labels, acquisition costs and model evaluations are all still absent.

## Acceptance and next execution boundary

| User criterion | Current status | Missing evidence |
|---|---|---|
| 25.1 fixed chemistry | Partial | Certified axial/product CIP, complete reactive maps, literature conditions |
| 25.2 competition network | Topology only | Computed paths, aggregation/inhibition, exchange timescales |
| 25.3 ensembles and method checks | Starts/pilot only | Complex/TS ensembles, correlated benchmarks, truncation errors |
| 25.4 validated barriers | Not run | Same-surface opt/frequency/IRC, solution free energies, two solvent PMF comparisons |
| 25.5 physical kinetics | Not run | Calibrated network, uncertainty and physical controls |
| 25.6 true external blind test | Split only | Unseen labels and cost-matched three-model/AL results |
| 25.7 convergence and validation | Not passed | ≤0.5 kcal/mol expansion change, repeat CIs, external configuration/error/coverage |

The immediate resource boundary is production electronic-structure/search/solvent capacity and independent experimental labels. OpenAI model-inference GPUs are not exposed to this session as a scientific job scheduler. The local implementation is useful groundwork, not completion of Phase25.
