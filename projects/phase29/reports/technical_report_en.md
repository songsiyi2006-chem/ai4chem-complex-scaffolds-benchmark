# Phase 29 technical report: public evidence and bounded local computation

This page records the initial delivery at commit `73b1ea0`. The [xTB matrix addendum](xtb_matrix_en.md) and [DFT crosscheck](dft_crosscheck_en.md) document subsequent calculations; initial counts and results below are retained for traceability.

Date: 2026-09-20. No new experiments, HPC calculations, or industrial validation are available.

## Outcome

The deliverable establishes reproducible source/identity audits, public-data arithmetic checks, molecular structure calculations, a semiempirical method-sensitivity diagnostic, a kinetic identifiability model, and segregated synthetic learning/manufacturing workflows. It does **not** demonstrate a same-substrate/same-partner C2/C4 switch, a validated mechanism, prospective chemical generalization, or manufacturing superiority.

The central correction is structural. The JACS model uses **4-phenylpyridine and benzaldehyde**, giving a secondary benzylic alcohol; C4 is already occupied. The Nature Communications model uses **2-phenylpyridine and 4-methoxyacetophenone**, giving a tertiary alcohol; the starting substrate's C2 is already occupied. These examples cannot constitute a matched switching pair. “Hydroxymethylation” must not be uniformly interpreted as installing –CH2OH.

## Sources and numerical checks

[JACS 2026](https://doi.org/10.1021/jacs.6c07265) has a publicly retrieved SI. Scheme S1 and the C1 characterization support the structure assignment. Conflicts between a 4-bromo drawing and 2-bromo preparation prose, and between 5 mg preparation loading and 2.5 mg TON discussion, remain explicit. An EXAFS coordination number of 4.3±0.2 and Zn–N distance of 2.00±0.01 Å support an average local environment, not an atom-resolved operating interface.

The [Nature Communications competitor](https://www.nature.com/articles/s41467-026-71858-2) provides public text, SI and coordinates. Its source spreadsheet contains coordinate records, not a reaction-yield dataset. Its Table S1 C1/C1′ footnote discrepancy is preserved. Product renumbering can call an alpha product C2 even when attachment occurred at C6 in the fixed starting-substrate numbering; the audit retains fixed parent maps.

Forty published condition records and fourteen thermochemical rows were transcribed with source locations. Trace and N.R. remain unquantified, not numeric zeros. No missing error bars or replicates were invented, and these sparse/confounded literature rows were not used for a yield-prediction model.

Recalculating common-reference Gibbs-energy differences from NC Table S3 gives TS1/TS2 barriers of **9.72295/12.03394 kcal/mol**, a 2.31099 kcal/mol difference. Subsequent barriers are 15.87392/21.14531 kcal/mol. `Gcorr+E=G` is consistent within tabulated precision. This is arithmetic reproduction of published DFT data, not a new optimization, frequency calculation, IRC, or interfacial calculation. No final product ratio was obtained by directly exponentiating a single barrier difference.

The NC model consumes 51.821 F/mol isolated product, corresponding to **3.859% FE at two electrons**, consistent with the reported 3.86%. JACS conventions give **8.876% at two electrons and 17.752% at four electrons**. Electron numbers and limiting-substrate versus product denominators are kept separate. NC standard and 3 mmol scale-up solution-volume STYs are approximately 0.550 and 0.237 g/L/h. The scale-up uses 80 h electrolysis and preparative TLC; it is not continuous-manufacturing validation. A JACS estimate using constant 3.01 V gives approximately 6.96 kWh/kg at the cell-only boundary; a matched NC voltage and full process inventory are missing.

## Executed molecular calculations

The computational panel contains twelve neutral pyridine structures spanning parent, alkyl, alkoxy, halogen, strongly withdrawing, carbonyl and aryl substituents. The common proposed partner is 4-methoxyacetophenone. Published scope, procurement, salt/speciation and analytical feasibility are not assumed for unverified members.

RDKit descriptor and ETKDGv3/MMFF94s workflows generated **48 structure records and 151 conformer attempts**, of which 150 converged. Each record retains its lowest converged conformer; the failed attempt remains visible. Thirty-six heavy-atom-mapped proposals conserve formula and map identity. They are not successful reactions. The parent pyridine's C2/C6 proposals are symmetry equivalents, so records are not counted as independent products. New carbinol stereochemistry is unspecified, and force-field energies only compare conformers of the same species.

A separate [quantum diagnostic](../results/quantum_preflight/README.md) executed twelve GFN1/GFN2 single points on N-protonated pyridine, 3-methoxy and 3-cyano analogues. Each method compared +1 closed-shell and neutral radical states at the same geometry in ALPB MeCN. All twelve logs report SCC convergence and normal termination, but also IEEE divide-by-zero/denormal flags, which are retained as numerical warnings. State-energy differences differ by approximately **0.736–1.021 eV between methods**. Independently equilibrated ALPB states do not establish a rigorous nonequilibrium vertical electron affinity, redox potential, free energy, or regioselectivity prediction.

## Model and learning boundaries

Five hypothetical ODE scenarios conserve pyridine equivalents and include relay state, surface occupancy, three product branches and secondary losses. Rates are assumed, not fitted. Preferential C4 destruction can produce s2≈0.919 while total regioisomer yield falls to ≈0.179, showing why a high site fraction alone is not useful selectivity. The 120-draw interval of approximately 0.283–0.890 reflects assumed lognormal rate uncertainty, not a posterior inferred from chemistry.

A central-difference sensitivity audit checks six rate parameters. Two endpoint observations have rank two; time points under one condition retain channel degeneracy. Adding a support intervention improves identifiability in this model. These diagnostics justify discriminating experiments, not mechanistic conclusions.

The learning benchmark is exclusively synthetic: 120 independent artificial families split 66/24/30. Ridge MAE is approximately 0.0457 versus 0.1546 for a constant baseline, 0.0717 for a shallow tree and 0.0598 for fixed-kernel GP. Synthetic conformal coverage is 90% over 30 test units. Hyperparameters are fixed; preprocessing uses training data only. An insufficient calibration sample produces an unbounded interval/abstention. None of these numbers measures chemical prediction performance.

Equal-cost acquisition curves evaluate random, diversity-per-cost and a distance-plus-site-prediction proxy. The full yield/risk/cost/information utility is only implemented and unit-tested; it was **not** the policy used by the budget benchmark. Real deployment always abstains because real independent labels are unavailable. No experimental savings are claimed.

## Deliverable acceptance

G0 is partial: key source structures are audited, but a common experimental domain and all panel procurement/analysis requirements are unresolved. G1 covers software, table checks and bounded molecular diagnostics. G2–G5 remain unfulfilled. Sixteen tests check identities, fixed-site eligibility, symmetry, leakage, mass balance, pulse integration, abstention, conformal edge cases and manufacturing units. See the [acceptance ledger](../results/acceptance.json), [data card](data_card.md), [mechanism matrix](mechanism_discrimination.md) and [industry requirements](industry_requirements.md).

The immediate research contribution is an auditable negative boundary and an executable discriminating design. It does not warrant a novelty-first or publication-readiness claim.
