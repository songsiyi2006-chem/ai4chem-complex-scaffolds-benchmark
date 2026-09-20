# Phase 28 technical report

Executed 2026-09-20 within the user's scope: **public data and bounded local computation; no new experiments or HPC**. The deliverable comprises traceable numerical reanalysis, structural-identifiability diagnostics, causal control baselines and explicitly uncalibrated cost scenarios. It does not establish a new electrode mechanism, industrial lifetime or commercial advantage.

## Published numerical evidence

The primary sources are the [Nature 2023 paper](https://www.nature.com/articles/s41586-023-05886-z) and [Nature Synthesis 2026 paper](https://www.nature.com/articles/s44160-026-01039-y). Their Figure 2 workbooks and relevant captions/SI passages were accessed. This is a bounded audit, not a page-by-page validation of every spectrum or an experimental reproduction. These are collaborative papers involving Guangxi Normal University and other institutions, not independent accomplishments of a single group. The 2026 squaramides must not be conflated with the 2023 molecular family.

Three numeric traces contain 196,496 source rows. Every pair and its Excel row are preserved in compressed CSV, with source URLs and SHA256. The 2023 constant-potential series ends at 120 h; its first/last-hour current-ordinate medians are 5.8175/5.6625 (−2.664%). The figure's current normalization remains unresolved, so no absolute charge is calculated. The 2023 800 mA/cm² series has median anode potentials 1.657/1.701 V versus NHE (+44 mV). The 2026 1000 mA/cm² series ends at 200.0566 h and has median potentials 1.406/1.423 V (+17 mV). These are descriptive windows, not newly measured lifetimes or a head-to-head comparison. Trace-specific iR handling remains unverified.

Repeated timestamps number 29,866, 11,043 and 115,861; unique times number 13,246, 7,707 and 18,773. The source does not establish why the repetitions occur. The sensitivity check first collapses each timestamp to its median, then evaluates endpoint windows, obtaining −2.686%, +44.25 mV and +17 mV. The plotted hourly 2.5–97.5% band describes within-trace records; it is not uncertainty across independent electrodes. No row-wise iid confidence interval, hazard distribution or lifetime extrapolation is reported.

In Nature 2023 Figure 2c, stored means cannot be reproduced by direct arithmetic over adjacent numeric columns. Because the caption describes ODC-based extrapolation, this is an unresolved transformation, not a demonstrated publication error. Published anode potentials are not integrated as full-cell power. NaOH-normalized published energy and Cl2-normalized local scenarios are not placed on a shared ranking. Nature Synthesis 2026 technoeconomic assumptions are recognized as scenarios, not quotations or measurements from a customer plant.

## Conservation and identifiability

The model tracks anchored inventory A, irreversible products D, soluble inventory S and removed material W with total inventory one. Piecewise constant-load intervals are propagated analytically; the maximum conservation error is 1.11e-16. All kinetic constants, voltage mappings, inventory calibration and fixed clearance are hypothetical.

With Q2=integral(load² dt), A=exp[-(k_ir+k_des)Q2]. Voltage observes only the sum, so the voltage sensitivity matrix has rank 1 and the partition profile has exactly zero chi-square range. Adding calibrated dissolved inventory and surface observations gives rank 2 for this two-parameter model when clearance and calibration are known. Adding A alone would not resolve the loss partition. If calibration or clearance becomes unknown, identifiability must be reassessed.

Four simplified candidates represent fixed-site loss, soluble mediation, carrier damage and reversible transport lag. Fixed-site and carrier explanations are exactly equivalent for voltage alone. Joint synthetic observations distinguish them within the specified model family. Constant-activity and empirical cumulative-load regression baselines were also fitted. No unsupported statistical residual correction was added.

Nine synthetic batches contain 72 points each, split by entire batch into four training, two validation and three blind challenges with start/stop and load changes. Training-only parameter estimates are k_ir=0.00288944 h^-1 and k_des=0.00151019 h^-1. Blind voltage RMSE values are 14.53, 8.89 and 6.25 mV, below the frozen 15 mV software benchmark threshold. No published trace calibrates these values. Passing the synthetic test does not validate real electrodes. No Bayesian posterior was sampled; sensitivity intervals are not posterior credible intervals, and the profile threshold is illustrative under the imposed error assumptions.

## Causal, equal-product comparison

Constant current, a simple price rule and four-step finite-grid MPC operate for 48 h with identical charge (1440 Ah), assumed qualified Cl2 (1.87610 kg), fixed FE (0.985) and assumed purity (0.998). All obey the fixed load, ramp and storage bounds, ending with zero inventory deviation. The shared initial physical buffer is two rated production hours, capacity is four hours, and physical stock equals two hours plus the reported deviation; negative deviations do not imply negative stock. MPC sees only current/past signals; future prices use persistence/trailing-mean forecasts. Guaranteed baseload availability is a model-domain requirement; the runner rejects lower availability.

| Synthetic strategy | Cell plus auxiliary kWh/t Cl2 | Scenario USD/t Cl2 | Model-predicted lifetime h |
|---|---:|---:|---:|
| Constant | 2611.29 | 442.43 | 79.26 |
| Rule | 2614.53 | 438.72 | 76.09 |
| MPC | 2611.92 | 444.01 | 77.05 |

The simple rule reduces nominal scenario cost by 0.84%; MPC increases it by 0.36%. The specified model has no recovery term. Jensen's inequality therefore predicts greater loss under varying load at the same mean current: waveform-specific lifetime benefit has not been built in. The negative MPC result is retained rather than tuned away.

Each strategy has 2000 matched assumption draws. Draws crossing the frozen A=0.7 failure threshold within the horizon (125 rule, 117 MPC) are retained as infeasible and excluded from equal-qualified-output benefit summaries. Conditional 2.5/50/97.5-percentile cost reductions are −0.445/+0.605/+1.420% for the rule and −0.937/−0.436/−0.062% for MPC. No feasible draw reaches the preregistered 5% decision threshold. These are scenario sensitivities, not a probability of industrial success.

Costs include electricity, auxiliary use, consumed layer life, support/membrane, maintenance, disposal and off-spec reserves. Layer use is prorated by consumed life, retaining residual value. No instantaneous zero-downtime replacement is invented. Real plant capital, financing, salt/water and separation inventories remain absent. Environmental output covers electricity only, not a complete LCA. Break-even layer-price/lifetime curves are provided; loading optimization is not claimed without loading data.

## Acceptance and reproduction

Eleven meaningful software checks pass, including conservation, Faraday units, right censoring, batch isolation, causal scheduling, fair output, constraints and unsupported-configuration rejection. Python 3.12.14, NumPy 2.4.6, SciPy 1.18.0 and Matplotlib 3.11.1 were reused with two CPU threads. No dependencies were installed. Run `python projects/phase28/run.py --self-test`, then use `--out NEW_DIRECTORY` to reproduce all local outputs from the extracted numerical inputs.

G0 is partial at the working-interface/measurement-convention level; G1 passes software conservation; G2 and G3 are synthetic only; G4 produces scenarios that do not meet the 5% target; G5 does not establish a real mechanism or industrial benefit. See [acceptance.json](../results/acceptance.json), [mechanism discrimination](mechanism_discrimination.md), [novelty matrix](novelty_matrix.md), [industry requirements](industry_requirements.md) and [resource plan](resource_plan.md).
