# Additional published-data calculations for Phase 28

This extension answers the request for more computation with **new analyses of published numeric data**, not more synthetic draws. The original `results/run_20260920/` remains unchanged. New inputs, frozen settings, code and outputs are separate.

## Inputs and statistical scope

The [Nature Synthesis 2026 Figure 2 source](https://www.nature.com/articles/s44160-026-01039-y/figures/2) provides 50 timepoint mean/SD records for cumulative Cl2 and Faradaic efficiency, plus 777 polarization coordinate pairs used here. Original sheet/row/column and workbook SHA256 are retained. The caption and SI printed p1/PDF p4 specify mean±SD from four independent experiments. Individual replicate values and cross-time pairing are unavailable.

The extension calculates 50 pointwise mean intervals, 48 time-contrast/correlation scenarios, 28 common-current interpolations and 72 startup/window sensitivity combinations. The original stability traces are not treated as independent electrode replicates.

## Small-sample intervals and missing covariance

Intervals use mean±t(0.975,3) SD/sqrt(4), with critical value3.18245. They assume independent, approximately normal replicates within each timepoint; the unavailable raw values prevent checking that assumption. These are **pointwise intervals, not a simultaneous time-series confidence band**.

At24h, cumulative Cl2 is0.338190mol with conditional95% mean interval[0.326924,0.349456]mol. The corresponding average rate is0.01409125mol/h,[0.01362184,0.01456066]. FE is99.79236%,[99.76278,99.82194]%. A through-origin fit to the cumulative means gives0.01391599mol/h and maximum residual0.00651228mol; this is descriptive, with no independence-based slope significance test. The zero-amount origin is fixed, not an inferred zero-variance population.

The24h−0h FE difference is−0.06296 percentage points. An unverified independent-groups Welch calculation gives[−0.11799,−0.00793]. If the same four runs were followed through time, the paired standard error instead depends on unknown correlationrho: sqrt((s0²+s24²−2rho s0s24)/4). Assumed pairedrho=0 gives[−0.12852,+0.00260],rho=−1 gives[−0.15105,+0.02513], andrho=0.5 gives[−0.11363,−0.01229]. Whether zero is excluded therefore depends on missing covariance. No definitive degradation significance is claimed. The24h−1h difference is only+0.00240 percentage points; failure to reject zero would still not establish equivalence without a prespecified equivalence margin.

## Polarization at common currents

Strictly increasing positive-current data are interpolated at10,50,100,300,500,800 and1000mA/cm². Extrapolation is forbidden. Linear versus PCHIP interpolation differs by at most0.04372mV; this measures interpolation sensitivity, not experimental accuracy.

The200h−100h coordinate difference is+25.480mV at100mA/cm²,−0.311mV at500mA/cm² and−5.209mV at1000mA/cm². The curves cross; the high-current negative difference cannot be generalized into an overall activity gain. Their durability iR processing remains unresolved, and the values are anodic potentials versus NHE, never whole-cell energy.

The2h source table/legend does not provide a separate numerical0h trace. Figure2a is therefore only a cross-panel proxy baseline, with unmatched electrode identity and correction status. At1000mA/cm², the200h-minus-proxy difference is−28.408mV against the iR-corrected curve or−71.995mV against the uncorrected curve. Both are disclosed as sensitivity cases, not a matched0/100/200h degradation series.

## Startup and endpoint-window robustness

Each repeated timestamp is collapsed to its ordinate median before applying six startup exclusions(0,1,5,10,20,30h after the first recorded time) and four endpoint-window widths(1,3,6,12h). All72 combinations are reported, avoiding selection of a favorable window.

The Nature2023 constant-current potential change spans−4 to+44.25mV across the full grid, and−4 to+9.25mV after excluding at least10h. Nature Synthesis2026 spans+3 to+17mV overall and+3 to+6mV after that exclusion. The2023 constant-potential current-ordinate relative change ranges from−2.686% to+0.762%; its normalization remains unresolved.

The original+44/+17mV endpoint summaries remain correct for their defined windows. The extension shows that they are strongly influenced by the early-time window and should not be interpreted as uniform degradation rates. The10h split is an analyst-defined grid summary, not an established physical steady-state boundary. Timing alone does not establish the physical cause of the initial drift. The sensitivity grid does not prove absence of degradation and supplies no electrode-population confidence interval.

## Reproduction and status

Run `python projects/phase28/run.py --self-test` and `python projects/phase28/code/additional_public_analysis.py --out NEW_DIRECTORY`. All18checks(11original+7new) pass. Three new plots were opened and visually reviewed. Exact input/output hashes are recorded in [input_output_manifest.json](../results/additional_public_analysis/input_output_manifest.json).

No new experiment, validated loss mechanism or industrial lifetime benefit is claimed. The contribution is additional source-grounded numerical evidence and a sharper account of statistical, startup and normalization limits.
