# Phase 29 addendum: 144 executed molecular sensitivity single points

Date: 2026-09-20. The extension executes **144 xTB single points** on twelve previously proposed substrates: two charge/spin states, two methods and three environments. All jobs satisfy SCC convergence and normal termination, while all also report IEEE floating-point flags. These remain flagged molecular diagnostics, not validated chemistry predictions.

The unaligned raw GFN1/GFN2 fixed-geometry state differences disagree by approximately **0.58–1.07 eV**, and some computed ordering reverses. The raw gap may include method-specific charge-state reference offsets and must not be treated entirely as predictive error. Ordering changes demonstrate that the molecular dependence is not identical across methods.

## Frozen design and meaning

The [configuration](../configs/xtb_matrix.json) was frozen at 12:50:39 UTC before execution. Ring N protonation is assumed for every substrate; alternative protonation, salt forms and experimental speciation are not established. States are the +1 closed-shell cation (UHF 0) and neutral radical (UHF 1), with explicit electron-parity checks.

P01/P04/P09 reuse byte-identical XYZ files from original Q01/Q02/Q03. Both SHA256 and canonical cation identity are checked. Nine additional single cation geometries use seeded ETKDGv3/MMFF94. Every charge state, method and environment uses the same nuclei for a given species. No geometry optimization at xTB, conformer ensemble, frequency or reaction-path calculation was performed.

The matrix uses GFN1/GFN2 in gas, ALPB acetonitrile and ALPB DMF, sequentially with two CPU threads, accuracy 0.1, 300 maximum SCC iterations and a 90-second per-job timeout. Define `Δ_model = E_model(neutral radical) − E_model(cation)`, converted using 27.211386245988 eV/Eh.

Gas outputs are the electronic model total-state energies. ALPB totals additionally include model solvation terms, as documented by the [official xTB documentation](https://xtb-docs.readthedocs.io/en/latest/gbsa.html). The default `gsolv` reference is retained without an extra 1 bar-to-1 M correction. Independently equilibrated charge-state continua do not provide a rigorous nonequilibrium-solvent vertical electron affinity. These outputs are neither full molecular Gibbs free energies nor redox potentials or activation barriers.

**Electron-reference limitation:** the [official tblite GFN2 tutorial](https://tblite.readthedocs.io/en/latest/tutorial/python/singlepoint.html) describes an empirical 4.846 V per-electron correction for ionization potentials obtained from charge-state energies, addressing free-electron self-interaction. The [xTB single-point documentation](https://xtb-docs.readthedocs.io/en/latest/sp.html) also distinguishes dedicated IP/EA workflows. No empirical shift is applied here; the GFN2 example's number is not assumed universal for GFN1 or every state in this study. Raw cross-method offsets alone do not measure predictive accuracy. A method-specific constant reference offset cancels from within-method rankings, substrate contrasts `Δ(Pi)−Δ(P01)`, and same-method solvent differences if it is common to those molecules/environments. Molecular dependence still requires independent comparison.

## Results and reproducibility

All 144 jobs yield **72 complete comparable state pairs**. Six MeCN differences on the three reused original geometries reproduce the previous preflight exactly at stored precision (maximum absolute discrepancy 0 eV). This is numerical reproducibility, not independent physical validation: twelve of the new single points rerun the previous configuration and must not be counted as independent molecules.

| Environment | Absolute unaligned raw GFN1/GFN2 difference / eV | Spearman ρ between methods | Reversed ordering among 66 species pairs |
|---|---:|---:|---:|
| Gas | 0.57998–1.00864 | 0.91608 | 6 |
| ALPB MeCN | 0.66915–1.04784 | 0.95105 | 4 |
| ALPB DMF | 0.70697–1.07360 | 0.95105 | 4 |

Within a method, the ALPB-minus-gas shift in the state difference spans **+1.72260 to +2.11740 eV**. This is model environmental sensitivity, not an observed solvent shift in electrode potential. Gas-to-MeCN comparison reverses 7/66 species pairs for each method. Within each method, MeCN and DMF orderings agree even though numerical values differ. The full [analysis](../results/xtb_matrix/analysis.json) retains all fifteen pairwise method/environment ordering comparisons; their ordering has no yield, activity or C2/C4 meaning.

![Executed sensitivity matrix](../results/xtb_matrix/method_environment_matrix.png)

All jobs retain stdout/stderr, commands, state identities, geometry and log hashes. None timed out or failed convergence. All 144 emit `IEEE_DIVIDE_BY_ZERO`; subsets additionally emit `IEEE_DENORMAL` or `IEEE_UNDERFLOW_FLAG`. The cause of these IEEE flags reported by the executable was not resolved, and convergence is not used to hide them.

## Interpretation and limits

Reference offsets must be distinguished from molecular/environmental dependence before independent calibration or electrochemical interpretation. The raw gaps are neither experimental error bars nor a direct measure of failure to predict chemical effects. There are no counterions, explicit solvent, acid/salt concentrations, catalyst surfaces, double layers, constant-potential treatment, coupling partners, rearomatization paths, secondary product losses or experimental labels. No reaction outcome is assigned from atomic charges.

A DFT comparison must match species, nuclei, charge, spin and gas environment. Gas DFT cannot directly validate independently equilibrated ALPB differences. A small-basis DFT calculation is itself an approximation, not experimental truth. The original three gas inputs are available for the independent DFT crosscheck. Scientific gates G2–G5 remain unmet.

## Files and rerun

See [state differences](../results/xtb_matrix/state_differences.csv), [method sensitivity](../results/xtb_matrix/method_sensitivity.csv), [environment sensitivity](../results/xtb_matrix/environment_sensitivity.csv), [execution manifest](../results/xtb_matrix/manifest.json) and [file hashes](../results/xtb_matrix/file_manifest.json).

```text
python projects/phase29/code/xtb_matrix.py --xtb <existing-xtb-executable> --out work/runs/phase29-xtb-matrix-new
```

Use the existing scientific environment and its xTB library path; nothing is installed. The output directory must be new. Eight additional tests pass and cover frozen scope, exact original geometry hashes, retained warnings, failed-SCC rejection, signs/units, pair comparability and ordering ties/reversals. See the final integrated validation record for the total including other calculation modules. The plotted figure was opened and visually checked for readable labels and clear environment separation.
