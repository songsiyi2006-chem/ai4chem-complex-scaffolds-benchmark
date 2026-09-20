# Molecular xTB preflight

EXECUTED_MOLECULAR_PREFLIGHT: hypothetical N-protonated solution species, not verified catalytic intermediates. Single MMFF94 geometry per species; no conformer ensemble, surface, electrolyte, counterion, explicit solvent, constant potential, thermal correction, or kinetic calculation. ALPB MeCN is independently equilibrated for each state. Method sensitivity is not experimental uncertainty. No redox-potential or C2/C4 prediction.

ΔE = E(neutral radical) − E(pyridinium cation), at identical nuclei within each method.
Negative values indicate a lower total electronic energy under this model; the free electron reference,
electrode reference and nonequilibrium solvent response needed for electrochemical interpretation are absent.

| Molecule | GFN1 ΔE / eV | GFN2 ΔE / eV | Absolute method difference / eV |
|---|---:|---:|---:|
| pyridinium | -8.8086 | -7.7924 | 1.0162 |
| 3-methoxypyridinium | -8.7787 | -7.7580 | 1.0207 |
| 3-cyanopyridinium | -9.2395 | -8.5037 | 0.7358 |

Raw input geometries, engine stdout/stderr, charge/spin settings and hashes are retained.
SCC convergence and normal termination: 12/12 jobs.
Numerical warnings: 12/12 jobs. This Windows build reports IEEE floating-point flags despite SCC convergence; the warnings are retained in each job manifest and stderr. These results are flagged diagnostics, not validated electrochemistry.
This is an optional diagnostic independent of the reaction-panel ranking and synthetic models.
