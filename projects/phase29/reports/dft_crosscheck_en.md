# Phase 29: fixed-geometry gas-phase PBE0 crosscheck

[Project](../README.md) · [中文](dft_crosscheck_zh.md) · [xTB matrix](xtb_matrix_en.md) · [Raw evidence](../results/dft_crosscheck/)

Eight Psi4 1.11 single points were executed on 2026-09-20; all eight SCFs converged. PBE0 is an independent approximate model, not experimental ground truth. These isolated molecules contain neither a catalytic surface nor a carbonyl reaction partner and cannot establish C2/C4 selectivity.

## Matched calculation

Q01 is hypothetical N-protonated pyridine, Q02 its 3-methoxy derivative and Q03 its 3-cyano derivative (matrix P01/P04/P09). Their original MMFF cation XYZ files were copied byte-for-byte. Each pair holds nuclei fixed and compares a +1 singlet RKS cation with a neutral doublet UKS radical: ΔE = E(radical) − E(cation).

All three molecules received PBE0/def2-SVP calculations; Q01 alone also received def2-TZVP. Density-fitted SCF used a 75×302 DFT grid, energy/density thresholds 1e-9/1e-7 and 150 maximum iterations. Jobs ran sequentially with two threads, 900 MB and a 300 s timeout each; total elapsed job time was approximately 396 s. [Configuration](../configs/dft_crosscheck.json), [execution manifest](../results/dft_crosscheck/manifest.json), [Psi4 DFT documentation](https://psi4.github.io/psi4docs/master/dft.html) and [SCF documentation](https://psi4.github.io/psi4docs/master/scf.html) make the settings inspectable.

| Molecule | Basis | ΔE / eV | Radical ⟨S²⟩ |
|---|---|---:|---:|
| Q01 | def2-SVP | −4.957325 | 0.776918 |
| Q02 | def2-SVP | −4.776458 | 0.770825 |
| Q03 | def2-SVP | −5.689793 | 0.775537 |
| Q01 | def2-TZVP | −4.997139 | 0.775941 |

Doublet ⟨S²⟩ exceeds the ideal 0.75 by approximately 0.021–0.027. It was evaluated from occupied alpha/beta orbitals and the AO overlap matrix; the formula and electron counts are retained in each result JSON. Tiny negative RKS values near zero are roundoff. Wavefunction stability was not checked. The Q01 basis change is −0.039814 eV; two non-diffuse bases on one molecule do not demonstrate basis convergence.

## Electron-reference caveat and centered contrasts

Raw GFN1 and GFN2 ΔE values differ from PBE0/def2-SVP by −5.86 to −5.59 and −4.89 to −4.81 eV, respectively. These are **unaligned charge-state reference offsets, not predictive-error estimates**. The official [tblite GFN2 tutorial](https://tblite.readthedocs.io/en/latest/tutorial/python/singlepoint.html) explains an empirical free-electron self-interaction correction for ionization differences; [xTB documentation](https://xtb-docs.readthedocs.io/en/latest/sp.html) separately describes IPEA reparameterization. We did not run `--vipea` or apply an unexplained universal shift across methods.

Within-method double differences ΔE(substituted)−ΔE(Q01) cancel a common additive electron-reference constant. Only matching gas-phase geometries enter this comparison:

| Centered effect / eV | PBE0/def2-SVP | GFN1 | GFN2 |
|---|---:|---:|---:|
| Q02 − Q01 | +0.180867 | +0.166602 | +0.162533 |
| Q03 − Q01 | −0.732468 | −0.475732 | −0.678319 |

Directions agree for both substituents. The cyano contrast differs from PBE0 by +0.256736 eV for GFN1 and +0.054149 eV for GFN2. Two contrasts do not establish general method accuracy or agreement with experiment. Centering removes constant offsets, not all model bias.

![Raw and centered gas-phase comparison](../results/dft_crosscheck/comparison/matched_model_comparison.png)

## Reproduction and limits

Use the existing chemistry Python and DLL environment documented in the project README; output directories must be new:

```powershell
python projects/phase29/code/dft_crosscheck.py --python C:/Users/HUIWEI/miniconda3/envs/phase7/python.exe --out work/runs/phase29-dft-new
python projects/phase29/code/compare_quantum_models.py --dft work/runs/phase29-dft-new --xtb-csv projects/phase29/results/xtb_matrix/state_differences.csv --out work/runs/phase29-dft-comparison-new
```

The outer Python is the chemistry environment; `--python` selects an existing Psi4 environment. No installation occurs. `--analyze NEW_COPY_OF_RESULT_DIR` regenerates derived energy tables without SCF execution; keep frozen evidence unchanged. Some native Psi4 system timing fields are nan; parent-process wall times are recorded separately and all energies are finite.

Standard GFN1/GFN2 open-shell treatment is not equivalent to UKS: `--uhf 1` sets the unpaired occupation and does not enable the separate `--spinpol --tblite` model, which was not used. See the [xTB spin-polarization documentation](https://xtb-docs.readthedocs.io/en/latest/spgfn.html). Matching electron counts and geometry does not match every spin approximation.

No DFT geometry optimization, conformer ensemble, diffuse-basis study, stability analysis, solvent/counterion, constant-potential interface, thermal free energy, transition state or IRC was evaluated. Experimental potentials, switching, catalytic mechanisms and manufacturing advantages remain unestablished. [Acceptance gates](../results/acceptance.json) are not upgraded by convergence.
