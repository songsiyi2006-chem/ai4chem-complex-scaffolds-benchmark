# AI4Chem — Complex Scaffolds Benchmark

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![RDKit](https://img.shields.io/badge/RDKit-2024.03%2B-38B2A3)
![License](https://img.shields.io/badge/License-MIT-yellow)
![Docs](https://img.shields.io/badge/Docs-EN%20%7C%20ZH-informational)
![ForceFields](https://img.shields.io/badge/Conformers-ETKDGv3%20%2B%20MMFF94%2BUFF-orange)
![MD](https://img.shields.io/badge/OpenMM-8.6%20GBSA%2FOBC2-9B59B6)

## Mission

**Phase 1 — Conformer & descriptor suite.** Stress-test frontier, fully open-source 3D conformer pipelines (**RDKit ETKDGv3 + MMFF94/UFF**) against **10 structurally complex, synthetically plausible, unindexed molecular entities** — macrocycles, cubane-type strain cages, beyond-Ro5 PROTAC prototypes, atropisomeric biaryls, perfluorinated cages, zwitterionic B–N dative systems, covalent warheads, helicenes, and sterically jammed peptoids — and quantify their physicochemical profiles, conformational energy landscapes, intramolecular H-bond networks, and **GNN featurization readiness** (PyG-loadable atom/bond tensors).

One target (M09) was caught by the pipeline with an **unkekulizable SMILES** — a genuine input-integrity failure — and a programmatically repaired aza-[5]-helicene reference (**M09R**) completes the case study.

**Phase 2 — Heavy dynamics suite.** Five next-generation med-chem targets (CRBN molecular glue, ADC Val-Cit-PAB linker, bicyclic disulfide peptidomimetic, allosteric covalent inhibitor core, VHL-mimetic hybrid) advance into **relaxed torsional barrier scanning** (36-point MMFF94 scans, ΔE‡ = 7.8–72.8 kcal/mol) and **production molecular dynamics** (OpenMM 8.6, Sage 2.1.0 valence + MMFF94 charges, GBSA/OBC2, 300 K, 200 ps) — labeling the fourth dimension that 2D GNNs cannot see. See [`DYNAMICS_REPORT_EN.md`](./DYNAMICS_REPORT_EN.md) / [`DYNAMICS_REPORT_ZH.md`](./DYNAMICS_REPORT_ZH.md).

## Quick Navigation

| Asset | Link |
|---|---|
| 🇬🇧 English technical report (Phase 1) | [`BENCHMARK_REPORT_EN.md`](./BENCHMARK_REPORT_EN.md) |
| 🇨🇳 中文技术报告（第一阶段） | [`BENCHMARK_REPORT_ZH.md`](./BENCHMARK_REPORT_ZH.md) |
| 🇬🇧 Dynamics report (Phase 2) | [`DYNAMICS_REPORT_EN.md`](./DYNAMICS_REPORT_EN.md) |
| 🇨🇳 动力学报告（第二阶段） | [`DYNAMICS_REPORT_ZH.md`](./DYNAMICS_REPORT_ZH.md) |
| Benchmark pipeline script (Phase 1) | [`molecule_benchmark.py`](./molecule_benchmark.py) |
| Figure generation script (Phase 1) | [`generate_assets.py`](./generate_assets.py) |
| Heavy dynamics pipeline (Phase 2) | [`run_heavy_dynamics_benchmark.py`](./run_heavy_dynamics_benchmark.py) |
| OpenFF→OpenMM parameterizer | [`export_openff_system.py`](./export_openff_system.py) |
| Workspace bootstrap script | [`setup_and_download.py`](./setup_and_download.py) |
| Machine-readable results (Phase 1 / 2) | [`bench_results/benchmark_results.json`](./bench_results/benchmark_results.json) / [`results_phase2/phase2_results.json`](./results_phase2/phase2_results.json) |
| Auto-generated run report | [`bench_results/benchmark_report.md`](./bench_results/benchmark_report.md) |
| 3D ensembles & minima | `bench_results/sdf/*_ensemble.sdf`, `bench_results/sdf/*_min.sdf` |
| PyG feature tensors | `bench_results/features/*.npz` |
| Phase 2 scans / trajectory | `results_phase2/T0*_torsion_scan.csv`, `results_phase2/T02_trajectory.dcd` |

## Figure Previews

### Fig. 1 — Structural Grid (10 targets + repaired reference)

![Molecular grid](./figures/fig1_molecular_grid.png)

### Fig. 2 — Chemical Space vs. Conventional Drug-Like Space

![Chemical space](./figures/fig2_chemical_space.png)

### Fig. 3 — Structural Complexity Radar

![Complexity radar](./figures/fig3_radar_complexity.png)

## Benchmark Summary Table

All values computed with RDKit 2026.03.5 (50 conformers/molecule, ETKDGv3 + MMFF94/UFF, all optimizations converged). E in kcal/mol; ΔE_ens = E_max − E_min over the optimized ensemble; Stereo a/u = assigned/unassigned chiral centers; d_min = shortest non-bonded heavy-atom distance (≥ 4 bonds apart) in the E_min conformer; IMHB = intramolecular H-bonds in the E_min conformer.

| # | Molecule | MW (Da) | cLogP | TPSA (Å²) | Fsp³ | RotB | Stereo a/u | MaxRing | FF | E_min | ΔE_ens | IMHB | d_min (Å) | Risk flags |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| M01 | Macrocyclic chameleon peptide | 529.6 | −0.72 | 165.2 | 0.538 | 5 | 4/0 | 13 | MMFF94 | 67.4 | **43.1** | 2 | 2.73 | MAC, ENTROPY |
| M02 | Azaspiro-cubane bioisostere | 366.3 | 1.78 | 32.8 | 0.611 | 2 | 0/**5** | 6 | MMFF94 | 106.3 | 8.6 | 0 | 2.78 | STRAIN, UST |
| M03 | PROTAC prototype (2 fragments) | 815.9 | 2.51 | 241.4 | 0.282 | 16 | 0/1 | 6 | MMFF94 | 29.2† | 12.2† | 1 | 2.75† | BIG, FLEX, FRAG, ENTROPY, UST |
| M04 | Atropisomeric biaryl | 445.5 | 4.02 | 101.9 | 0.375 | 7 | 2/0 | 6 | MMFF94 | 71.3 | 11.5 | 1 | 2.93 | ATRO (4/4), ENTROPY |
| M05 | Perfluorinated cage | 395.2 | 1.37 | 66.4 | **0.846** | 3 | 0/7 | 5 | MMFF94 | **−0.15** | 7.9 | 0 | 2.78 | STRAIN, FLUORO, UST |
| M06 | Oxetane polyketide mimic | 460.6 | 2.61 | 110.3 | 0.522 | 5 | **7/0** | 7 | MMFF94 | 38.0 | 27.4 | 1 | 2.76 | STRAIN, ENTROPY |
| M07 | B–N dative macrocycle | 342.2 | 1.60 | 44.5 | 0.095 | 1 | 0/0 | 6 | **UFF** | 78.6 | **0.3** | 0 | 3.02 | ZWIT, UFF |
| M08 | Bicyclo-acrylamide warhead | 360.7 | 2.96 | 58.6 | 0.333 | 4 | 0/2 | 6 | MMFF94 | 8.4 | 6.8 | 0 | 2.89 | STRAIN, UST |
| M09 | Hetero-[5]-helicene | — | — | — | — | — | — | — | — | — | — | — | — | **unkekulizable SMILES** |
| M09R | Aza-[5]-helicene (repaired ref.) | 279.3 | 5.69 | 12.9 | 0.000 | 0 | 0/0 | 6 | MMFF94 | 97.2 | **0.0** | 0 | **4.13** | — |
| M10 | Tetra-ortho peptoid core | 422.6 | 5.51 | 40.6 | 0.481 | 4 | 0/0 | 6 | MMFF94 | **124.0** | 15.1 | 0 | 3.14 | ENTROPY |

† M03 is a dot-disconnected two-fragment assembly as written; graph and 3D ensemble computed on the 39-heavy-atom major fragment.

## Phase 2 — Heavy Dynamics Suite (T01–T05)

36-point relaxed torsional scans (MMFF94, 10° grid) + 200 ps OpenMM MD of the most flexible target (T02). All five targets are **PAINS-clean**.

| ID | Target | MW (Da) | cLogP | TPSA (Å²) | Fsp³ | RotB | SAScore | QED | ΔE‡ scan (kcal/mol) |
|---|---|---|---|---|---|---|---|---|---|
| T01 | CRBN molecular glue (imide mimic) | 470.9 | 2.21 | 90.0 | 0.30 | 3 | 3.10 | 0.69 | **7.8** (hinge, min@50°) |
| T02 | ADC Val-Cit-PAB linker | 604.7 | 1.68 | 218.5 | 0.50 | 14 | 3.59 | 0.13 | **21.4** (hinge, min@220°) |
| T03 | Bicyclic disulfide peptidomimetic | 465.6 | −1.81 | 145.5 | 0.64 | 0 | **6.57** | 0.21 | **72.8** (S–S strain probe) |
| T04 | Covalent inhibitor core (switch-II) | 515.0 | 4.99 | 66.3 | 0.18 | 5 | 2.94 | 0.30 | **10.5** (hinge, min@260°) |
| T05 | VHL-mimetic proline hybrid | 545.7 | 3.10 | 135.4 | 0.35 | 7 | 3.60 | 0.42 | **17.3** (hinge, min@110°) |

**MD highlights (T02, GBSA/OBC2, 300 K, 200 ps, 200 frames):** RMSD plateau **2.99 Å** (hinge breathing 2.1–4.8 Å) · Rg **7.60 ± 0.38 Å** · PE flat at −6280 kJ/mol · intramolecular H-bonds **0/frame** (audited: solvent-exposed by design) · ~970 steps/s on CPU.

### Phase 2 Figure Previews

**Fig. P2-1 — Relaxed torsional scans (5 panels + barrier bar chart)**

![Torsion scans](./figures_phase2/fig1_torsion_scans.png)

**Fig. P2-2 — MD trajectories (A: RMSD, B: Rg, C: PE/KE/TE)**

![MD trajectories](./figures_phase2/fig2_md_trajectories.png)

**Fig. P2-3 — Med-chem radar (SAScore / QED / Fsp³ / TPSA / MW)**

![MedChem radar](./figures_phase2/fig3_medchem_radar.png)

## Reproduce

```bash
pip install rdkit pandas numpy scikit-learn matplotlib seaborn torch torch_geometric

python molecule_benchmark.py --workers 4 --conformers 50   # Phase 1 -> bench_results/
python generate_assets.py                                  # Phase 1 figures -> ./figures/

# Phase 2 (MD parameterization uses a conda env with openff-toolkit; see DYNAMICS_REPORT_EN.md §6)
python run_heavy_dynamics_benchmark.py                     # full suite -> results_phase2/ + figures_phase2/
python run_heavy_dynamics_benchmark.py --fig_only          # regenerate figures from saved results
```

## Repository Structure

```
.
├── molecule_benchmark.py           # Phase 1: per-molecule pipeline (parse → scaffold → GNN tensors → 3D ensemble → analysis)
├── generate_assets.py              # Phase 1: publication figures from benchmark_results.json
├── run_heavy_dynamics_benchmark.py # Phase 2: descriptors → torsion scans → OpenMM MD → figures
├── export_openff_system.py         # Phase 2: OpenFF Sage → OpenMM system.xml (MMFF94 charge surgery)
├── setup_and_download.py           # one-shot environment/dataset bootstrap + self-check
├── BENCHMARK_REPORT_EN.md          # Phase 1 English technical report
├── BENCHMARK_REPORT_ZH.md          # 第一阶段中文技术报告
├── DYNAMICS_REPORT_EN.md           # Phase 2 English dynamics report
├── DYNAMICS_REPORT_ZH.md           # 第二阶段动力学报告
├── figures/                        # Phase 1 fig1–fig3 (300 DPI PNG)
├── figures_phase2/                 # Phase 2 fig1–fig3 (300 DPI PNG)
├── bench_results/                  # Phase 1 outputs
│   ├── benchmark_results.json      # machine-readable records
│   ├── benchmark_report.md         # auto-generated run report
│   ├── sdf/                        # *_ensemble.sdf, *_min.sdf
│   └── features/                   # *.npz (PyG-loadable)
└── results_phase2/                 # Phase 2 outputs (scans, system.xml, DCD, metrics)
```

> This repository also preserves the earlier sibling project [`aqueous-solubility-ml-benchmark`](https://github.com/songsiyi2006-chem/aqueous-solubility-ml-benchmark) (ESOL/AqSolDB solubility modeling) in its history — see the initial commit.

## License

MIT — see [`LICENSE`](./LICENSE).
