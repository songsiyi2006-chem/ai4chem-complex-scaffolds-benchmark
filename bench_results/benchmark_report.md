# Hardcore Benchmark: 10 Novel Complex Molecules

*Generated 2026-09-02T13:43:26 | RDKit 2026.03.5 | 50 conformers/molecule, ETKDGv3 + MMFF94/UFF | wall time 27.3 s*

## Benchmark Table

| # | Molecule | MW (Da) | cLogP | TPSA (Å²) | Fsp³ | RotB | Stereo a/u | MaxRing | Conf | FF | E_min | ΔE_ens | IMHB | d_min (Å) | Risk |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| M01 | Macrocyclic chameleon peptidomimetic | 529.6 | -0.72 | 165.2 | 0.538 | 5 | 4/0 | 13 | 50 | MMFF94 | 67.4 | 43.1 | 2 | 2.73 | MAC, ENTROPY |
| M02 | Strained azaspiro-cubane bioisostere | 366.3 | 1.78 | 32.8 | 0.611 | 2 | 0/5 | 6 | 50 | MMFF94 | 106.3 | 8.6 | 0 | 2.78 | STRAIN, UST |
| M03 | Bivalent heterobifunctional degron (PROTAC prototype) | 815.9 | 2.51 | 241.4 | 0.282 | 16 | 0/1 | 6 | 50 | MMFF94 | 29.2 | 12.2 | 1 | 2.75 | BIG, FLEX, FRAG, ENTROPY, UST |
| M04 | Axially + centrally chiral atropisomeric biaryl | 445.5 | 4.02 | 101.9 | 0.375 | 7 | 2/0 | 6 | 50 | MMFF94 | 71.3 | 11.5 | 1 | 2.93 | ATRO, ENTROPY |
| M05 | Perfluorinated polycyclic cage vector | 395.2 | 1.37 | 66.4 | 0.846 | 3 | 0/7 | 5 | 50 | MMFF94 | -0.1 | 7.9 | 0 | 2.78 | STRAIN, FLUORO, UST |
| M06 | Oxetane-fused polyketide NP mimic | 460.6 | 2.61 | 110.3 | 0.522 | 5 | 7/0 | 7 | 50 | MMFF94 | 38.0 | 27.4 | 1 | 2.76 | STRAIN, ENTROPY |
| M07 | Zwitterionic B-N dative macrocycle | 342.2 | 1.60 | 44.5 | 0.095 | 1 | 0/0 | 6 | 50 | UFF | 78.6 | 0.3 | 0 | 3.02 | ZWIT, UFF |
| M08 | Strained bicyclo-acrylamide covalent warhead | 360.7 | 2.96 | 58.6 | 0.333 | 4 | 0/2 | 6 | 50 | MMFF94 | 8.4 | 6.8 | 0 | 2.89 | STRAIN, UST |
| M09 | Heteroatom-doped [5]-carbohelicene | parse FAILED: see analysis | - | - | - | - | - | - | - | - | - | - | - | - | - | |
| M09R | Aza-[5]-carbohelicene (repaired reference for M09) | 279.3 | 5.69 | 12.9 | 0.000 | 0 | 0/0 | 6 | 50 | MMFF94 | 97.2 | 0.0 | 0 | 4.13 | - |
| M10 | Tetra-ortho hindered peptoid core | 422.6 | 5.51 | 40.6 | 0.481 | 4 | 0/0 | 6 | 50 | MMFF94 | 124.0 | 15.1 | 0 | 3.14 | ENTROPY |

*E_min and ΔE_ens: total MMFF94/UFF steric energy in kcal/mol of the lowest-energy conformer and the ensemble spread (E_max − E_min). Energies are relative force-field quantities, not formation enthalpies. Stereo a/u = assigned/unassigned chiral centers; d_min = shortest non-bonded heavy-atom distance (≥4 bonds apart) in the E_min conformer.*

**Risk flag legend:** `MAC` macrocycle (>=12-ring): transannular effects invisible to 2D GNNs; `STRAIN` contains 3/4-membered rings: high angle strain, distance-geometry edge case; `ATRO` ortho-blocked biaryl axis (>=3/4 positions): atropisomerism, axial chirality lost in 2D; `FLUORO` >=6 fluorines: dense C-F electrostatics, LogP/vdW model extrapolation; `ZWIT` formal charges / zwitterion / dative bonding: OOD for typical pretrained 2D GNNs; `BIG` MW > 800 Da: outside Ro5 / typical GNN pretraining domain; `FLEX` >=15 rotatable bonds: conformational entropy, single-graph under-sampling; `FRAG` disconnected multi-fragment input: 3D run on largest fragment; `ENTROPY` ensemble dE >= 10 kcal/mol: extreme conformational polymorphism; `CLASH` non-bonded heavy-atom pair < 2.5 A in Emin conformer: steric clash; `UFF` MMFF94 params unavailable -> UFF fallback (lower accuracy); `UST` unassigned stereocenters detected

## GNN Feature Tensor Verification (PyG-ready)

| # | Nodes | Directed edges | Atom feat dim | Feature file | GCN smoke |
|---|---|---|---|---|---|
| M01 | 38 | 80 | 8 | `bench_results\features\M01.npz` | OK (torch 2.13.0+cpu, pyg 2.8.0.post1, out [38, 64]) |
| M02 | 26 | 64 | 8 | `bench_results\features\M02.npz` | OK (torch 2.13.0+cpu, pyg 2.8.0.post1, out [26, 64]) |
| M03 | 39 | 82 | 8 | `bench_results\features\M03.npz` | OK (torch 2.13.0+cpu, pyg 2.8.0.post1, out [39, 64]) |
| M04 | 32 | 66 | 8 | `bench_results\features\M04.npz` | OK (torch 2.13.0+cpu, pyg 2.8.0.post1, out [32, 64]) |
| M05 | 26 | 60 | 8 | `bench_results\features\M05.npz` | OK (torch 2.13.0+cpu, pyg 2.8.0.post1, out [26, 64]) |
| M06 | 32 | 70 | 8 | `bench_results\features\M06.npz` | OK (torch 2.13.0+cpu, pyg 2.8.0.post1, out [32, 64]) |
| M07 | 26 | 60 | 8 | `bench_results\features\M07.npz` | OK (torch 2.13.0+cpu, pyg 2.8.0.post1, out [26, 64]) |
| M08 | 24 | 52 | 8 | `bench_results\features\M08.npz` | OK (torch 2.13.0+cpu, pyg 2.8.0.post1, out [24, 64]) |
| M09 | - | - | - | - | n/a (parse failed) |
| M09R | 22 | 52 | 8 | `bench_results\features\M09R.npz` | OK (torch 2.13.0+cpu, pyg 2.8.0.post1, out [22, 64]) |
| M10 | 31 | 64 | 8 | `bench_results\features\M10.npz` | OK (torch 2.13.0+cpu, pyg 2.8.0.post1, out [31, 64]) |

## Per-Molecule Failure Analysis

### M01 — Macrocyclic chameleon peptidomimetic（大环变色龙拟肽）
- **SMILES:** `O=C1N[C@H](C(C)C)C(=O)N[C@@H](Cc2ccccc2)C(=O)N(C)[C@H](CC(=O)NC1)C(=O)N2CCC[C@H]2C(=O)O`
- **Bemis-Murcko scaffold (generic):** `CC1CCC(C)CC(CC2CCCCC2)C(C)CC(C(C)C2CCCC2)CC(C)CC1` (scaffold atom fraction 0.816)
- **Challenge:** Transannular H-bond networks, solvent-dependent conformational switching, extreme conformational entropy.
- **Parse status:** sanitized
- **IMHB in E_min conformer:** 2 — D9-H47···A0 (2.02 Å, 147.8°); D2-H38···A25 (2.33 Å, 136.0°)
- **Ensemble:** 50 confs optimized with MMFF94 (0 not fully converged); E ∈ [67.37, 110.42] kcal/mol
- **2D-GNN failure modes:** `MAC` macrocycle (>=12-ring): transannular effects invisible to 2D GNNs; `ENTROPY` ensemble dE >= 10 kcal/mol: extreme conformational polymorphism
- **Runtime:** 18.3 s

### M02 — Strained azaspiro-cubane bioisostere（稠合氮杂螺环-立方烷生物电子等排体）
- **SMILES:** `O=C(N1CC23C4C1C2C34)c1c(F)c(F)c(N5CC6(COC6)C5)c(F)c1F`
- **Bemis-Murcko scaffold (generic):** `CC(C1CCC(C2CC3(CCC3)C2)CC1)C1CC23C4C1C2C43` (scaffold atom fraction 0.846)
- **Challenge:** Extreme C-C-C angle distortion, bridgehead nitrogen strain, hard distance geometry.
- **Parse status:** sanitized
- **IMHB in E_min conformer:** 0 — none detected
- **Ensemble:** 50 confs optimized with MMFF94 (0 not fully converged); E ∈ [106.33, 114.92] kcal/mol
- **2D-GNN failure modes:** `STRAIN` contains 3/4-membered rings: high angle strain, distance-geometry edge case; `UST` unassigned stereocenters detected
- **Runtime:** 6.5 s

### M03 — Bivalent heterobifunctional degron (PROTAC prototype)（新型刚柔偶联PROTAC分子）
- **SMILES:** `O=C1c2ccccc2C(=O)N1C3CCC(=O)NC3=O.O=C(NCCOCCOCCOc4ccc(NC(=O)c5cnc(Nc6ccc(S(=O)(=O)C)cc6)nc5)cc4)C`
- **Bemis-Murcko scaffold (generic):** `CC(CC1CCCCC1)C1CCC(CC2CCCCC2)CC1` (scaffold atom fraction 0.564)
- **Challenge:** MW > 850 Da, 16+ rotatable linker bonds, inter-domain steric clashes on conformer collapse; input is a 2-fragment disconnected assembly.
- **Parse status:** sanitized
- **Note:** input is a 2-fragment disconnected assembly; bulk properties reported for the full assembly, graph/3D on the largest fragment (39 heavy atoms)
- **IMHB in E_min conformer:** 1 — D16-H54···A0 (1.82 Å, 149.7°)
- **Ensemble:** 50 confs optimized with MMFF94 (0 not fully converged); E ∈ [29.16, 41.33] kcal/mol
- **2D-GNN failure modes:** `BIG` MW > 800 Da: outside Ro5 / typical GNN pretraining domain; `FLEX` >=15 rotatable bonds: conformational entropy, single-graph under-sampling; `FRAG` disconnected multi-fragment input: 3D run on largest fragment; `ENTROPY` ensemble dE >= 10 kcal/mol: extreme conformational polymorphism; `UST` unassigned stereocenters detected
- **Runtime:** 16.3 s

### M04 — Axially + centrally chiral atropisomeric biaryl（轴手性与点手性结合的双芳基抑制剂）
- **SMILES:** `CC(=O)Oc1c(C)c(c2c(OC(=O)C)c(C)ccc2[C@@H](C)NC(=O)CF)ccc1[C@H](C)O`
- **Bemis-Murcko scaffold (generic):** `C1CCC(C2CCCCC2)CC1` (scaffold atom fraction 0.375)
- **Challenge:** High biaryl rotation barrier (>30 kcal/mol), perpendicular aromatic planes, non-planar conjugation.
- **Parse status:** sanitized
- **Atropisomer scan:** 1 rotatable biaryl axes, max ortho-blocking 4/4
- **IMHB in E_min conformer:** 1 — D31-H59···A2 (2.34 Å, 124.1°)
- **Ensemble:** 50 confs optimized with MMFF94 (0 not fully converged); E ∈ [71.34, 82.83] kcal/mol
- **2D-GNN failure modes:** `ATRO` ortho-blocked biaryl axis (>=3/4 positions): atropisomerism, axial chirality lost in 2D; `ENTROPY` ensemble dE >= 10 kcal/mol: extreme conformational polymorphism
- **Runtime:** 11.9 s

### M05 — Perfluorinated polycyclic cage vector（全氟代多环笼状分子）
- **SMILES:** `FC1(F)C2(F)C3(F)C1(F)C4(F)C2(F)C3(F)C4(C(=O)NC5(CC5)C(=O)O)(F)`
- **Bemis-Murcko scaffold (generic):** `CC(CC1CC1)C1C2C3CC4C3C1C42` (scaffold atom fraction 0.538)
- **Challenge:** Dense fluorine electrostatics, extreme hydrophobicity, abnormal surface/volume ratio, electron-deficient core.
- **Parse status:** sanitized
- **IMHB in E_min conformer:** 0 — none detected
- **Ensemble:** 50 confs optimized with MMFF94 (0 not fully converged); E ∈ [-0.15, 7.72] kcal/mol
- **2D-GNN failure modes:** `STRAIN` contains 3/4-membered rings: high angle strain, distance-geometry edge case; `FLUORO` >=6 fluorines: dense C-F electrostatics, LogP/vdW model extrapolation; `UST` unassigned stereocenters detected
- **Runtime:** 11.6 s

### M06 — Oxetane-fused polyketide NP mimic（含氧杂四元环的多手性中心天然产物类似物）
- **SMILES:** `C[C@H]1O[C@@]2(CO2)[C@@H](O)[C@H](C)[C@@H](OC(=O)c3ccccc3)[C@H]1C(=O)N[C@H](C)c4nc(C)cs4`
- **Bemis-Murcko scaffold (generic):** `CC(CC1CCC2(CCC1C(C)CCC1CCCC1)CC2)C1CCCCC1` (scaffold atom fraction 0.844)
- **Challenge:** 7 stereocenters, constrained oxygenated small ring, dense intramolecular steric repulsion.
- **Parse status:** sanitized
- **IMHB in E_min conformer:** 1 — D7-H39···A22 (1.87 Å, 149.4°)
- **Ensemble:** 50 confs optimized with MMFF94 (0 not fully converged); E ∈ [38.04, 65.47] kcal/mol
- **2D-GNN failure modes:** `STRAIN` contains 3/4-membered rings: high angle strain, distance-geometry edge case; `ENTROPY` ensemble dE >= 10 kcal/mol: extreme conformational polymorphism
- **Runtime:** 11.7 s

### M07 — Zwitterionic B-N dative macrocycle（含B-N配位内盐大环化合物）
- **SMILES:** `c1ccc2c(c1)[B-]3(c4ccccc42)OCCN[N+]3=Cc5cccc(O)c5`
- **Bemis-Murcko scaffold (generic):** `C1CCC(CC2CCCCC23C2CCCCC2C2CCCCC23)CC1` (scaffold atom fraction 0.962)
- **Challenge:** Dative B-N bond, formal B(-)/N(+) charges, tetrahedral borate coordination; MMFF lacks B params (UFF fallback).
- **Parse status:** sanitized
- **IMHB in E_min conformer:** 0 — none detected
- **Ensemble:** 50 confs optimized with UFF (0 not fully converged); E ∈ [78.56, 78.82] kcal/mol
- **2D-GNN failure modes:** `ZWIT` formal charges / zwitterion / dative bonding: OOD for typical pretrained 2D GNNs; `UFF` MMFF94 params unavailable -> UFF fallback (lower accuracy)
- **Runtime:** 2.9 s

### M08 — Strained bicyclo-acrylamide covalent warhead（张力双环丙烷共价弹头分子）
- **SMILES:** `C=CC(=O)N1CC2(CC12)C(=O)Nc3ccc(OC(F)(F)F)c(Cl)c3`
- **Bemis-Murcko scaffold (generic):** `CC(CC1CCCCC1)C12CCC1C2` (scaffold atom fraction 0.583)
- **Challenge:** Electrophilic warhead alignment, localized ring strain, reactive dihedral vectoring for cysteine trapping.
- **Parse status:** sanitized
- **IMHB in E_min conformer:** 0 — none detected
- **Ensemble:** 50 confs optimized with MMFF94 (0 not fully converged); E ∈ [8.43, 15.27] kcal/mol
- **2D-GNN failure modes:** `STRAIN` contains 3/4-membered rings: high angle strain, distance-geometry edge case; `UST` unassigned stereocenters detected
- **Runtime:** 1.8 s

### M09 — Heteroatom-doped [5]-carbohelicene（杂原子掺杂五螺烯手性发光分子）
- **SMILES:** `c1cc2c(s1)c3ccc4c(c3c2)c5ccc6ncccc6c5c7cccnc47`
- **Challenge:** Helical non-planar aromatic distortion, overlapping terminal rings, optical asymmetry.
- **Parse status:** sanitization failed at Kekulize: Can't kekulize mol.  Unkekulized atoms: 0 1 2 3 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
- **Note:** PARSE FAILURE: sanitization failed at Kekulize: Can't kekulize mol.  Unkekulized atoms: 0 1 2 3 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
- **Runtime:** - s

### M09R — Aza-[5]-carbohelicene (repaired reference for M09)（氮杂[5]螺烯（M09修复参照））
- **SMILES:** `c1ccc2c(c1)ccc1c2ccc2c3ccncc3ccc21`
- **Bemis-Murcko scaffold (generic):** `C1CCC2C(C1)CCC1C2CCC2C3CCCCC3CCC21` (scaffold atom fraction 1.0)
- **Challenge:** Reference structure: the M09 SMILES as supplied is unkekulizable (invalid ring-fusion aromaticity); this program-generated aza-doped [5]-helicene (C21H13N) demonstrates the intended helical topology and failure modes.
- **Parse status:** sanitized
- **IMHB in E_min conformer:** 0 — none detected
- **Ensemble:** 50 confs optimized with MMFF94 (0 not fully converged); E ∈ [97.19, 97.19] kcal/mol
- **Runtime:** 1.4 s

### M10 — Tetra-ortho hindered peptoid core（四邻位超拥挤拟肽骨架）
- **SMILES:** `CC1=C(C)C(=C(C)C(=C1C)N(C)C(=O)CN(C)C(=O)c2c(C)c(C)c(C)c(C)c2C)C`
- **Bemis-Murcko scaffold (generic):** `CC(CCC(C)C1CCCCC1)CC1CCCCC1` (scaffold atom fraction 0.613)
- **Challenge:** Severe steric clash blocking amide cis/trans planarization; hindered-rotation energy barriers.
- **Parse status:** sanitized
- **IMHB in E_min conformer:** 0 — none detected
- **Ensemble:** 50 confs optimized with MMFF94 (0 not fully converged); E ∈ [123.99, 139.09] kcal/mol
- **2D-GNN failure modes:** `ENTROPY` ensemble dE >= 10 kcal/mol: extreme conformational polymorphism
- **Runtime:** 6.7 s

## Artifacts
- `sdf/<ID>_ensemble.sdf` — full optimized ensemble (E_kcal_per_mol per conf)
- `sdf/<ID>_min.sdf` — lowest-energy conformer (ΔE = 0 reference)
- `features/<ID>.npz` — atom/bond feature tensors, PyG-loadable
- `benchmark_results.json` — complete machine-readable records
