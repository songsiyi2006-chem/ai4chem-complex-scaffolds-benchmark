# Stress-Testing Frontier 3D Conformer Pipelines on Unindexed, Structurally Complex Molecular Entities

**A 10-Molecule Hard-Core Benchmark of RDKit ETKDGv3 + MMFF94/UFF Conformational Search, Physicochemical Profiling, and GNN Featurization Readiness**

*Technical Report — generated 2026-09-02 · RDKit 2026.03.5 · PyTorch 2.13.0+cpu · PyTorch Geometric 2.8.0.post1 · all numbers are computed values serialized in `bench_results/benchmark_results.json`*

---

## 1. Executive Summary & Objective

This benchmark stress-tests a fully open-source, physics-based molecular modeling stack — **ETKDGv3 experimental-torsion distance geometry followed by MMFF94 (UFF fallback) relaxation** — against ten structurally exotic, synthetically plausible, deliberately *unindexed* scaffold classes: a conformationally labile macrocyclic peptidomimetic, a fused azaspiro/polycyclic cubane bioisostere, a beyond-Ro5 heterobifunctional degron (PROTAC prototype), a fully ortho-blocked atropisomeric biaryl, a perfluorinated polycyclic cage, an oxetane/epoxide-bearing polyketide mimic with seven stereocenters, a zwitterionic B–N dative macrocycle, a strained bicyclo-acrylamide covalent warhead, a heteroatom-doped [5]-carbohelicene, and a tetra-ortho-sterically-hindered peptoid core.

**Headline outcomes:**

- **9 of 10 targets completed the full pipeline** (parse → sanitize → stereo perception → Bemis–Murcko scaffolding → GNN tensor extraction → 50-conformer ETKDGv3 ensemble → MMFF94/UFF optimization → lowest-energy selection → IMHB/steric geometry analysis); 500 conformer optimizations converged with **zero non-converged structures**, in **27.3 s** wall time on 4 worker processes.
- **1 of 10 (M09, the helicene) failed at the SMILES level**: the supplied aromatic ring-closure pattern is *unkekulizable* — a genuine input-integrity defect caught before any compute was wasted. A repaired aza-[5]-helicene reference (**M09R**, C₂₁H₁₃N, programmatically generated and validated) completes the intended helical case study (ΔE_ensemble = 0.0 kcal/mol; bay-region C···C contact 4.13 Å confirming the helical pitch).
- The set deliberately spans **unconventional drug-like space**: MW 279–816 Da, cLogP −0.72 to 5.69, TPSA 12.9–241.4 Å², Fsp³ 0.000–0.846; ensemble energy spreads (ΔE_ens, a conformational-entropy proxy) range from **0.0 to 43.1 kcal/mol** — the latter quantifying the "chameleon" polymorphism of the M01 macrocycle.

![Figure 1 — 2D structural grid of the benchmark set](./figures/fig1_molecular_grid.png)

*Figure 1. The benchmark set. Ten target molecules plus the M09R repaired reference. M09 is shown as a parse-failure placeholder: its SMILES cannot be kekulized (Section 4).*

---

## 2. Chemical Space & Structural Diversity

![Figure 2 — chemical space of the benchmark set](./figures/fig2_chemical_space.png)

*Figure 2. MW vs. cLogP bubble map (bubble area ∝ TPSA, color = Fsp³). The dashed rectangle marks the Lipinski Ro5 reference zone. Four of ten entries violate at least one Ro5 bound; the set deliberately occupies beyond-Ro5 and fragment-like corners simultaneously.*

**Topology inventory.** The set contains one 13-membered macrocycle (M01), one fused azaspiro bicyclic cage with 4/5-membered rings (M02), one disconnected two-fragment assembly as written (M03), a 4/4-ortho-blocked biaryl axis (M04), an eight-fold fluorinated cage (M05), an epoxide-fused 7-membered oxacycle with seven fully specified stereocenters (M06), a B⁻/N⁺ zwitterionic dative macrocycle (M07), a cyclopropane-fused pyrrolidine warhead (M08), a rigid aza-helicene (M09R), and a pentamethyl-aryl hindered bis-amide (M10). Bemis–Murcko scaffold atom fractions span **0.375 (M04) to 1.000 (M09R)** — i.e., from decorated small-scaffold molecules to whole-molecule scaffolds — which matters directly for scaffold-split generalization estimates.

**Ring strain (M02, M05, M08).** M02's fused bicyclo[2.1.0]-pentane-type core forces C–C–C angles toward 60–90°; its E_min = 106.3 kcal/mol is among the highest in the set, and its five bridgehead/spiro stereocenters are *unassigned in the input SMILES* (UST flag) — the conformational ensemble therefore averages over unspecified diastereomers. M05's perfluorinated cage is the most sp³-dense entry (Fsp³ = 0.846) yet is geometrically well-behaved (E_min ≈ −0.15 kcal/mol, the set minimum). M08's cyclopropane-fused warhead contributes ΔE_ens = 6.8 kcal/mol of acrylamide-orientation polymorphism.

**Flexibility (M03).** At 815.9 Da (as written: a dot-disconnected two-fragment assembly, bulk properties computed on the full assembly; graph and 3D ensemble on the 39-heavy-atom major fragment) with **16 rotatable bonds**, M03 exemplifies beyond-Ro5 degrader space; its lowest-energy conformer folds into an intramolecular H-bond (1.82 Å, 149.7°), a collapse mode that any single-conformer representation will miss.

**Transannular effects (M01).** The 13-membered cyclic-tetrapeptide macrocycle exhibits **2 intramolecular H-bonds in its E_min conformer** (donor–H···acceptor 2.02 Å at 147.8° — a strong, nearly linear contact — plus 2.33 Å at 136.0°) and the **largest ensemble spread of the set (ΔE_ens = 43.1 kcal/mol)**, quantifying its "conformational chameleon" character: multiple disjoint H-bond networks compete across the ensemble.

---

## 3. Physics-Based Conformation & Energy Landscape

**Method.** Each molecule (largest fragment if disconnected) was protonated (`AddHs`), embedded with **ETKDGv3** (experimental torsion preferences + basic knowledge; `useSmallRingTorsions=True`, `useMacrocycleTorsions=True`, seed-controlled, no RMS pruning) to exactly **50 conformers**, then batch-optimized with **MMFF94** (`MMFFOptimizeMoleculeConfs`, maxIters = 2000; all converged, flag = 0). Where MMFF94 parameterization is impossible, the pipeline falls back to **UFF** — which occurred exactly once, for M07.

**MMFF94 vs. UFF parameterization limits (M07).** The B–N dative macrocycle carries formal B⁻/N⁺ charges and a four-coordinate borate; MMFF94 has no boron atom types, so `MMFFHasAllMoleculeParams` correctly returned False and UFF (which parameterizes boron with a universal geometric rule set) took over. Two consequences deserve emphasis:

1. **Absolute UFF energies are not comparable** to MMFF94 totals across molecules (UFF's diagonal-harmonic approximation lacks MMFF's charged dstn/charge–charge electrostatics), so M07's E_min = 78.6 kcal/mol must not be ranked against MMFF entries. Within-molecule *relative* comparisons remain meaningful.
2. The dative B–N bond **rigidifies the macrocycle dramatically**: M07's ensemble collapses to ΔE_ens = 0.3 kcal/mol — two orders of magnitude tighter than the non-dative M01 macrocycle (43.1 kcal/mol). Coordination chemistry here acts as a conformational lock, an effect a 2D topology model cannot express.

**Energy-landscape reading.** E_min spans −0.15 (M05) to 124.0 kcal/mol (M10) — the latter being the direct steric-crowding cost of four ortho methyls flanking two tertiary amides. ΔE_ens stratifies the set into rigid (M07: 0.3; M09R: 0.0 — a 50-conformer ensemble degenerate to a single helical basin), moderately flexible (M02: 8.6; M05: 7.9; M08: 6.8), and conformationally polymorphic (M10: 15.1; M04: 11.5; M03: 12.2; M06: 27.4; M01: 43.1). Minimum non-bonded heavy-atom distances (d_min ≥ 2.73 Å for all entries; M01's 2.73 Å coincides with its strong transannular H-bond donor–acceptor pair) show **no catastrophic steric clashes** survive in the relaxed minima.

**Intramolecular H-bond inventory (E_min conformers; criteria d(H···A) ≤ 2.5 Å, ∠(D–H···A) ≥ 120°, ≥ 3 bonds apart):** M01: 2 · M03: 1 (1.82 Å) · M04: 1 (2.34 Å) · M06: 1 (1.87 Å) · others: 0.

---

## 4. Edge-Case Investigation: The M09 / M09R Case Study

**Symptom.** `Chem.MolFromSmiles("c1cc2c(s1)c3ccc4c(c3c2)c5ccc6ncccc6c5c7cccnc47")` terminates with `Can't kekulize mol. Unkekulized atoms: 0 1 2 3 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27` — 27 of 28 aromatic atoms resist assignment of an alternating single/double-bond (Kekulé) pattern.

**Forensic analysis.** Tracing the ring-closure digits reveals the defect: the closure tagged `c2` (opened at atom 3, closed inside the branch of atom 9) does not complete an ortho-fused six-ring as a helicene requires; combined with the `c3` closure it instead carves out a **five-membered all-carbon ring fused between the thiophene and a benzene unit**, while the terminal `c47` double-closure spans an unintended long cycle. The resulting fused graph is **non-alternant with no perfect matching** in its π-bond graph — its Kekulé number is zero — so no valid electron-pair assignment exists, regardless of toolkit. This is a *structure-authoring* defect (ring-digit bookkeeping), not a toolkit limitation; the supplied molecule, as written, does not exist as a closed-shell aromatic species.

**Repair protocol (M09R).** Rather than hand-editing closure digits (error-prone, as demonstrated), the repair was programmatic:

1. Start from a *validated* parent [5]-carbohelicene, `c1ccc2c(c1)ccc3c2ccc4c3ccc5c4cccc5` (C₂₂H₁₄, five ortho-fused six-rings — confirmed by RDKit ring perception);
2. Perform a single programmatic aromatic CH→N substitution at a non-fusion, degree-2 carbon (preserving pyridine-type two-coordinate N), re-sanitize, and round-trip the canonical SMILES;
3. Retain the first chemically valid isomer: **M09R = `c1ccc2c(c1)ccc1c2ccc2c3ccncc3ccc21` (C₂₁H₁₃N)**.

**M09R results validate the intended failure mode.** The repaired helicene runs the full pipeline cleanly: ΔE_ens = 0.0 kcal/mol (an entirely rigid helix), Fsp³ = 0.000, cLogP = 5.69, TPSA = 12.9 Å², and a **bay-region terminal-ring C···C distance of 4.13 Å** — the geometric signature of overlapping helical termini. Two further lessons generalize: (i) aromatic SMILES from generative sources should always pass a *kekulization round-trip* gate before entering any dataset; (ii) helicity is expressed only through coordinates — no 2D graph invariant in the pipeline distinguishes M09R's P/M helices from its mirror image.

---

## 5. Summary Benchmark Table

All values computed (see `bench_results/benchmark_results.json`); E in kcal/mol; ΔE_ens = E_max − E_min over 50 optimized conformers. Stereo a/u = assigned/unassigned chiral centers. d_min = shortest non-bonded heavy-atom distance (≥ 4 bonds apart) in the E_min conformer.

| # | Molecule | MW (Da) | cLogP | TPSA (Å²) | Fsp³ | RotB | Stereo a/u | MaxRing | FF | E_min | ΔE_ens | IMHB | d_min (Å) | Risk flags |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| M01 | Macrocyclic chameleon peptide | 529.6 | −0.72 | 165.2 | 0.538 | 5 | 4/0 | 13 | MMFF94 | 67.4 | **43.1** | 2 | 2.73 | MAC, ENTROPY |
| M02 | Azaspiro-cubane bioisostere | 366.3 | 1.78 | 32.8 | 0.611 | 2 | 0/**5** | 6 | MMFF94 | 106.3 | 8.6 | 0 | 2.78 | STRAIN, UST |
| M03 | PROTAC prototype (2 fragments) | 815.9 | 2.51 | 241.4 | 0.282 | 16 | 0/1 | 6 | MMFF94 | 29.2† | 12.2† | 1 | 2.75† | BIG, FLEX, FRAG, ENTROPY, UST |
| M04 | Atropisomeric biaryl | 445.5 | 4.02 | 101.9 | 0.375 | 7 | 2/0 | 6 | MMFF94 | 71.3 | 11.5 | 1 | 2.93 | **ATRO (4/4)**, ENTROPY |
| M05 | Perfluorinated cage | 395.2 | 1.37 | 66.4 | **0.846** | 3 | 0/7 | 5 | MMFF94 | **−0.15** | 7.9 | 0 | 2.78 | STRAIN, FLUORO, UST |
| M06 | Oxetane polyketide mimic | 460.6 | 2.61 | 110.3 | 0.522 | 5 | **7/0** | 7 | MMFF94 | 38.0 | 27.4 | 1 | 2.76 | STRAIN, ENTROPY |
| M07 | B–N dative macrocycle | 342.2 | 1.60 | 44.5 | 0.095 | 1 | 0/0 | 6 | **UFF** | 78.6 | **0.3** | 0 | 3.02 | ZWIT, UFF |
| M08 | Bicyclo-acrylamide warhead | 360.7 | 2.96 | 58.6 | 0.333 | 4 | 0/2 | 6 | MMFF94 | 8.4 | 6.8 | 0 | 2.89 | STRAIN, UST |
| M09 | Hetero-[5]-helicene | — | — | — | — | — | — | — | — | — | — | — | — | **unkekulizable SMILES** |
| M09R | Aza-[5]-helicene (repaired ref.) | 279.3 | 5.69 | 12.9 | 0.000 | 0 | 0/0 | 6 | MMFF94 | 97.2 | **0.0** | 0 | **4.13** | — |
| M10 | Tetra-ortho peptoid core | 422.6 | 5.51 | 40.6 | 0.481 | 4 | 0/0 | 6 | MMFF94 | **124.0** | 15.1 | 0 | 3.14 | ENTROPY |

† M03 is a dot-disconnected two-fragment assembly as written; graph and 3D ensemble computed on the 39-heavy-atom major fragment, bulk properties on the full assembly.

**GNN featurization readiness** (integer-encoded 8-dim atom features / 4-dim bond features, bidirectional edge index; every graph passed a `GCNConv` forward smoke test under PyG 2.8.0.post1):

| ID | Nodes (heavy atoms) | Directed edges | ID | Nodes | Directed edges |
|---|---|---|---|---|---|
| M01 | 38 | 80 | M06 | 32 | 70 |
| M02 | 26 | 64 | M07 | 26 | 60 |
| M03 | 39 | 82 | M08 | 24 | 52 |
| M04 | 32 | 66 | M09R | 22 | 52 |
| M05 | 26 | 60 | M10 | 31 | 64 |

**Interpretation highlights.** (i) ΔE_ens is the single most informative stability axis in this set: it cleanly separates rigid coordination-locked (M07) and helically over-determined (M09R) scaffolds from entropic macrocyclic (M01) and hindered-rotamer (M10) systems. (ii) UST flags on M02/M05/M08 mean their *inputs* under-specify stereochemistry (5, 7, and 2 unassigned centers respectively) — any ML training on such inputs mixes diastereomer classes unless stereoisomers are enumerated or the unspecified atoms are made explicit. (iii) M10's E_min of 124.0 kcal/mol is the set's steric-crowding ceiling; combined with ΔE_ens = 15.1 it quantifies the amide-rotation frustration that defines the "tetra-ortho" peptoid design.

---

## 6. Implications for AI4Science & Geometric GNNs

![Figure 3 — complexity radar of representative scaffolds](./figures/fig3_radar_complexity.png)

*Figure 3. Normalized complexity fingerprints (each axis scaled to the benchmark-set maximum) for the four most structurally divergent scaffolds: macrocycle (M01), cubane-type cage (M02), PROTAC (M03), perfluorinated cage (M05).*

**Why 2D message-passing GNNs fail on this set.** A 2D GNN consumes only the molecular graph; by construction it is invariant to everything this benchmark shows to be decisive:

- **Conformational polymorphism.** M01 spans 43.1 kcal/mol of interconverting conformers with competing transannular H-bond networks; a single-graph 2D embedding necessarily averages ("smears") chameleon behavior into one point estimate with no uncertainty signal.
- **Axial chirality.** M04's biaryl axis carries 4/4 ortho blocking — an atropisomeric stereogenic element that is *not* an atom-level chiral tag. Standard 2D featurizations (including the chiral-tag channel used here) cannot represent it; the two atropisomers are one graph.
- **Strain energetics.** Cubane-like M02 and the M05 cage live in a narrow valid-geometry manifold (60–90° C–C–C angles); 2D models trained on relaxed drug-like space extrapolate blindly here, and even 3D pipelines need `useSmallRingTorsions`-style priors to embed them at all.
- **Charge/dative bonding.** M07's B⁻/N⁺ four-coordinate borate is out-of-distribution for atom-type vocabularies of most pretrained 2D models, and its chemistry (UFF-not-MMFF) exposes label-provenance issues: any force-field-derived property in a training set must carry its FF tag.
- **Helicity.** M09R's P/M enantiomers are graph-identical; only coordinates discriminate them — and only an equivariant model can keep that information through the layers.

**What to use instead.** These are precisely the failure modes that SE(3)/E(3)-equivariant architectures were designed to absorb: **MACE-style higher-order irreducible-representation message passing** (and e3nn-built equivalents) acting on conformer ensembles preserves directional/rotational information; PaiNN/TorchMD-Net vector-atom designs offer a lighter-weight route. Concretely, this benchmark motivates four design rules: (1) feed *ensembles* (or their amortized distribution), not single conformers, and treat the ensemble spread (our ΔE_ens) as an explicit aleatoric-uncertainty feature — M01 vs M07 should not look equally "certain"; (2) enforce stereochemical completeness at ingestion (the UST flags here are a data-quality gate, not a nuisance); (3) gate generative/aromatic inputs through kekulization round-trips (the M09 lesson) before they reach any loader; and (4) respect evaluation discipline: with Bemis–Murcko scaffold fractions from 0.375 to 1.000, random splits would overestimate generalization on this set — scaffold-based (or scaffold-generic) splits are mandatory for any model claiming performance on exotic topologies.

---

## Reproducibility

```bash
python molecule_benchmark.py --workers 4 --conformers 50   # full pipeline → bench_results/
python generate_assets.py                                  # figures → ./figures/ (300 DPI)
```

Artifacts: `bench_results/benchmark_results.json` (machine-readable records) · `bench_results/benchmark_report.md` (auto-generated run report) · `bench_results/sdf/{*_ensemble.sdf, *_min.sdf}` · `bench_results/features/*.npz` (PyG-loadable tensors) · `figures/fig1–fig3` (this report).

*Energies are total MMFF94/UFF steric energies — relative, force-field-internal quantities, not formation enthalpies; cross-force-field (M07 vs. rest) absolute comparisons are invalid.*
