# DYNAMICS REPORT — Phase 2 Heavy Compute Suite

**3D biophysical simulation, relaxed torsional barrier scanning & OpenMM molecular dynamics**

*Generated from `results_phase2/phase2_results.json` | RDKit 2026.03.5 (Python 3.14.6) · OpenMM 8.6 (CPU) · OpenFF Sage 2.1.0 + MMFF94 charges · conda env `phase2ff` for parameterization*

---

## 1. Executive Summary

Phase 2 extends this repository's static benchmark (Phase 1: 2D descriptors + conformer ensembles) into the **fourth dimension**. Five next-generation medicinal-chemistry targets — a CRBN-imide molecular glue, an ADC Val-Cit-PAB linker, a bicyclic disulfide peptidomimetic, an allosteric covalent inhibitor core, and a VHL-mimetic proline hybrid — were subjected to:

1. **Stage 1** — full descriptor profiling (MW, cLogP, TPSA, Fsp³, RotB, HBD/HBA) plus synthesizability (SAScore), drug-likeness (QED) and PAINS A/B/C alerts. **All five targets are PAINS-clean.**
2. **Stage 2** — 36-point *relaxed* torsional scans (0→360°, 10° grid, MMFF94 with ±5° flat-bottom torsion constraints and full re-minimization at every grid point), yielding rotational barriers ΔE‡ from **7.8 to 72.8 kcal/mol**.
3. **Stage 3** — a **200 ps production MD** trajectory of the most flexible target (**T02**, 14 rotatable bonds) in GBSA/OBC2 implicit solvent at 300 K, with RMSD, radius of gyration, energy and intramolecular H-bond analyses (200 frames, ~970 steps/s on CPU).

**Headline findings:**

| Finding | Number | Interpretation |
|---|---|---|
| Lowest rotational barrier (T01, glutarimide hinge) | **7.8 kcal/mol** | freely equilibrating at 300 K (t½ ≈ 50 ns) |
| Highest acyclic barrier (T02, Val-Cit amide/urethane hinge) | **21.4 kcal/mol** | conformationally locked on the second timescale (t½ ≈ 7 min) |
| Ring-strain probe (T03, disulfide dihedral) | **72.8 kcal/mol** | effectively a *covalent* constraint — the disulfide is not rotatable at any physiological temperature |
| T02 MD stability | RMSD plateau 2.99 Å, Rg 7.60 ± 0.38 Å, PE flat at −6280 kJ/mol | folded, thermally stable hinge without runaway unfolding in 200 ps |
| T02 intramolecular H-bonds | **0 per frame** (audited) | solvent-exposed by design — correct behavior for a cathepsin-B-cleavable linker |

---

## 2. Targets & Stage-1 Medicinal-Chemistry Descriptors

| ID | Target | SMILES (input) | MW (Da) | cLogP | TPSA (Å²) | Fsp³ | RotB | HBD/HBA | SAScore | QED | PAINS |
|---|---|---|---|---|---|---|---|---|---|---|---|
| T01 | Molecular glue degrader (CRBN imide mimic) | `O=C1NC(=O)C(N2C(=O)c3c(N4CCN(c5cccc(Cl)c5)CC4)c(F)ccc3C2=O)CC1` | 470.89 | 2.21 | 90.03 | 0.304 | 3 | 1/7 | 3.10 | 0.692 | 0 |
| T02 | ADC cathepsin-B linker (Val-Cit-PAB mimic) | `CC(C)[C@@H](NC(=O)[C@H](CCCNC(=O)N)NC(=O)OCc1ccc(NC(=O)OCC2=CCN(CC2)C(=O)C)cc1)C(=O)O` | 604.66 | 1.68 | 218.49 | 0.500 | 14 | 6/11 | 3.59 | 0.134 | 0 |
| T03 | Bicyclic disulfide-constrained peptidomimetic | `O=C1N[C@H]2CSSC[C@@H]3NC(=O)[C@H](CC(=O)N3)NC(=O)[C@H](CSSC2)NC1=O` | 465.60 | −1.81 | 145.50 | 0.643 | 0 | 6/8 | 6.57 | 0.206 | 0 |
| T04 | Allosteric covalent inhibitor core (switch-II pocket) | `C=CC(=O)N1CCN(CC1)c2nc(Nc3ccc(F)c(Cl)c3)nc4c2c(C#C)c(c5ccccc5)n4C` | 514.99 | 4.99 | 66.29 | 0.179 | 5 | 1/8 | 2.94 | 0.300 | 0 |
| T05 | Non-thalidomide E3 binder (VHL-mimetic proline hybrid) | `CC(C)(C)[C@H](NC(=O)[C@@H]1C[C@@H](O)CN1C(=O)c2ccc(c3nccs3)cc2)C(=O)NCc4ccc(C#N)cc4` | 545.67 | 3.10 | 135.42 | 0.345 | 7 | 4/8 | 3.60 | 0.418 | 0 |

**Reading the descriptor space.** The five targets deliberately span orthogonal design axes: T01 is a compact, Lipinski-compliant glue (QED 0.69 — the only "classically drug-like" member); T02 sits **far outside** Ro5 (TPSA 218 Å², MW 605, 14 RotB, QED 0.13) exactly as ADC linkers must; T03 is a zero-RotB, H-bond-dense bicyclic peptide (Fsp³ 0.64, cLogP −1.8 — beyond-Ro5 (bRo5) chameleon territory); T04 is a lipophilic, flat covalent core (Fsp³ 0.18, cLogP 5.0); T05 is a tenuously Ro5-compliant E3 ligand. The radar fingerprint of this spread is shown in `./figures_phase2/fig3_medchem_radar.png`.

**PAINS:** zero alerts across all three filter families for all five molecules — none of these designs survives by luck; they are *intended* to be clean.

---

## 3. Relaxed Torsional Barrier Scans (Stage 2)

![torsion scans](./figures_phase2/fig1_torsion_scans.png)

**Protocol.** For each target, the dominant "hinge" torsion was selected automatically (acyclic single bond maximizing fragment-split size; for T03, which has *zero* rotatable bonds, the disulfide C–S–S–C dihedral was used as a ring-strain probe). Each of 36 grid angles θ ∈ {0°, 10°, …, 350°} was enforced with an MMFF94 flat-bottom torsion constraint (window θ±5°, k = 1000), and the *whole molecule* was re-minimized at every point — a **relaxed** scan, not a rigid one. Energies are reported relative to the global scan minimum.

| ID | Torsion (atoms) | Scan mode | θ(E_min) | θ(barrier) | ΔE‡ (kcal/mol) | Eyring t½ @300 K* |
|---|---|---|---|---|---|---|
| T01 | 9-10-11-12 | hinge (rotatable) | 50° | 200° | **7.80** | ~50 ns |
| T02 | 15-16-18-19 | hinge (rotatable) | 220° | 80° | **21.35** | ~7 min |
| T03 | 21-22-23-24 | disulfide S-S strain probe | 260° | 220° | **72.82** | effectively permanent |
| T04 | 6-7-10-11 | hinge (rotatable) | 260° | 40° | **10.53** | ~5 µs |
| T05 | 1-4-5-6 | hinge (rotatable) | 110° | 10° | **17.30** | ~0.5 s |

\* Crude single-barrier Eyring estimates, k = (k_BT/h)·exp(−ΔE‡/RT) with RT = 0.596 kcal/mol; real kinetics are multi-dimensional and entropy-laden. Grid scans also *lower-bound* true maxima (peaks may fall between 10° samples).

### 3.1 Rotational barriers and bioactive-conformation selection

- **T01 (7.8 kcal/mol)** — the glutarimide-piperidone hinge of the CRBN mimic equilibrates *nanoseconds* after binding or synthesis. Its bioactive conformation is selected enthalpically by the protein (CRBN surface) rather than pre-payable by the molecule; no atropisomeric or conformational-isomer liability. The scan shows a broad, shallow double-minimum basin (50° ↔ 260° separated by <8 kcal/mol), i.e., genuine conformational freedom that a docking engine must sample rather than assume.
- **T04 (10.5 kcal/mol)** and **T05 (17.3 kcal/mol)** — intermediate regimes: microsecond-to-second half-lives. For the switch-II-pocket covalent core, the barrier separates the acrylamide **s-cis/s-trans** populations that vector the β-carbon toward the catalytic cysteine; a 10.5 kcal/mol barrier means both rotamers populate on folding timescales but not within a single binding encounter — conformational selection, not induced fit, will dominate recognition. For the VHL-hybrid, the 17.3 kcal/mol barrier localizes the tert-leucine–proline amide geometry; only ~1 in 10⁶ molecules at any instant sits at the barrier geometry, so the *pre-organized* rotamer determines which crystal conformer is observable.
- **T02 (21.4 kcal/mol)** — the single torsion that dominates the Val-Cit-PAB linker's folding (the Cit-urethane hinge) is **locked for minutes**. This is a designed feature: the linker must hold the payload rigid through plasma transit, then be cleaved in one defined place. A 2D model that treats all 14 rotatable bonds as equivalent samplers fundamentally cannot express "one locked hinge among fourteen free bonds."
- **T03 (72.8 kcal/mol)** — not a rotation at all: the C–S–S–C dihedral in the bicyclic dual-disulfide framework would need to *break the ring* to rotate. The 72.8 kcal/mol figure quantifies the conformational pre-organization: the peptidomimetic is effectively a rigid body whose H-bond donors/acceptors are pinned in space — the design intent of disulfide bicyclization, now with a number attached.

### 3.2 What this buys the benchmark

Every ΔE‡ above is a **label that no 2D graph model can see** but every downstream consumer (docking, free-energy, generative design) implicitly depends on. Section 5 turns this into the AI-facing argument.

---

## 4. OpenMM Molecular Dynamics of T02 (Stage 3)

![MD trajectories](./figures_phase2/fig2_md_trajectories.png)

### 4.1 System construction (the hard part, documented honestly)

The OpenFF ecosystem on PyPI is currently in a yanked-transition state, and AM1-BCC charge assignment via NAGL fails in this offline environment (`BadFileSuffixError` — model files not fetchable). The working recipe, encoded in `export_openff_system.py` and run in a dedicated conda env (`phase2ff`: openff-toolkit 0.19.0, openff-interchange 0.5.4):

1. ETKDGv3 + MMFF94 starting geometry (RDKit, main env).
2. `Interchange.from_smirnoff` with **Sage 2.1.0** valence + vdW after *deregistering* the Electrostatics/ToolkitAM1BCC/ChargeIncrementModel handlers (NAGL unavailable).
3. **Charge surgery:** RDKit MMFF94 partial charges (validated: Σq = formal charge = 0 within 10⁻⁶ e) written into the serialized `NonbondedForce`, with 1-4 Coulomb exceptions scaled ×1/1.2 to match the MMFF convention.

The result is a hybrid MMFF94-charge / Sage-valence system — a documented, reproducible approximation appropriate for a *comparative benchmark*, not for absolute binding free energies.

### 4.2 Simulation protocol

| Parameter | Value |
|---|---|
| Target | T02 (chosen automatically: max RotB = 14) |
| Solvent model | GBSA/OBC2 implicit (Bondi radii + OBC screening scales) |
| Integrator | LangevinMiddle, T = 300 K, γ = 1 ps⁻¹, Δt = 2 fs |
| Particles / constraints | 83 / 40 (X–H bond constraints) |
| Equilibration | 5 000 steps (10 ps) |
| Production | 100 000 steps = **200 ps**, DCD frame every 500 steps (200 frames) |
| Throughput | 969.7 steps/s (CPU platform) |
| Temperature realized | 301.5 K mean |

### 4.3 Results and interpretation

- **Thermal stability.** Potential energy is flat at −6280 ± ~40 kJ/mol after equilibration with no drift over 200 ps; total energy (panel C of fig 2) conserves as expected for the Langevin middle scheme. The system is thermally stable — no parameterization pathologies, no non-physical strain growth.
- **Global fold.** Backbone RMSD plateaus at **2.99 Å** (max excursion 4.81 Å) relative to the MMFF-minimized start. For a 605 Da, 14-RotB linker this is the signature of a *hinge-breathing* polymer, not an unfolding event: RMSD oscillates between extended (~4.8 Å) and compact (~2.1 Å) states on ~20–40 ps timescales.
- **Compactness.** Radius of gyration **Rg = 7.60 ± 0.38 Å** (range 5.83–8.06 Å). The ±0.38 Å breathing around a 7.6 Å mean is exactly the extended ⇄ collapsed interconversion an ADC linker must support (compact in plasma to evade hydrolysis, extended in the cathepsin pocket for cleavage access).
- **Intramolecular H-bonds: 0.000 per frame.** This null result was *audited*, not assumed: every heavy-atom D···A pair below 3.5 Å in the trajectory is either a geminal amide contact (N–H···O=C on the same residue, 2-bond separation, D-H···A angles 23–27° — chemically excluded) or fails the 120° angular criterion (closest genuine candidate: 2.63 Å at 104.7°). The linker is solvent-exposed by design; in implicit solvent it has no partner to compensate desolvation, so the census is expected to be zero. **Caveat:** frame sampling every 1 ps can miss sub-picosecond flickering H-bonds; a hydrogen-bond *lifetime* study would need 10–20 fs output.
- **Physical reasonableness checklist.** 40/40 constraints satisfied throughout; no particle loss; DCD frame count = 200/200; final-frame PDB visually compact (`results_phase2/T02_final.pdb`).

---

## 5. Implications for Generative AI & 3D Equivariant Force Fields

**For 2D/fragment generative models.** Phase 1 showed these molecules break 2D GNN assumptions (macrocycles, zwitterions, atropisomerism). Phase 2 sharpens the knife:

- **Property priors mis-rank designed space.** QED scores T02 at 0.134 and T03 at 0.206 — "poor" — yet both are *archetypes of clinically validated design patterns* (approved ADC linkers; disulfide-bicyclic peptide drugs). A generative model maximizing a QED-like prior will actively *avoid* correct solutions in the bRo5 regime. Torsion barriers and Rg-fluctuation signatures are orthogonal coordinates that the same models ignore entirely.
- **One number, fourteen bonds.** T02's ΔE‡ = 21.4 kcal/mol on a *specific* atom quartet is a per-bond, spatially-resolved label. 2D models see RotB = 14 (a count). The count and the barrier are not correlated across this set (T01: 3 RotB, 7.8; T02: 14 RotB, 21.4; T03: 0 RotB, 72.8) — so barrier prediction *requires* 3D structure, and barrier *data* is what this pipeline mass-produces (5 relaxed scans = 180 constrained minimizations per run).

**For 3D equivariant force fields (MACE / Allegro / e3nn family).**

- **Domain of applicability.** MACE-MP-0/Allegro-class models are trained chiefly on neutral, unstrained, drug-like or inorganic periodic systems. This suite's edges — a 72.8 kcal/mol strained disulfide dihedral, T03's Fsp³ = 0.64 peptide density, T04's acrylamide warhead (reactive electronic structure) — are outside that distribution. The relaxed-scan curves here are ready-made **stress-test inputs**: run any learned FF on the same 36-point protocol and compare curve shape, barrier height and minimum-angle placement against the MMFF94 reference (agreeing with MMFF94 is not the point; *failing gracefully* is).
- **A curriculum, not a ground truth.** MMFF94 labels are themselves an approximation (no polarization, no explicit solvent). The productive use is **transfer-learning curriculum**: MMFF-relaxed scans as cheap pre-training PES data for TorsionNet-style heads, with the final fine-tune on DFT (ωB97X-D/def2-TZVP would be the natural next phase) on a *curated subset* — the five scans here tell you exactly which subset matters (the barrier regions, ±30° around each maximum).
- **MD as the integration test.** The T02 trajectory demonstrates the full simulation stack working on a bRo5 molecule. Equivariant FFs promise MD at DFT quality; the benchmark's contribution is the *diagnostic protocol* — 200 ps at 300 K, then check energy drift, RMSD plateau behavior, Rg breathing amplitude, and H-bond census. Those four numbers are how a learned FF's failure modes become visible before it is trusted.

---

## 6. Reproducibility, Artifacts & Limitations

### Commands

```bash
# full pipeline (descriptors → scans → MD → figures), ~6-8 min on CPU
python run_heavy_dynamics_benchmark.py

# regenerate figures only from saved results (no re-run of MD)
python run_heavy_dynamics_benchmark.py --fig_only

# knobs
--scan_points 36 --equil_steps 5000 --md_steps 100000 --report_interval 500
--md_target auto|T02|... --skip_scan --skip_md --auto_shutdown --shutdown_delay 60
```

### Artifacts (`results_phase2/`)

| File | Content |
|---|---|
| `phase2_results.json` | complete machine-readable record (all stages) |
| `T0*_torsion_scan.csv` | 36-point angle/energy tables |
| `T0*_scan_min.sdf` | lowest-energy scan conformer per target |
| `T02_system.xml` | OpenMM-serialized system (Sage valence + MMFF94 charges) |
| `T02_start.pdb` / `T02_final.pdb` | MD endpoints |
| `T02_trajectory.dcd` | 200-frame production trajectory |
| `T02_md_metrics.csv` | per-frame RMSD/Rg/energies/T |

### Environment matrix

- Main env: Python 3.14.6, RDKit 2026.03.5, matplotlib, ASE 3.29, NumPy.
- Parameterization env: conda `phase2ff` (Python 3.12, openff-toolkit 0.19.0, openff-interchange 0.5.4) — invoked by the script via `--ff_env` subprocess.
- MD: OpenMM 8.6, CPU platform.

### Limitations (read before citing numbers)

1. **MMFF94 torsion profiles** are classical, gas-phase, no polarization; barriers are semiquantitative (±2–3 kcal/mol typical vs. DFT for drug-like torsions; worse for strained/electronic edge cases like T03/T04).
2. **Charge hybridization** (MMFF94 charges on Sage valence) sacrifices internal consistency; absolute PE values (−6280 kJ/mol) are not comparable across force fields — only *trends within this trajectory*.
3. **Implicit solvent** (OBC2) lacks explicit water competition and desolvation structure; the H-bond census is correspondingly conservative.
4. **200 ps × 1 molecule** samples hinge breathing but not rare events (no full cyclization/decyclization of the piperidine tail was observed; expect multi-ns timescales).
5. Eyring half-lives in §3 are single-barrier, single-pathway estimates — treat as order-of-magnitude.
6. The 10° scan grid lower-bounds barrier heights; no transition-state refinement (e.g., eigenvector-following) was performed.

---

*Phase 2 of the ai4chem-complex-scaffolds-benchmark. Phase 1 (10-molecule conformer/descriptor suite) is documented in `BENCHMARK_REPORT_EN.md`.*
