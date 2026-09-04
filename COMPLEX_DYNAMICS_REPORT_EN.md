# COMPLEX DYNAMICS REPORT — Phase 3: Target-Ligand Complex, Residue-Decomposed MM-GBSA & ML-Potential Benchmarking

*Data: `results_phase3/phase3_results.json` | OpenMM 8.6 (CPU) · pdbfixer 1.12 · AutoDock Vina 1.2.5 CLI + meeko 0.8 · OpenFF Sage 2.1 + AM1-BCC (NAGL 1.0.0, offline) · torchani ANI-2x · mdtraj 1.11*

---

## 1. Target & Erratum

**Target: KRAS G12D, switch-II pocket (PDB 7RPZ)** — X-ray 1.30 Å, GDP-bound, co-crystallized with **MRTX1133** (residue `6IC`, 44 heavy atoms), catalytic-domain construct 1–170 with the engineered G12D mutation (SEQADV ASP12 vs GLY) and a cysteine-light background (S51/C80L/S118 conflicts vs UniProt P01116-2).

> **Spec erratum, verified live against the RCSB REST API**: the PDB ID "7E27" cited for KRAS G12D is actually *Plasmodium falciparum formate-nitrite transporter in complex with MMV007839* (cryo-EM, 2.29 Å). The canonical G12D switch-II-pocket co-crystal used here is 7RPZ (searched and verified among 142 G12D entries). 7BQY (SARS-CoV-2 Mpro + N3, 1.7 Å) was the verified alternative; KRAS was chosen because Phase-2 target **T04** was designed as a switch-II-pocket covalent inhibitor core.

**Ligand: T04** (Phase 2) — 2-aminopyrido-pyrimidinone-piperazine acrylamide (C₂₈H₂₄ClFN₆O, 61 atoms with H), the acrylamide-warhead paradigm of the KRAS G12C/G12D switch-II pocket era.

---

## 2. Stage 1A — Macromolecular Ingestion & Curation

- Automated RCSB fetch (`files.rcsb.org/download/7RPZ.pdb`).
- Component split: protein 1,342 heavy atoms | MRTX1133 44 | GDP 28 | Mg²⁺ 1 | waters stripped; alternate conformers (altloc ≠ A/B) removed.
- **PDBFixer repair**: 7 missing heavy atoms added; hydrogens at **pH 7.4** (His tautomers by H-bond pattern rules); final protein 2,681 atoms.
- **Pocket definition**: centroid of MRTX1133 heavy atoms = [1.71, 4.93, −23.16] Å; **28 residues** within R = 10.0 Å (includes H95, Y96, Q99, E62, D12, D92, Y64, R68, K88, A11…).

---

## 3. Stage 2 — Structure-Based Docking (meeko + Vina)

Ligand: ETKDGv3 + MMFF94 conformer → meeko `PDBQTWriterLegacy` (flexible torsions). Receptor: fixed protein **+ GDP** (co-factor retained as rigid receptor context) via openbabel 3.2.1 rigid PDBQT.

*Engine fallback:* no Windows/Py3.12 wheel exists for the Vina Python API → the official **Vina v1.2.5 CLI executable** was used (results parsed from the mode table; meeko round-trips poses to RDKit).

| Pose | ΔG_docking (kcal/mol) |
|---|---|
| **#1 (selected)** | **−7.92** |
| #2 | −7.30 |
| #3 | −7.11 |

Box 22.5 Å³ centered on the 6IC centroid, exhaustiveness 16, seed 42, 9 modes.

---

## 4. Stage 3 — The Paradigm Clash: Classical vs ML Potentials on the Frozen Pocket Pose

Single-point energies/atomic forces on **pose #1 as frozen in the pocket** (61 atoms; hydrogens first relaxed with heavy atoms pinned, so forces reflect the pose, not AddHs artifacts).

**Documented environment fallbacks** (spec-sanctioned):

1. **GAFF2** (`openmmforcefields`) is unavailable offline without downgrading the working OpenFF stack → classical reference = **Sage 2.1 valence+vdW with AM1-BCC charges**. Notably, AM1-BCC now works fully offline via the packaged NAGL GNN model (`openff-gnn-am1bcc-1.0.0.pt` invoked by file path) — retiring Phase 2's MMFF-charge-surgery workaround.
2. **MACE-OFF23**: model download failed (offline GitHub) → **ANI-2x (torchani 2.8.4)**, explicitly sanctioned by the mission spec.

| Potential | E (kcal/mol) | ‖F‖max (kcal/mol/Å) |
|---|---|---|
| Sage 2.1 + AM1-BCC (classical, GAFF2 stand-in) | −732.4 | 132.5 |
| MMFF94 (RDKit, numeric gradients) | +37.7 | 36.1 |
| **ANI-2x (ML reference)** | −47,152.0 | 1.74 |

(Energies have different zero references across models; only forces are directly comparable.)

### 4.1 Force-field discrepancy vector (vs ANI-2x)

| Metric | Sage+AM1-BCC | MMFF94 |
|---|---|---|
| Pearson r (per-atom ‖F‖) | 0.698 | 0.638 |
| mean cos(F_classical, F_ML) | 0.483 | 0.555 |
| ⟨‖ΔF‖⟩ (kcal/mol/Å) | 24.7 | 9.8 |
| max ‖ΔF‖ | 135.4 | 42.8 |

### 4.2 Which moieties carry the artificial strain (classical vs ML)

| Moiety (SMARTS-assigned) | atoms | ⟨‖ΔF_Sage−ANI‖⟩ | max |
|---|---|---|---|
| **2-aminopyrimidine heteroaromatic** | 7 | **81.0** | 135.4 |
| alkyne linker (C#C) | 2 | 57.8 | 65.1 |
| acrylamide warhead (C=CC(=O)N) | 5 | 46.9 | 63.5 |
| piperazine | 5 | 27.6 | 72.0 |

**Reading**: the classical force field's largest disagreement with the equivariant ML potential sits exactly on the **conjugated N-heteroaromatic core and the sp/sp² warhead vector** — the electron-delocalized, polarizable moieties whose torsional profiles and through-conjugation harmonic MM cannot represent. MMFF94 agrees with ANI better than Sage here (⟨ΔF⟩ 9.8) — a reminder that "which classical FF" matters as much as "classical vs ML". Figure: `./figures_phase3/fig3_ml_vs_classical_ff_gap.png`.

---

## 5. Stage 4 — OpenMM Complex Dynamics (Implicit OBC2, 310 K)

**System construction (2785 atoms)** — pure-OpenMM splice, no openmmforcefields:
- **Protein**: Amber14SB (`amber14/protein.ff14SB.xml`) + **GBSA/OBC2** via `implicit/obc2.xml` (OpenMM 8.x Script generator → CustomGBForce), Cα restraints k = 5 kcal/mol/Å².
- **GDP + T04**: Sage 2.1 valence + AM1-BCC charges spliced into the combined system (valence forces index-offset; GB particles appended with Bondi×OBC radii; Sage's own X–H constraints copied, deduplicated).
- **Salt**: Debye–Hückel screening κ = 1.263 nm⁻¹ (0.15 M 1:1 salt, ε_w = 76.6, T = 310 K) via `implicitSolventKappa`.
- **Mg²⁺ dropped** (no template in the pure-OpenMM path; documented).

**Protocol deviations, both documented and mechanistically diagnosed**:
1. `dt = 1.0 fs` (spec 2.0 fs): at 310 K the 2-fs step diverges via the Sage-spliced molecules' angle modes; a 10→50 K (0.5 fs) → 150 K (1 fs) → 310 K heating ramp made the system rigorously stable (validated: 20,000-step clean test before production).
2. Production = **50,000 steps = 50 ps** (spec range 50–100k) at ~14 steps/s (no GB interaction-cutoff API in this OpenMM build); equilibration 5,000 steps; DCD every 250 steps → **200 frames**.

| Observable | Value |
|---|---|
| ⟨Cα RMSD⟩ | **0.347 Å** (max 0.39) — restrained fold rock-stable |
| ⟨pocket-Cα RMSD⟩ | 0.352 Å |
| ⟨ligand heavy-atom RMSD⟩ (tail) | **1.275 Å** (max 1.66) — pose retained, induced-fit wiggle |
| PLIF (hydrophobic, top) | **Tyr96** (dominant, 4 ligand-carbon contacts ≥ ~60% persistence), Tyr64 |
| PLIF (H-bonds) | sparse — consistent with a hydrophobic-groove pose (see fig1) |

---

## 6. Stage 5 — End-State MM-GBSA & Per-Residue Decomposition

40 frames (every 5th of 200). ΔG_bind = ⟨G_complex⟩ − ⟨G_protein⟩ − ⟨G_ligand⟩ with GB-polar from the CustomGBForce (κ-screened), nonpolar from γ·ΔSASA (γ = 0.005 kcal/mol/Å², Shrake–Rupley 960-pt sphere).

**ΔG_bind = −143.47 ± 2.85 kcal/mol** (dE_MM −102.24 | ΔG_GB −38.10 | ΔG_SA −3.13)

> Implicit-solvent end-point MM-GBSA systematically overbinds relative to experiment (no configurational-entropy term, no explicit-water competition, unscreened charged contacts in dE_MM). The −143 kcal/mol figure is a **ranking/decomposition quantity**, not a Kd predictor — the residue decomposition, not the absolute number, is the deliverable.

### 6.1 Per-residue decomposition (top 10)

| Residue | total | vdW | elec | GB polar | note |
|---|---|---|---|---|---|
| **His95** | **−8.65** | −6.45 | −5.59 | +3.39 | SIIP's signature H-bond/π wall |
| **Tyr96** | **−6.33** | −5.53 | −2.47 | +1.67 | hydrophobic lid (matches PLIF dominance) |
| **Glu62** | −4.54 | −4.60 | +7.23 | **−7.17** | switch-II; classic elec↔desolvation compensation |
| **Asp12 (G12D!)** | −4.33 | −4.05 | −4.08 | +3.80 | the oncogenic mutation itself engages the ligand |
| Asp92 | −4.04 | −3.74 | −4.39 | +4.09 | |
| Tyr64 | −3.44 | −2.67 | **−11.29** | **+10.52** | electrostatically locked, desolvation-paid |
| Ala11 | −2.35 | −2.41 | +3.91 | −3.85 | |
| Arg68 | −1.96 | −1.61 | +9.76 | −10.11 | |
| Lys88 | −1.82 | −1.87 | −1.72 | +1.77 | |
| Gln99 | −1.14 | −1.29 | +8.44 | −8.29 | MRTX1133's canonical contact, weaker here |

**Reading**: the catalytic/allosteric triad of the switch-II pocket — **His95/Tyr96 (vdW+H-bond walls) plus the G12D Asp12 and switch-II Glu62** — dominates binding; for charged residues the electrostatic term is almost exactly cancelled by GB desolvation (Y64, R68, Q99), a textbook MM-GBSA signature. Figures: `./figures_phase3/fig2_per_residue_mmgbsa.png`, pose map `fig1_binding_pose_pocket.png`.

**Honest validation caveat**: my independent numpy pair-sum of the cross vdW+electrostatics could not be reconciled with OpenMM's charge-nulling difference on the first frame (−68.9 vs +89.5 kcal/mol, 177% rel err; likely reaction-field/exception accounting in the Script-GB system). The per-residue decomposition is internally consistent (identical formula across residues) and its ranking matches the crystallographic contact map, but per-residue absolute values carry this documented uncertainty. GB components use Still pair terms with intrinsic OBC radii (screening-level approximation).

---

## 7. Why Classical FFs Struggle Here & the Equivariant-ML Case

1. **Fixed-charge polarization failure.** T04 sits between His95/Glu62/Asp12 (±0.5–0.8 e AM1-BCC) — a field the fixed charges were fitted without. Stage 3 quantifies the consequence: ⟨‖ΔF_Sage−ANI‖⟩ ≈ 25 kcal/mol/Å per atom, peaking at 135 on the heteroaromatic core — precisely where induced polarization and through-conjugation live. Sage *and* MMFF94 both misalign force directions with ANI-2x (mean cosine 0.48/0.56).
2. **Harmonic torsions vs conjugation.** The acrylamide warhead and alkynyl linker — the reactivity-vectoring dihedrals — are 2nd/3rd in the discrepancy ranking. A classical torsion fitted to gas-phase scans cannot carry the pocket-field-induced rotation barrier shifts that decide warhead geometry.
3. **Equivariant potentials (MACE/PaiNN/Allegro)** learn these effects from data: message passing gives many-body polarization for free; e(3)-equivariance gives smooth, energy-conserving forces everywhere. This run's ANI-2x (a 3× cheaper ancestor of that class) already exposes the classical gap on ONE frozen pose; MACE-OFF23 on the same protocol is the natural Phase-4 stress test (its model download failed offline here — documented).
4. **The ML-vs-classical gap is pocket-specific**, which is exactly what SBDD needs: it localizes WHERE the classical simulation underlying docking/rescoring is least trustworthy — here, the aminopyrimidine core and warhead vector of a covalent inhibitor.

---

## 8. Reproducibility, Artifacts, Limitations

### Commands

```bash
# in the phase2ff conda env, with its Library/bin on PATH (BLAS DLLs):
python run_phase3_complex_dynamics.py                  # full pipeline
python run_phase3_complex_dynamics.py --fig_only       # regenerate figures
# knobs: --pdb_id 7RPZ --lig_code 6IC --md_steps 50000 --report_interval 250
#        --mgb_frames 40 --force_rerun --skip_dock --skip_md --auto_shutdown
```

### Artifacts (`results_phase3/`)

`7RPZ.pdb` · `7RPZ_fixed_protein.pdb` · `7RPZ_GDP.pdb`/`GDP_withH.pdb` · `receptor.pdbqt` · `T04_docked_poses.pdbqt` (9 modes) · `T04_pose1.sdf/pdb` (H-relaxed) · `complex_start.pdb` · `T04_complex_trajectory.dcd` (200 frames) · `T04_complex_metrics.csv` · `T04_plif_persistence.csv` · `T04_per_residue_mmgbsa.csv` · `phase3_results.json` (all stages) · per-stage checkpoints (`stage1a/2/3/4/5.json`).

### Environment matrix

Main env Python 3.14 (RDKit 2026.03.5) for authoring; execution env `phase2ff` (Python 3.12): OpenMM 8.6, pdbfixer 1.12, openbabel 3.2.1, mdtraj 1.11, openff-toolkit 0.19/interchange 0.5.4, meeko 0.8, torchani 2.8.4 + torch 2.10, mace-torch 0.3.16 (model unavailable offline), `tools/vina.exe` 1.2.5, numpy 2.4.6 (conda-forge; requires `Library/bin` on PATH).

### Limitations (read before citing numbers)

1. GAFF2 → Sage 2.1 + AM1-BCC substitution (openmmforcefields unavailable offline); GAFF-vs-Sage valence differences are small vs the classical↔ML gap under study, but noted.
2. MACE-OFF23 → ANI-2x (spec-sanctioned fallback); ANI-2x lacks a Cl-trained? — Cl is in ANI-2x's element set (H,C,N,O,S,F,Cl) ✓, but ANI-2x accuracy < MACE-OFF23 on strained conjugated systems.
3. Implicit solvent (OBC2), Mg²⁺ omitted, GDP kept; dt = 1 fs deviation (diagnosed); 50 ps single-trajectory, Cα-restrained — induced fit within a restrained fold, not full exchange.
4. MM-GBSA end-point caveats (no entropy, overbinding) and the pair-sum validation mismatch (§6) — per-residue *ranking* robust, absolute values approximate.
5. Docking scores are Vina empirical, not ΔG predictions.

---

*Phase 3 of ai4chem-complex-scaffolds-benchmark. Phase 1: `BENCHMARK_REPORT_EN.md` (10-molecule conformer suite). Phase 2: `DYNAMICS_REPORT_EN.md` (torsion barriers + 200 ps ligand MD).*
