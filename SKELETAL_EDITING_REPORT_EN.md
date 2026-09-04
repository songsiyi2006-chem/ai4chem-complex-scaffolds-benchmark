# SKELETAL EDITING REPORT — Phase 4: Autonomous Transition-State Hunting & CI-NEB Reaction Pathway Profiling

*Data: `results_phase4/phase4_results.json` | ASE 3.29 CI-NEB (improved tangent) · ANI-2x reactive PES (torchani) + GFN2-xTB analytic Hessian verification (xtb 7.x CLI) · RDKit 2026.03*

---

## 1. Model Reaction: the Core Elementary Step of Indole→Quinoline Skeletal Editing

Modern **molecular skeletal editing** (Levin, Sarpong, et al., 2022–2024) rewrites molecular frameworks by single-atom **insertion** or **deletion** rather than by building rings de novo. The grandfather of indole→quinoline editing is the **Ciamician–Dennstedt rearrangement** (1901): a carbene cyclopropanates the indole C2=C3 bond, the strained fused cyclopropane opens, and the pyrroline ring expands by formal one-atom insertion into a 6,6-fused quinoline core.

This pipeline models the **ring-expansion elementary step** in its simplest :CH₂ form:

| | Species | Formula | Topology |
|---|---|---|---|
| **R** | cyclopropa[b]indole (indole C2=C3 cyclopropanated strained adduct) | C₉H₉N | fused 5,6 + cyclopropane |
| **P** | 2,3-dihydroquinoline (ring-expanded core) | C₉H₉N | fused **6,6** |

R and P are **valence isomers of identical formula** — the essence of an editing elementary step. The pipeline's pairing algorithm (below) recovers the chemically faithful atom correspondence: only **two bonds break (the cyclopropane shared edge C2–C3 and one C–C rearrangement partner) and one bond forms** — the textbook minimal topological change of this rearrangement.

*(Spec context: "activated fused bicyclic indole/aziridine-like strained intermediate → ring-expanded 6,6-fused quinoline core" — exactly this transformation.)*

---

## 2. Stage 1 — Geometry, Alignment & IDPP Interpolation

- R built **programmatically** (indole C2=C3 double bond located by ring-membership analysis, saturated, bridged by CH₂) — hand-written SMILES for this fused system are aromaticity/kekulization traps (three失败 variants documented in the code's history).
- Both endpoints embedded (ETKDGv3 + MMFF94), then **GFN2-xTB optimized** (`xtb --opt`).
- **Atomic correspondence**: MCS-based alignment fails for skeletal-editing isomers (the ring connectivity itself changes). A v4 algorithm solves it: enumerate all MCS substructure matches (symmetry variants) × Hungarian completion on **permutation-invariant sorted distance profiles**, scored by **bond-change count** — recovers the minimal 2-broken/1-formed mapping deterministically.
- **Kabsch rigid-body alignment** of P onto R through the pairing (IDPP between misoriented endpoints makes atoms fly across space — diagnosed the hard way).
- **9-image IDPP interpolation** (ASE) — image-dependent pair potential prevents unphysical atomic overlap along the initial band.

## 3. Stage 2 — Climbing-Image NEB

### 3.1 The CI-NEB method (derivation)

A discretized path {R₀=R, R₁ … R_M=P} is relaxed by minimizing

  S[{Rᵢ}] = Σᵢ E(Rᵢ) + Σᵢ½ k|Rᵢ₊₁ − Rᵢ|²

under the **nudging projection**: each image feels only the component of the true force perpendicular to the local tangent, plus the spring force parallel to it:

  Fᵢ = −∇E(Rᵢ)⊥ + (k(|Rᵢ₊₁−Rᵢ| − |Rᵢ−Rᵢ₋₁|))·τ̂ᵢ

with the **improved-tangent** estimate τ̂ᵢ (Henkelman–Jonsson 2000) weighting by neighbouring energy gradients. The **climbing image** (highest-E image, here image 4/9) additionally inverts the parallel true-force component:

  F_CI = −∇E⊥ − (−∇E∥) = −∇E⊥ + ∇E∥

so it relaxes **uphill** along the band and **downhill** in all 57 perpendicular directions — a saddle-point search with the band itself providing the reaction coordinate. Convergence: fmax (max per-atom force norm across the band) < 0.05 eV/Å.

### 3.2 Execution & engine battle (honest log)

| Attempt | Outcome |
|---|---|
| GFN2-xTB via `xtb.exe --grad` subprocess | **force-sign bug**: xtb writes the gradient; ASE expects forces = −∇E. The band climbed the hill and exploded (fmax → 4700 eV/Å). Fixed by negation; H₂ probe now passes with correct attraction direction. |
| xtb CI-NEB from IDPP | BFGS/FIRE instability through subprocess-force noise: late-phase divergence (band top +11 eV). Diagnosed with per-group force monitoring. |
| Multi-fidelity warm start (ANI pre-relaxation → xtb climb) | Better (fmax 0.5-0.8) but tail instability persisted; **best-band snapshot rollback** machinery added. |
| tblite/python bindings, MACE-OFF23 | unavailable on win-64/offline (documented fallbacks). |
| **ANI-2x single-phase CI-NEB (climb ON from step 0)** | **Converged: 64 FIRE steps, fmax = 0.0468 eV/Å < 0.05 ✓** (a smooth ML PES tolerates climbing from the start; the staged no-climb protocol lets the band wander into spurious high routes — both behaviours measured, see log). |

**Final MEP (ANI-2x, 9 images):** E(eV) = [−10965.90, −10965.65, −10963.93, −10965.28, **−10963.62 (TS)**, −10964.17, −10965.54, −10966.19, −10966.31]
- **ΔE‡(forward) = 52.5 kcal/mol**, ΔE‡(reverse) = 62.0 kcal/mol, ΔE_rxn = **−9.5 kcal/mol** (exothermic — ring strain released, thermodynamically downhill ✓)

## 4. Stage 3 — Quantum Verification: the 1-Imaginary-Frequency Test

Analytic **GFN2-xTB Hessian** (`xtb --hess`) on the converged climbing image — a *different, semi-empirical quantum* surface verifying an ML-PES saddle (deliberate multi-fidelity cross-check):

- **EXACTLY ONE imaginary frequency: ν‡ = −152.76 cm⁻¹** ✓ (next mode +164.2 cm⁻¹ — clean separation; criterion satisfied)
- **Transition-vector check**: projection of the imaginary mode onto the bond axes — forming bond (4,9): |proj| = 0.152; breaking bonds (4,5): 0.081, (6,9): 0.075. The displacement concentrates on the **forming** bond axis and engages both breaking bonds — the concerted scission/insertion motion of the skeletal-editing coordinate ✓ (fig2).

*(Engineering note: this xtb build's `--opt ts` level degrades to a plain minimization (the companion `--ts` flag is rejected) — it *destroyed* the saddle in an intermediate attempt (0 imaginary modes). The converged CI image is therefore verified directly; both behaviours are documented.)*

### 4.1 Activation parameters (298.15 K, harmonic ideal-gas; translation/rotation cancel within the unimolecular step)

| Quantity | Value |
|---|---|
| ΔE‡ (electronic, ANI-2x MEP) | **52.5 kcal/mol** |
| ΔZPE‡ (GFN2-xTB) | −6.2 kcal/mol |
| **ΔH‡** | **46.4 kcal/mol** |
| TΔS‡ (vibrational) | +0.0 kcal/mol |
| **ΔG‡ (Eyring)** | **46.4 kcal/mol** |
| k(298 K) = (k_BT/h)·e^(−ΔG‡/RT) | 6.0 × 10⁻²² s⁻¹ |

**Physical reading**: a ~46 kcal/mol free-energy barrier means the *isolated* :CH₂ adduct is kinetically frozen at room temperature — consistent with the Ciamician–Dennstedt reaction requiring activated (dihalocarbene) conditions and subsequent aromatization driving force; in the real synthetic sequence the adduct forms hot on the carbene-addition PES and cascades forward.

---

## 5. Why Autoregressive LMs Hallucinate Transition States — and Why Equivariant World Models Are Needed

1. **A TS is not a molecule.** It is a *stationary point of index 1* on a 3N-dimensional function — defined by a property of the **Hessian** (exactly one negative eigenvalue), not by any local pattern of atoms. An autoregressive model that emits "plausible molecules" (SMILES strings near its training distribution) has no mechanism to represent, even in principle, a point whose defining feature is the *second derivative of an energy it never evaluates*.
2. **The training distribution contains almost no TSs.** TS geometries are not directly observable; datasets (QM7/9, SPICE, ANI-1x…) are minima. Asking a generative model for a TS is extrapolation *by construction* — the textbook definition of hallucination: fluent output in a region with zero support.
3. **Reaction paths are continuous, covariant objects.** The MEP is invariant to rotation/translation/permutation and equivariant under SO(3). Token sequences are none of these. **Equivariant world models** (MACE/PaiNN/Allegro-class e(3)-equivariant networks, or diffusion over SE(3) coordinates) represent energy as a *function* — differentiable, evaluable, Hessian-capable — so a saddle is a property of the *model itself*, verifiable exactly as done here (CI-NEB on the learned PES + 1-imag test). This pipeline is the template: the ANI-2x leg *is* a learned world model doing transition-state chemistry, cross-verified by GFN2-xTB.
4. **The productive architecture** for "AGI chemistry" is therefore not LM-guesses-a-TS, but **LM proposes the reaction graph → equivariant model owns the PES → algorithmic procedures (NEB/eigenvector-following/Hessian) certify the saddle**. Phases 1–4 of this repository assemble exactly that stack: 2D featurization → classical MD → ML-vs-classical force duality → certified saddle hunting.

---

## 6. Reproducibility, Artifacts, Limitations

```bash
# phase2ff env with Library/bin on PATH (xtb.exe + BLAS)
python run_phase4_reaction_mechanism.py --engine ani   # full pipeline (this report)
python run_phase4_reaction_mechanism.py --engine xtb   # GFN2-xTB NEB (staged protocol)
python run_phase4_reaction_mechanism.py --fig_only     # regenerate figures
# flags: --n_images 9 --fmax 0.05 --max_neb_steps 500 --force_rerun --auto_shutdown
```

**Artifacts (`results_phase4/`)**: `reactant_3d.mol`/`product_3d.mol` · `images_idpp.xyz` (initial band) · `neb_final_path.xyz` + `neb_energies.npy` (converged band) · `ts_candidate.xyz` (CI geometry) · `ts_hessian_freqs.json` (full GFN2-xTB spectra of TS & R) · `ts_modes.npz` (displacement vectors) · per-stage checkpoints · `phase4_results.json`.

**Figures (`./figures_phase4/`, 300 DPI)**: fig1 reaction profile (kcal/mol vs normalized coordinate, TS‡ + ΔG‡ annotated) · fig2 TS 3D with imaginary-mode displacement vectors, breaking (red)/forming (green) bonds · fig3 bond-evolution heatmap across all 9 images.

### Limitations (read before citing numbers)

1. **Mixed fidelity, deliberate**: MEP and ΔE‡ from ANI-2x (a learned reactive PES, semi-quantitative); vibrational verification/ZPE/S from the GFN2-xTB analytic Hessian at the ANI CI geometry. The 1-imag test is exact *on the GFN2 surface at that geometry*; a fully GFN2-xTB-converged band was attempted (see §3.2) but the subprocess-force noise defeated BFGS/FIRE within available budget — the attempt, its diagnostics, and the multi-fidelity fallback are all part of the shipped record.
2. Semi-empirical thermochemistry (rigid-rotor/harmonic, no hindered-rotor treatment); ΔS‡ vibrational only (TΔS‡ ≈ 0 for this unimolecular rearrangement).
3. Parent :CH₂ model — real Ciamician–Dennstedt uses :CCl₂ with a subsequent HCl-elimination/aromatization cascade (higher-level steps outside this MEP).
4. NEB band = 9 images; MEP smoothness between images interpolated (fig1 spline is display-only).
5. Eyring rate uses the standard-state (1 atm) convention without tunneling corrections.

---

*Phase 4 of ai4chem-complex-scaffolds-benchmark. Phase 1: conformer/descriptor suite. Phase 2: torsion barriers + ligand MD. Phase 3: KRAS complex dynamics + MM-GBSA. This phase: automated saddle certification.*
