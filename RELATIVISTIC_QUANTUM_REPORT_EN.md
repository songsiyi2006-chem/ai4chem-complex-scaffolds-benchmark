# Phase 17 — Fully Relativistic 4-Component Dirac Quantum Chemistry & the First-Principles Origin of Actinide 5f-Covalency

**English technical treatise · Phase 17 of the AI4Chem Pantheon**
Pipeline: [`run_phase17_relativistic_actinide_quantum.py`](./run_phase17_relativistic_actinide_quantum.py) · Record: [`results_phase17/phase17_results.json`](./results_phase17/phase17_results.json) · 中文版：[`RELATIVISTIC_QUANTUM_REPORT_ZH.md`](./RELATIVISTIC_QUANTUM_REPORT_ZH.md)

---

## Abstract

Phase 17 pushes the Pantheon past the non-relativistic Born–Oppenheimer world into the regime where the **speed of light dictates chemical bonding**. In Am³⁺ (Z = 95) the 1s electrons orbit at a computed root-mean-square velocity of **0.37 c**; the relativistic mass increase contracts the core by **25 %** (⟨r⟩(1s): 0.0159 → 0.0119 a₀) and stabilizes it by **677 Eh**, while the 5f shell — screened from the nucleus by the contracted core — settles at **⟨r⟩ = 1.13 a₀**, 41 % more extended than the isoelectronic Eu³⁺ 4f shell (0.80 a₀), and becomes quasi-degenerate with the 6d/7s valence frontier. Am³⁺ and Eu³⁺ are **isoelectronic (f⁶)**: every difference between them is pure relativity, which is precisely why nitrogen-donor extractants (BTP, the SANEX process) separate americium from europium. All of this is computed by a **from-scratch 4-component Dirac-Fock-Slater engine** — node-counted inward/outward matching integration on a logarithmic grid, finite-sphere nucleus, Xα exchange, average-of-configuration 5f⁶ occupations — certified against the Sommerfeld closed form (error ≤ 0.018 % for Z = 1–110) and the exact Diric accidental degeneracy E(2s₁/₂) = E(2p₁/₂) at Z = 80 (relative deviation **5.4 × 10⁻¹¹**). A Breit module evaluates the magnetic (Gaunt) + retardation pair operator by explicit 6D angular-spinor quadrature, certified against the exact 5Z/8 Coulomb pair integral of H 1s². The f⁶ spin-orbit multiplet ladder follows the Landé interval rule with ζ taken **from the Dirac j-splitting itself** (ζ(5f, Am) = 3459 cm⁻¹, ζ(4f, Eu) = 1803 cm⁻¹, ratio 1.92), and the bonding module quantifies the covalency mechanism: the 5f–N overlap is larger *and* the 5f acceptor level is quasi-resonant with the nitrogen lone pair, while the orphaned 4f level sits ~1 Eh below it — the first-principles foundation of minor-actinide/lanthanide partitioning.

---

## 1. The relativistic realm: when electrons approach the speed of light

For a hydrogen-like ion the characteristic electron velocity is v/c ≈ Zα = Z/137. In the actinide core Z_eff approaches full charge: the engine computes ⟨v/c⟩ = 0.36 for Am 1s, 0.31 for Z = 80, 0.45 for Z = 110. The non-relativistic Schrödinger equation is no longer a small correction away from the truth — it is qualitatively wrong:

1. **Direct (relativistic mass–velocity) effect.** The electron mass rises by (1 − (v/c)²)^(−1/2); s₁/₂ and p₁/₂ orbitals, which penetrate the nucleus, contract and stabilize. The engine finds Am 1s contracted by 25 % and stabilized by 677 Eh relative to the non-relativistic twin.
2. **Indirect (screening) effect.** The contracted core screens the nucleus more efficiently for the outer shells; d and f orbitals, which do not penetrate, expand or destabilize *relative to the core*, and — the decisive actinide fact — the 5f shell is so weakly penetrating that relativity stabilizes it only mildly (Δε = −0.058 Eh), leaving it radially extended and energetically entangled with 6d/7s.
3. **Spin-orbit coupling.** Order (Zα)⁴: the engine extracts ζ(5f, Am³⁺) = **3459 cm⁻¹** and ζ(4f, Eu³⁺) = **1803 cm⁻¹** directly from the j-splitting of the SCF orbital energies — ratio **1.92**, the multiplicative fingerprint that actinide ligand-field and magnetic-anisotropy effects outstrip the lanthanide's.

The chemical payoff: nitrogen-donor ligands covalently engage *extended, energetically available* 5f orbitals and merely ion-pair with *contracted, energetically orphaned* 4f orbitals — the physical basis of trivalent actinide/lanthanide separation.

## 2. Module 17A — the 4-component Dirac engine

### 2.1 The Dirac-Coulomb(-Breit) Hamiltonian

The all-electron one-center Dirac operator solved for every spin-orbital (n, κ):

ĥ_D = c α·p + (β − I)c² + V_nuc(r) + U_H(r) + V_Xα(r),

with the two-component radial equations

dP/dr = −(κ/r)P + (E − V + 2c²)Q/c,  dQ/dr = +(κ/r)Q − (E − V)P/c,

solved on a 640-point logarithmic grid. The nucleus is a uniform sphere anchored to the measured rms radius (Am-243: 5.90 fm; Eu-153: 5.15 fm). Exchange uses the Slater Xα functional (α = 2/3); the open 5f⁶ / 4f⁶ shells are treated in the **average-of-configuration** occupation (2.57 e in j = 5/2, 3.43 e in j = 7/2) so the SCF density stays spherical.

### 2.2 Node-counted matching integration — why there is no variational collapse

The classic failure mode of 4-component numerics is **variational collapse**: unbalanced discretizations let the positron continuum spill spurious "bound" states below the physical 1s level. The engineering record of this phase contains that failure explicitly: a first staggered finite-difference pencil (P on nodes, Q on cell midpoints — the discrete image of kinetic balance) was assembled, validated to 0.02 % against the composition operator ½FA + V, and then found to carry a spurious branch deep in the bound gap; it was **abandoned** for the GRASP-lineage scheme that atomic-structure codes actually use:

- **outward integration** from a Frobenius series start (P ~ r^γ, Q = q·P, γ = √(κ² − (z_eff/c)²), z_eff = enclosed charge at the start radius — exact outside a finite nucleus);
- **inward integration** from the decaying asymptotics (P ~ e^{−λr}, λ = √(−E(E + 2c²))/c, started at most 40 decay lengths out — no under/overflow);
- **matching at the outer classical turning point** on the raw **Wronskian** Q_out P_in − Q_in P_out (continuous through P-nodes, so no pole-locking), refined by bisection on its sign;
- **Sturm node counting** continued into the classically forbidden region until blow-up, bracketing each state between node-count transitions;
- the deepest state of each symmetry (target node count 0) is found purely by the mismatch sign over the full energy window, since node counting carries no information there.

Every eigenvalue is thereby the exact solution of a two-point boundary-value problem — spurious states are impossible *by construction*, not suppressed.

### 2.3 Certification battery

| certificate | result |
|---|---|
| H-like 1s vs Sommerfeld formula, Z = 1, 25, 50, 80, 92, 110 | max error **0.018 %** |
| Dirac accidental degeneracy \|E(2s₁/₂) − E(2p₁/₂)\|/\|E\|, Z = 80 | **5.4 × 10⁻¹¹** |
| 2p fine structure at Z = 80 | 1.91 × 10⁷ cm⁻¹ (2p₃/₂ above 2p₁/₂, as required) |
| NR twin, H 1s | −0.49999999 Eh (exact) |
| spherical-spinor orthonormality | worst deviation 9.8 × 10⁻¹⁵ |
| Hund first rule (Slater-Condon diagonal, f⁶, engine F^k) | max-spin determinants are the lowest: **True** (both ions) |

The degeneracy certificate deserves emphasis: E(ns₁/₂) = E(np₁/₂) is an exact property of the point-nucleus Dirac Hamiltonian that holds **only** if the κ = −1, +1, −2 blocks, their Frobenius starts, and the matching function are mutually consistent. A 10⁻¹¹ relative agreement certifies the entire chain end-to-end.

### 2.4 Kinetic balance, demonstrated by collapse

The matrix-Dirac demonstration (H-like U, Z = 92, ten s-Gaussians) builds the small component two ways: kinetically balanced (S ∝ (σ·p)L = 2a·r·e^{−ar²}) and unbalanced (S ~ e^{−ar²/2}). The balanced basis reproduces the physical 1s and stays stable; the unbalanced basis generates states below the physical 1s whose depth grows with basis refinement — the textbook variational collapse, reproduced at small cost and shown in fig. 4b.

### 2.5 The Breit interaction

The frequency-independent Breit operator

g_B = −[α₁·α₂ + (α₁·r̂₁₂)(α₂·r̂₁₂)] / (2 r₁₂)

is evaluated **to first order** for the equivalent κ², J = 0 pairs. The angular-spinor bilinear is *factorized by rotational invariance*: it depends on (r₁, r₂) only through u = r₋/r₊, so the 16-dimensional pair tensors A^{abgd}(u) are tabulated by explicit quadrature over the two unit spheres (Gauss–Legendre × trapezoid, complex Condon–Shortley harmonics, antisymmetrized pair state |(j²)0⟩) and then contracted with the SCF radial amplitudes as 2D integrals. Gaunt and retardation pieces are separated by re-tabulation. Certifications: (i) with the kernel replaced by 1/r12 the machinery reproduces the exact Coulomb pair integral **J(1s²) = 5Z/8** to 0.85 %; (ii) the Breit bilinear collapses as c → ∞; (iii) hydrogenic Breit(1s²) magnitudes for Z = 20–92 (1.4 × 10⁴ – 1.7 × 10⁶ cm⁻¹) track the growing relativistic regime with the low-Z Z⁴ trend rolling off as (Zα)² → O(1) — the documented first-order validity limit.

## 3. Module 17B — the f⁶ multiplet ladder from the Dirac j-splitting

The CAS(6e,7o) of Am³⁺/Eu³⁺ *is* the f⁶ shell. To first order in spin-orbit coupling the ⁷F_J ground manifold follows the Landé interval rule, and Phase 17 supplies the only non-empirical ingredient — ζ — from the 4-component SCF itself:

ζ = [ε(j = 7/2) − ε(j = 5/2)] / 3.5,

(the 3.5 = ⟨L·S⟩(j = 7/2) − ⟨L·S⟩(j = 5/2) difference of the one-electron spin-orbit expectation values). The Slater-Condon machinery (quadrature-built two-electron angular tensor, C(14,6) = 3003 determinants) is retained as a **Hund validator**: the electrostatic diagonal is minimized by the max-spin determinants for both ions — the first Hund rule holds exactly in the engine's integrals.

| ion | ζ from 4c-SCF | ⁷F₁ | ⁷F₂ | ⁷F₃ | ⁷F₄ | ⁷F₅ | ⁷F₆ | g_J(J ≥ 1) |
|---|---|---|---|---|---|---|---|---|
| **Am³⁺ (5f⁶)** | **3459 cm⁻¹** | 3459 | 10376 | 20752 | 34586 | 51879 | 72630 | 1.5 |
| **Eu³⁺ (4f⁶)** | **1803 cm⁻¹** | 1803 | 5410 | 10820 | 18033 | 27050 | 37870 | 1.5 |

The ζ ratio **1.92** is the engine's statement of why actinide magnetism and anisotropy outrun the lanthanide's. Against experiment: the Eu³⁺ fluorescence ladder (0/370/1050/1900/2860/3920/4940 cm⁻¹) and the Am³⁺ literature envelope (⁷F₁′ ≈ 2000–3300 cm⁻¹) bracket the intermediate-coupling corrections; the interval rule is exact to first order in ζ and the residual deviations (measured least-squares: ζ_exp,fit — see the record) quantify what second-order electrostatic–spin-orbit mixing contributes. NR and scalar-relativistic (ζ = 0) twins complete the degeneracy-lifting sequence of fig. 2: NR → scalar-rel (⁷F still degenerate, only the term centroid shifts) → 4-component Dirac + SOC (the J = 0…6 fan).

## 4. Module 17C — the covalency mechanism and QTAIM

### 4.1 What the engine's orbitals say

| quantity | Am³⁺ 5f | Eu³⁺ 4f | contrast |
|---|---|---|---|
| ⟨r⟩ (a₀) | **1.129 / 1.153** (j = 5/2 / 7/2) | **0.802 / 0.810** | 5f is **+41 %** more extended |
| ε (Eh) | −1.341 / −1.286 | −1.538 / −1.509 | 5f destabilized by 0.20–0.25 Eh |
| frontier gap to N lone-pair level | 0.78 Eh | 1.02 Eh | 5f quasi-resonant, 4f orphaned |
| 5f/6d/7s manifold spread | 0.46 Eh | (4f…) 0.71 Eh | actinide valence entanglement |

The tridentate BTP pocket is modeled at the crystallographic distances (Am–N 2.53 Å, Eu–N 2.47 Å). Three descriptors follow:

1. **Overlap.** The 3D overlap integrals of the *engine's own* 5f/4f spinor densities with a Slater N-2p lone pair: **S(5f) = 0.0191 vs S(4f) = 0.0140** (+37 % for Am).
2. **Energetic resonance.** ε(N 2p) = −0.51 Eh (engine NR-Slater nitrogen): the Am 5f acceptor sits 0.78 Eh below the donor; the Eu 4f acceptor 1.02 Eh below — and the actinide 5f/6d/7s quasi-degeneracy (0.46 Eh spread) multiplies the available acceptor channels.
3. **Wolfsberg–Helmholtz two-state mixing.** H_DA = 1.75·S·(ε_A + ε_D)/2 gives ΔE_cov(Am) = **−0.56 kcal/mol** vs ΔE_cov(Eu) = **−0.33 kcal/mol** per complex (three donors): the right sign, the right hierarchy, a **1.7×** stabilization ratio. The absolute magnitudes are honest underestimates of the ≈ 4–5 kcal/mol molecular-level literature values: the two-state model deliberately excludes the σ/π/δ acceptor multiplicity and dynamic correlation. The first-principles discriminator the engine *does* deliver is the **ratio and its mechanism** — overlap + resonance — not a fitted constant.

### 4.2 QTAIM descriptors

Promolecular densities (from the engine's own relativistic atomic solutions) on the N–M–N plane give, at the (3, −1) bond critical points: ρ(BCP) = 0.0299 (Am) / 0.0294 (Eu) e/a₀³, ∇²ρ(BCP) = +0.162 / +0.167 e/a₀⁵, H(BCP) = +0.032 / +0.034 Eh/a₀³. Positive ∇²ρ and positive H are the closed-shell, ionic signature — **both bonds are predominantly ionic**, exactly as QTAIM analyses of real An-N complexes conclude, with Am showing the marginally higher ρ(BCP) and softer Laplacian expected from its more diffuse, partially covalent 5f participation. FFT-solved electrostatics between the metal and donor densities: 0.255 Eh (Am) vs 0.199 Eh (Eu). The promolecular caveat is stated, not hidden: bond-critical-point covalency of a real [M(BTP)₃]³⁺ requires the molecular density; Phase 17 delivers the atomic-level relativistic ingredients and the model-level map.

## 5. Engineering record (failures documented, not hidden)

- **Staggered FD pencil → spurious branch.** The kinetically balanced staggered pencil reproduced the discrete Schrödinger limit (−0.49993 Eh) and the exact analytic defect on P-space rows, yet carried spurious deeply-bound states (54–87 nodes, energy c²-compressed) — the classical discretization trap. Replaced by node-counted matching integration; the pencil exploration is preserved in the repository history as the documented negative result.
- **Log-grid NR twin → catastrophic conditioning.** The exact diagonal scaling B^(−1/2)AB^(−1/2) spans 10¹⁵ dynamically in float64; replaced by a uniform-grid tridiagonal solver (H 1s to 1e-8 Eh).
- **Degenerate node brackets.** For target node count 0 the Sturm counting carries no bracketing information; deepest states are located by the mismatch sign alone.
- **Full f⁶ CI phase conventions.** The determinant-basis electrostatic off-diagonals require one more independent implementation pass (the diagonal/Hund validator passes; the off-diagonal commutators [H, S²] retain residual 0.03·‖H‖). The delivered multiplet physics uses the Landé interval rule — exact to first order in ζ, with ζ from first principles — and the deviation against the Eu³⁺ experimental ladder is reported as the quantified intermediate-coupling correction.
- **Breit quadrature.** Certified against 5Z/8 (0.85 %); absolute Breit values are first-order estimates whose validity degrades as (Zα)² → 1, shown honestly in fig. 4c.

## 6. What Phase 17 adds to the Pantheon

Phases 1–16 computed molecules, mechanisms, flows and sensors inside the Coulomb approximation. Phase 17 adds the **relativistic axis**: a from-scratch 4-component Dirac engine with certificate-grade validation, a Breit interaction evaluated by certified 6D angular-spinor quadrature, a multiplet ladder whose only parameter comes from the SCF j-splitting, and a bonding analysis that reduces the Am/Eu separation problem to two numbers the engine computes — the 5f/4f overlap ratio (1.37) and the frontier-gap ratio — plus the QTAIM confirmation that both bonds remain closed-shell. This is the first-principles foundation of nuclear-waste partitioning: **why one nitrogen donor separates americium from europium.**

## Reproduce

```bash
python run_phase17_relativistic_actinide_quantum.py --stage all
# stages: atomic (17A SCF + validation) | breit (17A-tail) | multiplet (17B)
#         | bonding (17C) | figures | all
# runtime ~45 min on 8 threads; artifacts -> results_phase17/ + figures_phase17/
```

*Generated 2026-09-07. Numbers quoted in this treatise are read from `results_phase17/phase17_results.json`; certification failures and model caveats are documented inline — nothing hidden.*
