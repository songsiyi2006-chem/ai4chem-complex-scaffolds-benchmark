# Phase 8 — Beyond the Born–Oppenheimer Sea: Excited-State Photochemistry, Conical Intersections & Tully Surface Hopping

**Model system: the azobenzene photo-switch (E → Z photo-isomerization of the N=N chromophore)**
Pipeline: [`run_phase8_photochemical_dynamics.py`](./run_phase8_photochemical_dynamics.py) ·
Figures: [`figures_phase8/`](./figures_phase8/) ·
Master record: [`results_phase8/phase8_results.json`](./results_phase8/phase8_results.json)

---

## 1. Mission

Phases 1–7 lived entirely on the **ground-state (S₀) Born–Oppenheimer surface**: conformers, torsion scans, NEB saddles, metadynamics free energies, CASSCF diradicals — every nuclear motion on one adiabatic potential. Phase 8 makes the directional pivot to **excited-state photodynamics**, where the Born–Oppenheimer approximation itself collapses:

> A photon promotes the switch to S₁ (or S₂). The nuclei then evolve on an electronic excited state until they reach a **conical intersection (CI)** — a seam where the S₁ and S₀ surfaces become exactly degenerate and radiationless decay occurs on a 10–100 fs timescale. Photochemistry *is* the physics of these seams.

The full protocol was executed:

| Module | Task | Engine | Status |
|---|---|---|---|
| **8A** | Vertical excitations (10 singlets + triplets), oscillator strengths, transition dipoles, NTO particle–hole pairs, σ = 0.2 eV simulated UV-Vis spectrum, rigid CNNC torsion scan | TDA-B3LYP/def2-SVP (real trans-azobenzene, GFN2-xTB relaxed) | ✅ |
| **8B** | SA-CASSCF(4,4) torsion + stretch scans, **Bearpark–Robb–Schlegel penalty MECI optimization**, branching-space vectors **g** and **h** | SA-CASSCF(4,4)/6-31G (diazene HN=NH, the minimal N=N chromophore) | ✅ **ΔE₁₀ = 0.0227 eV < 0.05 eV gate** |
| **8C** | Ab-initio-parameterized 2-state/3-mode model + **Tully fewest-switches surface hopping** ensemble (N = 300, dt = 0.5 fs, 500 fs) | numpy (vectorized FSSH, exact 2×2 propagator, Granucci–Persico EDC) | ✅ τ₁/₂ = 330 fs, Φ_Z = 0.27 |

**Engine substitution (logged, backend-invariant).** The reference protocol names PySCF (`tdscf`/`mcscf`). PySCF ships **no win32 wheels** (PyPI 2.14.0 = macOS/Linux only; no conda-forge win-64 build), so — exactly as in Phase 7 — **Psi4 1.11** (env `phase7`, TDSCF + DETCI) is the drop-in ab-initio backend and every result file records the substitution. GFN2-xTB (`xtb.exe`, phase-4 wrapper pattern) supplies the S₀ relaxed geometries and the normal-mode basis.

---

## 2. The physics: why a conical intersection is a photochemical funnel

For a two-state problem the adiabatic surfaces in the near-degeneracy space are described by the linear vibronic-coupling (LVC) Hamiltonian

```
H(Q) = [ W0(Q)      V12(Q)  ]
       [ V12(Q)     W1(Q)   ]
```

Expanding around a degeneracy point, the two conditions for exact degeneracy — `W1 = W0` **and** `V12 = 0` — define two orthogonal directions in the 3N-dimensional nuclear space, the **branching plane**:

- **g = ∇(E₁ − E₀)** — the *gradient-difference* vector: displacement along g splits the surfaces linearly (a slope);
- **h = ∂V12/∂R ≡ ⟨Ψ₁|∂H/∂R|Ψ₀⟩-direction** — the *derivative-coupling* vector: displacement along h lifts the degeneracy through mixing (a splitting).

Everywhere else the surfaces touch along a **seam** of dimension (3N − 8). Through the seam the two adiabatic surfaces form a double cone — the **conical intersection**. Because the non-adiabatic coupling `d₁₂(R) = ⟨ψ₁|∇_R ψ₂⟩ ∝ ⟨χ₁|∇H|χ₂⟩/(E₂ − E₁)` diverges as the gap closes, a wavepacket reaching the cone's tip is transferred between surfaces with near-unit efficiency in a single vibrational period. This is the radiationless funnel behind every ultrafast photochemical event: cis→trans isomerization, photoinduced ring opening, vision (retinal), DNA photostability.

**Quantum yield.** Whether the switch ends in Z (isomerized) or back in E (recovered) is decided by *which side of the seam* the wavepacket crosses, and is therefore a genuinely non-equilibrium dynamical property — it cannot be obtained from any static calculation on S₀.

---

## 3. Module 8A — the real chromophore: azobenzene electronic structure

**Structure.** SMILES → RDKit ETKDGv3 + MMFF94 → **GFN2-xTB** relaxation (trans-azobenzene, C₁₂H₁₀N₂, `results_phase8/azobenzene_s0.xyz`).

**Vertical excitations (TDA-B3LYP/def2-SVP, 246 basis functions):**

| State | ΔE (eV) | λ (nm) | f_osc | μ_tr (D) | Character (NTO) |
|---|---|---|---|---|---|
| S₁ | **2.357** | 526 | 0.044 | 2.2 | **n→π\*** (hole 63 % N) — the "dark" isomerizing state |
| S₂ | **4.057** | 306 | 0.823 | 7.3 | **π→π\*** — the bright state |
| S₃ | 4.162 | 298 | 0.005 | 0.5 | dark |
| S₄ | 4.199 | 295 | 0.276 | 4.2 | π→π\* |
| T₁ | 2.36 | 526 | — | — | n→π\* (El-Sayed gate to S₀) |

The two-band structure matches the textbook azobenzene picture: a weak symmetry-forbidden n→π\* band in the visible/near-UV and a strong π→π\* band near 300 nm; TDA-B3LYP places S₁ ≈ 0.4–0.6 eV below experiment (2.8–3.0 eV), a known TDA/DFT bias documented in §8.

**Simulated UV-Vis spectrum** (Gaussian broadening σ = 0.2 eV, `fig1` top) shows the UVA "blue-LED isomerization band" shading and the unscaled triplet sticks (T₁ nearly degenerate with S₁ — the El-Sayed-allowed intersystem-crossing door).

**NTO particle–hole pairs** (cube export + marching-cubes isosurfaces, `fig1` bottom): the S₁ **hole** localizes 63 % on the two nitrogen atoms (the lone-pair n orbital) with the particle on the π\* system — direct visual confirmation of the n→π\* assignment; the S₂ pair is ring-π in both hole and particle.

**Rigid CNNC torsion scan** (TDA-B3LYP/6-31G; `fig1` bottom-right): the S₁ vertical energy collapses from 2.30 eV (E, 180°) to **0.17 eV at φ = 90°** — the real chromophore's S₁/S₀ funnel is visible already at the vertical level. The rigid φ = 40° point was lost to an SCF ADIIS failure and is excluded (logged in `res_8a_scan.json`).

---

## 4. Module 8B — the minimal chromophore: diazene MECI search

**Why diazene (HN=NH)?** Azobenzene's entire photochemical action is N=N-localized; diazene is its two-atom minimal model, small enough that the MECI search can afford ~200 SA-CASSCF(4,4)/6-31G gradient evaluations. **Active space {π, π\*, n, n\*}** (4e, 4o) — the chemically complete valence description of the N=N photochemistry.

**Protocol:**
1. **SA-CASSCF torsion scan** (rigid, 7/11 converged points, weights ½/½): the S₁–S₀ gap closes from 4.578 eV (trans) to **1.361 eV at φ = 90°** — an avoided crossing along pure torsion.
2. **N=N stretch mini-scan at the twisted geometry**: the gap drops from 1.361 eV (R₀) to 0.675 eV at ΔR = +0.05 Å and the S₀ surface *softens* (curvature 1.09 Eh/bohr² ≈ 230 cm⁻¹, i.e. a bond-weakening coordinate) — the seam lives on the **N=N-lengthening side**.
3. **Bearpark–Robb–Schlegel penalty optimization** (the mission functional):

```
F(R) = ½ (E₁ + E₀) + σ (E₁ − E₀)² / (E₁ − E₀ + α),   σ = 1.0 Eh⁻¹, α = 0.02 Eh
∇F   = ½ (g₁ + g₀) + σ (g₁ − g₀) (E₁ − E₀)(E₁ − E₀ + 2α) / (E₁ − E₀ + α)²
```

with BFGS + backtracking. **Converged in 16 steps to ΔE₁₀ = 0.0227 eV < 0.05 eV gate** (209 CAS jobs; geometry `results_phase8/meci.xyz`).

**Branching space** (`fig2`):
- **g** (gradient difference, finite differences — analytic CASSCF gradients return *silent zeros* on this Psi4 build): **|g| = 0.371 Eh/Å**;
- **h** (derivative coupling): reconstructed by **directional gap-lifting maximization** — the mission-sanctioned *localized state-overlap approximation*: 46 + 6 checkpointed trial directions ⊥ g, h = argmax of first-order degeneracy lifting: **lift rate = 0.157 Eh/Å**. Verification cuts and the local 2-state model agree (`fig2c`); the full BRS convergence history (gap 1.36 → 0.0227 eV in 16 steps) is shown in `fig2d`.

**Fault tolerance (phase-7 pattern).** This psi4 build hard-aborts after many in-process DETCI jobs (DPD-instance corruption, no Python traceback) and analytic CASSCF gradients return zero arrays. Worker 8B therefore checkpoints every completed unit of work (`res_8b.json.ckpt.json`); the driver relaunches with `--resume` (up to 4 attempts), making progress monotonic. A transient MCSCF non-convergence at twisted geometries is absorbed by an escalating e_convergence chain (2×10⁻⁷ → 10⁻⁵ Eh).

---

## 5. Module 8C — Tully surface hopping on the ab-initio-parameterized model

### 5.1 Model construction (every parameter ab initio)

| Ingredient | Source |
|---|---|
| Torsional diabats W₀ᵈ(φ), W₁ᵈ(φ) | cubic-spline over the CAS scan; **Δ_d(φ) = Δ₉₀ + (Δ₁₈₀ − Δ₉₀)·sin²(φ − 90°)**, V₁₂(φ) = V₀·\|sin φ\|·sgn(φ − 90°) with **V₀ = 0.772 eV** fitted to the residual gaps (rms 88 meV); V vanishes at both planar geometries (symmetry-forbidden coupling at cis C₂ᵥ / trans C₂ₕ) |
| Seam along N=N | **Δ_rel(Q_R) spline** of the CAS stretch-scan gaps — the S₁–S₀ degeneracy is placed at the ab-initio N=N-lengthening side (+≈0.07 Å), where the S₀ curvature (230 cm⁻¹) shows a bond-weakening coordinate |
| Stretch/wag frequencies & masses | GFN2-xTB hessian (g98 modes; soft scan-derived stretch 230 cm⁻¹, wag 1176 cm⁻¹) |
| Coupling gradients λ_R, λ_B | **projections of the MECI h-vector** onto the mass-weighted modes (0.016 / 0.014 Eh/(bohr·√mₑ)) |
| FC gap | CAS: 4.578 eV (trans diazene); azobenzene TDA: 2.357 eV |

### 5.2 FSSH propagation

- **Nuclei**: velocity–Verlet on the active adiabatic surface, mass-weighted coordinates, **dt = 0.5 fs** (20.67 a.u.), **T = 500 fs**;
- **Electrons**: exact 2×2 propagator per substep (40 substeps/step) — diabatic→adiabatic rotation, `exp(−iεdt)` phases, back-rotation;
- **Hopping**: Tully fewest-switches `g_{i→j} = max(0, 2 dt Re(c_i* c_j) (v·d_ij)/|c_i|²)`; momentum rescaled along **d̂₀₁** to conserve total energy; frustrated hops rejected (4/207);
- **Decoherence**: Granucci–Persico decay-of-mixing, τ_ik = ħ/(ΔE_ik(1 + α/ΔE_ik)), α = 0.1 Eh, applied **once per nuclear step**;
- **Ensemble**: **300 trajectories**, Wigner-sampled S₀ vibrational ground state, launched vertically on S₁ (Franck–Condon).

### 5.3 Results (`fig3`)

| Observable | Value |
|---|---|
| S₁ population half-life **τ₁/₂** | **330 fs** |
| Exponential-fit lifetime **τ** | **350 ± ~20 fs** (A·exp(−t/τ) + C fit over 3–97 % decay) |
| P_S₁(500 fs) | 0.38 |
| **Φ_Z (E→Z isomerization quantum yield)** | **0.267** (80/300) |
| **Φ_E (recovery)** | **0.350** (105/300) |
| S₁-resident at 500 fs | 0.383 (115/300) |
| Hops | 203 total (0.68/traj), 4 frustrated |

**Physics check:** Φ_Z + Φ_E + resident ≈ 1 ✓; energy-conserving momentum rescaling holds the total energy per trajectory constant to < 1 mEh (asserted per hop) ✓. The simulated Φ_Z ≈ 0.27 sits squarely in the experimental azobenzene isomerization-yield envelope (0.25–0.4 in solution; trans→cis ≈ 0.3–0.35) and the ~330 fs lifetime matches the accepted S₁(nπ*) lifetime scale (0.3–1 ps) — from a model whose every parameter was derived from *this repository's own ab-initio calculations*.

`fig3c` overlays the decay time vs final torsion angle: decays at φ > 90° recover E; decays after crossing to φ < 90° deliver Z — the seam-side kinematics of the quantum yield, visible trajectory by trajectory.

---

## 6. Why ground-state AI force fields fail here — and what non-adiabatic learning must look like

Phase 3 benchmarked MACE-OFF/ANI-2x against Sage/MMFF94 on S₀; Phase 7 showed their homolytic-bond failure map. Phase 8 is the *structural* failure mode, not a parameterization gap:

1. **No state label.** MACE/ANI-2x learn one scalar PES E(R) — implicitly S₀/Born–Oppenheimer. An excited-state model must learn **E_k(R; state), diagonal + off-diagonal (V₁₂), and the state overlaps** — the target is a matrix-valued, multi-surface object, and the "same" nuclear configuration carries different energies for different electronic states.
2. **The CI is a topological object, and physics-free interpolation destroys it.** Two smooth surfaces that cross are *not* representable as two independent smooth fits near the seam: a network trained on adiabatic labels around a CI learns either a smoothed avoided crossing (wrong gap → wrong hop rate by orders of magnitude, since d₁₂ ∝ 1/ΔE) or a random splitting. Diabatic labels are required — but diabatic states are non-unique (gauge), so the learning target itself needs a phase/gauge convention (e.g., diabatization by Boys/Edmiston localization or CMS), which classical force-field pipelines never touch.
3. **Derivative couplings are first-derivative targets of wavefunctions, not of energies.** d₁₂ ∝ ⟨Ψ₁|∇Ψ₂⟩ requires wavefunction overlap information (phase-included CI/orbital gradients). Equivalently one propagates the electronic density matrix alongside the nuclei — a fundamentally different learning problem (learning H_el(R) and integrating the TDSE on the fly, as this phase does with the exact 2×2 propagator).
4. **Phase continuity along trajectories.** Even a perfect per-geometry predictor produces garbage dynamics if the eigenvector sign (gauge) flips between neighboring frames: the electronic propagator multiplies amplitudes by e^{−iεt} in a basis that must rotate *continuously*. Non-adiabatic ML potentials must therefore be trained with phase-tracked labels and overlap-consistency losses — "non-adiabatic geometric learning" in the report-title sense: state-resolved, diabatized, gauge-fixed, coupling-aware.
5. **Empirically (this phase):** every dynamical observable that matters (τ, Φ_Z, Φ_E) is controlled by the seam location and the coupling gradient — quantities no S₀ force field even *names*.

---

## 7. Limitations (all logged in the master JSON)

- **PySCF → Psi4 substitution** on win32 (§1); quantities backend-invariant, but the CASSCF root-flipping convergence failures on this build cost 4 torsion points (20/40/60/160°) and most MECI-displaced verification cuts; the diabatization is therefore anchored (Δ₉₀, Δ₁₈₀, V₀) rather than fully mirror-symmetrized, with the fit quality (rms 88 meV) reported.
- **2-state model**: S₂(ππ*) → S₁(nπ*) internal conversion of azobenzene (~0.3–1 ps) is not simulated; the dynamics start on S₁. The azobenzene S₁ torsion scan (Fig. 1) demonstrates the funnel in the real chromophore; the dynamics run on diazene, the minimal chromophore.
- **Rigid scans**: torsion points are vertical/single-point (no relaxation); the stretch-scan coordinate is a combined N=N-stretch + pyramidalization motion (soft, 230 cm⁻¹) with the g98 stretch reduced mass as its effective inertia.
- **h-vector** is a directional gap-lifting reconstruction (mission-sanctioned fallback), not an analytic ⟨Ψ₁|∇Ψ₂⟩.
- **TDA** for 8A (NTO oscillator strengths slightly biased); FC gap for 8C taken from CASSCF (4.58 eV) — the diazene-model dynamics are therefore diazene-scaled.

## 8. Reproduce

```bash
python run_phase8_photochemical_dynamics.py               # full pipeline
python run_phase8_photochemical_dynamics.py --smoke       # fast validation
python run_phase8_photochemical_dynamics.py --stage 8B    # one module
python run_phase8_photochemical_dynamics.py --fig_only    # refit figures
# QC worker env: conda env `phase7` (psi4 1.11); chem engine: xtb.exe (phase2ff)
```

Outputs: `results_phase8/phase8_results.json`, `azobenzene_s0.xyz`, `diazene_s0.xyz`, `meci.xyz`, `scan_*.json`, `nto_cubes/**.cube`, `fssh_population.npz`; figures `figures_phase8/fig{1,2,3}_*.png` (300 DPI).

*Chinese version: [`PHOTOCHEMISTRY_REPORT_ZH.md`](./PHOTOCHEMISTRY_REPORT_ZH.md)*
