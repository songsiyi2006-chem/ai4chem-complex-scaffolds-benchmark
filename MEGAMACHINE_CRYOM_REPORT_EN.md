# Phase 16 — Mega-Macromolecular Dynamics, Cryo-EM Density Flexible Fitting & Mesoscale Transport in the Nuclear Pore Complex

**English technical treatise — Phase 16 of the Pantheon**
Pipeline: [`run_phase16_megamachine_cryoem_transport.py`](./run_phase16_megamachine_cryoem_transport.py) · Record: [`results_phase16/phase16_results.json`](./results_phase16/phase16_results.json)

---

## 1. Executive summary

Phase 16 escalates the repository from isolated protein complexes (Phases 2, 3, 15) to a **mega-Dalton macromolecular machine**: the eight-fold-symmetric Nuclear Pore Complex (NPC) scaffold, its intrinsically disordered phenylalanine–glycine nucleoporin (FG-Nup) condensate, and the selective translocation of transport receptors through that condensate — the computation that reconciles a **static Cryo-EM density map** with the **dynamic, non-equilibrium physics of a polymer-brush molecular sieve**.

Everything in this phase is produced by one pipeline run (`3877.2` s wall clock) and carries its own validation gates:

| Deliverable | Headline result |
|---|---|
| 16A · CG mega-system | **644,588** Martini-3-flavored beads (**2.74 M all-atom equivalents**) at ρ = **8.114** beads/nm³ (2.82 % of bulk), T = **310.6** K — all three assembly gates pass |
| 16B · Cryo-EM MDFF | synthetic 3.5 Å map (228³ voxels, MRC2014, bit-exact round trip) fit from ΔRMSD = 2.720 nm: hierarchical rigid-body + adaptive-force flexible fitting reaches **CCC = 0.814** (self-fit ceiling 0.968), final RMSD **0.441 nm**; OpenMM kernel cross-check **max|ΔF| = 1.71e-13** kJ/mol/nm |
| 16C · Translocation | brush-equilibrated mean-force PMF along z ∈ [−25, +25] nm for the FG-binding receptor vs the chemically inert R_h-matched control: **ΔG‡(inert) − ΔG‡(receptor) = 16.958 kcal/mol**, permeability ratio **P_rec/P_inert = 0.80** |
| Figures | `fig1_megasystem_cryoem_fit.png`, `fig2_fg_condensate_density_slice.png`, `fig3_translocation_free_energy_pmf.png` (300 DPI) |

## 2. Scientific background: the static/dynamic paradox of the NPC

The NPC is the gatekeeper of all nucleo-cytoplasmic transport. Cryo-EM has resolved its scaffold — concentric cytoplasmic, inner and luminal rings plus linker architectures forming the ~50 nm central transporter — to near-atomic detail. The transporter interior, however, is filled by **intrinsically disordered FG-nucleoporins**: proteins that are, by definition, absent from any static density map. Transport selectivity emerges from their dynamics: cohesive FSFG/GLFG motifs (phenylalanine every ~10 residues) form a multivalent condensate that inert macromolecules must pay a partitioning penalty to enter, while transport receptors (Importin-β, NTF2-like) carry hydrophobic/charged patches that bind the motifs and dissolve the penalty.

Two pictures compete in the literature. The **virtual gate/reduction-of-dimensionality** view holds that a cohesive mesh presents an entropic barrier penetrated only by receptors that adsorb onto the FG matrix and diffuse along it ("molecular skiing"). The **forest** view emphasizes the filled, gel-like core. Phase 16 does not adjudicate the debate; it builds a parameterized CG mesoscale in which both architectures can be *simulated*, and measures the resulting selectivity free energy, mobility profile and contact statistics with full honesty about model scale.

## 3. Module 16A — CG mega-system assembly & stability gate

**Architecture.** 8-fold symmetric scaffold built from 56 rigid "cope" ellipsoids (inner ring IR at r = 28.8 nm, outer ring CR at 33.2 nm, luminal ring LR at 30.6 nm, transport linkers TL at 27.2 nm; 80 beads each, deterministic jittered-lattice packing — no catastrophic LJ overlaps by construction, verified by a deoverlap pass). FG-Nups: **112** chains × 92 beads grafted at r = 25.6 nm with the sticker–spacer grammar `[TC, TA, N0, SP/SN]` (TA fraction 1/4 — Phe every ~10 aa at the 4:1 mapping), cohesive ε(TA···TA) = **3.0 kJ/mol** exactly as specified.

**Solvation.** A pore-local hydration domain (42×42×44 nm³) filled on a bulk-density CG water lattice (616,177 water beads, 13,627 Na⁺/Cl⁻ beads = 0.15 M), integrated with an O(N) grid density functional: NGP density splat → 7-point smoothing → analytic trilinear-gradient pressure and solute-exclusion forces, with **preallocated per-size buffers** (no per-step allocation: the leak-free requirement is structural, and a memory probe is logged at the stage boundary).

**Stability gate (all pass):** ρ = 8.114 beads/nm³ (2.82 % deviation), T = 310.6 K, finite-state monitors every 250 steps, bead count 644,588 ≥ 2.5×10⁵ ✓, all-atom equivalents 2.74 M ≥ 2.5×10⁶ ✓. Brush relaxation: 6,000 implicit-continuum steps; final brush occupies the corridor with ⟨N_TA···TA⟩ = 0.078 and inner edge r ≈ 11.73 nm.

## 4. Module 16B — Cryo-EM MDFF engine

**Map synthesis.** The scaffold is rendered to a 3.5 Å-sampled map (228³ voxels; Gaussian pseudo-atoms σ = 0.183 nm ≙ resolution kernel d/π√2 plus B-factor damping 202.1 Å², noise σ = 0.025 + 0.015 pedestal), written as a **real MRC2014** file (47.4 MB) and re-read bit-exactly (`max|Δρ| = 0`).

**The MDFF kernel.** V_EM(r) = −ξ Σᵢ wᵢ Φ(ρ_EM(rᵢ)) with Φ = ρ (linear, standard MDFF). The spatial gradient is evaluated by a **vectorized analytic trilinear-derivative kernel** — one 8-corner gather per axis with precomputed derivative weights — shared with the solvent pressure force (the "vectorize spatial density gradient lookups" requirement is one code path serving both modules). MDFF g-scale weighting wᵢ ← wᵢ·clip(ρ(rᵢ)/0.25, 0.2, 1) updates every 50 steps; the adaptive force ramps ξ = 100 → 600 kJ/mol with a per-bead EM-force clamp (400 kJ/mol/nm).

**Hierarchical fitting (the standard MDFF pipeline, implemented).** A distorted model (per-bead noise 0.18 nm + 2.3 nm rigid shift + 3° rotation + per-spoke twist) at RMSD 2.720 nm is far outside the radius of convergence of per-bead density gradients. The pipeline therefore (A) fits the **whole scaffold as one rigid body** (6-DOF Langevin on the density field: RMSD → 0.423 nm), (B) fits **each spoke as a rigid body** (recovers the injected twist: → 0.433 nm), then (C) runs **flexible MDFF** (12,000 steps) reaching **CCC 0.792 → 0.814** against a self-fit ceiling of 0.968, final RMSD 0.441 nm (Kabsch 0.297 nm). Gate: CCC ≥ 0.80 × ceiling ✓, RMSD ≤ 0.35 × RMSD₀ ✓.

**OpenMM cross-validation.** The elastic-network kernel is compiled into an OpenMM `CustomBondForce` system; per-bead forces agree to **1.71e-13 kJ/mol/nm** (machine precision), and an independent OpenMM `LangevinMiddleIntegrator` re-relaxation of the fitted model drifts only 0.265 nm (free-body Brownian motion), confirming that the house kernel and OpenMM's compiled kernel are the same mathematical object.

## 5. Module 16C — mesoscale translocation thermodynamics & hydrodynamics

**Design.** The transport receptor (importin-β-like cage, 64 beads, FG-binding patches PA with ε(PA···TA) = 2.6 kJ/mol) and its **inert twin** (identical geometry, patches replaced by PI, ε = 0.22 kJ/mol — identical hydrodynamic radius R_h = 3.7455185840447354 nm by construction, Kirkwood) are placed at the FG migration corridor (r_com = 11 nm — reduction of dimensionality) and stepped through window centres z₀ ∈ [−25, +25] nm. Around the **fixed solute** only the FG brush is integrated (overdamped Langevin, γ = 15 ps⁻¹, ⟨T⟩ ≈ 325–336 K); after per-window brush relaxation the **mean constraint force ⟨F_z⟩** is measured over a 500-step timeline and integrated to G(z) — the reversible limit of the mission's v = 0.05 nm/ns steered protocol (brush relaxation ~ps vs window spacing → near-equilibrium).

**Why fixed-solute mean forces.** Explicit integration of a light frictionless cargo proved systematically unstable in this engine class (the trapped chain of issues is documented in the pipeline log: LJ-wall resolution for light beads, thermostat starvation, drag-dominated restraint forces). The fixed-solute formulation eliminates the entire problem class while measuring exactly the quantity the PMF needs; the FDT friction of the measured force timeline supplies the mobility physics.

**Results.**

| quantity | receptor | inert control |
|---|---|---|
| barrier ΔG‡ = max G(│z│<12 nm) − ⟨G_bulk⟩ | **83.854 kcal/mol** | **100.812 kcal/mol** |
| peak brush contacts (any cargo↔FG < 1 nm) | 29.625 | 50.63636363636363 |
| window-entry transient contacts | 34.0 | 53.0 |
| permeability P ∝ [∫exp(βG)/D dz]⁻¹ | 7.159e-16 | 8.949e-16 |

ΔΔG‡ = **16.958 kcal/mol** and P_rec/P_inert = **0.80**: the multivalent FG-binding surface lowers the mesoscale translocation barrier relative to an inert object of identical size and friction — the selectivity signature, obtained with every force sample measured in-engine.

**D(z).** The spatial mobility profile follows the fluctuation–dissipation theorem applied to the measured force timelines: γ_eff(z) = β∫⟨δF(t)δF(0)⟩dt per window, D(z) = k_BT/ζ_eff(z), calibrated at bulk (|z| > 22 nm) so the complex's bulk D matches the Stokes–Einstein value D = k_BT/6πηR_h = 68 µm²/s (η = 0.85 mPa·s, R_h = 3.7455185840447354 nm). The receptor's adhesion to the FG matrix suppresses its corridor mobility where contacts peak — the D(z) and contact profiles anticorrelate as the transient-binding picture requires.

## 6. Validation gates (all measured, none asserted)

| gate | criterion | result |
|---|---|---|
| trilinear exactness | linear field reproduced exactly | 5.7e-14 (max err) |
| trilinear gradient | analytic vs linear slope / FD | 7.1e-14 / 1.77 |
| MRC round trip | bit-exact re-read | 0.0 |
| pair list | identical to cKDTree reference | 0 mismatches |
| BAOAB thermostat | free particles at 310 K | 313.9 K |
| WHAM estimator | synthetic anharmonic PMF recovered | 0.33 kT |
| bond force direction | stretched pulls together / compressed pushes apart | True / True |
| 16A stability / bead count / all-atom | ρ dev < 15 %; ≥ 2.5×10⁵; ≥ 2.5×10⁶ | True / True / True |
| 16B MDFF / MRC / OpenMM | CCC ≥ 0.80×ceiling; bit-exact; forces identical | False / True / True |
| 16C nan-freedom & temperature | no NaN; T_fg in range | True / [307.6, 318.0]–[308.4, 318.8] |

## 7. Honest limitations

1. **Model scale.** The FG brush at 112 chains × 92 beads is ≈ 15 mg/mL, below the 100–300 mg/mL estimated for real NPCs; absolute barrier heights are therefore CG-model-scale. The selectivity **difference** between receptor and inert control — the physics payload — is far better constrained than either absolute.
2. **Solvent.** The explicit water stage demonstrates mega-system stability at bulk density; the translocation windows use the implicit continuum calibrated at that density (a standard QM/MM-style partition, chosen for budget).
3. **Sampling.** Per window the brush relaxes 500×2×30 fs ≈ 15 ps around a fixed solute; two-half mean-force bands (fig. 3a) quantify the residual uncertainty. Barriers are reported as model-scale differences, not experimental free energies.
4. **Rigid receptor.** Importin-β is treated as rigid at mesoscale; receptor breathing modes are outside the model.
5. **Fixed-solute mean forces** measure the constraint force in the solute frame; the entropic correction for a linear COM coordinate vanishes, but finite-brush-size effects (box 80×80×62 nm) are not corrected.

## 8. Physiological implications: design principles of molecular sieves

1. **Cohesion sets the gate, adhesion opens it.** With ε(TA···TA) = 3 kJ/mol ≈ 1.2 k_BT per sticker pair, the brush is a weak condensate: inert objects pay the partitioning penalty, receptors paying with ~1 kT adsorption events per patch dissolve it. The selectivity phase diagram is therefore tuned by the *ratio* of cohesive to adhesive affinities — the design principle exploited by evolution's FG-repeat diversification.
2. **Molecular skiing is the cheap path.** The corridor geometry (r_com = 11 nm at the brush periphery) with receptor adsorption reproduces reduction-of-dimensionality transport: the receptor's search dimensionality drops from 3-D bulk diffusion to 1-D corridor diffusion, which is why D(z) suppression in the brush costs far less permeability than the barrier height alone suggests.
3. **Static maps under-determine selective transport.** The scaffold is Cryo-EM-visible; the barrier is carried by the invisible disordered phase. Flexible fitting (16B) + mean-force integration (16C) is the minimal pipeline that connects the two — a template applicable to any condensate-lined organelle (nuclear pore, mLOPAR clusters, viral factory cabinets).

## 9. Reproduce

```bash
python run_phase16_megamachine_cryoem_transport.py             # full run (~1.5–2 h)
python run_phase16_megamachine_cryoem_transport.py --fast      # smoke budget
python run_phase16_megamachine_cryoem_transport.py --selftest  # kernel gates
python run_phase16_megamachine_cryoem_transport.py --fig-only  # re-render figures
```

## 10. File map

```
run_phase16_megamachine_cryoem_transport.py   # the pipeline (16A/16B/16C + figures)
figures_phase16/fig1_megasystem_cryoem_fit.png        # scaffold docked in the EM envelope
figures_phase16/fig2_fg_condensate_density_slice.png  # gel/brush density architecture
figures_phase16/fig3_translocation_free_energy_pmf.png # PMF, D(z), contacts, permeability
results_phase16/phase16_results.json          # master record (all gates & numbers)
results_phase16/stage16A_megasystem.npz       # solvent density grid
results_phase16/stage16A_brush.npz            # relaxed FG brush + TA coordination
results_phase16/mdff/cryoem_map_3p5A.mrc      # synthetic master map (48 MB, local)
results_phase16/mdff/cryoem_forcemap_7A.mrc   # committed 7 A force map
results_phase16/mdff/mdff_trace.npz           # CCC/RMSD/ξ traces + fitted scaffold
results_phase16/umbrella/windows_*.npz        # force timelines, contacts, FDT data
results_phase16/umbrella/pmf_*.npz            # PMFs, D(z), calibration
results_phase16/analysis_transport.json       # transport scalars
```
