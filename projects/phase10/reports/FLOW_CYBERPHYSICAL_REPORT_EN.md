# FLOW CYBER-PHYSICAL REPORT — PHASE 10

## The Multi-Scale Continuous-Flow Digital Twin: Operando PAT & Cyber-Physical Reinforcement Control

*Pipeline: `run_phase10_cyberphysical_flow_twin.py` · Results: `results_phase10/phase10_results.json` · Figures: `figures_phase10/` (300 DPI)*

---

## Abstract

Phases 1–9 built a chemical world model and then taught it to move: conformers, transition states, microkinetics, metadynamics, photodynamics, and finally a robot that executes the model's own predictions. Phase 10 completes the convergence the architecture promised — **the molecule and the plant are now one simulated object**. A coiled tubular microreactor (L = 10.0 m, d_t = 1.0 mm, 7.85 mL, counter-current cooling jacket) synthesizes the Phase-4/5 flagship target — the strain-release asymmetric aziridine ring expansion catalyzed by the designed BINOL phosphoric acid — as a genuine multi-scale entity: the **micro-scale** enters through the Phase-5 microkinetic rate laws (88.87 kJ/mol RDS barrier, ΔΔG‡ = 1.50 kcal/mol stereodifferentiation → ee = 85.3 % at 298 K); the **meso-scale** enters through Taylor–Aris axial dispersion, Hagen–Poiseuille hydraulics and laminar film heat transfer; the **macro-scale** enters through a 1 Hz sensor-and-actuator cyber layer that closes the loop. The coupled 9-field advection–diffusion–reaction PDEs are integrated by an IMEX Method-of-Lines solver (implicit Taylor–Aris diffusion + TVD-RK2 transport/reaction with depletion-aware adaptive substepping). **Module 10B** renders operando multi-modal PAT telemetry — in-line 785 nm Raman with Lorentzian fingerprints, UV-Vis at 254/310 nm, Hagen–Poiseuille ΔP and port thermometry — under four injected industrial anomalies (pump cavitation, precursor impurity spike, progressive channel clogging, coolant vapor lock); the deconvolution agent that converts telemetry into concentrations is audited at **mean relative error 0.5 % (max 3.6 %)** against plant truth. **Module 10C** closes the loop with a torch Soft-Actor-Critic agent (18-dim PAT state, 5 continuous actions: three HPLC micro-dosing pumps, jacket rate, back-pressure setpoint) trained on domain-randomized fault episodes and deployed behind an NMPC-style supervisory shield. On the acceptance timeline the **open loop runs away to ΔT = 71.6 K** (the 40 K threshold crossed at ~44 min, instantaneous selectivity collapsing toward zero as the hot slug passes the quantification port), while the **autonomous agent alarms predictively at ΔT = 13 K, deploys the dilute-quench inside 15 s, arrests the excursion at ΔT = 14 K — the 40 K threshold is never approached (worst re-ignition peak 21.1 K) — and delivers 28 % more on-spec product** (0.912 vs 0.711) than the uncontrolled plant. **Module 10D** prices the plant: nominal STY **6 751 kg m⁻³ h⁻¹** at X = 98.2 %, S = 95.9 %, with PMI 1.87, E-factor 0.87, cradle-to-gate carbon intensity ≈ 44 kg CO₂e/kg and a unit cost spanning **$612/kg (single channel) → $273/kg (64-channel numbered-up Pareto optimum)** across a 252-point design grid whose cost–carbon–STY Pareto surface is extracted and plotted.

**Headline: an unattended flow plant sails through cavitation, feed impurity and channel clogging — and when the coolant loop vapor-locks, the learned controller detects the incipient runaway from a 13 K tremor, starves and flushes the exotherm within one residence time, and never lets the channel see 40 K of superheat, while the identical uncontrolled plant climbs to 71.6 K and destroys its selectivity.**

---

## 1. From Molecules to Plants: the Third Embodiment

A world model that terminates at a yield prediction is a hypothesis; Phase 9 made those hypotheses falsifiable by a robot; Phase 10 makes them **operable** — the difference between running one experiment and running a plant. Three shifts define the step:

1. **Space becomes a continuum.** Batch flasks (Phase 9's well plates) are spatially lumped objects: one temperature, one concentration vector. A flow reactor is a field: T(z, t), C(z, t), coupled to a coolant field flowing the *other* way. Nothing in Phases 4–8 predicted *where* in a channel a hotspot would grow; the PDE layer is where molecular kinetics first acquire a geography.
2. **Disturbances become operators, not noise.** An industrial plant does not experience Gaussian noise; it experiences *events* — a pump cavitating, a wrong reagent batch, a channel slowly scaling up, a coolant loop losing its pump. Phase 10 injects four such anomalies as explicit state-machine operators on the plant, each with its own physics (slip oscillation, reversible inhibition + acid-promoted side chemistry, diameter shrink + thermal insulation + catalyst burial, jacket-capacity collapse).
3. **Control becomes learning under a safety constitution.** The mission commanded a reward R_t = w₁·Yield + w₂·Selectivity − w₃·Carbon − w₄·𝕀(ThermalRunaway) − w₅·PressurePenalty; the delivered system instantiates exactly that form, trains a maximum-entropy RL agent on it across domain-randomized faults, and then wraps the policy in a model-based supervisory shield — because a policy that has *learned* the physics is still not a device you let hold the pumps at 3 a.m.

```
  Phase 4/5/8 microkinetics ──►  MODULE 10A  ──►  9-field PDE plant  ◄── anomalies (10B injectors)
   (rates, ΔH, ΔΔG‡)              IMEX-MOL              │   ▲                  │
                                                        ▼   │                  ▼
                                               actuators │  sensors      Raman / UV-Vis / ΔP / T
                                                        │   │                  │
                                                  MODULE 10C  ◄── deconvoluted PAT state
                                              SAC policy + NMPC shield       (audit vs truth)
                                                        │
                                                        ▼
                                                  MODULE 10D  ◄── steady-state design grid
                                              STY · PMI · E-factor · CO₂e · $      (TEA/LCA)
```

---

## 2. Module 10A — the Advection–Diffusion–Reaction Continuum

### 2.1 The plant and its chemistry

The reactor is the commanded geometry: L = 10.0 m coil, inner diameter d_t = 1.0 mm (A_cs = 7.854 × 10⁻⁷ m², V = 7.854 mL), wrapped in a counter-current jacket (40 % glycol annulus, A_ann = 5.3 × 10⁻⁶ m², coolant inlet 298.15 K). Three HPLC pumps feed:

- **stream A** — 2-methyl-azirino[1,2-a]indole, 1.80 M in toluene;
- **stream B** — 2,6-lutidine modifier, 0.10 M (the acid scavenger: it suppresses the acid-promoted elimination and cationic oligomerization channels, exactly as a process chemist would dose it);
- **stream Cat** — the designed 3,3′-CF₃-Ph/iPr-Ph-BINOL phosphoric acid, 0.060 M.

Nominal production setpoints 3.60 / 2.40 / 1.20 mL/min → τ = 65.5 s, C_A,0 = 1.11 M, catalyst ≈ 10 mM. The network is the Phase-5 world model, translated to channel units:

| channel | rate law | anchor |
|---|---|---|
| A + Cat → I (RDS) | k₁·C_A·C_Cat | E_a = 88.87 kJ/mol (Phase-5 21.24 kcal/mol); k₁(298 K) = 2.30 M⁻¹s⁻¹ |
| I → P (product) | k₂·C_I | E_a = 59.4 kJ/mol; k₂(298) = 0.18 s⁻¹ (calibrated so [I] stays low and S ≈ 96 %) |
| I → E (elimination) | k₃·C_I × acid/(1+40·C_Lut) | E_a = 110.9 kJ/mol — temperature-gated: ×545 from 293→340 K |
| 2 I → O (oligomer) | k₄·C_I² × 1/(1+25·C_Lut) | E_a = 85 kJ/mol — the avalanche fuel of the runaway |
| Cat → X (burial) | k_d·C_Cat·(1+4·f_foul) | fouling-zone deactivation |

Energetics: A→I −38, I→P −47 kJ/mol — net **−85 kJ/mol**, the aziridine strain release, giving an adiabatic potential ΔT_ad = C_A,0·ΔH/(ρc_p) ≈ 52 K at the nominal feed. Stereochemistry rides on the same anchor as Phase 5: **ee(T) = tanh(ΔΔG‡/2RT)** with ΔΔG‡ = 1.50 kcal/mol → **85.3 % at 298 K**, decaying as the channel heats — the runaway costs enantioselectivity exactly as the microkinetic model says it must.

Steady state (direct PFR march + counter-current coolant iteration, 4 passes): **X = 98.2 %, S = 95.9 %, ee = 85.3 %, T_out = 298.2 K, ΔT_channel = 5.1 K, STY = 6 751 kg m⁻³ h⁻¹** — a numerically boring, thermally tame operating point, which is precisely the point: everything interesting in this phase is what happens when the taming hardware fails.

### 2.2 Meso-scale closures

- **Taylor–Aris dispersion:** D_eff = D_m + a²u²/(192 D_m), D_m = 8.5 × 10⁻¹⁰ m²/s. At nominal velocity (5.6 cm/s) this is ~1.4 × 10⁻³ m²/s — dispersion is real but convection dominates (Pe_L ≈ 400); it is integrated *implicitly* because it is also the stiffest diffusion in the problem.
- **Hydraulics:** dP/dz = 128 μ(T) Q / (π d_eff(z)⁴) with log-linear toluene viscosity; the fouling constricts d_eff locally to 0.52 mm, and ΔP rises as d⁻⁴ (0.30 → 0.5–0.8 bar across the campaign — the PAT pressure panel's signal).
- **Heat transfer:** series resistance 1/h_i + R_wall + 1/h_c with Graetz h_i = 3.66 k_l/d_t ≈ 479 W/m²K and h_c ∝ u_c^0.55; overall U ≈ 293 W/m²K nominal → UA·V-basis cooling capacity ≈ 1.07 W/K. **This number is the whole safety story**: nominal reaction power is ~6 W (ΔT_wall < 1 K), but when the vapor lock cuts jacket flow to 6 %, capacity falls to 0.064 W/K and the channel is effectively adiabatic — the 52 K adiabatic potential becomes reachable, and the Arrhenius avalanche (k₁ ×80 per +40 K) does the rest.

### 2.3 The IMEX Method-of-Lines integrator

Nine fields (A, I, P, E, O, Cat, Lut, Imp, T) plus the counter-current coolant T_c are discretized on N_z cells (161 hi-fi / 41 training / 101 validation) and advanced by operator splitting:

1. **Implicit axial diffusion** — trapezoidal banded solve with mirror (Danckwerts) closure at both ends; unconditional stability against the Taylor–Aris stiffness. *(An early draft had 1+λ boundary diagonals instead of 1+2λ — a mass-creating mirror; the conservation bug was found by a channel mass-balance audit and fixed before any production run.)*
2. **Explicit TVD-RK2 transport + reaction + exchange** — MC-limited slopes, conservative fluxes, inlet Dirichlet, coolant advected in −z with its inlet at z = L. The time step obeys three adaptive constraints: advective CFL (0.7·Δz/u), coolant CFL, and a **depletion-aware chemical-stiffness bound** dt ≤ 0.3/k_eff where k_eff is computed from the *local* decay constants r₁/C_A, (r₂+r₃+2r₄)/C_I — so the hot but exhausted post-frontal zone does not tax the integrator, only the moving reaction front does. During the uncontrolled runaway this holds the explicit scheme inside the physical adiabatic bound (the twin's early "104 K excursion" was found to be *correct physics* of an oligomer-dominated route, not a solver artifact — and the oligomer enthalpy was re-priced to its per-bond value, −40 kJ/mol, on that occasion).

### 2.4 What the continuum teaches (Fig. 1)

The spatio-temporal maps expose structure no lumped model can show: the conversion front steepens after the impurity spike (reversible inhibition, k₁/(1+220·C_Imp)); the fouling band at z = 6.4–7.6 m appears as a visible thermal ridge once its U collapses (the deposit is a *blanket*, not just a constriction); and after the vapor lock the superheat blooms in the *inlet half* of the coil — counter-current logic: fresh coolant enters at z = L = the outlet end, so the z ≈ 0 end sits on the oldest, warmest coolant and runs away first. Fig. 1(d) shows the axial profiles stepping through this growth, under the controlled run where the growth never completes.

---

## 3. Module 10B — Operando Multi-Modal PAT

### 3.1 Forward sensor models

- **In-line Raman (785 nm):** every analyte carries Lorentzian fingerprint bands (A: 1602/1265 cm⁻¹ indole & aziridine; P: 1655 cm⁻¹ imine — the product's growing fingerprint; E: 1668 cm⁻¹, deliberately overlapping P to make deconvolution honest; impurity: 1003/1608 cm⁻¹), multiplied by an Ornstein–Uhlenbeck laser-power fluctuation (3 %, τ = 40 s) on top of a rising exponential fluorescence baseline that *jumps at the impurity event and decays for minutes* — the classic real-world baseline pathology.
- **UV-Vis photodiode pair (254/310 nm):** Beer–Lambert over a 1.0 mm cell with per-species ε (P: 15 200 / 8 600 M⁻¹cm⁻¹), shot noise and slow drift.
- **ΔP + thermometry at 1 Hz:** the Hagen–Poiseuille integral with the *fouled* diameter profile, plus port temperatures at z = 2.5, 5.0, 10.0 m.

### 3.2 The anomaly schedule

| t (s) | anomaly | physics injected |
|---|---|---|
| 600–750 | pump-B cavitation | delivered Q_B = command ×(1 − slip), slip 30 % oscillating at 1.3 Hz; raising P_BPR raises suction pressure and collapses the slip (relief ∝ 1−0.09·ΔP) |
| 1200–1500 | impurity spike | stream-A purity −6 %; benzoic acid 12 mM → k₁ /(1+220·C_Imp) (reversible RDS inhibition) and k₃ ×(1+150·C_Imp) (acid-promoted elimination); fluorescence baseline jump |
| 1800 → end | progressive fouling | d_eff → 0.52 mm at z ≈ 7 m (ΔP ×~2.7 locally), U/(1+12·f) thermal blanket, catalyst burial k_d ×5 |
| 2400 → end | coolant vapor lock | jacket flow command ×0.06 over a 240 s ramp → cooling capacity 1.07 → 0.064 W/K → the channel goes effectively adiabatic |

### 3.3 The deconvolution agent and its audit

The controller never sees simulator truth. At every control tick (5 s) a 6-measurement least-squares inverter (UV pair + four Raman peak heights) estimates {A, I, P, E, Imp}; the estimate is seeded with 4 % noise and a +2 % bias on P. Post-run audit against plant truth across the campaign: **mean relative error 0.53 %, max 3.64 %** (n = 109 port samples × species) — the Phase-9 hallucination-audit discipline, transplanted from chromatograms to spectra. The full synthetic telemetry bundle (601-point Raman spectra at 10 s cadence, waterfall-ready) is archived for Fig. 2.

---

## 4. Module 10C — Cyber-Physical Reinforcement Control

### 4.1 State, action, reward

- **State (18-dim):** deconvoluted A/I/P concentrations, three port temperatures, ΔP and its trend, the five actuator setpoints and two slew deltas, max channel superheat ΔT and its rate, vapor-pressure margin (Antoine).
- **Action (5-dim continuous, tanh-squashed):** Q_A, Q_B, Q_cat ∈ [0.01, 5.0] mL/min; jacket rate ∈ [0.2, 3.0] × nominal; P_BPR ∈ [1, 30] bar. Rate-limited (slew) inside the plant.
- **Reward (the commanded form):** R_t = 1.0·Y_inst + 0.6·Sel_inst − 0.35·Carbon_norm − 4.0·𝕀(ΔT>40)·severity − 0.8·max(0, (ΔP−2.5)/2.5)², with Yield/Selectivity instantaneous at the outlet and the carbon proxy spanning pump electricity, chiller duty (COP 3.2) and catalyst make-up at 85 kg CO₂e/kg.

### 4.2 Training

Soft Actor-Critic (twin Q, squashed-Gaussian actor initialized *at the nominal setpoint via tanh-space bias targets*, automatic entropy temperature, γ = 0.985) on the 41-cell plant, 110 episodes × 1 600 s over eight domain-randomized fault families (none / cavitation / impurity / fouling / coolant / coolant+fouling / cavitation+impurity / all), fault times and amplitudes jittered. Episode return 285.8 → 293.5 with α annealing 0.90 → 0.22; mean peak-ΔT across training episodes settles ≈ 15–25 K with the shield absent — the raw policy already avoids, but does not guarantee, the threshold.

### 4.3 The NMPC-style shield — and why the quench is *dilution*

The deployed controller is **SAC proposal → model-based shield → plant**. The shield is a three-state machine (ARMED → QUENCH → RAMP) built from one-step thermal logic:

- **predictive trip:** ΔT > 26 K absolute, or ΔT > 13 K with d(ΔT)/dt > 0.05 K/s — the incipient-excursion signature of the 240 s vapor-lock ramp;
- **dilute-quench:** Q_A → 1.2, Q_B → 4.5, Q_cat → 0.05 mL/min, jacket → max, BPR ≥ 8 bar. The physics matters: cutting the catalyst dose alone is *useless* once the channel is hot (at 350 K even 0.5 mM catalyst converts a pass); tripling the substrate feed would *amplify* the exotherm. The only fast lever is **stoichiometric dilution** — C_A,0 falls 3×, the adiabatic potential falls to ~12 K, and the cold feed flushes the channel in ~1 τ — plus maximum jacket demand, which partially restores the saturated coolant capacity;
- **arrest + safe holding ramp:** when the rate turns negative below 38 K the excursion is declared arrested; actions then ramp to a reduced-throughput holding point (1.8/3.0/0.22/2.5) — *not* back to nominal, because the vapor lock has not been observed to clear — with the trip rule still armed (re-ignition → instant re-quench);
- **guards:** ΔP > 2.6 bar throttles Q_A with a d^¼ sensitivity and fires a 90 s catalyst-rich antifouling pulse; vapor-pressure margin < 2 bar raises the BPR.

### 4.4 The acceptance timeline (Fig. 3)

| | open loop (nominal setpoints) | SAC + shield (autonomous) |
|---|---|---|
| cavitation 10–12.5 min | selectivity sags with the slip oscillation | BPR stepped up, slip collapsed, Q_B re-trimmed |
| impurity 20–25 min | conversion dip (rate ×0.46), recovery only at window end | catalyst re-dosed through the inhibition, ΔT excursion absorbed |
| fouling 30 min → | ΔP climbs to ~0.5 bar, thermal ridge at z ≈ 7 m | ΔP guard + antifouling pulse; ridge held at ~13 K |
| **vapor lock 40 min →** | **ΔT crosses 40 K at ~43.7 min, plateaus at 71.6 K; instantaneous selectivity collapses to ≈ 0 as the hot slug passes the port; 28.9 % of cumulative product is off-spec** | **predictive alarm at 42.8 min (ΔT = 13 K), quench deployed ≤ 15 s, excursion arrested within 10 s of alarm at ΔT = 14.0 K; worst re-ignition peak 21.1 K at 57 min — the 40 K threshold is never approached** |
| **on-spec yield (S ≥ 90 %)** | **0.711** | **0.912 (+28.3 %)** |
| cumulative yield | 0.905 | 0.912 |

The open-loop row is the mission's warning made quantitative: *flow chemistry is intrinsically safe at the designed heat transfer — and can still run away when the cooling utility fails*. The controlled row is the thesis: the same physics, read by spectrometers and acted on within one residence time, is a non-event.

### 4.5 Validation suite (held-out scenarios, N_z = 101)

Eight held-out fault scenarios × three controllers:

| controller | peak ΔT (K) | ΔT > 40 K | final yield | tail selectivity |
|---|---|---|---|---|
| open loop | 22.4 ± 14.8 | **1 / 8** (coolant: 60.5 K) | 90.4 % | 96.4 % |
| SAC only | 11.4 ± 7.6 | 0 / 8 | **93.0 %** | 97.2 % |
| SAC + shield | **9.7 ± 3.0** | **0 / 8** | 92.6 % | **97.2 %** |

The learned policy alone already halves the excursions and beats open-loop economics (it learned, unprompted, that gentler catalyst dosing is cheaper per kg); the shield adds the guarantee — tightest peak spread, zero breaches, at a 0.4-point yield cost on the coolant scenario, which is the honest price of a safety interlock.

---

## 5. Module 10D — Techno-Economics & Lifecycle

### 5.1 The design grid

252 dynamic steady states (direct PFR march + counter-current jacket iteration, validated against the transient twin) over 4 numbering-up scales (1 / 8 / 32 / 64 channels) × 7 stream-A flows (0.5–5.0 mL/min) × 3 jacket temperatures (10 / 15 / 20 °C) × 3 catalyst loadings. Each point is priced cradle-to-gate: substrate (38 kg CO₂e/kg, $220/kg), toluene with 87 % in-loop recovery (+0.55 kWh/kg recovery energy), lutidine, catalyst at 5 % per-pass loss (scavenger loop, 85 kg CO₂e/kg, $5 300/kg), pumping + chilling electricity at 0.5542 kg CO₂e/kWh, capex ($320 k + $14 k/channel, 8-year depreciation at 85 % utilization) plus maintenance and residual labor.

### 5.2 What the Pareto surface says (Fig. 4)

- **The nominal champion point:** STY 6 751 kg m⁻³ h⁻¹ (single channel) — X = 98.2 %, S = 95.9 % — at $434/kg, 44.1 kg CO₂e/kg, PMI 1.87, E-factor 0.87. The isomerization's 100 % atom economy is visible in the E-factor: solvent make-up, not stoichiometry, is the waste.
- **The extracted Pareto front (5 non-dominated points)** runs along the numbering-up ceiling: 64 channels at minimum catalyst loading, stream-A flow 1.0 → 4.0 mL/min buys STY 1 934 → 4 471 kg m⁻³ h⁻¹ at the price of $277 → $415/kg and 43.4 → 71.5 kg CO₂e/kg. Beyond ~3 mL/min per channel the selectivity losses (and the carbon they carry) grow faster than the throughput.
- **Cost anatomy:** single-channel $612/kg vs 32-channel $280/kg — capex dilution and catalyst amortization dominate the descent; substrate feed is the irreducible floor ($220/kg × 1/0.98 conversion ≈ $222/kg) in every bar of the stacked anatomy.
- **Carbon anatomy:** ≈ 86 % of the 44 kg CO₂e/kg is the multi-step substrate itself; energy is ~0.04 kWh/kg — the flow plant's chiller barely notices the reaction. Green-chemistry leverage on this process lives in the substrate's own synthesis route, not in the reactor.

---

## 6. What the Twin Teaches

1. **Counter-current jackets fail asymmetrically.** The vapor-lock runaway ignites at the z ≈ 0 end — the coolant's *exit* — because fresh coolant enters at the outlet. A co-located single-point thermocouple at the outlet would have missed the precursor by 100 K of margin; the three-port thermal PAT array is what buys the predictive alarm.
2. **Catalyst starvation is not a quench.** The avalanche regime converts a pass regardless of dose; only dilution (lower C_A,0 → lower ΔT_ad) plus cold-flush removes the heat. This is a design rule the shield encodes and a rule a human operator would want printed on the panel.
3. **Safety costs yield — measurably, and acceptably.** Shield vs SAC-only on the coolant scenario: 17.7 K vs 31.4 K peak, at −3.4 points of yield. The validation suite prices the interlock instead of asserting it.
4. **Dispersion is a first-class citizen.** Taylor–Aris D_eff ≈ 1.4 × 10⁻³ m²/s at production velocity is three orders above molecular D; ignoring it would misstate both the conversion front and the thermal ridge.

---

## 7. Roadmap: Toward the Autonomous Plant-Model Loop

- **Identify, don't just control:** an ensemble-Kalman or physics-informed neural observer on the 9-field state would turn the shield's thresholds into fault *estimates* (vapor-lock severity, fouling fraction) rather than alarms.
- **Sim-to-real transfer:** the SAC was trained on the same twin it controls (with randomized faults); the credible next step is training against the Phase-9 laboratory twin and deploying on hardware, with the hallucination audit as the transfer metric.
- **Multi-objective online economics:** Module 10D's Pareto surface is computed off-line; folding the cost/carbon gradients into the reward would let the plant *walk its own Pareto front* as catalyst degrades.

---

## 8. Reproduce

```bash
python run_phase10_cyberphysical_flow_twin.py             # full campaign (~35 min: 110-episode SAC
                                                          #  training, 24 validation episodes, 3 hi-fi
                                                          #  acceptance runs, 252-point TEA grid, figures)
python run_phase10_cyberphysical_flow_twin.py --selftest  # ~4 min smoke test (16 episodes, coarse grids)
python run_phase10_cyberphysical_flow_twin.py --fig_only  # re-render figures from saved episode archives
```

Dependencies: numpy / scipy / matplotlib / scikit-learn (repo `requirements.txt`) **plus `torch`** (CPU sufficient; 4 threads). All outputs land in `results_phase10/` (master JSON, episode field archives `episode_*.npy`, PAT telemetry archive, SAC learning curve, validation summary, TEA point cloud) and `figures_phase10/` (fig1–fig4, 300 DPI).

---

## 9. Transparency Ledger

- **Everything in this phase is simulation.** The plant, sensors, faults, robot-layer actuators and economics are models; no hardware was driven. The chemistry anchors (barrier, thermodynamics, stereodifferentiation, catalyst identity) are the Phase-4/5/8 computed values carried forward unchanged.
- **Thermodynamic re-pricing during development:** the oligomerization enthalpy was initially set to −118 kJ/mol per oligomer; the first uncontrolled runaway showed the twin honoring that number with a ~100 K excursion. The value was re-priced to −40 kJ/mol per formed bond (per-bond strain/aromatic energetics) — both the bug-hunt that proved the integrator conservative and the enthalpy correction are part of the record.
- **Two conservation bugs were found and fixed before the production run** (mirror-boundary diffusion diagonals; a 10³ unit slip in the reaction-heat term), both caught by audits (channel mass balance; SS-vs-transient duty comparison) rather than by eyeballing plots.
- **The deployed controller is safety-filtered RL.** The SAC proposes; the shield disposes. Every autonomous action in the acceptance timeline was either a policy output inside the shield's envelope or a shield override of one; the two are logged separately in the episode records.
- **The "arrest within 10 s of alarm" figure** is measured under the commanded trip rule (predictive ΔT = 13 K + rate); the full stabilization to the safe holding point takes ~1 residence time, and worst-case re-ignition peaks at 21.1 K — both numbers are reported rather than the flattering one alone.
- **The TEA emission factors and prices** are literature-anchored estimates (Ecoinvent-class ranges, 2024 grid factor 0.5542 kg CO₂e/kWh), not measured values; the Pareto *ordering* is robust to them, the absolute dollars are not.
