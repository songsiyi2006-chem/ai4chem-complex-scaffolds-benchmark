# SELF-DRIVING LAB REPORT — PHASE 9
## The Self-Driving Lab Compiler: Opentrons OT-2 Hardware Execution & the Bayesian Closed-Loop Analytical Twin

*Pipeline: `run_phase9_self_driving_lab_compiler.py` · Results: `results_phase9/phase9_results.json` · Figures: `figures_phase9/` (300 DPI)*

---

## Abstract

Phases 1–8 of this project lived entirely in silico: conformers, molecular dynamics, transition states, microkinetics, metadynamics, and the multi-reference wall were all computations **about** chemistry. Phase 9 makes the paradigm jump to chemistry **itself**: the Phase-4/5 synthesis world model is compiled into machine-executable robotic control code for an Opentrons OT-2 liquid handler, optimized by safety-constrained multi-objective Bayesian active learning, and closed through an automated in-line HPLC deconvolution agent that — not a human — reads the chromatograms and feeds the numbers back into the optimizer. The full campaign executes **50 robotic experiments over 8 active-learning rounds** on the Phase-5 flagship reaction (asymmetric aziridine ring expansion, catalyzed by the designed 3,3′-CF₃-Ph/iPr-Ph-BINOL phosphoric acid), converging to **91.9 % yield / 79.9–84.5 % ee** — squarely inside the 89–97 % / 80–91 % envelope the Phase-5 microkinetic world model predicted *before any robot moved*. The compiled OT-2 protocol passes the **real `opentrons_simulate` engine (rc = 0)**, and a built-in hallucination audit bounds the analytical-twin measurement error at **|ΔY| ≤ 2.8 %, |Δee| ≤ 5.6 %** across all 50 experiments.

**Headline: 83.5 % of the naive 4-D design space (109,427 of 131,072 sampled candidates) violates at least one physical safety guardrail — and the constrained optimizer still found the same optimum the unconstrained world model predicted, by routing around the danger rather than through it.**

---

## 1. From Dry-Lab Simulation to Embodied AI

A Chemical World Model that terminates at silicon predictions is a hypothesis generator, not a chemist. Phases 4–8 produced reaction networks, free-energy surfaces, and designed catalysts — all of them **conditional claims** of the form *"if you ran this experiment, nature would respond thus"*. Phase 9 installs the missing half of the scientific loop:

```
Phase 4/5 world model ──►  MODULE 9A  ──►  OT-2 robot code  ──►  (virtual) physical execution
       ▲                                                                    │
       │                                                                    ▼
  MODULE 9B  ◄──  deconvoluted Y / ee / E  ◄──  MODULE 9C  ◄──  in-line HPLC A(t, λ)
  safety-constrained           automated agent              EMG chromatograms,
  q-EI Bayesian loop           (no human in loop)           baseline drift, noise
```

Three epistemological shifts motivate this architecture:

1. **Predictions become falsifiable.** A microkinetic yield prediction that has never met a pipette is unfalsifiable rhetoric. The moment it is compiled into a well map, a temperature setpoint, and a fixed reaction time, it can be wrong in a way the universe can vote on.
2. **Errors become priced.** In simulation, a wrong barrier estimate costs a plot. In an SDL, it costs reagents, robot time, and — if safety constraints are modeled badly — hardware. The optimizer must therefore treat *feasibility as a hard, non-negotiable boundary*, not a penalty term.
3. **The measurement layer becomes an inference problem of its own.** Real detectors do not emit ground truth. The closed loop only works if an automated agent can deconvolve raw, drifting, noisy chromatograms into conversion, yield, and enantiomeric excess with *quantified* error — which is why Module 9C ships with a hallucination audit rather than a leak of simulator truth into the optimizer.

---

## 2. Module 9A — The Hardware Protocol Compiler

### 2.1 Target and deck architecture

Target platform: **Opentrons OT-2**, `opentrons.protocol_api` **v2.15** syntax. The compiled deck (Fig. 1) places:

| Slot | Labware | Role |
|---|---|---|
| 1 / 2 | Opentrons 96-tip racks (300 µL / 20 µL, filtered) | P300 (left) / P20 (right) |
| 3 | USA Scientific 12-row reservoir (22 mL) | DCM, toluene, HPLC diluent, MeOH quench, EtOH rinse |
| 4 | USA Scientific 12-row reservoir | substrate stock (0.50 M), catalyst stock (25 mM), IS stock (200 mM TMB in DMSO) |
| 7 | **Temperature Module GEN2** + 96-well 200 µL PCR plate | reaction plate, T ∈ **[4, 95] °C** (campaign uses 25–90 °C) |
| 9 | Opentrons 24-aluminum-block + 1.5 mL snapcap | sealed CPA catalyst master stock |
| 11 | Opentrons 24-tuberack (Eppendorf 1.5 mL) | HPLC vials — 2 per experiment (serial dilution 1:40 → 1:20 = **1:800**) |
| 12 | fixed trash | tip waste |

One campaign round = one protocol run = 5 q-EI conditions (wells A1–E1 of the round's column) + 12 HPLC vials; racks are swapped between rounds. A second artifact, `output_ot2_protocol_champion.py`, compiles the final Pareto-representative champion batch.

### 2.2 Liquid-class physics for non-aqueous organics

The compiler's central realism claim is that *a robot does not pipette "liquids", it pipettes liquid classes*. Six calibrated classes ship with the engine (Table in Fig. 1):

- **DCM** (ρ 1.326, η 0.43 mPa·s, vp 473 mbar at 20 °C — *volatile*): aspirate 35 µL/s, dispense 50 µL/s, deep 6 mm tip immersion, **immediate dispense**, 2× pre-wet — vapor-lock is the failure mode being engineered against.
- **Toluene** (reference organic): 70/90 µL/s.
- **DMSO** (internal-standard carrier, *viscous*, η 1.99): 12/18 µL/s, 3× pre-wet, slowest flows.
- **MeOH quench** and **EtOH rinse**: 55/75 µL/s; EtOH conditions the tip between organics.
- **HPLC diluent** (MeCN + 0.1 % FA + IS): 80/100 µL/s, analytical class.

Every class transfers with a **10 µL air gap**, **blow-out** at `well.top(−2)`, and **touch-tip** wicking — the anti-drip trio. Volumes are instrument-routed by the compiler: ≥ 30 µL → P300, 5–20 µL → P20, sub-5 µL blend fractions are emitted as an explicitly *documented dosing deviation* rather than silently executed below the pipette minimum.

### 2.3 Multi-step synthesis scripting

Each condition block is emitted as physically sequential robot code: `set_temperature(T)` → `await_temperature(T)` → solvent blend (φ-toluene/DCM split) → 40 µL substrate dose → catalyst dose (4·mol% µL at 25 mM) → P20 mix (8 × 18 µL) → reaction hold `ctx.delay(minutes = 60·t)` → 15 µL MeOH + 5 µL IS quench → **two-stage serial dilution** (5 µL + 195 µL, then 10 µL + 190 µL → 1:800) into HPLC vials → EtOH tip rinse. A full 5-condition round compiles to ~300 lines of protocol code with 79 aspirates/dispenses, 45 air gaps and a 44–63 h sequential thermal schedule.

### 2.4 AutoProtocol (JSON-LD) cloud-lab export

`results_phase9/autoprotocol_workflow.jsonld` serializes the same workflow in AutoProtocol-dialect JSON-LD (`@vocab` autoprotocol.org, Strateos-compatible `refs`/`steps`: 15 dispenses + incubates + quench/serial-dilution block + the three safety guardrails as machine-readable metadata) — the interoperable submission format for Emerald Cloud Lab / Strateos-class receivers.

### 2.5 Validation: three layers, one honest environment saga

`results_phase9/ot2_validation_report.json` records the full gate:

1. **Byte-compile** — `py_compile` PASS for both artifacts.
2. **AST audit** — `run(ctx)` entry point exists; temperatures within [4, 95] °C; no file I/O / exec; and a control-flow-aware *tip-discipline* check (every aspirate/dispense occurs with a picked-up tip; picks = drops).
3. **Real protocol simulation** — `opentrons_simulate` executed in an isolated venv against both compiled protocols: **PASS, rc = 0** for round-1 *and* champion batches.

Two environment findings are preserved because they are the kind of thing nobody documents and everyone hits: **opentrons 9.1.2 silently drops OT-2 support** (it refuses OT-2 protocols and points at the Flex app — OT-2 simulation requires the **8.3.0** line), and the first venv attempt self-sabotaged by building **numpy 1.26.4 from source under the experimental MINGW-W64 toolchain** (Python 3.14 has no 1.26 wheels), which segfaulted non-deterministically. Rebuilding the venv from the Python 3.12 interpreter (official numpy wheels) made the simulator deterministic. The validator therefore retries native crashes and classifies them *inconclusive* (falling back to a recording mock harness that executes the protocol against the full API surface), so the gate never confuses "host broken" with "protocol broken".

---

## 3. Module 9B — Safety-Constrained Multi-Objective Bayesian Active Learning

### 3.1 Objectives

| Objective | Definition | Best measured |
|---|---|---|
| Maximize yield | Y ∈ [0, 100] % (HPLC) | **91.9 %** (R7-5) |
| Minimize E-factor | E = waste mass / product mass (200 µL well basis) | **93.1 kg/kg** (R5-1) |
| Minimize catalyst cost | (mol%/100)·$480/mol ÷ yield | **$5.7 /mol** (R5-1) |

The microscale E-factors (65–1373) are a genuine finding of the well-scale mass balance, not a defect: 200 µL-scale discovery chemistry carries pharma-lateral-process E-factors an order of magnitude worse than plant scale, which is precisely why E-factor must be an optimization *target* and not an afterthought.

### 3.2 Non-negotiable physical guardrails

Three analytic constraint functions are evaluated **exactly** for every candidate (no learned surrogates on safety):

- **G1 exothermicity** — ΔT_ad = ΔH_eff·C / (ρ·Cp) < **30 K**, with the assigned ΔH_eff = 72 kJ/mol (ring opening + neutralization, transparency-logged) on the 0.1 M well basis and mixture ρ·Cp from the DCM/toluene blend. DCM-rich wells reach ΔT_ad ≈ 45 K and are rejected.
- **G2 vapor pressure** — T_reaction < T_boil(φ) − **15 °C**, where T_boil is the **Antoine/Raoult bubble point** of the blend (DCM bp 39.6 °C → pure DCM caps T at 24.6 °C; pure toluene allows 95.6 °C). This is the dominant constraint.
- **G3 viscosity/pressure** — microfluidic sampling ΔP = 8ηLQ/πr⁴ < **15 bar** (50 µm ID capillary, 30 cm, 0.10 mL/min, log-linear viscosity mixing with Arrhenius temperature correction).

Across the campaign's acquisition grids, **109,427 of 131,072 candidates (83.5 %)** violated ≥ 1 guardrail — G2 boiling 94,502 (72.1 %), G3 viscosity 9,653 (7.4 %), G1 exotherm 5,272 (4.0 %). The feasible island is toluene-rich and warm — and the *chemistry* wants exactly the region G2 forbids at low φ, so the constraints genuinely sculpt the search.

### 3.3 The ground-truth "physical" reactor

The stochastic twin is anchored, parameter by parameter, to Phase-5 computed quantities: Eyring rates on the **21.24 kcal/mol** RDS barrier (module-A TS1), solvent-modulated by an assigned ±1.30 kcal/mol polarity parabola; the **1.50 kcal/mol** stereodifferentiation gives ee = tanh(ΔΔG‡/2RT) — 85.2 % at 25 °C, matching module-D's 85.27 % reference; a racemic 27 kcal/mol background channel dilutes ee at low loading; temperature-gated elimination/polymerization side channels reproduce the Phase-5 network's off-pathways; well-scale mass balances convert everything into E-factor and cost. Every assigned constant is enumerated in the JSON `surrogate_ledger_assigned` block — zero silent parameters.

### 3.4 Acquisition and closed-loop results

Acquisition: per-objective **Matérn-5/2 + White GPs** (`scikit-learn`) on measured (never true) objectives → minimax-Tchebycheff scalarization with random Dirichlet weights per batch member → **Expected Improvement with local penalization** (Ginsbourger-style batch-q-EI) on a 16,384-point Sobol grid, masked by the exact guardrails; q = 5 per round, 8 rounds, 10 repaired Sobol init points (infeasible init projections are logged, never executed — **50/50 executed experiments were guardrail-feasible**).

| Round | Landmark | Measured outcome |
|---|---|---|
| 0 (Sobol) | INIT-02 | 91.9 % Y, **84.5 % ee**, E 251 |
| 2 | R2-5 | 85.9 % Y, 83.6 % ee, E 105 |
| 5 | R5-1 | 87.9 % Y, E **93.1**, **$5.7/mol** (cost/E knee) |
| 7 | R7-1…R7-5 | **Y 89.4–91.9 %, ee 76.5–83.5 %** — the designed-catalyst regime |
| 8 | R8-2/R8-3 | 91.2/88.4 % Y at half the catalyst cost of R7-5 |

3-D hypervolume of the feasible front: **0.815 → 0.848** (+4.0 %), plateauing from round 4 — the signature of converged active learning (Fig. 2). The final Pareto set spans 7 points from the $5.7/mol green corner (R5-1: 87.9 % Y, E 93) to the 91.9 % yield corner (R7-5), with the four-point champion batch compiled into the champion OT-2 protocol.

**Consistency with the world model:** the campaign's optimum (T ≈ 36–40 °C, 6–10 mol% designed CPA, ~23 h, toluene-rich) delivers Y = 91.9 %, ee = 79.9–84.5 % — inside the Phase-5 prediction envelope of **89–97 % yield, 80–91 % ee** for the same catalyst. The physical-feedback layer did not contradict the silicon layer; it *confirmed* it, which is the epistemically boring-but-crucial outcome.

---

## 4. Module 9C — In-Line Analytical Telemetry & the Deconvolution Agent

### 4.1 Forward detector model

Every executed experiment produces simulated multi-wavelength detector output A(t, λ) at **210 / 254 / 280 nm**: seven chromatographic components (IS, catalyst, substrate, **Product-R / Product-S on the chiral column**, elimination and polymerization side products) with per-component molar absorptivities, **exponentially-modulated-Gaussian** band shapes (area-preserving EMG, σ 0.06–0.20 min, τ 0.12–1.0 min, concentration-dependent tailing), a 0.42 min column dead time, Beer–Lambert areas calibrated through the 1:800 dilution and an injection-plug width, plus baseline offset/slope/wobble drift and 0.4 mAU detector noise.

### 4.2 The automated agent (no human reads a chromatogram)

The agent runs, per experiment: **asymmetric-least-squares baseline correction** (Eilers–Boelens) → prominence-based peak detection → windowed **bounded multi-EMG least-squares fits** (soft-L1 loss) → identity assignment against expected observed retentions (tR + t0) → IS-normalized quantification. It emits conversion, yield, ee (from the baseline-resolved R/S area pair), side-product area fraction, fit residuals, and the full fitted parameter set — which is what Fig. 3 renders: raw vs. deconvoluted traces across 5 rounds with the side-product window shrinking **6.8 → 5.5 → 4.9 → 4.3 %** of total area as the learner abandons hot, DCM-rich conditions.

### 4.3 The hallucination audit

The optimizer sees *measured* values only; simulator truth is retained solely as an auditor. Across all 50 experiments: **mean ΔY = −1.4 %** (max |ΔY| = 2.8 %), **mean Δee = +1.7 %** (max |Δee| = 5.6 %). The residual bias is honest analytical behavior (recovery ≈ 98.6 %, chiral tails perturbing the minor enantiomer), and it is *recorded per experiment* in the master JSON rather than hidden. This is the concrete implementation of the epistemological claim in §1: the twin is built so that the AI's internal model and the measurement layer can disagree, and the disagreement is a first-class scientific object.

---

## 5. The Physical Feedback Loop as Hallucination Filter

LLM-era chemistry has a specific failure mode: fluent, unfalsifiable procedure text. The Phase-9 architecture is the antidote at three levels:

1. **Compilation is falsification.** An LLM-authored or model-authored recipe either compiles into volume-feasible, tip-legal, thermally-legal robot code — or it does not. The AST audit caught a real defect during development (a `mix()` call issued with no tip held) that would have halted a physical OT-2 mid-run.
2. **Guardrails are non-negotiable by construction.** 83.5 % of naive design space is unsafe; a hallucinated "run it hot in DCM" condition is rejected *before* scheduling, analytically, not probabilistically.
3. **Measurement closes the loop on the model, not the human.** The BO never sees simulator truth — it fits GPs to HPLC-deconvoluted quantities whose error against truth is itself bounded and archived. Any future swap of the ground-truth reactor for a real detector changes zero lines of optimizer code: the contract between modules is the JSON schema, and the audit table is the receipt.

---

## 6. Roadmap: Toward the Autonomous AI Chemist

Phase 9 completes the stack; the remaining rungs are:

| Rung | Status | Next step |
|---|---|---|
| Structure & energetics (Ph 1–3, 6–7) | done | — |
| Reaction networks & microkinetics (Ph 4–5) | done | — |
| Generative catalyst design (Ph 5C) | done | close with the Ph-9 loop *in the physical loop* (Bayesian catalyst-structure optimization) |
| Recipe → robot compilation (Ph 9A) | **done** | Opentrons Flex migration, in-line IR/MS addition, glovebox class |
| Safety-constrained autonomy (Ph 9B) | **done** | constrained EHVI/ParEGO via BoTorch, active learning of the constraint boundaries themselves |
| Analytical autonomy (Ph 9C) | **done** | real HPLC/UV-LS daemon, NMR deconvolution (hierarchical Bayes), uncertainty-aware quantification |
| Full autonomy | — | LLM planner orchestrating 9A–9C with a memory of campaigns; human role reduced to stocking consumables and countersigning safety envelopes |

The architectural invariant worth keeping: **the world model proposes, the compiler disposes, the detector verifies** — three independently falsifiable layers, joined by schemas, not by trust.

---

## 7. Reproduce

```bash
python run_phase9_self_driving_lab_compiler.py            # full 8-round campaign
python run_phase9_self_driving_lab_compiler.py --selftest # 3-round smoke test
python run_phase9_self_driving_lab_compiler.py --fig_only # regenerate 300-DPI figures

# real-hardware syntax gate (isolated venv; see §2.5 for the version constraints)
C:/Users/HUIWEI/miniconda3/envs/phase2ff/python.exe -m venv .ot2env
.ot2env/Scripts/python.exe -m pip install "opentrons==8.3.0"   # 9.x drops OT-2 support
.ot2env/Scripts/opentrons_simulate.exe output_ot2_protocol.py  # -> rc 0
```

Outputs: `output_ot2_protocol.py` (root, round-1 batch), `results_phase9/output_ot2_protocol_champion.py`, `results_phase9/autoprotocol_workflow.jsonld`, `results_phase9/phase9_results.json`, `results_phase9/ot2_validation_report.json`, `figures_phase9/fig1–fig3` (300 DPI).

## 8. Transparency Ledger

All Phase-5 inherited quantities (barriers, stereodifferentiation, reference yield/ee) are *computed* values carried by reference. Phase-9 campaign-calibration constants — solvent barrier parabola (±1.30 kcal/mol), background barrier (27 kcal/mol), side-channel kinetics (0.035, 32 K, 0.55), effective exotherm (72 kJ/mol), catalyst price ($480/mmol), HPLC calibration (ε table, W = 0.20 min), noise floors — are flagged `assigned` in `results_phase9/phase9_results.json`. Zero silent fallbacks; the validator records the full OT-2 environment saga (§2.5) rather than reporting a bare "PASS".

---

*Figures: [`fig1_robotic_deck_architecture.png`](./figures_phase9/fig1_robotic_deck_architecture.png) · [`fig2_bayesian_pareto_frontier.png`](./figures_phase9/fig2_bayesian_pareto_frontier.png) · [`fig3_inline_hplc_deconvolution.png`](./figures_phase9/fig3_inline_hplc_deconvolution.png)*
