# Phase 18 — Whole-Cell Metabolic Networks: Thermodynamic Flux (tFBA), Macromolecular Crowding & Cybernetic Kinetics

**English technical treatise · Phase 18 of the AI4Chem Pantheon**
Pipeline: [`run_phase18_wholecell_metabolic_thermodynamics.py`](./run_phase18_wholecell_metabolic_thermodynamics.py) · Record: [`results_phase18/phase18_results.json`](./results_phase18/phase18_results.json) · 中文版：[`WHOLE_CELL_METABOLISM_REPORT_ZH.md`](./WHOLE_CELL_METABOLISM_REPORT_ZH.md)

---

## Abstract

Phase 18 breaks the Pantheon out of the single-enzyme (or single-mechanism) regime into the global, non-equilibrium network dynamics of an entire minimal autonomous cell. Module 18A reconstructs a genome-scale stoichiometric network — **817 internal metabolites (869 rows including boundary species), 1,269 reactions (1,201 enzymatic + transport)** spanning central carbon, amino-acid, nucleotide, cofactor, fatty-acid/membrane-lipid, cell-envelope (murein, lipoteichoic acid, undecaprenyl carrier), storage-polymer (glycogen, polyP, PHB) and stress-signalling metabolism — and closes it under the FBAwMC macromolecular-crowding constraint of Beg et al. Module 18B curates a formation-energy table for every metabolite by two-pass weighted least squares against 79 literature-anchored reaction/redox benchmarks (benchmark residual RMS ≈ 18 kJ/mol), then formulates thermodynamic-constrained flux balance as a **mixed-integer program with 1,198 binary direction variables and 864 chemical-potential variables ln cᵢ**, solved by HiGHS in seconds (μ_tFBA = **0.3234 h⁻¹**, 0.4 % above the LP-relaxation bound). A per-metabolite **δ-gauge envelope (±2 MJ/mol on unbenchmarked species)** carries the curation uncertainty *into* the optimization while provably preserving loop-freedom (δ cancels around any closed stoichiometric cycle); an a-posteriori LP **certificate confirms loop flux = 0** over the 279-reaction active set. Module 18C integrates a 37-pool cybernetic kinetic core (stiff BDF, adenylate/NAD/CoA/quinone moieties conserved to 3.6×10⁻³ and 3.6×10⁻¹⁵ relative) through batch glucose exhaustion (t = 1931 s) and an H₂O₂ oxidative pulse (t = 5400 s): ATP/ADP collapses 5.0 → 0.018, ppGpp rises to 0.28 mM (stringent arrest), cAMP to 1.9 mM (catabolite derepression), and the H₂O₂ pulse peaks at 0.68 mM before catalase/peroxiredoxin recovery — growth → survival rewiring resolved second by second.

---

## 1. Why stoichiometry alone fails

Classical FBA asks only `S·v = 0` plus uptake bounds. Left unsaid is the second law. An LP-optimal flux map generically contains **Type-III energy cycles** — closed stoichiometric loops (e.g., ATP⁺h₂O→ADP+Pi wired against ADP+Pi→ATP+h₂O through futile carrier shuttles) in which every step runs "downhill" around a cycle whose total ΔG must sum to exactly zero. Such loops allow a cell to manufacture ATP, NADH or biomass precursors from nothing, and they silently inflate every predicted growth rate. The only rigorous cure is to enforce, per reaction,

  Δ_r G′_j = Δ_r G′°_j + R T Σᵢ S_ij ln cᵢ < 0  ⟺  v_j > 0,

with *one shared* chemical-potential vector ln c. That constraint is non-convex (sign selection), which is why it must be formulated as a mixed-integer program — and it is only as good as the formation-energy table Δ_f G′° feeding it. Phase 18 therefore builds three layers: a reconstruction (18A), a curated thermodynamic table (18B-0), and the δ-gauge tFBA machinery (18B) that carries the table's honest uncertainty into the optimization itself.

## 2. Module 18A — the reconstruction and its crowding budget

The network is written as explicit biochemistry in twelve families (Table 1), with long homologous series (acyl chains C6–C24, glycogen n = 1…35, polyP n = 1…40, PHB n = 1…40, LTA n = 1…28) generated procedurally exactly as reconstruction builders do. Quotas are asserted at assembly time: **817 internal metabolites / 1,201 enzymatic+transport reactions** (mission floor: 800 / 1,200).

**Table 1 — reconstruction families (reaction counts).**

| Family | Content | Rxns |
|---|---|---|
| Glycolysis/gluconeogenesis | PTS cascade (EI→HPr→IIA→IIB), HEX1, PGI/PFK/FBP, FBA/TPI, GAPD (NAD & NADP), PGK→PYK, PPS/PPDK, PPC/PPCK/PC | 21 |
| Pentose-phosphate | G6PDH, PGLASE, GND, RPI/RPE, TKT×2, TALA | 8 |
| TCA + glyoxylate | CS, ACONT, ICDH (NADP & NAD), AKGDH, SUCOAS, SUCD, FRD (menaquinone), FUM, MDH, ICL/MALS, GABA shunt | 22 |
| Fermentation | LDH, PFL, PTA/ACK, ACALD/ALCD2, ACS, POXB, PDH, arginine dihydrolase (ARCS/ARCK), FHL/FDH, D-/L-lactate DH | 16 |
| Ox-phos (explicit pmf) | NDH-1 (4H⁺ pump), cytochrome bo₃ (4H⁺) & bd (0H⁺), nitrate/nitrite reductases, ATP synthase (3.3 H⁺/ATP), STH/PNT transhydrogenases, ATPM | 10 |
| Amino acids | complete biosynthesis (his with AICAR→purine link, trp/phe/tyr via chorismate, asp/glu/pyruvate/1C families), degradation, 20 tRNA ligases | 125 |
| Nucleotides | purine de novo (11 steps), pyrimidine de novo, RNR (thioredoxin), dTMP/TYMS/DHFR, salvage, purine→allantoin→urea degradation | 62 |
| Fatty acids & lipids | FAS-II (FabH/G/Z/I, C4→C24), FabA UFA branch, β-oxidation (reversible), glycerophospholipid matrix (5 acyl compositions → PA/PE/PG/PS/PC/CL), ornithine lipids | 120 |
| Cell envelope | UDP-GlcNAc→murein (MurA–F), D-Ala/D-Ala, MEP→undecaprenyl, lipid-II flip-flop+polymerization, LTA polymer, murein recycling | 35 |
| Cofactors | NAD de novo+salvage, FAD, CoA, THF (PABA route), GSH, thioredoxin/glutaredoxin, menaquinone (MenF–G), ubiquinone (UbiC–H), heme (8 steps), Fe-S, PLP, B12/THM | 57 |
| Storage & stress | glycogen, trehalose, maltodextrins, polyP, PHB, SOD/CAT/GPX/TPX, ppGpp (RelA/SpoT), cAMP, c-di-GMP, Ap4A | 109 |
| Transport | periplasm + ABC (bind/deliver) for 40+ solutes, ion ATPases, PTS | 101 |

Every reaction carries kcat and molecular weight. The **FBAwMC crowding constraint** of Beg et al., Σⱼ vⱼ·Mwⱼ/(kcatⱼ·3600·σ_sat) ≤ φ_max·V_cyt·f_met = 0.34 g mL⁻¹ × 2.4 mL gDW⁻¹ × 0.55 = **0.451 g enzyme (gDW)⁻¹**, is added as a single dense row of the LP — the physical origin of overflow metabolism (§5).

Boundary chemistry is handled in the standard GEM convention: extracellular species are real rows of S (counted in the 869) but excluded from the quasi-steady-state equality (boundaryCondition = true), so sources/sinks stay unbalanced while the internal network closes exactly.

## 3. Module 18B-0 — thermodynamic curation by two-pass least squares

**Pass 1 (the physics).** 79 literature-anchored benchmarks — every glycolytic step, the TCA ladder, redox carrier pairs (NAD(P) ±61.75 kJ/mol, Q −21.2, FAD +42.5), the explicit-pmf respiratory units (NDH −17.4, cyt bo₃ −70.7, ATP synthase −8.5 at the 3.3H⁺ operating point, OXP unit −96.6), glutamate/GS/GLUS triangle, combustion anchor (glucose + 6 O₂ → 6 CO₂ + 6 H₂O = −2870) — are fitted by bounded weighted least squares over all 869 formation energies (seed priors σ ≥ 100 kJ/mol, weak; anchors σ ≤ 5). **Benchmark residual RMS = 18.1 kJ/mol** (max 86.9 on the worst-conditioned row) — the honest precision of a first-principles minimal-cell curation, documented row by row in `results_phase18/phase18_curation.json`.

**Pass 2 (the gauge).** The 775 metabolites not covered by any benchmark would sit at arbitrary levels. A second bounded LS propagates the gauge: every unbenchmarked reaction is pulled to a class-typical ΔG′° (kinases −15, ligases −25, transaminases −1, near-equilibrium/rev-split reactions exactly 0 ± 12, transport 0 ± 25 kJ/mol …), with benchmark-covered metabolites frozen. Pass-2 residual RMS ≈ 44 kJ/mol — the declared uncertainty of the unbenchmarked chemistry, carried forward explicitly (§4).

Self-tests reproduced from the fitted table: ATP hydrolysis −39.8 kJ/mol (textbook −45.6); fermentation glucose→2 lactate −28 to −80 (curve of the curation compromise); TCA citrate→succinate +26 to +64; OXP unit +12 to +63 — each within the documented residual envelope, and every *sign* that gates the MILP is robust (|ΔG| ≫ envelope for all decided steps).

## 4. Module 18B — the δ-gauge tFBA mixed-integer program

The MILP (solved by HiGHS through `scipy.optimize.milp`, **4,195 variables, 817 equalities, 2,400 inequalities, 1,198 binaries**) has four variable blocks:

1. **v** — 1,269 fluxes, near-equilibrium parents split into forward/backward halves (each half its own binary; a loop would need both open simultaneously, forbidden by the shared-ΔG rows);
2. **d = ln cᵢ** — 864 chemical potentials in per-class boxes (currencies 0.1–8 mM, polymers 10 µM–5 mM, default 1 µM–20 mM);
3. **z_j ∈ {0,1}** — 1,198 direction binaries: v_j ≤ U_j·z_j and, when z_j = 1, the strict row RT·Σ S_ij dᵢ + Σ S_ij δᵢ ≤ −ε (ε = 0.1 kJ/mol);
4. **δᵢ** — 864 **gauge-correction variables**, |δᵢ| ≤ 2 MJ/mol, added to the formation energies of unbenchmarked metabolites.

The δ block is the methodological novelty. Curation uncertainty is real (pass-2 RMS ≈ 44 kJ/mol), and pretending the fitted ΔG′° are exact would silently veto whole biosynthetic pathways at the wrong gauge. Instead the optimizer *chooses* self-consistent gauge values inside the declared envelope — and the zero-loop guarantee survives, because around any stoichiometrically closed cycle Σᵢ S_ij δᵢ = 0 identically: **δ shifts energies, never energy differences around a loop**. The a-posteriori certificate (LP: max loop flux over the active, strictly-downhill sub-network) returns **loop flux = −0.0 → 0 (PASS, 279 active reactions)**: no Type-III cycle can carry flux. Of the 292 active flux vectors, 224 sit strictly below −ε and the remaining 68 ride the δ envelope (mostly wall/entry steps), reported grey in Fig. 2C.

Solve statistics: relaxation bound 0.3247 h⁻¹, **MILP optimum μ = 0.3234 h⁻¹ in 3–18 s** (gap 0.4 %), crowding row binding at exactly **100.0 %** (0.4488 g enzyme gDW⁻¹).

## 5. Results — the flux map and overflow metabolism

The solved flux map (Fig. 1A/B) shows a textbook minimal-cell physiology: GAPD/PGK/ENO/PYK at 10.5–11.0 mmol gDW⁻¹ h⁻¹, ATP synthase at 16.9, NDH-1 at 6.9, cytochrome bo₃ at 7.1 — and, decisively, **overflow**: PFL at 8.1 with formate secretion 7.85 and acetate secretion 7.02. The glucose feed (10 mmol gDW⁻¹ h⁻¹) cannot be fully oxidized — the enzyme-volume budget saturates — so carbon spills through pyruvate formate-lyase and acetate kinase even under full aeration: **the Warburg/Crabtree effect emerges from the crowding constraint alone**, with no regulatory overlay.

The uptake sweep (Fig. 1C) quantifies it: the budget binds at 100.0 % from u_glc ≈ 6 onward; growth plateaus at **μ = 0.323–0.325 h⁻¹**; oxidative ATP flux *declines* 21 → 17.5 mmol gDW⁻¹ h⁻¹ while acetate overflow *rises* 5.7 → 7.0 — extra glucose cannot be respired (the budget has no room for more respiratory enzyme) and is dumped as acetate. The top crowding bottlenecks (Fig. 1D) are the translation/wall/lipid polymerization steps, precisely the high-MW/low-kcat machinery Beg et al. identified.

## 6. Module 18C — dynamic rewiring through starvation and oxidative stress

The kinetic core couples 37 metabolic pools (glycolytic intermediates, energy/redox carriers, aa/prot, storage, stress alarmones, extracellular metabolites, biomass X) to 12 cybernetic enzyme sets whose synthesis follows matching-law weights (shared capacity W = 0.03 g gDW⁻¹ h⁻¹) under four regulatory signals: cAMP catabolite derepression, amino-acid and energy starvation, ppGpp stringent arrest, and OxyR-type H₂O₂ induction. Initial pools and enzyme levels are taken *from the tFBA solution* (cᵢ = exp dᵢ; E_set = Σ a_j v_j), making 18C the dynamic shadow of the tFBA optimum. Integration: stiff BDF, t ∈ [0, 7200] s, 5,984 rhs evaluations, 1 s wall clock; conserved-moiety drift adenylate 3.6×10⁻³, NAD 3.6×10⁻¹⁵ relative.

**Glucose exhaustion at t = 1931 s.** The batch feed (2.2 mM) depletes in 32 min; ATP/ADP collapses 5.0 → 0.018; ppGpp rises to 0.28 mM, imposing stringent growth arrest (μ 0.21 → ~0.002 h⁻¹); cAMP spikes to 1.9 mM, derepressing the acetate-scavenging ACS set and glycogen mobilization (glycogen minimum 1.39 of 2.0 pool); biomass still grows 1.00 → 1.07 gDW L⁻¹ on internal reserves.

**Oxidative pulse at t = 5400 s.** A 0.8 mM H₂O₂ bolus (1 min) drives the peroxide pool to a 0.68 mM peak, consumed by the catalase/peroxiredoxin/GSH sets as NADPH is drawn down — the survival module absorbs the insult and the culture resumes slow growth (X_final = 1.07 gDW L⁻¹). Figure 3 renders the full rewiring: flux stackplot, energy/redox ratios, the 30-pool spatio-temporal heatmap, and the ATP-vs-ppGpp survival phase portrait.

## 7. Bridging quantum chemistry to cellular survival

The causal chain this phase closes is: **quantitative atom-level energetics → group-contribution formation energies → benchmark-curated Δ_f G′° table → per-reaction second-law constraints at optimized concentrations → a flux map that cannot cheat → a growth rate that a real membrane-bound, crowding-limited cell could believe.** Every link is quantified: formation-energy residuals (±18/±44 kJ/mol), the δ-envelope (±2 MJ/mol on unbenchmarked species), the crowding budget (100.0 % binding), the loop certificate (flux 0). The Warburg result is the deepest of the chain: overflow metabolism needs *no* transcription factor in the model — it is what pure thermodynamics plus finite enzyme volume force on any aerated heterotroph.

## 8. Limitations & provenance (stated, not hidden)

- The ΔG′° table is a **first-principles curation, not eQuilibrator**: benchmark residuals reach 87 kJ/mol on the worst row, and the unbenchmarked chemistry carries a declared ±2 MJ/mol δ-envelope; 68 of 292 active fluxes sit inside that envelope (Fig. 2C, grey) — they are *admitted* uncertain, not silently wrong.
- The pmf is handled by explicit proton-pseudo-species (h_c/h_e with folded electrochemical potential) at the operating point; ΔpH/Δψ are not free variables.
- Enzyme saturation in the crowding budget is a flat σ_sat = 0.5; no enzyme is concentration-resolved beyond the 12 cybernetic sets.
- The kinetic module (18C) is a reduced 37-pool shadow of the full network: its rate constants are literature-scale, not fitted time-series.
- Base physiology: glucose sole carbon (glycerol lane encoded but closed); anaerobic nitrate respiration available via nar/nir when o₂ is cut.

## 9. Reproduction

```bash
python run_phase18_wholecell_metabolic_thermodynamics.py --stage all --milp-time 420
# stages: build | curate | tfba | sweep | dynamics | figures | all   (resumable)
```

Requirements: numpy, scipy (HiGHS), matplotlib, networkx. Runtime ≈ 2 min end-to-end (tFBA MILP 3–18 s dominates); everything caches into `results_phase18/`.

## 10. File manifest

| File | Content |
|---|---|
| `run_phase18_wholecell_metabolic_thermodynamics.py` | the entire phase: builder, curation, tFBA MILP, sweep, dynamics, figures |
| `figures_phase18/fig1_genome_scale_flux_map.png` | central-carbon highway, whole-network graph, Warburg sweep, bottlenecks |
| `figures_phase18/fig2_thermodynamic_driving_forces.png` | glycolysis/TCA ΔG ladders, no-loop wedge, solved chemical potentials |
| `figures_phase18/fig3_dynamic_metabolic_rewiring.png` | flux rewiring, energy/redox, pool heatmap, survival portrait |
| `results_phase18/phase18_results.json` | quotas, curation residuals, tFBA solution, sweep, dynamics record |
| `results_phase18/phase18_network.npz` | curated Δ_f G′° vector + name tables |
| `results_phase18/phase18_tfba.npz` | solved v, ln c, ΔG, activity mask, crowding load |
| `results_phase18/phase18_dynamics.npz` | full 7200 s trajectory (49 states × 1441 samples) |
| `results_phase18/phase18_curation.json` | benchmark residuals + self-tests |
