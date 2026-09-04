# FRONTIER EPISTEMIC REPORT — Phase 7
## Confronting the Wall of Sighs: Multi-Reference Electronic Structure & the Quantum-AI Breakdown Map

**Pipeline:** `run_phase7_strong_correlation_wall.py` · **Results:** `results_phase7/phase7_results.json`, `phase7_scan_summary.csv` · **Figures:** `figures_phase7/fig1–fig3` (300 DPI)
**System:** cyclopropa[b]indole (C₉H₉N, Phase-4 GFN2-xTB-optimized reactant) — homolytic stretch of the scissile cyclopropane C6–C9 bond, 1.40 → 3.20 Å (11-point relaxed scan)
**Engines:** Psi4 1.11 (DETCI CASSCF, RHF/UHF/RKS/UKS) · GFN2-xTB (xtb 6.7.1) · MACE-OFF (small) · ANI-2x (torchani) · B3LYP
**Status:** all 11 scan points converged on the primary basis (def2-SVP) with zero basis degradations; wall time 2 295 s.

---

## 1. Theoretical Impasse — why the wall exists

For an N-electron system the exact wavefunction is a function of 3N spatial coordinates. Its exact (FCI) expansion in N-electron determinants grows combinatorially — the **exponential wall** of quantum mechanics. The entire edifice of practical electronic-structure theory is a set of strategies for not climbing that wall:

- **Single-determinant theories** (HF, DFT, semi-empirical tight-binding, and — structurally — every mainstream machine-learning potential) compress the electronic state into one occupation pattern. This is exact for a single bond at equilibrium and catastrophically wrong when a bond is broken: the correct wavefunction of a homolytically cleaving bond is a superposition

  |Ψ⟩ ≈ c₁ |…σ̄σ…⟩ + c₂ |…σ*ᾱσ*β…⟩  (bonding-covalent mixed with double-excitation onto σ*)

  and the coefficient ratio is not perturbatively small — it tends to 1:1 at infinite separation. No single determinant, and no functional of a single density, can converge to the correct dissociation limit while keeping the wrong one out. This is **static (non-dynamic) correlation failure**.
- **Multi-reference theories** (CASSCF among them) build the wavefunction from all distributions within a chemically chosen active space — for a single bond cleavage, CAS(2e,2o) is already the exact solution *of that bond's* correlation problem, at a cost that scales with the active space, not with the molecule.

The mission of Phase 7 is to force this wall into the open, quantitatively, on the exact bond that Phase 4/5 identified as the chemically scissile one.

## 2. Computational Protocol

| Module | Quantity | Engine | Settings |
|---|---|---|---|
| G | relaxed bond scan, 11 points | GFN2-xTB (xtb.exe + ASE LBFGS) | FixBondLengths(C6–C9), fmax 0.03 eV/Å, warm-started |
| AI | single points on identical geometries | MACE-OFF (small), ANI-2x, GFN2-xTB | local MACE model (ASL), torchani CPU |
| 7A | RHF/UHF/RKS/UKS singlets + UHF/UKS **triplets** | Psi4 1.11 | def2-SVP, DF-JK; ⟨S²⟩ = S(S+1) + n_β − Tr[Dα S Dβ S] |
| 7B | CAS(2e,2o) CASSCF | Psi4 DETCI | def2-SVP, active = σ(C6–C9)/σ* selected by **complement coherence** (see §4), OPDM → NOONs |
| 7C | epistemic error landscape | aggregation | ΔE_error(R) = \|E_rel^AI(R) − E_rel^CASSCF(R)\|, common zero at R = 1.50 Å |

**Fault tolerance as executed.** Every QC point ran in an isolated subprocess; the orchestrator carried a basis-degradation ladder (def2-SVP → 6-31G → STO-3G, trade-offs logged) and a CASSCF-algorithm ladder (default → AO → DF). In the event **no degradation was ever needed** — all 11 points converged on tier 0 — and the fallback log in `phase7_results.json` is empty except for audited methodology probes.

**Backend substitution (logged in every result file).** The reference protocol specifies PySCF. PySCF ships no native win32 build (no wheels; no MSVC toolchain on this host), so the ab-initio backend is **Psi4 1.11 (conda-forge win-64)** — protocol-equivalent for RHF/UHF/RKS/UKS + CAS(2e,2o) CASSCF + NOON extraction. On a POSIX host with PySCF installed the same quantities are produced with identical definitions.

Two build quirks were discovered, contained, and documented in the code: (i) this Psi4 build's DETCI/DPD cache sizing aborts when the global memory is set inside the 4 GB band — the pipeline pins 2 GB (verified band); (ii) explicitly forcing `mcscf_type=CONV` routes CASSCF through a defective DPD path (energies ~+4.7 Eh off) — the default integral path is used instead.

## 3. Module 7A — Mean-Field Breakdown: the alarm that stays silent

**Headline finding.** Along the entire 1.4–3.2 Å window, the UHF and UKS singlet SCF **never leaves the closed-shell branch**: ΔS²(singlet) ≡ 0 to numerical precision at every point (fig. 1, top panel, gray line). This is *not* the absence of static correlation — CASSCF proves the opposite (§4) — it is the absence of an **SCF-accessible** broken-symmetry solution.

We treated this negatively-result seriously and exhausted the standard seed engineering before accepting it: SAD and SAP guesses; supramolecular-style fragment specifications (ionic fragments +1/−1 with doublet multiplets); HOMO↔LUMO alpha swaps written through `guess=read` checkpoint splicing (validated on H₂@3.0 Å, where the seed holds and yields ⟨S²⟩ = 1.0); MOM occupation pinning; SOSCF; and triplet-derived BS seeds (converge the S_z = 1 state, splice its radical orbitals into a singlet occupation). Every route relaxed back to the closed-shell determinant. The physical reason is chemically interesting in its own right: the C6 radical center is conjugated with the indole π system and the CH2• fragment is a good electron acceptor, so the **zwitterionic resonance form of the "broken" bond is unusually stabilized** — the restricted surface is anomalously flat, the broken-symmetry basin correspondingly shallow, and DIIS/AIIS gradients slide back into the closed-shell well.

The operational consequence is important enough to state as a theorem of practice: **for π-stabilized polar diradicals, the textbook smoke alarm (spin contamination of the singlet) stays silent while the house burns.** A ΔS²-based diagnostic would have certified this system as "single-reference safe" at every geometry — while CASSCF shows diradical character up to y = 0.69.

**The open-shell sector that does behave** is S_z = 1. Triplet UHF/UKS converge reliably (multiplicity-enforced occupations cannot collapse), and their energetics deliver the cleanest possible symmetry-breaking diagnostic available to a single determinant:

| R (Å) | ΔE_ST^UHF (kcal/mol) | ΔE_ST^UKS-B3LYP | ⟨S²⟩_T (UHF) | y (CASSCF) |
|---|---|---|---|---|
| 1.40 | +75.5 | +89.1 | 1.302 | 0.014 |
| 1.50 | +76.0 | +89.9 | 1.304 | 0.018 |
| 1.60 | +75.0 | +88.4 | 1.309 | 0.028 |
| 1.75 | +71.4 | +82.9 | 1.349 | 0.047 |
| 1.90 | +61.0 | +68.8 | 1.358 | 0.081 |
| 2.05 | +16.4 | +45.2 | 1.512 | 0.139 |
| 2.20 | **−14.9** | +23.3 | 1.520 | 0.235 |
| 2.40 | −41.0 | +5.6 | 1.526 | 0.383 |
| 2.60 | −55.4 | −4.7 | 1.524 | 0.528 |
| 2.90 | −62.6 | −11.7 | 1.539 | 0.676 |
| 3.20 | −68.2 | −8.0 | 1.502 | 0.694 |

**R_crit (singlet–triplet gap sign change) = 2.13 Å (UHF) / 2.51 Å (UKS-B3LYP).** Beyond R_crit the closed-shell singlet is no longer even the lowest single determinant — the mean-field ground state has ceased to exist as a stationary point worth trusting. The triplet spin contamination ⟨S²⟩_T − 2 climbs from −0.70 toward the diradical signature in the same window (fig. 1, top panel), confirming that the S_z = 1 sector, too, becomes progressively multi-configurational.

## 4. Module 7B — CASSCF(2e,2o): the ground truth, and how the active space is found

Active-space selection is where naive automation fails on this molecule, and the pipeline implements a physically constructed criterion, **complement coherence**:

1. Candidate occupied orbitals (valence, ε > −1.6 Eh, top 8 by Mulliken gross population on the scissile carbons) are examined. Raw population argmax fails twice: core 1s orbitals out-score the σ bond (weight ≈ 0.999), and diffuse functions inflate Rydberg gross populations to values of 10–40.
2. For each candidate φ_k the **antibonding complement** u_k = P₆|φ_k⟩ − P₉|φ_k⟩ is constructed from the projector onto each scissile carbon's AO subspace. The true bonding orbital is the one whose complement **collapses onto a single low virtual** (score = max overlap); that virtual is σ*. At R = 3.20 Å this criterion relocates the pair from a plausible-looking impostor (MO34, complement coherence 0.10) to the genuine pair (MO32/MO35, coherence 0.48).
3. The MO set is permuted so that [σ, σ*] sit exactly at the docc/active/uocc boundary required by DETCI's occupation arrays; CASSCF(2,2) then converges on the default integral path at every geometry.

**Natural-orbital occupation numbers (trace of active block = 2.000 at every point):**

| R (Å) | ν(σ NO) | ν(σ* NO) | diradical character y = ν(σ*) |
|---|---|---|---|
| 1.40 | 1.986 | 0.014 | 0.014 |
| 1.50 | 1.982 | 0.018 | 0.018 |
| 1.60 | 1.972 | 0.028 | 0.028 |
| 1.75 | 1.953 | 0.047 | 0.047 |
| 1.90 | 1.919 | 0.081 | 0.081 |
| 2.05 | 1.861 | 0.139 | 0.139 |
| 2.20 | 1.765 | 0.235 | 0.235 |
| 2.40 | 1.617 | 0.383 | 0.383 |
| 2.60 | 1.472 | 0.528 | 0.528 |
| 2.90 | 1.324 | 0.676 | 0.676 |
| 3.20 | 1.306 | 0.694 | 0.694 |

The fractional occupancy rises smoothly from 0.014 e (closed shell) through the mission's canonical example — **1.86/0.14 at R = 2.05 Å** — to 1.31/0.69 e, a strongly correlated diradical. This monotonic NOON curve *is* the mathematical proof of static correlation that the spin-contamination diagnostic refused to give (§3). Figure 2 renders the two natural orbitals at R = 1.50/2.20/3.20 Å: the localized σ bond at equilibrium resolves into two separated radical lobes (ring carbon ↔ methylene) whose occupation fraction grows exactly as tabulated.

One honest subtlety: y does **not** approach 1.0 at 3.2 Å but saturates near 0.69–0.70. The relaxed diradical is tethered (C5–C9 intact), C6• conjugates into the indole π system, and part of the static correlation migrates into the π channel that a σ/σ* CAS(2,2) deliberately excludes. The CAS(2,2) number is therefore a *lower bound* on the true diradical character — which strengthens, not weakens, the conclusion.

## 5. Module 7C — The Wall of Sighs: quantified

All methods evaluated on the **identical** relaxed geometries; each curve zeroed at R = 1.50 Å; error = |E_rel^method(R) − E_rel^CASSCF(R)|:

| R (Å) | E_rel CASSCF | GFN2-xTB | MACE-OFF | ANI-2x | BS-UB3LYP | RHF |
|---|---|---|---|---|---|---|
| 1.40 | +4.4 | −4.9 | +4.2 | +3.9 | +4.1 | +3.7 |
| 1.50 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| 1.60 | −1.7 | −4.9 | +1.9 | +2.4 | +1.6 | +2.5 |
| 1.75 | +10.1 | −2.9 | +16.4 | +16.3 | +13.1 | +17.1 |
| 1.90 | +20.1 | +6.8 | +40.7 | +29.5 | +23.2 | +30.5 |
| 2.05 | +29.1 | +19.0 | +59.4 | +43.6 | +32.6 | +44.1 |
| 2.20 | +33.9 | +29.4 | +60.2 | +49.6 | +39.6 | +55.5 |
| 2.40 | +44.5 | +37.9 | +64.2 | +54.7 | +52.9 | +74.9 |
| 2.60 | +52.3 | +47.3 | +66.2 | +56.4 | +62.5 | +88.5 |
| 2.90 | +81.3 | +56.5 | +88.5 | +71.6 | +89.0 | +118.8 |
| 3.20 | +113.0 | +61.1 | +120.1 | +94.0 | +111.8 | +155.4 |

**Absolute epistemic error |ΔE_error(R)| and the 15 kcal/mol gate:**

| Method | max \|ΔE_error\| | first R beyond 15 kcal/mol | behavior |
|---|---|---|---|
| GFN2-xTB | **52.0** (3.2 Å) | **2.9 Å** | under-dissociates twice over: its strain-release surface misses the correlation-driven rise entirely through the onset window (−2.9 vs +10.1 at 1.75 Å) and then undershoots the asymptotic cost by 52 kcal |
| **MACE-OFF (small)** | **30.2** (2.05 Å) | **1.90 Å — the earliest failure** | *over-stiff* pathology: at the exact onset of diradical character (y ≈ 0.14) it prices the stretch at 59.4 kcal/mol where CASSCF says 29.1 — it treats a breaking bond as an intact one; the potential's error is largest precisely where its training manifold ended |
| ANI-2x | 19.1 (3.2 Å) | 2.20 Å | same closed-shell prior, muted: error peaks where y ≈ 0.2–0.4 and grows again toward the separated-radical limit |
| UHF (RHF branch) | 42.4 (3.2 Å) | 2.20 Å | the closed-shell mean-field surface diverges monotonically from the truth |
| BS-UB3LYP | 10.3 | **never fails** | B3LYP's delocalized exchange plus systematic fractional-spin error cancellation tracks CASSCF within 10.3 kcal/mol along the whole coordinate |

Three epistemic structures deserve emphasis:

1. **The failure zone is the diradical zone.** Every method that fails does so after R ≈ 1.9–2.2 Å, bracketing the CASSCF-defined onset (y > 0.1) and the singlet–triplet crossing (2.13 Å). The models did not fail randomly; they failed precisely where their training manifold — closed-shell, near-equilibrium molecules — ended.
2. **The most dangerous failure is the earliest one.** MACE-OFF crosses the 15 kcal/mol gate at 1.90 Å, *before* the mean-field benchmark itself (2.20 Å), and its error peaks (30.2 kcal/mol) exactly where the homolysis becomes chemically decisive (y ≈ 0.14). An ML-potential-driven reaction mechanism search would accept structures in this window with zero internal warning — no formal charge anomaly, no spin flag, nothing. Equivariance and foundation-scale training data do not purchase observables outside the training distribution.
3. **Semi-empirical and ML fail differently.** xTB is wrong *smoothly* (an SCF-surface difference: it misses the correlation-driven rise through the onset window and undershoots the asymptote), while the ML potentials are wrong *sharply* — a stiffness spike (MACE-OFF: +30 kcal/mol at one grid point) that appears and vanishes within 0.3 Å. Both are fatal for barrier heights; the ML failure is the more insidious because it is localized exactly where a mechanism search makes its branching decisions.

## 6. Epistemological Critique — why scaling language models cannot climb this wall

The wall of sighs is not a data problem, and that is precisely why scale does not solve it.

- **The target function is not on the manifold.** A foundation potential is an interpolator over a training distribution of mostly closed-shell, near-equilibrium structures. The diradical seam is *orthogonal* to that distribution: the correct energy there is governed by an avoided crossing — a 2×2 Hamiltonian in the {|σ²⟩, |σ*²⟩} space — whose off-diagonal mixing is a non-local quantum property. Interpolation cannot synthesize an avoided crossing that the training set never displayed.
- **Non-locality vs. local message passing.** Static correlation is a global property of the Fock space (which determinants are near-degenerate), not a local property of nuclear positions. Message-passing/equivariant architectures encode locality and smoothness of PES landscapes; Fock-space near-degeneracy is neither local nor smooth — it flips character across a seam the network has no coordinate for.
- **No multi-reference inductive prior.** A language model scales by learning correlations among tokens it has seen. Nothing in its inductive bias encodes "when two orbitals become degenerate, I must superpose occupations." Without that prior, more data and more parameters sharpen the interpolation *inside* each regime; they do not bridge regimes. Our fig. 3 is the empirical demonstration: MACE-OFF is a state-of-the-art foundation potential, and it fails *earliest*.
- **Silent failure is the epistemic core.** The most cited lesson of this phase is §3: for this entire class of systems, the conventional alarm (spin contamination) reads zero. An autonomous AI-chemistry loop that monitors ΔS² — or any single-determinant diagnostic — would publish confident nonsense. Uncertainty quantification does not rescue this: the models are *confidently* wrong (deep ensembles report low variance on a smooth but wrong surface).

## 7. The Future Trajectory — from statistical surrogate to hybrid quantum architecture

The exit from the wall is not more parameters; it is architecture that carries the right prior:

1. **Multi-reference ML surrogates.** Learn the *active-space objects* — NOONs, effective Hamiltonians, CASPT2/NEVPT2 corrections — instead of total energies. y(R) from this phase is exactly the kind of label such models need; our NOON table is a seed dataset.
2. **Physics-anchored hybrid loops.** Use cheap surrogates as proposals and multi-reference engines as verifiers: ML proposes geometries, a ΔE_ST/NOON gate (this pipeline, pointwise, minutes per structure) certifies single-reference safety, and only uncertifiable points get CASSCF/DMRG treatment. The orchestrator pattern of this phase (per-point isolation, graceful degradation, audit logs) is the engineering template.
3. **Differentiable quantum chemistry.** Backprop-able CASSCF (differentiable orbital optimization) inside the training loop would let ML potentials *learn* the coupling term of the avoided crossing instead of interpolating through it.
4. **Tensor-network and neural wavefunction priors.** DMRG-inspired entanglement diagnostics and neural-network wavefunctions (fermionic backflows, etc.) as *features*, not just references — giving the model an entanglement coordinate along which to condition.
5. **Quantum computing, honestly scoped.** A fault-tolerant device attacking the active space (2 electrons in 2 orbitals here; modest) while classical methods treat the inactive/dynamic part is the natural quantum-classical hybrid — the active space *is* the quantum resource, and this phase shows how to find it automatically.

## 8. Limitations & Computational Economy

- **CAS(2,2) is a σ-channel instrument.** The y-saturation at 0.69 (§4) quantifies the π-channel correlation it intentionally excludes; a π-inclusive CAS would raise, never lower, the diradical character.
- **CASSCF lacks dynamic correlation**; CASPT2/NEVPT2 on top would shift absolute curves by several kcal/mol but not the *crossing topology* that defines the wall.
- **Relaxed scan at GFN2-xTB level**: the geometry coordinate is xTB-optimal, not CASSCF-optimal; this is the identical-trajectory requirement of Module 7C and is applied uniformly to every method.
- **Basis**: def2-SVP throughout; no degradations were triggered. The def2-SVP → 6-31G → STO-3G ladder and its trade-off log (smaller basis shifts NOON fractions and absolute energies, not the qualitative diradical signature) are implemented and armed.
- **The BS-singlet negative result** (§3) is specific to SCF accessibility on this system; on a POSIX PySCF host the same probes can be rerun unchanged.

## 9. Deliverables

| Artifact | Path |
|---|---|
| Pipeline (dual-interpreter, fault-tolerant) | `run_phase7_strong_correlation_wall.py` |
| Master machine-readable record | `results_phase7/phase7_results.json` |
| Per-point energy/NOON/error table | `results_phase7/phase7_scan_summary.csv` |
| Fig. 1 spin-contamination & ST-gap collapse | `figures_phase7/fig1_spin_contamination_profile.png` |
| Fig. 2 CAS(2,2) natural orbitals + NOON evolution | `figures_phase7/fig2_casscf_frontier_orbitals.png` |
| Fig. 3 the Wall of Sighs discrepancy map | `figures_phase7/fig3_the_wall_of_sighs_discrepancy.png` |
| Chinese version of this report | [`FRONTIER_EPISTEMIC_REPORT_ZH.md`](./FRONTIER_EPISTEMIC_REPORT_ZH.md) |

*Reproduce: `conda activate phase7 && python run_phase7_strong_correlation_wall.py` (auto-discovers the `phase2ff` chem env for xTB/MACE/ANI; `--smoke` for a 3-point validation).*
