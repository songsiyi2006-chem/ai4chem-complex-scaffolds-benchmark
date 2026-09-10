# Phase 25 — Local calculation progress and evidence review

**Snapshot: 2026-09-10T16:10:30.622155+00:00. Scientific acceptance remains NOT PASSED; product selectivity remains undetermined.** This update extends the initial `93b9ab1` release with actual optimizations, Hessians, complete ternary structures and a running local DFT queue. Phase26/27 have no new production calculations in this update.

## Completed and running work

- **33 xTB pilot jobs** are recorded, including controls and retained failed minimum candidates; **25** pass the preliminary local-minimum audit. A successful engine exit alone is insufficient.
- Eight complete **117-atom** C01/imine/dimethyl-Hantzsch-ester starts cover E/Z, both approach-face signs and two rigid-placement starts per combination. **7** complete ternary candidates currently pass the xTB-level checks. The catalyst is not truncated. Face signs are construction labels, never assigned R/S product channels.
- The isolated I01-E PBE-D3BJ/def2-SVP optimization converged in **19.76 min**, with energy **-594.909268679091 Eh**, maximum Cartesian gradient **4.534e-06 Eh/bohr**, and RMS gradient **1.407e-06 Eh/bohr**. E geometry is retained. A converged optimization is not a frequency-validated minimum.
- The same-surface native DFT Hessian has **41/163 completed atomic gradient jobs** in the frozen public checkpoint snapshot. **No complete DFT frequencies or DFT Gibbs energy are claimed.** Live output is excluded from the release manifest until separately reviewed.
- The prior four isolated-imine DFT single points remain historical evidence. The background queue adopts the current Hessian process, then runs Z optimization/frequencies and omegaB97X-D single points on the PBE-optimized E/Z geometries. Queued calculations are not counted as completed.

## DFT job ledger at the snapshot

| Job | Method | Stage | Status | Energy (Eh) | Final E/Z geometry |
|---|---|---|---|---:|---|
|pbe_E_freq_01|pbe-d3bj|freq|FAILED|pending|not assessed|
|pbe_E_freq_02|pbe-d3bj|freq|STARTED|pending|not assessed|
|pbe_E_opt_01|pbe-d3bj|opt|OPTIMIZED_FREQUENCY_PENDING|-594.909268679091|E|
|pbe_Z_opt_01|pbe-d3bj|opt|OPTIMIZED_FREQUENCY_PENDING|-594.905102976489|Z|

Both E and Z optimizations have now converged. Z retains Z geometry; its maximum gradient is 3.661e-06 Eh/bohr and wall time 26.54 min. The gas-phase PBE-D3BJ electronic E(Z)-E(E) difference is 2.614 kcal/mol. This is not a Gibbs-energy difference or catalytic barrier. It must not be directly attributed to a functional effect against the solvated xTB table, since both method and environment differ.

Z optimization was additionally started when free local memory increased; the finite queue adopts or skips that same job to avoid duplicate work. A STARTED/FAILED frequency row has no accepted frequency result.

## Methods and acceptance rules

New tight-binding jobs use xTB 6.7.1, GFN2-xTB, ALPB(toluene) or ALPB(CH2Cl2), charge 0, zero unpaired electrons, accuracy 0.2, and `vtight` optimization. Imaginary-mode remediation uses `extreme`. Optimization and Hessian use the same method and solvent. Runs are isolated, commands and raw output retained, and the Hessian follows a confirmed optimization using `--hess --strict`. The [xTB manual](https://xtb-docs.readthedocs.io/en/latest/hessian.html) explains the Hessian and imaginary-mode handling.

The independent audit requires the complete 3N frequency block, agreement between repeated printouts, six projected translation/rotation modes, and all internal modes above 1 cm^-1. This last threshold is a numerical screen, not a demonstrated frequency error bound. The audit also checks atom ordering, isolated-imine E/Z identity, and, for complete ternaries, catalyst-axis sign and expected heavy-bond distances. Soft modes and all rejected candidates remain visible. The screen is not a substitute for a chemically reviewed DFT stationary point or IRC.

The xTB thermochemistry summary can report zero imaginary frequencies after its default -20 cm^-1 cutoff even when a printed mode is negative. The initial I02-E structures had negative modes. Both signs of the engine-generated displacement were reoptimized; successful replacement candidates were retained while the others stayed rejected. Negative-mode records: **imine_dcm_repair_01/I02_E_dcm_plus: -31.66 cm^-1, imine_repair_01/I02_E_plus: -10.27 cm^-1, imines_dcm_01/I02_E: -23.21 cm^-1, imines_toluene_01/I02_E: -18.26 cm^-1, ternary_toluene_01/C01_I01_Z_face+1_start0: -4.88 cm^-1**. No absolute-value conversion hides these modes. Any negative-mode ternary candidate listed here remains unresolved and is excluded from the accepted-candidate table; the imine remediation does not certify that ternary.

## Solvent and low-frequency diagnostics

The table reports **single-local-minimum E/Z differences at the xTB level**, in kcal/mol. It is not a conformer-ensemble result, activation barrier, population prediction or reaction ee. The best audited candidate by G50 is selected per E/Z species from this limited search. Standard-state constants cancel in each equal-composition 1-to-1 comparison; absolute 1 M solution free energies and binding free energies are not certified.

| Imine | Solvent | G(Z)-G(E), cutoff 50 | G(Z)-G(E), cutoff 100 |
|---|---|---:|---:|
|I01|ALPB(toluene)|0.135|0.198|
|I02|ALPB(toluene)|0.386|0.475|
|I03|ALPB(toluene)|-0.761|-0.857|
|I04|ALPB(toluene)|0.493|0.441|
|I01|ALPB(ch2cl2)|0.285|0.355|
|I02|ALPB(ch2cl2)|0.386|0.413|
|I03|ALPB(ch2cl2)|-0.466|-0.575|
|I04|ALPB(ch2cl2)|0.428|0.410|

![E/Z cutoff sensitivity](results_phase25_27/phase25/figures/imine_EZ_sensitivity.png)

The 50/100 cm^-1 comparison reuses the actual saved Hessian in **native `xtb thermo`**, retaining its geometry-dependent inertia. All **25/25** reference replays agree with the original G correction within 1e-6 Eh. These are method-sensitivity contrasts, not confidence intervals. No fixed-inertia approximation is used in this reported table. The defaults are documented in [xcontrol](https://github.com/grimme-lab/xtb/blob/v6.7.1/man/xcontrol.7.adoc); the [thermodynamic implementation](https://github.com/grimme-lab/xtb/blob/v6.7.1/src/thermo.f90) uses the molecular rotational inertia. ALPB's default `gsolv` reference is retained and identified using the [reference-state documentation](https://xtb-docs.readthedocs.io/en/latest/gbsa.html#reference-states).

## Complete-system structures and controls

| Start | Final E/Z | Proton contact state | Relative G50 (kcal/mol) | Lowest mode (cm^-1) |
|---|---|---|---:|---:|
|E_face-1_start0|E|N_bound_candidate|0.000|12.14|
|E_face-1_start1|E|N_bound_candidate|2.189|6.64|
|E_face+1_start0|E|O_bound_candidate|6.231|6.27|
|E_face+1_start1|E|N_bound_candidate|4.287|5.41|
|Z_face-1_start0|Z|O_bound_candidate|7.083|7.65|
|Z_face-1_start1|Z|N_bound_candidate|4.035|8.78|
|Z_face+1_start1|Z|O_bound_candidate|19.963|7.46|

Relative G50 refers to the lowest audited candidate in this finite set. **These are minima, not competing transition states; converting this table to ee would be invalid.** Fixed-mapping heavy-atom RMSDs are provided, but search counts are not statistical degeneracies and no ensemble weights are assigned. Absolute catalyst axial CIP and product CIP remain pending.

![Complete optimized ternary and contacts](results_phase25_27/phase25/figures/ternary_minimum.png)

For E_face-1_start0, O...H is 1.629 A and N-H is 1.048 A; the designated donor C-H remains 1.094 A. At this approximate level, the geometry supports an iminium/phosphate ion-pair candidate without completed hydride transfer. Other starts include O-bound neutral hydrogen-bond candidates. These observations motivate separate protonation-state paths; they do not establish the experimental mechanism. The exact reflected full complex has a same-geometry xTB energy difference of **0.000e+00 kcal/mol at printed precision**. This parity check is not a kinetic ee-inversion test. Achiral and uncatalyzed rate controls remain uncomputed.

## Resources, failure recovery and unfinished criteria

This work uses the existing local CPU environment. Psi4 requests 500 MB and 2 threads for the completed optimization, then 3 threads for the Hessian/queue; xTB uses one thread. These are requested resources, not peak-memory guarantees. Recorded finished xTB jobs consumed **4318.9 allocated core-seconds** (wall time times assigned threads), including unsuccessful minimum candidates; active jobs and unrelated campaigns are excluded. This is not a cost-matched active-learning comparison.

The first DFT frequency attempt failed while decoding UTF-8 native output with Windows GBK, after its first gradient had completed. The failed output is retained. `python -X utf8` repairs that interface failure without modifying the installed engine. Exact QCSchema inputs and successful AtomicResults are hash-bound and checkpointed. A real small native-gradient test verifies fresh execution, identical cache replay and rejection of corrupted bindings. The published snapshot contains only completed, checked checkpoints. The local finite queue continues while the computer/processes remain running; it does not purchase compute, use OpenAI inference GPUs, publish automatically, or guarantee eventual convergence.

There are still no accepted DFT R/S catalytic transition states or IRCs, no correlated-wavefunction validation of selectivity-determining structures, no explicit-solvent free-energy replicas, no verified Curtin-Hammett separation, and no external blind outcome labels. Therefore no rates, conversion, signed ee, confidence coverage or 0.5 kcal/mol ensemble-convergence acceptance is reported. Phase26 constant-potential sampling and Phase27 multireference dynamics remain pending.

## Reproduction and evidence

[Progress data](results_phase25_27/phase25/pilot_progress.json), [raw xTB jobs](results_phase25_27/phase25/xtb_pilot), [independent audit](results_phase25_27/phase25/xtb_pilot/validation.json), [ternary diagnostics](results_phase25_27/phase25/ternary_diagnostics.json), [native thermochemistry](results_phase25_27/phase25/thermo_sensitivity/results.json), [DFT optimization](results_phase25_27/phase25/stationary_pilot/pbe_E_opt_01/result.json), [frozen native gradients](results_phase25_27/phase25/hessian_checkpoint_snapshot/manifest.json), [execution snapshot](results_phase25_27/phase25/local_execution_snapshot.json), [publication omissions](results_phase25_27/publication_omissions.json), [commands and limits](phase25_27/README.md).
