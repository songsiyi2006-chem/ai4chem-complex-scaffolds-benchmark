# Phase 21: cavity QED polaritonic chemistry

## Scope

An assigned harmonic Pauli-Fierz model demonstrates collective optical splitting,
dark modes and a local reaction bottleneck. The 2.012 eV symmetric quartic is a
schematic scale resembling the earlier skeletal-reaction barrier, not a recovered
Phase 4/5 PES or dipole function. No electronic-structure or experimental data are
fitted. A vacuum equilibrium Hamiltonian needs no pump; measuring FTIR requires
a weak probe. Those statements refer to different operations.

## 21A: dipole gauge, basis and units

In energy units (hbar=1), H=sum E0(b_i†b_i+1/2)+Ec(a†a+1/2)
+g0(a+a†)sum(b_i+b_i†)+(g0²/Ec)[sum(b_i+b_i†)]².
This includes counterrotating and all cross-molecule dipole self-energy terms.
The projected squared coordinate is evaluated analytically; squaring a truncated
coordinate matrix would incorrectly remove the top-level boundary contribution.
Occupations are photons 0..3 and vibrations 0..4. A full N=2 tensor-product
calculation verifies the bright/dark reduction; N-fold full tensor products are
not materialized for large N. Identical harmonic molecules reduce exactly to one
bright oscillator with g=g0 sqrt(N), plus N-1 uncoupled dark oscillators.
The larger 12-photon-state, 14-vibration-state bright basis checks convergence.
The reference normal modes use the full quadratic PF Hessian, not the rotating-wave
approximation. Rabi scaling is checked at resonance and fixed single-molecule g0.
The input coupling is assigned; no actual cavity volume, transition dipole, or
absolute SI vacuum field is inferred. g0=lambda mu01 sqrt(Ec/2) explains its
meaning but does not furnish a material calibration.

## 21B: what the model says about reactions

V(q)=B[1-(q/q0)²]², with effective mass chosen to give E0 at the well;
the imaginary barrier energy is E0/sqrt(2), consistent with this same quartic.
The cavity potential is 1/2[Ec Q+sum a_i q_i]², a_i=2g0 sqrt(E0/Ec)
in mass-weighted coordinates. Minimizing Q leaves the bare molecular PES unchanged.
DSE omission would generate an artificial static barrier change or instability.
The 2D PES output is an N=1 slice; it is not a many-body reaction pathway.

At the TS, only one molecular curvature becomes negative; spectator wells remain
unchanged. Exactly one unstable normal mode is excluded. Quantum harmonic TST uses
k=(kBT/h) exp(-B/kBT) prod_min[2sinh(Ek/2kBT)] /
prod_TS_stable[2sinh(Ek/2kBT)]. ZPE is already in these partition functions and
is not added again as a separate barrier correction. The stable quantum free
energy and changed normal-mode dividing surface can change the rate without a
relaxed classical PES barrier change. The classical limit equals the coupled
unstable frequency divided by the bare barrier frequency, relative to outside.
The resulting profiles are calculated, not multiplied by an assigned resonance
Lorentzian. A sharp minimum at the reactant frequency is not guaranteed.
Scans cover both reactant and barrier frequencies. No tunneling, solvent friction,
exact reactive trajectories or quantum recrossing calculation is included, so the
absolute rates are illustrative harmonic-TST estimates, especially for this high barrier.

The dilution scan holds collective g fixed, hence g0=g/sqrt(N). A local molecule's
bright-state participation is 1/N, while the photon sees the sum of dipoles.
This explains dilution within this harmonic equilibrium model; it does not resolve
every experimental collective-chemistry observation. The fixed-g0 Rabi scan and
fixed-g collective-local scan vary different physical parameters and must not be
conflated. Frequency scans hold g fixed, not cavity volume/dipole coupling fixed.

## 21C/D: spectroscopy and geometry

The passive two-port response inverts K-E²I-i E diag(gamma,kappa), retaining DSE
and counterrotating physics even for g/E0>0.1. kappa=Ec/Q, Q=200 by default,
gamma=3 meV. Symmetric ports give T, R and true absorptance A=1-R-T.
The reported splitting-versus-linewidth test is a chosen resolvability criterion,
not a universal definition of strong coupling. Spectral input-output damping is
phenomenological, not a microscopic ultrastrong-coupling bath master equation.
The angular cavity energy is Ec(0)/sqrt(1-sin²(theta)/n_eff²). This is a
Fabry-Perot angular-dispersion proxy, not a calibrated ATR prism/multilayer model.
An energy-angle heatmap and weak-probe FTIR spectra demonstrate the polariton
branches. Ideal identical dark modes are counted in the eigenvalue spectrum but
have zero direct cavity oscillator strength; no false dark absorption peak is added.
Disorder, direct molecular illumination or symmetry breaking would be needed to
make them visible in a chosen instrument. All traces are simulated.

## Interpretation and validation

Boundary conditions modify hybrid normal modes and equilibrium fluctuations in
this model. Optical anticrossing alone does not establish a chemical rate change.
Tests check Fock convergence, exact sqrt(N) scaling, dark degeneracy, DSE stability,
the relaxed-PES cancellation, classical determinants, uncoupled rate recovery and
optical passivity. Reported errors are model-internal validation, not experimental
accuracy. Molecular anharmonicity, realistic dipoles, losses in rate dynamics and
repeated chemical measurements are needed before attributing catalysis to vacuum
fields in an actual material. Earlier phases are preserved, not revalidated here.

## Computed results / 实际计算结果

| Quantity / 项目 | Value / 数值 |
|---|---:|
| Fundamental / 分子振动能量 (eV) | 0.18 |
| Collective coupling / 集体耦合 (eV) | 0.018 |
| LP / UP (eV) | 0.16289776 / 0.19889776 |
| Rabi splitting / 劈裂 (meV) | 36.000000 |
| gamma + kappa (meV) | 3.900000 |
| Resolved strong coupling / 满足所用线宽判据 | True |
| N=1 resonant quantum rate ratio / 共振处量子速率比 | 0.98780387 |
| N=64 resonant quantum rate ratio / 共振处量子速率比 | 0.99981325 |
| N=1 swept ratio range / 扫描范围 | 0.98375770 to 0.99439723 |
| Truncated Fock mode error / 小基组误差 (eV) | 1.968e-07 |

## Run / 运行

```bash
python -m pip install -r requirements_phase21.txt
python run_phase21_cavity_qed_polaritonic_chemistry.py --self-test
python run_phase21_cavity_qed_polaritonic_chemistry.py --output-dir .
```

`results_phase21/summary.json` includes all assigned parameters, environment versions,
numerical errors and tolerances. CSV files preserve modes, scaling, rate sweeps,
local dilution, the 2D PES, FTIR and the angle-resolved transmission grid.
`results_phase21/sha256.json` records generated file and script checksums.
Only numpy, scipy and matplotlib are required. No network calls occur during a run.

![Levels](figures_phase21/fig1_polaritonic_pauli_fierz_energy_levels.png)
![Dispersion](figures_phase21/fig2_vacuum_rabi_dispersion_anticrossing.png)
![Rates](figures_phase21/fig3_cavity_reaction_rate_resonance_profile.png)
![FTIR](figures_phase21/fig4_simulated_ftir_polariton_spectra.png)

## References / 参考文献

- [Li, Mandal and Huo: Cavity frequency-dependent theory for vibrational polariton chemistry](https://doi.org/10.1038/s41467-021-21610-9)
- [Wang, Flick and Yelin: Chemical reactivity under collective vibrational strong coupling](https://arxiv.org/abs/2206.08937)
- [Lindoy, Mandal and Reichman: Investigating the collective nature of cavity-modified chemical kinetics](https://doi.org/10.1515/nanoph-2024-0026)
