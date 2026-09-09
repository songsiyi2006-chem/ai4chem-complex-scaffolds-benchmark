# Phase 22: molecular spin qubits, clock transitions and dynamical decoupling

## Scope and molecular model

This is an assigned S=1/2, I=7/2 coordination-complex analogue, not a chemical
structure calculation or fitted vanadium/lanthanide material. Electronic/nuclear
Zeeman, anisotropic hyperfine and nuclear quadrupole tensors are included in Hz.
The general tensor implementation contains ZFS; for S=1/2 its traceless part
has no physical splitting and D is correctly set to zero. Nuclear Zeeman sign
and radian/cycle conversions are explicit. Bz spans 0..1 T. Sparse Kronecker
operators build the small central Hilbert space, then dense Hermitian diagonalization
is appropriate. Transition slopes use Hellmann-Feynman matrix elements;
curvature includes virtual transitions through all other central eigenstates.
Interior derivative roots are refined and filtered by microwave oscillator strength.
Near-degenerate nondegenerate-perturbation failures are rejected.

## What clock protection means

At a CT, dnu/dBz vanishes; curvature generally does not. This protects against
first-order longitudinal field noise, not arbitrary magnetic fluctuations,
hyperfine/strain noise, transverse coupling, spin flips or phonons. The two chosen
levels are a noise-insensitive operating pair, NOT a rigorously decoherence-free
subspace or an implemented error-correcting code. The report tests a local Taylor
expansion and does not prove macroscopic immunity. Sorted adiabatic states label
the field map; this is not a diabatic pulse-tracking calculation.

## Bath dynamics and approximation hierarchy

Thirty protons are sampled uniformly by volume in a 3.5..12 Angstrom shell with
2 Angstrom pair exclusion and a fixed seed. This is synthetic solvent geometry,
not crystallography. Nuclear dipolar interactions retain the full Cartesian
tensor, including nonsecular terms. Proton Zeeman and dipolar couplings are in Hz.
The central conditional energies use E_a(B+beta)=E_a(B)+E'_a beta+E''_a beta²/2,
where beta is the longitudinal dipolar field operator of the protons. The effective
local-field model applies to the total central axial moment; it is not a full
central-spin flip Hamiltonian. Off-diagonal central transitions and relaxation
are excluded. Quadratic terms are retained so a CT does not spuriously switch
all decoherence off. Their cross-nuclear terms enter pair clusters.

The bath initial state is maximally mixed (high nuclear-spin temperature).
Each one/two-spin cluster evolves exactly under two conditional Hamiltonians,
L_C=Tr(U1†U0)/dimension. Ensemble CCE-2 multiplies singles and irreducible pairs,
L=prod L_i prod[L_ij/(L_i L_j)]. It retains finite-memory quantum evolution,
but is not the exact 30-spin density matrix: that would have 2^60 complex entries,
about 18.4 exabytes in complex128. CCE-1/CCE-2 differences and a three-spin exact
reference are diagnostics, not proof of global CCE convergence at clock points.
Outside-cluster fields and higher clusters are omitted. Small denominators or
|L|>1 beyond tolerance invalidate the remaining curve; no clipping hides failure.

Independent assigned quasistatic technical Bz noise (20 microtesla rms) is averaged
analytically through first and second derivatives for Ramsey/FID. Balanced ideal
echo pulses cancel this static term. An assigned T1=10 ms contributes exp(-t/2T1)
to all curves; it is not a spin-phonon calculation or a measured molecular T1.

## Control and lifetimes

The projected rotating-wave drive derives its Rabi rate from the actual selected
Sx matrix element. Gauss-Hermite quadrature averages detuning and 0.5% amplitude
noise for average single-qubit gate fidelity, (|Tr(Utarget†U)|²+2)/6. This includes
neither leakage to other molecular levels nor bath evolution during the pulse.
It must not be quoted as experimental gate fidelity. Bloch trajectories are ideal
two-level rotations. No 99.9% threshold is enforced.

CPMG and UDD use ideal instantaneous pi pulses, with n=1,4,16,64. Conditional
propagators implement pulse toggling; the repeated CPMG block is exponentiated
exactly, while nonuniform UDD intervals are explicitly multiplied. Filter-function
numerators are exported, not used as a surrogate for the CCE calculation. These
are frequency-selective filters, not an absolute Nyquist cutoff. Finite pulse width,
pulse accumulation errors and hardware bandwidth are outside this comparison.
The gate duration is reported separately and is not inserted into ideal DD.
Lifetimes are the FIRST 1/e crossing on the finite grid, not fitted exponential
T2 values. Revivals are possible. Right-censored and CCE-unresolved curves do not
receive invented lifetimes or a claimed >100x gain. A gain ratio is only reported
when both FID and CPMG-64 crossings are resolved at the SAME field/noise settings.

## EPR, chemistry and perspective

The field-dependent lifetime map uses CCE-2 for FID and Hahn echo. The denser
spectroscopic map uses CCE-1 explicitly to keep the calculation tractable; it
captures conditional nuclear precession but omits pair correlations in that map.
Mean subtraction and a Hann window precede the real FFT. The axis is modulation
frequency versus total evolution time 2tau, not a unique isotope/structure assignment.
Normalized echo quadrature is Re L; conversion to electron Sy units multiplies
by the projected readout matrix element. A synthetic Ramsey detuning is 100 kHz.
No Bruker hardware response, experimental ESEEM data, noise calibration or ligand
structure inversion is claimed. The CT need not show the largest ESEEM amplitude.

Coordination chemistry offers precise molecular identity and tunable ligands,
but does not guarantee uniform crystal/solvent environments or capabilities
unachievable by semiconductor defects. Deuteration lowers gyromagnetic coupling
but introduces I=1 quadrupoles; 12C enrichment removes 13C spins, while other
nuclei, phonons and disorder remain. These are experimental design directions,
not simulations performed here or evidence of room-temperature processors.
Molecular clock-transition experiments motivate this benchmark but their fitted
parameters and measured coherence times are not transplanted into its outputs.

## Computed values / 计算结果

- Selected CT / 选定时钟跃迁: **0.0027208157 T**, states 3 and 11 (ascending energies).
- Transition / 跃迁频率: **1.18649392 GHz**.
- Axial slope / 轴向一阶导数: **0 Hz/T**.
- Curvature / 二阶导数: **6.78178e+12 Hz/T²**.
- pi pulse / π 脉冲: **6.002730 microseconds**.
- Projected noisy gate fidelity / 投影含噪声门保真度: **99.943046%**.
- CPMG-64/FID crossing ratio / 阈值交点增益: **8.85831669775934** (null = not resolvable).

| Sequence / 序列 | Lifetime (seconds) / 寿命 | Status / 状态 |
|---|---:|---|
| FID | 7.822288533831443e-05 | crossing |
| Hahn | 0.00373053789404695 | crossing |
| CPMG1 | 0.00373053789404695 | crossing |
| CPMG4 | 0.004996799101433063 | crossing |
| CPMG16 | 0.00011541593097462193 | crossing |
| CPMG64 | 0.0006929230913393051 | crossing |
| UDD1 | 0.003730537894047006 | crossing |
| UDD4 | 0.004882192833990317 | crossing |
| UDD16 | 0.0034304853590400933 | crossing |
| UDD64 | 0.000365455700450016 | crossing |

## Reproduction / 复现

```bash
python -m pip install -r requirements_phase22.txt
python run_phase22_molecular_spin_qubits.py --self-test
python run_phase22_molecular_spin_qubits.py
```

All parameters and numerical checks are in `results_phase22/summary.json`.
CSV files contain the bath geometry, CT candidates, field levels, full complex
coherences (including UDD-1/4/16/64), filter functions, ESEEM and lifetime statuses.
`sha256.json` records output and script hashes. Figures are 300 DPI.
The default run uses a fixed seed and only numpy, scipy and matplotlib.

![Levels](figures_phase22/fig1_breit_rabi_clock_transitions.png)
![Control](figures_phase22/fig2_rabi_nutation_and_bloch_sphere.png)
![Coherence](figures_phase22/fig3_dynamical_decoupling_coherence_gain.png)
![EPR](figures_phase22/fig4_pulsed_epr_eseem_twin.png)

## References / 参考文献

- [Shiddiq et al. (2016), Enhancing coherence in molecular spin qubits via atomic clock transitions](https://doi.org/10.1038/nature16984)
- [Yang and Liu, finite-spin-bath cluster-correlation expansion](https://arxiv.org/abs/0806.0098)
- [Yang and Liu, ensemble CCE dynamics](https://arxiv.org/abs/0902.3055)
- [Cluster-correlation expansion for studying decoherence of clock transitions in spin baths](https://arxiv.org/abs/2007.00412)
