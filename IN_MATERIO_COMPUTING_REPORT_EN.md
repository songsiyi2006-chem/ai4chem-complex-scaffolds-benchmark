# Phase 23 — In-materio computing: a reproducible effective-model benchmark

## Scope and results

This is a classical numerical simulation of a hypothetical volatile redox junction,
not a fabricated transition-metal bis-terpyridine/polyoxometalate device, a quantum
transport calculation, or experimental evidence of AGI. Coordination chemistry can
motivate tunable redox and ion-relaxation properties; the material-specific parameters
here are assigned. No atomistic electronic structure or solvent/temperature dependence
is computed. Physical reservoirs may reduce data movement for selected tasks; neither
an insurmountable silicon barrier nor a universal advantage over digital computers follows.

| Quantity / 指标 | Computed value / 计算值 |
|---|---:|
| One-step test NMSE (mean x,y,z) | 0.25545867 |
| One-step test NMSE x / y / z | [0.10287173447306878, 0.26713650512751624, 0.39636777835630116] |
| Ten-step test NMSE (mean x,y,z) | 0.98501813 |
| One-step persistence baseline, x only | 0.010697335 |
| Device dissipation J / input sample | 2.9845378e-09 |
| Device + wire + source-resistor J / sample | 1.3640139e-08 |
| Full system J/op | Not calculated / 未计算 |
| High-frequency log(area)-log(f) slope | -0.999997 |
| Hysteresis area refinement relative error | 1.71153e-05 |
| Positive / negative STDP sign tests | True / True |
| KCL max residual (A) | 9.98575e-14 |
| Halved network-step max state difference | 0.000129068 |
| Halved plasticity-step max normalized-window difference | 1.60583e-07 |

## A. Correcting and deriving the device model

The requested current G(w) V sinh(alpha V) is even in V and gives I V < 0 for V < 0,
which is not a passive two-terminal dissipative junction. Keeping G in siemens, we use

    G(w) = G_low + (G_high-G_low) w
    I(V,w) = G(w) sinh(alpha V) / alpha
    q(V) = sign(V) max(abs(V)-V_threshold, 0)
    dw/dt = eta [1-(2w-1)^(2p)] sinh(q/V0) - (w-w_eq)/tau

The antisymmetric exponential drift is Butler-Volmer-inspired: the difference of
forward/backward exponential rates produces a sinh term. It is not a quantitatively
derived Butler-Volmer kinetic model: V0 is an effective scale, and the state rates
are not tied to electronic current by a Faradaic conservation law. w_eq=0.5,
p=1 and tau=10 ms. Thermal noise, reaction free energies and chemical heat are absent.
Relaxation toward an interior equilibrium removes the absorbing zero-state boundary
of the original window model. At w=0, dw/dt=+w_eq/tau; at w=1 it is negative.
These inward directions preserve the interval. RK4 subdivides held-voltage steps
according to a rate bound and raises on bound violations; it never clips states.

I(0,w)=0 and I V>=0 exactly. A pinched curve is consistent with memristive dynamics
but not by itself a unique experimental identification of a real memristor.
Periodic-orbit shooting removes transient initialization. The unsigned two-lobe
area is reported rather than a whole-loop signed area that can cancel.
At bounded sinusoidal drive, the state swing becomes O(1/f); at infinite frequency
the frozen-state I(V) is single-valued. Numerical slopes test the finite-frequency trend,
not the experimental universality of the model. Parameters are in summary.json.

## B. Pulse-dependent plasticity, with an explicit retention limitation

Pre/post waveforms have a +0.35 V, 1 ms head and a -0.15 V exponential tail
with an 8 ms decay, truncated at 51 ms. Junction voltage is V_post-V_pre.
The 101 delays span -50 to +50 ms. The same state equation is integrated for every
delay; no signed exponential learning curve is inserted into the dynamics.
Reported Delta W=int[G_pair(t)-G0]dt has units S s, not S.
Its dimensionless normalization is Delta W/(G0 tau), not Delta G/G0.
All pulses and at least 69 ms of final relaxation fit inside the common 230 ms
integration window. Tail fits use 3-40 ms and report log-space fit quality, not
experimental fidelity. Endpoint states and peak deviations are separately exported.
Reversing electrode orientation reverses the order-dependent window.

The PPF ratio is G_after_second/G_after_first for positive 1 ms pulses; it is not
a ratio of baseline-subtracted increments. At zero voltage the state returns to
w_eq. Therefore these results represent volatile STP/STDP-like updates, NOT
durable LTP/LTD. A second slow/nonvolatile state and retention measurements would
be needed to claim long-term synaptic learning. A pulse-shaped curve is not such proof.

## C. Circuit and reservoir learning

Twenty-five junctions shunt a 5x5 resistive grid (40 edges) to ground. Each node has
a finite source resistor and a 2 pF parasitic capacitor. Fixed random voltage gains,
biases, conductances, rates and 3-30 ms relaxation times represent assigned hardware
heterogeneity. A scalar, training-normalized Lorenz x signal is distributed by this
fixed mask; no delayed input, polynomial expansion or engineered software features
are used. This topology is a mesh with crossbar-like parasitics, not a lithographically
resolved row/column crossbar or an atomistic molecular network.

Each substep solves by Newton iteration:

    (C/h)(v-v_previous) + G_source(v-u) + G_wire L v + I(v,w) = 0.

The Jacobian is positive definite. Molecular states then advance by RK4 at held
nodal voltage. This is first-order operator splitting with backward-Euler capacitors,
not a globally fourth-order circuit solver. A halved substep comparison and maximum
KCL residual are exported. The fading-memory check compares two initial states
under the same constant source input; it does not establish global echo-state stability.

Lorenz parameters are sigma=10, rho=28, beta=8/3. A dimensionless 0.02 interval maps
to one assigned physical 0.5 ms sample; this clock conversion is not an experimentally
measured speed. The first 1500 generated samples are discarded. There are 3600 inputs,
with washout 200, train end 2200, validation end 2800 and 800 held-out test samples.
Ten-sample gaps exclude cross-boundary target leakage. Input normalization and readout
centering/scaling use TRAIN ONLY. Nine ridge penalties are compared on validation only;
no reservoir parameters or readout weights are tuned using test targets.

Only an affine linear readout is fitted by SVD ridge regression; centering is folded
into a single 25x3 weight matrix plus a bias. Targets are x,y,z at +1 and +10 steps.
NMSE is MSE divided by each TEST coordinate's variance; their mean is also reported.
The x-only persistence and input-only linear baselines are exported separately.
All predictions are input-driven (teacher-forced), NOT autonomous chaotic rollouts.
A low one-step x error can be obtained by persistence; three-coordinate and ten-step
scores must not be conflated with that easier task. There are no uncertainty bars
from independent chaotic realizations or device-noise ensembles in this deterministic run.

## Energy accounting, not a sub-femtojoule assumption

An operation is one input sample of the entire 25-node network, not one transistor
event. Device int V I dt, wire loss and source-resistor loss are summed separately.
Supplied source energy, capacitor stored-energy change and backward-Euler numerical
loss are recorded and checked in a discrete energy balance. The sampling integral
uses the same right-endpoint held-state solution as KCL. The residual verifies this
discretization, not exact continuous-time energy accuracy.
No ADC, DAC, input-mask distribution, output readout, training, sensing noise,
reset, fabrication or chemical reaction heat is included. Hence even the circuit
subtotal is not total system energy, and dividing by an arbitrary number of internal
events to manufacture a sub-fJ/op number is not permitted. The actual computed
threshold flags are retained, including failures.

## D. Small-signal electrochemical impedance

At each DC bias (-0.25, 0, +0.25 V), solve F(w,V)=0. Linearizing gives

    Y_j(omega) = I_V + I_w F_V/(i omega-F_w)
    Z_W = R_W tanh(sqrt(i omega tau_D))/sqrt(i omega tau_D)
    Z = R_s + [i omega C_geom + (1/Y_j + Z_W)^(-1)]^(-1).

The scan uses f=0.01-1e6 Hz and omega=2 pi f, resolving the request's frequency/angular
frequency ambiguity. Nominal 10 mV phasor currents are saved. The finite-length,
transmissive Warburg branch has a finite DC resistance and an intermediate 45-degree
asymptote; it is not a blocking-boundary divergent pseudocapacitor.
R_s=500 ohm, C_geom=2 nF, R_W=30 kilohm and tau_D=0.5 s are explicitly assigned
film/electrolyte parameters, not the network's per-node parasitics. The single redox
state cannot generate a diffusion continuum by itself, so the Warburg branch is
an added equivalent-circuit assumption. Frozen-state R_ct=1/I_V is not the full DC
differential resistance. The code records stationary stability and a 10 mV secant
linearity diagnostic, not a nonlinear harmonic simulation or inverse deconvolution.
Nyquist/Bode curves are synthetic and do not identify a ligand or validate a potentiostat.
The biased state-response term can produce an inductive-sign loop; not every arc
is a purely capacitive semicircle and no fitted experimental component separation is claimed.

## Outlook and reproducibility

Useful next experiments are pulse/retention and temperature sweeps on a specified
redox film, simultaneous electrical and chemical state measurements, independently
measured line/capacitor parasitics, and noise-inclusive held-out edge-sensing benchmarks.
Molecular-electronic sensor interfaces are a plausible research application;
sub-fJ AGI coprocessors remain a proposal, not an output of this benchmark.

Run `python -m pip install -r requirements_phase23.txt`, then the script with
`--self-test` and without it. `--output-dir PATH` selects a dedicated output folder.
The script has no Git/network side effects; publishing is a separate reviewed action.
Results include CSV traces, readout matrices, split indices, all validation trials,
summary/threshold flags, source/data SHA-256 hashes and four 300 DPI PNGs.

## References / 参考文献

1. Joglekar & Wolf, *The elusive memristor: properties of basic electrical circuits*:
   https://arxiv.org/abs/0807.3994 — window-model background, not parameter calibration.
2. Du et al., *Reservoir computing using dynamic memristors for temporal information processing*:
   https://www.nature.com/articles/s41467-017-02337-y — physical short-term memory and readout-only learning.
3. Zhong et al., *Dynamic memristor-based reservoir computing for high-efficiency temporal signal processing*:
   https://www.nature.com/articles/s41467-020-20692-1 — experimental reservoir work; its measured accuracy is not ours.
4. Goswami et al., *Decision trees within a molecular memristor*:
   https://www.nature.com/articles/s41586-021-03748-0 — molecular redox-state computing; a different device/task.
5. Gamry Instruments, *Basics of electrochemical impedance spectroscopy*:
   https://www.gamry.com/application-notes/EIS/basics-of-electrochemical-impedance-spectroscopy/
   — small-signal equivalent circuits and diffusion impedance.
