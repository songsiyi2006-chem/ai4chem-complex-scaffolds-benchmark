#!/usr/bin/env python3
"""Phase 21: dipole-self-energy-complete harmonic Pauli-Fierz benchmark.

Python >=3.10; numpy, scipy, matplotlib. No downloads or earlier-phase runtime.
Run --self-test, then run with --output-dir PATH (default: script directory).
The reactive PES and dipole slopes are assigned model parameters, not reconstructed
Phase 4/5 molecular electronic structure. Spectra are weak-probe simulations;
equilibrium rate calculations have no drive. No rate-resonance target is imposed.
"""
from __future__ import annotations
import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import platform
import unittest

import numpy as np
import scipy
from scipy.linalg import eigh, eigvalsh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm

KB=8.617333262145e-5  # eV/K
HBAR=6.582119569e-16  # eV s
PLANCK=2*np.pi*HBAR
EV_CM=8065.543937
REFS=[
    ('Li, Mandal and Huo: Cavity frequency-dependent theory for vibrational polariton chemistry',
     'https://doi.org/10.1038/s41467-021-21610-9'),
    ('Wang, Flick and Yelin: Chemical reactivity under collective vibrational strong coupling',
     'https://arxiv.org/abs/2206.08937'),
    ('Lindoy, Mandal and Reichman: Investigating the collective nature of cavity-modified chemical kinetics',
     'https://doi.org/10.1515/nanoph-2024-0026'),
]

@dataclass(frozen=True)
class Config:
    vibration_ev: float=.18
    collective_g_ev: float=.018
    quality_factor: float=200.
    gamma_ev: float=.003
    temperature_k: float=300.
    barrier_ev: float=2.012  # illustrative quartic, about 46.4 kcal/mol
    minimum_q: float=1.0  # dimensionless reaction coordinate
    photon_states: int=4  # occupations 0..3
    vibration_states: int=5  # occupations 0..4
    energy_points: int=601
    angle_points: int=141
    cavity_normal_ev: float=.16
    effective_index: float=1.3

    @property
    def barrier_frequency_ev(self):
        # Symmetric quartic has |V''(TS)|=V''(minimum)/2.
        return self.vibration_ev/np.sqrt(2)

    def validate(self):
        if not all(np.isfinite(v) for v in asdict(self).values()):
            raise ValueError('Parameters must be finite')
        if min(self.vibration_ev,self.gamma_ev,self.temperature_k,self.barrier_ev,self.minimum_q,self.cavity_normal_ev)<=0:
            raise ValueError('Energies, temperature and geometry must be positive')
        if self.collective_g_ev<0 or not 50<=self.quality_factor<=500 or self.effective_index<=1:
            raise ValueError('g>=0, Q in [50,500], effective index>1 required')
        if self.photon_states<4 or self.vibration_states<5 or self.energy_points<101 or self.angle_points<31:
            raise ValueError('Insufficient Fock or spectral sampling')


def coordinate(n):
    """Projected b+b† and projected (b+b†)^2, including boundary correction.

    Squaring the truncated coordinate directly loses n at the top boundary.
    """
    x=np.diag(np.sqrt(np.arange(1,n)),1)+np.diag(np.sqrt(np.arange(1,n)),-1)
    x2=np.diag(2*np.arange(n)+1.)
    if n>2:
        v=np.sqrt(np.arange(1,n-1)*np.arange(2,n))
        x2+=np.diag(v,2)+np.diag(v,-2)
    return x,x2


def embed(matrix,index,dims):
    a=np.ones((1,1))
    for k,n in enumerate(dims): a=np.kron(a,matrix if k==index else np.eye(n))
    return a


def fock_hamiltonian(w0,wc,g0,nmol=1,nph=4,nvib=5):
    """Full product basis for small N. H=H0+g0 Xc sum Xi+g0²/wc(sum Xi)².

    Counterrotating terms and all intermolecular DSE terms are retained.
    General large-N harmonic ensemble is handled by exact bright/dark reduction.
    """
    if nmol not in (1,2): raise ValueError('Full Fock basis limited to N=1 or 2; use bright/dark reduction for large N')
    dims=[nvib]*nmol+[nph];dim=int(np.prod(dims))
    h=np.zeros((dim,dim))
    x,x2=coordinate(nvib); xc,_=coordinate(nph)
    xs=[embed(x,i,dims) for i in range(nmol)]
    for i in range(nmol):
        h+=embed(np.diag(w0*(np.arange(nvib)+.5)),i,dims)
        h+=(g0*g0/wc)*embed(x2,i,dims)
    h+=embed(np.diag(wc*(np.arange(nph)+.5)),nmol,dims)
    sumx=sum(xs)
    h+=g0*sumx@embed(xc,nmol,dims)
    for i in range(nmol):
        for j in range(i): h+=(2*g0*g0/wc)*(xs[i]@xs[j])
    return h


def stiffness(w0,wc,g0,nmol=1,transition=False):
    """Mass-weighted Hessian, units eV² with hbar=1.

    V=V_mol+1/2(wc Q+sum a_i q_i)^2; a_i=2g0 sqrt(w0/wc).
    Local TS changes ONLY molecular coordinate 0 to negative curvature.
    The dipole slope is held fixed across the local barrier.
    """
    a=np.full(nmol,2*g0*np.sqrt(w0/wc))
    bare=np.full(nmol,w0*w0)
    if transition: bare[0]=-w0*w0/2
    k=np.zeros((nmol+1,nmol+1))
    k[:nmol,:nmol]=np.diag(bare)+np.outer(a,a)
    k[:nmol,-1]=wc*a;k[-1,:nmol]=wc*a;k[-1,-1]=wc*wc
    return k


def polaritons(w0,wc,g):
    return np.sqrt(eigvalsh(stiffness(w0,wc,g)))


def log_partition_denominator(energies,temp):
    x=np.asarray(energies)/(2*KB*temp)
    return float(np.sum(x+np.log(-np.expm1(-2*x))))  # log[2 sinh(x)]


def rates(cfg,wc,nmol,g):
    """Quantum harmonic TST at coupled normal-mode dividing surface.

    k=(kBT/h) exp(-B/kBT) prod[2sinh(Emin/2kBT)]/prod_stable[2sinh(ETS/2kBT)].
    One negative TS mode is excluded. No separate ZPE correction is added,
    since these oscillator partition functions already contain ZPE.
    Classical limit ratio = |omega_unstable,cav|/omega_barrier.
    This does not include tunneling, dissipative recrossing or exact dynamics.
    """
    g0=g/np.sqrt(nmol)
    km=stiffness(cfg.vibration_ev,wc,g0,nmol)
    kt=stiffness(cfg.vibration_ev,wc,g0,nmol,True)
    vm=eigvalsh(km);vt=eigvalsh(kt)
    if np.any(vm<=0) or np.count_nonzero(vt<0)!=1: raise RuntimeError('Invalid minimum or TS inertia')
    em=np.sqrt(vm);et=np.sqrt(vt[vt>0])
    lq=log_partition_denominator(em,cfg.temperature_k)-log_partition_denominator(et,cfg.temperature_k)
    outside=log_partition_denominator([cfg.vibration_ev],cfg.temperature_k)
    logk=np.log(KB*cfg.temperature_k/PLANCK)-cfg.barrier_ev/(KB*cfg.temperature_k)+lq
    ratio=float(np.exp(lq-outside))
    classical=float(np.sqrt(-vt[0])/cfg.barrier_frequency_ev)
    # DSE completion of square restores the original PES after minimizing Q.
    return ratio,classical,float(np.exp(logk))


def transmission(cfg,energy,wc,g):
    """Vectorized damped coupled-oscillator two-port response, beyond RWA.

    It is an angle-dispersive Fabry-Perot input-output proxy, not a prism ATR
    transfer-matrix calculation. kappa and gamma are energy linewidths (eV).
    Port fractions 1/2,1/2; absorption is vibrational damping, not 1-T.
    """
    e=np.asarray(energy);wc=np.asarray(wc)
    w0=cfg.vibration_ev;kappa=wc/cfg.quality_factor
    amm=w0*w0+4*g*g*w0/wc-e*e-1j*e*cfg.gamma_ev
    acc=wc*wc-e*e-1j*e*kappa
    off=2*g*np.sqrt(w0*wc)
    chi=amm/(amm*acc-off*off)
    t=1j*e*kappa*chi
    r=1+1j*e*kappa*chi
    T=abs(t)**2;R=abs(r)**2;A=1-T-R
    if np.min(A)<-1e-10 or np.max(T)>1+1e-10: raise RuntimeError('Nonpassive optical response')
    return T,R,np.maximum(A,0)


def pes(cfg,q,photon_q,wc,g):
    """Dimensionless quartic molecular q mapped to mass-weighted coordinate.

    Effective mass M=V''(well)/w0² fixes the prescribed fundamental.
    photon_q units eV^(-1/2) in hbar=1 mass-weighted convention.
    """
    v=cfg.barrier_ev*(1-(np.asarray(q)/cfg.minimum_q)**2)**2
    mass=8*cfg.barrier_ev/(cfg.minimum_q**2*cfg.vibration_ev**2)
    a=2*g*np.sqrt(cfg.vibration_ev/wc)
    return v+.5*(wc*np.asarray(photon_q)+a*np.sqrt(mass)*np.asarray(q))**2


def validate(cfg):
    w=cfg.vibration_ev;g=cfg.collective_g_ev
    analytic=polaritons(w,w,g)
    small=eigvalsh(fock_hamiltonian(w,w,g,nph=cfg.photon_states,nvib=cfg.vibration_states))
    large=eigvalsh(fock_hamiltonian(w,w,g,nph=12,nvib=14))
    # Bright LP/UP are the two lowest excited levels for default g/w=0.1.
    # General mode selection uses photon oscillator transition matrix elements.
    def mode_gaps(nph,nv):
        vals,vec=eigh(fock_hamiltonian(w,w,g,nph=nph,nvib=nv))
        xc,_=coordinate(nph); op=np.kron(np.eye(nv),xc)
        strength=abs(vec.T@(op@vec[:,0]))**2
        nearest=[]
        for expected in analytic:
            candidate=np.where(strength>1e-5)[0]
            index=candidate[np.argmin(abs((vals[candidate]-vals[0])-expected))]
            nearest.append(vals[index]-vals[0])
        return np.array(nearest)
    low=mode_gaps(cfg.photon_states,cfg.vibration_states);high=mode_gaps(12,14)
    ensemble=eigvalsh(fock_hamiltonian(w,w,g/np.sqrt(2),nmol=2))
    expected=np.sort(np.r_[analytic,w])
    ensemble_error=float(np.max(abs(ensemble[1:4]-ensemble[0]-expected)))
    counts=np.array([1,4,9,16,25]);g0=.0015
    gaps=np.array([np.diff(polaritons(w,w,g0*np.sqrt(n)))[0] for n in counts])
    q=np.linspace(-1,1,51);mass=8*cfg.barrier_ev/(cfg.minimum_q**2*w*w)
    optimum=-2*g*np.sqrt(w/w)*np.sqrt(mass)*q/w
    bare=cfg.barrier_ev*(1-q*q)**2
    checks={
        'fock_low_mode_error_ev':float(np.max(abs(low-analytic))),
        'fock_high_mode_error_ev':float(np.max(abs(high-analytic))),
        'fock_ground_convergence_ev':float(abs(small[0]-large[0])),
        'N2_product_basis_mode_error_ev':ensemble_error,
        'rabi_sqrtN_error_ev':float(np.max(abs(gaps-2*g0*np.sqrt(counts)))),
        'DSE_relaxed_PES_error_ev':float(np.max(abs(pes(cfg,q,optimum,w,g)-bare))),
        'uncoupled_quantum_rate_error':abs(rates(cfg,w,8,0)[0]-1),
        'uncoupled_classical_rate_error':abs(rates(cfg,w,8,0)[1]-1),
    }
    tolerance={k:1e-9 for k in checks}
    tolerance.update(fock_low_mode_error_ev=5e-4,fock_ground_convergence_ev=5e-5,N2_product_basis_mode_error_ev=5e-4,fock_high_mode_error_ev=1e-8)
    if any(not np.isfinite(v) or v>tolerance[k] for k,v in checks.items()):
        raise RuntimeError(f'Validation failed: {checks}')
    return {'status':'passed','errors':checks,'tolerances':tolerance},analytic,low,high,counts,gaps


class Tests(unittest.TestCase):
    def test_dse_boundary(self):
        x,x2=coordinate(5)
        self.assertAlmostEqual(x2[-1,-1]-(x@x)[-1,-1],5.)
    def test_fock_and_scaling(self): self.assertEqual(validate(Config())[0]['status'],'passed')
    def test_dark_modes(self):
        cfg=Config();w=cfg.vibration_ev
        vals=np.sqrt(eigvalsh(stiffness(w,w,cfg.collective_g_ev/4,16)))
        self.assertEqual(np.sum(abs(vals-w)<1e-10),15)
    def test_passivity(self):
        for quality in (50,200,500):
            cfg=Config(quality_factor=quality)
            for g in (0,.018,.04):
                t,r,a=transmission(cfg,np.linspace(.03,.4,4001),.18,g)
                self.assertTrue(np.allclose(t+r+a,1,atol=1e-10))
        t,r,a=transmission(Config(),.18,.18,0)
        self.assertAlmostEqual(float(t),1);self.assertAlmostEqual(float(a),0)
    def test_classical_determinant(self):
        cfg=Config();w=cfg.vibration_ev;wc=.13;n=6;g0=.018/np.sqrt(n)
        km=stiffness(w,wc,g0,n);kt=stiffness(w,wc,g0,n,True)
        self.assertAlmostEqual(float(np.linalg.det(km)/(wc**2*w**(2*n))),1,places=10)
        self.assertEqual(np.sum(eigvalsh(kt)<0),1)
    def test_rate_null_and_dilution(self):
        cfg=Config()
        self.assertAlmostEqual(rates(cfg,.11,16,0)[0],1,places=12)
        self.assertLess(abs(rates(cfg,.18,64,.018)[0]-1),abs(rates(cfg,.18,1,.018)[0]-1))
    def test_no_dse_instability(self):
        w=.18;g=.12
        k=stiffness(w,w,g);self.assertTrue(np.all(eigvalsh(k)>0))
        k[0,0]-=4*g*g
        self.assertTrue(np.any(eigvalsh(k)<0))
    def test_bad_config(self):
        for cfg in (Config(quality_factor=20),Config(collective_g_ev=-1),Config(temperature_k=0)):
            with self.assertRaises(ValueError):cfg.validate()


def figures(cfg,data,out):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    def save(fig,name): fig.savefig(out/name,dpi=300,facecolor='white');plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(10,4.8),layout='constrained')
    for val,label,color in zip(data['analytic'],('LP','UP'),('tab:blue','tab:red')):
        ax[0].hlines(val,0,.8,color=color,lw=3,label=label)
    ax[0].hlines(cfg.vibration_ev,1,1.8,color='black',lw=3,label='N-1 dark modes (N=16)')
    ax[0].set(xlim=(-.1,1.9),ylabel='Excitation energy (eV)',title='21A | PF modes with dipole self-energy',xticks=[.4,1.4],xticklabels=['Bright + cavity','Dark reservoir'])
    ax[0].legend(fontsize=8)
    ax[1].plot(np.sqrt(data['counts']),data['gaps']*1000,'o',label='PF diagonalization')
    ax[1].plot(np.sqrt(data['counts']),.003*np.sqrt(data['counts'])*1000,'--',label='2 g0 sqrt(N)')
    ax[1].set(xlabel='sqrt(N), fixed g0=1.5 meV',ylabel='Resonant splitting (meV)',title='Collective scaling check');ax[1].legend()
    save(fig,'fig1_polaritonic_pauli_fierz_energy_levels.png')
    fig,ax=plt.subplots(figsize=(9,5.7),layout='constrained')
    mesh=ax.pcolormesh(data['angles'],data['energy'],data['map_T'],shading='auto',norm=PowerNorm(.4,vmin=0,vmax=1),cmap='magma')
    for j,label in enumerate(('LP','UP')):ax.plot(data['angles'],data['branches'][:,j],color='cyan',ls='--',lw=1,label=label)
    ax.plot(data['angles'],data['wc_angles'],'w:',label='Bare cavity');ax.axhline(cfg.vibration_ev,color='white',ls='-.',lw=1)
    ax.set(xlabel='External angle (degrees)',ylabel='Probe energy (eV)',ylim=(.135,.24),title='21C | Fabry-Perot angular transmission proxy')
    ax.legend(fontsize=8);fig.colorbar(mesh,ax=ax,label='Transmission (power-law color scale, gamma=0.4)')
    save(fig,'fig2_vacuum_rabi_dispersion_anticrossing.png')
    fig,ax=plt.subplots(1,2,figsize=(11,5),layout='constrained')
    for n,curve in data['rate_curves'].items():ax[0].plot(data['sweep']-cfg.vibration_ev,curve[:,0],label=f'N={n}, fixed collective g')
    ax[0].plot(data['sweep']-cfg.vibration_ev,data['rate_curves'][1][:,1],ls=':',label='N=1 classical limit')
    ax[0].axvline(0,color='black',ls='--',lw=.7);ax[0].axvline(cfg.barrier_frequency_ev-cfg.vibration_ev,color='grey',ls='-.',lw=.7)
    ax[0].set(xlabel='Detuning Ec - E0 (eV)',ylabel='k(cavity) / k(outside)',title='21B | Harmonic TST; no forced sharp resonance');ax[0].legend(fontsize=8)
    ax[1].loglog(data['dilution_N'],abs(data['dilution']-1),'o-',label='Computed local rate change')
    ax[1].loglog(data['dilution_N'],abs(data['dilution'][0]-1)/data['dilution_N'],'--',label='1/N guide')
    ax[1].set(xlabel='N at fixed collective g',ylabel='|rate ratio - 1|',title='Collective vs local');ax[1].legend(fontsize=8)
    save(fig,'fig3_cavity_reaction_rate_resonance_profile.png')
    fig,ax=plt.subplots(1,2,figsize=(11,5),layout='constrained')
    for g,t,r,a in data['ftir']:
        label=f'g/E0={g/cfg.vibration_ev:.3f}'
        ax[0].plot(data['energy']*EV_CM,a,label=label)
        ax[1].plot(data['energy']*EV_CM,t,label=label)
    ax[0].set(xlabel='Wavenumber (cm^-1)',ylabel='Absorptance (1-R-T)',title='21D | Weak-probe FTIR simulation')
    ax[1].set(xlabel='Wavenumber (cm^-1)',ylabel='Transmission',title='DSE and counterrotating response retained')
    for a in ax:a.legend(fontsize=8);a.axvline(cfg.vibration_ev*EV_CM,color='grey',ls=':',lw=.8)
    fig.suptitle('Ideal dark modes have zero direct oscillator strength; no fabricated dark absorption peak')
    save(fig,'fig4_simulated_ftir_polariton_spectra.png')


def write_reports(cfg,summary,root):
    refs='\n'.join(f'- [{t}]({u})' for t,u in REFS)
    table=f"""
## Computed results / 实际计算结果

| Quantity / 项目 | Value / 数值 |
|---|---:|
| Fundamental / 分子振动能量 (eV) | {cfg.vibration_ev} |
| Collective coupling / 集体耦合 (eV) | {cfg.collective_g_ev} |
| LP / UP (eV) | {summary['LP_ev']:.8f} / {summary['UP_ev']:.8f} |
| Rabi splitting / 劈裂 (meV) | {summary['rabi_mev']:.6f} |
| gamma + kappa (meV) | {summary['linewidth_sum_mev']:.6f} |
| Resolved strong coupling / 满足所用线宽判据 | {summary['resolved_strong_coupling']} |
| N=1 resonant quantum rate ratio / 共振处量子速率比 | {summary['rate_ratio_N1']:.8f} |
| N=64 resonant quantum rate ratio / 共振处量子速率比 | {summary['rate_ratio_N64']:.8f} |
| N=1 swept ratio range / 扫描范围 | {summary['rate_min']:.8f} to {summary['rate_max']:.8f} |
| Truncated Fock mode error / 小基组误差 (eV) | {summary['validation']['errors']['fock_low_mode_error_ev']:.4g} |

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

{refs}
"""
    en="""# Phase 21: cavity QED polaritonic chemistry

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
"""
    zh="""# 第二十一阶段：腔量子电动力学与极化激元化学

## 范围

使用指定参数的谐振 Pauli-Fierz 有效模型，研究集体光谱劈裂、暗态和局部反应瓶颈。
2.012 eV 对称四次势仅借用此前骨架反应能垒的量级，不是从 Phase 4/5 重新提取的真实势能面，
也没有拟合实际分子的偶极函数。平衡真空场计算没有外部泵浦；测量 FTIR 则需要弱探测光。
两者并不矛盾，本阶段没有实测光谱或化学反应数据。

## 21A：哈密顿量与数值基组

使用能量单位，H=sum E0(b_i†b_i+1/2)+Ec(a†a+1/2)
+g0(a+a†)sum(b_i+b_i†)+(g0²/Ec)[sum(b_i+b_i†)]²。
完整保留反旋转项和不同分子间的偶极自能交叉项。平方坐标算符直接投影，避免先截断再平方
造成最高能级边界项丢失。光子占据数为 0 至 3，振动占据数为 0 至 4。
N=2 的完整张量积对角化用于验证；大量相同谐振分子则精确约化为耦合强度 g=g0 sqrt(N)
的亮振子和 N-1 个暗振子。额外采用 12 个光子态、14 个振动态检验亮态基组收敛。
完整二次型正常模提供不依赖旋波近似的参照，并在共振条件、固定 g0 时检验 sqrt(N) 劈裂。
g0=lambda mu01 sqrt(Ec/2) 说明有效参数的物理含义；未指定真实腔体积和跃迁偶极，
因此不推断绝对真空场强，也不声称材料层面的从头算精度。

## 21B：局部反应与速率

势能面 V(q)=B[1-(q/q0)²]²，通过有效质量使势阱振动能量为 E0，
同一势对应的虚频绝对值为 E0/sqrt(2)。光子势为 1/2[Ec Q+sum a_i q_i]²，
a_i=2g0 sqrt(E0/Ec)，采用质量加权坐标。对 Q 极小化后恢复原始分子势能面；
省略自能才会产生虚假的静态势垒变化甚至失稳。导出的二维势能面是 N=1 的示意切片。

过渡态仅改变一个分子的曲率，其余保持势阱状态。剔除唯一不稳定模后计算谐振量子 TST：
k=(kBT/h) exp(-B/kBT) prod_min[2sinh(Ek/2kBT)]/prod_TS_stable[2sinh(Ek/2kBT)]。
零点能已包含在振子配分函数中，不另加一次零点能垒修正。
即使弛豫后的经典势垒不变，稳定模自由能和正常模分割面的变化仍可改变这个模型的速率。
其经典极限速率比等于耦合后不稳定频率与裸势垒频率之比。
扫描同时覆盖势阱实频与势垒虚频绝对值，不人为乘上共振峰函数。
此模型不保证在分子基频处出现尖锐抑制，也不包含隧穿、溶剂摩擦、真实反应轨迹或量子再穿越。
尤其在较高势垒下，绝对速率仅是谐振 TST 的示意估计。

局部稀释扫描固定集体 g，因此单分子 g0=g/sqrt(N)。光场耦合的是总偶极，
而单个断键坐标在亮态中的参与权重为 1/N。这解释本模型中的集体与局部差异，
不等于解决所有实验体系的集体反应问题。固定 g0 的劈裂扫描与固定 g 的稀释扫描
改变的是不同参数；扫腔频时也固定 g，不视为固定真实腔体积的实验。

## 21C/D：光谱与仪器代理模型

反演 K-E²I-i E diag(gamma,kappa) 得到被动双端口光学响应，保留自能及反旋转项，
包括 g/E0>0.1 的扫描。默认 Q=200，kappa=Ec/Q，分子线宽 gamma=3 meV。
使用对称端口，同时输出透射 T、反射 R 和真实吸收率 A=1-R-T，检验能量守恒。
劈裂大于所用线宽之和是本基准的分辨判据，不冒充唯一强耦合定义；
光谱阻尼为经验模型，不是超强耦合微观热浴主方程。

角度色散采用 Ec(theta)=Ec(0)/sqrt(1-sin²(theta)/n_eff²)，
属于 Fabry-Perot 腔角分辨代理模型，未计算真实 ATR 棱镜、多层薄膜或完整传输矩阵。
相同分子理想暗态没有直接光腔振子强度：在能级图计数暗态，但不伪造 FTIR 暗态峰。
要使暗态可见，需要无序、直接分子照明或其他对称性破缺机制及对应仪器建模。

## 解释与验证

本阶段说明边界条件在有效模型中可以改变混合正常模与平衡涨落，但光谱反交叉本身
不证明化学反应速率改变。测试覆盖 Fock 收敛、sqrt(N)、暗态简并、自能稳定性、
弛豫势能面抵消、经典行列式、零耦合速率恢复及光谱被动性。
实际残差和容差写入 JSON；这些检查验证程序与模型内部一致性，不能取代真实分子
非谐性、偶极标定、耗散反应动力学和重复实验。保留既有阶段，不宣称重新验证此前结论。
"""
    for lang,body in [('EN',en),('ZH',zh)]:
        (root/f'POLARITON_CHEMISTRY_REPORT_{lang}.md').write_text(body+table,encoding='utf-8',newline='\n')


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
    ap.add_argument('--self-test',action='store_true')
    args=ap.parse_args(argv)
    if args.self_test:
        result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests))
        return 0 if result.wasSuccessful() else 1
    cfg=Config();cfg.validate()
    validation,analytic,low,high,counts,gaps=validate(cfg)
    print('PF basis, DSE, collective scaling and null checks passed.',flush=True)
    energy=np.linspace(.06,.31,cfg.energy_points);angles=np.linspace(0,70,cfg.angle_points)
    wcangles=cfg.cavity_normal_ev/np.sqrt(1-(np.sin(np.deg2rad(angles))/cfg.effective_index)**2)
    mapT,_,_=transmission(cfg,energy[:,None],wcangles[None,:],cfg.collective_g_ev)
    branches=np.array([polaritons(cfg.vibration_ev,wc,cfg.collective_g_ev) for wc in wcangles])
    sweep=np.linspace(.06,.30,161)
    curves={n:np.array([rates(cfg,wc,n,cfg.collective_g_ev) for wc in sweep]) for n in (1,4,16,64)}
    dn=np.array([1,2,4,8,16,32,64,128])
    dilution=np.array([rates(cfg,cfg.vibration_ev,int(n),cfg.collective_g_ev)[0] for n in dn])
    ftir=[(g,*transmission(cfg,energy,cfg.vibration_ev,g)) for g in (.002,.009,.022,.04)]
    summary={'parameters':asdict(cfg),'validation':validation,'LP_ev':float(analytic[0]),'UP_ev':float(analytic[1]),
        'rabi_mev':float(1000*np.diff(analytic)[0]),'linewidth_sum_mev':1000*(cfg.gamma_ev+cfg.vibration_ev/cfg.quality_factor),
        'resolved_strong_coupling':bool(np.diff(analytic)[0]>cfg.gamma_ev+cfg.vibration_ev/cfg.quality_factor),
        'rate_ratio_N1':float(dilution[0]),'rate_ratio_N64':float(dilution[-2]),
        'rate_min':float(curves[1][:,0].min()),'rate_max':float(curves[1][:,0].max()),
        'rate_min_cavity_ev':float(sweep[np.argmin(curves[1][:,0])]),
        'provenance':'Assigned harmonic PF model, schematic quartic PES, weak-probe optical proxy; no experimental data or ab initio PES.',
        'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__}}
    root=args.output_dir.resolve();res=root/'results_phase21';figdir=root/'figures_phase21'
    res.mkdir(parents=True,exist_ok=True);figdir.mkdir(parents=True,exist_ok=True)
    def csv(name,arr,header):
        with (res/name).open('w',encoding='utf-8',newline='\n') as handle:np.savetxt(handle,arr,delimiter=',',header=header,comments='',fmt='%.12g')
    csv('fock_convergence.csv',np.column_stack((analytic,low,high)),'normal_mode_ev,Fock_4x5_ev,Fock_12x14_ev')
    csv('collective_scaling.csv',np.column_stack((counts,gaps,2*.0015*np.sqrt(counts))),'N,computed_split_ev,2g0sqrtN_ev')
    csv('angular_modes.csv',np.column_stack((angles,wcangles,branches)),'angle_deg,cavity_ev,LP_ev,UP_ev')
    ee,aa=np.meshgrid(energy,angles,indexing='ij')
    csv('angular_transmission.csv',np.column_stack((aa.ravel(),ee.ravel(),mapT.ravel())),'angle_deg,probe_ev,transmission')
    csv('reaction_rates.csv',np.column_stack((sweep,*[curves[n][:,j] for n in curves for j in range(3)])),
        'cavity_ev,'+','.join(f'N{n}_{k}' for n in curves for k in ('quantum_ratio','classical_ratio','quantum_rate_s-1')))
    csv('local_dilution.csv',np.column_stack((dn,dilution)),'N,quantum_rate_ratio_fixed_collective_g')
    csv('ftir.csv',np.column_stack((energy,*[arr for g,t,r,a in ftir for arr in (t,r,a)])),
        'probe_ev,'+','.join(f'g{g:.3f}_{s}' for g,*_ in ftir for s in ('T','R','A')))
    q=np.linspace(-1.3,1.3,101);Q=np.linspace(-18,18,121);qq,pp=np.meshgrid(q,Q)
    csv('reactive_pes_2d.csv',np.column_stack((qq.ravel(),pp.ravel(),pes(cfg,qq,pp,cfg.vibration_ev,cfg.collective_g_ev).ravel())),
        'dimensionless_reaction_q,photon_massweighted_Q,potential_ev')
    (res/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    data=dict(analytic=analytic,counts=counts,gaps=gaps,energy=energy,angles=angles,map_T=mapT,branches=branches,wc_angles=wcangles,sweep=sweep,rate_curves=curves,dilution_N=dn,dilution=dilution,ftir=ftir)
    figures(cfg,data,figdir);write_reports(cfg,summary,root)
    files=sorted([*res.glob('*.csv'),res/'summary.json',*figdir.glob('*.png'),*root.glob('POLARITON_CHEMISTRY_REPORT_*.md')])
    hashes={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    hashes['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (res/'sha256.json').write_text(json.dumps(hashes,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:summary[k] for k in ('rabi_mev','rate_ratio_N1','rate_ratio_N64','rate_min_cavity_ev')},indent=2))
    return 0

if __name__=='__main__': raise SystemExit(main())
