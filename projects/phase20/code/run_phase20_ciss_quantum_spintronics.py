#!/usr/bin/env python3
"""Phase 20: auditable effective-model CISS transport and illustrative twins.

Requires Python >=3.10, numpy >=1.23, scipy >=1.9, matplotlib >=3.6.
Run: python run_phase20_ciss_quantum_spintronics.py --output-dir .
Test: python run_phase20_ciss_quantum_spintronics.py --self-test

No target polarization or catalytic benefit is imposed. A coherent single
orbital chain with nonmagnetic contacts is a required null control. A second
orbital path makes the device block-tridiagonal with 4x4 site blocks and allows
spin/orbital entanglement. This is a phenomenological Pauli/Rashba tight-binding
model, not an ab initio Dirac calculation or experimental validation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import unittest
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np
import scipy
from scipy import sparse
from scipy.sparse.linalg import splu
from scipy.special import expit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

I2 = np.eye(2, dtype=complex)
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
KB_EV = 8.617333262145e-5
G0_PER_SPIN = 3.874045864931824e-5  # e^2/h, A/V; energy integrals are in eV
P_FLOOR = 1e-10
REFERENCES = (
    ("Varela et al., one-dimensional Rashba transport and the coherent null control",
     "https://arxiv.org/abs/2301.02156"),
    ("Evers et al., Theory of Chirality Induced Spin Selectivity: Progress and Challenges",
     "https://arxiv.org/abs/2108.09998"),
    ("Spin polarization and OER in a monolayer chiral COF model catalyst (2025)",
     "https://doi.org/10.1021/jacs.5c09729"),
)


@dataclass(frozen=True)
class Config:
    sites: int = 24
    orbitals: int = 2
    turns: float = 2.5
    radius_angstrom: float = 4.0
    pitch_per_turn_angstrom: float = 5.4
    hopping_ev: float = 0.55
    rung_ev: float = 0.28
    orbital_offset_ev: float = 0.16
    soc_ev: float = 0.03
    lead_hopping_ev: float = 1.5
    contact_ev: float = 0.65
    lead_exchange_ev: float = 0.0
    eta_ev: float = 1e-10  # numerical retarded regulator, NOT a dephasing bath
    energy_points: int = 401
    temperature_k: float = 300.0
    tip_polarization: float = 0.4
    oer_alpha: float = 0.5
    oer_j0_ma_cm2: float = 1e-3
    oer_barrier_gain_ev: float = 0.10  # assigned sensitivity, not fitted chemistry
    mchd_coefficient_per_t: float = 0.002  # assigned optical sensitivity
    seed: int = 2026

    def validate(self):
        if self.sites < 4 or self.orbitals not in (1, 2) or self.energy_points < 21:
            raise ValueError("sites >=4, orbitals in {1,2}, energy-points >=21 required")
        if not all(np.isfinite(v) for v in asdict(self).values()):
            raise ValueError("All parameters must be finite")
        if min(self.hopping_ev, self.lead_hopping_ev, self.contact_ev,
               self.radius_angstrom, self.pitch_per_turn_angstrom,
               self.turns, self.temperature_k, self.eta_ev, self.oer_j0_ma_cm2) <= 0:
            raise ValueError("Geometry, hopping, contact, temperature and regulator must be positive")
        if self.soc_ev < 0 or not 0 < self.oer_alpha <= 1:
            raise ValueError("SOC >=0 and 0<alpha<=1 required")
        if abs(self.tip_polarization) > 1 or self.oer_barrier_gain_ev < 0:
            raise ValueError("Unphysical tip polarization or OER gain")
        if self.lead_exchange_ev != 0:
            raise ValueError("This benchmark's mirror/null tests require nonmagnetic leads")


def helix(cfg: Config, chirality: int):
    """Reflection y -> -y maps + to - at fixed transport direction z.

    Coordinates are Angstrom. Radial vectors represent normalized grad(V);
    cross(grad(V), tangent) is an axial SOC direction, not a magnetic field in T.
    """
    if chirality not in (-1, 1):
        raise ValueError("chirality must be +1 or -1")
    phi = np.linspace(0, 2*np.pi*cfg.turns, cfg.sites)
    xyz = np.zeros((cfg.sites, cfg.orbitals, 3))
    for a in range(cfg.orbitals):
        angle = chirality*(phi + a*1.15)
        radius = cfg.radius_angstrom*(1 + 0.15*a)
        xyz[:, a] = np.column_stack((radius*np.cos(angle), radius*np.sin(angle),
                                    cfg.pitch_per_turn_angstrom*phi/(2*np.pi)))
    tangent = np.diff(xyz, axis=0)
    tangent /= np.linalg.norm(tangent, axis=2, keepdims=True)
    radial = (xyz[:-1] + xyz[1:])/2
    radial[:, :, 2] = 0
    radial /= np.linalg.norm(radial, axis=2, keepdims=True)
    soc_axis = np.cross(radial, tangent)
    return xyz, radial, tangent, soc_axis


def hamiltonian(cfg: Config, chirality: int):
    """Hermitian, time-reversal invariant Pauli Hamiltonian in CSC format.

    Basis (site, orbital, spin-z). H[i,i+1] = -t I + i lambda n.sigma.
    The bond SOC sign is a convention fixed by this directed hopping.
    Scalar longitudinal modulation and rung coupling are mirror-even.
    """
    _, _, _, axis = helix(cfg, chirality)
    b = 2*cfg.orbitals
    h = sparse.lil_matrix((b*cfg.sites, b*cfg.sites), dtype=complex)
    for i in range(cfg.sites):
        for a in range(cfg.orbitals):
            j = i*b + 2*a
            onsite = (a-0.5)*cfg.orbital_offset_ev + 0.08*np.cos(2*np.pi*i/(cfg.sites-1))
            h[j:j+2, j:j+2] = onsite*I2
            if i < cfg.sites-1:
                bond = -cfg.hopping_ev*(1-0.12*a)*I2 + 1j*cfg.soc_ev*np.einsum('k,kij->ij', axis[i, a], PAULI)
                h[j:j+2, j+b:j+b+2] = bond
                h[j+b:j+b+2, j:j+2] = bond.conj().T
        if cfg.orbitals == 2:
            j = i*b
            h[j:j+2, j+2:j+4] = -cfg.rung_ev*I2
            h[j+2:j+4, j:j+2] = -cfg.rung_ev*I2
    return h.tocsc()


def surface_green(energy, hopping):
    """Exact retarded surface GF of a semi-infinite chain, units 1/eV.

    Select decaying root outside the band and negative-imaginary root inside;
    taking the principal sqrt alone on the negative-energy side is incorrect.
    """
    e = float(energy)
    if abs(e) <= 2*hopping:
        return (e - 1j*np.sqrt(max(0.0, 4*hopping*hopping-e*e)))/(2*hopping*hopping)
    return (e - np.sign(e)*np.sqrt(e*e-4*hopping*hopping))/(2*hopping*hopping)


def scattering(cfg, h, energy):
    """Sparse LU selected-column solve; no explicit dense matrix inverse.

    Transmission rows = outgoing RIGHT spin/orbital, columns = incident LEFT.
    This convention fixes the physical meaning of T_up/down and P at the drain.
    """
    b = 2*cfg.orbitals
    dim = h.shape[0]
    sigma = cfg.contact_ev**2*surface_green(energy, cfg.lead_hopping_ev)
    gamma = max(0.0, -2*sigma.imag)
    d = np.zeros(dim, complex)
    d[:b] += sigma
    d[-b:] += sigma
    a = (energy+1j*cfg.eta_ev)*sparse.eye(dim, format='csc') - h - sparse.diags(d, format='csc')
    rhs = np.zeros((dim, b), complex)
    rhs[:b] = np.eye(b)
    sol = splu(a).solve(rhs)
    residual = float(np.max(np.abs(a@sol-rhs)))
    t = -1j*gamma*sol[-b:]
    r = np.eye(b)-1j*gamma*sol[:b]
    rho = t@t.conj().T
    total = float(np.trace(rho).real)
    components = np.array([np.trace(np.kron(np.eye(cfg.orbitals), s)@rho).real for s in PAULI])
    p = components/total if total > P_FLOOR else np.full(3, np.nan)
    spin = np.abs(t.reshape(cfg.orbitals, 2, cfg.orbitals, 2))**2
    # Sum orbital indices, preserve outgoing/incident spin indices.
    spin_matrix = spin.sum(axis=(0, 2))
    deficit = float(b-np.trace(r@r.conj().T).real-total)
    return dict(total=total, up=float(spin_matrix[0].sum()), down=float(spin_matrix[1].sum()),
                polarization=p, spin_matrix=spin_matrix, residual=residual, deficit=deficit)


def spectrum(cfg, chirality, energies):
    h = hamiltonian(cfg, chirality)
    rows = [scattering(cfg, h, e) for e in energies]
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


def integrate(y, x):
    # Support NumPy 1.23 through 2.x without deprecated np.trapz on newer releases.
    if hasattr(np, 'trapezoid'):
        return np.trapezoid(y, x=x, axis=0)
    return np.trapz(y, x=x, axis=0)


def transport_average(cfg, energies, spec):
    f = expit(-energies/(KB_EV*cfg.temperature_k))
    w = f*(1-f)/(KB_EV*cfg.temperature_k)
    numerator = integrate((spec['up']-spec['down'])*w, energies)
    denom = integrate(spec['total']*w, energies)
    if denom <= P_FLOOR:
        raise ValueError('Insufficient thermally averaged transmission for OER/twin coupling')
    return float(numerator/denom)


def oer(cfg, polarization):
    """Assigned kinetic sensitivity. Not a calculation of OER intermediates.

    j0(P)=j0 exp(gain P^2/kBT); alpha unchanged, so the Tafel slope is unchanged.
    Delta eta = gain P^2/alpha (volts). Spin filtering alone proves neither a
    180 mV shift nor a slope change or elimination of singlet oxygen.
    """
    b = np.log(10)*KB_EV*cfg.temperature_k/cfg.oer_alpha
    shift = cfg.oer_barrier_gain_ev*polarization**2/cfg.oer_alpha
    logj = np.linspace(-3, 2, 251)
    base = b*(logj-np.log10(cfg.oer_j0_ma_cm2))
    return dict(logj=logj, achiral_eta=base, chiral_eta=base-shift,
                tafel_v_dec=float(b), shift_v=float(shift), p_used=float(polarization))


def twins(cfg, energies, spec):
    """Weak-analyzer Landauer SP-STM and independent phenomenological MChD.

    Tip is a downstream spin analyzer, not a self-consistent magnetic junction.
    Optical coefficients are assigned, NOT inferred from electronic transmission.
    """
    voltage = np.linspace(-0.30, 0.30, 121)
    currents = []
    for v in voltage:
        window = expit((v/2-energies)/(KB_EV*cfg.temperature_k))-expit((-v/2-energies)/(KB_EV*cfg.temperature_k))
        spin_part = cfg.tip_polarization*(spec['up']-spec['down'])
        currents.append([G0_PER_SPIN*integrate((spec['total']+sign*spin_part)*window/2, energies) for sign in (1, -1)])
    currents = np.array(currents)
    denominator = currents.sum(axis=1)
    asym = np.divide(currents[:,0]-currents[:,1], denominator,
                     out=np.full(len(voltage), np.nan), where=np.abs(denominator)>1e-18)
    wavelength = np.linspace(250, 550, 301)
    absorption = 0.06+0.7*np.exp(-0.5*((wavelength-360)/35)**2)
    # B=+/-1 T and k along +z. Chirality chi multiplies this signed template.
    delta = 2*cfg.mchd_coefficient_per_t*absorption
    rng = np.random.default_rng(cfg.seed)
    measured = delta+rng.normal(0, 0.0001, len(delta))
    return dict(voltage=voltage, currents=currents, asymmetry=asym,
                wavelength=wavelength, absorption=absorption, delta=delta,
                noisy_delta=measured)


def audit(cfg):
    """Numerical/physical invariants, independent dense reference, grid stability."""
    probe = np.array([-0.73, -0.29, 0.07, 0.31, 0.68])
    d = spectrum(cfg, 1, probe)
    l = spectrum(cfg, -1, probe)
    zero = spectrum(replace(cfg, soc_ev=0), 1, probe)
    single = spectrum(replace(cfg, orbitals=1), 1, probe)
    h = hamiltonian(cfg, 1).toarray()
    tr = np.kron(np.eye(cfg.sites*cfg.orbitals), 1j*PAULI[1])
    mirror = np.kron(np.eye(cfg.sites*cfg.orbitals), PAULI[1])
    hl = hamiltonian(cfg, -1).toarray()
    small = replace(cfg, sites=6)
    hs = hamiltonian(small, 1).toarray()
    b = 2*small.orbitals
    sig = cfg.contact_ev**2*surface_green(0.13,cfg.lead_hopping_ev)
    a = (0.13+1j*cfg.eta_ev)*np.eye(len(hs))-hs
    a[:b,:b] -= sig*np.eye(b); a[-b:,-b:] -= sig*np.eye(b)
    dense = np.linalg.inv(a)[-b:,:b]*(-2*sig.imag)
    dense_total = float(np.sum(np.abs(dense)**2))
    sparse_total = scattering(small, sparse.csc_matrix(hs), 0.13)['total']
    fine = spectrum(replace(cfg, eta_ev=cfg.eta_ev/10), 1, probe)
    metrics = {
        'hermiticity_max': float(np.max(np.abs(h-h.conj().T))),
        'time_reversal_max': float(np.max(np.abs(tr@h.conj()@tr.conj().T-h))),
        'mirror_hamiltonian_max': float(np.max(np.abs(mirror@h@mirror-hl))),
        'mirror_total_max': float(np.max(np.abs(d['total']-l['total']))),
        'mirror_vector_max': float(np.nanmax(np.abs(l['polarization']-d['polarization']*[-1,1,-1]))),
        'zero_soc_p_max': float(np.nanmax(np.abs(zero['polarization']))),
        'single_channel_p_max': float(np.nanmax(np.abs(single['polarization']))),
        'sparse_dense_transmission_max': abs(dense_total-sparse_total),
        'linear_solve_residual_max': float(np.max(d['residual'])),
        'flux_deficit_max': float(np.max(np.abs(d['deficit']))),
        'eta_transmission_change_max': float(np.max(np.abs(d['total']-fine['total']))),
        'spin_sum_error_max': float(np.max(np.abs(d['up']+d['down']-d['total']))),
    }
    tolerance = {k: 1e-7 for k in metrics}
    tolerance['single_channel_p_max'] = 2e-6
    tolerance['flux_deficit_max'] = 2e-5
    tolerance['eta_transmission_change_max'] = 2e-5
    failures = {k:v for k,v in metrics.items() if not np.isfinite(v) or v>tolerance[k]}
    if failures:
        raise RuntimeError(f'Validation failed: {failures}')
    if np.min(d['spin_matrix']) < -1e-12 or np.max(d['total']) > b+1e-7:
        raise RuntimeError('Transmission outside channel bounds')
    if np.nanmax(np.linalg.norm(d['polarization'], axis=1)) > 1+1e-7:
        raise RuntimeError('Polarization norm exceeds 1')
    return dict(status='passed', errors=metrics, tolerances=tolerance,
                checks='Hermiticity, TRS, mirror, SOC=0, single path, dense reference, flux, eta, positivity, channel and spin bounds')


class PhysicsTests(unittest.TestCase):
    def test_invariants(self):
        for orbitals in (1,2):
            for soc in (0,0.01,0.05):
                self.assertEqual(audit(Config(sites=8, orbitals=orbitals, soc_ev=soc))['status'], 'passed')

    def test_surface_branch(self):
        for e in (-10,-3,-1,0,1,3,10):
            g=surface_green(e,1)
            self.assertLessEqual(g.imag,0)
            self.assertAlmostEqual(abs(g*(e-g)-1),0,places=11)
        self.assertLess(abs(surface_green(100,1)-0.01),2e-6)

    def test_closed_band(self):
        cfg=Config(sites=8)
        s=scattering(cfg,hamiltonian(cfg,1),5)
        self.assertEqual(s['total'],0)
        self.assertTrue(np.isnan(s['polarization']).all())

    def test_perfect_uniform_wire(self):
        cfg=Config(sites=8,orbitals=1,lead_hopping_ev=1,contact_ev=1)
        chain=sparse.diags([-np.ones(7),-np.ones(7)],[-1,1],shape=(8,8),format='csc')
        h=sparse.kron(chain,I2,format='csc')
        for e in (-1.7,0,1.7):
            s=scattering(cfg,h,e)
            self.assertAlmostEqual(s['total'],2,places=7)
            self.assertLess(abs(s['deficit']),1e-7)

    def test_oer_null_and_evenness(self):
        cfg=Config()
        self.assertEqual(oer(cfg,0)['shift_v'],0)
        self.assertEqual(oer(cfg,0.3)['shift_v'],oer(cfg,-0.3)['shift_v'])
        q=oer(cfg,0.3)
        slope=np.polyfit(q['logj'],q['chiral_eta'],1)[0]
        self.assertAlmostEqual(slope,q['tafel_v_dec'],places=12)

    def test_twins_zero_and_reproducible(self):
        cfg=Config(sites=8,soc_ev=0)
        e=np.linspace(-1,1,101)
        s=spectrum(cfg,1,e)
        a=twins(cfg,e,s); b=twins(cfg,e,s)
        self.assertTrue(np.array_equal(a['noisy_delta'],b['noisy_delta']))
        self.assertLess(np.nanmax(np.abs(a['asymmetry'])),1e-12)
        self.assertTrue(np.isnan(a['asymmetry'][60]))
        self.assertTrue(np.allclose(a['currents'],-a['currents'][::-1],atol=1e-15))

    def test_validation(self):
        for cfg in (Config(sites=0),Config(soc_ev=-1),Config(temperature_k=0),Config(lead_exchange_ev=.1)):
            with self.assertRaises(ValueError): cfg.validate()


def figure_outputs(cfg, energy, d, l, single, kinetics, optical, directory):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300})
    def save(fig,name):
        fig.savefig(directory/name,dpi=300,facecolor='white');plt.close(fig)
    xyz,radial,tangent,axis=helix(cfg,1)
    fig=plt.figure(figsize=(9,6),layout='constrained')
    ax=fig.add_subplot(111,projection='3d')
    for a in range(cfg.orbitals): ax.plot(*xyz[:,a].T,label=f'Orbital path {a+1}')
    loc=(xyz[:-1,0]+xyz[1:,0])/2
    for v,color,label in ((radial[:,0],'tab:orange','grad V direction'),(axis[:,0],'tab:green','SOC effective-field direction'),(tangent[:,0],'tab:red','Electron tangent')):
        ax.quiver(*loc[::3].T,*v[::3].T,length=1.7,color=color,label=label)
    ax.set(xlabel='x (Angstrom)',ylabel='y (Angstrom)',zlabel='z (Angstrom)',title='20A | Helical geometry and normalized SOC vectors')
    ax.legend(loc='upper left',fontsize=8)
    save(fig,'fig1_helical_spin_orbit_vector.png')
    fig,axs=plt.subplots(3,1,figsize=(9,9),sharex=True,layout='constrained')
    axs[0].plot(energy,d['up'],label='Right spin up');axs[0].plot(energy,d['down'],'--',label='Right spin down')
    axs[0].set(ylabel='Transmission',title='20B | Effective two-path coherent NEGF model');axs[0].legend()
    for s,label,style in ((d,'D (+)','-'),(l,'L (-)','--')):
        axs[1].plot(energy,100*s['polarization'][:,2],style,label=label)
    axs[1].set(ylabel='Drain Pz (%)');axs[1].legend();axs[1].axhline(0,color='black',lw=.5)
    axs[2].plot(energy,100*single['polarization'][:,2],color='black',label='One orbital: numerical null')
    axs[2].set(xlabel='E - EF (eV)',ylabel='Control Pz (%)',ylim=(-.001,.001));axs[2].legend()
    save(fig,'fig2_spin_transmission_polarization.png')
    fig,axs=plt.subplots(1,2,figsize=(10,4.8),layout='constrained')
    axs[0].plot(kinetics['achiral_eta'],kinetics['logj'],label='P=0 reference')
    axs[0].plot(kinetics['chiral_eta'],kinetics['logj'],'--',label='Computed thermal Pz')
    axs[0].set(xlabel='Overpotential (V)',ylabel='log10 J (mA / cm2)',title='20C | Assigned kinetic sensitivity');axs[0].legend()
    p=np.linspace(0,1,201)
    axs[1].plot(p,1000*cfg.oer_barrier_gain_ev*p*p/cfg.oer_alpha)
    axs[1].scatter([abs(kinetics['p_used'])],[1000*kinetics['shift_v']],color='black',label='Computed Pz')
    axs[1].set(xlabel='Assumed |Pz|',ylabel='Conditional overpotential shift (mV)',title='Not a predicted OER mechanism');axs[1].legend()
    fig.suptitle(f"Same Tafel slope: {1000*kinetics['tafel_v_dec']:.2f} mV/dec; no imposed 180 mV gain")
    save(fig,'fig3_oer_overpotential_tafel.png')
    fig,axs=plt.subplots(2,2,figsize=(10,8),layout='constrained')
    for i,label in enumerate(('Tip +z','Tip -z')):axs[0,0].plot(optical['voltage'],1e6*optical['currents'][:,i],label=label)
    axs[0,0].set(xlabel='Bias (V)',ylabel='Current (microampere)',title='20D | Weak-analyzer SP-STM');axs[0,0].legend()
    axs[0,1].plot(optical['voltage'],optical['asymmetry'])
    axs[0,1].set(xlabel='Bias (V)',ylabel='Current asymmetry',title='Undefined at zero current')
    w=optical['wavelength'];a=optical['absorption'];delta=optical['delta']
    axs[1,0].plot(w,a+delta/2,label='B=+1 T');axs[1,0].plot(w,a-delta/2,'--',label='B=-1 T')
    axs[1,0].set(xlabel='Wavelength (nm)',ylabel='Absorbance (assigned)',title='Phenomenological MChD');axs[1,0].legend()
    axs[1,1].plot(w,delta,label='D ideal');axs[1,1].plot(w,-delta,'--',label='L ideal')
    axs[1,1].scatter(w[::6],optical['noisy_delta'][::6],s=6,color='black',alpha=.5,label='Synthetic noise')
    axs[1,1].set(xlabel='Wavelength (nm)',ylabel='A(+B) - A(-B)',title='Odd in chirality and k dot B');axs[1,1].legend()
    save(fig,'fig4_mchd_spstm_spectroscopic_twin.png')


def reports(cfg, summary, root):
    refs='\n'.join(f'- [{title}]({url})' for title,url in REFERENCES)
    common=f"""
## Numerical record / 数值记录

| Quantity / 项目 | Value / 数值 |
|---|---:|
| Sites / 格点数 | {cfg.sites} |
| Orbital paths / 轨道通道数 | {cfg.orbitals} |
| SOC (meV) | {1000*cfg.soc_ev:.6g} |
| Energy window / 能窗 | -1 to +1 eV |
| Energy points / 能量点数 | {cfg.energy_points} |
| Sampled maximum absolute Pz (%) / 采样点最大纵向极化绝对值 | {summary['peak_abs_pz_percent']:.8g} |
| Thermal Pz (%) / 热加权纵向极化 | {100*summary['thermal_pz']:.8g} |
| Single-channel max absolute P (%) / 单通道极化残差 | {summary['single_channel_max_p_percent']:.8g} |
| Conditional OER shift (mV) / 条件性过电位变化 | {summary['oer_shift_mv']:.8g} |
| Tafel slope, both (mV/dec) / 两条曲线共同斜率 | {summary['tafel_mv_dec']:.8g} |
| Coarse/fine thermal Pz difference / 能网格加密差值 | {summary['energy_grid_pz_difference']:.8g} |

## Reproduction / 复现

```bash
python -m pip install -r requirements_phase20.txt
python run_phase20_ciss_quantum_spintronics.py --self-test
python run_phase20_ciss_quantum_spintronics.py --output-dir .
```

SOC scan / SOC 参数扫描：`--soc-mev 10`, `--soc-mev 30`, `--soc-mev 50`.
Use separate `--output-dir` paths for scans / 扫描时使用独立输出目录。
Only numpy, scipy, matplotlib are required / 仅需这三个第三方依赖。
Run from any directory; default output is beside the script / 默认输出位于脚本目录。

## Artifacts / 产物

- `results_phase20/summary.json`: parameters, calculated values, environment and validation tolerances.
- `results_phase20/transport.csv`: D/L spin transmission, polarization vectors and one-path control.
- `results_phase20/helix.csv`: geometry and normalized gradient, tangent and SOC directions.
- `results_phase20/oer.csv`: conditional kinetic curves, not measured currents.
- `results_phase20/oer_sensitivity.csv`: assigned polarization scan for the second OER panel.
- `results_phase20/spstm.csv`, `results_phase20/mchd.csv`: simulated analytical signals.
- `results_phase20/sha256.json`: generated-file checksums and script checksum.
- Four PNG figures in `figures_phase20/`, 300 DPI, generated from the exported arrays.

![Helix](figures_phase20/fig1_helical_spin_orbit_vector.png)
![Transport](figures_phase20/fig2_spin_transmission_polarization.png)
![OER](figures_phase20/fig3_oer_overpotential_tafel.png)
![Twins](figures_phase20/fig4_mchd_spstm_spectroscopic_twin.png)

## References / 参考文献

{refs}
"""
    en="""# Phase 20: CISS quantum transport, catalytic sensitivity and analytical twins

## Scope and result interpretation

This is an effective-model computational benchmark. It is not experimental evidence,
an atomistically parameterized material, or a prediction of near-perfect filtering.
All reported values below are generated by the script. A single coherent orbital
path with nonmagnetic leads gives the required null result. The two-path model is
an explicit extension, not a hidden spin-dependent loss term.

## 20A: Hamiltonian and relativistic interpretation

Each helix is r(phi)=(R cos(chi phi), R sin(chi phi), p phi), with chi=+/-1
and p=pitch/(2 pi). Each longitudinal slice has two spinful orbital paths.
Nearest-slice hopping is -t I + i lambda n.sigma, with n=radial cross tangent;
the reverse block is its Hermitian conjugate. Real rung hopping closes loops.
The resulting sparse matrix is block-tridiagonal (4 by 4 blocks in the main run).
The radial electrostatic-gradient direction and directed tangent are explicitly
computed from geometry. The SOC axis is an axial vector; it is not a measured field.

This implements a Rashba-like Pauli SOC effective hopping. Atomic/interface SOC
and electrostatic gradients are subsumed in the assigned 10-50 meV coupling.
It does not derive that coupling from a scalar potential, effective mass or c,
and does not separately model bulk Dresselhaus SOC. It is not a four-component
relativistic quantum chemistry calculation. Angstrom is used for geometry and eV
for all Hamiltonian entries. Geometry changes at fixed hopping are model changes,
not spatial convergence of a continuous Schrödinger equation.

## 20B: Contacts, NEGF, spin and symmetry

Contacts are identical semi-infinite nonmagnetic chains, one per orbital/spin.
Their exact surface Green function uses the retarded decaying branch;
Sigma=t_contact^2 g_surface and Gamma=-2 Im Sigma. Sparse LU solves only contact
columns of (E+i eta-H-Sigma_L-Sigma_R)G=I. The regulator eta is numerical and
does not stand for decoherence or an unaccounted absorbing reservoir.
The outgoing-right transmission amplitude is sqrt(Gamma_R) G_RL sqrt(Gamma_L).
T_up/down sum over all incident spins and orbital channels. The full three-vector
P is Tr[(I_orb tensor sigma)tt†]/Tr[tt†]. Transmission is dimensionless.
Polarization below a transmission floor is undefined, stored as NaN in CSV.

Single-path time-reversal symmetry and unitarity prevent filtering of unpolarized
incident electrons; spin precession alone is insufficient. Two orbital paths
permit orbital/spin entanglement and interference, without guaranteeing large P.
The model is recomputed independently for each enantiomer. Under reflection
y -> -y, an axial vector transforms as (Px,Py,Pz) -> (-Px,Py,-Pz).
Thus longitudinal Pz reverses, whereas reversal of the entire vector is NOT a
general symmetry theorem. H_L=U H_D U† with U=I tensor sigma_y provides the
exact matrix identity for this model. Contacts and current direction are held fixed.

## 20C: OER is a conditional sensitivity study

The thermal Pz averages spin current divided by charge current using -df/dE,
not an average of polarization ratios. Assigned exchange current j0 and transfer
coefficient alpha define a one-electron effective anodic Tafel branch:
log10 j = log10 j0 + [alpha eta + gain Pz^2]/(kBT ln 10).
The gain is an assigned barrier-response coefficient (eV); no OER energetics were
computed. Delta eta=gain Pz^2/alpha is therefore conditional. The transfer
coefficient stays fixed and both slopes are identical. No 180 mV benefit or
slope reduction is imposed. The sensitivity panel shows what stronger assigned
polarization would do without presenting it as a measured or NEGF-derived result.
Very low-overpotential points are an extrapolation of the anodic Tafel branch,
not a full Butler-Volmer net current. No equilibrium water-splitting voltage is
added because the horizontal axis is overpotential, not absolute electrode voltage.

OER has several intermediates and catalyst-dependent steps. Triplet oxygen does
not justify assigning all losses to spin flips. Neither singlet-oxygen prohibition
nor elimination of radical damage follows from this model. Experimental CISS/OER
reports motivate research, but do not calibrate our arbitrary scaffold.

## 20D: analytical simulations and biological outlook

SP-STM integrates finite-temperature Landauer current over the exported energy
window with a weak tip analyzer, (T +/- p_tip(T_up-T_down))/2. The Hamiltonian
does not change under tip reversal; finite-bias electrostatic feedback and magnetic
contact backaction are omitted. A_STM is undefined at zero current. MChD is a
separate assigned optical response, A_chi(B)=A0(1+chi c k_hat.B). Reversing chirality,
field or propagation direction reverses the ideal difference. Added Gaussian noise
has a fixed seed and is explicitly synthetic. Optical transitions are not computed
from the transport Hamiltonian, and these traces are not validation measurements.

Chiral biomolecules motivate questions about spin in electron transfer, but this
script includes no mitochondrial proteins, oxygen intermediates, ROS dynamics or
biological measurements. Protection against oxidative damage remains outside its
evidence. Spin-controlled catalysis is a research direction requiring independent
electrochemical, structural and spin measurements, not a demonstrated application here.

## Validation and limits

The run stops on failed Hermiticity, time reversal, mirror covariance, positivity,
channel bounds, sparse/dense agreement, current conservation or null controls.
It checks sensitivity to a tenfold smaller retarded regulator and computes a
refined energy-grid thermal average. The JSON lists actual errors and tolerances.
No tests establish a universal CISS mechanism, material accuracy or production
readiness. Increasing site count changes the model at fixed hopping; fitting to
real materials, inelastic transport, interface chemistry and OER microkinetics
are future work. Earlier phases are not revalidated by this Phase 20 addition.
"""
    zh="""# 第二十阶段：CISS 量子输运、催化敏感性与分析信号仿真

## 范围与结论边界

本阶段是有效模型计算基准。下表数值由脚本实际计算；不把近乎完美的自旋过滤、
180 mV 过电位下降或真实实验验证作为预设结论。单轨道相干输运是必须通过的零极化对照；
双轨道模型是明确增加的物理自由度，不通过隐蔽的自旋吸收项制造极化。

## 20A：螺旋结构与自旋轨道耦合

使用 r(phi)=(R cos(chi phi), R sin(chi phi), p phi)，chi=±1，p=每圈螺距/(2 pi)。
每个纵向格点含两个带自旋的轨道路径；相邻格点跃迁为 -t I+i lambda n·sigma，
n 由径向静电势梯度方向与电子运动切向叉乘得到，反向块取厄米共轭，实数横向跃迁连接两条路径。
主模型是 4×4 分块的稀疏块三对角矩阵。几何单位为埃，能量统一使用 eV。

这是 Pauli/Rashba 型有效紧束缚模型，10–50 meV 耦合为指定参数，代表未显式计算的
原子及界面 SOC。没有从实际势场、有效质量或光速推导该数值，也没有独立加入体相
Dresselhaus 耦合。图中的有效场箭头仅表示方向，不是以特斯拉计量的真实磁场。
本模型不属于四分量相对论从头算；固定跃迁时改变格点数也不等同于连续方程网格收敛。

## 20B：NEGF、自旋透射与对称性

左右电极为非磁性半无限链，各轨道和自旋通道相同。用满足推迟边界条件的表面格林函数
计算自能 Sigma=t_contact² g 和展宽 Gamma=-2 Im Sigma；稀疏 LU 只求接触所需列，
不显式求整个稠密逆矩阵。微小 eta 是数值正则项，不代表退相干、非弹性散射或漏电浴。

透射矩阵采用右端出射、左端入射的指标约定；T_up/down 对全部入射自旋及轨道求和。
由出射密度矩阵 tt† 计算三个方向的极化，极低透射时将极化标为未定义。
单轨道、非磁性、双端相干系统在时间反演和幺正约束下不能仅靠自旋进动产生净过滤。
增加轨道路径允许自旋与轨道纠缠及路径干涉，但不保证高极化。

分别重建 D/L 哈密顿量并独立计算。采用 y→-y 的镜像时，自旋作为轴矢量满足
(Px,Py,Pz)→(-Px,Py,-Pz)：沿电流方向的 Pz 反号，整个矢量并非无条件反号。
脚本验证 H_L=U H_D U†（U=I⊗sigma_y）、纵向极化反号和总透射不变。
电极和输运方向固定，因此不将翻转电流或电极磁化混同于手性翻转。

## 20C：OER 条件性动力学敏感性

先用费米分布导数加权，求热自旋流与热电荷流之比。再指定交换电流、转移系数 alpha
和能垒响应参数 gain，使用 log10 j=log10 j0+[alpha eta+gain Pz²]/(kBT ln10)。
由此得到的过电位变化 gain Pz²/alpha 是模型条件性结果，不是计算出的 OER 反应能垒。
alpha 保持不变，两条 Tafel 曲线斜率相同；脚本不强制生成斜率下降或至少 180 mV 的改善。
附图另给假定极化大小的敏感性扫描，并与实际输运计算点明确区分。
低过电位部分只是阳极 Tafel 支的外推，不是完整 Butler–Volmer 净电流；横轴为过电位，
不冒充含平衡电位的绝对电极电压。

OER 涉及多个催化剂相关中间体，三重态氧的存在不能证明所有损失都源于自旋翻转。
本阶段未计算吸附自由能、过渡态或 ROS 生成，也不能证明单重态氧被完全禁止。
已有实验研究支持继续探索自旋与 OER 的联系，但并不为这里的任意螺旋模型提供校准。

## 20D：分析信号与生物学展望

SP-STM 使用弱磁针尖分析器，将 (T±p_tip(T_up-T_down))/2 代入有限温度 Landauer 积分。
针尖反转不反馈修改哈密顿量，未包含有限偏压下电势重排或磁性接触反作用；
零电流处电流不对称度为未定义，而不是人为置为可观测信号。
MChD 是独立的经验吸收模型 A_chi(B)=A0(1+chi c k_hat·B)，系数为指定值，
理想信号对手性、磁场、光传播方向分别反号。高斯噪声固定随机种子并标为合成噪声。
光学跃迁没有从输运哈密顿量求出，所有图均为计算或仿真信号，不能视为实测验证。

生物分子的手性使电子转移中的自旋效应值得研究，但本脚本没有线粒体蛋白、氧中间体、
ROS 动力学或生物实验。不能据此声称生命利用 CISS 防止氧化损伤。
面向清洁能源的自旋催化仍需独立电化学、结构及自旋测量，并结合催化机制建模。

## 验证与后续工作

运行时检查厄米性、时间反演、镜像关系、无 SOC 及单轨道对照、透射非负与通道上限、
稀疏求解和稠密参考一致性、概率守恒、eta 减小十倍的敏感性，并计算能量网格加密后的
热极化变化。失败时停止输出，实际残差与容差写入 JSON。验证证明的是模型内部一致性，
不等于真实材料精度、通用 CISS 理论或工业成熟度。
后续可研究界面参数标定、非弹性输运与 OER 微观动力学；本阶段不宣称重新验证了其他阶段。
"""
    for language,body in (('EN',en),('ZH',zh)):
        (root/f'CISS_QUANTUM_SPINTRONICS_REPORT_{language}.md').write_text(body+common,encoding='utf-8',newline='\n')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--sites',type=int,default=24)
    parser.add_argument('--energy-points',type=int,default=401)
    parser.add_argument('--soc-mev',type=float,default=30)
    parser.add_argument('--self-test',action='store_true')
    args=parser.parse_args(argv)
    if args.self_test:
        suite=unittest.defaultTestLoader.loadTestsFromTestCase(PhysicsTests)
        return 0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1
    if not 10 <= args.soc_mev <= 50:
        parser.error('--soc-mev must be in [10,50]; zero-SOC control is automatic')
    cfg=Config(sites=args.sites,energy_points=args.energy_points,soc_ev=args.soc_mev/1000)
    cfg.validate()
    print('20A/B: solving D/L, single-channel and zero-SOC controls...',flush=True)
    validation=audit(cfg)
    energy=np.linspace(-1,1,cfg.energy_points)
    d=spectrum(cfg,1,energy);l=spectrum(cfg,-1,energy)
    single=spectrum(replace(cfg,orbitals=1),1,energy)
    zero=spectrum(replace(cfg,soc_ev=0),1,energy)
    p=transport_average(cfg,energy,d)
    fine_energy=np.linspace(-1,1,2*cfg.energy_points-1)
    refined=spectrum(cfg,1,fine_energy)
    refined_p=transport_average(cfg,fine_energy,refined)
    grid_error=abs(p-refined_p)
    if grid_error>5e-4:
        raise RuntimeError(f'Thermal polarization grid error {grid_error:g}; increase --energy-points')
    # Whole-energy-window checks, in addition to audit's independent small reference.
    full_errors={
        'mirror_vector':float(np.nanmax(abs(l['polarization']-d['polarization']*[-1,1,-1]))),
        'single_null':float(np.nanmax(abs(single['polarization']))),
        'zero_soc_null':float(np.nanmax(abs(zero['polarization']))),
        'flux':float(np.max(abs(d['deficit']))),
        'residual':float(np.max(d['residual']))}
    if not all(np.isfinite(v) and v<2e-5 for v in full_errors.values()):
        raise RuntimeError(f'Full-window validation failed: {full_errors}')
    validation['full_window_errors']=full_errors
    kinetics=oer(cfg,p); optical=twins(cfg,energy,d)
    root=args.output_dir.resolve(); result=root/'results_phase20';figures=root/'figures_phase20'
    result.mkdir(parents=True,exist_ok=True);figures.mkdir(parents=True,exist_ok=True)
    summary=dict(model='effective two-path coherent Pauli/Rashba NEGF; assigned OER and optical sensitivities',
                 parameters=asdict(cfg),validation=validation,
                 peak_abs_pz_percent=float(100*np.nanmax(abs(d['polarization'][:,2]))),
                 thermal_pz=p,energy_grid_pz_difference=grid_error,
                 single_channel_max_p_percent=float(100*np.nanmax(abs(single['polarization']))),
                 oer_shift_mv=1000*kinetics['shift_v'],tafel_mv_dec=1000*kinetics['tafel_v_dec'],
                 environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
                 parameter_provenance='All geometry, hopping, SOC, contact, OER and optical parameters are assigned; no fitted experimental dataset.',
                 references=[dict(title=t,url=u) for t,u in REFERENCES])
    (result/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    def csv(name,columns,header):
        with (result/name).open('w',encoding='utf-8',newline='\n') as handle:
            np.savetxt(handle,np.column_stack(columns),delimiter=',',header=header,comments='',fmt='%.12g')
    csv('transport.csv',[energy,d['up'],d['down'],l['up'],l['down'],*d['polarization'].T,*l['polarization'].T,single['polarization'][:,2],zero['polarization'][:,2],d['deficit']],
        'energy_ev,D_Tup,D_Tdown,L_Tup,L_Tdown,D_Px,D_Py,D_Pz,L_Px,L_Py,L_Pz,single_Pz,zero_SOC_Pz,flux_deficit')
    csv('spin_resolved.csv',[energy,*d['spin_matrix'].reshape(-1,4).T,*l['spin_matrix'].reshape(-1,4).T],
        'energy_ev,D_up_up,D_up_down,D_down_up,D_down_down,L_up_up,L_up_down,L_down_up,L_down_down')
    xyz,radial,tangent,socaxis=helix(cfg,1)
    # Last node has no outgoing bond: direction entries are explicitly undefined.
    padded=[np.concatenate((v,np.full((1,cfg.orbitals,3),np.nan)),axis=0).reshape(-1,3) for v in (radial,tangent,socaxis)]
    csv('helix.csv',[np.repeat(np.arange(cfg.sites),cfg.orbitals),np.tile(np.arange(cfg.orbitals),cfg.sites),*xyz.reshape(-1,3).T,*padded[0].T,*padded[1].T,*padded[2].T],
        'site,orbital,x_angstrom,y_angstrom,z_angstrom,grad_x,grad_y,grad_z,tangent_x,tangent_y,tangent_z,soc_axis_x,soc_axis_y,soc_axis_z')
    csv('oer.csv',[kinetics['logj'],kinetics['achiral_eta'],kinetics['chiral_eta']], 'log10_j_ma_cm2,achiral_eta_v,conditional_chiral_eta_v')
    sensitivity_p=np.linspace(0,1,201)
    csv('oer_sensitivity.csv',[sensitivity_p,1000*cfg.oer_barrier_gain_ev*sensitivity_p**2/cfg.oer_alpha],
        'assigned_abs_pz,conditional_shift_mv')
    csv('spstm.csv',[optical['voltage'],*optical['currents'].T,optical['asymmetry']], 'bias_v,tip_plus_current_a,tip_minus_current_a,asymmetry')
    csv('mchd.csv',[optical['wavelength'],optical['absorption'],optical['delta'],-optical['delta'],optical['noisy_delta']], 'wavelength_nm,baseline_absorbance,D_delta_absorbance,L_delta_absorbance,synthetic_D_delta')
    print('20C/D: generating conditional kinetics, analytical twins and 300 DPI figures...',flush=True)
    figure_outputs(cfg,energy,d,l,single,kinetics,optical,figures)
    reports(cfg,summary,root)
    paths=sorted([*result.glob('*.csv'),result/'summary.json',*figures.glob('fig*.png'),*root.glob('CISS_QUANTUM_SPINTRONICS_REPORT_*.md')])
    manifest={str(p.relative_to(root)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    manifest['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (result/'sha256.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:summary[k] for k in ('peak_abs_pz_percent','thermal_pz','oer_shift_mv','energy_grid_pz_difference')},indent=2))
    print(f'Validated outputs: {root}',flush=True)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
