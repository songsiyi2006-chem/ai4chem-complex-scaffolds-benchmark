#!/usr/bin/env python3
"""Phase 22: molecular clock transitions, CCE-2 and ideal dynamical decoupling.

Standalone Python >=3.10; numpy, scipy and matplotlib. No network access at run time.
Run --self-test then run normally. Outputs default to the script's directory.
Assigned S=1/2, I=7/2 effective complex; not a fitted vanadium compound.
Thirty bath spins are treated with ensemble CCE-2, NOT a 2**30 density matrix.
Finite-pulse errors are studied for a projected gate; DD pulses are instantaneous.
No target fidelity, coherence gain or room-temperature claim is imposed.
"""
from __future__ import annotations
import argparse
from dataclasses import dataclass, asdict
import hashlib
import itertools
import json
from pathlib import Path
import platform
import unittest

import numpy as np
import scipy
from scipy import sparse
from scipy.linalg import eigh, expm
from scipy.optimize import brentq
from numpy.polynomial.hermite import hermgauss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

MU_B_H=13.99624555e9  # Hz/T
MU_N_H=7.622593229e6  # Hz/T
GAMMA_H=267.52218744e6  # rad/s/T, 1H
HBAR=1.054571817e-34
PLANCK=2*np.pi*HBAR
MU0_4PI=1e-7
PI2=2*np.pi
REFERENCES=[
 ('Shiddiq et al. (2016), Enhancing coherence in molecular spin qubits via atomic clock transitions', 'https://doi.org/10.1038/nature16984'),
 ('Yang and Liu, finite-spin-bath cluster-correlation expansion', 'https://arxiv.org/abs/0806.0098'),
 ('Yang and Liu, ensemble CCE dynamics', 'https://arxiv.org/abs/0902.3055'),
 ('Cluster-correlation expansion for studying decoherence of clock transitions in spin baths', 'https://arxiv.org/abs/2007.00412'),
]

@dataclass(frozen=True)
class Config:
    S: float=.5
    I: float=3.5
    g: tuple=(1.98,1.99,2.00)
    A_hz: tuple=(280e6,300e6,320e6)
    D_hz: tuple=(0.,0.,0.)  # S=1/2 has no genuine traceless ZFS
    P_hz: tuple=(-.1e6,-.1e6,.2e6)
    nuclear_g: float=1.471
    bath_spins: int=30
    rmin_angstrom: float=3.5
    rmax_angstrom: float=12.
    min_pair_angstrom: float=2.
    static_sigma_t: float=20e-6  # independent quasistatic technical field noise
    T1_s: float=.01  # assigned relaxation floor; not a spin-phonon calculation
    drive_hz: float=100e3
    amplitude_sigma: float=.005
    seed: int=2026
    time_points: int=121
    max_time_s: float=.04


def spin(j):
    m=np.arange(j,-j-1,-1);n=len(m)
    plus=np.zeros((n,n),complex)
    for col in range(1,n): plus[col-1,col]=np.sqrt(j*(j+1)-m[col]*(m[col]+1))
    return np.array([(plus+plus.T.conj())/2,(plus-plus.T.conj())/(2j),np.diag(m)])


class Molecule:
    def __init__(self,cfg):
        self.cfg=cfg;s=spin(cfg.S);i=spin(cfg.I)
        self.s=np.array([sparse.kron(x,sparse.eye(len(i[0]))).toarray() for x in s])
        self.i=np.array([sparse.kron(sparse.eye(len(s[0])),x).toarray() for x in i])
        self.hstatic=sum(cfg.D_hz[k]*self.s[k]@self.s[k]+cfg.A_hz[k]*self.s[k]@self.i[k]+cfg.P_hz[k]*self.i[k]@self.i[k] for k in range(3))
        self.db=MU_B_H*cfg.g[2]*self.s[2]-cfg.nuclear_g*MU_N_H*self.i[2]

    def hamiltonian(self,field):
        b=np.asarray(field)
        if b.ndim==0: return self.hstatic+float(b)*self.db
        if b.shape!=(3,): raise ValueError('Field must be scalar Bz or three components')
        return self.hstatic+sum(b[k]*(MU_B_H*self.cfg.g[k]*self.s[k]-self.cfg.nuclear_g*MU_N_H*self.i[k]) for k in range(3))

    def levels(self,b,curvature=False):
        e,v=eigh(self.hamiltonian(b));d=v.conj().T@self.db@v
        slopes=np.real(np.diag(d))
        if not curvature:return e,v,slopes
        delta=e[:,None]-e[None,:];np.fill_diagonal(delta,np.inf)
        if np.min(abs(delta))<1.: raise ValueError('Near-degenerate levels invalidate nondegenerate curvature')
        curves=2*np.sum(abs(d)**2/delta,axis=1).real
        return e,v,slopes,curves

    def transition(self,b,pair):
        e,v,d,c=self.levels(b,True);lo,hi=pair
        return dict(frequency_hz=float(e[hi]-e[lo]),slope_hz_t=float(d[hi]-d[lo]),curvature_hz_t2=float(c[hi]-c[lo]),
                    slopes=np.array([d[lo],d[hi]]),curves=np.array([c[lo],c[hi]]),
                    drive_element=float(abs(v[:,lo].conj()@self.s[0]@v[:,hi])),
                    readout_element=float(abs(v[:,lo].conj()@self.s[1]@v[:,hi])))

    def clock_transitions(self):
        # Dense search below 0.15 T and full 0..1 T plotted; high-field search retained.
        fields=np.unique(np.r_[np.linspace(.0001,.15,401),np.linspace(.15,1,171)])
        levels=[self.levels(b) for b in fields]
        energy=np.array([r[0] for r in levels]);d=np.array([r[2] for r in levels]);found=[]
        for lo,hi in itertools.combinations(range(energy.shape[1]),2):
            f=energy[:,hi]-energy[:,lo];grad=d[:,hi]-d[:,lo]
            crosses=np.where((grad[:-1]*grad[1:]<0)&(f[:-1]>50e6)&(f[:-1]<3e9))[0]
            for k in crosses:
                root=brentq(lambda b:self.levels(b)[2][hi]-self.levels(b)[2][lo],fields[k],fields[k+1],xtol=1e-14)
                try:
                    tr=self.transition(root,(lo,hi))
                except ValueError:
                    continue  # True crossings can jump sorted derivatives; not clock roots.
                if abs(tr['slope_hz_t'])<1. and tr['drive_element']>.03 and 100e6<tr['frequency_hz']<3e9 and root>.002:
                    found.append(dict(field_t=float(root),lo=lo,hi=hi,**{k:v for k,v in tr.items() if k not in ('slopes','curves')}))
        if not found:raise RuntimeError('No drive-accessible interior clock transition found')
        found.sort(key=lambda r:(-r['drive_element'],r['field_t']))
        return found


def bath_positions(cfg):
    rng=np.random.default_rng(cfg.seed);points=[]
    for _ in range(100000):
        u=rng.normal(size=3);u/=np.linalg.norm(u)
        r=(rng.uniform(cfg.rmin_angstrom**3,cfg.rmax_angstrom**3))**(1/3)
        p=u*r
        if all(np.linalg.norm(p-q)>=cfg.min_pair_angstrom for q in points):points.append(p)
        if len(points)==cfg.bath_spins:return np.array(points)
    raise RuntimeError('Bath placement failed')


def embed(op,k,n):
    a=np.ones((1,1),complex)
    for j in range(n):a=np.kron(a,op if j==k else np.eye(2))
    return a


class Bath:
    def __init__(self,cfg,points):
        self.cfg=cfg;self.points=points;self.n=len(points);self.proton=spin(.5)
        meters=points*1e-10;r=np.linalg.norm(meters,axis=1);u=meters/r[:,None]
        self.field=MU0_4PI*GAMMA_H*HBAR/r[:,None]**3*(3*u[:,2,None]*u-np.array([0,0,1]))
        self.pairs=list(itertools.combinations(range(self.n),2))
        self.templates={}

    def cluster(self,indices,b,tr,interactions=True):
        n=len(indices);ops=[[embed(s,k,n) for s in self.proton] for k in range(n)]
        hb=np.zeros((2**n,2**n),complex);beta=hb.copy()
        for k,i in enumerate(indices):
            hb-=GAMMA_H/PI2*b*ops[k][2]
            beta+=sum(self.field[i,a]*ops[k][a] for a in range(3))
        if interactions:
            for k,l in itertools.combinations(range(n),2):
                r=(self.points[indices[k]]-self.points[indices[l]])*1e-10;dist=np.linalg.norm(r);u=r/dist
                pref=MU0_4PI*(GAMMA_H*HBAR)**2/(PLANCK*dist**3) # Hz, not rad/s
                hb+=pref*(sum(ops[k][a]@ops[l][a] for a in range(3))-3*sum(u[a]*ops[k][a] for a in range(3))@sum(u[a]*ops[l][a] for a in range(3)))
        return np.array([hb+tr['slopes'][a]*beta+.5*tr['curves'][a]*(beta@beta) for a in (0,1)])

    def prepared(self,b,tr,order=2):
        groups=[[(i,) for i in range(self.n)]]
        if order>=2:groups.append(self.pairs)
        return [np.linalg.eigh(np.array([self.cluster(c,b,tr) for c in group])) for group in groups]


def pulse_fractions(sequence,n):
    if sequence=='FID':return np.array([])
    if n<1:raise ValueError('Pulse count must be positive')
    if sequence in ('CPMG','Hahn'):return (np.arange(n)+.5)/n
    if sequence=='UDD':return np.sin(np.pi*np.arange(1,n+1)/(2*n+2))**2
    raise ValueError('Unknown pulse sequence')


def cluster_coherence(prepared,times,sequence='FID',n=0):
    """Exact conditional cluster propagators, maximally mixed cluster bath.
    Arrays include all clusters of a fixed size; temporal propagation is batched.
    L=Tr(U1†U0)/dim. Frequencies in Hz require -i 2 pi H t.
    """
    values,vectors=prepared;dim=values.shape[-1];times=np.atleast_1d(times)
    cache={}
    def propagator(state,fraction):
        key=(state,round(float(fraction),15))
        if key not in cache:
            e=values[:,state];v=vectors[:,state]
            phase=np.exp(-1j*PI2*e[:,None,:]*times[None,:,None]*fraction)
            cache[key]=(v[:,None,:,:]*phase[:,:,None,:])@v.conj().transpose(0,2,1)[:,None,:,:]
        return cache[key]
    if sequence=='FID':u0=propagator(0,1);u1=propagator(1,1)
    elif sequence in ('CPMG','Hahn') and n==1:
        u0=propagator(1,.5)@propagator(0,.5);u1=propagator(0,.5)@propagator(1,.5)
    elif sequence=='CPMG' and n%2==0:
        def branch(s):
            half=propagator(s,.5/n);other=propagator(1-s,1/n);full=propagator(s,1/n)
            return half@np.linalg.matrix_power(other@full,n//2-1)@other@half
        u0=branch(0);u1=branch(1)
    else:
        durations=np.diff(np.r_[0,pulse_fractions(sequence,n),1]);u0=np.eye(dim);u1=np.eye(dim)
        for k,f in enumerate(durations):
            u0=propagator(k%2,f)@u0;u1=propagator(1-k%2,f)@u1
    return np.sum(u1.conj()*u0,axis=(-1,-2))/dim


def cce(prepared,pairs,times,sequence='FID',n=0):
    single=cluster_coherence(prepared[0],times,sequence,n)
    logs=np.sum(np.log(single.astype(complex)),axis=0)
    denominator_ok=np.min(abs(single),axis=0)>1e-7
    if len(prepared)>1:
        pair=cluster_coherence(prepared[1],times,sequence,n)
        for k,(i,j) in enumerate(pairs):logs+=np.log(pair[k].astype(complex))-np.log(single[i].astype(complex))-np.log(single[j].astype(complex))
    raw=np.exp(np.minimum(logs.real,700)+1j*logs.imag)
    valid=denominator_ok&np.isfinite(raw)&(abs(raw)<=1+1e-5)
    # CCE truncation is not guaranteed completely positive: never silently clip.
    if np.any(~valid):valid[np.flatnonzero(~valid)[0]:]=False
    result=np.where(valid,raw,np.nan+1j*np.nan)
    return result,valid


def static_noise(times,tr,sigma):
    """Exact Gaussian average of dnu=d b + (c/2)b², b~N(0,sigma²)."""
    z=1+1j*PI2*tr['curvature_hz_t2']*sigma*sigma*times
    return z**(-.5)*np.exp(-.5*(PI2*tr['slope_hz_t']*sigma*times)**2/z)


def envelope(cfg,prepared,pairs,times,tr,sequence,n):
    value,valid=cce(prepared,pairs,times,sequence,n)
    if sequence=='FID':value*=static_noise(times,tr,cfg.static_sigma_t)
    value*=np.exp(-times/(2*cfg.T1_s))
    return value,valid


def lifetime(times,values):
    a=abs(values);valid=np.isfinite(a);hits=np.flatnonzero(valid&(a<=np.exp(-1)))
    if len(hits):
        k=hits[0]
        if k==0:return dict(seconds=float(times[0]),status='left_censored')
        t=float(np.interp(np.exp(-1),[a[k],a[k-1]],[times[k],times[k-1]]))
        return dict(seconds=t,status='crossing')
    if not valid.all():return dict(seconds=None,status='CCE_unresolved',last_valid_s=float(times[np.flatnonzero(valid)[-1]]) if valid.any() else 0.)
    return dict(seconds=None,status='right_censored',lower_bound_s=float(times[-1]))


def gate_metrics(cfg,tr):
    rabihz=2*cfg.drive_hz*tr['drive_element'];tpi=1/(2*rabihz)
    nodes,weights=hermgauss(15);fidelity=0.
    sx=np.array([[0,1],[1,0]],complex);sz=np.diag([1.,-1.]);target=-1j*sx
    for i,x in enumerate(nodes):
        b=np.sqrt(2)*cfg.static_sigma_t*x
        det=tr['slope_hz_t']*b+.5*tr['curvature_hz_t2']*b*b
        for j,y in enumerate(nodes):
            amp=1+np.sqrt(2)*cfg.amplitude_sigma*y
            u=expm(-1j*np.pi*tpi*(rabihz*amp*sx+det*sz))
            fidelity+=weights[i]*weights[j]/np.pi*(abs(np.trace(target.conj().T@u))**2+2)/6
    return dict(rabi_hz=float(rabihz),pi_pulse_s=float(tpi),projected_average_gate_fidelity=float(fidelity),
                scope='Two-level rotating-wave gate with assigned static detuning/amplitude noise; no leakage or full pulse/bath dynamics')


def filter_amplitude(z,sequence,n):
    pulses=pulse_fractions(sequence,n)
    # Dimensionless numerator of i omega integral y(t) exp(i omega t) dt.
    ans=1+(-1)**(n+1)*np.exp(1j*z) if sequence!='FID' else 1-np.exp(1j*z)
    for j,f in enumerate(pulses,1):ans+=2*(-1)**j*np.exp(1j*z*f)
    return abs(ans)**2


def validation(mol,ct,bath):
    b=ct['field_t'];pair=(ct['lo'],ct['hi']);tr=mol.transition(b,pair);eps=1e-6
    finite=(mol.transition(b+eps,pair)['frequency_hz']-mol.transition(b-eps,pair)['frequency_hz'])/(2*eps)
    second=(mol.transition(b+eps,pair)['frequency_hz']-2*tr['frequency_hz']+mol.transition(b-eps,pair)['frequency_hz'])/eps**2
    checks=dict(hermiticity_hz=float(np.max(abs(mol.hamiltonian(b)-mol.hamiltonian(b).conj().T))),
                clock_gradient_hz_t=abs(tr['slope_hz_t']),gradient_FD_error_hz_t=float(abs(finite-tr['slope_hz_t'])),
                curvature_relative_error=float(abs(second/tr['curvature_hz_t2']-1)))
    limits=dict(hermiticity_hz=1e-6,clock_gradient_hz_t=1.,gradient_FD_error_hz_t=2e4,curvature_relative_error=1e-3)
    if any(not np.isfinite(v) or v>limits[k] for k,v in checks.items()):raise RuntimeError(f'CT validation failed: {checks}')
    # Small bath exact density-matrix reference versus CCE-2; not a global CCE convergence proof.
    small=Bath(bath.cfg,bath.points[:3]);prep=small.prepared(b,tr)
    times=np.array([0.,1e-7,1e-6])
    approx,_=cce(prep,small.pairs,times,'Hahn',1)
    exact=cluster_coherence(np.linalg.eigh(small.cluster((0,1,2),b,tr)[None]),times,'Hahn',1)[0]
    error=float(np.max(abs(approx-exact)))
    return dict(status='passed',errors=checks,tolerances=limits,three_spin_CCE2_exact_difference=error,
                three_spin_comparison_times_s=times.tolist(),CCE_scope='Order 2 only; 30-spin convergence not established')


class Tests(unittest.TestCase):
    def test_spin_algebra(self):
        for j in (.5,1,3.5):
            s=spin(j);self.assertTrue(np.allclose(s[0]@s[1]-s[1]@s[0],1j*s[2]));self.assertTrue(np.allclose(sum(x@x for x in s),j*(j+1)*np.eye(len(s[0]))))
    def test_clock(self):
        m=Molecule(Config());cts=m.clock_transitions();self.assertGreater(len(cts),0)
        t=m.transition(cts[0]['field_t'],(cts[0]['lo'],cts[0]['hi']));self.assertLess(abs(t['slope_hz_t']),1.)
    def test_bath(self):
        cfg=Config();p=bath_positions(cfg);self.assertEqual(len(p),30);self.assertTrue(np.array_equal(p,bath_positions(cfg)))
        self.assertTrue(np.all(np.linalg.norm(p,axis=1)>=3.5));self.assertTrue(np.all(np.linalg.norm(p,axis=1)<=12))
    def test_propagation(self):
        rng=np.random.default_rng(8);a=rng.normal(size=(2,2,2));h=a+a.transpose(0,2,1)
        prep=np.linalg.eigh(h[None]);t=.017
        for seq,n in [('FID',0),('Hahn',1),('CPMG',4),('UDD',4),('CPMG',64)]:
            durations=np.diff(np.r_[0,pulse_fractions(seq,n),1]);u0=np.eye(2,dtype=complex);u1=u0.copy()
            for k,f in enumerate(durations):u0=expm(-1j*PI2*h[k%2]*t*f)@u0;u1=expm(-1j*PI2*h[1-k%2]*t*f)@u1
            expected=np.trace(u1.conj().T@u0)/2
            self.assertAlmostEqual(abs(cluster_coherence(prep,[t],seq,n)[0,0]-expected),0,places=10)
    def test_identical_branches(self):
        cfg=Config(bath_spins=3);bath=Bath(cfg,bath_positions(cfg));tr=dict(slopes=np.array([1e9,1e9]),curves=np.array([1e10,1e10]))
        prep=bath.prepared(.02,tr);v,ok=cce(prep,bath.pairs,np.linspace(0,1e-4,10),'CPMG',16)
        self.assertTrue(ok.all());self.assertTrue(np.allclose(v,1,atol=1e-10))
    def test_lifetimes(self):
        t=np.linspace(0,2,101);self.assertAlmostEqual(lifetime(t,np.exp(-t))['seconds'],1)
        self.assertEqual(lifetime(t,np.ones(101))['status'],'right_censored')
        self.assertEqual(lifetime(t,np.r_[1,np.full(100,np.nan)])['status'],'CCE_unresolved')
    def test_filter_dc(self):
        for seq in ('CPMG','UDD'):
            for n in (1,4,16,64):self.assertLess(filter_amplitude(0,seq,n),1e-20)
    def test_static_and_gate(self):
        cfg=Config(static_sigma_t=0,amplitude_sigma=0);tr=dict(slope_hz_t=0.,curvature_hz_t2=1e12,drive_element=.5)
        self.assertTrue(np.allclose(static_noise(np.linspace(0,.01,10),tr,0),1))
        self.assertAlmostEqual(gate_metrics(cfg,tr)['projected_average_gate_fidelity'],1,places=12)


def make_figures(cfg,data,out):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    def save(fig,name):fig.savefig(out/name,dpi=300,facecolor='white');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    axs[0].plot(data['fields'],data['energies']/1e9,lw=.65);axs[0].axvline(data['ct']['field_t'],ls='--',color='black')
    axs[0].set(xlabel='Bz (T)',ylabel='Energy / h (GHz)',title='22A | S=1/2, I=7/2 assigned spin Hamiltonian')
    axs[1].semilogy(data['local_fields'],np.maximum(abs(data['local_slopes'])/1e9,1e-10))
    axs[1].axvline(data['ct']['field_t'],ls='--',color='black',label=f"CT = {data['ct']['field_t']:.7f} T")
    axs[1].set(xlabel='Bz (T)',ylabel='|d nu / d Bz| (GHz/T)',title='First-order axial protection; curvature is finite');axs[1].legend(fontsize=8)
    save(fig,'fig1_breit_rabi_clock_transitions.png')
    fig=plt.figure(figsize=(11,5),layout='constrained');ax=fig.add_subplot(121);bloch=fig.add_subplot(122,projection='3d')
    for power in (.5,1,2):ax.plot(data['rabi_times']*1e6,np.sin(np.pi*data['gate']['rabi_hz']*power*data['rabi_times'])**2,label=f'Drive amplitude x{power}')
    ax.set(xlabel='Pulse time (microsecond)',ylabel='Excited population',title='22C | Projected resonant Rabi nutation');ax.legend()
    u,v=np.meshgrid(np.linspace(0,2*np.pi,25),np.linspace(0,np.pi,15));bloch.plot_wireframe(np.cos(u)*np.sin(v),np.sin(u)*np.sin(v),np.cos(v),color='grey',alpha=.2,lw=.4)
    th=np.linspace(0,np.pi,101);bloch.plot(np.zeros_like(th),-np.sin(th),np.cos(th),color='tab:blue',lw=2)
    bloch.scatter([0,0],[ -1,0],[0,-1],color='red');bloch.set(xlabel='x',ylabel='y',zlabel='z',title='Ideal pi/2 and pi rotations');bloch.set_box_aspect((1,1,1))
    save(fig,'fig2_rabi_nutation_and_bloch_sphere.png')
    fig,axs=plt.subplots(1,2,figsize=(11,5),layout='constrained')
    for key in ('FID','Hahn','CPMG4','CPMG16','CPMG64'):axs[0].semilogx(data['times']*1e6,abs(data['curves'][key]),label=key)
    axs[0].axhline(np.exp(-1),color='grey',ls=':');axs[0].set(xlabel='Total evolution time (microsecond)',ylabel='|L(t)|',ylim=(0,1.05),title='22B/C | CT, CCE-2 + assigned T1 floor');axs[0].legend(fontsize=8)
    for n in (1,4,16,64):axs[1].loglog(data['filter_z'],np.maximum(filter_amplitude(data['filter_z'],'CPMG',n),1e-12),label=f'CPMG-{n}')
    axs[1].set(xlabel='Angular frequency x total time',ylabel='Dimensionless filter numerator',title='Ideal instantaneous-pulse filter functions');axs[1].legend(fontsize=8)
    save(fig,'fig3_dynamical_decoupling_coherence_gain.png')
    fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for k in (0,len(data['epr_fields'])//2,-1):axs[0,0].plot(data['epr_times']*1e6,data['epr_echo'][k],label=f"B={data['epr_fields'][k]:.4f} T")
    axs[0,0].set(xlabel='Total echo time 2 tau (microsecond)',ylabel='Normalized echo quadrature',title='22D | CCE-1 ESEEM preview');axs[0,0].legend(fontsize=8)
    mesh=axs[0,1].pcolormesh(data['fft_freq']/1e6,data['epr_fields'],data['fft'],shading='auto',cmap='magma')
    axs[0,1].set(xlabel='Modulation frequency (MHz)',ylabel='Bz (T)',title='FT of mean-subtracted Hahn modulation');fig.colorbar(mesh,ax=axs[0,1],label='Amplitude (a.u.)')
    for key,label in [('field_fid','FID T2*'),('field_hahn','Hahn T2')]:
        vals=[r['seconds'] if r['status']=='crossing' else np.nan for r in data[key]]
        axs[1,0].semilogy(data['lifetime_fields'],np.array(vals)*1e6,'o-',label=label)
    axs[1,0].axvline(data['ct']['field_t'],ls=':',color='black');axs[1,0].set(xlabel='Bz (T)',ylabel='First 1/e crossing (microsecond)',title='CCE-2 field map; missing = unresolved/censored');axs[1,0].legend(fontsize=8)
    axs[1,1].plot(data['ramsey_times']*1e6,data['ramsey_signal']);axs[1,1].set(xlabel='Ramsey delay (microsecond)',ylabel='Readout probability',title='Assigned 100 kHz Ramsey detuning')
    fig.suptitle('Synthetic signals; no Bruker calibration or structural inversion')
    save(fig,'fig4_pulsed_epr_eseem_twin.png')


def reports(cfg,summary,root):
    ct=summary['clock_transition'];gate=summary['gate']
    rows='\n'.join(f"| {key} | {v['seconds'] if v['seconds'] is not None else 'not resolved'} | {v['status']} |" for key,v in summary['lifetimes'].items())
    shared=f"""
## Computed values / 计算结果

- Selected CT / 选定时钟跃迁: **{ct['field_t']:.10f} T**, states {ct['lo']} and {ct['hi']} (ascending energies).
- Transition / 跃迁频率: **{ct['frequency_hz']/1e9:.8f} GHz**.
- Axial slope / 轴向一阶导数: **{ct['slope_hz_t']:.6g} Hz/T**.
- Curvature / 二阶导数: **{ct['curvature_hz_t2']:.6g} Hz/T²**.
- pi pulse / π 脉冲: **{gate['pi_pulse_s']*1e6:.6f} microseconds**.
- Projected noisy gate fidelity / 投影含噪声门保真度: **{100*gate['projected_average_gate_fidelity']:.6f}%**.
- CPMG-64/FID crossing ratio / 阈值交点增益: **{summary['cpmg64_gain']}** (null = not resolvable).

| Sequence / 序列 | Lifetime (seconds) / 寿命 | Status / 状态 |
|---|---:|---|
{rows}

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

"""+'\n'.join(f'- [{title}]({url})' for title,url in REFERENCES)+'\n'
    en="""# Phase 22: molecular spin qubits, clock transitions and dynamical decoupling

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
"""
    zh="""# 第二十二阶段：分子自旋量子比特、时钟跃迁与动力学解耦

## 模型范围

采用指定参数的 S=1/2、I=7/2 配位化合物有效模型，不是已拟合的钒或镧系分子。
包含电子/核塞曼项、各向异性超精细与核四极矩张量，统一用 Hz。
程序支持 ZFS 项，但 S=1/2 没有真正的无迹零场劈裂，所以 D 取零。
稀疏 Kronecker 构造中心自旋算符，小维度矩阵采用厄米对角化。
扫描 0 至 1 T，通过 Hellmann-Feynman 导数找根，再以微波矩阵元筛选可驱动跃迁。
二阶导数包含其他中心能级的虚跃迁，并与有限差分核验。

## 时钟保护的含义

CT 的 dnu/dBz=0 只消除纵向磁场的一阶灵敏度，二阶曲率通常非零。
不能据此声称对所有磁噪声、应变、超精细噪声、电子翻转或声子免疫。
所选两能级是较不敏感的工作点，不是严格无退相干子空间，也没有实现量子纠错编码。
场扫描采用按能量排序的绝热本征态，不是跨场非绝热脉冲追踪。

## 核浴与 CCE

固定随机种子，在 3.5 至 12 埃球壳内按体积均匀放置 30 个质子，最近间距 2 埃。
坐标为合成环境，不是真实晶格。保留完整核偶极张量、非长期近似项和质子塞曼作用，
显式区分角频率与 Hz。用中心能级的 E_a(B+beta)=E_a+E'_a beta+E''_a beta²/2
形成两组条件核浴哈密顿量。beta 为质子的纵向偶极场算符；二阶交叉项保留于成对簇中，
避免在 CT 人为把所有退相干关闭。模型仅处理有效纵向场和纯退相干，不含完整中心自旋翻转。

高核自旋温度下初态为最大混合态，每个单核和双核簇进行精确条件传播，
L_C=Tr(U1†U0)/dim；CCE-2 用单簇乘积和 L_ij/(L_i L_j) 的不可约成对关联组合结果。
它保留有限记忆演化，但不是精确的 30 核全密度矩阵：后者有 2^60 个复数元素，
complex128 约需 18.4 EB 内存。提供 CCE-1 对比和三核精确参考，不能据此宣称
30 核在 CT 附近已经达到高阶收敛；未包含簇外均场和更高阶簇。
分母过小或截断导致 |L| 超出物理范围时标记后续无效，不通过裁剪掩盖失败。

另加入独立的指定准静态技术磁噪声（20 微特斯拉 rms），FID 使用包含曲率的解析高斯平均，
理想平衡回波抵消该静态项。所有曲线另乘指定 T1=10 ms 对应的 exp(-t/2T1)。
这不是计算或测量所得的自旋声子弛豫时间，也不是从核浴中拟合出来的参数。

## 驱动、解耦和寿命

由实际所选 Sx 矩阵元计算投影旋波模型的 Rabi 频率与 π 脉冲时间。
用高斯求积平均失谐及 0.5% 幅度噪声，通过 (|Tr(Utarget†U)|²+2)/6 计算平均门保真度。
未计入向其他分子能级泄漏和脉冲期间核浴演化，因此不是可直接用于实验的门精度。
Bloch 球为理想二能级旋转，不预设保真度必须超过 99.9%。

CPMG 和 UDD 均实现 1、4、16、64 个理想瞬时 π 脉冲，通过条件哈密顿量交替传播，
CPMG 重复块用矩阵整数幂加速，UDD 非等间隔逐段相乘。导出滤波函数，但不以
滤波示意曲线替代核浴动力学；其作用不是一个绝对 Nyquist 截止频率。
有限脉宽、累计控制误差和硬件带宽未纳入 DD，单门脉宽单独报告。

寿命定义为有限时间网格的首次 1/e 交点，允许非指数衰减与复苏。
超出观测窗记为右删失，CCE 失效记为未确定，不虚构 T2。
只有在同一场和同一噪声配置下 FID 与 CPMG-64 均有可确定交点，才报告增益；
不强制得到超过两数量级延长。

## EPR 信号与化学展望

相干寿命随磁场图使用 CCE-2；较密集 ESEEM 光谱图明确采用 CCE-1，包含条件核进动，
不含该图中的成对核关联。Hahn 信号去均值、乘 Hann 窗后进行实 FFT；频率共轭变量
是总回波时间 2tau，不是未经验证的同位素或配位结构鉴定。
归一化回波正交分量为 Re L，换算成电子 Sy 单位需乘投影读出矩阵元。
Ramsey 的 100 kHz 失谐为指定设置。未模拟 Bruker 仪器响应或使用真实 EPR 测量。

配位化学可调控分子组成和配体，但不保证晶体/溶剂环境完全一致，不能断言这些能力
是半导体缺陷体系无法实现的。氘代减小旋磁比，却引入 I=1 四极矩核；12C 富集可减少
13C 核自旋，其他核、声子和无序仍存在。这些是实验设计方向，不是本脚本已完成的
材料优化，也不是室温分子量子处理器的证据。参考实验仅用于研究背景。
"""
    for lang,body in [('EN',en),('ZH',zh)]:
        (root/f'MOLECULAR_SPIN_QUBIT_REPORT_{lang}.md').write_text(body+shared,encoding='utf-8',newline='\n')


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--self-test',action='store_true');ap.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
    args=ap.parse_args(argv)
    if args.self_test:return 0 if unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests)).wasSuccessful() else 1
    cfg=Config();mol=Molecule(cfg);cts=mol.clock_transitions();ct=cts[0];pair=(ct['lo'],ct['hi']);b=ct['field_t'];tr=mol.transition(b,pair)
    points=bath_positions(cfg);bath=Bath(cfg,points);checks=validation(mol,ct,bath);gate=gate_metrics(cfg,tr)
    print(f'Clock found: B={b:.9f} T, transition={tr["frequency_hz"]/1e9:.6f} GHz',flush=True)
    times=np.r_[0,np.geomspace(1e-8,cfg.max_time_s,cfg.time_points)]
    prep=bath.prepared(b,tr);curves={};lifetimes={};validity={}
    sequences=[('FID','FID',0),('Hahn','Hahn',1)]+[(f'{s}{n}',s,n) for s in ('CPMG','UDD') for n in (1,4,16,64)]
    for key,seq,n in sequences:
        curves[key],valid=envelope(cfg,prep,bath.pairs,times,tr,seq,n);lifetimes[key]=lifetime(times,curves[key]);validity[key]=int(np.sum(valid))
        print(f'{key}: {lifetimes[key]}',flush=True)
    cce1={key:envelope(cfg,prep[:1],bath.pairs,times,tr,seq,n)[0] for key,seq,n in sequences}
    gain=None
    if lifetimes['FID']['status']=='crossing' and lifetimes['CPMG64']['status']=='crossing':gain=lifetimes['CPMG64']['seconds']/lifetimes['FID']['seconds']
    fields=np.linspace(0,1,401);energies=np.array([mol.levels(f)[0] for f in fields])
    local_fields=np.sort(np.r_[np.linspace(max(.0001,b-.008),min(.999,b+.008),201),b]);local_slopes=np.array([mol.transition(f,pair)['slope_hz_t'] for f in local_fields])
    lifetime_fields=np.sort(np.r_[np.linspace(max(.0001,b-.006),b+.006,8),b]);field_fid=[];field_hahn=[]
    for f in lifetime_fields:
        t=mol.transition(f,pair);p=bath.prepared(f,t)
        for name,seq,n in [(field_fid,'FID',0),(field_hahn,'Hahn',1)]:name.append(lifetime(times,envelope(cfg,p,bath.pairs,times,t,seq,n)[0]))
    print('Field lifetime map complete; generating CCE-1 spectroscopic grid.',flush=True)
    epr_fields=np.linspace(max(.0001,b-.01),b+.01,17)
    # Sample at >8 times the largest proton Larmor frequency. FFT axis uses total echo time.
    dt=1/(8*GAMMA_H/PI2*max(epr_fields));epr_times=np.arange(512)*dt;epr_echo=[]
    for f in epr_fields:
        t=mol.transition(f,pair);p=bath.prepared(f,t,1);v,_=envelope(cfg,p,bath.pairs,epr_times,t,'Hahn',1);epr_echo.append(v.real)
    epr_echo=np.array(epr_echo)
    if not np.isfinite(epr_echo).all():raise RuntimeError('Spectroscopic grid left CCE validity domain')
    fft=abs(np.fft.rfft((epr_echo-epr_echo.mean(axis=1,keepdims=True))*np.hanning(len(epr_times)),axis=1))*2/np.sum(np.hanning(len(epr_times)))
    fft_freq=np.fft.rfftfreq(len(epr_times),dt)
    ramsey_times=np.linspace(0,min(cfg.max_time_s,80e-6),401);ramsey_coh,_=envelope(cfg,prep,bath.pairs,ramsey_times,tr,'FID',0)
    ramsey_signal=.5*(1+(ramsey_coh*np.exp(1j*PI2*100e3*ramsey_times)).real)
    root=args.output_dir.resolve();res=root/'results_phase22';fig=root/'figures_phase22';res.mkdir(parents=True,exist_ok=True);fig.mkdir(parents=True,exist_ok=True)
    data=dict(ct=ct,fields=fields,energies=energies,local_fields=local_fields,local_slopes=local_slopes,gate=gate,rabi_times=np.linspace(0,2*gate['pi_pulse_s'],301),times=times,curves=curves,filter_z=np.geomspace(.01,1000,1001),
              epr_fields=epr_fields,epr_times=epr_times,epr_echo=epr_echo,fft_freq=fft_freq,fft=fft,lifetime_fields=lifetime_fields,field_fid=field_fid,field_hahn=field_hahn,ramsey_times=ramsey_times,ramsey_signal=ramsey_signal)
    summary=dict(parameters=asdict(cfg),clock_transition=ct,gate=gate,validation=checks,lifetimes=lifetimes,cpmg64_gain=gain,valid_time_points=validity,
                 field_lifetimes=[dict(field_t=float(f),FID=a,Hahn=c) for f,a,c in zip(lifetime_fields,field_fid,field_hahn)],
                 environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
                 provenance='Assigned effective spin model and synthetic 30-proton geometry; CCE-2 truncated bath, ideal DD, CCE-1 spectroscopy; not experimental data')
    def csv(name,columns,header):
        with (res/name).open('w',encoding='utf-8',newline='\n') as handle:np.savetxt(handle,np.column_stack(columns),delimiter=',',header=header,comments='',fmt='%.12g')
    csv('bath_geometry.csv',[np.arange(30),*points.T,*bath.field.T],'nucleus,x_angstrom,y_angstrom,z_angstrom,beta_x_t,beta_y_t,beta_z_t')
    csv('breit_rabi.csv',[fields,*energies.T],'B_t,'+','.join(f'E{i}_hz' for i in range(energies.shape[1])))
    csv('clock_candidates.csv',[[r[k] for r in cts] for k in ('field_t','lo','hi','frequency_hz','slope_hz_t','curvature_hz_t2','drive_element')],'B_t,lo,hi,nu_hz,slope_hz_t,curvature_hz_t2,Sx_element')
    csv('clock_gradient.csv',[local_fields,local_slopes],'B_t,gradient_hz_t')
    csv('coherence.csv',[times,*[a for v in curves.values() for a in (v.real,v.imag,abs(v))]],'time_s,'+','.join(f'{k}_{v}' for k in curves for v in ('real','imag','abs')))
    csv('cce1_comparison.csv',[times,*[abs(v) for v in cce1.values()]],'time_s,'+','.join(cce1))
    csv('filter_functions.csv',[data['filter_z'],*[filter_amplitude(data['filter_z'],s,n) for s in ('CPMG','UDD') for n in (1,4,16,64)]],'omega_t,'+','.join(f'{s}{n}' for s in ('CPMG','UDD') for n in (1,4,16,64)))
    csv('rabi.csv',[data['rabi_times'],*[np.sin(np.pi*gate['rabi_hz']*p*data['rabi_times'])**2 for p in (.5,1,2)]],'time_s,amplitude_half,amplitude_one,amplitude_two')
    csv('ramsey.csv',[ramsey_times,ramsey_signal],'time_s,readout_probability')
    f,t=np.meshgrid(epr_fields,epr_times,indexing='ij');csv('hahn_eseem.csv',[f.ravel(),t.ravel(),epr_echo.ravel()],'B_t,total_echo_time_s,normalized_quadrature')
    f,w=np.meshgrid(epr_fields,fft_freq,indexing='ij');csv('eseem_fft.csv',[f.ravel(),w.ravel(),fft.ravel()],'B_t,modulation_frequency_hz,amplitude')
    (res/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    make_figures(cfg,data,fig);reports(cfg,summary,root)
    paths=sorted([*res.glob('*.csv'),res/'summary.json',*fig.glob('*.png'),*root.glob('MOLECULAR_SPIN_QUBIT_REPORT_*.md')])
    hashes={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths};hashes['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (res/'sha256.json').write_text(json.dumps(hashes,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(f'Completed: {root}; measured crossing gain={gain}',flush=True)
    return 0

if __name__=='__main__':raise SystemExit(main())
