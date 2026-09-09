#!/usr/bin/env python3
"""Phase 23: reproducible effective-model in-materio computing benchmark.

Run: python run_phase23_in_materio_neuromorphic_computing.py [--self-test]
Dependencies: numpy, scipy.linalg, matplotlib. No downloaded datasets or ML library.
All model parameters are assigned, not fitted to a particular molecular complex.
No experimental STDP, sub-fJ system energy or AGI performance is asserted.
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
from scipy.linalg import solve, svd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


@dataclass(frozen=True)
class Config:
    seed: int = 2026
    glow_s: float = 1e-6
    ghigh_s: float = 20e-6
    alpha_per_v: float = 3.0
    eta_per_s: float = 8.0
    v0_v: float = .18
    threshold_v: float = .03
    tau_s: float = .010
    equilibrium: float = .5
    window_p: int = 1
    nodes: int = 25
    source_g_s: float = 100e-6
    wire_g_s: float = 8e-6
    node_c_f: float = 2e-12
    sample_s: float = .0005
    substeps: int = 4
    lorenz_step: float = .02
    samples: int = 3600
    washout: int = 200
    train_end: int = 2200
    validation_end: int = 2800


def trapz(y, x, axis=-1):
    """Trapezoidal integral compatible with NumPy 1.23 and NumPy 2.x."""
    y = np.asarray(y)
    a = np.take(y, np.arange(y.shape[axis]-1), axis=axis)
    b = np.take(y, np.arange(1, y.shape[axis]), axis=axis)
    shape = [1]*y.ndim
    shape[axis] = -1
    return np.sum((a+b)*.5*np.diff(x).reshape(shape), axis=axis)


class Junction:
    def __init__(self, cfg, scale=1., tau_scale=1., rate_scale=1.):
        self.cfg = cfg
        self.scale = np.asarray(scale)
        self.tau = cfg.tau_s*np.asarray(tau_scale)
        self.eta = cfg.eta_per_s*np.asarray(rate_scale)

    def conductance(self, w):
        c = self.cfg
        return self.scale*(c.glow_s+(c.ghigh_s-c.glow_s)*w)

    def current(self, v, w):
        # Odd, passive corrected law: I V >= 0; small-V slope is G(w).
        a = self.cfg.alpha_per_v
        return self.conductance(w)*np.sinh(a*np.asarray(v))/a

    def drive(self, v):
        # Symmetric dead band, not a one-sided threshold that drives at zero bias.
        q = np.sign(v)*np.maximum(np.abs(v)-self.cfg.threshold_v, 0.)
        return self.eta*np.sinh(q/self.cfg.v0_v)

    def rhs(self, w, v):
        window = 1-(2*w-1)**(2*self.cfg.window_p)
        return self.drive(v)*window-(w-self.cfg.equilibrium)/self.tau

    def advance(self, w, v, dt):
        """RK4 for a held voltage; step subdivision bounds stiffness, no clipping."""
        stiffness = np.max(4*self.cfg.window_p*np.abs(self.drive(v))+1/self.tau)
        count = max(1, int(np.ceil(dt*stiffness/.15)))
        h = dt/count
        for _ in range(count):
            k1 = self.rhs(w, v)
            k2 = self.rhs(w+h*k1/2, v)
            k3 = self.rhs(w+h*k2/2, v)
            k4 = self.rhs(w+h*k3, v)
            w = w+h*(k1+2*k2+2*k3+k4)/6
        if np.any(~np.isfinite(w)) or np.any(w < -1e-12) or np.any(w > 1+1e-12):
            raise FloatingPointError('State left [0,1]; no state clipping is permitted')
        return w


def sinusoidal_cycle(j, frequency, initial, points=800):
    n = max(points, int(np.ceil(1/(frequency*j.cfg.tau_s/30))))
    n += (-n) % 4
    t = np.linspace(0, 1/frequency, n+1)
    v = .4*np.sin(2*np.pi*frequency*t)
    vm = .4*np.sin(2*np.pi*frequency*(t[:-1]+np.diff(t)/2))
    w = np.empty(n+1); w[0] = initial
    for k, h in enumerate(np.diff(t)):
        # Nonautonomous RK4 evaluates the sine at the stages, not a DC surrogate.
        q = w[k]
        k1 = j.rhs(q, v[k]); k2 = j.rhs(q+h*k1/2, vm[k])
        k3 = j.rhs(q+h*k2/2, vm[k]); k4 = j.rhs(q+h*k3, v[k+1])
        w[k+1] = q+h*(k1+2*k2+2*k3+k4)/6
    if np.min(w) < 0 or np.max(w) > 1:
        raise FloatingPointError('Sinusoidal RK4 state bound failure')
    return t, v, w, j.current(v, w)


def hysteresis(cfg):
    j = Junction(cfg); loops = []; records = []
    for f in np.logspace(0, 5, 9):
        # Shooting for the periodic orbit avoids comparing arbitrary transients.
        lo, hi = 0., 1.
        for _ in range(25):
            mid = (lo+hi)/2
            data = sinusoidal_cycle(j, f, mid)
            if data[2][-1] > mid: lo = mid
            else: hi = mid
        t, v, w, current = sinusoidal_cycle(j, f, (lo+hi)/2)
        half = (len(t)-1)//2
        # Sum unsigned lobe areas; a signed whole-cycle integral can cancel.
        area = abs(trapz(current[:half+1], v[:half+1]))+abs(trapz(current[half:], v[half:]))
        records.append([f, area, abs(w[-1]-w[0]), float(np.min(v*current)), w.min(), w.max()])
        loops.append(np.column_stack([np.full(len(t),f),t,v,w,current]))
    records = np.array(records)
    high_slope = float(np.polyfit(np.log10(records[-3:,0]), np.log10(records[-3:,1]), 1)[0])
    # A second resolution for an actual area convergence diagnostic.
    f = records[3,0]; initial = loops[3][0,3]
    t,v,w,i = sinusoidal_cycle(j,f,initial,1600); half=(len(t)-1)//2
    refined = abs(trapz(i[:half+1],v[:half+1]))+abs(trapz(i[half:],v[half:]))
    return loops, records, dict(high_frequency_log_slope=high_slope,
        area_refinement_relative_error=float(abs(refined-records[3,1])/refined),
        area_100khz_over_max=float(records[-1,1]/records[:,1].max()))


def spike(age):
    """A specified biphasic electrical waveform, not a fitted biological spike."""
    age = np.asarray(age)
    return np.where((age>=0)&(age<.001), .35,
                    np.where((age>=.001)&(age<.051), -.15*np.exp(-(age-.001)/.008), 0.))


def plasticity(cfg, dt=5e-5):
    j = Junction(cfg); delays = np.linspace(-.05,.05,101)
    times = np.arange(0,.23+dt/2,dt)
    w = np.full(len(delays),cfg.equilibrium)
    integrated = np.zeros_like(w); g0=float(j.conductance(cfg.equilibrium))
    peaks = np.zeros_like(w); last = np.zeros_like(w); trace=[]
    trace_index=int(np.argmin(abs(delays-.01)))
    for t in times[:-1]:
        # Positive voltage is post minus pre. Reversing electrodes reverses STDP.
        v = spike(t+dt/2-(.06+delays))-spike(t+dt/2-.06)
        new = j.advance(w,v,dt)
        integrated += .5*(j.conductance(w)+j.conductance(new)-2*g0)*dt
        peaks = np.maximum(peaks,np.abs(new-cfg.equilibrium))
        trace.append([t+dt,float(spike(t+dt/2-.06)),
            float(spike(t+dt/2-(.06+delays[trace_index]))),
            float(v[trace_index]),float(new[trace_index])])
        last = new; w = new
    normalized = integrated/(g0*cfg.tau_s)
    fits = []
    for sign in (1,-1):
        mask = (delays*sign>=.003-1e-12)&(delays*sign<=.04+1e-12)&(normalized*sign>1e-9)
        if np.sum(mask)>=3:
            x = np.abs(delays[mask]); y=np.log(abs(normalized[mask]))
            coeff = np.polyfit(x,y,1); pred=np.polyval(coeff,x)
            fits.append(dict(sign=sign,amplitude=float(np.exp(coeff[1])),
                decay_s=float(-1/coeff[0]) if coeff[0]<0 else None,
                log_space_r2=float(1-np.sum((y-pred)**2)/np.sum((y-y.mean())**2))))
        else: fits.append(dict(sign=sign,amplitude=None,decay_s=None,log_space_r2=None))
    # Two positive rectangular programming pulses, sampled at their ends.
    intervals=np.geomspace(.002,.1,35); q=np.full(len(intervals),cfg.equilibrium)
    first=np.zeros_like(q); second=np.zeros_like(q)
    for t in np.arange(0,.103,dt):
        mid=t+dt/2
        on=((mid>=0)&(mid<.001))|((mid>=intervals)&(mid<intervals+.001))
        q=j.advance(q,np.where(on,.35,0.),dt)
        if t<.001<=t+dt+1e-12: first[:]=j.conductance(q)
        pick=(t<intervals+.001)&(intervals+.001<=t+dt+1e-12)
        second[pick]=j.conductance(q)[pick]
    if np.any(first==0) or np.any(second==0):raise RuntimeError('PPF sampling failed')
    # PPF ratio refers to total small-bias conductance, not zero-baseline increments.
    return dict(delays=delays,area_siemens_s=integrated,normalized=normalized,
        final_w=last,peak_change=peaks,fits=fits,intervals=intervals,ppf=second/first,trace=np.array(trace),
        all_positive_delays_potentiate=bool(np.all(normalized[delays>0]>0)),
        all_negative_delays_depress=bool(np.all(normalized[delays<0]<0)))


class Mesh:
    """25 nonlinear shunts, 40 resistive edges, finite source R and node C.

    Backward-Euler nodal charge equation + operator-split RK4 molecular state.
    It is a 5x5 resistor mesh, NOT a full transistor-addressed crossbar layout.
    """
    def __init__(self,cfg):
        if cfg.nodes!=25:raise ValueError('The topology requires exactly 25 nodes')
        self.cfg=cfg; rng=np.random.default_rng(cfg.seed)
        self.scale=rng.uniform(.7,1.3,25)
        self.tau_scale=np.geomspace(.3,3.,25); rng.shuffle(self.tau_scale)
        self.rate_scale=rng.uniform(.65,1.35,25)
        self.j=Junction(cfg,self.scale,self.tau_scale,self.rate_scale)
        self.mask=rng.uniform(-1,1,25); self.bias=rng.uniform(-.07,.07,25)
        self.edges=[]; self.lap=np.zeros((25,25))
        for row in range(5):
            for col in range(5):
                n=5*row+col
                for m in ([n+1] if col<4 else [])+([n+5] if row<4 else []):
                    self.edges.append((n,m)); self.lap[n,n]+=1; self.lap[m,m]+=1
                    self.lap[n,m]-=1; self.lap[m,n]-=1

    def nodal(self,w,u,previous,dt):
        c=self.cfg; cap=c.node_c_f/dt
        base=c.wire_g_s*self.lap+np.eye(25)*(c.source_g_s+cap)
        rhs=c.source_g_s*u+cap*previous; v=previous.copy()
        for _ in range(15):
            residual=base@v+self.j.current(v,w)-rhs
            if np.max(abs(residual))<1e-13:break
            jac=base+np.diag(self.j.conductance(w)*np.cosh(c.alpha_per_v*v))
            v-=solve(jac,residual,assume_a='pos',check_finite=False)
        residual=base@v+self.j.current(v,w)-rhs
        if np.max(abs(residual))>=1e-12:raise RuntimeError('KCL Newton solve did not converge')
        return v,float(np.max(abs(residual)))

    def run(self,inputs,substeps=None,initial=.5):
        c=self.cfg; sub=c.substeps if substeps is None else substeps
        h=c.sample_s/sub; w=np.full(25,initial); v=np.zeros(25)
        states=[]; volts=[]; energy=[]; residuals=[]
        for value in inputs:
            u=self.bias+.27*self.mask*value
            totals=np.zeros(6); residual=0.
            for _ in range(sub):
                vn,r=self.nodal(w,u,v,h); residual=max(residual,r)
                # Energy balance is evaluated at the same held-state KCL solution.
                device=float(np.sum(vn*self.j.current(vn,w)))*h
                wires=c.wire_g_s*float(vn@self.lap@vn)*h
                sources=c.source_g_s*float(np.sum((u-vn)**2))*h
                supplied=c.source_g_s*float(np.sum(u*(u-vn)))*h
                dc=.5*c.node_c_f*float(np.sum(vn**2-v**2))
                be_loss=.5*c.node_c_f*float(np.sum((vn-v)**2))
                totals+=np.array([device,wires,sources,supplied,dc,be_loss])
                w=self.j.advance(w,vn,h); v=vn
            states.append(w.copy()); volts.append(v.copy()); energy.append(totals); residuals.append(residual)
        return np.array(states),np.array(volts),np.array(energy),np.array(residuals)


def lorenz(count,step):
    def f(q):
        x,y,z=q; return np.array([10*(y-x),x*(28-z)-y,x*y-(8/3)*z])
    q=np.array([1.,1.,1.]); out=[]; dt=step/4
    for n in range(count+1500):
        for _ in range(4):
            k1=f(q); k2=f(q+dt*k1/2); k3=f(q+dt*k2/2); k4=f(q+dt*k3)
            q+=dt*(k1+2*k2+2*k3+k4)/6
        if n>=1500:out.append(q.copy())
    return np.array(out)


def ridge_fit(x,y,lam):
    # Training-only centering, no polynomial, delay embedding, or hidden features.
    xm=x.mean(axis=0); xs=x.std(axis=0); xs=np.where(xs>1e-12,xs,1.)
    ym=y.mean(axis=0); u,s,vt=svd((x-xm)/xs,full_matrices=False,check_finite=False)
    coef=(vt.T*(s/(s*s+len(x)*lam)))@(u.T@(y-ym))
    physical=coef/xs[:,None]; intercept=ym-xm@physical
    return physical,intercept


def nmse(y,pred):
    return np.mean((y-pred)**2,axis=0)/np.var(y,axis=0)


def reservoir(cfg):
    signal=lorenz(cfg.samples+11,cfg.lorenz_step)
    avg=signal[cfg.washout:cfg.train_end,0].mean()
    scale=signal[cfg.washout:cfg.train_end,0].std()
    inputs=(signal[:cfg.samples,0]-avg)/scale
    mesh=Mesh(cfg); states,volts,energy,residual=mesh.run(inputs)
    train=np.arange(cfg.washout,cfg.train_end-10)
    valid=np.arange(cfg.train_end,cfg.validation_end-10)
    test=np.arange(cfg.validation_end,cfg.samples)
    # Ten-sample gaps prevent future labels leaking across split boundaries.
    predictions={}; scores={}; weights={}; ridge_trials=[]
    for horizon in (1,10):
        target=signal[np.arange(cfg.samples)+horizon]
        best=None
        for lam in np.logspace(-10,-2,9):
            coef,intercept=ridge_fit(states[train],target[train],lam)
            score=float(np.mean(nmse(target[valid],states[valid]@coef+intercept)))
            ridge_trials.append([horizon,lam,score])
            if best is None or score<best[0]:best=(score,float(lam),coef,intercept)
        _,lam,coef,intercept=best
        prediction=states@coef+intercept
        baseline_coef,baseline_intercept=ridge_fit(inputs[train,None],target[train],1e-10)
        input_only=inputs[:,None]@baseline_coef+baseline_intercept
        scores[str(horizon)]=dict(ridge_lambda=lam,validation_nmse=float(best[0]),
            test_nmse_xyz=nmse(target[test],prediction[test]).tolist(),
            test_mean_nmse=float(np.mean(nmse(target[test],prediction[test]))),
            input_only_nmse_xyz=nmse(target[test],input_only[test]).tolist(),
            persistence_x_nmse=float(nmse(target[test,0],signal[test,0])),
            below_001=bool(np.mean(nmse(target[test],prediction[test]))<.01))
        predictions[horizon]=(target,prediction,input_only)
        weights[str(horizon)]=dict(matrix=coef.tolist(),intercept=intercept.tolist())
    # Numerical resolution and fading-memory checks are diagnostics, not retraining.
    fine=mesh.run(inputs[:160],substeps=2*cfg.substeps)[0]
    convergence=float(np.max(abs(states[:160]-fine)))
    a=mesh.run(np.zeros(500),initial=.2)[0]; b=mesh.run(np.zeros(500),initial=.8)[0]
    distance=np.sqrt(np.mean((a-b)**2,axis=1))
    balance=energy[:,3]-np.sum(energy[:,:3],axis=1)-energy[:,4]-energy[:,5]
    circuit=float(np.mean(np.sum(energy[test,:3],axis=1)))
    return dict(signal=signal,inputs=inputs,states=states,volts=volts,energy=energy,
        residual=residual,predictions=predictions,weights=weights,scores=scores,
        ridge_trials=np.array(ridge_trials),test=test,mesh=mesh,fading=distance,
        diagnostics=dict(kcl_max_residual_a=float(residual.max()),
            energy_balance_max_error_j=float(abs(balance).max()),
            state_half_step_max_error=convergence,
            fading_distance_final=float(distance[-1]),
            state_min=float(states.min()),state_max=float(states.max())),
        energy_summary=dict(operation='one 25-node input sample, excluding readout and external electronics',
            test_device_j_per_sample=float(energy[test,0].mean()),
            test_wire_j_per_sample=float(energy[test,1].mean()),
            test_source_resistor_j_per_sample=float(energy[test,2].mean()),
            test_circuit_dissipation_j_per_sample=circuit,
            test_source_supplied_j_per_sample=float(energy[test,3].mean()),
            test_BE_artificial_loss_j_per_sample=float(energy[test,5].mean()),
            circuit_sub_fj=bool(circuit<1e-15),system_total_j_per_operation=None),
        split=dict(washout=cfg.washout,train_first=int(train[0]),train_last=int(train[-1]),
            validation_first=int(valid[0]),validation_last=int(valid[-1]),
            test_first=int(test[0]),test_last=int(test[-1]),
            input_train_mean=float(avg),input_train_std=float(scale)))


def equilibrium_state(j,bias):
    lo,hi=0.,1.
    for _ in range(60):
        mid=(lo+hi)/2
        if j.rhs(mid,bias)>0:lo=mid
        else:hi=mid
    return (lo+hi)/2


def impedance(cfg):
    j=Junction(cfg); frequencies=np.logspace(-2,6,481); omega=2*np.pi*frequencies
    rs=500.; capacitance=2e-9; rw=30e3; diffusion_s=.5
    root=np.sqrt(1j*omega*diffusion_s)
    # Finite-length transmissive Warburg: finite DC limit, 45-degree asymptote.
    # Stable complex tanh on Re(root)>0; avoid platform-specific overflow.
    decay=np.exp(-2*root)
    zw=rw*(-np.expm1(-2*root))/(1+decay)/root
    rows=[]; parameters=[]
    for bias in (-.25,0.,.25):
        w=equilibrium_state(j,bias); a=cfg.alpha_per_v
        iv=float(j.conductance(w)*np.cosh(a*bias))
        iw=(cfg.ghigh_s-cfg.glow_s)*np.sinh(a*bias)/a
        fw=-4*cfg.window_p*(2*w-1)**(2*cfg.window_p-1)*float(j.drive(bias))-1/cfg.tau_s
        fv=0. if abs(bias)<=cfg.threshold_v else float(j.eta)/cfg.v0_v*np.cosh((abs(bias)-cfg.threshold_v)/cfg.v0_v)*(1-(2*w-1)**(2*cfg.window_p))
        if fw>=0:raise RuntimeError('Unstable DC operating point cannot be used for stationary EIS')
        ystate=iv+iw*fv/(1j*omega-fw)
        z=rs+1/(1j*omega*capacitance+1/(1/ystate+zw))
        current_phasor=.01/z
        rows.append(np.column_stack([np.full(len(omega),bias),frequencies,z.real,z.imag,current_phasor.real,current_phasor.imag]))
        # Finite 10 mV symmetric change checks nonlinearity of the static I-V;
        # not a substitute for nonlinear AC harmonics or experimental calibration.
        plus=equilibrium_state(j,bias+.01); minus=equilibrium_state(j,bias-.01)
        finite=(j.current(bias+.01,plus)-j.current(bias-.01,minus))/.02
        linear=iv-iw*fv/fw
        parameters.append(dict(bias_v=bias,state=w,Rct_frozen_ohm=1/iv,
            state_relaxation_s=-1/fw,static_10mv_secant_relative_error=float(abs(finite-linear)/abs(linear))))
    return rows,dict(series_r_ohm=rs,geometric_c_f=capacitance,warburg_r_ohm=rw,
        diffusion_time_s=diffusion_s,ac_amplitude_v=.01,bias_points=parameters,
        scope='linearized molecular admittance plus assigned finite diffusion and capacitance; not fitted spectroscopy')


def figures(cfg,loops,hyst,plastic,rc,eis,destination):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'savefig.dpi':300,'figure.dpi':110})
    def save(fig,name):
        fig.savefig(destination/name,dpi=300,facecolor='white');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
    for k in (0,2,4,6,8):
        q=loops[k];axs[0].plot(q[:,2],q[:,4]*1e6,label=f'{q[0,0]:.0f} Hz')
    axs[0].set(xlabel='Voltage (V)',ylabel='Current (microampere)',title='A  Passive, periodic memristive junction')
    axs[0].axhline(0,color='.8',lw=.6);axs[0].axvline(0,color='.8',lw=.6);axs[0].legend(fontsize=8)
    axs[1].loglog(hyst[:,0],hyst[:,1],'-o')
    axs[1].set(xlabel='Frequency (Hz)',ylabel='Sum of unsigned lobe areas (A V)',title='B  Finite-frequency hysteresis collapse')
    save(fig,'fig1_molecular_pinched_hysteresis_loops.png')
    fig,axs=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
    x=plastic['delays'];y=plastic['normalized'];axs[0].plot(x*1e3,y,label='Integrated electrical response')
    for fit in plastic['fits']:
        if fit['decay_s'] is not None:
            xx=np.linspace(.003,.04,100);sg=fit['sign']
            axs[0].plot(sg*xx*1e3,sg*fit['amplitude']*np.exp(-xx/fit['decay_s']),'--',label=f"{sg:+d} tail fit: {fit['decay_s']*1e3:.1f} ms")
    axs[0].axhline(0,color='.6',lw=.6);axs[0].axvline(0,color='.6',lw=.6)
    axs[0].set(xlabel='Post minus pre delay (ms)',ylabel='Integral delta G dt / (G0 tau)',title='A  Volatile STDP-like response');axs[0].legend(fontsize=8)
    axs[1].semilogx(plastic['intervals']*1e3,plastic['ppf'],'-o',ms=3)
    axs[1].set(xlabel='Pulse onset separation (ms)',ylabel='G after pulse 2 / G after pulse 1',title='B  Paired-pulse facilitation')
    save(fig,'fig2_synaptic_stdp_hebbian_window.png')
    fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained');a=axs.ravel()
    test=rc['test'];section=test[:350];target,pred,_=rc['predictions'][1];t=section*cfg.lorenz_step
    for ax,dim in zip(a[:2],(0,2)):
        ax.plot(t,target[section,dim],color='black',label='Lorenz target')
        ax.plot(t,pred[section,dim],color='tab:orange',ls='--',label='Readout prediction')
        ax.set(xlabel='Lorenz time (dimensionless)',ylabel=('x','y','z')[dim],title=f"{'A' if dim==0 else 'B'}  Held-out, one-step {('x','y','z')[dim]} prediction")
        ax.legend(fontsize=8)
    image=a[2].imshow(rc['states'][section].T,aspect='auto',origin='lower',extent=[t[0],t[-1],0,25],vmin=0,vmax=1,cmap='viridis')
    a[2].set(xlabel='Lorenz time (dimensionless)',ylabel='Physical node',title='C  Simulated molecular states');fig.colorbar(image,ax=a[2],label='Redox state w')
    a[3].semilogy(np.arange(1,len(rc['fading'])+1)*cfg.sample_s*1e3,np.maximum(rc['fading'],1e-16))
    a[3].set(xlabel='Physical time (ms)',ylabel='RMS state separation',title='D  Fading memory under identical DC input')
    save(fig,'fig3_in_materio_chaotic_reservoir_tracking.png')
    fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    for row in eis:
        z=row[:,2]+1j*row[:,3];label=f'{row[0,0]:+.2f} V'
        axs[0].plot(z.real/1e3,-z.imag/1e3,label=label)
        axs[1].semilogx(row[:,1],np.angle(z,deg=True),label=label)
    axs[0].set(xlabel="Z real (kilohm)",ylabel='-Z imaginary (kilohm)',title='A  Linearized junction + assigned diffusion branch')
    axs[0].set_aspect('equal',adjustable='datalim');axs[0].legend()
    axs[1].set(xlabel='Frequency (Hz)',ylabel='Impedance phase (degrees)',title='B  10 mV nominal small-signal Bode phase');axs[1].legend()
    save(fig,'fig4_electrochemical_impedance_nyquist.png')


def reports(cfg,summary,root):
    s=summary;r=s['reservoir_scores'];e=s['energy'];p=s['plasticity'];h=s['hysteresis']
    common=f"""
| Quantity / 指标 | Computed value / 计算值 |
|---|---:|
| One-step test NMSE (mean x,y,z) | {r['1']['test_mean_nmse']:.8g} |
| One-step test NMSE x / y / z | {r['1']['test_nmse_xyz']} |
| Ten-step test NMSE (mean x,y,z) | {r['10']['test_mean_nmse']:.8g} |
| One-step persistence baseline, x only | {r['1']['persistence_x_nmse']:.8g} |
| Device dissipation J / input sample | {e['test_device_j_per_sample']:.8g} |
| Device + wire + source-resistor J / sample | {e['test_circuit_dissipation_j_per_sample']:.8g} |
| Full system J/op | Not calculated / 未计算 |
| High-frequency log(area)-log(f) slope | {h['high_frequency_log_slope']:.6g} |
| Hysteresis area refinement relative error | {h['area_refinement_relative_error']:.6g} |
| Positive / negative STDP sign tests | {p['all_positive_delays_potentiate']} / {p['all_negative_delays_depress']} |
| KCL max residual (A) | {s['diagnostics']['kcl_max_residual_a']:.6g} |
| Halved network-step max state difference | {s['diagnostics']['state_half_step_max_error']:.6g} |
| Halved plasticity-step max normalized-window difference | {p['half_step_max_window_difference']:.6g} |
"""
    refs="""
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
"""
    en=r"""# Phase 23 — In-materio computing: a reproducible effective-model benchmark

## Scope and results

This is a classical numerical simulation of a hypothetical volatile redox junction,
not a fabricated transition-metal bis-terpyridine/polyoxometalate device, a quantum
transport calculation, or experimental evidence of AGI. Coordination chemistry can
motivate tunable redox and ion-relaxation properties; the material-specific parameters
here are assigned. No atomistic electronic structure or solvent/temperature dependence
is computed. Physical reservoirs may reduce data movement for selected tasks; neither
an insurmountable silicon barrier nor a universal advantage over digital computers follows.
"""+common+r"""
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
"""+refs
    zh=r"""# 第23阶段：材料内计算与分子忆阻器——可复现的有效模型基准

## 定位与实际结果

本阶段是经典数值模拟，不是已制备的双三联吡啶/多金属氧酸盐器件，不是量子输运，
也不是实现 AGI 的证据。配位化学可以提供可调的氧化还原和离子弛豫性质，但本程序
的参数是指定值，未由某一真实配合物、第一性原理或实验拟合得到。
材料内计算有望减少部分任务的数据搬运；“硅路线不可逾越”及全面能效优势不能据此推出。
"""+common+r"""
## 23A：为什么必须修改原公式

原式 I=G(w)V sinh(alpha V) 是电压的偶函数，负电压时 I V<0，不能用作这里的被动耗散器件。
采用 I=G(w)sinh(alpha V)/alpha，确保小电压电导为 G，且任意电压 I V>=0。
状态方程为：

    q = sign(V) max(|V|-V_threshold,0)
    dw/dt = eta [1-(2w-1)^(2p)] sinh(q/V0) - (w-w_eq)/tau

保留 Joglekar 窗口，p=1、w_eq=0.5、tau=10 ms。正反指数反应速率相减产生 sinh，
所以这是受 Butler-Volmer 启发的模型，而非定量推导的界面动力学；电流和状态变化
没有通过法拉第电荷守恒关联，V0 是有效参数。没有计算热涨落、反应自由能或化学热。
原模型在 w=0 会锁死；向内部平衡态恢复后，两个边界的导数都指向区间内部。
RK4 根据速率上界细分步长，越界时报错，不用裁剪状态掩盖数值不稳定。

频率覆盖 1 Hz—100 kHz，通过周期轨道射击法排除任意初态瞬变。
报告两瓣面积绝对值之和，避免有向积分相互抵消。高频时状态摆幅为 O(1/f)，
趋于冻结状态的单值 I-V 曲线；实际输出有限频率斜率和加密网格误差。
过原点的滞回与忆阻动力学相容，但不是鉴定真实忆阻器的唯一实验证据。

## 23B：脉冲塑性与长期保持的边界

每个脉冲由 +0.35 V、1 ms 的头部，以及初始 -0.15 V、8 ms 衰减常数的尾部构成，
尾部在 51 ms 截断。器件电压定义为 V_post-V_pre，交换电极会交换学习窗口方向。
扫描 -50 至 +50 ms 共101个延迟，统一积分同一个器件方程，不直接植入正负指数学习规则。

Delta W=int[G_pair(t)-G0]dt 的量纲为 S·s，不是电导；图上无量纲量为 Delta W/(G0 tau)，
不能误标为 Delta G/G0。统一 230 ms 窗口包含所有脉冲及至少69 ms的最终恢复。
3—40 ms尾部进行对数空间指数拟合，并保存拟合质量；这不是“实验 Hebbian 保真度”。
另外导出最终状态和最大状态偏移。PPF 定义为第二脉冲后电导/第一脉冲后电导，
不是扣除基线后的增量之比。零偏压下状态恢复到0.5，因此只能称挥发性 STP/STDP-like，
不能声称长期 LTP/LTD；需要额外慢变量、非易失态和保持实验才能研究长期学习。

## 23C：25节点物理网络与预测协议

5×5节点通过40条电阻边连接，每个节点有接地非线性忆阻器、有限输入电阻和2 pF寄生电容。
固定随机输入增益、偏置、电导和速率，以及3—30 ms弛豫时间代表指定的器件差异。
输入是仅用训练段标准化的 Lorenz x，经过固定电压掩码分配到节点，不加入延迟嵌入、
多项式或手工特征。拓扑是带类交叉阵列寄生项的电阻网格，不是完整行列选通交叉阵列。

每个子步用 Newton 法和正定线性求解器求解：

    (C/h)(v-v_previous)+G_source(v-u)+G_wire L v+I(v,w)=0

然后在保持的节点电压下推进分子状态。整体是后向Euler电容加算子分裂，
全局精度不是四阶。导出KCL残差和子步减半差异。两个初始状态在相同直流输入下
收敛用于展示衰退记忆，但不是任意输入下的全局回声状态稳定性证明。

Lorenz参数为10、28、8/3，每0.02无量纲时间映射为0.5 ms物理采样，是指定时钟，
不是实测计算速度。丢弃最初1500个生成样本后保留3600个输入；洗去200个样本，
训练截止2200，验证截止2800，最后800个样本完全留作测试。
边界预留10个样本，防止未来标签跨分区泄漏。标准化、读出中心化和缩放只用训练集。
只在验证集比较9个岭回归系数；测试集不参与器件参数、读出或超参数选择。

通过SVD训练一个25×3线性矩阵及偏置，预测未来1步和10步的x、y、z。
NMSE分别按测试段每个坐标的方差归一化，同时报告三坐标平均值。
另提供只使用当前输入的线性基线和x坐标的持续性基线。
这是持续输入真实x的驱动式预测，不是自主闭环长期生成混沌。
单步x预测容易因信号平滑而获得低误差，不等于三维或较长步长的预测性能。
本次为确定性单种子结果，未计算跨器件噪声或多个独立混沌样本的置信区间。

## 能耗：必须说明每次操作是什么

一次操作定义为整个25节点网络接收一个输入样本，不随意拆成内部事件以缩小数值。
分别积分器件VI耗散、线阻损耗和输入电阻损耗，并保存电源提供能量、电容储能变化、
后向Euler人工数值损耗，验证同一离散电路的能量平衡。
右端点求积与KCL中保持状态一致；离散守恒残差不是连续时间能耗精确性的证明。
尚未计入ADC、DAC、掩码分配、读出、训练、复位、噪声裕度、制造及化学反应热，
所以报告的是电路子总量，系统总能耗为空；不得标成已实现亚飞焦耳系统。
精度和能耗是否达到原定阈值均由程序计算并保存，未通过也保留。

## 23D：阻抗分析模型

分别求 -0.25、0、+0.25 V下的稳定平衡态，对状态方程F和电流线性化：

    Y_j = I_V + I_w F_V/(i omega-F_w)
    Z_W = R_W tanh(sqrt(i omega tau_D))/sqrt(i omega tau_D)
    Z = R_s + [i omega C_geom + (1/Y_j + Z_W)^(-1)]^(-1)

明确扫描的是f=0.01—1e6 Hz，角频率为2πf。保存名义10 mV交流扰动的电流相量。
Warburg采用有限长度、透射边界，低频电阻有限，并非阻塞边界下发散的赝电容。
R_s=500欧姆、C_geom=2 nF、R_W=30千欧、tau_D=0.5 s都是指定的薄膜/电解质参数，
与网络的单节点寄生电容不同。单一红氧状态无法产生连续扩散谱，故Warburg是额外
等效电路假设，不是从前述状态ODE直接推导出来的扩散。
冻结状态R_ct=1/I_V与完整直流微分电阻不同。输出稳定性检查和10 mV静态割线误差，
并未完成非线性交流谐波仿真、参数反演或实验谱峰解卷积。图谱也不能验证配位环境。
偏压下的状态响应项可能产生感性符号的小环，不能把所有弧线都解释为纯电容半圆。

## 实验路线与复现

下一步应选择明确组成的氧化还原薄膜，测量脉冲、保持、温度依赖和电化学状态，
独立测定线阻/电容，再开展含噪声的化学传感边缘计算测试。
分子电子传感接口是值得验证的方向；亚飞焦耳AGI协处理器仍是展望，不是本阶段成果。

安装 requirements_phase23.txt，运行脚本的 --self-test，再正常运行。
支持 --output-dir 指定输出目录；脚本本身不联网、不执行Git、不自动改写README。
交付包含全部数据、岭回归验证记录、读出矩阵、数据划分、摘要、阈值判断、SHA-256校验
及四张300 DPI图。提交和发布由独立的人工可审查流程完成。
"""+refs
    for lang,text in [('EN',en),('ZH',zh)]:
        (root/f'IN_MATERIO_COMPUTING_REPORT_{lang}.md').write_text(text,encoding='utf-8',newline='\n')


class Tests(unittest.TestCase):
    def test_passivity_and_pinch(self):
        j=Junction(Config());v=np.linspace(-.8,.8,201)
        for w in (0,.5,1):
            self.assertTrue(np.all(v*j.current(v,w)>=0));self.assertEqual(j.current(0,w),0)
    def test_bounds_and_relaxation(self):
        c=Config();j=Junction(c);w=np.array([0.,.2,.8,1.])
        self.assertTrue(np.all(j.rhs(0.,np.array([-.8,.8]))>0))
        self.assertTrue(np.all(j.rhs(1.,np.array([-.8,.8]))<0))
        final=j.advance(w,0.,.03)
        np.testing.assert_allclose(final,.5+(w-.5)*np.exp(-3),atol=5e-7)
        for v in (-.8,.8):
            q=j.advance(w,v,.01);self.assertTrue(np.all((q>=0)&(q<=1)))
    def test_mesh_energy_and_kcl(self):
        m=Mesh(Config());state,v,e,r=m.run(np.array([0.,1.,-1.,.2]))
        self.assertLess(r.max(),1e-12)
        np.testing.assert_allclose(e[:,3],e[:,:3].sum(axis=1)+e[:,4]+e[:,5],atol=1e-16,rtol=1e-7)
        self.assertTrue(np.all(e[:,:3]>=-1e-20));self.assertEqual(len(m.edges),40)
    def test_mesh_zero_control(self):
        m=Mesh(Config());m.bias[:]=0
        s,v,e,r=m.run(np.zeros(5))
        np.testing.assert_allclose(s,.5);np.testing.assert_allclose(e,0,atol=1e-25)
    def test_ridge_known_map(self):
        rng=np.random.default_rng(2);x=rng.normal(size=(100,4));a=rng.normal(size=(4,3));y=x@a+2
        b,c=ridge_fit(x,y,1e-12);np.testing.assert_allclose(x@b+c,y,atol=1e-9)
    def test_impedance(self):
        rows,p=impedance(Config())
        for row in rows:
            self.assertTrue(np.isfinite(row).all());self.assertTrue(np.all(row[:,2]>0))
            self.assertLess(abs(row[-1,2]-p['series_r_ohm']),1.)
    def test_reproducibility(self):
        m=Mesh(Config());n=Mesh(Config());np.testing.assert_array_equal(m.mask,n.mask)
        np.testing.assert_array_equal(lorenz(5,.02),lorenz(5,.02))
    def test_plasticity_zero_delay_and_signs(self):
        p=plasticity(Config())
        self.assertAlmostEqual(float(p['normalized'][50]),0.,places=12)
        self.assertTrue(p['all_positive_delays_potentiate'])
        self.assertTrue(p['all_negative_delays_depress'])
        self.assertGreater(float(p['ppf'][0]),1.)
    def test_nodal_symmetry_and_dc(self):
        m=Mesh(Config());w=np.full(25,.5);u=np.full(25,.2)
        v,r=m.nodal(w,u,np.zeros(25),.001)
        self.assertLess(r,1e-12)
        self.assertTrue(np.all(v>0));self.assertTrue(np.all(v<.2))
        np.testing.assert_allclose(m.lap@np.ones(25),0,atol=1e-14)
    def test_quadrature(self):
        self.assertAlmostEqual(float(trapz(np.array([0.,1.,2.]),np.array([0.,1.,2.]))),2.)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
    args=parser.parse_args(argv)
    if args.self_test:
        result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests))
        return 0 if result.wasSuccessful() else 1
    cfg=Config();root=args.output_dir.resolve();res=root/'results_phase23';fig=root/'figures_phase23'
    res.mkdir(parents=True,exist_ok=True);fig.mkdir(parents=True,exist_ok=True)
    print('23A: computing periodic hysteresis and refinement check',flush=True)
    loops,hr,hd=hysteresis(cfg)
    print('23B: propagating paired-pulse plasticity and half-step comparison',flush=True)
    plastic=plasticity(cfg);plastic_fine=plasticity(cfg,dt=2.5e-5)
    print('23C: solving 25-node KCL reservoir and held-out ridge prediction',flush=True);rc=reservoir(cfg)
    print('23D: linearized impedance and figure generation',flush=True);eis,ep=impedance(cfg)
    ps={k:plastic[k] for k in ('fits','all_positive_delays_potentiate','all_negative_delays_depress')}
    ps['half_step_max_window_difference']=float(np.max(abs(plastic['normalized']-plastic_fine['normalized'])))
    ps['half_step_max_ppf_difference']=float(np.max(abs(plastic['ppf']-plastic_fine['ppf'])))
    if rc['diagnostics']['state_half_step_max_error']>.002:
        raise RuntimeError('Network state refinement exceeds declared 0.002 tolerance')
    if ps['half_step_max_window_difference']>.001:
        raise RuntimeError('STDP refinement exceeds declared 0.001 tolerance')
    summary=dict(parameters=asdict(cfg),hysteresis=hd,plasticity=ps,reservoir_scores=rc['scores'],
        energy=rc['energy_summary'],diagnostics=rc['diagnostics'],split=rc['split'],eis=ep,
        provenance='Assigned classical volatile junction model; synthetic Lorenz data; no experimental device or AGI claim',
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__))
    def csv(name,data,header):
        with (res/name).open('w',encoding='utf-8',newline='\n') as f:
            np.savetxt(f,np.asarray(data),delimiter=',',header=header,comments='',fmt='%.12g')
    csv('hysteresis.csv',np.concatenate(loops),'frequency_hz,time_s,voltage_v,state,current_a')
    csv('hysteresis_area.csv',hr,'frequency_hz,unsigned_area_a_v,periodic_state_residual,min_power_w,min_state,max_state')
    csv('stdp.csv',np.column_stack([plastic[k] for k in ('delays','area_siemens_s','normalized','final_w','peak_change')]),'delay_s,delta_weight_siemens_s,delta_weight_over_g0_tau,final_state,peak_state_change')
    csv('ppf.csv',np.column_stack([plastic['intervals'],plastic['ppf']]),'interval_s,second_over_first_conductance')
    csv('paired_spike_trace.csv',plastic['trace'],'state_time_s,pre_midpoint_v,post_midpoint_v,junction_midpoint_v,state')
    time=np.arange(cfg.samples)*cfg.sample_s
    csv('reservoir_states.csv',np.column_stack([time,rc['states']]),'time_s,'+','.join(f'w{k}' for k in range(25)))
    csv('node_voltages.csv',np.column_stack([time,rc['volts']]),'time_s,'+','.join(f'v{k}' for k in range(25)))
    csv('energy.csv',np.column_stack([time,rc['energy'],rc['residual']]),'time_s,device_j,wire_j,source_resistor_j,supplied_j,capacitor_change_j,BE_numerical_loss_j,kcl_residual_a')
    for horizon,(target,pred,baseline) in rc['predictions'].items():
        split=np.full(cfg.samples,-1)
        split[cfg.washout:cfg.train_end-10]=0
        split[cfg.train_end:cfg.validation_end-10]=1
        split[cfg.validation_end:]=2
        csv(f'prediction_h{horizon}.csv',np.column_stack([time,rc['inputs'],split,target,pred,baseline]),'time_s,input_normalized,partition_minus1unused_0train_1val_2test,target_x,target_y,target_z,pred_x,pred_y,pred_z,input_only_x,input_only_y,input_only_z')
    csv('ridge_validation.csv',rc['ridge_trials'],'horizon_samples,lambda,validation_mean_nmse')
    csv('fading_memory.csv',np.column_stack([np.arange(1,501)*cfg.sample_s,rc['fading']]),'time_s,rms_state_distance')
    m=rc['mesh'];csv('mesh_parameters.csv',np.column_stack([np.arange(25),m.mask,m.bias,m.scale,m.tau_scale*cfg.tau_s,m.rate_scale*cfg.eta_per_s]),'node,input_mask,bias_v,conductance_scale,tau_s,eta_per_s')
    csv('mesh_edges.csv',m.edges,'node_a,node_b')
    csv('impedance.csv',np.concatenate(eis),'bias_v,frequency_hz,z_real_ohm,z_imag_ohm,i_real_a,i_imag_a')
    for name,data in [('summary.json',summary),('readout.json',rc['weights'])]:
        (res/name).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    figures(cfg,loops,hr,plastic,rc,eis,fig);reports(cfg,summary,root)
    paths=sorted([*res.glob('*.csv'),res/'summary.json',res/'readout.json',*fig.glob('*.png'),*root.glob('IN_MATERIO_COMPUTING_REPORT_*.md')])
    hashes={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    hashes['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (res/'sha256.json').write_text(json.dumps(hashes,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(scores=rc['scores'],energy=rc['energy_summary'],diagnostics=rc['diagnostics'],plasticity=ps),indent=2),flush=True)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
