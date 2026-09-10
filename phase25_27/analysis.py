"""Numerical analysis primitives; production callers must supply validated data.

Free energies in kcal/mol for solution kinetics; eV for electrodes.
ee = (R-S)/(R+S); effective ddG = G_S^eff-G_R^eff = RT ln(R/S).
An effective ddG inferred from yields is not automatically a TS difference.
"""
from dataclasses import dataclass
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import logsumexp
from scipy.stats import norm, t as student_t

R_KCAL=0.00198720425864083
KB_EV=8.617333262145e-5
KB_OVER_H=2.0836619123327574e10


def positive_temperature(T):
    if not np.isfinite(T) or T<=0: raise ValueError('Temperature must be finite and positive')


def ensemble_free_energy(energies, degeneracies, T=298.15):
    positive_temperature(T)
    g=np.asarray(energies,dtype=float); d=np.asarray(degeneracies,dtype=float)
    if g.ndim!=1 or not len(g) or g.shape!=d.shape or not np.isfinite(g).all() or not np.isfinite(d).all() or np.any(d<=0):
        raise ValueError('Finite energies and positive physical degeneracies required')
    # Count each physically distinct basin once; optimizer hit counts are NOT d.
    lw=np.log(d)-g/(R_KCAL*T); z=logsumexp(lw)
    return float(-R_KCAL*T*z),np.exp(lw-z)


def signed_selectivity(R,S,T=298.15):
    positive_temperature(T)
    if not np.isfinite([R,S]).all() or min(R,S)<0: raise ValueError('Nonnegative finite yields required')
    if R+S==0: return dict(ee=None,effective_ddG_kcal=None,major='selectivity_undetermined')
    ee=(R-S)/(R+S)
    ddg=None if min(R,S)==0 else R_KCAL*T*math.log(R/S)
    return dict(ee=ee,effective_ddG_kcal=ddg,major='R' if ee>0 else ('S' if ee<0 else 'racemic'))


def selectivity_interval(ee_draws):
    x=np.asarray(ee_draws,dtype=float)
    if x.ndim!=1 or len(x)<2 or not np.isfinite(x).all() or np.max(np.abs(x))>1: raise ValueError('Finite ee draws in [-1,1] required')
    lo,mid,hi=np.quantile(x,[0.025,0.5,0.975])
    return dict(median=float(mid),CI95=[float(lo),float(hi)],major='R' if lo>0 else ('S' if hi<0 else 'selectivity_undetermined'))


def correlated_energy_draws(mean,covariance,seed,n=2000):
    mean=np.asarray(mean,dtype=float); cov=np.asarray(covariance,dtype=float)
    if cov.shape!=(len(mean),len(mean)) or not np.isfinite(mean).all() or not np.isfinite(cov).all() or not np.allclose(cov,cov.T): raise ValueError('Invalid covariance')
    if np.linalg.eigvalsh(cov).min()<-1e-10: raise ValueError('Covariance must be positive semidefinite')
    return np.random.default_rng(seed).multivariate_normal(mean,cov,size=n,check_valid='raise')


def standard_state_correction(T=298.15,pressure_Pa=101325,concentration_mol_L=1):
    positive_temperature(T)
    if pressure_Pa<=0 or concentration_mol_L<=0: raise ValueError('Positive standard state required')
    return R_KCAL*T*math.log(concentration_mol_L*1000*8.314462618*T/pressure_Pa)


@dataclass
class Network:
    """Closed homogeneous mass-action network with reversible shared TS energies.

    nu_r/nu_p: reactions x species. All standard state G use the same moiety reference.
    TS G is for the complete composition of each elementary reaction.
    Mixing entropy arises from activities; it must not also be added to G0.
    No implicit BEP law, barrier clipping, steady-state assumption or fast-exchange reduction.
    """
    species: list
    nu_r: np.ndarray
    nu_p: np.ndarray
    G0: np.ndarray
    Gts: np.ndarray
    T: float=298.15
    c0: float=1.0

    def __post_init__(self):
        positive_temperature(self.T)
        self.nu_r=np.asarray(self.nu_r,float); self.nu_p=np.asarray(self.nu_p,float)
        self.G0=np.asarray(self.G0,float); self.Gts=np.asarray(self.Gts,float)
        if self.nu_r.shape!=self.nu_p.shape or self.nu_r.ndim!=2 or self.nu_r.shape[1]!=len(self.species): raise ValueError('Invalid stoichiometry shape')
        if self.G0.shape!=(len(self.species),) or self.Gts.shape!=(len(self.nu_r),): raise ValueError('Missing free energies')
        if not np.isfinite(self.c0) or self.c0<=0 or not all(np.isfinite(x).all() for x in [self.nu_r,self.nu_p,self.G0,self.Gts]): raise ValueError('Missing/invalid production data')
        if np.any(self.nu_r<0) or np.any(self.nu_p<0): raise ValueError('Negative stoichiometry')
        self.S=(self.nu_p-self.nu_r).T
        self.barriers_forward=self.Gts-self.nu_r@self.G0
        self.barriers_reverse=self.Gts-self.nu_p@self.G0
        if min(self.barriers_forward.min(),self.barriers_reverse.min())<0: raise ValueError('Submerged TS requires an explicit capture/diffusion model')
        self.kf=KB_OVER_H*self.T*np.exp(-self.barriers_forward/(R_KCAL*self.T))
        self.kr=KB_OVER_H*self.T*np.exp(-self.barriers_reverse/(R_KCAL*self.T))

    def flux(self,c):
        c=np.asarray(c,float)
        if c.shape!=(len(self.species),) or not np.isfinite(c).all() or np.any(c < -1e-10): raise ValueError('Invalid concentrations')
        a=np.maximum(c,0)/self.c0
        return self.c0*(self.kf*np.prod(a[None,:]**self.nu_r,axis=1)-self.kr*np.prod(a[None,:]**self.nu_p,axis=1))

    def integrate(self,c_initial,times,conservation_vectors):
        ts=np.asarray(times,float); c=np.asarray(c_initial,float); L=np.asarray(conservation_vectors,float)
        if ts.ndim!=1 or len(ts)<2 or not np.isfinite(ts).all() or not np.isfinite(c).all() or np.any(np.diff(ts)<=0) or ts[0]<0 or np.any(c<0): raise ValueError('Invalid integration schedule')
        if not np.allclose(L@self.S,0,atol=1e-12): raise ValueError('Network violates declared conservation laws')
        sol=solve_ivp(lambda time,y:self.S@self.flux(y),(ts[0],ts[-1]),c,method='BDF',t_eval=ts,rtol=1e-8,atol=1e-12)
        if not sol.success: raise RuntimeError(sol.message)
        if np.min(sol.y)<-1e-9: raise RuntimeError('Negative population')
        if not np.allclose(L@sol.y,(L@c)[:,None],atol=1e-8,rtol=1e-7): raise RuntimeError('Conservation drift')
        return sol


def curtin_hammett_diagnostic(Q,reaction_hazards,threshold=100):
    """Conservative spectral mixing/escape separation, not an unconditional proof.

    Q is a connected reversible conformer CTMC generator with ROW sums zero.
    Irreducibility/reversibility and concentration-dependent hazards need physical validation.
    Repeat across temperatures, concentrations and uncertain energy draws.
    """
    Q=np.asarray(Q,float); h=np.asarray(reaction_hazards,float)
    if Q.shape!=(len(h),len(h)) or len(h)<2 or not np.isfinite(Q).all() or not np.isfinite(h).all() or np.any(h<0): raise ValueError('Invalid generator')
    off=Q-np.diag(np.diag(Q))
    if np.any(off<-1e-12) or not np.allclose(Q.sum(1),0): raise ValueError('Invalid generator')
    ev=np.linalg.eigvals(Q)
    if np.max(np.abs(ev.imag))>1e-8 or np.sum(np.abs(ev)<1e-10)!=1: return dict(supported=False,reason='disconnected or non-reversible spectrum')
    A=Q.T.copy(); A[-1,:]=1; b=np.zeros(len(h));b[-1]=1
    pi=np.linalg.solve(A,b)
    stationary_flux=pi[:,None]*Q
    if np.any(pi<=0) or not np.allclose(stationary_flux,stationary_flux.T,rtol=1e-7,atol=max(1e-12,abs(Q).max()*1e-10)):
        return dict(supported=False,reason='stationary process fails detailed balance')
    gap=float(np.min(-ev.real[ev.real<-1e-10])); hazard=float(h.max())
    ratio=None if hazard==0 else gap/hazard
    return dict(mixing_rate_s=gap,max_escape_rate_s=hazard,separation_ratio=ratio,supported=hazard>0 and ratio>=threshold,
                limitation='spectral necessary diagnostic; compare full transient network before CH reduction')


def electrode_grand_energy(G_eV,N_e,N_e_neutral,mu_e_eV):
    """Omega=G-mu_e*(N_e-N0); same nuclear/ionic composition and gauge required.

    Variable atom counts additionally require reservoir chemical potentials.
    Positive cell charge in units of e equals N0-Ne, not Ne-N0.
    """
    values=np.asarray([G_eV,N_e,N_e_neutral,mu_e_eV],float)
    if not np.isfinite(values).all(): raise ValueError('Missing path-point data')
    return float(G_eV-mu_e_eV*(N_e-N_e_neutral))


def audit_potential_path(points,target_U,tolerance=0.05):
    if len(points)<3: raise ValueError('Need reactant, transition region and product')
    required=['G_eV','N_e','N_e_neutral','mu_e_eV','actual_U_V_SHE','total_cell_charge_e','electrode_charge_partition_e']
    for point in points:
        if any(k not in point or point[k] is None or not np.isfinite(point[k]) for k in required): raise ValueError('Missing constant-potential evidence')
        if abs(point['actual_U_V_SHE']-target_U)>tolerance: raise ValueError('Actual potential out of tolerance')
        if abs(point['total_cell_charge_e']-(point['N_e_neutral']-point['N_e']))>1e-6: raise ValueError('Charge/electron-number sign mismatch')
    omega=[electrode_grand_energy(p['G_eV'],p['N_e'],p['N_e_neutral'],p['mu_e_eV']) for p in points]
    return dict(Omega_eV=omega,max_potential_deviation_V=max(abs(p['actual_U_V_SHE']-target_U) for p in points),
        free_energy_barrier_eV=None,note='pointwise energy audit is not an explicit-solvent activation PMF')


def independent_mean_ci(values):
    x=np.asarray(values,float)
    if len(x)<3 or not np.isfinite(x).all(): raise ValueError('At least 3 independent estimates needed')
    half=float(student_t.ppf(0.975,len(x)-1)*x.std(ddof=1)/math.sqrt(len(x)))
    return dict(mean=float(x.mean()),CI95=[float(x.mean()-half),float(x.mean()+half)],halfwidth=half,n=len(x))


def competing_risk_incidence(times,events):
    """Aalen-Johansen competing-risk cumulative incidence with right censoring.

    Each trajectory contributes one final *absorbing* event or 'censored'.
    First bond rupture and recombination are intermediate states, not necessarily terminal events.
    Independent/non-informative censoring required; no asymptotic quantum-yield inference.
    """
    times=np.asarray(times,float); events=np.asarray(events,str)
    if times.shape!=events.shape or not len(times) or not np.isfinite(times).all() or np.any(times<0): raise ValueError('Invalid trajectories')
    causes=sorted(set(events)-{'censored'}); F={k:0. for k in causes}; S=1.; history=[]
    for time in np.unique(times):
        risk=int(np.sum(times>=time)); mask=times==time
        d={k:int(np.sum(mask&(events==k))) for k in causes}
        previous=S
        for k in causes: F[k]+=previous*d[k]/risk
        S*=1-sum(d.values())/risk
        history.append(dict(time=float(time),survival=float(S),CIF=dict(F),at_risk=risk,censored=int(np.sum(mask&(events=='censored')))))
    return history


def wilson_interval(successes,total,confidence=.95):
    if total<1 or successes<0 or successes>total or int(total)!=total or int(successes)!=successes: raise ValueError('Invalid binomial counts')
    z=norm.ppf((1+confidence)/2); p=successes/total; denom=1+z*z/total
    center=(p+z*z/(2*total))/denom; half=z*math.sqrt(p*(1-p)/total+z*z/(4*total*total))/denom
    return dict(CI=[float(center-half),float(center+half)],halfwidth=float(half),n=int(total))


def validate_stationary_path(record):
    """Strict metadata prefilter. It does not replace inspection of raw QM output."""
    fields=['opt_surface','freq_surface','irc_surface','imaginary_frequencies_cm1','imaginary_mode_is_reaction','irc_forward_endpoint','irc_reverse_endpoint','expected_endpoints','raw_log_sha256','geometry_sha256']
    if any(k not in record or record[k] is None for k in fields): raise ValueError('Incomplete path evidence')
    if len({record[k] for k in ['opt_surface','freq_surface','irc_surface']})!=1: raise ValueError('Inconsistent potential energy surfaces')
    if len(record['imaginary_frequencies_cm1'])!=1 or record['imaginary_frequencies_cm1'][0]>=0 or record['imaginary_mode_is_reaction'] is not True: raise ValueError('Not a verified first-order reaction saddle')
    if sorted([record['irc_forward_endpoint'],record['irc_reverse_endpoint']])!=sorted(record['expected_endpoints']): raise ValueError('Wrong IRC endpoints')
    return True
