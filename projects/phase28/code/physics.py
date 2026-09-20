"""Conservative, deliberately uncalibrated Phase 28 scenario equations."""
from __future__ import annotations
import numpy as np

FARADAY_C_MOL = 96485.33212
CL2_KG_MOL = 0.07090


def simulate(load, dt=1.0, kir=0.003, kdes=0.0015, clearance=0.12):
    """State fractions [anchored, irreversible, soluble, removed], initial total=1.

    Rates have h^-1 at unit load; load is dimensionless j/(0.3 A cm^-2).
    dA/dt=-(kir+kdes)load^2 A; dD/dt=kir load^2 A;
    dS/dt=kdes load^2 A-clearance*S; W receives clearance*S.
    Exact interval propagation conserves the molecular inventory. Removed material
    remains in W; S is a normalized inventory, NOT an experimentally measured mol/L.
    """
    load = np.asarray(load, dtype=float)
    if np.any(load < 0) or min(dt, kir, kdes, clearance) < 0 or dt == 0:
        raise ValueError("Nonnegative currents/rates and positive dt required")
    states = np.zeros((len(load) + 1, 4)); states[0, 0] = 1
    for n, x in enumerate(load):
        a, d, s, _ = states[n]
        k = (kir + kdes) * x*x
        anew = a*np.exp(-k*dt)
        dnew = d + (kir/(kir+kdes)*(a-anew) if kir+kdes else 0)
        if abs(clearance-k) < 1e-12:
            integral = dt*np.exp(-k*dt)
        else:
            integral = (np.exp(-k*dt)-np.exp(-clearance*dt))/(clearance-k)
        snew = s*np.exp(-clearance*dt) + kdes*x*x*a*integral
        states[n+1] = (anew, dnew, snew, 1-anew-dnew-snew)
    return states


def observations(load, kir=0.003, kdes=0.0015, clearance=0.12):
    states = simulate(load, kir=kir, kdes=kdes, clearance=clearance)
    a = states[1:, 0]
    voltage = 3.0 + 0.12*np.asarray(load) + 0.8*(1-a)
    return np.column_stack((voltage, a, states[1:, 2]))


def candidate_observations(name, load, rate, amplitude=0.8):
    """Minimal distinguishable observational candidates, not established chemistry.

    H1 fixed-site irreversible loss; H2 soluble species additionally active;
    H3 carrier damage with stable organic surface; H4 reversible transfer lag.
    """
    load = np.asarray(load)
    if name in ("H1", "H3"):
        states = simulate(load, kir=rate, kdes=0)
        a = states[1:, 0]
        v = 3 + 0.12*load + amplitude*(1-a)
        return np.column_stack((v, a if name == "H1" else np.ones(len(a)), np.zeros(len(a))))
    if name == "H2":
        states = simulate(load, kir=0, kdes=rate)
        a, s = states[1:, 0], states[1:, 2]
        return np.column_stack((3+0.12*load+amplitude*(1-a-2*s), a, s))
    if name == "H4":
        m = 0.0; out = []
        for x in load:
            m = x*x + (m-x*x)*np.exp(-rate)
            out.append((3+0.12*x+amplitude*m, 1.0, 0.0))
        return np.asarray(out)
    raise ValueError(name)


def product_kg(load, dt_h=1, area_cm2=100, j_ref=0.3, fe=0.985):
    if not 0 <= fe <= 1: raise ValueError("FE must be in [0,1]")
    return np.asarray(load)*j_ref*area_cm2*dt_h*3600*fe*CL2_KG_MOL/(2*FARADAY_C_MOL)


def kaplan_meier(times, events):
    """Right-censored survival estimator; tied failures precede censor removals."""
    times, events = np.asarray(times), np.asarray(events, bool)
    if len(times) != len(events) or np.any(times < 0): raise ValueError("Invalid survival records")
    risk, survival, rows = len(times), 1.0, []
    for t in np.unique(times):
        deaths = int(np.sum((times == t) & events))
        censored = int(np.sum((times == t) & ~events))
        survival *= 1-deaths/risk
        rows.append((float(t), risk, deaths, censored, survival))
        risk -= deaths+censored
    return rows


def no_batch_leakage(records):
    seen = {}
    for record in records:
        batch, split = record["batch_id"], record["split"]
        if batch in seen and seen[batch] != split: raise ValueError("Batch leakage")
        seen[batch] = split
    return True
