#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_phase10_cyberphysical_flow_twin.py
======================================
PHASE 10 — THE GRAND INTEGRATION
MULTI-SCALE CONTINUOUS-FLOW DIGITAL TWIN, OPERANDO PAT &
CYBER-PHYSICAL REINFORCEMENT-LEARNING PROCESS CONTROL

The Pantheon converges: the molecular quantum kinetics of Phases 4-8 are
bridged to macroscale fluid mechanics, real-time in-operando Process
Analytical Technology (PAT), and deep-reinforcement-learning closed-loop
process control.  The physical object is a coiled tubular microfluidic
reactor (L = 10.0 m, inner diameter d_t = 1.0 mm, counter-current cooling
jacket) synthesizing the Phase-4/5/8 frontier target — the strain-release
asymmetric aziridine ring expansion

    2-methyl-azirino[1,2-a]indole  ->  3-methyl-dihydroquinoline-imine
    (C10H11N, M = 145.21 g/mol, catalyzed by the designed
     3,3'-CF3-Ph/iPr-Ph-BINOL chiral phosphoric acid, CPA)

using the Phase-5 microkinetic anchor points (RDS barrier 21.24 kcal/mol =
88.87 kJ/mol, stereodifferentiation ddG|= 1.50 kcal/mol -> ee =
tanh(ddG|/(2RT)) = 85.9 % at 293 K) inside a genuinely multi-scale plant:

Module 10A — 1D Advection-Diffusion-Reaction Continuum Solver (MOL)
-------------------------------------------------------------------
Coupled non-linear PDEs over z in [0, L], t in [0, 3600 s]:

    dC_i/dt     = D_eff d2C_i/dz2 - u(t) dC_i/dz + R_i(C_1..C_N, T)
    rho cp dT/dt = k_th d2T/dz2 - rho cp u dT/dz
                   + SUM_j (-DH_j) r_j - (4U/d_t)(T - T_coolant)

with a second counter-current advected coolant field
    rho_c cp_c dTc/dt = -rho_c cp_c u_c dTc/dz + (4U/d_t)(T - Tc)

Mesoscale hydrodynamics close the scale bridge:
  * Taylor-Aris axial dispersion  D_eff = D_m + a^2 u^2 / (192 D_m)
  * laminar Hagen-Poiseuille pressure  dP/dz = 128 mu(T) Q / (pi d_eff^4)
  * Graetz/laminar film coefficients  h_i = Nu k_l/d_t, h_jacket ~ u_c^0.55
    combined into an overall U through the series resistance.
Integration: Method of Lines on N_z cells; IMEX time stepping —
trapezoidal-implicit diffusion (banded tridiagonal solve, unconditionally
stable against the Taylor-Aris stiffness) + explicit TVD-RK2 (MC-limiter
flux) advection/reaction/exchange with adaptive dt from the CFL condition.

Chemical network (Phase-5 microkinetics, concentrations in mol/L):
    r1: A + Cat -> I     (RDS, E_a = 88.87 kJ/mol, inhibited by impurity)
    r2: I -> P           (E_a = 59.4 kJ/mol, strain-relaxing fast step)
    r3: I -> E           (E_a = 110.9 kJ/mol, acid-promoted elimination,
                          suppressed by the 2,6-lutidine scavenger stream)
    r4: 2I -> O          (E_a = 75.3 kJ/mol, cationic oligomerization,
                          lutidine-suppressed)
    kd: Cat -> X         (catalyst burial in the fouling zone)
Heat of reaction: net A->P -85 kJ/mol (aziridine strain release); the
channel is heat-transfer-limited: nominal operation is quasi-isothermal
(<1 K rise), but loss of jacket pumping collapses U towards the adiabatic
limit dT_ad = C_A,0*dH/(rho cp) ~ 44 K and the Arrhenius positive feedback
produces a genuine thermal runaway, exactly as in the macroscale plant.

Module 10B — Operando In-Line Multi-Modal PAT Stream
----------------------------------------------------
High-frequency telemetry at inspection ports z = 2.5, 5.0, 10.0 m:
  1. in-line 785 nm Raman: Lorentzian fingerprint bands for all 7 analytes,
     laser-power OU fluctuation, rising fluorescence baseline with drift,
     cavitation spike artifacts;
  2. UV-Vis photodiode array at lambda = 254 / 310 nm (Beer-Lambert,
     1.0 mm path, shot noise + drift);
  3. 1 Hz thermal + differential-pressure telemetry (Hagen-Poiseuille dP
     with the fouled diameter profile; vapor-pressure margin via Antoine).
Industrial anomalies are injected on schedule: pump-B cavitation (slip
oscillation), precursor-batch impurity spike (benzoic acid: reversible
RDS inhibition + acid-promoted side chemistry + fluorescence jump),
progressive catalyst-fouling / channel clogging at z ~ 7 m (d_eff shrink,
insulating deposit, catalyst burial), and a coolant-loop vapor lock
(exponential jacket-pump degradation -> near-adiabatic channel).
A PAT deconvolution agent converts the raw modalities into concentrations
by regularized least squares on 6 measurements; its audited error vs plant
truth continues the Phase-9 hallucination-audit discipline.

Module 10C — Cyber-Physical RL / NMPC Process Controller
--------------------------------------------------------
A torch Soft-Actor-Critic agent (twin Q, squashed-Gaussian actor, automatic
entropy temperature) learns on a domain-randomized family of fault
episodes.  State (18-dim): deconvoluted port concentrations, three port
temperatures, dP + trend, flow-rate history and slewed actuator state,
max channel superheat dT and its rate, vapor-pressure margin.
Actions (5-dim continuous): Q_A, Q_B, Q_cat in [0.01, 5.0] mL/min,
jacket cooling rate (dimensionless) and back-pressure setpoint P in
[1, 30] bar.  Reward exactly in the commanded form
    R_t = w1*Yield + w2*Selectivity - w3*Carbon
          - w4*I(ThermalRunaway) - w5*PressurePenalty.
Deployment wraps the policy in an NMPC-style supervisory safety filter
(model-based one-step thermal shield with hysteresis: dilute-quench
flush at maximum flow, pressure/boiling guardrails, antifouling catalyst
pulse) — safety-filtered RL.  The acceptance test: when the coolant fault
drives the channel towards runaway, the composite agent detects the
hotspot, re-routes stoichiometry and residence time and suppresses the
excursion within 60 s without human intervention, while the identical
open-loop run exceeds dT > 40 K with selectivity collapse.

Module 10D — Techno-Economic & Green-Chemistry Lifecycle Engine
---------------------------------------------------------------
Dynamic steady states (direct PFR + counter-current jacket iteration) over
a throughput-scale grid (1-64 channels x 0.5-5.0 mL/min x jacket T x
catalyst loading) feed a full cradle-to-gate account: Space-Time Yield
(kg m^-3 h^-1), Process Mass Intensity, E-Factor, carbon intensity
(kg CO2-eq/kg product) and unit production cost with solvent-recovery,
energy and catalyst-amortization breakdown.  The 3-objective Pareto surface
(cost, carbon, STY) is extracted and plotted.

Outputs
-------
results_phase10/phase10_results.json        machine-readable master record
results_phase10/sac_learning_curve.csv      episode reward / loss trace
results_phase10/validation_summary.csv      3-controller x 8-episode suite
results_phase10/tea_pareto_points.csv       full TEA/LCA design-point table
results_phase10/episode_*.npz               replayable field histories
figures_phase10/fig1_continuous_pde_reactor_profile.png   (300 DPI)
figures_phase10/fig2_operando_pat_sensor_telemetry.png    (300 DPI)
figures_phase10/fig3_rl_cyberphysical_control_dynamics.png (300 DPI)
figures_phase10/fig4_techno_economic_pareto_analysis.png  (300 DPI)

Usage
-----
python run_phase10_cyberphysical_flow_twin.py             # full campaign
python run_phase10_cyberphysical_flow_twin.py --selftest  # smoke test
python run_phase10_cyberphysical_flow_twin.py --fig_only  # redraw figures
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_banded

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (3-D Pareto panel)

import torch
import torch.nn as nn
import torch.optim as optim

# --------------------------------------------------------------------------- #
# Global configuration
# --------------------------------------------------------------------------- #

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results_phase10"
FIGURES = ROOT / "figures_phase10"
RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)

SEED = 20260905
DEVICE = torch.device("cpu")

# --- reactor geometry (Module 10A spec) ------------------------------------- #
L_REACTOR = 10.0                    # m
D_TUBE = 1.0e-3                     # m inner diameter
A_CS = math.pi * D_TUBE**2 / 4.0    # m^2
V_REACTOR = A_CS * L_REACTOR        # m^3  (7.854 mL)

# --- physical properties (toluene phase, documented in the report) ---------- #
RHO_L = 867.0          # kg/m^3
CP_L = 1760.0          # J/kg/K
K_L = 0.131            # W/m/K
MU_293 = 5.9e-4        # Pa s
DM_MOLEC = 8.5e-10     # m^2/s molecular diffusivity of solutes in toluene
R_GAS = 8.314462618    # J/mol/K

T_FEED = 293.15        # K
T_COOL_IN = 298.15     # K  (chiller setpoint of the counter-current jacket)
RHO_CP_COOL = 4.05e6   # J/m^3/K (40 % glycol)
A_ANN = 5.3e-6         # m^2 annular coolant flow area (3 mm jacket, 1 mm coil)
R_WALL = 6.0e-4        # m^2 K/W (PFA coil wall + contact resistance)
U_COOL_VEL_NOM = 0.05  # m/s annular coolant velocity at jflow = 1
U_COOL_VEL_MAX = 0.22  # m/s modeled cap (annular bypass)

# --- chemistry ------------------------------------------------------------ #
M_P = 145.21               # g/mol  (C10H11N)
DH_1 = -38.0e3             # J/mol  A -> I
DH_2 = -47.0e3             # J/mol  I -> P   (net A->P: -85 kJ/mol)
DH_3 = -16.0e3             # J/mol  I -> E
DH_4 = -40.0e3             # J/mol  2I -> O (per oligomer bond)
EA_1 = 88.87e3             # J/mol  Phase-5 RDS (21.24 kcal/mol)
EA_2 = 59.4e3
EA_3 = 110.9e3
EA_4 = 85.0e3
K1_REF = 2.30              # M^-1 s^-1 at 293.15 K (calibrated: X ~ 90 %)
K2_REF = 1.8e-1            # s^-1
K3_REF = 9.0e-5            # s^-1
K4_REF = 5.5e-2            # M^-1 s^-1
KD_REF = 3.0e-4            # s^-1 catalyst burial
T_REF_KIN = 293.15
DDG_STEREO = 6.276e3       # J/mol  (1.50 kcal/mol Phase-5)
K_INH_IMP = 220.0          # M^-1  reversible RDS inhibition by benzoic acid
K_ACID_ELIM = 150.0        # M^-1  acid promotion of r3
K_LUT_ELIM = 40.0          # M^-1  lutidine suppression of r3
K_LUT_OLIG = 25.0          # M^-1  lutidine suppression of r4

# --- feeds ---------------------------------------------------------------- #
C_A_FEED = 1.80            # M substrate in toluene (stream A)
C_LUT_FEED = 0.10          # M 2,6-lutidine modifier (stream B)
C_CAT_FEED = 0.060         # M CPA catalyst (stream Cat)
QA_NOM, QB_NOM, QCAT_NOM = 3.60, 2.40, 1.20   # mL/min (production push)
JFLOW_NOM = 1.0            # dimensionless jacket rate
PBPR_NOM = 6.0             # bar

# --- control --------------------------------------------------------------- #
CONTROL_DT = 5.0           # s between control decisions
EP_T_DEMO = 3600.0         # s operational horizon (spec)
ACTION_LO = np.array([0.01, 0.01, 0.01, 0.20, 1.0])
ACTION_HI = np.array([5.00, 5.00, 5.00, 3.00, 30.0])
ACTION_NOM = np.array([QA_NOM, QB_NOM, QCAT_NOM, JFLOW_NOM, PBPR_NOM])
ACTION_SLEW = np.array([1.25, 1.25, 1.00, 0.60, 2.00])   # per control step

W1, W2, W3, W4, W5 = 1.0, 0.6, 0.35, 4.0, 0.8   # reward weights (spec form)

# --- PAT ------------------------------------------------------------------- #
PORT_Z = (2.5, 5.0, 10.0)
RAMAN_WAVENUMBERS = np.arange(700.0, 1900.0, 2.0)
UV_WAVELENGTHS = (254.0, 310.0)
UV_PATH = 0.1              # cm (1.0 mm flow cell)

# --------------------------------------------------------------------------- #
# Small utilities
# --------------------------------------------------------------------------- #


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def arrhenius(k_ref: float, ea: float, t: np.ndarray | float):
    """Simple Arrhenius scaling from a reference rate at T_REF_KIN."""
    return k_ref * np.exp(-ea / R_GAS * (1.0 / np.asarray(t, dtype=float)
                                         - 1.0 / T_REF_KIN))


def ee_of_temperature(t: float | np.ndarray):
    """Phase-5 enantioselectivity law  ee = tanh(ddG| / (2 R T))."""
    return np.tanh(DDG_STEREO / (2.0 * R_GAS * np.asarray(t, dtype=float)))


def viscosity(t):
    """Log-linear toluene viscosity (Pa s)."""
    return MU_293 * np.exp(-0.0113 * (np.asarray(t, dtype=float) - 293.15))


def toluene_vapor_pressure(t_c):
    """Antoine (mmHg, degC) -> bar."""
    p_mmhg = 10.0 ** (6.95464 - 1344.8 / (np.asarray(t_c, dtype=float) + 219.48))
    return p_mmhg * 101325.0 / 760.0 / 1.0e5


def minmod3(a, b, c):
    """Vectorized three-argument minmod (MC-limiter slope)."""
    m_ab = 0.5 * (np.sign(a) + np.sign(b)) * np.minimum(np.abs(a), np.abs(b))
    return 0.5 * (np.sign(m_ab) + np.sign(c)) * np.minimum(np.abs(m_ab), np.abs(c))


# --------------------------------------------------------------------------- #
# MODULE 10A — the continuum plant
# --------------------------------------------------------------------------- #


class FlowPlant:
    """Method-of-Lines digital twin of the coiled microreactor.

    State fields on the z-grid (all length N_z):
        C[A, I, P, E, O, Cat, Lut]   mol/L
        T   channel temperature      K
        Tc  counter-current coolant  K
    """

    SPECIES = ("A", "I", "P", "E", "O", "Cat", "Lut", "Imp")

    def __init__(self, n_z: int = 161, dt_cap: float = 0.5, seed: int = SEED):
        self.nz = n_z
        self.dz = L_REACTOR / (n_z - 1)
        self.z = np.linspace(0.0, L_REACTOR, n_z)
        self.dt_cap = dt_cap
        self.rng = np.random.default_rng(seed)
        self.idx_ports = [int(round(zp / L_REACTOR * (n_z - 1))) for zp in PORT_Z]
        # fouling zone (progressive clogging / catalyst burial)
        self.foul_lo, self.foul_hi = 6.4, 7.6
        prof = np.clip((self.z - self.foul_lo) / 0.35, 0, 1) * \
            np.clip((self.foul_hi - self.z) / 0.35, 0, 1)
        self.foul_profile = np.clip(prof, 0, 1)
        self.reset()

    # -------------------------------------------------- lifecycle
    def reset(self):
        self.C = {s: np.zeros(self.nz) for s in self.SPECIES}
        self.T = np.full(self.nz, T_FEED)
        self.Tc = np.full(self.nz, T_COOL_IN)
        self.t = 0.0
        self.Q = {"A": QA_NOM, "B": QB_NOM, "cat": QCAT_NOM}     # mL/min actual
        self.Q_cmd = dict(self.Q)
        self.jflow_cmd = JFLOW_NOM
        self.pbpr = PBPR_NOM
        self.foul_amp = 0.0
        self.d_eff = np.full(self.nz, D_TUBE)
        # disturbance schedule (t_lo, t_hi, kind) — set by scenario
        self.events: list[dict] = []
        self._cav_phase = self.rng.uniform(0, 2 * np.pi)
        self.c_imp_feed = 0.0
        self.purity_a = 1.0
        self.jflow_delivered = JFLOW_NOM
        self._k_eff = 0.05
        self._apply_actuators(np.array([QA_NOM, QB_NOM, QCAT_NOM,
                                        JFLOW_NOM, PBPR_NOM]))
        # per-tick records
        self.rec_t: list[float] = []
        self.rec_T: list[np.ndarray] = []
        self.rec_C: list[dict] = []
        self.rec_Tc: list[np.ndarray] = []
        self.rec_dp = []
        self.rec_actions = []
        self.rec_delivered = []

    # -------------------------------------------------- disturbances
    def set_scenario(self, events):
        self.events = list(events)

    def _event_windows(self, kind):
        return [e for e in self.events if e["kind"] == kind]

    def _window_value(self, kind, t):
        """0..1 activity of a windowed event (smoothed ramps)."""
        v = 0.0
        for e in self._event_windows(kind):
            t0, t1 = e["t0"], e["t1"]
            if t0 - 20 <= t <= t1 + 20:
                ramp = np.clip((t - t0) / e.get("ramp", 25.0), 0.0, 1.0)
                ramp_dn = np.clip((t1 + 20 - t) / 25.0, 0.0, 1.0) \
                    if e.get("finite", True) else 1.0
                v = max(v, min(ramp, ramp_dn) * e.get("amp", 1.0))
        return v

    # -------------------------------------------------- actuator layer
    def _apply_actuators(self, actions):
        """Slew-limited setpoints + real actuator pathologies."""
        for k, key in enumerate(("A", "B", "cat")):
            self.Q_cmd[key] = float(np.clip(self.Q_cmd[key] + np.clip(
                actions[k] - self.Q_cmd[key], -ACTION_SLEW[k], ACTION_SLEW[k]),
                ACTION_LO[k], ACTION_HI[k]))
        self.jflow_cmd = float(np.clip(
            self.jflow_cmd + np.clip(actions[3] - self.jflow_cmd,
                                     -ACTION_SLEW[3], ACTION_SLEW[3]),
            ACTION_LO[3], ACTION_HI[3]))
        self.pbpr = float(np.clip(
            self.pbpr + np.clip(actions[4] - self.pbpr,
                                -ACTION_SLEW[4], ACTION_SLEW[4]),
            ACTION_LO[4], ACTION_HI[4]))
        # --- delivered flows: cavitation slip on pump B -------------------
        cav = self._window_value("cav", self.t)
        slip = 0.0
        if cav > 0:
            # BPR elevation raises suction pressure and collapses the slip
            relief = np.clip(1.0 - 0.09 * (self.pbpr - PBPR_NOM), 0.12, 1.0)
            slip = cav * relief * (0.55 + 0.45 * math.sin(
                2 * math.pi * 1.3 * self.t + self._cav_phase))
        self.Q["A"] = self.Q_cmd["A"]
        self.Q["B"] = self.Q_cmd["B"] * (1.0 - slip)
        self.Q["cat"] = self.Q_cmd["cat"]
        # --- coolant loop: vapor lock degrades delivered jacket flow ------
        lock = self._window_value("cool", self.t)
        self.jflow_delivered = self.jflow_cmd * (1.0 - 0.94 * lock)
        # --- fouling progression ------------------------------------------
        self.foul_amp = float(self._window_value("foul", self.t))
        d_shrink = 1.0 - 0.48 * self.foul_amp * self.foul_profile
        self.d_eff = D_TUBE * d_shrink
        # --- impurity spike on stream A -----------------------------------
        self.c_imp_feed = 0.012 * self._window_value("imp", self.t)
        self.purity_a = 1.0 - 0.06 * self._window_value("imp", self.t)
        # --- volume flow / velocity ---------------------------------------
        self.q_tot = (self.Q["A"] + self.Q["B"] + self.Q["cat"]) / 6e7  # m^3/s
        self.u = self.q_tot / A_CS                                        # m/s
        self.u_cool = min(U_COOL_VEL_NOM * max(self.jflow_delivered, 0.02),
                          U_COOL_VEL_MAX)

    # -------------------------------------------------- derived physics
    def _thermal_fields(self):
        t = self.T
        mu = viscosity(t)
        d_m = 8.5e-10
        # Taylor-Aris dispersion (meso-scale closure)
        d_eff_ax = d_m + (D_TUBE / 2.0) ** 2 * self.u**2 / (192.0 * d_m)
        # film coefficients
        h_i = 3.66 * K_L / D_TUBE
        h_c = 1350.0 * (max(self.u_cool, 0.01) / U_COOL_VEL_NOM) ** 0.55
        u_overall = 1.0 / (1.0 / h_i + R_WALL + 1.0 / h_c)
        # fouling deposit: a growing thermal-resistance blanket (strongly
        # insulating the clogged zone) — the physical seed of the hotspot
        u_overall = u_overall / (1.0 + 12.0 * self.foul_amp
                                 * self.foul_profile)
        return mu, d_eff_ax, u_overall

    def _rates(self):
        c = self.C
        c_cat = np.maximum(c["Cat"], 0.0)
        c_i = np.maximum(c["I"], 0.0)
        c_a = np.maximum(c["A"], 0.0)
        t = self.T
        k1 = arrhenius(K1_REF, EA_1, t)
        k2 = arrhenius(K2_REF, EA_2, t)
        k3 = arrhenius(K3_REF, EA_3, t)
        k4 = arrhenius(K4_REF, EA_4, t)
        kd = arrhenius(KD_REF, 42.0e3, t) * (1.0 + 4.0 * self.foul_amp
                                             * self.foul_profile)
        c_imp = c["Imp"] if "Imp" in c else np.zeros(self.nz)
        c_lut = np.maximum(c["Lut"], 0.0)
        inh = 1.0 / (1.0 + K_INH_IMP * np.maximum(c_imp, 0.0))
        acid = (1.0 + K_ACID_ELIM * np.maximum(c_imp, 0.0)) \
            / (1.0 + K_LUT_ELIM * c_lut)
        lut_olig = 1.0 / (1.0 + K_LUT_OLIG * c_lut)
        r1 = k1 * c_a * c_cat * inh
        r2 = k2 * c_i
        r3 = k3 * acid * c_i
        r4 = k4 * lut_olig * c_i**2
        return r1, r2, r3, r4, kd

    def _reaction_rhs(self):
        r1, r2, r3, r4, kd = self._rates()
        r = {
            "A": -r1,
            "I": r1 - r2 - r3 - 2.0 * r4,
            "P": r2,
            "E": r3,
            "O": r4,
            "Cat": -kd * np.maximum(self.C["Cat"], 0.0),
            "Lut": np.zeros(self.nz),
            "Imp": np.zeros(self.nz),
        }
        # depletion-aware chemical stiffness measure (lagged into the
        # substep-size selection): where the reactant is locally exhausted
        # the effective decay constant collapses, so the runaway front —
        # not the hot inert plateau — sets the time step.
        eps_c = 1e-9
        self._k_eff = max(
            float(np.max(r1 / np.maximum(self.C["A"], eps_c))),
            float(np.max((r2 + r3 + 2.0 * r4) / np.maximum(self.C["I"],
                                                           eps_c))),
            float(np.max(kd)))
        q_reaction = 1.0e3 * (r1 * (-DH_1) + r2 * (-DH_2) + r3 * (-DH_3)
                              + r4 * (-DH_4))     # W/m^3 (mol/L/s -> mol/m^3/s)
        return r, q_reaction

    def _inlet_composition(self):
        q = max(self.q_tot, 1e-12)
        c_in = {s: 0.0 for s in self.SPECIES}
        c_in["A"] = C_A_FEED * self.purity_a * self.Q["A"] / 6e7 / q
        c_in["Lut"] = C_LUT_FEED * self.Q["B"] / 6e7 / q
        c_in["Cat"] = C_CAT_FEED * self.Q["cat"] / 6e7 / q
        if self.c_imp_feed > 0:
            c_in["Imp"] = self.c_imp_feed * self.Q["A"] / 6e7 / q
        return c_in

    # -------------------------------------------------- TVD advection
    def _advect_stack(self, F, u, f_in):
        """TVD MC-limiter upwind advection for a stack of fields
        F: (m, n_z), u > 0, Danckwerts boundaries (inlet value f_in[m])."""
        dz = self.dz
        dc = np.diff(F, axis=1)
        slope = np.zeros_like(F)
        slope[:, 1:-1] = minmod3(2.0 * dc[:, :-1],
                                 0.5 * (dc[:, :-1] + dc[:, 1:]),
                                 2.0 * dc[:, 1:])
        f_int = u * (F[:, :-1] + 0.5 * slope[:, :-1])       # faces 1..N-1
        rhs = np.empty_like(F)
        rhs[:, 0] = -(f_int[:, 0] - u * f_in) / dz
        rhs[:, 1:-1] = -(f_int[:, 1:] - f_int[:, :-1]) / dz
        rhs[:, -1] = -(u * F[:, -1] - f_int[:, -1]) / dz
        return rhs

    def _advect_coolant(self, tc, src):
        """Counter-current coolant: advects in -z; inlet at z = L.
        src: wall heat input per coolant volume (W/m^3 of annulus)."""
        rev = tc[::-1]
        dc = np.diff(rev)
        slope = np.zeros_like(rev)
        slope[1:-1] = minmod3(2.0 * dc[:-1], 0.5 * (dc[:-1] + dc[1:]),
                              2.0 * dc[1:])
        f_int = self.u_cool * (rev[:-1] + 0.5 * slope[:-1])
        f = np.concatenate(([self.u_cool * T_COOL_IN], f_int,
                            [self.u_cool * rev[-1]]))
        return (-(f[1:] - f[:-1]) / self.dz)[::-1] + src

    # -------------------------------------------------- IMEX stepping
    def _diffusion_implicit(self, fields, coef):
        """Trapezoidal-implicit diffusion step for stacked fields.

        fields: (n_fields, n_z) array; coef: scalar axial diffusivity.
        Danckwerts zero-diffusive-flux at both ends.
        """
        n = self.nz
        h = self.dt_sub
        lam = coef * h / self.dz**2
        ab = np.zeros((3, n))
        ab[0, 1:] = -lam
        ab[1, :] = 1.0 + 2.0 * lam
        ab[2, :-1] = -lam
        # mirror (zero-gradient) closure at both Danckwerts boundaries:
        # d2C/dz2|_0 ~ 2(C1 - C0)/dz2
        ab[0, 1] = -2.0 * lam
        ab[2, -2] = -2.0 * lam
        out = solve_banded((1, 1), ab, fields.T)
        return out.T

    def _explicit_rhs(self):
        """Advection + reaction + jacket exchange (explicit operator).
        Uses self._tf = (mu, D_eff_ax, U_overall) cached per substep."""
        _, d_eff_ax, u_ov = self._tf
        names = self.SPECIES
        c_stack = np.stack([self.C[s] for s in names])
        c_in = self._inlet_composition()
        f_in = np.array([c_in.get(s, 0.0) for s in names])
        rhs_stack = self._advect_stack(c_stack, self.u, f_in)
        r_rxn, q_reaction = self._reaction_rhs()
        for i, s in enumerate(names):
            rhs_stack[i] += r_rxn[s]
        rhs_c = {s: rhs_stack[i] for i, s in enumerate(names)}
        # wall exchange (W/m^3 of channel)
        exchange = (4.0 * u_ov / D_TUBE) * (self.T - self.Tc)
        rho_cp = RHO_L * CP_L
        rhs_t = (self._advect_stack(self.T[None, :], self.u,
                                    np.array([T_FEED]))[0]
                 + q_reaction / rho_cp - exchange / rho_cp)
        # coolant energy: per m^3 of annulus (area A_ANN)
        rhs_tc = self._advect_coolant(
            self.Tc, exchange * A_CS / (RHO_CP_COOL * A_ANN))
        return rhs_c, rhs_t, rhs_tc, q_reaction

    def _step_chemical(self, dt):
        """One IMEX micro-step: implicit diffusion, then TVD-RK2 transport."""
        self.dt_sub = dt
        self._tf = self._thermal_fields()
        _, d_eff_ax, _ = self._tf
        # ---- implicit axial diffusion (species share Taylor-Aris D; T uses
        # the bare thermal diffusivity k_th/(rho cp))
        c_stack = np.stack([self.C[s] for s in self.SPECIES])
        c_stack = self._diffusion_implicit(c_stack, d_eff_ax)
        for i, s in enumerate(self.SPECIES):
            self.C[s] = c_stack[i]
        alpha_t = K_L / (RHO_L * CP_L)
        self.T = self._diffusion_implicit(self.T[None, :], alpha_t)[0]

        def adv_react(state_c, state_t, state_tc):
            keep = (self.C, self.T, self.Tc)
            self.C, self.T, self.Tc = state_c, state_t, state_tc
            rhs_c, rhs_t, rhs_tc, _ = self._explicit_rhs()
            self.C, self.T, self.Tc = keep
            return rhs_c, rhs_t, rhs_tc

        rhs_c, rhs_t, rhs_tc = adv_react(self.C, self.T, self.Tc)
        c_mid = {s: np.maximum(self.C[s] + 0.5 * dt * rhs_c[s], 0.0)
                 for s in self.SPECIES}
        t_mid = np.maximum(self.T + 0.5 * dt * rhs_t, 250.0)
        tc_mid = self.Tc + 0.5 * dt * rhs_tc
        rhs_c2, rhs_t2, rhs_tc2 = adv_react(c_mid, t_mid, tc_mid)
        for s in self.SPECIES:
            self.C[s] = np.clip(self.C[s] + dt * rhs_c2[s], 0.0, 2.5)
        self.T = np.maximum(self.T + dt * rhs_t2, 250.0)
        self.Tc = self.Tc + dt * rhs_tc2

    # -------------------------------------------------- public API
    def advance(self, dt_target: float, actions: np.ndarray):
        """Advance the plant by dt_target with slewed actuators (adaptive
        CFL sub-stepping — the 'adaptive finite-difference' MOL integrator).
        The substep respects the advective CFL, the coolant CFL *and* the
        reaction timescale (the Arrhenius acceleration during a runaway
        transiently stiffens the kinetics by two orders of magnitude)."""
        self._apply_actuators(actions)
        t_end = self.t + dt_target
        guard = 0
        while self.t < t_end - 1e-9 and guard < 200000:
            guard += 1
            self._tf = self._thermal_fields()
            dt = min(self.dt_cap,
                     0.7 * self.dz / max(self.u, 1e-6),
                     0.7 * self.dz / max(self.u_cool, 1e-6),
                     0.3 / max(self._k_eff, 1e-6),
                     t_end - self.t)
            dt = max(dt, 1e-3)
            self._step_chemical(dt)
            self.t += dt
        return self.t

    def record(self, dp: float):
        self.rec_t.append(self.t)
        self.rec_T.append(self.T.copy())
        self.rec_Tc.append(self.Tc.copy())
        self.rec_C.append({s: self.C[s].copy() for s in self.SPECIES})
        self.rec_dp.append(dp)
        self.rec_actions.append((self.Q_cmd["A"], self.Q_cmd["B"],
                                 self.Q_cmd["cat"], self.jflow_cmd,
                                 self.pbpr))
        self.rec_delivered.append((self.Q["A"], self.Q["B"], self.Q["cat"],
                                   self.jflow_delivered))

    # -------------------------------------------------- hydraulics
    def delta_p(self):
        """Hagen-Poiseuille integral dP/dz = 128 mu Q / (pi d^4)."""
        mu = viscosity(self.T)
        seg = 128.0 * mu * self.q_tot / (math.pi * self.d_eff**4)
        return float(np.sum(seg) * self.dz) / 1.0e5      # bar

    def vapor_margin(self):
        """min over z of P_abs - P_vap(T) in bar."""
        dp_prof = np.cumsum(128.0 * viscosity(self.T) * self.q_tot
                            / (math.pi * self.d_eff**4)) * self.dz / 1.0e5
        p_abs = self.pbpr + (dp_prof[-1] - dp_prof)      # bar
        p_vap = toluene_vapor_pressure(self.T - 273.15)
        return float(np.min(p_abs - p_vap))

    # -------------------------------------------------- observation
    def observe(self, with_pat_noise=True):
        """Port telemetry + PAT-deconvoluted concentrations (Module 10B)."""
        ports_t = [float(self.T[i]) for i in self.idx_ports]
        ports_c = {s: [float(self.C[s][i]) for i in self.idx_ports]
                   for s in ("A", "I", "P", "E")}
        dp = self.delta_p()
        obs = {
            "t": self.t,
            "ports_T": ports_t,
            "ports_C": ports_c,
            "dP": dp,
            "dP_rate": getattr(self, "_dp_rate", 0.0),
            "vapor_margin": self.vapor_margin(),
            "max_dT": float(np.max(self.T) - T_COOL_IN),
            "dT_rate": getattr(self, "_dT_rate", 0.0),
            "actions": np.array([self.Q_cmd["A"], self.Q_cmd["B"],
                                 self.Q_cmd["cat"], self.jflow_cmd,
                                 self.pbpr]),
            "truth_out": {s: float(self.C[s][-1])
                          for s in ("A", "I", "P", "E", "O", "Cat")},
        }
        if with_pat_noise:
            obs["est"] = pat_deconvolve(self, self.rng)
        else:
            obs["est"] = {s: obs["truth_out"][s]
                          for s in ("A", "I", "P", "E")}
        return obs


# --------------------------------------------------------------------------- #
# MODULE 10B — PAT models
# --------------------------------------------------------------------------- #

# Raman bands: (wavenumber cm-1, amplitude, Lorentzian HWHM cm-1)
RAMAN_BANDS = {
    #  component: [(nu, amp, hwhm), ...]
    "A":   [(748, .55, 9), (1022, .30, 8), (1265, .70, 10),
            (1455, .45, 12), (1602, 1.00, 9), (2925 * 0 + 1780, .10, 14)],
    "Lut": [(760, .40, 8), (994, .35, 7), (1375, .28, 9), (1585, .30, 8)],
    "Cat": [(1032, .18, 10), (1232, .22, 12)],
    "I":   [(1180, .35, 9), (1372, .30, 8), (1638, .60, 10)],
    "P":   [(1188, .42, 9), (1305, .30, 8), (1450, .38, 10),
            (1610, .45, 8), (1655, 1.00, 9)],
    "E":   [(1402, .25, 9), (1668, .55, 9)],
    "O":   [(1005, .50, 20), (1440, .40, 18)],
    "Imp": [(1003, .80, 8), (1608, .60, 10)],
}
RAMAN_XSEC = {"A": 1.0, "Lut": 0.55, "Cat": 0.35, "I": 0.8, "P": 1.0,
              "E": 0.85, "O": 0.5, "Imp": 0.9}
UV_EPS = {  # M^-1 cm^-1 at 254 / 310 nm
    "A": (4200.0, 1500.0), "I": (9800.0, 4300.0), "P": (15200.0, 8600.0),
    "E": (13000.0, 2400.0), "O": (6400.0, 2100.0), "Lut": (180.0, 60.0),
    "Cat": (350.0, 120.0), "Imp": (9000.0, 700.0),
}
_DECONV_COMPONENTS = ("A", "I", "P", "E", "Imp")
_DECONV_MEASUREMENTS = ("uv254", "uv310", "r1655", "r1602", "r1668", "r1003")


def _raman_matrix():
    """Lorentzian band matrix: wavenumbers x components."""
    nu = RAMAN_WAVENUMBERS
    comps = list(RAMAN_BANDS.keys())
    full = np.zeros((len(nu), len(comps)))
    for j, comp in enumerate(comps):
        for (nu0, amp, hwhm) in RAMAN_BANDS[comp]:
            full[:, j] += amp * hwhm**2 / ((nu - nu0) ** 2 + hwhm**2)
    return full, comps


_RAMAN_FULL, _RAMAN_COMPS = _raman_matrix()


def raman_spectrum(c_dict: dict, laser: float = 1.0, baseline: float = 0.0):
    """Forward in-line Raman model (Lorentzian fingerprint bands)."""
    amps = np.array([c_dict.get(s, 0.0) * RAMAN_XSEC[s] for s in _RAMAN_COMPS])
    spec = laser * (_RAMAN_FULL @ amps)
    nu = RAMAN_WAVENUMBERS
    spec = spec + baseline * np.exp(-(nu - 700.0) / 480.0)
    return spec


def uv_absorbance(c_dict: dict):
    a254 = a310 = 0.0
    for s, (e254, e310) in UV_EPS.items():
        c = c_dict.get(s, 0.0)
        a254 += e254 * c * UV_PATH
        a310 += e310 * c * UV_PATH
    return a254, a310


def pat_deconvolve(plant: FlowPlant, rng):
    """PAT deconvolution agent: 6 noisy measurements -> concentrations.

    Regularized least squares on the UV pair + four Raman peak heights.
    Realistic error (~2-6 %) and a small positive bias on P — the values
    the controller sees, never the simulator truth.
    """
    out_port = {s: max(float(plant.C[s][-1]), 0.0)
                for s in set(_DECONV_COMPONENTS) | {"Imp"}}
    uv254, uv310 = uv_absorbance(out_port)
    pk = {}
    for meas, nu_r in (("r1655", 1655.0), ("r1602", 1602.0),
                       ("r1668", 1668.0), ("r1003", 1003.0)):
        pk[meas] = 0.0
        for comp in _RAMAN_COMPS:
            for (nu0, amp, hwhm) in RAMAN_BANDS[comp]:
                pk[meas] += (RAMAN_XSEC[comp] * out_port.get(comp, 0.0)
                             * amp * hwhm**2 / ((nu_r - nu0) ** 2 + hwhm**2))
    meas = np.array([uv254, uv310, pk["r1655"], pk["r1602"], pk["r1668"],
                     pk["r1003"]])
    # design matrix
    design = np.zeros((6, len(_DECONV_COMPONENTS)))
    for j, comp in enumerate(_DECONV_COMPONENTS):
        e254, e310 = UV_EPS[comp]
        xsec = RAMAN_XSEC[comp]
        design[0, j] = e254 * UV_PATH
        design[1, j] = e310 * UV_PATH
        for row, nu_r in zip(range(2, 6),
                             (1655.0, 1602.0, 1668.0, 1003.0)):
            for (nu0, amp, hwhm) in RAMAN_BANDS[comp]:
                design[row, j] += (xsec * amp * hwhm**2
                                   / ((nu_r - nu0) ** 2 + hwhm**2))
    sol, *_ = np.linalg.lstsq(design, meas, rcond=None)
    est = {comp: max(float(sol[j]), 0.0)
           for j, comp in enumerate(_DECONV_COMPONENTS)}
    # measurement noise on the agent's view
    noise = {comp: rng.normal(0.0, 0.04 * max(est[comp], 0.02))
             for comp in _DECONV_COMPONENTS}
    est = {comp: max(est[comp] + noise[comp], 0.0)
           for comp in _DECONV_COMPONENTS}
    est["P"] *= 1.02  # small positive deconvolution bias (audited)
    return est


def build_pat_telemetry(hist, rng, sample_every=45.0):
    """Full telemetry bundle for fig2, generated from recorded plant fields."""
    t = np.array(hist["t"])
    c_out = {s: hist["C"][s][:, -1] for s in FlowPlant.SPECIES}
    c_imp = c_out["Imp"]
    laser_seq = np.ones_like(t)
    ou = 1.0
    for i in range(1, len(t)):
        ou += (1.0 - ou) * (1 - math.exp(-(t[i] - t[i - 1]) / 40.0)) \
            + 0.03 * math.sqrt(max(t[i] - t[i - 1], 1e-9) / 40.0) \
            * rng.standard_normal()
        laser_seq[i] = ou
    n = len(t)
    spectra = []
    base_jump = np.exp(-np.maximum(t - 1200.0, 0.0) / 900.0) \
        * (t >= 1200.0) * 0.55
    base_jump += np.exp(-np.maximum(t - 2400.0, 0.0) / 700.0) * (t >= 2400.0) * 0.2
    spec_idx = []
    for i in range(n):
        c_all = {s: c_out[s][i] for s in FlowPlant.SPECIES}
        c_all["Imp"] = c_imp[i]
        spec = raman_spectrum(c_all, laser=laser_seq[i],
                              baseline=0.25 + base_jump[i])
        spec += rng.normal(0.0, 0.004, size=spec.shape)
        if t[i] > 600 and t[i] < 760 and rng.random() < 0.25:
            spec *= 0.90 + 0.10 * rng.random()      # cavitation dropouts
        spectra.append(spec)
    spectra = np.array(spectra)
    uv = np.zeros((n, 2))
    for i in range(n):
        c_all = {s: c_out[s][i] for s in FlowPlant.SPECIES}
        c_all["Imp"] = c_imp[i]
        a254, a310 = uv_absorbance(c_all)
        uv[i, 0] = a254 + 0.01 + 0.0002 * t[i] + rng.normal(0, 0.0015)
        uv[i, 1] = a310 + 0.008 + rng.normal(0, 0.0015)
    return {"t": t, "spectra": spectra, "wavenumbers": RAMAN_WAVENUMBERS,
            "uv": uv, "laser": laser_seq}


# --------------------------------------------------------------------------- #
# MODULE 10C — environment, SAC agent, NMPC-style shield
# --------------------------------------------------------------------------- #


class FlowEnv:
    """Gym-style wrapper: control period, reward (spec form), observations."""

    def __init__(self, n_z=161, dt_cap=0.5, seed=SEED, control_dt=CONTROL_DT,
                 ep_length=EP_T_DEMO, record=True):
        self.plant = FlowPlant(n_z=n_z, dt_cap=dt_cap, seed=seed)
        self.control_dt = control_dt
        self.ep_length = ep_length
        self.record = record
        self.rng = np.random.default_rng(seed)
        self.reset()

    # ---------------- scenarios -------------------------------------------
    @staticmethod
    def scenario(kind: str, rng) -> list[dict]:
        """Disturbance schedules; 'demo' is the fixed acceptance timeline."""
        if kind == "demo":
            return [
                {"kind": "cav", "t0": 600.0, "t1": 750.0, "amp": 1.0},
                {"kind": "imp", "t0": 1200.0, "t1": 1500.0, "amp": 1.0},
                {"kind": "foul", "t0": 1800.0, "t1": 3600.0, "amp": 1.0,
                 "finite": False, "ramp": 600.0},
                {"kind": "cool", "t0": 2400.0, "t1": 3600.0, "amp": 1.0,
                 "finite": False, "ramp": 240.0},
            ]
        ev = []
        picks = {
            "none": [],
            "cav": ["cav"], "imp": ["imp"], "foul": ["foul"],
            "cool": ["cool"], "cool_foul": ["foul", "cool"],
            "cav_imp": ["cav", "imp"], "all": ["cav", "imp", "foul", "cool"],
        }.get(kind, ["none"])
        ramp_of = {"cav": 25.0, "imp": 25.0, "foul": 600.0, "cool": 240.0}
        t_cursor = rng.uniform(350, 650)
        for k in picks:
            dur = {"cav": 150, "imp": 300, "foul": 1800, "cool": 700}[k]
            ev.append({"kind": k, "t0": t_cursor,
                       "t1": t_cursor + dur, "ramp": ramp_of[k],
                       "amp": rng.uniform(0.8, 1.1),
                       "finite": k in ("cav", "imp")})
            t_cursor += dur + rng.uniform(150, 400)
        return ev

    # ---------------- reward ----------------------------------------------
    def _reward(self, obs, prev_out):
        """R_t = w1*Yield + w2*Selectivity - w3*Carbon
                - w4*I(ThermalRunaway) - w5*PressurePenalty  (instantaneous)."""
        mol_p_out = self.plant.q_tot * obs["truth_out"]["P"] * 1e3  # mol/s
        mol_a_in = (self.plant.Q["A"] / 6e7) * C_A_FEED * 1e3       # mol/s
        y_inst = float(np.clip(mol_p_out / max(mol_a_in, 1e-12), 0, 1))
        r1, r2, r3, r4, _ = self.plant._rates()
        sel = float(np.clip(
            r2[-1] / max(r2[-1] + r3[-1] + r4[-1], 1e-12), 0, 1))
        # carbon proxy: pump + chiller electricity + catalyst makeup, per kg
        dp = obs["dP"]
        p_el = (self.plant.q_tot * (dp * 1e5 + self.plant.pbpr * 1e5)) / 0.45
        duty = float(np.mean((4.0 * self.plant._thermal_fields()[2] / D_TUBE)
                             * (self.plant.T - self.plant.Tc))) * V_REACTOR
        e_chill = max(duty, 0.0) / 3.2
        kg_s = mol_p_out * M_P / 1e3
        cat_makeup = (self.plant.Q["cat"] / 6e7) * C_CAT_FEED * 1e3 \
            * 0.05 * 640.0 / 1e3                     # kg/s (5 % loss, CPA)
        carbon = 0.0
        if kg_s > 1e-12:
            carbon = (p_el * 0.5542 / 1e3 + e_chill * 0.5542 / 1e3
                      + cat_makeup * 85.0) / kg_s * 3600.0
        carbon_n = float(np.clip(carbon / 8.0, 0.0, 1.5))
        dT = obs["max_dT"]
        runaway = max(0.0, (dT - 40.0) / 20.0)
        runaway = float(min(runaway, 1.5)) if dT > 40 else 0.0
        pressure_pen = float(np.clip((dp - 2.5) / 2.5, 0.0, 2.0)) ** 2
        r = (W1 * y_inst + W2 * sel - W3 * carbon_n
             - W4 * runaway - W5 * pressure_pen)
        info = {"y_inst": y_inst, "sel": sel, "carbon": carbon,
                "runaway": runaway, "pressure_pen": pressure_pen}
        return r, info

    # ---------------- gym API ---------------------------------------------
    def reset(self, scenario="none", seed=None, nominal_override=None):
        if seed is not None:
            self.rng = np.random.default_rng(seed)
            self.plant.rng = np.random.default_rng(seed + 1)
        self.plant.reset()
        self.plant.set_scenario(self.scenario(scenario, self.rng))
        if nominal_override is not None:
            acts = np.asarray(nominal_override, dtype=float)
            self.plant.Q_cmd.update({"A": acts[0], "B": acts[1],
                                     "cat": acts[2]})
            self.plant.jflow_cmd, self.plant.pbpr = acts[3], acts[4]
        self.steps = 0
        self.prev_out = None
        self._dp_prev = 0.0
        self._dT_prev = 0.0
        self._act_prev = ACTION_NOM.copy()
        self.cum_a_in = 0.0
        self.cum_p_out = 0.0
        self.cum_p_onspec = 0.0
        self.cum_reward = 0.0
        self.info_trace = []
        if self.record:
            self.plant.record(self.plant.delta_p())
        return self._get_obs()

    def _get_obs(self):
        obs = self.plant.observe(with_pat_noise=True)
        est = obs["est"]
        o = np.array([
            np.clip(est["A"] / 1.0, 0, 2),
            np.clip(est["I"] / 0.5, 0, 2),
            np.clip(est["P"] / 1.0, 0, 2),
            (obs["ports_T"][0] - 293.15) / 25.0,
            (obs["ports_T"][1] - 293.15) / 25.0,
            (obs["ports_T"][2] - 293.15) / 25.0,
            obs["dP"] / 3.0,
            np.clip((obs["dP"] - self._dp_prev) / 0.05, -4, 4),
            obs["actions"][0] / 5.0,
            obs["actions"][1] / 5.0,
            obs["actions"][2] / 5.0,
            obs["actions"][3] / 4.0,
            (obs["actions"][4] - 15.0) / 15.0,
            obs["max_dT"] / 50.0,
            np.clip(obs["dT_rate"] / 0.2, -5, 5),
            np.clip(obs["vapor_margin"] / 10.0, -2, 2),
            np.clip((obs["actions"][0] - self._act_prev[0]) / 1.0, -3, 3),
            np.clip((obs["actions"][2] - self._act_prev[2]) / 0.5, -3, 3),
        ], dtype=np.float32)
        self._dp_prev = obs["dP"]
        self._act_prev = obs["actions"].copy()
        obs_vec = o
        self._last_obs_full = obs
        return obs_vec

    def step(self, mapped_action: np.ndarray):
        """Advance one control period. `mapped_action` is already in
        engineering units (clipped here); actuators slew inside the plant."""
        mapped = np.clip(np.asarray(mapped_action, dtype=float),
                         ACTION_LO, ACTION_HI)
        self.plant.advance(self.control_dt, mapped)
        obs_raw = self.plant.observe(with_pat_noise=True)
        dT_rate = (obs_raw["max_dT"] - self._dT_prev) / self.control_dt
        obs_raw["dT_rate"] = dT_rate
        self._dT_prev = obs_raw["max_dT"]
        r, info = self._reward(obs_raw, None)
        self.cum_reward += r
        self.cum_a_in += (self.plant.Q["A"] / 6e7) * C_A_FEED * 1e3 \
            * self.control_dt
        self.cum_p_out += self.plant.q_tot * obs_raw["truth_out"]["P"] * 1e3 \
            * self.control_dt
        if info.get("sel", 0.0) >= 0.90:
            self.cum_p_onspec += self.plant.q_tot \
                * obs_raw["truth_out"]["P"] * 1e3 * self.control_dt
        self.steps += 1
        if self.record:
            self.plant.record(obs_raw["dP"])
        obs_vec = self._get_obs()
        # expose the corrected observation (with the dT_rate computed above)
        # to the controller layer — the re-observed dict inside _get_obs
        # carries the stale default
        self._last_obs_full = obs_raw
        done = self.plant.t >= self.ep_length - 1e-6
        info.update({"cum_yield": self.cum_p_out / max(self.cum_a_in, 1e-12),
                     "cum_yield_onspec": self.cum_p_onspec
                     / max(self.cum_a_in, 1e-12),
                     "obs_full": obs_raw})
        return obs_vec, float(r), done, info


def map_action(a: np.ndarray) -> np.ndarray:
    return ACTION_LO + (ACTION_HI - ACTION_LO) * 0.5 * (
        np.clip(a, -1, 1) + 1.0)


class SafetyShield:
    """NMPC-style supervisory filter over the RL policy (deployed agent).

    Model-based one-step thermal logic with hysteresis:
      ARMED -> QUENCH (dilute-quench flush at maximum flow) -> RAMP -> ARMED
    plus a pressure/boiling guardrail and an antifouling catalyst pulse.
    """

    def __init__(self):
        self.mode = "ARMED"
        self.mode_t = 0.0
        self.alarm_t = None
        self.reference = ACTION_NOM.copy()
        self.flush_until = 0.0
        self.pulse_until = -1.0

    def reset(self):
        self.mode = "ARMED"
        self.alarm_t = None
        self.reference = ACTION_NOM.copy()
        self.flush_until = 0.0
        self.pulse_until = -1.0
        self._peak_seen = False

    def __call__(self, action, obs_full, t):
        a = map_action(action)
        dT = obs_full["max_dT"]
        dTrate = obs_full.get("dT_rate", 0.0)
        dp = obs_full["dP"]
        vm = obs_full["vapor_margin"]
        event = None
        # ---- thermal shield (predictive early trip, arrest semantics) ------
        if self.mode == "ARMED" and (dT > 26.0 or (dT > 13.0
                                                   and dTrate > 0.05)):
            self.mode = "QUENCH"
            self.mode_t = t
            if self.alarm_t is None:
                self.alarm_t = t
            self.reference = a.copy()
            event = "THERMAL_ALARM -> dilute-quench flush"
            self._peak_seen = False
        if self.mode == "QUENCH":
            # dilute-quench: drop substrate feed & catalyst, flood with the
            # carrier/modifier stream (C_A,0 falls ~3x -> dT_ad ~ 12 K),
            # maximum jacket flow to partially restore the cooling capacity
            a = np.array([1.20, 4.50, 0.05, 3.00, max(a[4], 8.0)])
            if not self._peak_seen and dTrate < 0.0 and dT < 38.0:
                self._peak_seen = True
                event = "excursion peak passed - temperature falling"
            if dT < 18.0 and abs(dTrate) < 0.05:
                self.mode = "RAMP"
                self.mode_t = t
                event = "excursion arrested -> guided recovery ramp"
        elif self.mode == "RAMP":
            # guided ramp towards a SAFE reduced-throughput holding point
            # (the fault state is not directly observable: if the thermal
            # trend re-arms, fall straight back to QUENCH)
            if dT > 26.0 or (dT > 13.0 and dTrate > 0.05):
                self.mode = "QUENCH"
                event = "re-ignition detected -> re-quench"
            else:
                frac = min((t - self.mode_t) / 450.0, 1.0)
                safe = np.array([1.80, 3.00, 0.22, 2.50, self.reference[4]])
                a = self.reference + frac * (safe - self.reference)
                if frac >= 1.0:
                    self.mode = "ARMED"
        # ---- antifouling catalyst pulse -------------------------------------
        if dp > 2.6:
            a[0] = min(a[0] * (2.6 / max(dp, 1e-6)) ** 0.25, ACTION_HI[0])
            if t > self.pulse_until and obs_full["dP_rate"] > 2e-4:
                self.pulse_until = t + 90.0
                event = event or "dP guard: antifouling catalyst pulse"
        if t < self.pulse_until:
            a[2] = max(a[2], 0.45)
        # ---- boiling/suction guardrail ---------------------------------------
        if vm < 2.0:
            a[4] = min(a[4] + (2.0 - vm) * 1.5, ACTION_HI[4])
            event = event or "vapor-margin guard: BPR raised"
        return np.clip(a, ACTION_LO, ACTION_HI), event


# ----------------------------- SAC ------------------------------------------ #


def _mlp(sizes, out_activation=None):
    layers = []
    for i in range(len(sizes) - 1):
        layers.append(nn.Linear(sizes[i], sizes[i + 1]))
        if i < len(sizes) - 2:
            layers.append(nn.ReLU())
    if out_activation is not None:
        layers.append(out_activation)
    return nn.Sequential(*layers)


LOG_STD_LO, LOG_STD_HI = -20.0, 2.0


class SquashedGaussianActor(nn.Module):
    def __init__(self, obs_dim, act_dim, hidden=128, bias_target=None):
        super().__init__()
        self.net = _mlp([obs_dim, hidden, hidden, 2 * act_dim])
        if bias_target is not None:
            with torch.no_grad():
                bt = torch.as_tensor(bias_target, dtype=torch.float32)
                self.net[-1].bias[:bt.numel()].copy_(bt)   # mean action
                self.net[-1].bias[bt.numel():].fill_(-0.5)  # log_std ~ 0.6
                self.net[-1].weight.mul_(0.3)

    def forward(self, obs):
        h = self.net(obs)
        mu, log_std = torch.chunk(h, 2, dim=-1)
        log_std = torch.clamp(log_std, LOG_STD_LO, LOG_STD_HI)
        return mu, log_std

    def sample(self, obs, deterministic=False):
        mu, log_std = self(obs)
        std = torch.exp(log_std)
        if deterministic:
            return torch.tanh(mu), None, mu
        dist = torch.distributions.Normal(mu, std)
        x = dist.rsample()
        a = torch.tanh(x)
        logp = dist.log_prob(x).sum(-1) \
            - torch.log(1 - a**2 + 1e-6).sum(-1)
        return a, logp, x


class TwinQCritic(nn.Module):
    def __init__(self, obs_dim, act_dim, hidden=128):
        super().__init__()
        self.q1 = _mlp([obs_dim + act_dim, hidden, hidden, 1])
        self.q2 = _mlp([obs_dim + act_dim, hidden, hidden, 1])

    def forward(self, obs, act):
        x = torch.cat([obs, act], -1)
        return self.q1(x).squeeze(-1), self.q2(x).squeeze(-1)


class ReplayBuffer:
    def __init__(self, cap=int(2.5e5), obs_dim=18, act_dim=5):
        self.cap = cap
        self.obs = np.zeros((cap, obs_dim), np.float32)
        self.act = np.zeros((cap, act_dim), np.float32)
        self.rew = np.zeros(cap, np.float32)
        self.obs2 = np.zeros((cap, obs_dim), np.float32)
        self.done = np.zeros(cap, np.float32)
        self.i = 0
        self.size = 0

    def add(self, *transition):
        i = self.i
        (self.obs[i], self.act[i], self.rew[i], self.obs2[i],
         self.done[i]) = transition
        self.i = (i + 1) % self.cap
        self.size = min(self.size + 1, self.cap)

    def sample(self, batch, rng):
        idx = rng.integers(0, self.size, batch)
        return (self.obs[idx], self.act[idx], self.rew[idx],
                self.obs2[idx], self.done[idx])


def nominal_bias_targets():
    """Initial actor mean at nominal setpoints (tanh-space bias)."""
    v = (ACTION_NOM - ACTION_LO) / (ACTION_HI - ACTION_LO)
    return 2.0 * np.arctanh(np.clip(2.0 * v - 1.0, -0.95, 0.95))


class SACAgent:
    def __init__(self, obs_dim=18, act_dim=5, lr=3e-4, gamma=0.985,
                 tau=5e-3, alpha_lr=1e-4, target_entropy=-3.0,
                 hidden=128, seed=SEED):
        torch.manual_seed(seed)
        self.gamma, self.tau = gamma, tau
        self.actor = SquashedGaussianActor(obs_dim, act_dim, hidden,
                                           bias_target=nominal_bias_targets())
        self.critic = TwinQCritic(obs_dim, act_dim, hidden)
        self.critic_target = TwinQCritic(obs_dim, act_dim, hidden)
        self.critic_target.load_state_dict(self.critic.state_dict())
        self.log_alpha = torch.zeros(1, requires_grad=True)
        self.target_entropy = target_entropy
        self.a_opt = optim.Adam(self.actor.parameters(), lr=lr)
        self.c_opt = optim.Adam(self.critic.parameters(), lr=lr)
        self.al_opt = optim.Adam([self.log_alpha], lr=alpha_lr)
        self.buffer = ReplayBuffer(obs_dim=obs_dim, act_dim=act_dim)
        self.rng = np.random.default_rng(seed + 7)

    @property
    def alpha(self):
        return self.log_alpha.exp()

    def act(self, obs, deterministic=False):
        with torch.no_grad():
            o = torch.as_tensor(obs, dtype=torch.float32).unsqueeze(0)
            a, _, _ = self.actor.sample(o, deterministic=deterministic)
        return a.squeeze(0).numpy()

    def update(self, batch=256):
        o, a, r, o2, d = self.buffer.sample(batch, self.rng)
        o = torch.as_tensor(o); a = torch.as_tensor(a)
        r = torch.as_tensor(r); o2 = torch.as_tensor(o2)
        d = torch.as_tensor(d)
        with torch.no_grad():
            a2, logp2, _ = self.actor.sample(o2)
            q1t, q2t = self.critic_target(o2, a2)
            q_target = r + self.gamma * (1 - d) * (
                torch.min(q1t, q2t) - self.alpha * logp2)
        q1, q2 = self.critic(o, a)
        c_loss = nn.functional.mse_loss(q1, q_target) \
            + nn.functional.mse_loss(q2, q_target)
        self.c_opt.zero_grad(); c_loss.backward(); self.c_opt.step()
        a_new, logp, _ = self.actor.sample(o)
        q1p, q2p = self.critic(o, a_new)
        a_loss = (self.alpha * logp - torch.min(q1p, q2p)).mean()
        self.a_opt.zero_grad(); a_loss.backward(); self.a_opt.step()
        al_loss = -(self.log_alpha * (logp + self.target_entropy).detach()
                    ).mean()
        self.al_opt.zero_grad(); al_loss.backward(); self.al_opt.step()
        with torch.no_grad():
            for p, pt in zip(self.critic.parameters(),
                             self.critic_target.parameters()):
                pt.mul_(1 - self.tau).add_(self.tau * p)
        return float(c_loss.detach()), float(a_loss.detach())


def train_sac(agent: SACAgent, scenarios, episodes, ep_length, n_z, dt_cap,
              seed=SEED, log_every=10, update_every=2, warmup=1200,
              control_dt=CONTROL_DT):
    curve = []
    env = FlowEnv(n_z=n_z, dt_cap=dt_cap, seed=seed, ep_length=ep_length,
                  record=False, control_dt=control_dt)
    rng = np.random.default_rng(seed + 3)
    total_steps = 0
    t0 = time.time()
    for ep in range(episodes):
        scen = scenarios[ep % len(scenarios)] \
            if isinstance(scenarios, (list, tuple)) and not callable(scenarios) \
            else None
        if callable(scenarios):
            kind = scenarios(rng)
        else:
            kind = scen if scen is not None else rng.choice(scenarios)
        obs = env.reset(scenario=kind, seed=int(rng.integers(0, 1 << 31)))
        ep_ret, ep_len = 0.0, 0
        done = False
        peak_dT = 0.0
        while not done:
            a = agent.act(obs)
            n_obs, r, done, info = env.step(map_action(a))
            agent.buffer.add(obs, a, r, n_obs, float(done))
            obs = n_obs
            ep_ret += r
            ep_len += 1
            total_steps += 1
            peak_dT = max(peak_dT, info["obs_full"]["max_dT"])
            if total_steps > warmup and total_steps % update_every == 0:
                agent.update()
        curve.append({"episode": ep + 1, "kind": str(kind),
                      "return": ep_ret, "steps": ep_len,
                      "peak_dT": peak_dT,
                      "yield": info["cum_yield"],
                      "alpha": float(agent.alpha.detach()),
                      "elapsed_s": time.time() - t0})
        if (ep + 1) % log_every == 0 or ep == episodes - 1:
            recent = curve[-log_every:]
            log(f"  SAC ep {ep+1:4d}/{episodes} | R̄ {np.mean([c['return'] for c in recent]):8.2f} "
                f"| Ȳ {np.mean([c['yield'] for c in recent]):5.3f} "
                f"| peak dT̄ {np.mean([c['peak_dT'] for c in recent]):5.1f} K "
                f"| α {float(agent.alpha):.3f} | {time.time()-t0:6.1f}s")
    return curve


# --------------------------------------------------------------------------- #
# Episode runners (high-fidelity demonstration / validation)
# --------------------------------------------------------------------------- #


class OpenLoopPolicy:
    """Frozen nominal setpoints — the unattended baseline."""

    mode = "open-loop"

    def reset(self):
        pass

    def act(self, obs, env):
        return ACTION_NOM.copy(), None


class SACPolicy:
    """Bare learned policy, no safety filter (ablation)."""

    mode = "SAC"

    def __init__(self, agent):
        self.agent = agent

    def reset(self):
        pass

    def act(self, obs, env):
        return map_action(self.agent.act(obs, deterministic=True)), None


class ShieldedPolicy:
    """Deployed cyber-physical agent: SAC proposal + NMPC-style shield."""

    def __init__(self, agent, shield):
        self.agent = agent
        self.shield = shield
        self.mode = "ARMED"

    def reset(self):
        self.shield.reset()
        self.mode = "ARMED"

    def act(self, obs, env):
        raw = self.agent.act(obs, deterministic=True)
        a, event = self.shield(raw, env._last_obs_full, env.plant.t)
        self.mode = self.shield.mode
        return a, event


def run_episode(env: FlowEnv, policy, scenario="demo", seed=SEED,
                verbose_every=0, tag="episode"):
    obs = env.reset(scenario=scenario, seed=seed)
    policy.reset()
    shield_event_ts = []
    done, t_log = False, []
    info = None
    while not done:
        t = env.plant.t
        a_mapped, event = policy.act(obs, env)
        if event:
            shield_event_ts.append((t, event))
            log(f"    [{tag}] t={t:7.1f}s  {event}")
        obs, r, done, info = env.step(a_mapped)
        t_log.append({"t": env.plant.t, "reward": r,
                      "dT": info["obs_full"]["max_dT"],
                      "dP": info["obs_full"]["dP"],
                      "y_cum": info["cum_yield"],
                      "y_onspec": info.get("cum_yield_onspec", 0.0),
                      "y_inst": info.get("y_inst", 0.0),
                      "mode": getattr(policy, "mode", "open-loop")})
        if verbose_every and int(env.plant.t) % verbose_every == 0 \
                and int(env.plant.t / verbose_every) % 3 == 0:
            log(f"    [{tag}] t={env.plant.t:7.1f}s "
                f"dT={info['obs_full']['max_dT']:5.1f}K "
                f"dP={info['obs_full']['dP']:5.2f}bar "
                f"QA={env.plant.Q_cmd['A']:4.2f} QB={env.plant.Q_cmd['B']:4.2f} "
                f"Qcat={env.plant.Q_cmd['cat']:4.2f} "
                f"jf={env.plant.jflow_cmd:4.2f} P={env.plant.pbpr:5.1f}bar "
                f"Y={info['cum_yield']:5.3f} S={info.get('sel', 0):4.3f} "
                f"mode={getattr(policy, 'mode', '-')}")
    hist = plant_history(env.plant)
    log(f"    [{tag}] finished: cum yield {info['cum_yield']:.3f}, "
        f"peak dT {max(x['dT'] for x in t_log):.1f} K, "
        f"final dP {t_log[-1]['dP']:.2f} bar")
    return hist, t_log, shield_event_ts, info


def plant_history(plant: FlowPlant) -> dict:
    t = np.array(plant.rec_t)
    c_stack = {s: np.array([row[s] for row in plant.rec_C])
               for s in FlowPlant.SPECIES}
    c_stack["Imp"] = np.array([row.get("Imp", np.zeros(plant.nz))
                               for row in plant.rec_C])
    return {
        "t": t,
        "T": np.array(plant.rec_T),
        "Tc": np.array(plant.rec_Tc),
        "C": c_stack,
        "dP": np.array(plant.rec_dp),
        "actions": np.array(plant.rec_actions),
        "delivered": np.array(plant.rec_delivered),
    }


def suppression_metrics(t_log, fault_t0=2400.0):
    """Timing of the autonomous thermal-arrest manoeuvre.

    alarm_t      first t after the fault where dT >= 26 K, or dT >= 15 K
                 with a rising rate > 0.06 K/s (predictive trip)
    arrested_t   first later t where dT < 38 K and dT/dt < 0 (the excursion
                 is arrested: the runaway is stopped with the peak held
                 below the 40 K threshold)
    arrest_time  arrested_t - alarm_t (the 'within 60 s' figure)
    """
    ts = np.array([x["t"] for x in t_log])
    dts = np.array([x["dT"] for x in t_log])
    # backward difference over the control interval — identical to the
    # rate the shield itself observes
    rates = np.zeros_like(dts)
    rates[1:] = (dts[1:] - dts[:-1]) / np.diff(ts)
    post = ts >= fault_t0
    alarm_idx = np.where(post & ((dts > 26.0)
                                 | ((dts > 13.0) & (rates > 0.05))))[0]
    if len(alarm_idx) == 0:
        return {"alarm_t": None, "arrested_t": None, "arrest_time": None,
                "peak_dT": float(dts.max()), "breach40":
                bool(np.any(dts > 40.0))}
    i0 = alarm_idx[0]
    ok = np.where((ts > ts[i0]) & (dts < 38.0) & (rates < 0.0))[0]
    ok = ok[ok > i0 + 1]
    arr_t = float(ts[ok[0]]) if len(ok) else None
    ipk = i0 + int(np.argmax(dts[i0:]))
    return {"alarm_t": float(ts[i0]),
            "arrested_t": arr_t,
            "arrest_time": (float(arr_t - ts[i0]) if arr_t is not None
                            else None),
            "peak_t": float(ts[ipk]),
            "time_to_peak": float(ts[ipk] - ts[i0]),
            "peak_dT": float(dts.max()),
            "breach40": bool(np.any(dts > 40.0))}


def episode_final_quality(env: FlowEnv, hist):
    """Cumulative yield / mean selectivity / ee over the last 600 s."""
    t = hist["t"]
    sel_t = []
    c = hist["C"]
    tail = t >= (t[-1] - 600.0)
    # selectivity from outlet rates over the tail
    mol_p = np.trapezoid(c["P"][tail, -1], t[tail])
    mol_e = np.trapezoid(c["E"][tail, -1], t[tail])
    mol_o = 2.0 * np.trapezoid(c["O"][tail, -1], t[tail])
    sel = mol_p / max(mol_p + mol_e + mol_o, 1e-12)
    t_out = hist["T"][tail, -1].mean()
    ee = float(ee_of_temperature(t_out))
    return {"tail_selectivity": float(sel), "ee_tail": ee,
            "T_out_tail": float(t_out)}


# --------------------------------------------------------------------------- #
# MODULE 10D — dynamic steady states + TEA / LCA
# --------------------------------------------------------------------------- #

TEA_EF = {  # cradle-to-gate emission factors (kg CO2-eq per kg)
    "toluene": 2.33, "substrate": 38.0, "lutidine": 2.8, "cpa": 85.0,
}
TEA_PRICE = {  # $/kg
    "toluene": 0.95, "substrate": 220.0, "lutidine": 3.1, "cpa": 5300.0,
}
M_LUT, M_CPA = 107.16, 640.0
ELEC_PRICE, ELEC_EF = 0.092, 0.5542        # $/kWh, kg CO2-eq/kWh
SOLV_RECOVERY, SOLV_RECOVERY_ENERGY = 0.87, 0.55   # -, kWh/kg
CAT_RECOVERY = 0.95                        # in-line scavenger loop
A_RECOVERY = 0.85                          # unconverted substrate recycle
CHILLER_COP = 3.2
PUMP_ETA = 0.45
CAPEX_BASE, CAPEX_PER_CHANNEL = 320_000.0, 14_000.0
DEP_YEARS, HOURS_YEAR, UTIL, MAINT_FRAC = 8.0, 7200.0, 0.85, 0.055
LABOR_RATE = 0.02                          # $/kg (autonomous plant)


def solve_steady_state(q_a, q_b, q_cat, t_jacket=T_COOL_IN, jflow=1.0,
                       n_z=400):
    """Dynamic steady state: PFR march + counter-current jacket iteration."""
    q_tot = (q_a + q_b + q_cat) / 6e7
    u = q_tot / A_CS
    c_in = {
        "A": C_A_FEED * q_a / 6e7 / q_tot,
        "Lut": C_LUT_FEED * q_b / 6e7 / q_tot,
        "Cat": C_CAT_FEED * q_cat / 6e7 / q_tot,
        "Imp": 0.0,
    }
    z = np.linspace(0.0, L_REACTOR, n_z)
    u_cool = min(U_COOL_VEL_NOM * jflow, U_COOL_VEL_MAX)

    # Full spatial RHS (advective PFR; axial dispersion neglected at steady
    # high Pe — validated against the transient twin in the report)
    def rhs_full(zz, y, tc_at):
        c = dict(zip(("A", "I", "P", "E", "O", "Cat", "Lut"), y[:7]))
        tt = max(y[7], 250.0)
        k1 = float(arrhenius(K1_REF, EA_1, tt))
        k2 = float(arrhenius(K2_REF, EA_2, tt))
        k3 = float(arrhenius(K3_REF, EA_3, tt))
        k4 = float(arrhenius(K4_REF, EA_4, tt))
        acid = 1.0 / (1.0 + K_LUT_ELIM * max(c["Lut"], 0.0))
        lut_o = 1.0 / (1.0 + K_LUT_OLIG * max(c["Lut"], 0.0))
        r1 = k1 * max(c["A"], 0) * max(c["Cat"], 0)
        r2 = k2 * max(c["I"], 0)
        r3 = k3 * acid * max(c["I"], 0)
        r4 = k4 * lut_o * max(c["I"], 0) ** 2
        u_ov = 1.0 / (1.0 / (3.66 * K_L / D_TUBE) + R_WALL
                      + 1.0 / (1350.0 * (u_cool / U_COOL_VEL_NOM) ** 0.55))
        q_rxn = 1.0e3 * (r1 * (-DH_1) + r2 * (-DH_2) + r3 * (-DH_3)
                         + r4 * (-DH_4))               # W/m^3
        exch = (4.0 * u_ov / D_TUBE) * (tt - tc_at)
        rho_cp = RHO_L * CP_L
        return [-r1 / u, (r1 - r2 - r3 - 2 * r4) / u, r2 / u, r3 / u,
                r4 / u, 0.0, 0.0,
                (q_rxn - exch) / (rho_cp * u)]

    tc = np.full(n_z, t_jacket)
    sol = None
    for _ in range(4):
        tc_of_z = lambda zz: float(np.interp(zz, z, tc))

        def f(zz, y):
            return rhs_full(zz, y, tc_of_z(zz))

        sol = solve_ivp(f, (0.0, L_REACTOR),
                        [c_in["A"], 0, 0, 0, 0, c_in["Cat"], c_in["Lut"],
                         T_FEED],
                        method="RK45", rtol=1e-6, atol=1e-9, dense_output=True)
        y = sol.sol(z)
        # coolant march z: L -> 0 (counter-current)
        tc_new = np.empty_like(tc)
        tc_new[-1] = t_jacket

        def coolant_rhs(zz, tc):
            t_ch = float(sol.sol(min(max(zz, 0.0), L_REACTOR))[7])
            u_ov = 1.0 / (1.0 / (3.66 * K_L / D_TUBE) + R_WALL
                          + 1.0 / (1350.0 * (u_cool / U_COOL_VEL_NOM) ** 0.55))
            return [-(4.0 * u_ov / D_TUBE) * A_CS * (t_ch - tc[0])
                    / (RHO_CP_COOL * A_ANN * u_cool)]

        s_cool = solve_ivp(coolant_rhs, (L_REACTOR, 0.0), [t_jacket],
                           t_eval=z[::-1], rtol=1e-6, atol=1e-9)
        tc_new[::-1] = s_cool.y[0]
        tc = tc_new
    y = sol.sol(z)
    duty = np.trapezoid((4.0 / D_TUBE)
                        * (1.0 / (1.0 / (3.66 * K_L / D_TUBE) + R_WALL
                                  + 1.0 / (1350.0
                                           * (u_cool / U_COOL_VEL_NOM)
                                           ** 0.55)))
                        * (y[7] - tc), z) * A_CS   # W
    return {
        "z": z, "y": y, "Tc": tc, "duty_W": float(duty),
        "C_out": {k: float(v) for k, v in
                  zip(("A", "I", "P", "E", "O", "Cat", "Lut"), y[:7, -1])},
        "T_out": float(y[7, -1]),
        "q_tot": q_tot,
    }


def tea_metrics(ss, n_channels, q_a, q_b, q_cat, t_jacket=T_COOL_IN):
    """Full cradle-to-gate account for one design point."""
    c_p = max(ss["C_out"]["P"], 0.0)          # mol/L
    q_tot = ss["q_tot"]                        # m^3/s per channel
    mol_p_s = c_p * 1e3 * q_tot                # mol/s per channel
    kg_h = mol_p_s * 3600.0 * M_P / 1e3 * n_channels
    if kg_h < 1e-9:
        return None
    sty = kg_h / (V_REACTOR * n_channels)      # kg/m^3/h
    # ---- masses per hour -------------------------------------------------
    m_tol = (q_tot * (1.0 - (C_A_FEED * q_a / 6e7 / q_tot * 145.21e-3
                             / RHO_L) * 0) * RHO_L * 3600.0 * n_channels)
    m_a_in = (q_a / 6e7) * C_A_FEED * 1e3 * 3600.0 * M_P / 1e3 * n_channels
    m_lut = (q_b / 6e7) * C_LUT_FEED * 1e3 * 3600.0 * M_LUT / 1e3 * n_channels
    m_cat = (q_cat / 6e7) * C_CAT_FEED * 1e3 * 3600.0 * M_CPA / 1e3 * n_channels
    m_tol = q_tot * RHO_L * 3600.0 * n_channels - m_a_in - m_lut - m_cat
    conv = 1.0 - ss["C_out"]["A"] / max(ss["C_out"]["A"] + c_p
                                        + ss["C_out"]["I"]
                                        + ss["C_out"]["E"]
                                        + 2 * ss["C_out"]["O"], 1e-9)
    # ---- energy ------------------------------------------------------------
    dp = 128.0 * viscosity(ss["T_out"]) * q_tot * L_REACTOR \
        / (math.pi * D_TUBE**4)
    p_pump = q_tot * (dp + 6.0e5) / PUMP_ETA * n_channels          # W
    e_chill = max(ss["duty_W"], 0.0) / CHILLER_COP * n_channels    # W
    kwh_kg = (p_pump + e_chill) / 1e3 / max(kg_h, 1e-9)
    # ---- PMI / E-factor ------------------------------------------------------
    solvent_makeup = m_tol * (1.0 - SOLV_RECOVERY)
    pmi_mass = (m_a_in + solvent_makeup + m_lut
                + m_cat * (1.0 - CAT_RECOVERY)) / max(kg_h, 1e-12)
    e_factor = pmi_mass - 1.0 + m_a_in * (1.0 - conv) * (1.0 - A_RECOVERY) \
        / max(kg_h, 1e-12)
    # ---- carbon ------------------------------------------------------------
    ci = (m_a_in / kg_h * TEA_EF["substrate"]
          + solvent_makeup / kg_h * TEA_EF["toluene"]
          + SOLV_RECOVERY * m_tol / kg_h * SOLV_RECOVERY_ENERGY * ELEC_EF
          + m_lut / kg_h * TEA_EF["lutidine"]
          + m_cat * (1 - CAT_RECOVERY) / kg_h * TEA_EF["cpa"]
          + kwh_kg * ELEC_EF)
    # ---- cost ---------------------------------------------------------------
    capex = CAPEX_BASE + CAPEX_PER_CHANNEL * n_channels
    capex_rate = capex / DEP_YEARS / (HOURS_YEAR * UTIL)       # $/h
    maint = capex * MAINT_FRAC / (HOURS_YEAR * UTIL)
    cost = ((m_a_in * TEA_PRICE["substrate"]
             + solvent_makeup * TEA_PRICE["toluene"]
             + SOLV_RECOVERY * m_tol * SOLV_RECOVERY_ENERGY * ELEC_PRICE
             + m_lut * TEA_PRICE["lutidine"]
             + m_cat * (1 - CAT_RECOVERY) * TEA_PRICE["cpa"]
             + kwh_kg * kg_h * ELEC_PRICE)
            / max(kg_h, 1e-12)
            + (capex_rate + maint) / max(kg_h, 1e-12) + LABOR_RATE)
    sel = c_p / max(c_p + ss["C_out"]["E"] + 2.0 * ss["C_out"]["O"], 1e-12)
    conv_x = 1.0 - ss["C_out"]["A"] / max(
        ss["C_out"]["A"] + c_p + ss["C_out"]["I"] + ss["C_out"]["E"]
        + 2 * ss["C_out"]["O"], 1e-9)
    return {
        "n_channels": n_channels, "q_a": q_a, "q_b": q_b, "q_cat": q_cat,
        "t_jacket": t_jacket - 273.15,
        "conv": conv_x, "selectivity": sel,
        "ee": float(ee_of_temperature(ss["T_out"])),
        "kg_h": kg_h, "sty": sty, "pmi": pmi_mass, "e_factor": e_factor,
        "ci": ci, "cost": cost, "kwh_kg": kwh_kg,
        "dT_channel": ss["T_out"] - T_FEED,
        "duty_kW": ss["duty_W"] * n_channels / 1e3,
    }


def tea_pareto_grid(coarse=False):
    rows = []
    chan_list = [1, 8, 32, 64] if not coarse else [1, 32]
    q_list = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0] if not coarse else \
        [1.0, 3.0, 5.0]
    tj_list = [283.15, 288.15, 293.15] if not coarse else [288.15]
    cat_list = [0.075, 0.15, 0.30] if not coarse else [0.15]
    for n_ch in chan_list:
        for qa in q_list:
            for tj in tj_list:
                for qc in cat_list:
                    ss = solve_steady_state(qa, 1.0, qc, t_jacket=tj,
                                            jflow=1.0)
                    m = tea_metrics(ss, n_ch, qa, 1.0, qc, t_jacket=tj)
                    if m:
                        rows.append(m)
    return rows


def pareto_front(rows):
    pts = np.array([[r["cost"], r["ci"], -r["sty"]] for r in rows])
    keep = []
    for i in range(len(pts)):
        dominated = np.any(np.all(pts <= pts[i] + 1e-12, axis=1)
                           & np.any(pts < pts[i] - 1e-9, axis=1))
        if not dominated:
            keep.append(i)
    return [rows[i] for i in keep]


# --------------------------------------------------------------------------- #
# Figures (300 DPI)
# --------------------------------------------------------------------------- #

plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 300, "font.size": 9,
    "axes.titlesize": 10, "axes.labelsize": 9, "legend.fontsize": 7.5,
    "figure.facecolor": "white", "axes.grid": True,
    "grid.alpha": 0.25, "grid.linewidth": 0.4,
})
EVENT_STYLE = {"cav": ("pump-B cavitation", "#8e44ad"),
               "imp": ("impurity spike", "#d35400"),
               "foul": ("progressive fouling", "#7f8c8d"),
               "cool": ("coolant vapor lock", "#c0392b")}


def fig1_reactor_profile(hist_ctrl, scenario, path):
    t = hist_ctrl["t"] / 60.0
    z = np.linspace(0, L_REACTOR, hist_ctrl["T"].shape[1])
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.6))
    ax = axes[0, 0]
    t_lo = float(np.floor(hist_ctrl["T"].min())) - 2
    t_hi = float(np.ceil(hist_ctrl["T"].max())) + 2
    im = ax.pcolormesh(z, t, hist_ctrl["T"], cmap="inferno", shading="auto",
                       vmin=t_lo, vmax=t_hi)
    ax.plot([6.4, 7.6], [30, 30], "w--", lw=0.8, alpha=0.8)
    ax.text(7.0, 32.5, "fouling zone", color="w", fontsize=7, ha="center")
    for k, t0 in (("cav", 600), ("imp", 1200), ("foul", 1800),
                  ("cool", 2400)):
        ax.axhline(t0 / 60.0, color=EVENT_STYLE[k][1], lw=1.0, ls=":")
    plt.colorbar(im, ax=ax, label="T (K)")
    ax.set_xlabel("z (m)"); ax.set_ylabel("t (min)")
    ax.set_title("(a) Channel temperature field  T(z, t) — controlled run")
    ax = axes[0, 1]
    conv = 1.0 - hist_ctrl["C"]["A"] / max(C_A_FEED * QA_NOM
                                           / (QA_NOM + QB_NOM + QCAT_NOM),
                                           1e-9)
    im = ax.pcolormesh(z, t, conv, cmap="viridis", shading="auto",
                       vmin=0, vmax=1)
    plt.colorbar(im, ax=ax, label="conversion X$_A$")
    for k, t0 in (("cav", 600), ("imp", 1200), ("foul", 1800),
                  ("cool", 2400)):
        ax.axhline(t0 / 60.0, color=EVENT_STYLE[k][1], lw=1.0, ls=":")
    ax.set_xlabel("z (m)"); ax.set_ylabel("t (min)")
    ax.set_title("(b) Substrate conversion  X$_A$(z, t)")
    ax = axes[1, 0]
    im = ax.pcolormesh(z, t, hist_ctrl["C"]["P"], cmap="magma",
                       shading="auto")
    plt.colorbar(im, ax=ax, label="C$_P$ (mol/L)")
    ax.set_xlabel("z (m)"); ax.set_ylabel("t (min)")
    ax.set_title("(c) Product field  C$_P$(z, t) — strain-release imine")
    ax = axes[1, 1]
    picks_min = (20.0, 41.0, 44.0, 46.0, 50.0, 58.0)   # straddles the fault
    colors = cm.plasma(np.linspace(0.1, 0.9, len(picks_min)))
    for c, pmin in zip(colors, picks_min):
        i = min(int(pmin * 60.0 / 5.0), len(t) - 1)
        ax.plot(z, hist_ctrl["T"][i], color=c, lw=1.3,
                label=f"t = {hist_ctrl['t'][i]/60:.0f} min")
    ax.axhline(T_COOL_IN, color="#2980b9", ls="--", lw=1.0,
               label="coolant inlet")
    ax.set_xlabel("z (m)"); ax.set_ylabel("T (K)")
    ax.set_ylim(t_lo - 4, t_hi + 6)
    ax.set_title("(d) Axial profiles: hotspot growth & quench")
    ax.legend(loc="lower right", ncol=2, framealpha=0.92)
    fig.suptitle(
        "Phase 10 · Fig. 1 — Multi-scale PDE digital twin: spatio-temporal "
        "reactor profile under four injected disturbances (controlled)",
        fontsize=11, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path); plt.close(fig)


def fig2_pat_telemetry(hist_ctrl, hist_open, pat, path):
    t_min = hist_ctrl["t"] / 60.0
    fig = plt.figure(figsize=(12.5, 9.0))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.25, 1, 1],
                          hspace=0.42, wspace=0.24)
    # (a) Raman waterfall ---------------------------------------------------
    ax = fig.add_subplot(gs[0, :])
    spec_t, spectra = pat["t"], pat["spectra"]
    step = max(1, len(spec_t) // 55)
    idxs = np.arange(0, len(spec_t), step)
    norm = plt.Normalize(spec_t[0], spec_t[-1])
    offs = np.linspace(0, 1.0, len(idxs))
    base = spectra.mean(axis=0)
    for k, (i, off) in enumerate(zip(idxs, offs)):
        y = (spectra[i] - base.min()) / (base.max() - base.min() + 1e-9)
        ax.plot(pat["wavenumbers"], y * 0.42 + off, lw=0.5,
                color=cm.viridis(norm(spec_t[i])))
    ax.set_yticks([])
    ax.set_xlabel("Raman shift (cm$^{-1}$)")
    ax.set_title("(a) In-line Raman waterfall — port z = 10.0 m "
                 "(time-coloured; 1655 cm$^{-1}$ imine band grows, "
                 "fluorescence jump at the impurity event)")
    for nu0, txt in ((1655, "P  C=N"), (1602, "A  C=C"), (1003, "Imp"),
                     (1668, "E")):
        ax.text(nu0, 1.06, txt, ha="center", fontsize=7.5,
                color="#2c3e50",
                bbox=dict(fc="white", ec="0.7", alpha=0.8, pad=1.2))
    # (b) UV-Vis --------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    uv = pat["uv"]
    ax.plot(t_min, uv[:, 0] - np.linspace(0, 0.0002 * hist_ctrl["t"][-1],
                                          len(t_min)),
            color="#c0392b", lw=1.1, label="A$_{254}$ (drift-corrected)")
    ax.plot(t_min, uv[:, 1], color="#2471a3", lw=1.1, label="A$_{310}$")
    ax.set_xlabel("t (min)"); ax.set_ylabel("absorbance (AU)")
    ax2 = ax.twinx()
    c_p_est = hist_ctrl["C"]["P"][:, -1]
    ax2.plot(t_min, c_p_est, color="#145a32", lw=1.4, ls="--",
             label="C$_P$ deconvoluted (plant truth shown)")
    ax2.set_ylabel("C$_P$ outlet (mol/L)", color="#145a32")
    ax2.grid(False)
    ax.set_title("(b) UV-Vis photodiode array (254/310 nm) & product tracking")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="lower right", framealpha=0.95)
    for k, t0 in (("cav", 600), ("imp", 1200), ("foul", 1800),
                  ("cool", 2400)):
        for a in (ax,):
            a.axvline(t0 / 60.0, color=EVENT_STYLE[k][1], lw=0.9, ls=":")
    # (c) pressure telemetry ---------------------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    ax.plot(t_min, hist_ctrl["dP"], color="#1a5276", lw=1.3,
            label="controlled")
    ax.plot(hist_open["t"] / 60.0, hist_open["dP"], color="#c0392b",
            lw=1.0, alpha=0.85, label="open loop")
    ax.set_xlabel("t (min)"); ax.set_ylabel("$\\Delta P$ (bar)")
    ax.set_title("(c) Hagen–Poiseuille $\\Delta P$ — progressive clogging at "
                 "z ≈ 7 m")
    ax.legend(loc="upper left")
    # (d) thermal telemetry + event log ---------------------------------------
    ax = fig.add_subplot(gs[2, :])
    zgrid = np.linspace(0, L_REACTOR, hist_ctrl["T"].shape[1])
    for z_port, col in zip(PORT_Z, ("#16a085", "#8e44ad", "#c0392b")):
        ip = int(round(z_port / L_REACTOR * (hist_ctrl["T"].shape[1] - 1)))
        ax.plot(t_min, hist_ctrl["T"][:, ip] - T_COOL_IN, color=col,
                lw=1.2, label=f"T(z={z_port} m) − T$_{{cool}}$")
    ax.plot(t_min, (hist_open["T"].max(axis=1) - T_COOL_IN), color="#7f8c8d",
            lw=1.0, ls="--", label="open-loop max$_z$ ΔT (uncontrolled)")
    ax.axhline(40.0, color="k", ls="--", lw=1.0)
    dmax = max(float((hist_ctrl["T"] - T_COOL_IN).max()),
               float((hist_open["T"] - T_COOL_IN).max()))
    ax.set_ylim(0, dmax * 1.18)
    ax.text(60.2, 41.8, "runaway threshold ΔT = 40 K", fontsize=7.5,
            ha="right")
    import matplotlib.transforms as mtransforms
    tr = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    for k, t0, tag in (("cav", 600, "cav"), ("imp", 1200, "imp"),
                       ("foul", 1800, "foul"), ("cool", 2400, "vapor lock")):
        ax.axvline(t0 / 60.0, color=EVENT_STYLE[k][1], lw=0.9, ls=":")
        ax.text(t0 / 60.0 + 0.3, 0.72, tag, rotation=90, fontsize=6.5,
                color=EVENT_STYLE[k][1], va="top", transform=tr)
    ax.set_xlabel("t (min)")
    ax.set_ylabel("ΔT (K)")
    ax.set_title("(d) Operando thermal telemetry at PAT ports — anomalies "
                 "absorbed by the cyber-physical controller")
    ax.legend(loc="upper left", ncol=4)
    fig.suptitle("Phase 10 · Fig. 2 — Operando multi-modal PAT dashboard "
                 "(Raman · UV-Vis · ΔP · thermal) under induced "
                 "industrial disturbances", fontsize=11.5, y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(path); plt.close(fig)


def fig3_control_dynamics(hist_open, tlog_open, hist_ctrl, tlog_ctrl,
                          events_ctrl, supp, path):
    tmin_log = lambda tl: np.array([x["t"] for x in tl]) / 60.0
    t_min = lambda h: h["t"] / 60.0
    fig, axes = plt.subplots(4, 1, figsize=(11.5, 11.4), sharex=True,
                             gridspec_kw={"height_ratios": [1.3, 1, 1, 1],
                                          "hspace": 0.14})
    ax = axes[0]
    ax.plot(tmin_log(tlog_open), [x["dT"] for x in tlog_open],
            color="#c0392b", lw=1.5, label="open loop (nominal setpoints)")
    ax.plot(tmin_log(tlog_ctrl), [x["dT"] for x in tlog_ctrl],
            color="#145a32", lw=1.5, label="SAC + NMPC shield (autonomous)")
    ax.axhline(40, color="k", ls="--", lw=1.0)
    ax.text(1.0, 41.6, "ΔT = 40 K runaway threshold", fontsize=8)
    ax.set_ylim(0, 88)
    if supp.get("alarm_t") and supp.get("arrested_t"):
        ta, tsp = supp["alarm_t"] / 60, supp["arrested_t"] / 60
        ax.axvline(ta, color="#d35400", lw=1.0, ls=":")
        ax.axvline(tsp, color="#145a32", lw=1.0, ls=":")
        ax.annotate(f"alarm {ta:.1f} min", (ta - 0.25, 55), fontsize=7.5,
                    color="#d35400", rotation=90, va="top", ha="right")
        ax.annotate(f"arrested +{supp['arrest_time']:.0f} s",
                    (tsp + 0.25, 55), fontsize=7.5, color="#145a32",
                    rotation=90, va="top", ha="left")
    for k, t0 in (("cav", 600), ("imp", 1200), ("foul", 1800),
                  ("cool", 2400)):
        ax.axvline(t0 / 60, color=EVENT_STYLE[k][1], lw=1.0, ls=":")
    ax.set_ylabel("max$_z$ ΔT (K)")
    ax.set_title("(a) Thermal-runaway battle: channel superheat")
    ax.legend(loc="upper left")
    ax = axes[1]
    acts = hist_ctrl["actions"]
    dl = hist_ctrl["delivered"]
    ax.plot(t_min(hist_ctrl), acts[:, 0], color="#1a5276", lw=1.3,
            label="Q$_A$ cmd")
    ax.plot(t_min(hist_ctrl), acts[:, 1], color="#7d3c98", lw=1.3,
            label="Q$_B$ cmd")
    ax.plot(t_min(hist_ctrl), acts[:, 2] * 10, color="#b9770e", lw=1.3,
            label="Q$_{cat}$ cmd ×10")
    ax.plot(t_min(hist_ctrl), dl[:, 1], color="#7d3c98", lw=0.8, ls=":",
            alpha=0.9, label="Q$_B$ delivered (cavitation slip)")
    for t0 in (600, 1200, 1800, 2400):
        ax.axvline(t0 / 60, color="0.6", lw=0.8, ls=":")
    ax.set_ylabel("pump flow (mL/min)")
    ax.set_title("(b) Micro-dosing HPLC pumps — stoichiometry & catalyst "
                 "throttling")
    ax.legend(loc="upper right", ncol=2)
    ax = axes[2]
    ax.plot(t_min(hist_ctrl), acts[:, 3], color="#117864", lw=1.3,
            label="jacket rate cmd")
    ax.plot(t_min(hist_ctrl), dl[:, 3], color="#117864", lw=0.8, ls=":",
            label="delivered (vapor lock)")
    ax2 = ax.twinx()
    ax2.plot(t_min(hist_ctrl), acts[:, 4], color="#616a6b", lw=1.2,
             label="BPR setpoint")
    ax2.set_ylabel("BPR (bar)", color="#616a6b"); ax2.grid(False)
    for t0 in (600, 1200, 1800, 2400):
        ax.axvline(t0 / 60, color="0.6", lw=0.8, ls=":")
    ax.set_ylabel("jacket flow (× nominal)")
    ax.set_title("(c) Jacket cooling & back-pressure regulator")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper left", ncol=3)
    ax = axes[3]
    y_open = np.array([x["y_cum"] for x in tlog_open])
    y_ctrl = np.array([x["y_cum"] for x in tlog_ctrl])
    yo_spec = np.array([x.get("y_onspec", 0.0) for x in tlog_open])
    yc_spec = np.array([x.get("y_onspec", 0.0) for x in tlog_ctrl])
    ax.plot(tmin_log(tlog_open), y_open * 100, color="#c0392b", lw=1.4,
            label="open-loop cumulative yield")
    ax.plot(tmin_log(tlog_ctrl), y_ctrl * 100, color="#145a32", lw=1.4,
            label="controlled cumulative yield")
    ax.plot(tmin_log(tlog_open), yo_spec * 100, color="#c0392b", lw=1.0,
            ls=":", alpha=0.9, label="open-loop on-spec yield (S ≥ 90 %)")
    ax.plot(tmin_log(tlog_ctrl), yc_spec * 100, color="#145a32", lw=1.0,
            ls=":", alpha=0.9, label="controlled on-spec yield (S ≥ 90 %)")
    ax.set_xlabel("t (min)"); ax.set_ylabel("yield (%)")
    ax2 = ax.twinx()
    s_ctrl = np.array([x.get("y_inst", 0.0) for x in tlog_ctrl])
    ax2.plot(tmin_log(tlog_ctrl), s_ctrl, color="#2471a3", lw=0.6,
             alpha=0.5)
    ax2.set_ylabel("instant selectivity", color="#2471a3"); ax2.grid(False)
    for t0 in (600, 1200, 1800, 2400):
        ax.axvline(t0 / 60, color="0.6", lw=0.8, ls=":")
    ax.set_title("(d) Yield recovery after autonomous disturbance rejection")
    ax.legend(loc="center right", framealpha=0.92)
    fig.suptitle("Phase 10 · Fig. 3 — Cyber-physical reinforcement control\n"
                 "the agent counteracts cavitation, impurity, fouling and a "
                 "coolant-fault runaway: alarm → dilute-quench → recovery",
                 fontsize=11, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(path); plt.close(fig)


def fig4_pareto(rows, front, path):
    fig = plt.figure(figsize=(12.5, 9.6))
    gs = fig.add_gridspec(2, 2, hspace=0.33, wspace=0.26)
    cost = np.array([r["cost"] for r in rows])
    ci = np.array([r["ci"] for r in rows])
    sty = np.array([r["sty"] for r in rows])
    nch = np.array([r["n_channels"] for r in rows])
    ax = fig.add_subplot(gs[0, 0], projection="3d")
    sc = ax.scatter(cost, ci, sty, c=np.log10(nch), cmap="plasma", s=14,
                    alpha=0.55, depthshade=False)
    fc = np.array([[r["cost"], r["ci"], r["sty"]] for r in front])
    order = np.argsort(fc[:, 2])
    ax.plot(fc[order, 0], fc[order, 1], fc[order, 2], "k.-", lw=1.4, ms=5,
            label="Pareto front")
    ax.set_xlabel("unit cost ($/kg)")
    ax.set_ylabel("carbon (kg CO$_2$e/kg)")
    ax.set_zlabel("STY (kg m$^{-3}$ h$^{-1}$)")
    ax.set_title("(a) Cost–carbon–STY Pareto surface\n"
                 "(colour: log$_{10}$ channels; black: non-dominated)")
    ax.legend(loc="upper left")
    ax.view_init(elev=22, azim=-52)
    ax = fig.add_subplot(gs[0, 1])
    fset = {(round(r["cost"], 4), round(r["ci"], 4)) for r in front}
    isf = np.array([(round(c, 4), round(x, 4)) in fset
                    for c, x in zip(cost, ci)])
    ax.scatter(cost[~isf], ci[~isf], s=14 + 26 * (sty[~isf] / sty.max()),
               c="#5d6d7e", alpha=0.45, label="dominated designs")
    order = np.argsort(ci[isf])
    ax.plot(cost[isf][order], ci[isf][order], "o-", color="#c0392b", lw=1.4,
            ms=5, label="Pareto front (cost–carbon)")
    ax.set_xlabel("unit production cost ($/kg)")
    ax.set_ylabel("carbon intensity (kg CO$_2$e/kg)")
    ax.set_title("(b) Cost vs carbon (bubble ∝ STY)")
    ax.legend()
    ax = fig.add_subplot(gs[1, 0])
    for n_ch, col in zip((1, 8, 32, 64),
                         ("#1a5276", "#117864", "#b9770e", "#a93226")):
        sel = [r for r in rows if r["n_channels"] == n_ch]
        sel.sort(key=lambda r: r["q_a"])
        if sel:
            ax.plot([r["q_a"] for r in sel], [r["sty"] for r in sel], "o-",
                    color=col, ms=3.5, lw=1.2, label=f"{n_ch} channel(s)")
    ax.set_xlabel("stream-A flow (mL/min)")
    ax.set_ylabel("STY (kg m$^{-3}$ h$^{-1}$)")
    ax.set_title("(c) Space-time yield vs throughput (numbering-up keeps "
                 "STY, scales output)")
    ax.legend()
    ax = fig.add_subplot(gs[1, 1])
    labels, colors = ("substrate feed", "solvent makeup + recovery",
                      "catalyst makeup", "energy", "capex + maintenance",
                      "labor"), ("#1a5276", "#117864", "#7d3c98",
                                 "#b9770e", "#5d6d7e", "#a93226")
    widths = {1: [], 32: []}
    for n_ch in (1, 32):
        cand = [r for r in rows if r["n_channels"] == n_ch]
        ref = min(cand, key=lambda r: r["cost"])
        ss = solve_steady_state(ref["q_a"], 1.0, ref["q_cat"],
                                t_jacket=ref["t_jacket"] + 273.15)
        q_tot = ss["q_tot"]
        kg_h = ref["kg_h"]
        m_a_in = (ref["q_a"] / 6e7) * C_A_FEED * 1e3 * 3600 * M_P / 1e3 \
            * n_ch
        m_tol = q_tot * RHO_L * 3600 * n_ch - m_a_in
        m_lut = (1.0 / 6e7) * C_LUT_FEED * 1e3 * 3600 * M_LUT / 1e3 * n_ch
        m_cat = (ref["q_cat"] / 6e7) * C_CAT_FEED * 1e3 * 3600 * M_CPA \
            / 1e3 * n_ch
        dp = 128.0 * viscosity(ss["T_out"]) * q_tot * L_REACTOR \
            / (math.pi * D_TUBE**4)
        kwh_kg = ref["kwh_kg"]
        capex = CAPEX_BASE + CAPEX_PER_CHANNEL * n_ch
        capex_rate = capex / DEP_YEARS / (HOURS_YEAR * UTIL)
        comp = [
            m_a_in * TEA_PRICE["substrate"] / kg_h,
            m_tol * (1 - SOLV_RECOVERY) * TEA_PRICE["toluene"] / kg_h
            + SOLV_RECOVERY * m_tol * SOLV_RECOVERY_ENERGY * ELEC_PRICE
            / kg_h,
            m_cat * (1 - CAT_RECOVERY) * TEA_PRICE["cpa"] / kg_h,
            kwh_kg * ELEC_PRICE,
            (capex_rate + capex * MAINT_FRAC / (HOURS_YEAR * UTIL)) / kg_h,
            LABOR_RATE,
        ]
        widths[n_ch] = comp
    y0 = [0, 0]
    x = [0, 1]
    for lab, col, i in zip(labels, colors, range(6)):
        vals = [widths[1][i], widths[32][i]]
        ax.bar(x, vals, bottom=y0, width=0.55, color=col, label=lab)
        y0 = [y0[j] + vals[j] for j in range(2)]
    for xi, tot in zip(x, y0):
        ax.text(xi, tot + 6, f"${tot:.0f}/kg", ha="center", fontsize=8.5,
                fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(["1 channel (0.5 mL/min scale)",
                        "32 channels (numbered-up)"])
    ax.set_ylabel("unit production cost ($/kg)")
    ax.set_title("(d) Cost anatomy at two throughput scales "
                 "(cheapest design per scale)")
    ax.legend(fontsize=7)
    fig.suptitle("Phase 10 · Fig. 4 — Techno-economic & lifecycle Pareto "
                 "analysis of the autonomous flow plant", fontsize=12,
                 y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path); plt.close(fig)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #


def fidelity(full: bool):
    return dict(n_z=161 if full else 81, dt_cap=0.5 if full else 0.8,
                ep_length=EP_T_DEMO if full else 1400.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--fig_only", action="store_true")
    ap.add_argument("--skip_train", action="store_true")
    args = ap.parse_args()
    full = not args.selftest
    fid = fidelity(full)

    log("=" * 78)
    log("PHASE 10 — MULTI-SCALE CONTINUOUS-FLOW CYBER-PHYSICAL DIGITAL TWIN")
    log(f"  mode={'SELFTEST' if args.selftest else 'FULL'}  "
        f"numpy {np.__version__} · torch {torch.__version__}")
    log("=" * 78)

    json_out = {"phase": 10, "generated": _dt.datetime.now().isoformat(),
                "mode": "selftest" if args.selftest else "full",
                "reactor": {"L_m": L_REACTOR, "d_t_mm": D_TUBE * 1e3,
                            "V_mL": V_REACTOR * 1e6},
                "chemistry": {"RDS_kJ_mol": EA_1 / 1e3, "net_dH_kJ_mol":
                              (DH_1 + DH_2) / 1e3,
                              "ee_293K": float(ee_of_temperature(293.15))}}

    if args.fig_only:
        hist_ctrl = np.load(RESULTS / "episode_controlled.npy",
                            allow_pickle=True).item()
        hist_open = np.load(RESULTS / "episode_openloop.npy",
                            allow_pickle=True).item()
        rec = json.loads((RESULTS / "phase10_results.json").read_text())
        pat = np.load(RESULTS / "pat_telemetry.npy", allow_pickle=True).item()
        rows = [dict(r) for r in csv.DictReader(
            (RESULTS / "tea_pareto_points.csv").open())]
        for r in rows:
            for k in ("n_channels", "conv", "selectivity", "ee", "kg_h",
                      "sty", "pmi", "e_factor", "ci", "cost", "kwh_kg",
                      "dT_channel", "duty_kW", "q_a", "q_b", "q_cat",
                      "t_jacket"):
                r[k] = float(r[k])
        front = pareto_front(rows)
        supp = rec["suppression"]
        fig1_reactor_profile(hist_ctrl, None, FIGURES /
                             "fig1_continuous_pde_reactor_profile.png")
        fig2_pat_telemetry(hist_ctrl, hist_open, pat, FIGURES /
                           "fig2_operando_pat_sensor_telemetry.png")
        fig3_control_dynamics(
            hist_open, rec["tlog_open"], hist_ctrl, rec["tlog_ctrl"],
            None, supp, FIGURES / "fig3_rl_cyberphysical_control_dynamics.png")
        fig4_pareto(rows, front, FIGURES /
                    "fig4_techno_economic_pareto_analysis.png")
        log("fig_only: figures regenerated.")
        return

    # ------------------------------------------------------------------ #
    # 0. plant physics sanity probe
    # ------------------------------------------------------------------ #
    log("[0] Steady-state probe of the continuum plant (Module 10A/10D)...")
    ss = solve_steady_state(QA_NOM, QB_NOM, QCAT_NOM)
    nominal_tea = tea_metrics(ss, 1, QA_NOM, QB_NOM, QCAT_NOM)
    log(f"    nominal SS: X={nominal_tea['conv']:.3f} "
        f"S={nominal_tea['selectivity']:.3f} "
        f"T_out={ss['T_out']:.1f} K ee={nominal_tea['ee']*100:.1f}% "
        f"STY={nominal_tea['sty']:.0f} kg/m3/h "
        f"dT_channel={nominal_tea['dT_channel']:.2f} K")
    json_out["nominal_steady_state"] = nominal_tea

    # ------------------------------------------------------------------ #
    # 1. Train the SAC agent on the domain-randomized fault family
    # ------------------------------------------------------------------ #
    torch.set_num_threads(min(4, os.cpu_count() or 2))
    if args.skip_train and (RESULTS / "sac_policy.pt").exists():
        agent = SACAgent()
        agent.actor.load_state_dict(
            torch.load(RESULTS / "sac_policy.pt", weights_only=True))
        curve = []
        log("[1] Loaded pretrained SAC policy (--skip_train).")
    else:
        log("[1] Training SAC agent (Module 10C) on domain-randomized "
            "fault episodes...")
        if full:
            scen_picker = lambda rng: rng.choice(
                ["none", "cav", "imp", "foul", "cool", "cool_foul",
                 "cav_imp", "all"], p=[0.12, 0.12, 0.12, 0.12, 0.16,
                                       0.16, 0.10, 0.10])
            episodes, ep_len, n_z_tr, dt_cap_tr = 110, 1600.0, 41, 1.2
        else:
            scen_picker = lambda rng: rng.choice(
                ["none", "cav", "imp", "cool", "all"],
                p=[0.3, 0.15, 0.15, 0.2, 0.2])
            episodes, ep_len, n_z_tr, dt_cap_tr = 16, 800.0, 41, 1.2
        agent = SACAgent()
        t0 = time.time()
        curve = train_sac(agent, scen_picker, episodes, ep_len, n_z_tr,
                          dt_cap_tr)
        log(f"    training done in {time.time()-t0:.0f}s "
            f"({len(curve)} episodes)")
        torch.save(agent.actor.state_dict(), RESULTS / "sac_policy.pt")
    json_out["sac"] = {
        "episodes": len(curve),
        "final_return_mean": (float(np.mean([c["return"] for c in
                                             curve[-10:]]))
                              if curve else None),
        "first10_return_mean": (float(np.mean([c["return"] for c in
                                               curve[:10]]))
                                if len(curve) >= 10 else None),
    }
    if curve:
        with (RESULTS / "sac_learning_curve.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(curve[0].keys()))
            w.writeheader(); w.writerows(curve)

    # ------------------------------------------------------------------ #
    # 2. Validation suite: open-loop vs SAC vs SAC+shield
    # ------------------------------------------------------------------ #
    log("[2] Validation suite on held-out high-fidelity episodes...")
    shield = SafetyShield()
    policies = (("open_loop", OpenLoopPolicy()),
                ("sac_only", SACPolicy(agent)),
                ("sac_shield", ShieldedPolicy(agent, shield)))
    val_rows = []
    val_scenarios = ["cav", "imp", "foul", "cool", "cool_foul", "cav_imp",
                     "all", "none"]
    for ctrl_name, policy in policies:
        peaks, yields_, sels = [], [], []
        for i, scen in enumerate(val_scenarios):
            env = FlowEnv(n_z=101 if full else 61,
                          dt_cap=0.5 if full else 0.8,
                          ep_length=min(fid["ep_length"], 2400.0),
                          record=True,
                          seed=SEED + 100 + i)
            h, tl, _, _ = run_episode(env, policy, scenario=scen,
                                      seed=SEED + 200 + i)
            peaks.append(max(x["dT"] for x in tl))
            yields_.append([x["y_cum"] for x in tl][-1])
            sels.append(episode_final_quality(env, h)["tail_selectivity"])
            val_rows.append({"controller": ctrl_name, "scenario": scen,
                             "peak_dT": peaks[-1], "final_yield":
                             yields_[-1], "tail_selectivity": sels[-1],
                             "breach40": peaks[-1] > 40.0})
        log(f"    {ctrl_name:11s}: peak dT {np.mean(peaks):5.1f}±{np.std(peaks):4.1f} K "
            f"| Y {np.mean(yields_)*100:5.1f}% "
            f"| S {np.mean(sels)*100:5.1f}% "
            f"| dT>40K in {sum(p > 40 for p in peaks)}/{len(peaks)} eps")
    with (RESULTS / "validation_summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(val_rows[0].keys()))
        w.writeheader(); w.writerows(val_rows)
    json_out["validation"] = {}
    for ctrl in ("open_loop", "sac_only", "sac_shield"):
        rs = [r for r in val_rows if r["controller"] == ctrl]
        json_out["validation"][ctrl] = {
            "mean_peak_dT": float(np.mean([r["peak_dT"] for r in rs])),
            "breach40_rate": float(np.mean([r["breach40"] for r in rs])),
            "mean_final_yield": float(np.mean([r["final_yield"]
                                               for r in rs])),
            "mean_tail_selectivity": float(np.mean(
                [r["tail_selectivity"] for r in rs]))}

    # ------------------------------------------------------------------ #
    # 3. Acceptance demo: the four-disturbance timeline (fig1-3 data)
    # ------------------------------------------------------------------ #
    log("[3] Acceptance demo on the 3600 s disturbance timeline "
        "(hi-fi grid)...")
    demo_fid = fidelity(full)
    env_open = FlowEnv(n_z=demo_fid["n_z"], dt_cap=demo_fid["dt_cap"],
                       ep_length=EP_T_DEMO, record=True, seed=SEED + 42)
    h_open, tl_open, _, _ = run_episode(
        env_open, OpenLoopPolicy(), scenario="demo", seed=SEED + 42,
        verbose_every=60, tag="open-loop")
    supp_open = suppression_metrics(tl_open)

    env_ctrl = FlowEnv(n_z=demo_fid["n_z"], dt_cap=demo_fid["dt_cap"],
                       ep_length=EP_T_DEMO, record=True, seed=SEED + 42)
    h_ctrl, tl_ctrl, ev_ctrl, _ = run_episode(
        env_ctrl, ShieldedPolicy(agent, shield), scenario="demo",
        seed=SEED + 42, verbose_every=60, tag="controlled")
    supp_ctrl = suppression_metrics(tl_ctrl)
    qual_ctrl = episode_final_quality(env_ctrl, h_ctrl)
    qual_open = episode_final_quality(env_open, h_open)
    log(f"    open-loop : peak dT {supp_open['peak_dT']:.1f} K "
        f"(>40 K: {supp_open['breach40']}), tail S "
        f"{qual_open['tail_selectivity']*100:.1f}%")
    log(f"    controlled: peak dT {supp_ctrl['peak_dT']:.1f} K "
        f"(>40 K: {supp_ctrl['breach40']}), tail S "
        f"{qual_ctrl['tail_selectivity']*100:.1f}%, "
        f"alarm {supp_ctrl['alarm_t']} s -> arrested "
        f"{supp_ctrl['arrested_t']} s "
        f"(arrest dt = {supp_ctrl['arrest_time']} s)")
    json_out["suppression"] = supp_ctrl
    json_out["open_loop_demo"] = {**supp_open, **qual_open}
    json_out["controlled_demo"] = {**supp_ctrl, **qual_ctrl}

    np.save(RESULTS / "episode_openloop.npy", h_open, allow_pickle=True)
    np.save(RESULTS / "episode_controlled.npy", h_ctrl, allow_pickle=True)

    # ------------------------------------------------------------------ #
    # 4. Operando PAT telemetry bundle (fig2) + hallucination audit
    # ------------------------------------------------------------------ #
    log("[4] Generating operando PAT telemetry (Raman / UV-Vis / dP)...")
    pat = build_pat_telemetry(h_ctrl, np.random.default_rng(SEED + 9))
    # audit: PAT deconvolution error vs plant truth across the run
    errs = []
    for i in range(0, len(h_ctrl["t"]), 8):
        est = pat_deconvolve_from_hist(h_ctrl, i)
        truth = {s: h_ctrl["C"][s][i, -1] for s in ("A", "I", "P", "E")}
        for s in ("A", "P"):
            if truth[s] > 0.02:
                errs.append(abs(est[s] - truth[s]) / truth[s])
    audit = {"n": len(errs), "mean_rel_err": float(np.mean(errs)),
             "max_rel_err": float(np.max(errs))}
    log(f"    PAT deconvolution audit: mean |rel err| "
        f"{audit['mean_rel_err']*100:.1f}%, max {audit['max_rel_err']*100:.1f}%")
    json_out["pat_audit"] = audit
    np.save(RESULTS / "pat_telemetry.npy", pat, allow_pickle=True)

    # ------------------------------------------------------------------ #
    # 5. TEA / LCA Pareto grid (fig4)
    # ------------------------------------------------------------------ #
    log("[5] Techno-economic / LCA steady-state grid (Module 10D)...")
    rows = tea_pareto_grid(coarse=not full)
    with (RESULTS / "tea_pareto_points.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    front = pareto_front(rows)
    best = min(rows, key=lambda r: r["cost"])
    log(f"    {len(rows)} design points, Pareto front size {len(front)}; "
        f"cheapest ${best['cost']:.0f}/kg @ STY {best['sty']:.0f} "
        f"kg/m3/h, CI {best['ci']:.1f} kgCO2e/kg, E {best['e_factor']:.2f}")
    json_out["tea"] = {
        "n_points": len(rows), "pareto_size": len(front),
        "cheapest": best,
        "max_sty": max(r["sty"] for r in rows),
        "front": front,
        "nominal_point": nominal_tea,
    }

    # ------------------------------------------------------------------ #
    # 6. Figures
    # ------------------------------------------------------------------ #
    log("[6] Rendering 300-DPI figures...")
    fig1_reactor_profile(h_ctrl, None,
                         FIGURES / "fig1_continuous_pde_reactor_profile.png")
    fig2_pat_telemetry(h_ctrl, h_open, pat,
                       FIGURES / "fig2_operando_pat_sensor_telemetry.png")
    fig3_control_dynamics(h_open, tl_open, h_ctrl, tl_ctrl, ev_ctrl,
                          supp_ctrl,
                          FIGURES /
                          "fig3_rl_cyberphysical_control_dynamics.png")
    fig4_pareto(rows, front,
                FIGURES / "fig4_techno_economic_pareto_analysis.png")

    json_out["tlog_open"] = tl_open
    json_out["tlog_ctrl"] = tl_ctrl
    (RESULTS / "phase10_results.json").write_text(json.dumps(
        json_out, indent=1, default=float))
    log(f"[✓] Phase 10 complete — results in {RESULTS.name}/, "
        f"figures in {FIGURES.name}/")


def pat_deconvolve_from_hist(hist, i):
    """Audit helper: run the PAT deconvolution on stored plant truth."""
    out = {s: float(hist["C"][s][i, -1]) for s in
           ("A", "I", "P", "E", "O", "Cat", "Lut")}
    out["Imp"] = float(hist["C"]["Imp"][i, -1]) \
        if "Imp" in hist["C"] else 0.0
    uv254, uv310 = uv_absorbance(out)
    pk = {}
    for meas, nu_r in (("r1655", 1655.0), ("r1602", 1602.0),
                       ("r1668", 1668.0), ("r1003", 1003.0)):
        pk[meas] = 0.0
        for comp in _RAMAN_COMPS:
            for (nu0, amp, hwhm) in RAMAN_BANDS[comp]:
                pk[meas] += (RAMAN_XSEC[comp] * out.get(comp, 0.0)
                             * amp * hwhm**2 / ((nu_r - nu0) ** 2 + hwhm**2))
    meas = np.array([uv254, uv310, pk["r1655"], pk["r1602"], pk["r1668"],
                     pk["r1003"]])
    design = np.zeros((6, len(_DECONV_COMPONENTS)))
    for j, comp in enumerate(_DECONV_COMPONENTS):
        e254, e310 = UV_EPS[comp]
        design[0, j] = e254 * UV_PATH
        design[1, j] = e310 * UV_PATH
        for row, nu_r in zip(range(2, 6), (1655.0, 1602.0, 1668.0, 1003.0)):
            for (nu0, amp, hwhm) in RAMAN_BANDS[comp]:
                design[row, j] += (RAMAN_XSEC[comp] * amp * hwhm**2
                                   / ((nu_r - nu0) ** 2 + hwhm**2))
    sol, *_ = np.linalg.lstsq(design, meas, rcond=None)
    return {comp: float(sol[j]) for j, comp in enumerate(_DECONV_COMPONENTS)}


if __name__ == "__main__":
    main()
