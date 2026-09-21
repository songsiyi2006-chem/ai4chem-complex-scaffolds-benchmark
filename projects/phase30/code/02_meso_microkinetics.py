"""Single-coverage stiff ODE + calibrated-nowhere Stern/diffuse-layer closure."""
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from common import F, R, ELECTRONS, PRODUCTS, module

micro = module("01_micro_surface.py")


def edl(config, local_pH, current_A_m2):
    eps = 78.4 * 8.8541878128e-12
    debye = math.sqrt(eps * R * config.temperature_K /
                     (2 * F**2 * config.ionic_strength_mol_m3))
    diffuse_cap = eps / debye
    stern_cap = .2  # F/m2, assumption
    u_she = micro.rhe_to_she(config.potential_RHE_V, config.bulk_pH, config.temperature_K)
    psi = (u_she - (-.5)) * stern_cap / (stern_cap + diffuse_cap)
    # Alkalinity-generation/film/buffer balance + Boltzmann proton partition.
    predicted_pH = (config.bulk_pH + current_A_m2 /
                    (F * config.liquid_film_m_s * config.buffer_capacity_mol_m3_pH)
                    + F * psi / (math.log(10) * R * config.temperature_K))
    return psi, debye, predicted_pH


def state_at_pH(config, c, local_pH, validate_ode=False):
    psi, debye, _ = edl(config, local_pH, 0.)
    k, _ = micro.rates(config, psi, local_pH)
    a = k["ads"] * c / 33.
    b = a + k["CO"]
    theta = 0. if a == 0 else 2*a / (b + math.sqrt(b*b + 8*k["C2H4"]*a))
    rhs = lambda t, y: np.array([a - b*y[0] - 2*k["C2H4"]*y[0]**2])
    derivative = abs(float(rhs(0, [theta])[0]))
    ode = None
    if validate_ode:
        horizon = 40 / max(b + 4*k["C2H4"]*theta, 1e-12)
        sol = solve_ivp(rhs, [0, horizon], [0.], method="BDF", rtol=1e-9, atol=1e-11)
        if not sol.success or abs(sol.y[0, -1] - theta) > 1e-7 or sol.y.min() < -1e-9:
            raise ArithmeticError("ODE failed steady-state/nonnegativity crosscheck")
        ode = dict(method="BDF", success=True, nfev=sol.nfev,
                   final_coverage=float(sol.y[0, -1]), analytic_error=abs(float(sol.y[0, -1])-theta))
    flux = dict(CO=k["CO"]*theta, C2H4=k["C2H4"]*theta**2,
                HCOOH=k["HCOOH"]*c/33.*(1-theta), H2=k["H2"]*(1-theta))
    flux = {key: value*config.site_density_mol_m2 for key, value in flux.items()}
    partial = {key: F*ELECTRONS[key]*flux[key] for key in PRODUCTS}
    current = sum(partial.values())
    if not 0 <= theta <= 1 or derivative > 1e-7 * max(1., a):
        raise ArithmeticError("Coverage/site balance failed")
    return dict(theta_CO=theta, site_residual_s=derivative, local_pH=local_pH,
                diffuse_potential_V=psi, debye_length_m=debye, current_A_m2=current,
                flux_mol_m2_s=flux, partial_current_A_m2=partial,
                FE={key: value/current for key, value in partial.items()} if current else {},
                ode_validation=ode)


def solve(config, concentration_mol_m3, validate_ode=False):
    if not math.isfinite(concentration_mol_m3) or concentration_mol_m3 < 0:
        raise ValueError("CO2 concentration must be finite and nonnegative")
    def residual(pH):
        s = state_at_pH(config, concentration_mol_m3, pH)
        return pH - edl(config, pH, s["current_A_m2"])[2]
    lo, hi = residual(0.), residual(14.)
    if lo * hi > 0:
        raise ValueError("EDL/buffer solution outside pH 0..14; candidate infeasible")
    pH = brentq(residual, 0., 14., xtol=1e-10)
    result = state_at_pH(config, concentration_mol_m3, pH, validate_ode)
    result["pH_balance_residual"] = residual(pH)
    return result
