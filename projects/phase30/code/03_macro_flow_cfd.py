"""Conservative 1-D finite-volume surrogate, explicitly NOT a 3-D CFD solver.

Liquid is plug-flow; gas transfer is capped by remaining pure-CO2 feed.
Products leave separately; gas holdup, salt precipitation and flooding unresolved.
"""
import math
from scipy.optimize import brentq
from common import PRODUCTS, CARBON, ELECTRONS, F, module, LABEL

meso = module("02_meso_microkinetics.py")


def channel(config):
    config.validate()
    q = config.velocity_m_s * config.width_m * config.height_m
    area = config.width_m * config.length_m
    dx = config.length_m/config.cells
    area_cell = config.width_m*dx
    volume_cell = area_cell*config.height_m
    c, gas, temperature = config.co2_inlet_mol_m3, config.gas_feed_mol_s, config.temperature_K
    products = dict.fromkeys(PRODUCTS, 0.)
    carbonate = power = heat_generated = heat_removed = 0.
    rows = []
    for i in range(config.cells):
        local_config = config.with_changes(temperature_K=temperature)
        def evaluate(outlet):
            s = meso.solve(local_config, outlet)
            carbon = sum(CARBON[k]*s["flux_mol_m2_s"][k]*area_cell for k in PRODUCTS)
            loss = config.carbonate_loss_s * outlet * volume_cell
            transfer = min(gas, config.kla_s * volume_cell *
                           (config.co2_equilibrium_mol_m3-outlet))
            return q*(outlet-c) - transfer + carbon + loss, s, transfer, loss
        upper = config.co2_equilibrium_mol_m3
        # Zero is always a lower bracket for physical nonnegative inlet supply.
        outlet = brentq(lambda x: evaluate(x)[0], 0., upper, xtol=1e-11)
        residual, s, transfer, loss = evaluate(outlet)
        gas -= transfer
        if gas < -1e-15 or outlet < 0 or abs(residual) > 1e-12:
            raise ArithmeticError("Cell carbon balance / nonnegativity failure")
        j = s["current_A_m2"]
        voltage = 1.229 + .35 - config.potential_RHE_V + j*config.area_resistance_ohm_m2
        electric = j*voltage*area_cell
        heat = config.heat_fraction*electric
        # Implicit convection + isothermal wall heat exchange.
        flow_capacity = q*997.*4180.  # W/K
        ua = config.heat_transfer_W_m2_K*area_cell
        next_t = (flow_capacity*temperature + ua*config.temperature_K + heat)/(flow_capacity+ua)
        removed = ua*(next_t-config.temperature_K)
        if next_t >= config.max_temperature_K or j > config.max_current_A_m2:
            raise ValueError("Temperature/current safety envelope exceeded")
        for key in PRODUCTS:
            products[key] += s["flux_mol_m2_s"][key]*area_cell
        power += electric
        heat_generated += heat
        heat_removed += removed
        carbonate += loss
        rows.append(dict(cell=i, x_m=(i+1)*dx, dissolved_CO2_mol_m3=outlet,
                         gas_CO2_mol_s=gas, temperature_K=next_t, local_pH=s["local_pH"],
                         current_A_m2=j, voltage_V=voltage, theta_CO=s["theta_CO"],
                         carbon_residual_mol_s=residual))
        c, temperature = outlet, next_t
    carbon_in = config.gas_feed_mol_s + q*config.co2_inlet_mol_m3
    product_carbon = sum(CARBON[k]*products[k] for k in PRODUCTS)
    carbon_out = gas + q*c + carbonate + product_carbon
    current = F*sum(ELECTRONS[k]*products[k] for k in PRODUCTS)
    dh = 2*config.width_m*config.height_m/(config.width_m+config.height_m)
    reynolds = 997.*config.velocity_m_s*dh/.00089
    # Rectangular channel approximation for h<=w; single-liquid equivalent only.
    pressure_drop = (12*.00089*config.length_m*q /
                     (config.width_m*config.height_m**3*(1-.63*config.height_m/config.width_m)))
    if reynolds > 1000 or pressure_drop > config.max_pressure_drop_Pa:
        raise ValueError("Laminar/pressure closure outside supported envelope")
    heat_residual = heat_generated-heat_removed-q*997.*4180.*(temperature-config.temperature_K)
    carbon_residual = carbon_in-carbon_out
    if abs(carbon_residual) > max(1e-13, carbon_in*1e-7) or abs(heat_residual) > 1e-8:
        raise ArithmeticError("Global balance failed")
    return dict(evidence_type=LABEL, model="1D_conservative_gas_liquid_solid_boundary_surrogate",
                cfd_3d_executed=False, area_m2=area, liquid_flow_m3_s=q,
                products_mol_s=products, carbonate_mol_s=carbonate,
                carbon_in_mol_s=carbon_in, carbon_out_mol_s=carbon_out,
                carbon_residual_mol_s=carbon_residual, carbon_utilization=product_carbon/carbon_in,
                FE={k: F*ELECTRONS[k]*products[k]/current for k in PRODUCTS},
                current_A=current, mean_current_A_m2=current/area,
                power_W=power, effective_voltage_V=power/current,
                heat_generated_W=heat_generated, heat_removed_W=heat_removed,
                sensible_heat_W=q*997.*4180.*(temperature-config.temperature_K),
                heat_residual_W=heat_residual, max_temperature_K=max(r["temperature_K"] for r in rows),
                reynolds=reynolds, pressure_drop_Pa=pressure_drop,
                pump_power_W=pressure_drop*q/.6, cells=rows)
