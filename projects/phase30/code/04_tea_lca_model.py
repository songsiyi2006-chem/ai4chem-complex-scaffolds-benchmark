"""Discounted target-product cost and undiscounted screening GHG inventory.

All shared burdens assigned to one chosen sale product; no coproduct revenue,
CO2 uptake credit, subsidy, thermal integration or avoided-product credit.
"""
import math
from common import MW_KG, LABEL


def lifetime(config, hours=None):
    """Constant-load duty cycles: decay reduces production, not electricity.

    Cycle resets at replacement, downtime consumes available scheduled hours.
    Event-driven within annual boundaries; no replacement at project termination.
    """
    scheduled = 8760*config.capacity_factor
    horizon = scheduled*config.years if hours is None else hours
    k = config.degradation_per_operating_h
    life = -math.log(config.replacement_threshold)/k if k else math.inf
    if life < 1:
        raise ValueError("Unphysical replacement interval below 1 operating hour")
    elapsed = age = downtime_left = 0.
    records, annual = [], []
    for year in range(1, config.years+1):
        boundary = min(year*scheduled, horizon)
        operating = effective = downtime = 0.
        replacements = 0
        while elapsed < boundary - 1e-9:
            if downtime_left > 0:
                dt = min(downtime_left, boundary-elapsed)
                downtime += dt; elapsed += dt; downtime_left -= dt
                continue
            dt = min(life-age, boundary-elapsed)
            effective += (math.exp(-k*age)*(-math.expm1(-k*dt))/k) if k else dt
            operating += dt; elapsed += dt; age += dt
            records.append(dict(year=year, scheduled_h=elapsed,
                                activity=math.exp(-k*age), event="operating_segment_end"))
            if age >= life-1e-8 and elapsed < horizon-1e-8:
                replacements += 1
                age = 0.; downtime_left = config.replacement_downtime_h
                records.append(dict(year=year, scheduled_h=elapsed, activity=1., event="replacement"))
        annual.append(dict(year=year, operating_h=operating, effective_product_h=effective,
                           downtime_h=downtime, replacements=replacements))
        if elapsed >= horizon-1e-8:
            break
    return annual, records


def assess(config, flow, stacks=None):
    stacks = config.stacks if stacks is None else stacks
    if type(stacks) is not int or stacks < 1:
        raise ValueError("Positive integer stack count required")
    rated_kW = config.stack_kW*stacks
    replicas = rated_kW*1000/flow["power_W"]
    output_kg_h = flow["products_mol_s"][config.target]*replicas*MW_KG[config.target]*3600*config.recovery_fraction
    if output_kg_h <= 1e-12:
        raise ValueError("No target product; cost not defined")
    pump_kW = flow["pump_power_W"]*replicas/1000.
    # Purchased gas + dissolved feed equivalent, no recycle credits.
    feed_kg_h = flow["carbon_in_mol_s"]*replicas*.04401*3600
    capital = rated_kW*config.capex_USD_kW
    costs = dict(capex=capital, electricity=0., feed_CO2=0., replacements=0., fixed_opex=0.)
    emissions = dict(electricity=0., feed_CO2=0., embodied= rated_kW*config.embodied_kgCO2e_kW,
                     replacements=0.)
    discounted_kg = undiscounted_kg = consumed_kWh = 0.
    annual, trace = lifetime(config)
    for row in annual:
        product = output_kg_h*row["effective_product_h"]
        electricity = ((rated_kW+pump_kW)*row["operating_h"]
                       + config.separation_kWh_kg*product)
        feed = feed_kg_h*row["operating_h"]
        discount = (1+config.discount_rate)**row["year"]
        row.update(product_kg=product, electricity_kWh=electricity, purchased_CO2_kg=feed,
                   discount_factor=discount)
        costs["electricity"] += electricity*config.electricity_USD_kWh/discount
        costs["feed_CO2"] += feed*config.co2_USD_kg/discount
        costs["replacements"] += row["replacements"]*capital*config.replacement_fraction/discount
        costs["fixed_opex"] += capital*config.fixed_opex_fraction_year/discount
        emissions["electricity"] += electricity*config.grid_kgCO2e_kWh
        emissions["feed_CO2"] += feed*config.feed_kgCO2e_kg
        emissions["replacements"] += (row["replacements"]*rated_kW*config.embodied_kgCO2e_kW
                                      *config.replacement_fraction)
        discounted_kg += product/discount
        undiscounted_kg += product
        consumed_kWh += electricity
    return dict(evidence_type=LABEL, target=config.target, installed_power_kW=rated_kW, stacks=stacks,
                stack_kW=config.stack_kW, geometric_area_m2=replicas*flow["area_m2"],
                idealized_channel_equivalents=replicas, initial_recovered_kg_h=output_kg_h,
                functional_unit="1 kg recovered target at gate; all shared burdens to target",
                monetary_basis="constant illustrative USD; no market-price or inflation forecast",
                discounted_cost_USD=costs, discounted_product_kg=discounted_kg,
                LCO_product_USD_kg=sum(costs.values())/discounted_kg,
                LCOE_USD_kWh=None, LCOE_reason="electrolyzer consumes electricity, does not generate it",
                GHG_inventory_kgCO2e=emissions, product_lifetime_kg=undiscounted_kg,
                screening_GHG_kgCO2e_kg=sum(emissions.values())/undiscounted_kg,
                electricity_kWh_kg=consumed_kWh/undiscounted_kg,
                total_replacements=sum(r["replacements"] for r in annual),
                annual=annual, degradation_trace=trace, industrial_lifetime_validated=False,
                limitations=["10-stack linear replication excludes manifolds and scale-dependent failure",
                             "activity loss diverts lost target output to unspecified side products at fixed load",
                             "GHG screening omits full equipment/electrolyte/refrigeration LCI; not ISO-complete LCA",
                             "no uptake credit: feedstock carbon is not permanent CO2 removal"])
