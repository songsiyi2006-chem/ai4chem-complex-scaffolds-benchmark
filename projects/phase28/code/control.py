"""Causal finite-grid receding-horizon control and explicit scenario accounting."""
from itertools import product
import numpy as np
from physics import observations, product_kg


def feasible(x, previous, inventory, remaining, available=1.25):
    # inventory is deviation from a common initial buffer of 2 rated charge-hours.
    # Physical storage is 2+inventory, bounded [0,4]; terminal deviation is zero.
    new_inventory = inventory+x-1
    return (0.75-1e-10 <= x <= available+1e-10 and abs(x-previous) <= 0.25+1e-10
            and abs(new_inventory) <= 2+1e-10 and abs(new_inventory) <= 0.25*remaining+1e-10)


def choose_action(strategy, history_prices, history_available, previous, inventory, remaining):
    """Only current/past external signals accepted. Future test values unavailable.

    Forecast is current price in first slot and trailing mean thereafter. Availability
    forecast is persistence. Objective is a proxy; measured site losses are unknown.
    """
    current_price, available = history_prices[-1], history_available[-1]
    choices = [x for x in (0.75,1.0,1.25) if feasible(x,previous,inventory,remaining,available)]
    if strategy == "constant": return 1.0
    if strategy == "rule":
        reference = np.mean(history_prices[-12:])
        preference = 1.25 if current_price < 0.95*reference else (0.75 if current_price > 1.05*reference else 1.0)
        return min(choices, key=lambda x:(abs(x-preference),abs(inventory+x-1)))
    if strategy != "mpc": raise ValueError(strategy)
    horizon = min(4, remaining+1)
    prediction = [current_price]+[float(np.mean(history_prices[-12:]))]*(horizon-1)
    best, action = np.inf, choices[0]
    for sequence in product((0.75,1.0,1.25), repeat=horizon):
        inv, prev, cost, ok = inventory, previous, 0.0, True
        for i,x in enumerate(sequence):
            if not feasible(x,prev,inv,remaining-i,available): ok=False; break
            # 30 A reference module; 3 V intercept + 0.12 V load term.
            cost += 0.030*x*(3+0.12*x)*prediction[i]+0.002*x*x
            inv += x-1; prev=x
        if not ok: continue
        # Value inventory at forecast mean; finite buffer discourages extreme bias.
        cost += 0.0008*inv*inv - inv*0.030*3.12*np.mean(prediction)
        if cost < best: best,action=cost,sequence[0]
    return action


def schedule(strategy, prices, availability):
    """Return equal-production loads with the same initial/final 2 h product buffer."""
    if len(prices)!=len(availability) or len(prices)==0:
        raise ValueError("Equal-length nonempty external sequences required")
    if np.any(np.asarray(availability)<1) or np.any(np.asarray(availability)>1.25):
        raise ValueError("This benchmark requires guaranteed baseload availability between 1 and 1.25")
    load, inventory, previous = [], 0.0, 1.0
    for i in range(len(prices)):
        x = choose_action(strategy,prices[:i+1],availability[:i+1],previous,inventory,len(prices)-i-1)
        load.append(x); inventory += x-1; previous=x
    return np.asarray(load)


def account(load, prices, k=0.0045, layer_cost=0.2, electricity_scale=1.0, fixed_per_kg=0.09):
    """All money inputs are scenarios, USD. No observed plant cost is represented.

    Layer replacement is accrued by consumed log-activity life (salvage retained),
    NOT silently charged a full replacement to a partly used short-window batch.
    """
    voltage = observations(load,kir=2*k/3,kdes=k/3)[:,0]
    electricity = 0.030*np.asarray(load)*voltage # kWh per hour at dt=1 h
    auxiliary = np.full(len(load),0.006)  # synthetic 6 W process auxiliaries
    kg = float(product_kg(load).sum())
    life_fraction = k*float(np.sum(load*load))/-np.log(0.7)
    cost_electric = float(np.sum((electricity+auxiliary)*prices))*electricity_scale
    cost_layer = layer_cost*life_fraction
    # Explicit line items per kg: support/membrane, maintenance, disposal, offspec.
    fixed_parts = {"support_membrane_USD":0.04*kg, "maintenance_USD":0.025*kg,
                   "disposal_USD":0.015*kg,"offspec_reserve_USD":0.010*kg}
    scale = fixed_per_kg/0.09
    fixed_parts = {key:value*scale for key,value in fixed_parts.items()}
    result = {"qualified_Cl2_kg":kg,"charge_Ah":float(sum(load)*30),
              "electricity_kWh":float(sum(electricity)),"auxiliary_kWh":float(sum(auxiliary)),
              "total_kWh_per_tCl2":float(sum(electricity+auxiliary))/kg*1000,
              "electricity_cost_USD":cost_electric,"layer_consumption_USD":cost_layer,
              "life_fraction_consumed":life_fraction,"end_active_fraction":float(np.exp(-k*np.sum(load*load))),
              "predicted_life_h":float(-np.log(0.7)/(k*np.mean(load*load))),
              **fixed_parts}
    result["total_cost_USD_per_tCl2"] = (cost_electric+cost_layer+sum(fixed_parts.values()))/kg*1000
    return result
