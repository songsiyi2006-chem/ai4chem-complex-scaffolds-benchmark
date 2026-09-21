"""Three deterministic proposal roles with fixed-kernel GP acquisition.

Agents are auditable software policies, not LLMs or autonomous hardware actors.
Only queried simulator outputs enter the GP; the final Pareto set is observed-only.
"""
from itertools import product
import numpy as np
from scipy.special import ndtr
from common import module

macro = module("03_macro_flow_cfd.py")
tea = module("04_tea_lca_model.py")


def pareto_mask(values):
    a = np.asarray(values, dtype=float)  # all minimization
    if a.ndim != 2 or not np.isfinite(a).all():
        raise ValueError("Finite 2D objective array required")
    return np.array([not np.any(np.all(a <= row, axis=1) & np.any(a < row, axis=1)) for row in a])


def gp_predict(x, y, query):
    # Fixed smoothness/lengthscale and numerical jitter, not a fitted noise estimate.
    scale = max(float(np.std(y)), 1e-12)
    mean = float(np.mean(y))
    y = (y-mean)/scale
    kernel = lambda a,b: np.exp(-.5*np.sum(((a[:,None]-b[None,:])/.5)**2, axis=-1))
    k = kernel(x,x) + np.eye(len(x))*1e-7
    chol = np.linalg.cholesky(k)
    cross = kernel(x,query)
    mu = cross.T @ np.linalg.solve(chol.T, np.linalg.solve(chol,y))
    v = np.linalg.solve(chol,cross)
    sigma = np.sqrt(np.maximum(1-np.sum(v*v,axis=0), 1e-12))
    return mu*scale+mean, sigma*scale


def expected_improvement(best, mu, sigma):
    improvement = best-mu  # minimization
    z = improvement/np.maximum(sigma,1e-15)
    return improvement*ndtr(z)+sigma*np.exp(-z*z/2)/np.sqrt(2*np.pi)


def optimize(config):
    candidates = [dict(material=m, ligand_shift_eV=l, potential_RHE_V=u, velocity_m_s=v)
                  for m,l,u,v in product(("Cu","Ag"),(-.06,0.,.06),(-.5,-.65,-.8,-.95),(.025,.05,.1))]
    x = np.array([[int(c["material"]=="Ag"), c["ligand_shift_eV"],c["potential_RHE_V"],c["velocity_m_s"]]
                  for c in candidates], dtype=float)
    x = (x-x.min(axis=0))/np.ptp(x,axis=0)
    rng = np.random.default_rng(config.seed)
    initial = [int(rng.integers(len(candidates)))]
    while len(initial) < config.initial_samples:
        distances = np.min(np.sum((x[:,None]-x[initial][None,:])**2, axis=-1),axis=1)
        distances[initial] = -1
        initial.append(int(np.argmax(distances)))
    observations, seen, ledger = [], set(), []
    def observe(index, round_id, role, acquisition=None, train_ids=None):
        seen.add(index)
        row = dict(candidate_id=index, round=round_id, role=role, **candidates[index])
        try:
            c = config.with_changes(**candidates[index])
            flow = macro.channel(c)
            result = tea.assess(c, flow)
            row.update(feasible=True, cost_USD_kg=result["LCO_product_USD_kg"],
                       GHG_kgCO2e_kg=result["screening_GHG_kgCO2e_kg"],
                       productivity_kg_h=result["initial_recovered_kg_h"],
                       target_FE=flow["FE"][config.target],
                       carbon_residual_mol_s=flow["carbon_residual_mol_s"])
        except (ValueError, ArithmeticError) as exc:
            row.update(feasible=False, rejection=str(exc))
        observations.append(row)
        ledger.append(dict(candidate_id=index, role=role, round=round_id, acquisition=acquisition,
                           training_candidate_ids=list(train_ids or []), observation_available_after_proposal=True))
    for i in initial:
        observe(i,0,"maxmin_seed")
    for round_id in range(1,config.rounds+1):
        valid = [r for r in observations if r["feasible"]]
        training = [r["candidate_id"] for r in valid]
        batch = []
        for role in ("cost_EI","carbon_EI","uncertainty")[:config.batch_size]:
            available = sorted(set(range(len(candidates)))-seen-set(batch))
            if len(valid) < 2:
                # Explicit fallback: enough valid data do not yet exist for GP.
                idx, acquisition, actual_role = available[0], None, "untrained_exploration"
            else:
                key = "GHG_kgCO2e_kg" if role=="carbon_EI" else "cost_USD_kg"
                # Zero screening emissions is valid for a zero-burden scenario.
                y = np.log1p([r[key] for r in valid])
                mu, sd = gp_predict(x[training],y,x[available])
                score = sd if role=="uncertainty" else expected_improvement(min(y),mu,sd)
                pick = int(np.argmax(score))
                idx, acquisition, actual_role = available[pick], float(score[pick]), role
            batch.append(idx)
            ledger.append(dict(round=round_id,event="proposal",candidate_id=idx,role=actual_role,
                               acquisition=acquisition,training_candidate_ids=training.copy()))
        # Synchronous batch: every proposal precedes any observation this round.
        for idx in batch:
            proposal = next(r for r in reversed(ledger) if r.get("event")=="proposal" and r["candidate_id"]==idx)
            observe(idx,round_id,proposal["role"],proposal["acquisition"],training)
    valid = [r for r in observations if r["feasible"]]
    if not valid:
        raise RuntimeError("No feasible evaluated point; no Pareto front manufactured")
    mask = pareto_mask([[r["cost_USD_kg"],r["GHG_kgCO2e_kg"],-r["productivity_kg_h"]] for r in valid])
    front = [r for r, keep in zip(valid,mask) if keep]
    return dict(algorithm="fixed_kernel_GP_marginal_EI_multi_policy_batch",
                qEI=False, experimental_feedback=False, candidate_pool=len(candidates),
                evaluations=len(observations), observations=observations, proposal_ledger=ledger,
                pareto=front, objectives=["min cost", "min screening GHG", "max target kg/h"],
                uncertainty="GP conditional variance is acquisition-only, not chemistry confidence")
