from __future__ import annotations
import csv
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import least_squares, minimize_scalar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from physics import simulate, observations, candidate_observations, no_batch_leakage, kaplan_meier
from control import schedule, account


def write_csv(path, rows):
    rows = list(rows)
    if not rows: return
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def dump(path,obj):
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")


def validate_supported_config(cfg):
    """Reject edits to fixed benchmark conventions instead of silently ignoring them."""
    fixed={"dt_h":1.0,"clearance_per_h":0.12,"current_density_A_cm2":0.3,
           "geometric_area_cm2":100.0,"faradaic_efficiency":0.985,
           "failure_active_fraction":0.7,"minimum_load":0.75,"maximum_load":1.25,
           "max_ramp_load_per_hour":0.25,"storage_charge_hours":2.0,
           "nominal_layer_replacement_USD":0.2,"kir_per_h":0.003,"kdes_per_h":0.0015}
    for key,value in fixed.items():
        if cfg[key]!=value: raise ValueError(f"Unsupported fixed benchmark convention {key}: expected {value}")
    if cfg["assumed_product_purity"]<cfg["minimum_product_purity"]:
        raise ValueError("The fixed-purity scenario does not meet its quality gate")


def run(root: Path,out: Path):
    cfg=json.loads((root/"configs/local.json").read_text(encoding="utf-8"))
    validate_supported_config(cfg)
    rng=np.random.default_rng(cfg["seed"])
    figures=out/"figures"; figures.mkdir(parents=True)
    sigma=np.array([cfg["sigma_voltage_V"],cfg["sigma_surface_fraction"],cfg["sigma_solution_fraction"]])
    n=cfg["fit_duration_h"]; time=np.arange(n)+1
    batches=[]; rows=[]
    for split,key in [("train","fit_batches"),("validation","validation_batches"),("blind","blind_batches")]:
        for bindex,batch in enumerate(cfg[key]):
            load=1+0.2*np.sin(np.arange(n)/(3.0+bindex))
            if split=="blind":
                load[20:26]=0 # specified synthetic start-stop challenge
                load[40:50]=1.25
            scale=np.exp(rng.normal(0,0.08))
            clean=observations(load,kir=cfg["kir_per_h"]*scale,kdes=cfg["kdes_per_h"]*scale)
            observed=clean+rng.normal(0,sigma,size=clean.shape)
            batches.append(dict(batch_id=batch,split=split,load=load,y=observed,truth_scale=scale))
            for i in range(n):
                rows.append(dict(batch_id=batch,split=split,time_h=int(time[i]),load=float(load[i]),
                    voltage_V=float(observed[i,0]),surface_fraction=float(observed[i,1]),
                    soluble_inventory_fraction=float(observed[i,2]),evidence="SYNTHETIC_SOFTWARE_BENCHMARK"))
    no_batch_leakage(rows)
    write_csv(out/"synthetic_batches.csv",rows)
    dump(out/"synthetic_truth.json",[{"batch_id":b["batch_id"],"rate_multiplier":b["truth_scale"]} for b in batches])
    train=[b for b in batches if b["split"]=="train"]

    def residual(params, selected=train, channels=(0,1,2)):
        return np.concatenate([((observations(b["load"],kir=params[0],kdes=params[1])-b["y"])/sigma)[:,channels].ravel() for b in selected])
    fit=least_squares(residual,[.002,.002],bounds=(1e-8,.03),xtol=1e-11,ftol=1e-11,gtol=1e-11)
    # Independently fit identifiable total loss from voltage, then partition analytically.
    fv=minimize_scalar(lambda k:np.sum(residual([k,0],channels=(0,))**2),bounds=(1e-8,.02),method="bounded",options={"xatol":1e-12})
    total=float(fv.x)
    profile=[]
    for kd in np.linspace(0,total,81):
        voltage_sse=float(np.sum(residual([total-kd,kd],channels=(0,))**2))
        opt=minimize_scalar(lambda ki:np.sum(residual([ki,kd])**2),bounds=(1e-8,.015),method="bounded")
        profile.append(dict(kdes_per_h=float(kd),voltage_chi2=voltage_sse,joint_chi2=float(opt.fun),joint_best_kir_per_h=float(opt.x)))
    write_csv(out/"profile_likelihood.csv",profile)
    # Finite difference Jacobian of both paths at fitted point, standardized observations.
    jacobians={}
    for label,channels in [("voltage_only",(0,)),("joint_observations",(0,1,2))]:
        eps=1e-6; basis=np.eye(2)*eps
        jac=np.column_stack([(residual(fit.x+b,channels=channels)-residual(fit.x-b,channels=channels))/(2*eps) for b in basis])
        sv=np.linalg.svd(jac,compute_uv=False)
        correlation=float(np.dot(jac[:,0],jac[:,1])/(np.linalg.norm(jac[:,0])*np.linalg.norm(jac[:,1])))
        jacobians[label]={"singular_values":sv.tolist(),"rank_relative_1e-8":int(np.sum(sv>sv[0]*1e-8)),"column_correlation":correlation}
    validation=[]
    for b in batches:
        pred=observations(b["load"],kir=fit.x[0],kdes=fit.x[1])
        validation.append(dict(batch_id=b["batch_id"],split=b["split"],voltage_rmse_V=float(np.sqrt(np.mean((pred[:,0]-b["y"][:,0])**2))),
            joint_nrmse=float(np.sqrt(np.mean(((pred-b["y"])/sigma)**2)))))
    write_csv(out/"batch_predictions.csv",validation)
    # All candidates fitted on TRAIN only and scored with the same assumed channel errors.
    comparisons=[]
    for mechanism in ("H1","H2","H3","H4"):
        for mode,channels in [("voltage_only",(0,)),("joint",(0,1,2))]:
            def r(par,selected=train):
                return np.concatenate([((candidate_observations(mechanism,b["load"],par[0],par[1])-b["y"])/sigma)[:,channels].ravel() for b in selected])
            candidate=least_squares(r,[.004,.8],bounds=([1e-6,.001],[.5,2]))
            blind=[b for b in batches if b["split"]=="blind"]
            comparisons.append(dict(mechanism=mechanism,observations=mode,rate_per_h=float(candidate.x[0]),amplitude_V=float(candidate.x[1]),
                train_chi2=float(sum(r(candidate.x)**2)),blind_normalized_rmse=float(np.sqrt(np.mean(r(candidate.x,blind)**2)))))
    # Required empirical/no-degradation baselines, fitted as voltage functions only.
    for name in ("constant_activity","empirical_linear_charge"):
        X=[]; y=[]
        for b in train:
            exposure=np.cumsum(b["load"]**2)
            X.append(np.column_stack([np.ones(n),exposure]) if name=="empirical_linear_charge" else np.ones((n,1)))
            y.append(b["y"][:,0]-3-.12*b["load"])
        coefficients=np.linalg.lstsq(np.vstack(X),np.concatenate(y),rcond=None)[0]
        for b in [x for x in batches if x["split"]=="blind"]:
            bx=np.column_stack([np.ones(n),np.cumsum(b["load"]**2)]) if name=="empirical_linear_charge" else np.ones((n,1))
            pred=3+.12*b["load"]+bx@coefficients
            comparisons.append(dict(mechanism=name,observations="voltage_only:"+b["batch_id"],rate_per_h=0,amplitude_V=float(coefficients[0]),train_chi2=float(np.sum((np.vstack(X)@coefficients-np.concatenate(y))**2)/sigma[0]**2),blind_normalized_rmse=float(np.sqrt(np.mean((pred-b["y"][:,0])**2))/sigma[0])))
    write_csv(out/"mechanism_comparison.csv",comparisons)

    hours=np.arange(cfg["duration_h"])
    prices=.11+.035*np.sin(2*np.pi*(hours-6)/24)+.012*np.sin(hours*1.77)
    availability=np.where((hours%17)<3,1.0,1.25)
    records=[]; controls={}; summaries=[]
    for strategy in ("constant","rule","mpc"):
        load=schedule(strategy,prices,availability); controls[strategy]=load
        summary=account(load,prices); summary["strategy"]=strategy
        summary["mean_load"]=float(np.mean(load)); summary["mean_square_load"]=float(np.mean(load**2))
        summary["terminal_storage_charge_h"]=float(sum(load-1)); summary["max_storage_charge_h"]=float(max(abs(np.cumsum(load-1))))
        summary["max_ramp_load_per_h"]=float(max(abs(np.diff(np.r_[1,load]))))
        temp=25.0; temps=[]
        for x in load: temp=25+(temp-25)*.8+1.0*x*x; temps.append(temp)
        summary["max_temperature_C"]=float(max(temps)); summaries.append(summary)
        if max(temps)>cfg["max_temperature_C"]:
            raise ValueError("Scenario temperature exceeds the frozen admissible limit")
        for i in hours:
            records.append(dict(strategy=strategy,time_h=int(i),price_USD_kWh=float(prices[i]),availability_load=float(availability[i]),load=float(load[i]),storage_charge_h=float(np.sum(load[:i+1]-1)),temperature_C=float(temps[i]),purity_assumed=cfg["assumed_product_purity"],evidence="MODEL_INFERENCE_SYNTHETIC_SCENARIO"))
    write_csv(out/"control_timeseries.csv",records); write_csv(out/"fair_comparison.csv",summaries)
    uncertainty=[]; breakeven=[]
    for draw in range(cfg["uncertainty_draws"]):
        k=.0045*np.exp(rng.normal(0,.30)); layer=rng.uniform(.05,.50); escale=rng.uniform(.65,1.35)
        fixed=rng.uniform(.05,.18); ef=rng.uniform(.15,.65)
        baseline=account(controls["constant"],prices,k,layer,escale,fixed)
        for strategy in ("rule","mpc"):
            evaluated=account(controls[strategy],prices,k,layer,escale,fixed)
            feasible_draw=min(evaluated["end_active_fraction"],baseline["end_active_fraction"])>=cfg["failure_active_fraction"]
            reduction=1-evaluated["total_cost_USD_per_tCl2"]/baseline["total_cost_USD_per_tCl2"]
            uncertainty.append(dict(draw=draw,strategy=strategy,k_total_per_h=float(k),layer_USD=float(layer),electricity_scale=float(escale),fixed_USD_kg=float(fixed),emission_factor_kg_kWh=float(ef),window_feasible=feasible_draw,cost_reduction_fraction=(float(reduction) if feasible_draw else ""),predicted_life_h=evaluated["predicted_life_h"],electricity_only_kgCO2e_per_tCl2=ef*evaluated["total_kWh_per_tCl2"]))
    write_csv(out/"scenario_uncertainty.csv",uncertainty)
    # Break-even replacement value vs assumed decay/life, with all other costs fixed.
    for life in np.linspace(50,250,41):
        k=-np.log(.7)/life
        base=account(controls["constant"],prices,k,0)
        for strategy in ("rule","mpc"):
            candidate=account(controls[strategy],prices,k,0)
            delta_life=candidate["life_fraction_consumed"]-base["life_fraction_consumed"]
            saving=base["electricity_cost_USD"]-candidate["electricity_cost_USD"]
            breakeven.append(dict(strategy=strategy,assumed_constant_life_h=float(life),extra_life_fraction=delta_life,electricity_saving_USD=saving,maximum_layer_replacement_USD=(saving/delta_life if delta_life>1e-12 else "undefined")))
    write_csv(out/"breakeven.csv",breakeven)
    censored_times=np.minimum(np.array([45,68,90,120,180,75]),96)
    events=np.array([45,68,90,120,180,75])<=96
    write_csv(out/"synthetic_survival.csv",[dict(time_h=a,at_risk=b,events=c,censored=d,survival=e) for a,b,c,d,e in kaplan_meier(censored_times,events)])
    # Illustrations explicitly identify synthetic status on every panel.
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for key,label in [("voltage_chi2","Voltage only"),("joint_chi2","Voltage + surface + solution")]:
        values=np.array([p[key] for p in profile]); axes[0].plot([p["kdes_per_h"] for p in profile],values-values.min(),label=label)
    axes[0].axhline(3.84,color="gray",linestyle="--",label="chi-square 3.84 reference")
    axes[0].set(xlabel="Desorption rate / h$^{-1}$",ylabel="Profile delta chi-square",title="Synthetic: loss routes not identified by voltage")
    axes[0].legend(fontsize=7); axes[0].set_ylim(-1,40)
    blind=batches[-1]; prediction=observations(blind["load"],kir=fit.x[0],kdes=fit.x[1])
    axes[1].plot(time,blind["y"][:,0],".",ms=3,label="Synthetic held-out batch")
    axes[1].plot(time,prediction[:,0],label="Training-only physical fit")
    axes[1].set(xlabel="Time / h",ylabel="Cell voltage / V",title="Synthetic blind challenge; not experimental validation")
    axes[1].legend(fontsize=7); fig.tight_layout(); fig.savefig(figures/"identifiability.png",dpi=180); plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    for strategy,load in controls.items(): axes[0].step(hours,load,where="post",label=strategy)
    axes[0].set(xlabel="Time / h",ylabel="Normalized current",title="Synthetic equal-product schedules"); axes[0].legend(fontsize=8)
    for strategy in ("rule","mpc"):
        vals=[100*r["cost_reduction_fraction"] for r in uncertainty if r["strategy"]==strategy and r["window_feasible"]]
        axes[1].hist(vals,bins=35,alpha=.55,label=strategy)
    axes[1].axvline(5,color="black",linestyle="--"); axes[1].set(xlabel="Scenario cost reduction / %",ylabel="Scenario draws",title="Assumption sensitivity; not a posterior"); axes[1].legend(fontsize=8)
    for strategy in ("rule","mpc"):
        subset=[r for r in breakeven if r["strategy"]==strategy]
        axes[2].plot([r["assumed_constant_life_h"] for r in subset],[r["maximum_layer_replacement_USD"] for r in subset],label=strategy)
    axes[2].axhline(.2,color="gray",linestyle="--"); axes[2].set(xlabel="Assumed static lifetime / h",ylabel="Break-even layer price / USD per module",title="Synthetic 100 cm$^2$ module; no quote"); axes[2].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(figures/"control_economics.png",dpi=180); plt.close(fig)
    uq={}
    for strategy in ("rule","mpc"):
        subset=[r for r in uncertainty if r["strategy"]==strategy]
        valid=[r for r in subset if r["window_feasible"]]
        uq[strategy]={"feasible_draws":len(valid),"infeasible_draws":len(subset)-len(valid),"cost_reduction_fraction_q025_q50_q975_conditional_feasible":np.quantile([r["cost_reduction_fraction"] for r in valid],[.025,.5,.975]).tolist(),"fraction_feasible_draws_reaching_5_percent":float(np.mean([r["cost_reduction_fraction"]>=.05 for r in valid])),"lifetime_h_q025_q50_q975_all_draws":np.quantile([r["predicted_life_h"] for r in subset],[.025,.5,.975]).tolist()}
    observed_blind=[r for r in validation if r["split"]=="blind"]
    summary={"evidence":"EXECUTED_CALCULATION: synthetic benchmark and analytic implications only","configuration_sha256":hashlib.sha256((root/"configs/local.json").read_bytes()).hexdigest(),"seed":cfg["seed"],"software":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"matplotlib":matplotlib.__version__},"batch_count":len(batches),"point_count":len(rows),"fit_kir_per_h":float(fit.x[0]),"fit_kdes_per_h":float(fit.x[1]),"voltage_total_loss_per_h":total,"voltage_profile_chi2_range":float(np.ptp([p["voltage_chi2"] for p in profile])),"local_identifiability":jacobians,"synthetic_blind_pass":all(r["voltage_rmse_V"]<cfg["synthetic_blind_voltage_rmse_limit_V"] and r["joint_nrmse"]<cfg["synthetic_blind_joint_channel_nrmse_limit"] for r in observed_blind),"synthetic_blind_predictions":observed_blind,"controls":summaries,"uncertainty":uq,"posterior_predictive_check":"NOT_RUN: no Bayesian posterior; frequentist residual and held-out synthetic checks are provided","statistical_residual_correction":"NOT_ADDED: no external evidence supporting a correction","physical_mass_conservation_max_error":float(np.max(np.abs(simulate(controls["mpc"]).sum(axis=1)-1))),"scientific_mechanism_status":"NOT_IDENTIFIED_IN_REAL_ELECTRODES","industrial_cost_status":"NOT_VALIDATED"}
    dump(out/"summary.json",summary)
    return summary
