"""Assemble isolated evidence, figures and machine-readable acceptance."""
import csv
from dataclasses import asdict
import importlib.metadata
import platform
from pathlib import Path
import shutil
import sys
import time
import traceback
import numpy as np
from common import Config, dump, sha, module, PRODUCTS, LABEL

MICRO = module("01_micro_surface.py")
MESO = module("02_meso_microkinetics.py")
MACRO = module("03_macro_flow_cfd.py")
TEA = module("04_tea_lca_model.py")
AGENT = module("05_agent_closed_loop.py")
PROJECT = Path(__file__).resolve().parents[1]


def source_hashes():
    return {str(p.relative_to(PROJECT)).replace("\\","/"): sha(p)
            for sub in ("code","configs","tests") for p in sorted((PROJECT/sub).rglob("*"))
            if p.is_file() and p.suffix in (".py",".json")} | {"run.py":sha(PROJECT/"run.py")}


def csv_write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w",encoding="utf-8",newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def prepare(out, config):
    out = Path(out).resolve()
    repo = PROJECT.parents[1]
    if out == repo or out in repo.parents or (out.is_relative_to(repo) and not out.is_relative_to(repo/"work")):
        raise ValueError("Output must be a NEW directory under repo/work or outside the repository")
    out.mkdir(parents=True,exist_ok=False)
    hashes = source_hashes()
    for name in hashes:
        dest = out/"source"/name
        dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(PROJECT/name,dest)
    dump(out/"resolved_config.json",asdict(config))
    dump(out/"run_status.json",dict(status="PREPARED",source_sha256=hashes,
                                   config_sha256=sha(out/"resolved_config.json"),evidence_type=LABEL))
    return out


def make_figures(out, config, map_rows, flow, engineering, sensitivity):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,
                         "axes.spines.right":False,"savefig.dpi":300})
    folder = out/"figures"
    folder.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(12,3.6))
    ax.set(xlim=(0,10),ylim=(0,3)); ax.axis("off")
    for x, text in zip((.95,3.6,6.25,8.95),
                       ("MICRO\nBarrier contract\nDFT/AIMD pending",
                        "MESO\nCoverage ODE + EDL\nUncalibrated closure",
                        "MACRO\n1-D transport surrogate\n3-D CFD pending",
                        "ENGINEERING\n100 kW x 10 scenario\nTEA / screening GHG")):
        ax.text(x,1.75,text,ha="center",va="center",bbox=dict(boxstyle="square,pad=.6",fc="white",ec="#444444"))
    for x in (2.05,4.7,7.35):
        ax.annotate("",xy=(x+.47,1.75),xytext=(x,1.75),arrowprops=dict(arrowstyle="->",color="#315b80"))
    ax.annotate("",xy=(.95,.55),xytext=(8.95,.55),arrowprops=dict(arrowstyle="->",color="#315b80"))
    ax.text(5,.15,"Cost-EI / carbon-EI / exploration -> feasibility gate -> queried model only",
            ha="center")
    ax.set_title("Phase30 coupling architecture | software prototype; no industrial validation",pad=15)
    fig.tight_layout(); fig.savefig(folder/"multiscale_architecture.png"); plt.close(fig)
    us = sorted(set(r["potential_RHE_V"] for r in map_rows))
    phs = sorted(set(r["bulk_pH"] for r in map_rows))
    fig, axes = plt.subplots(1,3,figsize=(12,4.2),layout="constrained",sharey=True)
    cmap = LinearSegmentedColormap.from_list("phase30_blue",["#f3f5f7","#315b80"])
    cmap.set_bad("#cccccc")
    for ax,key in zip(axes,("C2H4","CO","HCOOH")):
        data = np.full((len(phs),len(us)),np.nan)
        for row in map_rows:
            if row["feasible"]:
                data[phs.index(row["bulk_pH"]),us.index(row["potential_RHE_V"])] = 100*row["FE_"+key]
        img = ax.pcolormesh(us,phs,data,vmin=0,vmax=100,cmap=cmap,shading="nearest")
        ax.set(title=key+" Faradaic efficiency",xlabel="Cathode potential (V vs bulk RHE)")
    axes[0].set_ylabel("Bulk pH")
    fig.colorbar(img,ax=axes,label="FE (%); denominator includes H2",shrink=.8)
    fig.suptitle("Uncalibrated Cu rate map | grey = outside supported pH envelope",fontsize=12)
    fig.savefig(folder/"microkinetics_rate_map.png"); plt.close(fig)
    fig, axes = plt.subplots(1,2,figsize=(11,4),layout="constrained")
    trace = engineering["degradation_trace"]
    axes[0].plot([0]+[r["scheduled_h"]/1000 for r in trace],[1]+[r["activity"] for r in trace],
                 color="#315b80",lw=1.5)
    axes[0].axhline(config.replacement_threshold,color="#555555",ls="--",label="replacement threshold")
    axes[0].set(xlabel="Scheduled operating window (10³ h)",ylabel="Assumed target activity",ylim=(0,1.05))
    axes[0].legend(frameon=False)
    axes[1].plot([r["degradation_per_operating_h"]*1e5 for r in sensitivity],
                 [r["LCO_product_USD_kg"] for r in sensitivity],marker="o",color="#315b80")
    axes[1].set(xlabel="Assumed decay constant (10⁻⁵ / operating h)",ylabel=f"Levelized {config.target} cost (USD/kg)")
    axes[1].set_ylim(bottom=0)
    fig.suptitle("Unvalidated lifetime scenarios | fixed-load operation; no confidence interval",fontsize=12)
    fig.savefig(folder/"stack_degradation_tea.png"); plt.close(fig)


def execute(out):
    import json
    out = Path(out).resolve()
    status = json.loads((out/"run_status.json").read_text(encoding="utf-8"))
    if status["status"] != "PREPARED":
        raise ValueError("Only an intact PREPARED directory can execute; no overwriting completed/failed runs")
    if status["source_sha256"] != source_hashes() or status["config_sha256"] != sha(out/"resolved_config.json"):
        raise ValueError("Prepared source/config changed")
    for name, checksum in status["source_sha256"].items():
        if sha(out/"source"/name) != checksum:
            raise ValueError("Prepared source snapshot changed")
    config = Config(**json.loads((out/"resolved_config.json").read_text(encoding="utf-8"))).validate()
    # Exclusive creation prevents two callers executing the same prepared workspace.
    with (out/"execution.lock").open("x",encoding="utf-8") as lock:
        lock.write("No automatic stale-lock recovery; use a new workspace.\n")
    started = time.monotonic()
    status.update(status="RUNNING",started_UTC=time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),
                  python=platform.python_version(), platform=platform.system(),
                  versions={n:importlib.metadata.version(n) for n in ("numpy","scipy","matplotlib")})
    dump(out/"run_status.json",status)
    try:
        raw = MICRO.packet(config)
        dump(out/"results/raw_intermediates.json",raw)
        print("01 micro: synthetic parameter packet; DFT/AIMD not executed",flush=True)
        meso = MESO.solve(config,config.co2_inlet_mol_m3,validate_ode=True)
        dump(out/"results/meso_reference.json",meso)
        flow = MACRO.channel(config)
        fine = MACRO.channel(config.with_changes(cells=config.cells*2))
        relative_grid_change = abs(flow["current_A"]-fine["current_A"])/fine["current_A"]
        if relative_grid_change > .1:
            raise ArithmeticError("Pilot mesh-refinement guard failed (>10% current change)")
        dump(out/"results/flow_reference.json",flow)
        csv_write(out/"results/channel_profile.csv",flow["cells"])
        print(f"02-03 balances: carbon={flow['carbon_residual_mol_s']:.3e} mol/s; mesh delta={relative_grid_change:.3g}",flush=True)
        single = TEA.assess(config,flow,stacks=1)
        engineering = TEA.assess(config,flow)
        dump(out/"reports/tea_lca_breakdown.json",dict(single_stack=single,stack_array=engineering))
        sensitivity = [dict(degradation_per_operating_h=k,
                            **{name:res[name] for name in ("LCO_product_USD_kg","total_replacements")})
                       for k in (0.,1e-5,2e-5,4e-5)
                       for res in [TEA.assess(config.with_changes(degradation_per_operating_h=k),flow)]]
        csv_write(out/"results/degradation_sensitivity.csv",sensitivity)
        optimization = AGENT.optimize(config)
        dump(out/"results/closed_loop.json",optimization)
        csv_write(out/"results/pareto_optimization.csv",optimization["pareto"])
        print(f"04-05 closed loop: {optimization['evaluations']} queried; {len(optimization['pareto'])} observed Pareto points",flush=True)
        map_rows = []
        for pH in np.linspace(5,9,9):
            for u in np.linspace(-.95,-.4,12):
                row = dict(potential_RHE_V=float(u),bulk_pH=float(pH),feasible=False,
                           **{"FE_"+p:None for p in PRODUCTS})
                try:
                    s = MESO.solve(config.with_changes(material="Cu",bulk_pH=float(pH),potential_RHE_V=float(u)),33.)
                    row.update(feasible=True,**{"FE_"+p:s["FE"][p] for p in PRODUCTS})
                except ValueError:
                    pass
                map_rows.append(row)
        csv_write(out/"results/rate_map.csv",map_rows)
        make_figures(out,config,map_rows,flow,engineering,sensitivity)
        acceptance = dict(software_pipeline_executed=True, evidence_type=LABEL,
                          science_status="NOT_SCIENTIFICALLY_VALIDATED",
                          dft_aimd=False,cfd_3d=False,operando_data=False,industrial_validation=False,
                          checks=dict(carbon_balance_mol_s=flow["carbon_residual_mol_s"],
                                      heat_balance_W=flow["heat_residual_W"],
                                      FE_sum=sum(flow["FE"].values()),
                                      ode_analytic_error=meso["ode_validation"]["analytic_error"],
                                      doubled_mesh_current_relative_change=relative_grid_change),
                          resource_scope="CPU-only bounded pilot; no external jobs or network calls",
                          missing=["constant-potential barriers and raw convergence", "calibrated EDL/speciation",
                                   "3D two-phase mesh validation", "measured FE/polarization/degradation",
                                   "complete life-cycle inventory", "prospective blind optimization comparison"])
        dump(out/"results/acceptance.json",acceptance)
        (out/"reports/demo_execution.md").write_text(
            "# Phase30 executed model scenario\n\nNOT SCIENTIFICALLY VALIDATED.\n\n"
            f"Target: {config.target}; installed array: {engineering['installed_power_kW']} kW.\n"
            f"Conditional model cost: {engineering['LCO_product_USD_kg']:.4f} USD/kg.\n"
            f"Conditional screening GHG: {engineering['screening_GHG_kgCO2e_kg']:.4f} kgCO2e/kg.\n\n"
            "These are uncalibrated assumptions, not DFT predictions, market estimates or measured hardware results.\n",
            encoding="utf-8")
        if source_hashes() != status["source_sha256"]:
            raise RuntimeError("Sources changed during execution; evidence not accepted")
        status.update(status="COMPLETED_DEMO",scientific_acceptance=False,elapsed_s=time.monotonic()-started)
    except BaseException as exc:
        status.update(status="FAILED",scientific_acceptance=False,error=str(exc),elapsed_s=time.monotonic()-started)
        (out/"failure.txt").write_text(traceback.format_exc(),encoding="utf-8")
        raise
    finally:
        dump(out/"run_status.json",status)
        files = {str(p.relative_to(out)).replace("\\","/"):dict(sha256=sha(p),bytes=p.stat().st_size)
                 for p in sorted(out.rglob("*")) if p.is_file() and p.name not in ("artifact_manifest.json",)}
        dump(out/"artifact_manifest.json",dict(schema_version="phase30.audit.v1",files=files,
                                              scope="all files in isolated run; manifest excludes itself"))
    return status


def publish_example(run):
    """Explicit local packaging of one hash-verified demo; never automatic Git writes."""
    import json
    run = Path(run).resolve()
    status = json.loads((run/"run_status.json").read_text(encoding="utf-8"))
    if status["status"] != "COMPLETED_DEMO" or source_hashes() != status["source_sha256"]:
        raise ValueError("Require a completed demo with current source hashes")
    manifest = json.loads((run/"artifact_manifest.json").read_text(encoding="utf-8"))
    for name, record in manifest["files"].items():
        src = (run/name).resolve()
        if not src.is_relative_to(run) or sha(src) != record["sha256"]:
            raise ValueError("Run artifact missing/changed or unsafe manifest")
    selected = [p for sub in ("results","reports","figures") for p in (run/sub).glob("*") if p.is_file()]
    if any((PROJECT/p.relative_to(run)).exists() for p in selected):
        raise ValueError("Published example paths already exist; do not overwrite evidence")
    for src in selected:
        dst = PROJECT/src.relative_to(run)
        dst.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(src,dst)
    dump(PROJECT/"results/example_provenance.json",
         dict(status=status, resolved_config=json.loads((run/"resolved_config.json").read_text(encoding="utf-8")),
              selected_files={str(p.relative_to(run)).replace("\\","/"):sha(p) for p in selected},
              excluded="Full input source snapshots, lock and runtime manifest remain in the isolated run"))
