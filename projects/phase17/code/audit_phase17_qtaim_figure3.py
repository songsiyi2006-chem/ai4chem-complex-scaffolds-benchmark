"""Bounded cached-density QTAIM-proxy/fig3 correction after the Breit audit."""
import argparse
import copy
import importlib.util
import json
import os
from pathlib import Path
import sys
from datetime import datetime, timezone
from audit_phase17_scaled_hydrogenic import sha, preserve, require_completed


def validate_merge(prior, candidate):
    """Accept only the explicit Breit-only delta, never silently merge other work."""
    if {k:v for k,v in prior.items() if k != "breit"} != {k:v for k,v in candidate.items() if k != "breit"}:
        raise RuntimeError("Input candidate contains changes outside Breit")
    allowed = {"coulomb_certification", "hydrogenic_grid", "hydrogenic_breit_cm"}
    if {k:v for k,v in prior["breit"].items() if k not in allowed} != {k:v for k,v in candidate["breit"].items() if k not in allowed}:
        raise RuntimeError("Input candidate contains unexpected Breit changes")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--attempt", required=True, type=Path)
    ap.add_argument("--source", required=True, type=Path)
    ap.add_argument("--expected-source-sha256", required=True)
    ap.add_argument("--breit-candidate", required=True, type=Path)
    args = ap.parse_args()
    attempt = args.attempt.resolve()
    status_path = attempt/"rerun_status.json"
    require_completed(json.loads(status_path.read_bytes()))
    source = args.source.read_bytes()
    if sha(source) != args.expected_source_sha256.lower():
        raise RuntimeError("Exact source SHA mismatch")
    prior_path = attempt/"results_phase17"/"phase17_results.json"
    prior_bytes = prior_path.read_bytes()
    prior = json.loads(prior_bytes)
    candidate_bytes = args.breit_candidate.read_bytes()
    candidate = json.loads(candidate_bytes)
    breit_audit_path = args.breit_candidate.parent/"audit.json"
    breit_audit_bytes = breit_audit_path.read_bytes()
    breit_audit = json.loads(breit_audit_bytes)
    if (breit_audit.get("status") != "completed" or
        breit_audit.get("original_json_sha256") != sha(prior_bytes) or
        breit_audit.get("corrected_json_sha256") != sha(candidate_bytes)):
        raise RuntimeError("Breit audit incomplete or candidate/prior hash mismatch")
    validate_merge(prior, candidate)
    plane_path = attempt/"results_phase17"/"phase17_qtaim_planes.npz"
    plane_bytes = plane_path.read_bytes()
    out = attempt/("qtaim_postaudit_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ"))
    out.mkdir()
    for name, data in (("verified_source.py", source),
                       ("original_results.prior.json", prior_bytes),
                       ("breit_candidate.prior.json", candidate_bytes),
                       ("planes.prior.npz", plane_bytes)):
        preserve(out/name, data)
    audit = dict(status="running", source_sha256=sha(source),
                 helper_sha256=sha(Path(__file__).read_bytes()),
                 utility_helper_sha256=sha(Path(__file__).with_name("audit_phase17_scaled_hydrogenic.py").read_bytes()),
                 original_json_sha256=sha(prior_bytes),
                 input_breit_candidate_sha256=sha(candidate_bytes),
                 input_breit_audit_sha256=sha(breit_audit_bytes),
                 input_planes_sha256=sha(plane_bytes),
                 changed_fields=["bonding.qtaim.Am", "bonding.qtaim.Eu", "bonding.model_scope"],
                 inherited_fields="Breit from hash-verified completed audit; all other fields unchanged",
                 computations="Cached fresh promolecular plane differentiation and fig3 only; no SCF/FFT/overlap/angular stages",
                 limitations=["N SCF metadata not persisted in original cache; convergence cannot be independently certified here",
                              "Linear radial interpolation and display-grid derivatives lack refinement certification",
                              "Coarsening comparison is diagnostic, not a fine-grid convergence pass",
                              "Promolecule minimum not stationary-point/Hessian-certified molecular QTAIM",
                              "WH mixing and density-weighted S_rms are proxies, not molecular EDA or signed overlap"])
    def record():
        (out/"audit.json").write_text(json.dumps(audit, indent=2, allow_nan=False), encoding="utf-8")
    record()
    try:
        for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
            os.environ[key] = "2"
        os.environ["MPLBACKEND"] = "Agg"
        spec = importlib.util.spec_from_file_location("phase17_qtaim_verified", out/"verified_source.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
        import numpy as np
        corrected = copy.deepcopy(candidate)
        corrected["bonding"]["model_scope"] = "; ".join(audit["limitations"])
        audit["ions"] = {}
        with np.load(out/"planes.prior.npz", allow_pickle=False) as cache:
            for ion in ("Am", "Eu"):
                xs, zs = cache[ion+"_xs"], cache[ion+"_zs"]
                rho = cache[ion+"_rho2d"]
                R = float(cache[ion+"_R"])
                if rho.shape != (len(xs), len(zs)) or not np.isfinite(rho).all() or np.any(rho <= 0):
                    raise RuntimeError("Invalid cached plane density")
                fine = mod.promolecular_axis_descriptor(rho, xs, zs, R)
                coarse = mod.promolecular_axis_descriptor(rho[::2, ::2], xs[::2], zs[::2], R)
                lap = mod._laplacian_plane(rho, xs[1]-xs[0], zs[1]-zs[0])
                if not np.isfinite(lap[1:-1, :-1]).all():
                    raise RuntimeError("Nonfinite interior corrected Laplacian")
                oldlap = cache[ion+"_lap2d"]
                minimum = np.unravel_index(np.argmin(oldlap), oldlap.shape)
                audit["ions"][ion] = dict(old=prior["bonding"]["qtaim"][ion], corrected=fine,
                    coarsened=coarse, lap_coarsening_absolute_change=abs(fine["lap_bcp"]-coarse["lap_bcp"]),
                    old_minimum_actual_x=float(xs[minimum[0]]), old_minimum_actual_z=float(zs[minimum[1]]),
                    corrected_lap_min=float(np.nanmin(lap)), corrected_lap_max=float(np.nanmax(lap)),
                    undefined_boundary_cells=int(np.isnan(lap).sum()),
                    core_display_mask_radius=float(3*max(xs[1]-xs[0], zs[1]-zs[0])))
                corrected["bonding"]["qtaim"][ion].update(fine)
                mod.PLANE_CACHE[ion] = dict(xs=xs, zs=zs, rho2d=rho, lap2d=lap, R=R, x_bcp=fine["x_bcp_a0"])
        # Everything except the enumerated bonding descriptors must survive verbatim.
        assert {k:v for k,v in corrected.items() if k != "bonding"} == {k:v for k,v in candidate.items() if k != "bonding"}
        assert {k:v for k,v in corrected["bonding"].items() if k not in ("qtaim", "model_scope")} == {k:v for k,v in candidate["bonding"].items() if k not in ("qtaim", "model_scope")}
        data = json.dumps(corrected, indent=2, allow_nan=False).encode()
        (out/"phase17_results.merged_corrected.json").write_bytes(data)
        np.savez_compressed(out/"phase17_qtaim_planes.corrected.npz", **{
            f"{ion}_{key}": value for ion, plane in mod.PLANE_CACHE.items() for key, value in plane.items()})
        mod.fig3_qtaim(corrected)
        require_completed(json.loads(status_path.read_bytes()))
        if prior_path.read_bytes() != prior_bytes:
            raise RuntimeError("Canonical JSON changed during audit")
        name = "fig3_qtaim_relativistic_covalency_map.png"
        target = attempt/"figures_phase17"/name
        oldfig = target.read_bytes()
        preserve(out/("fig3.prior."+sha(oldfig)+".png"), oldfig)
        newfig = (mod.FIG/name).read_bytes()
        target.write_bytes(newfig)
        audit.update(status="completed", corrected_json_sha256=sha(data),
                     prior_fig3_sha256=sha(oldfig), corrected_fig3_sha256=sha(newfig),
                     acceptance="Coordinate/operator corrections applied; numerical derivative convergence and molecular interpretation NOT certified")
        record()
        print(f"QTAIM-proxy audit complete: {out}", flush=True)
    except BaseException as exc:
        audit.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        record()
        raise


if __name__ == "__main__":
    main()
