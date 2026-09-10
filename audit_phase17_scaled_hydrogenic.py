"""Targeted, post-completion hydrogenic audit; never runs atomic or Gaunt stages."""
import argparse
import ast
import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import stat
import sys
from datetime import datetime, timezone


def sha(data):
    return hashlib.sha256(data).hexdigest()


def require_completed(status):
    if status.get("status") != "process_completed":
        raise RuntimeError("Main attempt has not completed; do not run concurrently")
    stages = {s["stage"]: s.get("exit_code") for s in status.get("stages", [])}
    if any(stages.get(s) != 0 for s in ("breit", "multiplet", "figures")):
        raise RuntimeError("All three main stages must have recorded exit code zero")


def preserve(path, data):
    """Exclusive content-addressed prior: verify, never replace, mark read-only."""
    if path.exists():
        if path.read_bytes() != data:
            raise RuntimeError(f"Immutable prior mismatch: {path}")
    else:
        with path.open("xb") as stream:
            stream.write(data)
    path.chmod(stat.S_IREAD)


def derive_coulomb(rows):
    one = [r for r in rows if r["Z"] == 1]
    if len(one) != 1 or not math.isfinite(one[0]["J_quadrature"]):
        raise ValueError("Exactly one finite Z=1 anchor required")
    anchor = one[0]["J_quadrature"]
    result = []
    for row in rows:
        z = row["Z"]
        if z <= 0:
            raise ValueError("Positive Z required")
        value, exact = z * anchor, 5 * z / 8
        result.append(dict(Z=z, J_quadrature=value, J_exact=exact,
                           rel_err=abs(value-exact)/exact,
                           grid="uniform Z*r", n_radial=4000,
                           provenance="derived from this attempt Z=1, not recomputed"))
    return result


def verify_anchor_source(old_bytes, new_bytes):
    names = {"breit_pair_energy", "breit_pair_tabulation", "angular_pair_contraction"}
    def definitions(data):
        return {node.name: ast.dump(node, include_attributes=False)
                for node in ast.parse(data).body
                if isinstance(node, ast.FunctionDef) and node.name in names}
    old, new = definitions(old_bytes), definitions(new_bytes)
    if set(old) != names or old != new:
        raise RuntimeError("Coulomb anchor requires identical angular and radial implementations")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--attempt", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--radial-rtol", type=float, default=1e-3,
                        help="3000 versus 6000 node relative-change gate; not angular accuracy")
    args = parser.parse_args()
    if not 0 < args.radial_rtol < 1:
        parser.error("radial-rtol must lie between zero and one")
    attempt = args.attempt.resolve()
    status_path = attempt / "rerun_status.json"
    status_bytes = status_path.read_bytes()
    status = json.loads(status_bytes)
    require_completed(status)
    source_bytes = args.source.read_bytes()
    if sha(source_bytes) != args.expected_source_sha256.lower():
        raise RuntimeError("Exact source SHA256 mismatch")
    prior_path = attempt / "results_phase17" / "phase17_results.json"
    prior_bytes = prior_path.read_bytes()
    prior = json.loads(prior_bytes)
    for key in ("validation", "breit", "bonding"):
        if key not in prior:
            raise RuntimeError(f"Missing completed result section: {key}")
    # Z=1 must come from the verified source snapshot with the corrected r</r> rule.
    old_source = attempt / args.source.name
    old_bytes = old_source.read_bytes()
    if sha(old_bytes) != status["source_sha256"][args.source.name]:
        raise RuntimeError("Attempt source differs from its recorded SHA")
    verify_anchor_source(old_bytes, source_bytes)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out = attempt / ("hydrogenic_postaudit_" + stamp)
    out.mkdir()  # New directory guarantees a fresh angular cache, not historical reuse.
    preserve(out / f"phase17_results.prior.{sha(prior_bytes)}.json", prior_bytes)
    preserve(out / "verified_source.py", source_bytes)
    audit = dict(status="running", source_sha256=sha(source_bytes),
                 helper_sha256=sha(Path(__file__).read_bytes()),
                 original_json_sha256=sha(prior_bytes),
                 original_source_sha256=sha(old_bytes),
                 main_status_sha256=sha(status_bytes), attempt=str(attempt),
                 prerequisite_artifacts=status.get("prerequisite_artifacts"),
                 recomputed=["full-default Breit angular table", "scaled hydrogenic Breit"],
                 not_recomputed=["atomic", "Gaunt", "Coulomb angular", "multiplet", "bonding"],
                 radial_rtol=args.radial_rtol,
                 angular_convergence="not established by this radial-only audit")
    def record():
        (out / "audit.json").write_text(json.dumps(audit, indent=2, allow_nan=False), encoding="utf-8")
    record()
    try:
        for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
            os.environ[key] = "2"
        os.environ["MPLBACKEND"] = "Agg"
        spec = importlib.util.spec_from_file_location("phase17_verified_audit", out / "verified_source.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
        import numpy as np
        if (mod.ANG_TH, mod.ANG_PH, mod.U_GRID) != (30, 36, 40):
            raise RuntimeError("Expected unchanged full angular defaults 30/36/40")
        tab = mod.cached_breit_pair_tabulation(-1, kernel="breit")
        audit["angular_grid"] = dict(theta=30, phi=36, radius_ratios=40)
        audit["cache_sha256"] = {p.name: sha(p.read_bytes()) for p in mod.RES.glob("*.npz")}
        corrected = copy.deepcopy(prior)
        br = corrected["breit"]
        br["coulomb_certification"] = derive_coulomb(prior["breit"]["coulomb_certification"])
        audit["coulomb_linearity"] = (
            "r_Z=r_1/Z, w_Z=w_1/Z, normalized P_Z=sqrt(Z)*P_1; "
            "r</r> and angular table unchanged, 1/r> scales by Z. Therefore "
            "the identical discrete quadrature gives J_Z=Z*J_1. Uses 4000 "
            "nodes with Z*r in [1e-5,40]. Inherits measured Z=1 angular error; "
            "not an independent high-Z certification.")
        audit["coulomb_rows"] = br["coulomb_certification"]
        diagnostics = []
        for z in br["hydrogenic_Z"]:
            energies = []
            for nodes in (3000, 6000):
                x, w = mod.hydrogenic_radial_grid(z, nodes)
                gamma = math.sqrt(1 - (z / mod.C_LIGHT)**2)
                p = x**gamma * np.exp(-z*x)
                q = -p * mod.C_LIGHT * (1-gamma) / z
                norm = math.sqrt(float(np.sum((p*p+q*q)*w)))
                p, q = p/norm, q/norm
                energy = complex(mod.breit_pair_energy(tab, p, q, x, w))
                if not (math.isfinite(energy.real) and math.isfinite(energy.imag)):
                    raise RuntimeError("Nonfinite hydrogenic Breit energy")
                if abs(energy.imag) > 1e-10 * max(abs(energy.real), 1e-15):
                    raise RuntimeError("Nonnegligible imaginary Breit energy")
                energies.append(energy.real)
            change = abs(energies[0]-energies[1])/max(abs(energies[1]), 1e-30)
            diagnostics.append(dict(Z=z, signed_energy_eh_3000=energies[0],
                                    signed_energy_eh_6000=energies[1],
                                    radial_relative_change=change,
                                    radial_gate_pass=change <= args.radial_rtol))
            audit["hydrogenic_rows"] = diagnostics
            record()
            print(f"Z={z}: E3000={energies[0]:.12g}, E6000={energies[1]:.12g}, relative change={change:.6g}", flush=True)
        br["hydrogenic_breit_cm"] = [abs(r["signed_energy_eh_3000"])*mod.HARTREE_CM for r in diagnostics]
        br["hydrogenic_grid"] = dict(coordinate="uniform Z*r", n_radial=3000,
                                     refinement_nodes=6000, audit_record=str(out / "audit.json"))
        # Enforce the deliberately narrow patch: all other sections/fields stay identical.
        allowed = {"coulomb_certification", "hydrogenic_breit_cm", "hydrogenic_grid"}
        assert {k:v for k,v in br.items() if k not in allowed} == {k:v for k,v in prior["breit"].items() if k not in allowed}
        assert {k:v for k,v in corrected.items() if k != "breit"} == {k:v for k,v in prior.items() if k != "breit"}
        candidate = json.dumps(corrected, indent=2, allow_nan=False).encode()
        (out / "phase17_results.corrected.json").write_bytes(candidate)
        audit["corrected_json_sha256"] = sha(candidate)
        if not all(r["radial_gate_pass"] for r in diagnostics):
            raise RuntimeError("Measured radial refinement gate failed; candidate retained, figure not published")
        require_completed(json.loads(status_path.read_bytes()))
        if prior_path.read_bytes() != prior_bytes:
            raise RuntimeError("Main results changed during audit; refusing figure publication")
        mod.fig4_validation(corrected, None)  # Only figure 4; module output is isolated.
        name = "fig4_relativistic_foundation_validation.png"
        fig_bytes = (mod.FIG / name).read_bytes()
        target = attempt / "figures_phase17" / name
        if target.exists():
            old_fig = target.read_bytes()
            preserve(out / ("prior." + sha(old_fig) + ".png"), old_fig)
            audit["prior_fig4_sha256"] = sha(old_fig)
        target.write_bytes(fig_bytes)
        audit.update(status="completed", fig4_sha256=sha(fig_bytes),
                     acceptance="radial refinement gate passed; overall scientific acceptance remains pending")
        record()
        print(f"Postaudit complete: {out}", flush=True)
    except BaseException as exc:
        audit.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        record()
        raise


if __name__ == "__main__":
    main()
