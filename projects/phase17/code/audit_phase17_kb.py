"""Small-matrix KB-only correction after the two recorded Phase17 post-audits."""
import argparse
import copy
import importlib.util
import json
import os
from pathlib import Path
import sys
from datetime import datetime, timezone
from audit_phase17_scaled_hydrogenic import sha, preserve, require_completed


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--attempt", type=Path, required=True)
    ap.add_argument("--candidate", type=Path, required=True)
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--expected-source-sha256", required=True)
    args = ap.parse_args()
    attempt = args.attempt.resolve()
    require_completed(json.loads((attempt/"rerun_status.json").read_bytes()))
    data, source = args.candidate.read_bytes(), args.source.read_bytes()
    upstream_bytes = (args.candidate.parent/"audit.json").read_bytes()
    upstream = json.loads(upstream_bytes)
    canonical = attempt/"results_phase17/phase17_results.json"
    original = canonical.read_bytes()
    if (sha(source) != args.expected_source_sha256.lower() or
        upstream.get("status") != "completed" or upstream.get("corrected_json_sha256") != sha(data) or
        upstream.get("original_json_sha256") != sha(original)):
        raise RuntimeError("Source/upstream/canonical provenance mismatch")
    out = attempt/("kb_postaudit_"+datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ"))
    out.mkdir()
    for name, contents in (("verified_source.py", source), ("candidate.prior.json", data),
                           ("canonical.prior.json", original)):
        preserve(out/name, contents)
    audit = dict(status="running", original_json_sha256=sha(original),
                 input_candidate_sha256=sha(data), upstream_audit_sha256=sha(upstream_bytes),
                 source_sha256=sha(source), helper_sha256=sha(Path(__file__).read_bytes()),
                 changed_fields=["certificates.kb_balanced_1s_err", "certificates.kb_refinement"],
                 scope="Gaussian matrix demo and fig4 only; no atomic/N/Breit/Gaunt/QTAIM recomputation",
                 criteria=dict(analytic_relative_error=1e-3, last_basis_relative_change=1e-3,
                               eigen_residual_relative=1e-8),
                 criteria_scope="Declared finite-basis diagnostic checks, not an existing SCF tolerance or proof of no spectral pollution")
    def record():
        (out/"audit.json").write_text(json.dumps(audit, indent=2, allow_nan=False), encoding="utf-8")
    record()
    try:
        for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
            os.environ[key] = "2"
        os.environ["MPLBACKEND"] = "Agg"
        spec = importlib.util.spec_from_file_location("phase17_kb_verified", out/"verified_source.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
        rows = []
        for n in (10, 16, 24, 32):
            _, exact, diag = mod.gaussian_kb_demo(n_basis=n, return_diagnostics=True)
            energy = diag["balanced"]["bound_1s_eh"]
            rows.append(dict(n_basis=n, exact_eh=float(exact), relative_error=abs(energy-exact)/abs(exact), **diag))
        last = rows[-1]
        change = abs(last["balanced"]["bound_1s_eh"]-rows[-2]["balanced"]["bound_1s_eh"])/abs(last["balanced"]["bound_1s_eh"])
        passed = (last["relative_error"] <= 1e-3 and change <= 1e-3 and
                  all(r["balanced"]["eigen_residual_relative"] <= 1e-8 for r in rows))
        audit.update(rows=rows, last_basis_relative_change=change, diagnostic_gate_pass=passed)
        corrected = copy.deepcopy(json.loads(data))
        corrected["certificates"]["kb_balanced_1s_err"] = last["relative_error"]
        corrected["certificates"]["kb_refinement"] = dict(
            rows=rows, last_basis_relative_change=change, diagnostic_gate_pass=passed,
            original_basis_size=10, corrected_default_basis_size=32,
            interpretation="Positive-total-energy bound branch of corrected radial RKB matrix; not global variational stability certification")
        assert {k:v for k,v in corrected.items() if k != "certificates"} == {k:v for k,v in json.loads(data).items() if k != "certificates"}
        result = json.dumps(corrected, indent=2, allow_nan=False).encode()
        (out/"phase17_results.merged_corrected.json").write_bytes(result)
        audit["corrected_json_sha256"] = sha(result)
        if not passed:
            raise RuntimeError("KB diagnostic gate failed; candidate retained, no figure published")
        mod.fig4_validation(corrected, None)
        if canonical.read_bytes() != original:
            raise RuntimeError("Canonical results changed")
        name = "fig4_relativistic_foundation_validation.png"
        target = attempt/"figures_phase17"/name
        old = target.read_bytes()
        preserve(out/("fig4.prior."+sha(old)+".png"), old)
        new = (mod.FIG/name).read_bytes()
        target.write_bytes(new)
        audit.update(status="completed", prior_fig4_sha256=sha(old), fig4_sha256=sha(new))
        record()
        print(f"KB post-audit complete: {out}", flush=True)
    except BaseException as exc:
        audit.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        record()
        raise


if __name__ == "__main__":
    main()
