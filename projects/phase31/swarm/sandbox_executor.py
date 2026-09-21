"""Trusted child-process isolation. NOT a security sandbox for hostile agent code."""
import json
import os
from pathlib import Path
import subprocess
import sys
from common import ROOT, EvidenceBlocked, load_module, write_json, sha256


class QuantumExecutionSandbox:
    def __init__(self, timeout_seconds=600):
        if not 0 < timeout_seconds <= 3600:
            raise ValueError("bounded timeout required")
        self.timeout = timeout_seconds

    def execute_with_healing(self, spec, out):
        engine = load_module("02_pyscf_qm_engine")
        engine.validate(spec)  # Method/input failures are not recoverable SCF failures.
        import importlib.util
        if importlib.util.find_spec("pyscf") is None:
            raise EvidenceBlocked("Missing PySCF: no attempts launched")
        out = Path(out)
        out.mkdir(parents=True, exist_ok=False)
        schedule = [
            dict(max_cycle=120, damp=0., level_shift=0., ediis=False),
            dict(max_cycle=180, damp=0.2, level_shift=0., ediis=True),
            dict(max_cycle=240, damp=0.3, level_shift=0.1, ediis=True),
            dict(max_cycle=300, damp=0.1, level_shift=0.2, ediis=False),
            dict(max_cycle=360, damp=0., level_shift=0.05, ediis=False),
        ]
        worker = ROOT / "code" / "02_pyscf_qm_engine.py"
        invariant = engine.fingerprint(spec)
        for index, numerical in enumerate(schedule):
            attempt = out / f"attempt_{index + 1:02}"
            attempt.mkdir()
            request = attempt / "request.json"
            write_json(request, dict(spec=spec, numerical=numerical))
            env = {k: v for k, v in os.environ.items() if k.upper() in {
                "PATH", "SYSTEMROOT", "WINDIR", "TEMP", "TMP", "LANG", "LD_LIBRARY_PATH"}}
            env.update(PYTHONPATH=str(ROOT) + os.pathsep + str(ROOT / "code"),
                       OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       PYTHONIOENCODING="utf-8")
            timed_out = False
            with (attempt / "stdout.log").open("wb") as stdout, (attempt / "stderr.log").open("wb") as stderr:
                try:
                    completed = subprocess.run(
                        [sys.executable, str(worker), "--spec", str(request), "--out", str(attempt)],
                        cwd=attempt, env=env, shell=False, stdout=stdout, stderr=stderr, timeout=self.timeout)
                    exit_code = completed.returncode
                except subprocess.TimeoutExpired:
                    timed_out, exit_code = True, None
            result = None
            if exit_code == 0 and (attempt / "result.json").is_file():
                result = json.loads((attempt / "result.json").read_text(encoding="utf-8"))
            accepted = bool(result and result.get("converged") is True and
                            result.get("method_fingerprint") == invariant)
            record = dict(attempt=index + 1, worker_sha256=sha256(worker), input_sha256=sha256(request),
                          numerical=numerical, method_fingerprint=invariant, exit_code=exit_code,
                          timed_out=timed_out, accepted=accepted)
            with (out / "error_healing_log.jsonl").open("a", encoding="utf-8") as stream:
                stream.write(json.dumps(record, allow_nan=False) + "\n")
            if accepted:
                import math
                if not math.isfinite(result.get("energy_hartree", float("nan"))):
                    raise ValueError("nonfinite energy despite convergence flag")
                return result
            # Only an explicit unconverged SCF result qualifies for numeric retry.
            if timed_out or exit_code != 0 or not result or result.get("converged") is not False:
                raise EvidenceBlocked("Non-SCF failure needs review; logs retained, no automatic code rewriting")
        raise EvidenceBlocked("Five numerical attempts exhausted; molecule is NOT blacklisted as chemically impossible")
