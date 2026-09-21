"""Phase31: audit-first local pilot and production-interface specification."""
import argparse
import json
import os
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "code"))
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", "--workspace", type=Path)
    p.add_argument("--mode", choices=("plan", "local-pilot", "production"), default="plan")
    p.add_argument("--token-cap", type=int, default=100_000_000)
    actions = p.add_mutually_exclusive_group()
    actions.add_argument("--self-test", action="store_true")
    actions.add_argument("--preflight", action="store_true")
    actions.add_argument("--fetch-reference", type=Path)
    actions.add_argument("--prepare-only", action="store_true")
    actions.add_argument("--resume", action="store_true")
    actions.add_argument("--worker", action="store_true")
    actions.add_argument("--finalize", action="store_true")
    actions.add_argument("--status", action="store_true")
    args = p.parse_args(argv)
    if args.self_test:
        result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.discover(str(ROOT / "tests")))
        return int(not result.wasSuccessful())
    from swarm.coordinator import preflight, prepare, work, finalize, validate_resume
    if args.preflight:
        print(json.dumps(preflight(), indent=2))
        return 0
    if args.fetch_reference:
        from reference_data import fetch
        print(json.dumps(fetch(args.fetch_reference), indent=2))
        return 0
    if args.mode == "production":
        p.error("Production blocked: unvalidated QM method, QM/MM forces, alchemical topology, sampling and TS/IRC evidence")
    if args.out is None:
        p.error("--out NEW_DIRECTORY required")
    if args.token_cap < 0:
        p.error("nonnegative token cap required")
    if args.resume or args.worker or args.finalize or args.status:
        out, frozen = validate_resume(args.out)
        # Resume always uses frozen mode and cap; override options are forbidden.
        options = [option.split("=", 1)[0] for option in (argv if argv is not None else sys.argv[1:])]
        if "--mode" in options or "--token-cap" in options:
            p.error("resume uses frozen mode and token cap; no overrides")
    else:
        out = prepare(args.out, args.mode, args.token_cap)
    if args.prepare_only:
        print(f"Prepared only: {out}")
        return 0
    if args.status:
        from swarm.memory_ledger import Ledger
        ledger = Ledger(out / "ledger.sqlite", frozen["token_cap"])
        try:
            print(json.dumps(ledger.snapshot(), indent=2))
        finally:
            ledger.close()
        return 0
    if not args.finalize:
        work(out)
    if not args.worker:
        summary = finalize(out)
        print(json.dumps(dict(output=str(out), production=False,
                              provider_tokens=summary["actual_provider_tokens"]), indent=2))
        return int(any(t["status"] == "FAILED" for t in summary["tasks"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
