"""Phase30 CLI: isolated uncalibrated demo, test, preflight, or input audit."""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"code"))
for key in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS"):
    os.environ[key] = "1"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out","--workspace",dest="out",type=Path,help="NEW isolated directory")
    parser.add_argument("--config",type=Path)
    parser.add_argument("--mode",choices=("demo","production"),default="demo")
    actions = parser.add_mutually_exclusive_group()
    actions.add_argument("--self-test",action="store_true")
    actions.add_argument("--preflight",action="store_true")
    actions.add_argument("--prepare-only",action="store_true")
    actions.add_argument("--resume-prepared",action="store_true")
    actions.add_argument("--validate-micro-input",type=Path)
    actions.add_argument("--publish-example",type=Path,help="Explicitly package one completed demo into this project")
    args = parser.parse_args(argv)
    if args.mode == "production":
        parser.error("Production disabled: calibrated DFT/EDL/3D CFD/operando/lifetime evidence gates are not satisfied")
    if args.preflight:
        versions = {}
        for name in ("numpy","scipy","matplotlib"):
            try: versions[name] = importlib.metadata.version(name)
            except importlib.metadata.PackageNotFoundError: versions[name] = None
        print(json.dumps(dict(packages=versions, network_needed=False, demo_only=True,
                              production_available=False, dft_engine_launched=False),indent=2))
        return int(any(v is None for v in versions.values()))
    if args.self_test:
        result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.discover(str(ROOT/"tests")))
        return int(not result.wasSuccessful())
    from common import load_config, validate_micro_artifact
    if args.validate_micro_input:
        print(json.dumps(validate_micro_artifact(args.validate_micro_input),indent=2))
        return 0
    from pipeline import prepare, execute, publish_example
    if args.publish_example:
        publish_example(args.publish_example)
        print("Packaged verified DEMO files only; no Git commit/push performed.")
        return 0
    if args.out is None:
        parser.error("Specify --out NEW_DIRECTORY or an inspection/test action")
    if args.resume_prepared:
        if args.config:
            parser.error("Prepared runs use their frozen config; --config cannot override it")
        out = args.out
    else:
        out = prepare(args.out,load_config(args.config))
    if args.prepare_only:
        print(f"Prepared only: {out}")
        return 0
    result = execute(out)
    print(json.dumps(dict(status=result["status"],evidence="UNCALIBRATED_DEMO",output=str(out)),indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
