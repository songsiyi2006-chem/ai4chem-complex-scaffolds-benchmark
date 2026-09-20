"""Run the bounded, isolated Phase28 analysis without modifying older phases."""
from pathlib import Path
import argparse
import json
import os
import sys
import unittest

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"code"))
for key in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS"):
    os.environ[key]="2"

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,help="New output directory; existing directories are refused")
    parser.add_argument("--self-test",action="store_true",help="Run conservation, causality, unit and censoring tests")
    args=parser.parse_args()
    if args.self_test:
        result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.discover(str(ROOT/"tests")))
        return 0 if result.wasSuccessful() else 1
    if args.out is None: parser.error("Specify --out NEW_DIRECTORY or --self-test")
    out=args.out.resolve()
    if out.exists(): parser.error("Output directory already exists; choose a new directory")
    out.mkdir(parents=True)
    from analysis import run
    summary=run(ROOT,out)
    # Published trace inputs, where present, are reanalyzed separately from synthetic data.
    from literature_analysis import analyze_extracted
    literature=analyze_extracted(ROOT,out)
    print(json.dumps({"out":str(out),"synthetic_blind_pass":summary["synthetic_blind_pass"],"literature":literature,"scientific_status":"PUBLIC_DATA_REANALYSIS_AND_UNCALIBRATED_LOCAL_BENCHMARK"},ensure_ascii=False))
    return 0

if __name__=="__main__": raise SystemExit(main())
