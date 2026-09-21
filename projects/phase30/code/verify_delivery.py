"""Execute tests and verify one immutable run. Writes local verification evidence."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from PIL import Image
from common import sha, dump
from pipeline import PROJECT, source_hashes


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run",required=True,type=Path)
    args=parser.parse_args()
    run=args.run.resolve()
    before=source_hashes()
    repo=PROJECT.parents[1]
    env=os.environ.copy()
    env.update(PYTHONIOENCODING="utf-8",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1")
    commands=[[sys.executable,str(PROJECT/"run.py"),"--self-test"],
              [sys.executable,"-m","unittest","discover","-s","tools","-p","test_native_phases.py","-v"]]
    tests=[]
    combined=[]
    for command in commands:
        proc=subprocess.run(command,cwd=repo,env=env,capture_output=True,text=True,encoding="utf-8",
                            errors="replace",timeout=60)
        text=proc.stdout+proc.stderr
        count=re.search(r"Ran (\d+) tests?",text)
        tests.append(dict(suite="phase30" if "--self-test" in command else "native_dispatch",
                          exit_code=proc.returncode,count=int(count.group(1)) if count else None))
        combined.append(text)
    records=json.loads((run/"artifact_manifest.json").read_text(encoding="utf-8"))["files"]
    mismatches=[]
    for name,record in records.items():
        p=(run/name).resolve()
        if not p.is_relative_to(run) or not p.is_file() or sha(p)!=record["sha256"]:
            mismatches.append(name)
    status=json.loads((run/"run_status.json").read_text(encoding="utf-8"))
    images=[]
    for path in sorted((run/"figures").glob("*.png")):
        with Image.open(path) as img:
            images.append(dict(name=path.name,width=img.width,height=img.height,
                               dpi=img.info.get("dpi")))
    broken=[]
    for doc in [PROJECT/"README.md",PROJECT/"schemas/CONTRACTS.md",PROJECT/"reports/phase30_technical_report.md"]:
        text=re.sub(r"```.*?```","",doc.read_text(encoding="utf-8"),flags=re.S)
        for target in re.findall(r"\]\(([^)\s]+)\)",text):
            if not re.match(r"[a-zA-Z][\w+.-]*:|#|/",target) and not (doc.parent/target.split("#")[0]).exists():
                broken.append(dict(document=doc.name,target=target))
    (PROJECT/"results").mkdir(exist_ok=True)
    # These are generated test outputs, not scientific observations.
    (PROJECT/"results/software_tests.txt").write_text("\n".join(combined),encoding="utf-8")
    # This report is written last; its own link is expected to be absent before first write.
    broken=[b for b in broken if b["target"]!="results/software_verification.json"]
    ok=(all(r["exit_code"]==0 for r in tests) and not mismatches and not broken and
        before==source_hashes()==status["source_sha256"] and status["status"]=="COMPLETED_DEMO" and
        len(images)==3 and all(abs(i["dpi"][0]-300)<1 for i in images))
    dump(PROJECT/"results/software_verification.json",
         dict(status="PASS" if ok else "FAIL", tests=tests, source_sha256=before,
              run_artifact_files_checked=len(records),hash_mismatches=mismatches,
              source_unchanged=before==source_hashes(),tested_run_matches_source=before==status["source_sha256"],
              broken_project_links=broken,figures=images,
              visual_qa="Images inspected manually during authoring; no color-vision simulation",
              scope="software, hashes, image metadata and Phase30 links only; not full repository or science acceptance"))
    print(json.dumps(dict(status="PASS" if ok else "FAIL",tests=tests,hash_mismatches=mismatches,broken_links=broken)))
    return int(not ok)


if __name__=="__main__": raise SystemExit(main())
