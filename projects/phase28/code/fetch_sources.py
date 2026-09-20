"""Fetch the two public figure-source workbooks and verify recorded SHA256."""
import argparse
import hashlib
import json
from pathlib import Path
import urllib.request

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--out",type=Path,required=True); a=p.parse_args()
    if a.out.exists(): p.error("Use a new directory; source files will not be overwritten")
    a.out.mkdir(parents=True)
    root=Path(__file__).resolve().parents[1]
    manifest=json.loads((root/"data/literature/source_manifest.json").read_text(encoding="utf-8"))
    for record in manifest:
        if record.get("reproduce_extraction") is not True: continue
        request=urllib.request.Request(record["url"],headers={"User-Agent":"Phase28-public-source-audit/1.0"})
        with urllib.request.urlopen(request,timeout=120) as response: data=response.read()
        checksum=hashlib.sha256(data).hexdigest()
        if checksum!=record["sha256"]: raise RuntimeError("Source hash changed: "+record["filename"])
        (a.out/record["filename"]).write_bytes(data)
        print(record["filename"],checksum)

if __name__=="__main__": main()
