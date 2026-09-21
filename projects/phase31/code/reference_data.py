"""Explicit retrieval of public structural reference, retaining raw bytes and hashes."""
from datetime import datetime, timezone
from pathlib import Path
import urllib.request
from common import write_json, sha256

URLS = {"6SIS.pdb": "https://files.rcsb.org/download/6SIS.pdb",
        "LFE_ideal.sdf": "https://files.rcsb.org/ligands/download/LFE_ideal.sdf"}


def fetch(out):
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    entries = []
    for name, url in URLS.items():
        request = urllib.request.Request(url, headers={"User-Agent": "Phase31-evidence-benchmark/1"})
        with urllib.request.urlopen(request, timeout=45) as response:
            data = response.read(10_000_001)
        if len(data) > 10_000_000 or not data:
            raise ValueError("unexpected reference size")
        (out / name).write_bytes(data)
        entries.append(dict(file=name, source=url, bytes=len(data), sha256=sha256(out / name)))
    text = (out / "6SIS.pdb").read_text(encoding="ascii")
    if "6SIS" not in text.splitlines()[0] or "LFE" not in text:
        raise ValueError("unexpected PDB contents")
    chains, ligands = {}, {}
    for line in text.splitlines():
        if line.startswith("ATOM  "):
            chains[line[21]] = chains.get(line[21], 0) + 1
        elif line.startswith("HETATM") and line[17:20] == "LFE":
            key = line[21:27].strip()
            ligands[key] = ligands.get(key, 0) + 1
    report = dict(retrieved_utc=datetime.now(timezone.utc).isoformat(), evidence="PUBLIC_STRUCTURE_NOT_PREPARED_SIMULATION",
                  entries=entries, polymer_atom_records_by_author_chain=chains, LFE_atom_records=ligands,
                  reference_target="human BRD4 BD2", E3_receptor="VHL (with Elongin B/C)",
                  resolution_A=3.5, ligand_class="macrocyclic heterobifunctional PROTAC, not a proven molecular glue",
                  intended_author_chains={"BRD4":"A", "ElonginB":"B", "ElonginC":"C", "VHL":"D"},
                  assembly_selection="Requires assembly/contacts/alternate-location review before preparation",
                  full_E3_E2_ubiquitin_system_present=False, covalent_site=None,
                  prepared_for_dynamics=False)
    write_json(out / "provenance.json", report)
    return report
