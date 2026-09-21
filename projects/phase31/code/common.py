"""Shared provenance and numbered module loading."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def load_module(name):
    key = "phase31_" + name
    if key in sys.modules:
        return sys.modules[key]
    spec = importlib.util.spec_from_file_location(key, ROOT / "code" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[key] = module
    spec.loader.exec_module(module)
    return module


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def source_hashes():
    paths = sorted(ROOT.rglob("*.py")) + sorted((ROOT / "swarm" / "prompts").glob("*.md"))
    return {p.relative_to(ROOT).as_posix(): sha256(p) for p in paths}


class EvidenceBlocked(RuntimeError):
    pass
