"""Compatibility entry point for Phase 26; see README.md."""
import importlib.util
from pathlib import Path
import sys

_path = Path(__file__).resolve().parents[2] / 'tools' / 'run_phase.py'
_spec = importlib.util.spec_from_file_location('ai4chem_runner', _path)
_runner = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_runner)

if __name__ == '__main__':
    raise SystemExit(_runner.main(['26', *sys.argv[1:]]))
