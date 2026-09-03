#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
setup_and_download.py
=====================
One-shot bootstrap for a molecular machine-learning workspace.

Pipeline (idempotent, safe to re-run):

  STEP 1  Auto-install dependencies (RDKit, PyTorch, PyTorch Geometric,
          DeepChem, pandas, scikit-learn, matplotlib) via pip.
  STEP 2  Create the standard workspace tree:
          ./data/{raw,processed}  ./models/checkpoints  ./src  ./logs
  STEP 3  Download the classic MoleculeNet benchmarks from the official
          DeepChem S3 bucket and store them locally:
            * ESOL     (delaney-processed.csv, 1128 molecules)   -> CSV + SDF
            * FreeSolv (SAMPL.csv,            642 molecules)    -> CSV + SDF
            * QM9      (qm9.csv,             ~130k molecules)   -> CSV
  STEP 4  Environment self-check: GPU probe, matmul micro-benchmark,
          RDKit parsing / sanitization / descriptor tests, an
          ETKDGv3 + MMFF94 conformer test, and a PyTorch-Geometric
          smoke test. Every result is serialized to ./results.json
          before exit (fault-tolerant execution protocol).
  STEP 5  Optional clean power-off: with --auto_shutdown the machine is
          shut down --shutdown_delay seconds (default 60) AFTER all
          self-checks pass. A failing self-check NEVER powers off.

Examples
--------
    python setup_and_download.py
    python setup_and_download.py --auto_shutdown                  # power off 60 s after success
    python setup_and_download.py --auto_shutdown --shutdown_delay 120
    python setup_and_download.py --skip_install --skip_download   # self-check only
    python setup_and_download.py --skip_qm9                       # skip the ~30 MB QM9 CSV
    python setup_and_download.py --sdf_3d                         # 3D SDF (ETKDGv3 + MMFF94)

Notes
-----
* Heavy libraries (rdkit / pandas / torch ...) are imported lazily inside
  functions: on a fresh interpreter the imports must happen AFTER STEP 1
  has installed the packages.
* Mainland-China users can accelerate pip with
  ``--pip_index_url https://pypi.tuna.tsinghua.edu.cn/simple`` and CUDA torch
  wheels with ``--torch_index_url https://download.pytorch.org/whl/cu121``.
* Cancel a scheduled shutdown with ``shutdown /a`` (Windows) or
  ``sudo shutdown -c`` (Linux).
"""

from __future__ import annotations

import argparse
import datetime as _dt
import importlib
import itertools
import json
import os
import platform
import subprocess
import sys
import time
import traceback
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

# Windows consoles frequently default to a legacy code page (cp936/cp437);
# reconfigure to UTF-8 so logging never dies on non-ASCII output.
if sys.platform.startswith("win"):
    for _stream in (sys.stdout, sys.stderr):
        try:
            _stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

RESULTS: Dict[str, Any] = {
    "timestamp_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
    "platform": f"{platform.system()} {platform.release()} ({platform.machine()})",
    "python": sys.version.split()[0],
    "packages": {},
    "workspace": {},
    "datasets": {},
    "gpu": {},
    "benchmark": {},
    "rdkit_tests": {},
    "pyg_smoke_test": {},
    "all_checks_passed": False,
}

BANNER = r"""
------------------------------------------------------------------------
 setup_and_download.py -- molecular-ML workspace bootstrap
 deps -> workspace -> MoleculeNet(ESOL / FreeSolv / QM9) -> self-check
------------------------------------------------------------------------
"""

REFERENCE_MOLECULES: Dict[str, str] = {
    "ethanol": "CCO",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "caffeine": "Cn1cnc2c(=O)n(C)c(=O)n2C",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
}

# (import_name, strict pip spec, looser fallback spec)
CORE_PACKAGES: List[Tuple[str, str, str]] = [
    ("numpy", "numpy>=1.24", "numpy"),
    ("pandas", "pandas>=2.0", "pandas"),
    ("sklearn", "scikit-learn>=1.3", "scikit-learn"),
    ("rdkit", "rdkit>=2023.9", "rdkit"),
    ("torch", "torch>=2.3.0", "torch"),
    ("torch_geometric", "torch_geometric>=2.5.0", "torch_geometric"),
]

DATASET_REGISTRY: Dict[str, Dict[str, Any]] = {
    "esol": {
        "url": "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv",
        "raw_name": "delaney-processed.csv",
        "smiles_col": "smiles",
        "target_hint": "measured log sol",
        "target_name": "esol_log_solubility",
        "write_sdf": True,
        "min_parse_rate": 0.98,
    },
    "freesolv": {
        "url": "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/SAMPL.csv",
        "raw_name": "SAMPL.csv",
        "smiles_col": "smiles",
        "target_hint": "expt",
        "target_name": "freesolv_dg_kcal_mol",
        "write_sdf": True,
        "min_parse_rate": 0.98,
    },
    "qm9": {
        "url": "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/qm9.csv",
        "raw_name": "qm9.csv",
        "smiles_col": "smiles",
        "target_hint": "u0",
        "target_name": "qm9_u0_hartree",
        "write_sdf": False,  # qm9.sdf (~430 MB, 3D) intentionally not pulled by default
        "min_parse_rate": 0.85,  # a few QM9 nitro SMILES are known to resist RDKit
    },
}


class _Tee:
    """Duplicate everything written to stdout into a log file stream."""

    def __init__(self, *streams: Any) -> None:
        self._streams = streams

    def write(self, data: str) -> int:
        for stream in self._streams:
            try:
                stream.write(data)
                stream.flush()
            except Exception:
                pass
        return len(data)

    def flush(self) -> None:
        for stream in self._streams:
            try:
                stream.flush()
            except Exception:
                pass


class _Fmt(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Help formatter that keeps the epilog layout and shows defaults."""


def _hr(title: str = "") -> None:
    """Print a horizontal section separator, optionally with a title."""
    if title:
        print(f"\n{'=' * 72}\n{title}\n{'=' * 72}")
    else:
        print("=" * 72)


def init_logging(logs_dir: Path) -> Path:
    """Mirror all stdout output into ``logs/setup_<timestamp>.log``.

    The file handle is intentionally never closed explicitly: the _Tee
    wrapper flushes on every write and the OS finalizes it at exit.

    Args:
        logs_dir: directory for the log file (created if needed).

    Returns:
        Path of the log file being written.
    """
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_path = logs_dir / f"setup_{_dt.datetime.now():%Y%m%d_%H%M%S}.log"
    fh = open(log_path, "w", encoding="utf-8")
    sys.stdout = _Tee(sys.stdout, fh)
    return log_path


def parse_args() -> argparse.Namespace:
    """Parse command-line options."""
    parser = argparse.ArgumentParser(
        prog="setup_and_download.py",
        description="Bootstrap a molecular-ML workspace: install deps, create the standard "
                    "directory tree, download MoleculeNet benchmarks (ESOL/FreeSolv/QM9), "
                    "run an environment self-check and (optionally) power off on success.",
        formatter_class=_Fmt,
        epilog="Examples:\n"
               "  python setup_and_download.py --auto_shutdown\n"
               "  python setup_and_download.py --skip_install --skip_download\n"
               "  python setup_and_download.py --pip_index_url https://pypi.tuna.tsinghua.edu.cn/simple\n",
    )
    parser.add_argument("--auto_shutdown", action="store_true",
                        help="schedule a clean system power-off --shutdown_delay seconds AFTER all self-checks pass")
    parser.add_argument("--shutdown_delay", type=int, default=60,
                        help="seconds to wait before power-off (Windows uses /t, Linux uses +minutes)")
    parser.add_argument("--skip_install", action="store_true", help="skip dependency installation")
    parser.add_argument("--skip_download", action="store_true", help="skip dataset download/validation")
    parser.add_argument("--skip_qm9", action="store_true", help="skip the large QM9 CSV (~30 MB)")
    parser.add_argument("--skip_deepchem", action="store_true", help="skip the optional DeepChem install")
    parser.add_argument("--sdf_3d", action="store_true",
                        help="embed ETKDGv3+MMFF94 3D conformers into the SDF exports (slower)")
    parser.add_argument("--pip_index_url", default=None, help="custom PyPI index / mirror URL")
    parser.add_argument("--torch_index_url", default=None,
                        help="custom torch wheel index, e.g. https://download.pytorch.org/whl/cu121")
    parser.add_argument("--data_root", default="data", help="root folder for datasets")
    return parser.parse_args()


# --------------------------------------------------------------------------- #
# STEP 1 -- dependencies
# --------------------------------------------------------------------------- #
def pip_install(spec: str, index_url: Optional[str] = None) -> bool:
    """Install one pip requirement spec with the current interpreter.

    Args:
        spec: pip requirement spec, e.g. ``"torch>=2.3.0"``.
        index_url: optional ``--index-url`` for mirrors / CUDA wheel repos.

    Returns:
        True if pip exited with code 0.
    """
    cmd = [sys.executable, "-m", "pip", "install", "--no-warn-script-location", spec]
    if index_url:
        cmd += ["--index-url", index_url]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, errors="replace", timeout=3600)
    except (subprocess.TimeoutExpired, OSError) as exc:
        print(f"      pip could not be executed: {exc}")
        return False
    if proc.returncode != 0:
        tail = "\n".join((proc.stderr or proc.stdout or "").strip().splitlines()[-12:])
        print(f"      pip exit code {proc.returncode}; last output:\n{tail}")
        return False
    return True


def ensure_package(import_name: str, pip_spec: str, *, fallback: Optional[str] = None,
                   index_url: Optional[str] = None, optional: bool = False) -> bool:
    """Import ``import_name``; if missing, pip-install ``pip_spec`` (then ``fallback``).

    Args:
        import_name: module name used by ``importlib``.
        pip_spec: strict pip requirement spec tried first.
        fallback: looser spec (usually unpinned) retried on failure.
        index_url: optional ``--index-url`` forwarded to pip.
        optional: if True, installation failure is tolerated (warning only).

    Returns:
        True if the package ends up importable.

    Raises:
        RuntimeError: for a required package that could not be installed.
    """
    try:
        mod = importlib.import_module(import_name)
    except ImportError:
        mod = None
    if mod is not None:
        version = getattr(mod, "__version__", "unknown")
        RESULTS["packages"][import_name] = {"status": "present", "version": version}
        print(f"[deps] {import_name:<18} present   (v{version})")
        return True

    for spec in [pip_spec] + ([fallback] if fallback else []):
        print(f"[deps] {import_name:<18} missing -> pip install '{spec}'")
        if pip_install(spec, index_url):
            importlib.invalidate_caches()
            try:
                mod = importlib.import_module(import_name)
            except ImportError:
                mod = None
            if mod is not None:
                version = getattr(mod, "__version__", "unknown")
                RESULTS["packages"][import_name] = {"status": "installed", "version": version}
                print(f"[deps] {import_name:<18} installed  (v{version})")
                return True
        print(f"[deps] {import_name:<18} still unavailable after '{spec}'")

    if optional:
        RESULTS["packages"][import_name] = {"status": "optional-miss"}
        print(f"[deps] {import_name:<18} OPTIONAL install failed -- continuing without it")
        return False
    RESULTS["packages"][import_name] = {"status": "failed"}
    raise RuntimeError(f"required package '{import_name}' could not be installed")


def install_dependencies(args: argparse.Namespace) -> None:
    """Install every core dependency plus the optional extras (DeepChem, matplotlib)."""
    index = args.pip_index_url
    for import_name, spec, fallback in CORE_PACKAGES:
        spec_index = args.torch_index_url if import_name == "torch" else index
        ensure_package(import_name, spec, fallback=fallback, index_url=spec_index)
    if args.skip_deepchem:
        RESULTS["packages"]["deepchem"] = {"status": "skipped"}
        print(f"[deps] {'deepchem':<18} skipped   (--skip_deepchem)")
    else:
        ensure_package("deepchem", "deepchem>=2.8.0", fallback="deepchem",
                       index_url=index, optional=True)
    ensure_package("matplotlib", "matplotlib>=3.7", fallback="matplotlib",
                   index_url=index, optional=True)


# --------------------------------------------------------------------------- #
# STEP 2 -- workspace
# --------------------------------------------------------------------------- #
def build_workspace(data_root: str) -> Dict[str, Path]:
    """Create the standard project directory tree.

    Args:
        data_root: root folder for datasets (``data/raw`` + ``data/processed``).

    Returns:
        Mapping of role name -> created Path.
    """
    _hr("STEP 2 -- workspace scaffolding")
    dirs = {
        "data_raw": Path(data_root) / "raw",
        "data_processed": Path(data_root) / "processed",
        "models": Path("models"),
        "model_checkpoints": Path("models") / "checkpoints",
        "src": Path("src"),
        "logs": Path("logs"),
    }
    for name, path in dirs.items():
        path.mkdir(parents=True, exist_ok=True)
        RESULTS["workspace"][name] = str(path)
    (Path("src") / "__init__.py").touch(exist_ok=True)
    for name, path in dirs.items():
        print(f"    {name:<16} -> {path}")
    return dirs


# --------------------------------------------------------------------------- #
# STEP 3 -- datasets
# --------------------------------------------------------------------------- #
def download_file(url: str, dest: Path, retries: int = 3, chunk: int = 1 << 20) -> Path:
    """Download ``url`` to ``dest`` with retries, resume-safe .part staging.

    Args:
        url: direct HTTP(S) source.
        dest: destination path (skipped if a non-empty file already exists).
        retries: number of attempts with exponential backoff.
        chunk: read-buffer size in bytes.

    Returns:
        The destination Path on success.

    Raises:
        RuntimeError: after all retries are exhausted.
    """
    if dest.exists() and dest.stat().st_size > 0:
        print(f"    {dest.name}: already on disk "
              f"({dest.stat().st_size / 1e6:.1f} MB) -- skipping download")
        return dest

    last_err: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (setup_and_download.py)"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                total = int(resp.headers.get("Content-Length", 0) or 0)
                tmp = dest.with_name(dest.name + ".part")
                done = 0
                t0 = time.time()
                last_report = 0.0
                with open(tmp, "wb") as fh:
                    while True:
                        block = resp.read(chunk)
                        if not block:
                            break
                        fh.write(block)
                        done += len(block)
                        if time.time() - last_report > 5.0:
                            last_report = time.time()
                            if total:
                                print(f"    {dest.name}: {done / 1e6:8.1f} / {total / 1e6:.1f} MB "
                                      f"({100 * done / total:5.1f} %)", flush=True)
                            else:
                                print(f"    {dest.name}: {done / 1e6:8.1f} MB downloaded", flush=True)
                os.replace(tmp, dest)
                elapsed = time.time() - t0
                print(f"    {dest.name}: finished {done / 1e6:.1f} MB in {elapsed:.1f} s "
                      f"({done / 1e6 / max(elapsed, 1e-9):.1f} MB/s)")
                return dest
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            last_err = exc
            wait = 2 ** attempt
            print(f"    attempt {attempt}/{retries} failed: {exc} -- retrying in {wait} s")
            time.sleep(wait)
    raise RuntimeError(f"could not download {url}: {last_err}")


def _find_column(columns: List[str], *hints: str) -> Optional[str]:
    """Find a dataframe column by exact-then-substring match of the hints.

    Args:
        columns: available column names.
        *hints: candidate names/fragments, in priority order.

    Returns:
        The matched column name, or None.
    """
    lowered = {c.lower(): c for c in columns}
    for hint in hints:
        if hint.lower() in lowered:
            return lowered[hint.lower()]
    for hint in hints:
        for low, original in lowered.items():
            if hint.lower() in low:
                return original
    return None


def rdkit_parse_rate(smiles: Iterable[str], max_n: int = 300) -> Tuple[float, int]:
    """Fraction of SMILES that RDKit can parse AND sanitize (MolFromSmiles != None).

    Args:
        smiles: iterable of SMILES strings.
        max_n: cap on the number of probed molecules (speed).

    Returns:
        (parse_rate, n_probed).
    """
    from rdkit import Chem  # local import: guaranteed installed by STEP 1

    n = ok = 0
    for smi in itertools.islice(iter(smiles), max_n):
        if not isinstance(smi, str) or not smi.strip():
            continue
        n += 1
        if Chem.MolFromSmiles(smi) is not None:
            ok += 1
    return (ok / n if n else 0.0), n


def write_sdf(df: "pd.DataFrame", target_name: str, sdf_path: Path, *, embed3d: bool = False) -> int:
    """Write molecules + target property to an MDL SDF file with RDKit.

    Args:
        df: DataFrame with ``smiles`` and ``target_{target_name}`` columns.
        target_name: dataset key used in the SD property name.
        sdf_path: output .sdf path.
        embed3d: if True, attach an ETKDGv3 conformer relaxed with MMFF94
            (UFF fallback); on embedding failure the 2D connection table
            is written instead.

    Returns:
        Number of molecules actually written.
    """
    import pandas as pd
    from rdkit import Chem
    from rdkit.Chem import AllChem

    n_written = 0
    writer = Chem.SDWriter(str(sdf_path))
    try:
        for smi, tgt in zip(df["smiles"], df[f"target_{target_name}"]):
            if not isinstance(smi, str):
                continue
            mol = Chem.MolFromSmiles(smi)  # None on parse OR sanitization failure
            if mol is None:
                continue
            try:
                Chem.SanitizeMol(mol)  # defensive re-sanitization (protocol rule 2)
            except Exception:
                continue
            mol.SetProp("_Name", f"{target_name}_{n_written:06d}")
            mol.SetProp(f"target_{target_name}",
                        "nan" if pd.isna(tgt) else f"{float(tgt):.6f}")
            if embed3d:
                mol_h = Chem.AddHs(mol)
                try:
                    params = AllChem.ETKDGv3()
                    params.randomSeed = 42
                    if AllChem.EmbedMolecule(mol_h, params) == 0:
                        if AllChem.MMFFHasAllMoleculeParams(mol_h):
                            AllChem.MMFFOptimizeMolecule(mol_h, maxIters=500)
                        else:
                            AllChem.UFFOptimizeMolecule(mol_h, maxIters=500)
                        mol = mol_h
                except Exception:
                    pass  # keep the 2D record rather than dropping the molecule
            writer.write(mol)
            n_written += 1
    finally:
        writer.close()
    return n_written


def process_dataset(name: str, meta: Dict[str, Any], raw_dir: Path, processed_dir: Path,
                    *, sdf_3d: bool) -> Dict[str, Any]:
    """Validate a raw MoleculeNet CSV and emit a standardized processed copy (+SDF).

    The processed CSV uses a uniform schema across all three datasets:
    ``smiles`` + ``target_<dataset_name>``, ready for downstream featurizers.

    Args:
        name: dataset key (esol / freesolv / qm9).
        meta: the DATASET_REGISTRY entry.
        raw_dir: folder containing the downloaded file.
        processed_dir: output folder for the standardized CSV/SDF.
        sdf_3d: forwarded to :func:`write_sdf`.

    Returns:
        Result record (also stored in RESULTS['datasets']).
    """
    import pandas as pd

    raw_path = raw_dir / meta["raw_name"]
    df = pd.read_csv(raw_path)
    smiles_col = meta["smiles_col"]
    if smiles_col not in df.columns:
        raise KeyError(f"column '{smiles_col}' not found; header={list(df.columns)[:10]}")

    target_col = _find_column(list(df.columns), meta["target_hint"])
    if target_col is None:
        raise KeyError(f"target column matching '{meta['target_hint']}' not found "
                       f"in {list(df.columns)[:10]}")

    rate, n_samp = rdkit_parse_rate(df[smiles_col].astype(str), max_n=300)

    out = pd.DataFrame({
        "smiles": df[smiles_col].astype(str),
        f"target_{meta['target_name']}": pd.to_numeric(df[target_col], errors="coerce"),
    })
    out_path = processed_dir / f"{name}.csv"
    out.to_csv(out_path, index=False)

    sdf_path: Optional[Path] = None
    n_sdf: Optional[int] = None
    if meta["write_sdf"]:
        sdf_path = processed_dir / f"{name}.sdf"
        n_sdf = write_sdf(out, meta["target_name"], sdf_path, embed3d=sdf_3d)

    passed = rate >= float(meta["min_parse_rate"])
    rec: Dict[str, Any] = {
        "status": "ok" if passed else "failed",
        "error": None if passed else
                 f"SMILES parse rate {rate:.3f} below threshold {meta['min_parse_rate']}",
        "url": meta["url"],
        "raw_file": str(raw_path),
        "raw_size_mb": round(raw_path.stat().st_size / 1e6, 2),
        "n_rows": int(len(df)),
        "target_column": target_col,
        "smiles_parse_rate_probed": round(rate, 4),
        "n_probed": n_samp,
        "processed_csv": str(out_path),
        "sdf": str(sdf_path) if sdf_path else None,
        "n_molecules_in_sdf": n_sdf,
    }
    if name == "qm9":
        rec["note"] = ("raw qm9.csv keeps ALL DeepChem columns; the ~430 MB qm9.sdf (3D) "
                       "is intentionally not downloaded by default")
    tag = "OK  " if passed else "FAIL"
    suffix = f"  sdf={n_sdf} mols" if n_sdf is not None else ""
    print(f"    [{tag}] {name:<8} rows={rec['n_rows']:>7}  size={rec['raw_size_mb']:>7.2f} MB  "
          f"parse_rate={rate:.3f}  target='{target_col}'{suffix}")
    return rec


# --------------------------------------------------------------------------- #
# STEP 4 -- self-check
# --------------------------------------------------------------------------- #
def gpu_report() -> Dict[str, Any]:
    """Probe PyTorch / CUDA / cuDNN availability and device inventory.

    A missing GPU is reported as information, NOT a failure (CPU-only is
    a valid configuration for these 2D benchmarks).

    Returns:
        GPU/CUDA status dict (also stored in RESULTS['gpu']).
    """
    info: Dict[str, Any] = {}
    try:
        import torch
    except ImportError as exc:
        info["error"] = f"torch not importable: {exc}"
        print(f"    GPU: cannot probe -- {info['error']}")
        return info

    info["torch_version"] = torch.__version__
    info["cuda_available"] = bool(torch.cuda.is_available())
    if info["cuda_available"]:
        info["device_count"] = torch.cuda.device_count()
        info["cuda_version"] = torch.version.cuda
        info["cudnn_version"] = torch.backends.cudnn.version()
        info["devices"] = []
        print(f"    GPU: {info['device_count']} CUDA device(s), "
              f"CUDA {info['cuda_version']}, cuDNN {info['cudnn_version']}")
        for i in range(info["device_count"]):
            props = torch.cuda.get_device_properties(i)
            device = {
                "index": i,
                "name": props.name,
                "total_memory_gb": round(props.total_memory / 1024 ** 3, 1),
                "compute_capability": f"{props.major}.{props.minor}",
            }
            info["devices"].append(device)
            print(f"      [{i}] {device['name']}  {device['total_memory_gb']} GB  "
                  f"sm_{device['compute_capability']}")
    else:
        info["note"] = "no CUDA device visible -- PyTorch runs on CPU (not a failure)"
        print("    GPU: none visible -- CPU mode (this does NOT fail the self-check)")
    return info


def matmul_benchmark() -> Dict[str, Any]:
    """Dense matmul micro-benchmark (TFLOP/s) on the fastest available device.

    Returns:
        Benchmark record dict.
    """
    try:
        import torch

        if torch.cuda.is_available():
            device, n, iters = "cuda", 4096, 20
            a = torch.randn(n, n, device="cuda")
            b = torch.randn(n, n, device="cuda")
            for _ in range(3):
                a @ b  # warm-up + kernel autotune
            torch.cuda.synchronize()
            t0 = time.perf_counter()
            for _ in range(iters):
                c = a @ b
            torch.cuda.synchronize()
            elapsed = time.perf_counter() - t0
            del a, b, c
            torch.cuda.empty_cache()
        else:
            device, n, iters = "cpu", 1024, 5
            a = torch.randn(n, n)
            b = torch.randn(n, n)
            a @ b  # warm-up
            t0 = time.perf_counter()
            for _ in range(iters):
                c = a @ b
            elapsed = time.perf_counter() - t0

        tflops = 2.0 * (n ** 3) * iters / elapsed / 1e12
        res = {"device": device, "matrix_n": n, "iters": iters,
               "seconds": round(elapsed, 3), "tflops": round(tflops, 2)}
        print(f"    matmul: {res['tflops']:.2f} TFLOP/s on {device} "
              f"(n={n}, {iters} iters, {res['seconds']} s)")
        return res
    except Exception as exc:
        print(f"    matmul benchmark failed: {exc}")
        return {"error": f"{type(exc).__name__}: {exc}"}


def rdkit_conformer_test() -> Dict[str, Any]:
    """ETKDGv3 embedding + MMFF94 relaxation of aspirin (protocol rule 3).

    Returns:
        Conformer test record dict.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    out: Dict[str, Any] = {"ok": False}
    try:
        mol = Chem.AddHs(Chem.MolFromSmiles(REFERENCE_MOLECULES["aspirin"]))
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        cid = AllChem.EmbedMolecule(mol, params)
        out["embed_etkdgv3_conformer_id"] = int(cid)
        if cid < 0:
            out["error"] = "ETKDGv3 embedding failed"
            return out
        if AllChem.MMFFHasAllMoleculeParams(mol):
            code = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            out["mmff94"] = "converged" if code == 0 else f"not_converged_code_{code}"
        else:
            out["mmff94"] = "params_unavailable"
        out["ok"] = True
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"
    return out


def rdkit_selfcheck() -> Tuple[Dict[str, Any], bool]:
    """Parse, sanitize and describe reference molecules; run the conformer test.

    Returns:
        (per-molecule results dict, all_ok flag).
    """
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

    results: Dict[str, Any] = {}
    all_ok = True
    print("    molecule parsing / sanitization / descriptors:")
    for name, smi in REFERENCE_MOLECULES.items():
        entry: Dict[str, Any] = {"smiles": smi, "ok": False}
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                entry["error"] = "MolFromSmiles returned None"
            else:
                Chem.SanitizeMol(mol)
                entry.update({
                    "ok": True,
                    "n_atoms": int(mol.GetNumAtoms()),
                    "n_bonds": int(mol.GetNumBonds()),
                    "canonical_smiles": Chem.MolToSmiles(mol),
                    "mw": round(Descriptors.MolWt(mol), 2),
                    "clogp": round(Crippen.MolLogP(mol), 2),
                    "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
                })
        except Exception as exc:
            entry["error"] = f"{type(exc).__name__}: {exc}"
        all_ok = all_ok and bool(entry["ok"])
        results[name] = entry
        status = "OK  " if entry["ok"] else "FAIL"
        print(f"      [{status}] {name:<9} atoms={entry.get('n_atoms', '?')!s:>3}  "
              f"MW={entry.get('mw', '?')!s:>7}  cLogP={entry.get('clogp', '?')!s:>6}  "
              f"TPSA={entry.get('tpsa', '?')!s:>6}")

    conf = rdkit_conformer_test()
    results["conformer_etkdgv3_mmff94"] = conf
    all_ok = all_ok and bool(conf.get("ok", False))
    conf_status = "OK  " if conf.get("ok") else "FAIL"
    print(f"      [{conf_status}] ETKDGv3 + MMFF94 conformer on aspirin "
          f"({conf.get('mmff94', conf.get('error', '?'))})")
    return results, bool(all_ok)


def pyg_smoke_test() -> Dict[str, Any]:
    """Build aspirin as a PyG graph and run one GCNConv forward pass.

    Returns:
        Smoke-test record dict.
    """
    out: Dict[str, Any] = {"ok": False}
    try:
        import torch
        import torch_geometric
        from rdkit import Chem
        from torch_geometric.data import Data
        from torch_geometric.nn import GCNConv

        mol = Chem.MolFromSmiles(REFERENCE_MOLECULES["aspirin"])
        x = torch.tensor([[atom.GetAtomicNum()] for atom in mol.GetAtoms()], dtype=torch.float)
        edge_index = torch.tensor(
            [[bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()] for bond in mol.GetBonds()],
            dtype=torch.long,
        ).t().contiguous()
        edge_index = torch.cat([edge_index, edge_index.flip(0)], dim=1)  # make undirected
        graph = Data(x=x, edge_index=edge_index)
        h = GCNConv(graph.num_node_features, 32)(graph.x, graph.edge_index)
        out.update({
            "ok": tuple(h.shape) == (graph.num_nodes, 32),
            "torch_geometric": torch_geometric.__version__,
            "n_nodes": int(graph.num_nodes),
            "n_directed_edges": int(graph.num_edges),
            "output_shape": list(h.shape),
        })
        print(f"    PyG smoke test: v{out['torch_geometric']} graph({out['n_nodes']} nodes, "
              f"{out['n_directed_edges']} directed edges) -> GCNConv out {out['output_shape']}")
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"
        print(f"    PyG smoke test FAILED: {out['error']}")
    return out


# --------------------------------------------------------------------------- #
# STEP 5 -- shutdown & results
# --------------------------------------------------------------------------- #
def schedule_shutdown(delay_s: int) -> Dict[str, Any]:
    """Schedule a clean, cancellable system power-off.

    Windows uses ``shutdown /s /t <delay>`` (seconds); POSIX systems use
    ``shutdown -h +<minutes>`` with a ``sudo -n`` retry for unprivileged
    shells (AutoDL / RunPod / AWS GPU boxes usually run as root).

    Args:
        delay_s: grace period in seconds before power-off.

    Returns:
        Record of what was scheduled and how to cancel it.
    """
    info: Dict[str, Any] = {"scheduled": False, "delay_seconds": int(delay_s)}
    minutes = max(1, int(round(delay_s / 60)))
    try:
        if sys.platform.startswith("win"):
            cmd = ["shutdown", "/s", "/t", str(int(delay_s)),
                   "/c", "setup_and_download.py: all checks passed -- powering off."]
            info["cancel_command"] = "shutdown /a"
        elif sys.platform == "darwin":
            cmd = ["sudo", "-n", "shutdown", "-h", f"+{minutes}"]
            info["cancel_command"] = "sudo killall shutdown"
        else:  # linux
            cmd = ["shutdown", "-h", f"+{minutes}"]
            info["cancel_command"] = "sudo shutdown -c"

        proc = subprocess.run(cmd, capture_output=True, text=True, errors="replace", timeout=30)
        info["command"] = " ".join(cmd)
        info["scheduled"] = proc.returncode == 0
        if not info["scheduled"] and not sys.platform.startswith("win"):
            sudo_cmd = ["sudo", "-n"] + cmd
            proc2 = subprocess.run(sudo_cmd, capture_output=True, text=True, errors="replace", timeout=30)
            info["command"] = " ".join(sudo_cmd)
            info["scheduled"] = proc2.returncode == 0
            if not info["scheduled"]:
                info["error"] = (proc2.stderr or proc2.stdout or "").strip()
        elif not info["scheduled"]:
            info["error"] = (proc.stderr or proc.stdout or "").strip()
    except Exception as exc:
        info["error"] = f"{type(exc).__name__}: {exc}"

    if info["scheduled"]:
        print("\n" + "!" * 62)
        print(f"  SYSTEM SHUTDOWN SCHEDULED -- power-off in ~{delay_s} s")
        print(f"  cancel with : {info['cancel_command']}")
        print("!" * 62 + "\n")
    else:
        print(f"\n[shutdown] could NOT be scheduled: {info.get('error', 'unknown reason')}")
    return info


def write_results_json(path: str = "results.json") -> Path:
    """Atomically serialize RESULTS to JSON (tmp file + os.replace).

    Args:
        path: output JSON path.

    Returns:
        The final JSON path.
    """
    out = Path(path)
    tmp = out.with_name(out.name + ".tmp")
    tmp.write_text(json.dumps(RESULTS, indent=2, ensure_ascii=False, default=str),
                   encoding="utf-8")
    os.replace(tmp, out)
    return out


def _print_final_summary() -> None:
    """Print the end-of-run summary block."""
    _hr("FINAL SUMMARY")
    print(f"  python            : {RESULTS['python']}  on  {RESULTS['platform']}")
    packages = ", ".join(f"{k}={v.get('version', v.get('status'))}"
                         for k, v in RESULTS["packages"].items())
    print(f"  packages          : {packages if packages else '(none recorded)'}")
    for name, rec in RESULTS["datasets"].items():
        print(f"  dataset {name:<8} : {rec.get('status', '-')}  "
              f"rows={rec.get('n_rows', '-')}  csv={rec.get('processed_csv', '-')}")
    gpu = RESULTS.get("gpu", {})
    if gpu.get("cuda_available"):
        print(f"  GPU               : CUDA available x{gpu.get('device_count')} "
              f"({', '.join(d['name'] for d in gpu.get('devices', []))})")
    else:
        print(f"  GPU               : {gpu.get('note', gpu.get('error', 'not probed'))}")
    print(f"  ALL CHECKS PASSED : {RESULTS['all_checks_passed']}")
    if "fatal_error" in RESULTS:
        print(f"  FATAL ERROR       : {RESULTS['fatal_error']}")
    print("  artifacts         : results.json | logs/setup_*.log | data/ | models/ | src/")


def main() -> int:
    """Run the full bootstrap pipeline; return the process exit code."""
    args = parse_args()
    exit_code = 0
    try:
        log_path = init_logging(Path("logs"))
        print(BANNER)
        print(f"started {RESULTS['timestamp_utc']} | platform: {RESULTS['platform']} "
              f"| python {RESULTS['python']}")
        print(f"argv: {' '.join(sys.argv)}")
        print(f"log : {log_path}")

        if not args.skip_install:
            _hr("STEP 1 -- dependency auto-install (pip)")
            install_dependencies(args)
        else:
            print("\nSTEP 1 skipped (--skip_install)")

        build_workspace(args.data_root)

        if not args.skip_download:
            _hr("STEP 3 -- download MoleculeNet datasets (official DeepChem S3 bucket)")
            raw_dir = Path(args.data_root) / "raw"
            processed_dir = Path(args.data_root) / "processed"
            for name, meta in DATASET_REGISTRY.items():
                if name == "qm9" and args.skip_qm9:
                    RESULTS["datasets"][name] = {"status": "skipped (--skip_qm9)"}
                    print("    [SKIP] qm9")
                    continue
                try:
                    download_file(meta["url"], raw_dir / meta["raw_name"])
                    RESULTS["datasets"][name] = process_dataset(
                        name, meta, raw_dir, processed_dir, sdf_3d=args.sdf_3d)
                except Exception as exc:
                    RESULTS["datasets"][name] = {"status": "failed",
                                                 "error": f"{type(exc).__name__}: {exc}"}
                    print(f"    [FAIL] {name}: {exc}")
        else:
            print("\nSTEP 3 skipped (--skip_download)")

        _hr("STEP 4 -- environment self-check")
        RESULTS["gpu"] = gpu_report()
        RESULTS["benchmark"]["matmul"] = matmul_benchmark()
        rdkit_results, rdkit_ok = rdkit_selfcheck()
        RESULTS["rdkit_tests"] = rdkit_results
        RESULTS["pyg_smoke_test"] = pyg_smoke_test()

        core_names = {"numpy", "pandas", "sklearn", "rdkit", "torch", "torch_geometric"}
        core_pkgs_ok = all(RESULTS["packages"].get(n, {}).get("status") in {"present", "installed"}
                           for n in core_names)
        datasets_ok = True
        if not args.skip_download:
            for rec in RESULTS["datasets"].values():
                status = str(rec.get("status", ""))
                if not status.startswith("skipped") and status != "ok":
                    datasets_ok = False
        pyg_ok = bool(RESULTS["pyg_smoke_test"].get("ok", False))
        RESULTS["all_checks_passed"] = bool(core_pkgs_ok and datasets_ok and rdkit_ok and pyg_ok)

        print(f"\n    core packages OK: {core_pkgs_ok} | datasets OK: {datasets_ok} | "
              f"RDKit OK: {rdkit_ok} | PyG OK: {pyg_ok}")
        print(f"    ALL CHECKS PASSED: {RESULTS['all_checks_passed']}")

    except Exception as exc:
        RESULTS["fatal_error"] = f"{type(exc).__name__}: {exc}"
        RESULTS["all_checks_passed"] = False
        print(f"\n[FATAL] {RESULTS['fatal_error']}")
        traceback.print_exc()
        exit_code = 1
    finally:
        # Shutdown policy: only on a fully passing run; results are always
        # serialized first (fault-tolerant execution protocol).
        if args.auto_shutdown:
            if RESULTS.get("all_checks_passed"):
                RESULTS["shutdown"] = schedule_shutdown(args.shutdown_delay)
            else:
                RESULTS["shutdown"] = {"scheduled": False,
                                       "reason": "self-check did not pass -- power-off suppressed"}
                print("\n[shutdown] --auto_shutdown given but checks FAILED -- power-off suppressed.")
        _print_final_summary()
        results_path = write_results_json()
        print(f"\nresults serialized to: {results_path}")
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
