"""Stdlib-only AST/mock tests. No engine imports, QC, or trajectories."""
import argparse
import ast
import hashlib
import json
import math
import os
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import traceback
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

SOURCE = Path(__file__).with_name("run_phase8_photochemical_dynamics.py")
TREE = ast.parse(SOURCE.read_text(encoding="utf-8"))


def load(*names, **extra):
    ns = dict(Path=Path, argparse=argparse, hashlib=hashlib, json=json, math=math,
              os=os, subprocess=subprocess, sys=sys, time=time, traceback=traceback,
              __file__=str(SOURCE), _log=Mock(), _warn=Mock(),
              np=SimpleNamespace(isfinite=math.isfinite))
    nodes = [n for n in TREE.body if isinstance(n, ast.Assign)
             and all(t.id.isupper() for target in n.targets for t in ast.walk(target)
                     if isinstance(t, ast.Name))]
    # Constants contain only stdlib operations; never execute module imports.
    nodes += [n for n in TREE.body if isinstance(n, ast.FunctionDef) and n.name in names]
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SOURCE), "exec"), ns)
    ns.update(extra)
    return ns


def valid_b():
    def rows(field, values):
        return [{field: v, "e0_eh": -1., "e1_eh": -.9999, "gap_eV": .0027} for v in values]
    return dict(basis="6-31g", phi_ci_deg=90.,
                torsion_scan=rows("phi", [180., 160., 140., 120., 100., 90., 80., 60., 40., 20., 0.]),
                stretch_scan_at_ci=rows("d", [-.1, -.05, 0., .05, .1]),
                branching_cuts={k: rows("t_ang", [-.2, -.1, 0., .1, .2]) for k in ("g", "h")},
                g_vector={"per_atom": [[1., 0., 0.]] * 4},
                h_vector={"per_atom": [[0., 1., 0.]] * 4},
                meci=dict(converged=True, gap_eV=.0027, gradF_norm=.01))


class PreflightTests(unittest.TestCase):
    def test_compiles_without_importing_engines(self):
        compile(TREE, str(SOURCE), "exec")

    def test_qc_validation_rejects_missing_nonfinite_and_false_convergence(self):
        validate = load("validate_qc_result", "_meci_gate")["validate_qc_result"]
        validate("8b", valid_b())
        for change in (lambda r: r["torsion_scan"].pop(),
                       lambda r: r["stretch_scan_at_ci"][0].update(e0_eh=float("nan")),
                       lambda r: r["meci"].update(gradF_norm=.05),
                       lambda r: r["meci"].update(converged=False),
                       lambda r: r.update(basis="sto-3g")):
            record = valid_b()
            change(record)
            with self.assertRaises(ValueError):
                validate("8b", record)
        for record in ({}, None, {"error": "worker failed"}):
            with self.assertRaises(ValueError):
                validate("8b", record)

    def test_8a_full_grid_and_states_required(self):
        ns = load("validate_qc_result")
        record = dict(basis="def2-svp", e_scf_eh=-1., states={
            spin: [dict(root=i, dE_eV=2., f_osc=.1) for i in range(1, 11)]
            for spin in ("singlets", "triplets")}, torsion_scan={"points": [
                dict(phi=p, e0_eh=-1., s1_eV=2., s2_eV=3.) for p in ns["SCAN_AZO_PHI"]]})
        ns["validate_qc_result"]("8a", record)
        record["states"]["triplets"].pop()
        with self.assertRaises(ValueError):
            ns["validate_qc_result"]("8a", record)

    def test_incomplete_scan_cannot_advance_checkpoint(self):
        fn = load("require_scan_points")["require_scan_points"]
        rows = valid_b()["stretch_scan_at_ci"]
        fn(rows, "d", [-.1, -.05, 0., .05, .1])
        with self.assertRaisesRegex(ValueError, "checkpoint retained"):
            fn(rows[:-1], "d", [-.1, -.05, 0., .05, .1])
        worker = next(n for n in TREE.body if isinstance(n, ast.FunctionDef) and n.name == "worker_8b")
        text = ast.unparse(worker)
        self.assertIn('ck.get(\'identity\') != identity', text)
        self.assertEqual(text.count("require_scan_points("), 3)

    def test_identity_changes_with_geometry_basis_and_smoke(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            geo = root / "g.xyz"
            geo.write_text("fixture")
            scan = root / "scan.json"
            scan.write_text(json.dumps({key: {"a": {"file": str(geo)}}
                                        for key in ("scan_torsion", "scan_stretch")}))
            args = SimpleNamespace(xyz=str(geo), scan_t_json=str(scan), scan_s_json=str(scan), smoke=False)
            fn = load("checkpoint_identity", "parse_json")["checkpoint_identity"]
            original = fn(args, "6-31g")
            self.assertNotEqual(original, fn(args, "sto-3g"))
            args.smoke = True
            self.assertNotEqual(original, fn(args, "6-31g"))
            args.smoke = False
            geo.write_text("changed fixture")
            self.assertNotEqual(original, fn(args, "6-31g"))

    def test_timeout_retains_both_streams_and_thread_dll_budget(self):
        with tempfile.TemporaryDirectory() as tmp:
            ns = load("run_worker", OUT=Path(tmp), wait_for_qc_memory=Mock())
            error = subprocess.TimeoutExpired("worker", 5, output=b"last point", stderr=b"failure")
            with patch.object(subprocess, "run", side_effect=error) as run, patch.dict(
                    os.environ, {"PYTHONPATH": "wrong", "PYTHONHOME": "wrong"}):
                with self.assertRaises(subprocess.TimeoutExpired):
                    ns["run_worker"]("8b", [], timeout=5, threads=1)
            ns["wait_for_qc_memory"].assert_called_once_with()
            env = run.call_args.kwargs["env"]
            self.assertNotIn("PYTHONHOME", env)
            self.assertNotIn("PYTHONPATH", env)
            self.assertEqual(env["OMP_NUM_THREADS"], "1")
            self.assertEqual(run.call_args.args[0][-2:], ["--threads", "1"])
            self.assertEqual(env["PATH"].split(os.pathsep)[0],
                             str(Path(ns["ENV_PY_QC"]).parent / "Library" / "bin"))
            self.assertIn("last pointfailure", (Path(tmp) / "worker_8b.log").read_text())

    def test_memory_hold_prevents_spawn(self):
        ns = load("run_worker", wait_for_qc_memory=Mock(side_effect=RuntimeError("hold")))
        with patch.object(subprocess, "run") as run:
            with self.assertRaisesRegex(RuntimeError, "hold"):
                ns["run_worker"]("8b", [])
            run.assert_not_called()

    def test_xtb_uses_unpaired_electrons_and_rejects_nonzero_exit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            geometry = root / "g.xyz"
            geometry.write_text("fixture")
            fn = load("run_xtb", shutil=shutil)["run_xtb"]
            with patch.object(subprocess, "run", return_value=SimpleNamespace(
                    returncode=2, stdout="partial", stderr="bad option")) as run:
                with self.assertRaisesRegex(RuntimeError, "xTB rc=2"):
                    fn("xtb.exe", geometry, root / "job", opt=True)
            command = run.call_args.args[0]
            self.assertNotIn("--mult", command)
            self.assertEqual(command[command.index("--uhf") + 1], "0")
            self.assertEqual((root / "job" / "xtb.stderr").read_text(), "bad option")

    def test_cli_thread_default_and_required_worker_options(self):
        with tempfile.TemporaryDirectory() as tmp:
            ns = load("main", OUT=Path(tmp), FIG=Path(tmp), worker_8a=Mock())
            with patch.dict(os.environ, {"OMP_NUM_THREADS": "1"}), patch.object(
                    sys, "argv", ["driver", "--worker", "8a", "--xyz", "x", "--res", "r"]):
                ns["main"]()
            self.assertEqual(ns["worker_8a"].call_args.args[0].threads, 1)
            for argv in (["driver", "--threads", "0"], ["driver", "--worker", "8b"]):
                with patch.object(sys, "argv", argv), patch("sys.stderr"), self.assertRaises(SystemExit) as error:
                    ns["main"]()
                self.assertEqual(error.exception.code, 2)

    def test_failed_rerun_removes_stale_module_and_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            master = root / "phase8_results.json"
            master.write_text(json.dumps({"all_stages_ok": True, "module_8b": valid_b()}))
            ns = load("main", "parse_json", "dump_json", OUT=root, FIG=root,
                      stage_8b=Mock(side_effect=RuntimeError("failed")))
            with patch.object(sys, "argv", ["driver", "--stage", "8B"]):
                self.assertEqual(ns["main"](), 1)
            saved = json.loads(master.read_text())
            self.assertNotIn("module_8b", saved)
            self.assertNotIn("all_stages_ok", saved)
            self.assertIn("8B", saved["errors"])

    def test_fig_only_failure_is_persisted_and_retry_clears_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            master = root / "phase8_results.json"
            master.write_text("{}")
            ns = load("main", "parse_json", "dump_json", OUT=root, FIG=root,
                      render_phase8_figures=Mock())
            renderer = Mock(side_effect=RuntimeError("figure failure"))
            ns["render_phase8_figures"] = renderer
            with patch.dict(sys.modules, {"phase8_figures": SimpleNamespace(render_all=renderer)}), patch.object(
                    sys, "argv", ["driver", "--fig_only"]):
                self.assertEqual(ns["main"](), 1)
                self.assertIn("fig", json.loads(master.read_text())["errors"])
                renderer.side_effect = None
                self.assertEqual(ns["main"](), 0)
            self.assertNotIn("fig", json.loads(master.read_text())["errors"])

    def test_renderer_silent_skip_is_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            ns = load("render_phase8_figures", FIG=Path(tmp))
            with patch.dict(sys.modules, {"phase8_figures": SimpleNamespace(render_all=Mock())}):
                with self.assertRaisesRegex(RuntimeError, "missing/not refreshed"):
                    ns["render_phase8_figures"]()

    def test_failed_cached_qc_blocks_dynamics(self):
        run = Mock()
        ns = load("stage_8c", "validate_qc_result", "_meci_gate", run_fssh=run)
        with self.assertRaises(ValueError):
            ns["stage_8c"](SimpleNamespace(smoke=False), {}, {})
        run.assert_not_called()

    def test_td_spin_and_cube_indices_match_installed_api(self):
        source = ast.unparse(TREE)
        self.assertNotIn("tdscf_singlet", source)
        self.assertIn("'tdscf_triplets': 'NONE' if singlet else 'ONLY'", source)
        self.assertIn("'cubeprop_orbitals': [nocc, nocc + 1]", source)
        self.assertIn("nt = singlet_ntos.get(root)", source)
        self.assertIn("Ca_view[:] = original_ca", source)


if __name__ == "__main__":
    unittest.main()
