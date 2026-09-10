"""Cheap launcher/prerequisite regressions; no molecular computations."""
import ast
import datetime
import hashlib
import json
import struct
import math
import os
from pathlib import Path
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch
import numpy as np

ROOT = Path(__file__).resolve().parent


def function(phase, name, **namespace):
    path = next(ROOT.glob(f"run_phase{phase}_*.py"))
    tree = ast.parse(path.read_text(encoding="utf-8"))
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
    namespace.update(Path=Path, os=os, time=time, __file__=str(path),
                     _log=lambda *a: None, wait_for_qc_memory=lambda: None)
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), namespace)
    return namespace[name]


class RerunTests(unittest.TestCase):
    def test_future_dcd_writer_records_actual_stride(self):
        from openmm import app, unit
        tree = ast.parse(next(ROOT.glob("run_phase6_*.py")).read_text(encoding="utf-8"))
        assignment = next(node for node in ast.walk(tree) if isinstance(node, ast.Assign)
                          and isinstance(node.value, ast.Call)
                          and isinstance(node.value.func, ast.Attribute)
                          and node.value.func.attr == "DCDFile")
        topology = app.Topology()
        residue = topology.addResidue("LIG", topology.addChain())
        topology.addAtom("C", app.element.carbon, residue)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "trajectory.dcd"
            ns = dict(app=app, unit=unit, sim=SimpleNamespace(topology=topology),
                      dcd_path=str(path), mode="wb", DT_FS=2., dcd_every=5000, start=0)
            exec(compile(ast.Module(body=[assignment], type_ignores=[]), "<DCD-writer>", "exec"), ns)
            dcd = ns["dcd"]
            dcd.writeModel([[0., 0., 0.]])
            dcd.writeModel([[.1, .1, .1]])
            dcd._file.close()
            result = function(6, "read_dcd_timing")(path)
            self.assertEqual(result["steps"], [5000, 10000])
            np.testing.assert_allclose(result["times_ps"], [10., 20.], atol=1e-5)

    def test_dcd_copy_corrects_timestamps_without_changing_coordinates(self):
        reader = function(6, "read_dcd_timing")
        correct = function(6, "correct_dcd_timing_copy", read_dcd_timing=reader,
                           hashlib=hashlib, json=json)
        def record(block):
            return struct.pack("<i", len(block)) + block + struct.pack("<i", len(block))
        header = bytearray(b"CORD" + struct.pack("<20i", 2, 0, 1, 1, *([0] * 16)))
        struct.pack_into("<f", header, 40, .002 / .04888821)
        original = record(header) + record(b"fixture") + record(struct.pack("<i", 1))
        original += b"".join(record(struct.pack("<f", value)) for value in (1., 2., 3., 4., 5., 6.))
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); source = root / "raw.dcd"; dest = root / "fixed.dcd"
            source.write_bytes(original)
            progress = root / "progress.csv"
            progress.write_text("step,time_ps\n5000,10\n10000,20\n")
            result = correct(source, dest, progress)
            self.assertEqual(source.read_bytes(), original)
            self.assertEqual(dest.read_bytes()[24:], original[24:])
            self.assertTrue(result["coordinate_bytes_identical"])
            np.testing.assert_allclose(reader(dest)["times_ps"], [10., 20.], atol=1e-5)
            self.assertEqual(reader(dest)["last_step"], 10000)
            with self.assertRaises(ValueError):
                correct(source, source, progress)

    def test_implicit_resume_requires_completed_leg_and_fingerprints_artifacts(self):
        fn = function(6, "load_completed_explicit", json=json, hashlib=hashlib,
                      _dt=datetime, DEPOSIT_STEPS_EXPL=500)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "results_phase6"
            (out / "cache").mkdir(parents=True)
            prior = {"explicit": {"steps": 1000, "target_steps": 1000, "n_hills": 2}}
            (out / "phase6_results.json").write_text(json.dumps(prior))
            for name in ("metadyn_state.json", "prod_explicit.chk", "progress_explicit.csv",
                         "traj_explicit.dcd", "solvated_equilibrated.pdb",
                         "cache/charges_reactant.txt", "cache/charges_product.txt",
                         "cache/morse_fit.json"):
                (out / name).write_text("fixture")
            (root / "rerun_status.json").write_text("fixture")
            (root / Path(fn.__globals__["__file__"]).name).write_text("original source")
            state = {"explicit": {"steps_done": 500, "hills": [[1., .15, .1]]}}
            with patch.object(Path, "cwd", return_value=root):
                with self.assertRaisesRegex(ValueError, "incomplete"):
                    fn(out, state, 1000)
                state["explicit"] = {"steps_done": 1000, "hills": [[1., .15, .1]] * 2}
                restored, provenance = fn(out, state, 1000)
            self.assertEqual(restored, prior)
            self.assertEqual(provenance["completed_explicit_steps"], 1000)
            for artifact in provenance["artifacts"]:
                self.assertEqual(artifact["sha256"],
                                 hashlib.sha256(Path(artifact["source"]).read_bytes()).hexdigest())

    def test_memory_gate_precedes_psi4_process_creation(self):
        sub = SimpleNamespace(run=Mock())
        for phase, name in ((7, "run_cmd"), (8, "run_worker")):
            fn = function(phase, name, subprocess=sub, ENV_PY_QC=str(ROOT / "python.exe"), OUT=ROOT)
            fn.__globals__["wait_for_qc_memory"] = Mock(side_effect=RuntimeError("memory hold"))
            with self.assertRaisesRegex(RuntimeError, "memory hold"):
                if phase == 7:
                    fn([str(ROOT / "python.exe"), "driver.py", "--worker", "7B"])
                else:
                    fn("8b", [])
        sub.run.assert_not_called()

    def test_spin_square_pure_singlet_doublet_and_triplet(self):
        fn = function(7, "_spin_square", np=np)
        for nalpha, nbeta, expected in ((1, 1, 0.), (2, 1, .75), (3, 1, 2.)):
            wf = SimpleNamespace(
                S=lambda: np.eye(3), nalpha=lambda: nalpha, nbeta=lambda: nbeta,
                Da_subset=lambda _: np.diag([1.] * nalpha + [0.] * (3-nalpha)),
                Db_subset=lambda _: np.diag([1.] * nbeta + [0.] * (3-nbeta)))
            self.assertAlmostEqual(fn(wf), expected)

    def test_phase6_keeps_small_hills_and_rejects_silent_saturation(self):
        meta = SimpleNamespace(getCollectiveVariableValues=lambda _: (.15, .2))
        ctx = SimpleNamespace(setParameter=Mock())
        fn = function(6, "deposit", math=math, N_MAX_HILLS=2,
                      sigma_constants=lambda: (.005, .08),
                      bias_eval=lambda *args: 10.)
        hills = []
        fn(meta, ctx, hills, 1., 1.)
        self.assertEqual(len(hills), 1)
        self.assertLess(hills[0][0], .001)
        fn(meta, ctx, hills, 1., 1.)
        with self.assertRaisesRegex(RuntimeError, "capacity exhausted"):
            fn(meta, ctx, hills, 1., 1.)
        self.assertEqual(len(hills), 2)

    def test_inverse_bfgs_satisfies_secant_and_positive_curvature(self):
        fn = function(8, "_inverse_bfgs_update", np=np)
        initial = np.eye(3) * .5
        step = np.array([.2, -.1, .3])
        delta_gradient = np.diag([2., 3., 4.]) @ step
        updated = fn(initial, step, delta_gradient)
        np.testing.assert_allclose(updated @ delta_gradient, step)
        np.testing.assert_allclose(updated, updated.T)
        self.assertTrue(np.all(np.linalg.eigvalsh(updated) > 0))
        np.testing.assert_array_equal(fn(initial, step, -step), initial)

    def test_meci_requires_finite_gap_and_gradient_with_original_gates(self):
        fn = function(8, "_meci_gate", np=np, MECI_GAP_EV_GATE=.05, MECI_GRAD_GATE=.05)
        self.assertTrue(fn(.01, .01))
        for gap, gradient in ((.01, .1), (.1, .01), (-.1, .01),
                              (.01, np.nan), (np.inf, .01), (.05, .01), (.01, .05)):
            self.assertFalse(fn(gap, gradient))

    def test_phase7_triplet_contamination_uses_spin_offset(self):
        path = next(ROOT.glob("run_phase7_*.py"))
        tree = ast.parse(path.read_text(encoding="utf-8"))
        loop = next(n for n in ast.walk(tree) if isinstance(n, ast.For)
                    and isinstance(n.target, ast.Tuple)
                    and [getattr(v, "id", None) for v in n.target.elts] == ["name", "ss"])
        seen = []
        namespace = dict(np=np, Rs=[1.5, 2.0, 3.0], s2_uhf_t=[2., 2.4, None],
                         s2_uks_t=[2.1, 2.2, 2.5], S2_GATE=.3, r_crit={},
                         interp_crossing=lambda x, y, gate: (seen.append(y), None))
        exec(compile(ast.Module(body=[loop], type_ignores=[]), str(path), "exec"), namespace)
        np.testing.assert_allclose(seen[0], [0., .4, np.nan], equal_nan=True)
        np.testing.assert_allclose(seen[1], [.1, .2, .5])

    def test_phase6_reports_all_missing_prerequisites_before_qm(self):
        fn = function(6, "stage0_assets")
        with tempfile.TemporaryDirectory() as tmp, patch("os.getcwd", return_value=tmp):
            with patch.object(Path, "is_file", return_value=False):
                with self.assertRaisesRegex(FileNotFoundError,
                                            "reactant_3d.mol.*product_3d.mol.*images_idpp.xyz"):
                    fn()

    def test_phase7_child_uses_own_dlls_and_no_parent_pythonpath(self):
        sub = SimpleNamespace(run=Mock(return_value=SimpleNamespace(returncode=0)))
        fn = function(7, "run_cmd", subprocess=sub)
        python = ROOT / "test_env" / "python.exe"
        with patch.dict(os.environ, {"PYTHONPATH": "wrong-packages", "PATH": "parent-bin"}):
            fn([python, "worker.py"], timeout=12)
        opts = sub.run.call_args.kwargs
        self.assertNotIn("PYTHONPATH", opts["env"])
        self.assertEqual(opts["env"]["PATH"].split(os.pathsep)[:2],
                         [str(python.parent / "Library" / "bin"), str(python.parent)])
        self.assertEqual(opts["timeout"], 12)

    def test_phase8_workers_keep_thread_budget_and_logs(self):
        sub = SimpleNamespace(run=Mock(return_value=SimpleNamespace(
            returncode=0, stdout="qc output", stderr="qc warning")))
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            py = out / "qc" / "python.exe"
            fn = function(8, "run_worker", subprocess=sub, ENV_PY_QC=str(py), OUT=out)
            with patch.dict(os.environ, {"PYTHONPATH": "wrong-packages"}):
                fn("8b", ["--res", "result.json"], threads=2)
            cmd = sub.run.call_args.args[0]
            env = sub.run.call_args.kwargs["env"]
            self.assertEqual(cmd[-2:], ["--threads", "2"])
            self.assertNotIn("--smoke", cmd)
            self.assertNotIn("PYTHONPATH", env)
            for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
                self.assertEqual(env[key], "2")
            log = (out / "worker_8b.log").read_text()
            self.assertIn("qc output", log)
            self.assertIn("qc warning", log)

    def test_phase8_worker_failure_preserves_diagnostics(self):
        sub = SimpleNamespace(run=Mock(return_value=SimpleNamespace(
            returncode=7, stdout="last successful point", stderr="native abort")))
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            fn = function(8, "run_worker", subprocess=sub, ENV_PY_QC=str(out / "python.exe"), OUT=out)
            with self.assertRaisesRegex(RuntimeError, "rc=7"):
                fn("8b", [])
            self.assertIn("last successful point", (out / "worker_8b.log").read_text())

    def test_phase8_dispatch_passes_threads_to_every_worker(self):
        dispatch = Mock(return_value={})
        for name, expected in (("stage_8a", 2), ("stage_8b", 1)):
            dispatch.reset_mock()
            fn = function(8, name, OUT=ROOT, ENGINE_NOTE="test", run_worker=dispatch,
                          parse_json=lambda _: {}, _warn=lambda *a: None)
            fn(SimpleNamespace(smoke=False, threads=2))
            self.assertEqual(dispatch.call_count, expected)
            for call in dispatch.call_args_list:
                self.assertEqual(call.kwargs["threads"], 2)


if __name__ == "__main__":
    unittest.main()
