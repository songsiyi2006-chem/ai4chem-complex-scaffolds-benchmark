"""Stage4 completion regressions using fake trajectories; never import/run MD."""
import ast
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch


SOURCE = Path(__file__).with_name("run_phase3_complex_dynamics.py")


class Stage4CompletionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name)
        self.args = SimpleNamespace(md_steps=100000, report_interval=500,
                                    equil_steps=5000)
        names = {"_stage4_expected_frames", "_stage4_validate_cache",
                 "stage4_complex_md"}
        tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
        nodes = [node for node in tree.body
                 if isinstance(node, ast.FunctionDef) and node.name in names]
        self.assertEqual(len(nodes), len(names))
        self.replay = Mock()
        self.ns = dict(Path=Path, json=json, RESULTS={}, _log=Mock(),
                       _replay_stage4_analysis=self.replay)
        exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SOURCE), "exec"),
             self.ns)
        self.dcd = self.out / "T04_complex_trajectory.dcd"
        self.dcd.write_bytes(b"incomplete trajectory: must not overwrite")
        (self.out / "complex_start.pdb").write_text("preserve topology")
        self.load = Mock()
        # Block all molecular engine imports even if a regression falls through.
        self.modules = patch.dict(sys.modules, {
            "mdtraj": SimpleNamespace(load=self.load), "openmm": None,
            "openmm.app": None, "rdkit": None,
        })
        self.modules.start()
        self.addCleanup(self.modules.stop)
        self.frames(200)

    def frames(self, count):
        self.load.return_value = SimpleNamespace(xyz=SimpleNamespace(shape=(count, 1, 3)))

    def cache(self, **overrides):
        rec = dict(production_steps=100000, equil_steps=5000,
                   report_interval=500, dt_fs=1.0, production_ps=100.0,
                   n_frames=200, md={"n_frames": 200})
        rec.update(overrides)
        (self.out / "stage4.json").write_text(json.dumps(rec), encoding="utf-8")
        return rec

    def run_stage(self, force=False):
        return self.ns["stage4_complex_md"](self.out, self.args, force)

    def assert_preserved_failure(self, pattern, force=False):
        before = {p.name: p.read_bytes() for p in self.out.iterdir()}
        with self.assertRaisesRegex((RuntimeError, ValueError), pattern):
            self.run_stage(force)
        self.assertEqual(before, {p.name: p.read_bytes() for p in self.out.iterdir()})
        self.replay.assert_not_called()
        self.assertNotIn("stage4_complex_md", self.ns["RESULTS"])

    def test_default_requires_200_frames(self):
        self.assertEqual(self.ns["_stage4_expected_frames"](self.args), 200)

    def test_incomplete_or_excess_dcd_is_preserved(self):
        for count in (0, 99, 100, 199, 201):
            with self.subTest(frames=count):
                self.frames(count)
                self.assert_preserved_failure("expected 200")

    def test_short_dcd_cannot_be_hidden_by_complete_cache(self):
        self.cache()
        self.frames(100)
        self.assert_preserved_failure("expected 200")

    def test_valid_cache_skips_replay(self):
        rec = self.cache()
        self.assertTrue(self.run_stage())
        self.assertEqual(self.ns["RESULTS"]["stage4_complex_md"], rec)
        self.load.assert_called_once()
        self.replay.assert_not_called()

    def test_inconsistent_cache_is_preserved(self):
        for overrides in (dict(production_steps=50000), dict(equil_steps=1000),
                          dict(report_interval=1000), dict(dt_fs=2.0),
                          dict(production_ps=50.0), dict(n_frames=100),
                          dict(md={"n_frames": 100}), dict(md=None)):
            with self.subTest(overrides=overrides):
                self.cache(**overrides)
                self.assert_preserved_failure("metadata mismatch")

    def test_legacy_cache_missing_interval_is_not_trusted(self):
        rec = self.cache()
        del rec["report_interval"]
        (self.out / "stage4.json").write_text(json.dumps(rec))
        self.assert_preserved_failure("report_interval")

    def test_malformed_cache_is_preserved(self):
        for value in ("{broken", "null", "[]"):
            with self.subTest(value=value):
                (self.out / "stage4.json").write_text(value)
                self.assert_preserved_failure(".")

    def test_unreadable_trajectory_is_preserved(self):
        self.load.side_effect = ValueError("truncated DCD")
        self.assert_preserved_failure("unreadable")

    def test_orphan_checkpoint_is_preserved(self):
        self.cache()
        self.dcd.unlink()
        self.assert_preserved_failure("no trajectory")

    def test_force_does_not_overwrite_existing_trajectory(self):
        self.frames(100)
        self.assert_preserved_failure("preserved even with --force_rerun", force=True)
        self.load.assert_not_called()

    def test_complete_dcd_recovery_writes_matching_metadata(self):
        summary = dict(n_frames=200, ca_rmsd_mean_A=1.0, lig_rmsd_mean_A=2.0)
        analyzer = SimpleNamespace(summary=Mock(return_value=summary), rows=[], plif={})
        self.replay.return_value = dict(analyzer=analyzer, n_frames=200)
        original = self.dcd.read_bytes()
        self.assertTrue(self.run_stage())
        self.replay.assert_called_once_with(self.out, self.args)
        rec = json.loads((self.out / "stage4.json").read_text())
        self.assertEqual(rec["production_steps"], 100000)
        self.assertEqual(rec["report_interval"], 500)
        self.assertEqual(rec["n_frames"], 200)
        self.assertEqual(original, self.dcd.read_bytes())
        self.replay.reset_mock()
        self.assertTrue(self.run_stage())
        self.replay.assert_not_called()

    def test_custom_steps_and_absolute_report_alignment(self):
        # Heating+equil ends at 9250; reports at 9500, 10000, 10500.
        self.args.equil_steps = 5250
        self.args.md_steps = 1250
        self.assertEqual(self.ns["_stage4_expected_frames"](self.args), 3)
        self.frames(2)
        self.assert_preserved_failure("expected 3")
        self.args.md_steps = 50000
        self.assertEqual(self.ns["_stage4_expected_frames"](self.args), 100)

    def test_invalid_configuration_is_rejected_without_writes(self):
        for key, value in (("md_steps", 0), ("md_steps", -1),
                           ("report_interval", 0), ("report_interval", -1),
                           ("equil_steps", -1), ("md_steps", 1)):
            with self.subTest(key=key, value=value):
                old = getattr(self.args, key)
                setattr(self.args, key, value)
                self.assert_preserved_failure(".")
                setattr(self.args, key, old)


if __name__ == "__main__":
    unittest.main()
