"""Single-thread, tiny Reference contexts only; no campaign MD or ML imports."""
import ast
import copy
import hashlib
import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import uuid


SOURCE = Path(__file__).with_name("run_phase3_complex_dynamics.py")


class CheckpointTests(unittest.TestCase):
    def setUp(self):
        import openmm as mm
        from openmm import app, unit
        self.mm, self.app, self.unit = mm, app, unit
        self.tmp = tempfile.TemporaryDirectory(prefix="phase3_checkpoint_test_")
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name)
        self.root = self.out / "stage4_checkpoints"
        self.root.mkdir()
        names = {"_md_file_identity", "_md_write_durable", "_md_identity",
                 "_MDCheckpointWriter", "_md_validate_checkpoint",
                 "_stage4_expected_frames", "stage4_complex_md"}
        tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
        nodes = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))
                 and n.name in names]
        self.assertEqual(len(nodes), len(names))
        self.ns = dict(Path=Path, hashlib=hashlib, json=json, os=os, uuid=uuid)
        exec(compile(ast.Module(nodes, type_ignores=[]), str(SOURCE), "exec"), self.ns)
        system = mm.System()
        system.addParticle(12)
        force = mm.CustomExternalForce("0.5*(x*x+y*y+z*z)")
        force.addParticle(0, [])
        system.addForce(force)
        topology = app.Topology()
        residue = topology.addResidue("ARG", topology.addChain())
        topology.addAtom("CA", app.element.carbon, residue)
        integ = mm.LangevinMiddleIntegrator(310*unit.kelvin, 1/unit.picosecond,
                                          1*unit.femtosecond)
        integ.setRandomNumberSeed(17)
        self.sim = app.Simulation(topology, system, integ, mm.Platform.getPlatformByName("Reference"))
        self.sim.context.setPositions([[.1, .2, .3]])
        self.sim.context.setVelocitiesToTemperature(310*unit.kelvin, 23)
        self.sim.step(3)
        self.dcd = app.DCDReporter(str(self.out / "T04_complex_trajectory.dcd"), 2)
        self.addCleanup(self.dcd._out.close)
        self.sim.reporters.append(self.dcd)
        self.analyzer = SimpleNamespace(frames=0, rows=[], plif={})
        self.args = SimpleNamespace(md_steps=6, equil_steps=5, report_interval=2,
                                    mgb_frames=40, pdb_id="toy", lig_code="toy")
        self.identity = self.ns["_md_identity"](self.sim, self.args,
                                               {"source": self.ns["_md_file_identity"](SOURCE)})
        self.writer = self.ns["_MDCheckpointWriter"](self.root, self.identity, 3, 6, 2)

    def advance(self):
        self.sim.step(2 - self.sim.currentStep % 2)
        state = self.sim.context.getState(getEnergy=True)
        self.analyzer.rows.append([0., 0., 0., 310., float(
            state.getPotentialEnergy().value_in_unit(self.unit.kilojoule_per_mole))])
        self.analyzer.frames += 1
        self.analyzer.plif[("hbond", 0, "A:ARG1")] = self.analyzer.frames

    def save_first(self):
        self.advance()
        self.writer.save(self.sim, self.dcd, self.analyzer)

    def validate(self, identity=None):
        return self.ns["_md_validate_checkpoint"](self.root,
                                                  self.identity if identity is None else identity)

    def snapshot(self):
        return {str(p.relative_to(self.out)): p.read_bytes()
                for p in self.out.rglob("*") if p.is_file()}

    def test_binary_restores_exact_reference_rng_and_step(self):
        self.save_first()
        progress = self.validate()
        self.assertEqual(progress["absolute_step"], 4)
        self.assertEqual(progress["production_steps_completed"], 1)
        self.assertEqual(progress["frame_absolute_steps"], [4])
        self.assertEqual(progress["rows"], self.analyzer.rows)
        self.assertEqual(progress["plif_counts"], [[["hbond", 0, "A:ARG1"], 1]])
        manifest = json.loads((self.root / "manifest.json").read_text())
        binary = (self.root / manifest["generation"] / "state.chk").read_bytes()
        self.sim.reporters.clear()
        self.sim.step(2)
        def values(sim):
            st = sim.context.getState(getPositions=True, getVelocities=True)
            return (str(st.getPositions()), str(st.getVelocities()), st.getTime(), st.getStepCount())
        expected = values(self.sim)
        mm = self.mm
        restored = self.app.Simulation(self.sim.topology,
            mm.XmlSerializer.deserialize(self.identity["system_xml"]),
            mm.XmlSerializer.deserialize(self.identity["integrator_xml"]),
            mm.Platform.getPlatformByName("Reference"))
        restored.context.loadCheckpoint(binary)  # diagnostic toy only, no application resume API
        self.assertEqual(restored.currentStep, 4)
        restored.step(2)
        self.assertEqual(values(restored), expected)

    def test_identity_mismatches_rejected_without_mutation(self):
        self.save_first()
        before = self.snapshot()
        for field in ("inputs", "config", "platform", "platform_properties", "host",
                      "openmm_version", "system_xml", "integrator_xml", "integrator_rng_seed"):
            with self.subTest(field=field):
                changed = copy.deepcopy(self.identity)
                changed[field] = "different"
                with self.assertRaisesRegex(RuntimeError, "identity mismatch"):
                    self.validate(changed)
        self.assertEqual(before, self.snapshot())

    def test_uncommitted_dcd_advance_is_rejected(self):
        self.save_first()
        self.advance()
        self.dcd._out.flush()
        with self.assertRaisesRegex(RuntimeError, "DCD changed"):
            self.validate()

    def test_manifest_replace_failure_preserves_prior_commit(self):
        self.save_first()
        previous = (self.root / "manifest.json").read_bytes()
        self.advance()
        with patch.object(os, "replace", side_effect=OSError("injected replace failure")):
            with self.assertRaisesRegex(OSError, "injected"):
                self.writer.save(self.sim, self.dcd, self.analyzer)
        self.assertEqual((self.root / "manifest.json").read_bytes(), previous)
        self.assertEqual(self.writer.last_step, 4)
        with self.assertRaisesRegex(RuntimeError, "DCD changed"):
            self.validate()

    def test_data_write_failure_cannot_publish_manifest(self):
        self.advance()
        original = self.ns["_md_write_durable"]
        def fault(path, data):
            if path.name == "progress.json":
                raise OSError("injected disk full")
            original(path, data)
        with patch.dict(self.ns, _md_write_durable=fault):
            with self.assertRaisesRegex(OSError, "disk full"):
                self.writer.save(self.sim, self.dcd, self.analyzer)
        self.assertFalse((self.root / "manifest.json").exists())

    def test_corrupt_binary_and_alignment_fail_closed(self):
        self.save_first()
        manifest = json.loads((self.root / "manifest.json").read_text())
        (self.root / manifest["generation"] / "state.chk").write_bytes(b"truncated")
        with self.assertRaisesRegex(RuntimeError, "integrity mismatch"):
            self.validate()
        self.advance()
        self.analyzer.frames -= 1
        with self.assertRaisesRegex(RuntimeError, "alignment mismatch"):
            self.writer.save(self.sim, self.dcd, self.analyzer)

    def test_flush_failure_and_nonfinite_metrics_never_publish(self):
        self.advance()
        with patch.object(os, "fsync", side_effect=OSError("injected flush failure")):
            with self.assertRaisesRegex(OSError, "flush failure"):
                self.writer.save(self.sim, self.dcd, self.analyzer)
        self.assertFalse((self.root / "manifest.json").exists())
        self.analyzer.rows[0][0] = float("nan")
        with self.assertRaises(ValueError):
            self.writer.save(self.sim, self.dcd, self.analyzer)
        self.assertFalse((self.root / "manifest.json").exists())

    def test_resume_and_force_preserve_partial_without_engine_work(self):
        for resume, force in ((True, False), (False, True), (False, False)):
            with self.subTest(resume=resume, force=force):
                # Remove only our empty toy DCD so partial-root handling is exercised.
                self.dcd._out.close()
                (self.out / "T04_complex_trajectory.dcd").unlink(missing_ok=True)
                self.args.resume_md = resume
                before = self.snapshot()
                with self.assertRaisesRegex(RuntimeError, "disabled|preserved"):
                    self.ns["stage4_complex_md"](self.out, self.args, force)
                self.assertEqual(before, self.snapshot())

    def test_report_boundaries_and_final_nonreport_chunk(self):
        for _ in range(3):
            self.advance()
            self.writer.save(self.sim, self.dcd, self.analyzer)
        self.sim.step(1)  # final step 9: no additional frame
        self.writer.save(self.sim, self.dcd, self.analyzer)
        progress = self.validate()
        self.assertEqual(progress["production_steps_completed"], 6)
        self.assertEqual(progress["frame_absolute_steps"], [4, 6, 8])
        self.assertEqual(progress["frames"], 3)
        self.assertEqual([x["absolute_step"] for x in progress["chunk_boundaries"]], [4, 6, 8, 9])
        self.assertEqual(len(list(self.root.glob("step_*"))), 4)


if __name__ == "__main__":
    unittest.main()
