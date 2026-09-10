"""Real OpenMM microchecks; never imports or launches the production entrypoint."""
import ast
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

import numpy as np
import openmm as mm
from openmm import app, unit

SOURCE = Path(__file__).with_name("run_phase15_quantum_biology_spin_allostery.py")
TREE = ast.parse(SOURCE.read_text(encoding="utf-8"))
NAMES = {"checkpoint_identity", "save_window_checkpoint", "load_window_checkpoint"}
NS = dict(Path=Path, json=json, os=os, sys=sys)
exec(compile(ast.Module(body=[n for n in TREE.body if isinstance(n, ast.FunctionDef)
                             and n.name in NAMES], type_ignores=[]), str(SOURCE), "exec"), NS)


def simulation(engine="Reference"):
    topology = app.Topology()
    residue = topology.addResidue("Ar", topology.addChain())
    topology.addAtom("Ar", app.element.argon, residue)
    system = mm.System()
    system.addParticle(40.)
    force = mm.CustomExternalForce("0.5*k*x*x")
    force.addGlobalParameter("k", 10.)
    force.addParticle(0, [])
    system.addForce(force)
    integrator = mm.LangevinMiddleIntegrator(300., 1., .002)
    integrator.setRandomNumberSeed(1515)
    sim = app.Simulation(topology, system, integrator,
                         mm.Platform.getPlatformByName(engine),
                         {"Threads": "1"} if engine == "CPU" else {})
    sim.context.setPositions([[.1, .2, .3]])
    sim.context.setVelocitiesToTemperature(300., 19)
    return sim


class CheckpointTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "window_checkpoint.zip"
        self.sim = simulation()
        self.identity = NS["checkpoint_identity"](self.sim, "FAD_oxid", {"steps": 20}, SOURCE)
        self.sim.step(10)
        self.progress = dict(phase="windows", executed_steps=10, requested_steps=20,
                             next_window=1, windows=[dict(samples_nm=[.21, .24])],
                             result=dict(latch_dist_nm=[.31, .32]))
        NS["save_window_checkpoint"](self.path, self.sim, self.identity, self.progress)

    def load(self, sim=None, identity=None, resume=True):
        return NS["load_window_checkpoint"](self.path, sim or self.sim,
                                            identity or self.identity, resume=resume)

    def test_real_roundtrip_continues_random_stream_and_samples(self):
        self.sim.context.setParameter("k", 17.)
        NS["save_window_checkpoint"](self.path, self.sim, self.identity, self.progress)
        self.sim.step(10)
        expected = self.sim.context.getState(getPositions=True, getVelocities=True)
        resumed = simulation()
        actual_identity = NS["checkpoint_identity"](resumed, "FAD_oxid", {"steps": 20}, SOURCE)
        self.assertEqual(actual_identity, self.identity)
        restored = self.load(resumed, actual_identity)
        self.assertEqual(restored, self.progress)
        self.assertEqual(resumed.currentStep, 10)
        self.assertEqual(resumed.context.getParameter("k"), 17.)
        resumed.step(10)
        actual = resumed.context.getState(getPositions=True, getVelocities=True)
        np.testing.assert_array_equal(expected.getPositions(asNumpy=True)._value,
                                      actual.getPositions(asNumpy=True)._value)
        np.testing.assert_array_equal(expected.getVelocities(asNumpy=True)._value,
                                      actual.getVelocities(asNumpy=True)._value)
        restored.update(phase="complete", executed_steps=20)
        NS["save_window_checkpoint"](self.path, resumed, actual_identity, restored)
        with zipfile.ZipFile(self.path) as archive:
            self.assertTrue(json.loads(archive.read("manifest.json"))["complete"])
        self.assertEqual(self.load(resumed)["executed_steps"], 20)

    def test_cpu_real_checkpoint(self):
        sim = simulation("CPU")
        identity = NS["checkpoint_identity"](sim, "cpu", {}, SOURCE)
        path = Path(self.temp.name) / "cpu.zip"
        sim.step(3)
        progress = dict(phase="complete", executed_steps=3, requested_steps=3)
        NS["save_window_checkpoint"](path, sim, identity, progress)
        other = simulation("CPU")
        other_identity = NS["checkpoint_identity"](other, "cpu", {}, SOURCE)
        NS["load_window_checkpoint"](path, other, other_identity, resume=True)
        self.assertEqual(other.currentStep, 3)
        other.step(2)
        self.assertEqual(other.currentStep, 5)

    def test_all_wrong_hashes_and_state_rejected_without_mutation(self):
        before = self.path.read_bytes()
        for key in ("source_hash", "config_hash", "system_hash", "integrator_hash", "version_hash", "state"):
            with self.subTest(key=key):
                bad = dict(self.identity, **{key: "wrong"})
                with self.assertRaisesRegex(ValueError, "Incompatible"):
                    self.load(identity=bad)
                with self.assertRaisesRegex(ValueError, "overwrite incompatible"):
                    NS["save_window_checkpoint"](self.path, self.sim, bad, self.progress)
                self.assertEqual(self.path.read_bytes(), before)

    def test_explicit_flag_required(self):
        with self.assertRaisesRegex(ValueError, "explicit"):
            self.load(resume=False)

    def test_atomic_replace_failure_keeps_prior_generation(self):
        before = self.path.read_bytes()
        self.sim.step(2)
        with patch("os.replace", side_effect=OSError("interrupted publication")):
            with self.assertRaises(OSError):
                NS["save_window_checkpoint"](self.path, self.sim, self.identity,
                                              dict(self.progress, executed_steps=12))
        self.assertEqual(self.path.read_bytes(), before)
        self.assertEqual(self.load()["executed_steps"], 10)

    def test_corrupt_samples_rejected(self):
        with zipfile.ZipFile(self.path) as archive:
            files = {name: archive.read(name) for name in archive.namelist()}
        files["samples.json"] = b'{}'
        with zipfile.ZipFile(self.path, "w") as archive:
            for name, contents in files.items():
                archive.writestr(name, contents)
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            self.load()

    def test_incomplete_cannot_be_marked_complete(self):
        with self.assertRaisesRegex(ValueError, "incomplete"):
            NS["save_window_checkpoint"](self.path, self.sim, self.identity,
                                          dict(self.progress, phase="complete"))

    def test_full_sampling_defaults_unchanged(self):
        config = next(n.value for n in TREE.body if isinstance(n, ast.Assign)
                      and any(isinstance(t, ast.Name) and t.id == "CONFIG" for t in n.targets))
        c = {k.arg: ast.literal_eval(k.value) for k in config.keywords}
        self.assertEqual((70 + c["MD_PRODUCTION_PS"] + c["UMB_N_WIN"] *
                          (c["UMB_SETTLE_PS"] + c["UMB_SAMPLE_PS"])) / (c["MD_DT_FS"] / 1000), 833000)


if __name__ == "__main__":
    unittest.main()
