"""CPU-only, five-step checkpoint tests. No scientific/full pipeline training."""
import copy
import random
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import torch

import run_phase12_hamiltonian_law_discovery as phase12


class InterruptedTraining(Exception):
    pass


class CheckpointTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.samples = np.random.default_rng(81).uniform(0.1, 1.0, (8, 3))

    def train(self, mode="hamiltonian", **kwargs):
        options = dict(n_steps=5, seed=17, pairs=(self.samples, self.samples * .9),
                       transverse_samples=self.samples * 1.2,
                       checkpoint_every=2, rhs_identity="-x with test RNG perturbation")
        options.update(kwargs)

        def rhs(x):
            # Exercise every persisted RNG in addition to Torch minibatch sampling.
            return -x + 1e-3 * (np.random.random() + random.random())

        return phase12.train_scalar_functional(rhs, mode, self.samples,
                                               self.samples * .8, **options)

    def assert_state_equal(self, a, b):
        if torch.is_tensor(a):
            self.assertTrue(torch.equal(a, b))
        elif isinstance(a, dict):
            self.assertEqual(a.keys(), b.keys())
            for key in a:
                self.assert_state_equal(a[key], b[key])
        elif isinstance(a, (tuple, list)):
            self.assertEqual(len(a), len(b))
            for x, y in zip(a, b):
                self.assert_state_equal(x, y)
        else:
            self.assertEqual(a, b)

    def test_interrupted_resume_equals_continuous_both_modes(self):
        for mode in ("hamiltonian", "lyapunov"):
            with self.subTest(mode=mode):
                full = self.root / (mode + "-full.pt")
                interrupted = self.root / (mode + "-interrupted.pt")
                random.seed(12)
                np.random.seed(13)
                self.train(mode, checkpoint_path=full)
                random.seed(12)
                np.random.seed(13)
                original_step = torch.optim.Adam.step
                steps = 0

                def stop_during_third_step(optimizer, *args, **kwargs):
                    nonlocal steps
                    result = original_step(optimizer, *args, **kwargs)
                    steps += 1
                    if steps == 3:  # weights updated; scheduler/save not yet run
                        raise InterruptedTraining()
                    return result

                with patch.object(torch.optim.Adam, "step", stop_during_third_step):
                    with self.assertRaises(InterruptedTraining):
                        self.train(mode, checkpoint_path=interrupted)
                self.assertEqual(torch.load(interrupted, weights_only=True)["step"], 2)
                random.seed(99)
                np.random.seed(98)
                torch.manual_seed(97)
                self.train(mode, checkpoint_path=interrupted, resume=True)
                expected = torch.load(full, weights_only=True)
                actual = torch.load(interrupted, weights_only=True)
                self.assert_state_equal(expected, actual)
                # Completed checkpoint must perform zero additional RHS evaluations.
                with patch.object(np.random, "random", side_effect=AssertionError("extra step")):
                    net = self.train(mode, checkpoint_path=interrupted, resume=True)
                self.assert_state_equal(expected["model"], net.state_dict())
                random.seed(12)
                np.random.seed(13)
                plain = self.train(mode)
                self.assert_state_equal(expected["model"], plain.state_dict())
                self.assertEqual(len(list(self.root.glob("*.tmp"))), 0)

    def test_mismatches_and_missing_states_rejected_without_overwrite(self):
        path = self.root / "checkpoint.pt"
        self.train(checkpoint_path=path)
        original = path.read_bytes()
        for change in ({"n_steps": 6}, {"lr": .001}, {"seed": 18},
                       {"rhs_identity": "different law"},
                       {"pairs": None}, {"transverse_samples": self.samples * 2}):
            with self.subTest(change=list(change)):
                with self.assertRaisesRegex(ValueError, "mismatch"):
                    self.train(checkpoint_path=path, resume=True, **change)
                self.assertEqual(path.read_bytes(), original)
        with patch.object(phase12, "SOURCE_SHA256", "changed"):
            with self.assertRaisesRegex(ValueError, "mismatch"):
                self.train(checkpoint_path=path, resume=True)
        original_samples = self.samples.copy()
        self.samples[0, 0] += .01
        with self.assertRaisesRegex(ValueError, "mismatch"):
            self.train(checkpoint_path=path, resume=True)
        self.samples = original_samples
        saved = torch.load(path, weights_only=True)
        for key, value in (("step", -1), ("step", 6), ("model", None),
                           ("optimizer", None), ("scheduler", None), ("torch_rng", None)):
            bad = copy.deepcopy(saved)
            if value is None:
                bad.pop(key)
            else:
                bad[key] = value
            with patch.object(torch, "load", return_value=bad):
                with self.assertRaises(ValueError):
                    self.train(checkpoint_path=path, resume=True)
        with self.assertRaises(FileExistsError):
            self.train(checkpoint_path=path)
        self.assertEqual(path.read_bytes(), original)

    def test_atomic_failure_preserves_old_checkpoint(self):
        path = self.root / "checkpoint.pt"
        phase12._atomic_checkpoint_write(path, lambda f: f.write(b"previous"))

        def fail(f):
            f.write(b"incomplete")
            raise OSError("simulated disk failure")

        with self.assertRaises(OSError):
            phase12._atomic_checkpoint_write(path, fail)
        self.assertEqual(path.read_bytes(), b"previous")
        with patch.object(phase12.os, "replace", side_effect=OSError("replace failed")):
            with self.assertRaises(OSError):
                phase12._atomic_checkpoint_write(path, lambda f: f.write(b"new"))
        self.assertEqual(path.read_bytes(), b"previous")
        self.assertEqual(list(self.root.iterdir()), [path])

    def test_explicit_resume_and_old_attempt_preservation(self):
        old = self.root / "old-attempt"
        old.mkdir()
        marker = old / "log.txt"
        marker.touch()
        with self.assertRaises(FileExistsError):
            phase12._prepare_checkpoint_dir(old, False, False)
        with self.assertRaisesRegex(ValueError, "without weights"):
            phase12._prepare_checkpoint_dir(old, True, False)
        self.assertTrue(marker.exists())
        fresh = phase12._prepare_checkpoint_dir(self.root / "new", False, False)
        with self.assertRaisesRegex(ValueError, "weights"):
            phase12._prepare_checkpoint_dir(fresh, True, False)
        self.train(checkpoint_path=fresh / "hamiltonian.pt")
        self.assertEqual(phase12._prepare_checkpoint_dir(fresh, True, False), fresh)
        with self.assertRaisesRegex(ValueError, "mismatch"):
            phase12._prepare_checkpoint_dir(fresh, True, True)
        with self.assertRaisesRegex(ValueError, "weights"):
            self.train(checkpoint_path=self.root / "missing.pt", resume=True)
        with self.assertRaisesRegex(ValueError, "resume requires"):
            self.train(resume=True)
        with self.assertRaisesRegex(ValueError, "positive"):
            self.train(checkpoint_every=0)


if __name__ == "__main__":
    unittest.main()
