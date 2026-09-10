"""Cheap integration probes; these are not production reruns."""
import ast
import math
import os
import json
import time
import tempfile
from types import SimpleNamespace
from pathlib import Path
import unittest
import numpy as np
import torch
from torch import nn, optim

ROOT = Path(__file__).resolve().parent


def definitions(phase, names, **extra):
    path = next(ROOT.glob(f"run_phase{phase}_*.py"))
    tree = ast.parse(path.read_text(encoding="utf-8"))
    nodes = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))
             and n.name in names]
    assert {n.name for n in nodes} == set(names)
    ns = dict(np=np, torch=torch, nn=nn, optim=optim, math=math, os=os,
              Path=Path, SEED=11, F=nn.functional, elapsed=lambda: "test")
    ns.update(extra)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), "exec"), ns)
    return ns


class RerunTests(unittest.TestCase):
    def test_plant_timestep_respects_short_remaining_interval(self):
        tree = ast.parse((ROOT / "run_phase10_cyberphysical_flow_twin.py").read_text(encoding="utf-8"))
        node = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "advance")
        ns = {"np": np}
        exec(compile(ast.Module(body=[node], type_ignores=[]), "<advance>", "exec"), ns)
        steps = []
        plant = SimpleNamespace(t=0., dt_cap=.5, dz=1., u=1., u_cool=1., _k_eff=1000.,
                                T=np.ones(2), Tc=np.ones(2), C={"A":np.ones(2)},
                                _apply_actuators=lambda a: None, _thermal_fields=lambda: None,
                                _step_chemical=steps.append)
        ns["advance"](plant, .0004, np.zeros(5))
        self.assertAlmostEqual(plant.t, .0004)
        self.assertLessEqual(max(steps), .0003)
        plant.dt_cap = 0.
        with self.assertRaises(FloatingPointError):
            ns["advance"](plant, 1., np.zeros(5))

    def setUp(self):
        torch.set_num_threads(2)
        self.dtype = torch.get_default_dtype()

    def tearDown(self):
        torch.set_default_dtype(self.dtype)

    def test_vmc_ad_and_parameter_gradient(self):
        torch.set_default_dtype(torch.float64)
        ns = definitions(11, ["FermiPauliNet", "antisymmetry_test", "laplacian_fd_check"])
        self.assertTrue(ns["antisymmetry_test"]("cpu")[0])
        self.assertTrue(ns["laplacian_fd_check"]("cpu")[0])
        net = ns["FermiPauliNet"]([2], [[0., 0., 0.]], 1, 1,
                                  n_det=2, hidden=8, h2_dim=4, depth=1)
        r = torch.randn(8, 2, 3)
        lap, gsq = net.kinetic_terms(r)
        self.assertTrue(torch.isfinite(lap + gsq).all())
        energy = -0.5 * (lap + gsq)
        loss = ((energy - energy.mean()).detach() * net.logabs(r)).mean()
        loss.backward()
        grads = [p.grad for p in net.parameters() if p.grad is not None]
        self.assertTrue(grads)
        self.assertTrue(all(torch.isfinite(g).all() for g in grads))

    def test_scalar_training_both_modes(self):
        ns = definitions(12, ["train_scalar_functional", "nn_value_and_grad", "scalar_audit_metrics"])
        X = np.random.default_rng(12).normal(size=(16, 3))
        for mode in ["hamiltonian", "lyapunov"]:
            net = ns["train_scalar_functional"](lambda x: -x, mode, X, X * .1,
                      n_steps=2, transverse_samples=X * 2, pairs=(X, X * .9))
            value, grad = ns["nn_value_and_grad"](net, X)
            self.assertTrue(np.isfinite(value).all() and np.isfinite(grad).all())
            audit = ns["scalar_audit_metrics"](value, np.sum(grad * -X, axis=1))
            self.assertFalse(audit["global_certificate"])

    def test_sac_update(self):
        torch.set_default_dtype(torch.float32)
        ns = definitions(10, ["_mlp", "SquashedGaussianActor", "TwinQCritic",
                             "ReplayBuffer", "SACAgent"],
                         LOG_STD_LO=-20., LOG_STD_HI=2.,
                         nominal_bias_targets=lambda: np.zeros(5))
        agent = ns["SACAgent"](hidden=8)
        for _ in range(8):
            agent.buffer.add(np.ones(18), np.zeros(5), 1., np.zeros(18), 0.)
        for _ in range(2):
            self.assertTrue(np.isfinite(agent.update(batch=4)).all())

    def test_bath_charge_order_and_translation(self):
        fn = definitions(13, ["bath_observables"])["bath_observables"]
        positions = np.array([[1.,0.,0.], [1.1,0.,0.], [1.,.1,0.],
                              [0.,0.,1.], [.1,0.,1.], [0.,.1,1.]])
        probe = np.array([0., 0., .4])
        a, d = probe + [0.,0.,.035], probe - [0.,0.,.035]
        q = np.array([-.834,.417,.417,-.834,.417,.417])
        expected = np.sum(q / np.linalg.norm(positions-a,axis=1)
                          - q / np.linalg.norm(positions-d,axis=1))
        result = fn(positions, probe, a, d)
        self.assertAlmostEqual(result[0], expected)
        shift = np.array([3., -2., 5.])
        np.testing.assert_allclose(result, fn(positions+shift,probe+shift,a+shift,d+shift), atol=1e-13)

    def test_backend_path(self):
        fn = definitions(13, ["backend_environment"])["backend_environment"]
        env = fn(Path("C:/env/phase7/python.exe"))
        self.assertTrue(env["PATH"].startswith(str(Path("C:/env/phase7/Library/bin"))))
        self.assertNotIn("PYTHONPATH", env)

    def test_partial_reference_failure_is_explicit(self):
        with tempfile.TemporaryDirectory() as folder:
            def energy(geometry, basis, method, log):
                if geometry == "bad":
                    raise RuntimeError("synthetic reference failure")
                return -1.1
            lit = {"H2_diss": {"hf": None, "ccsdt": -1., "fci": -1.,
                                "source": "explicit fallback"}}
            ns = definitions(11, ["compute_references"], RESULTS=Path(folder),
                             PHASE7_PY=Path(__file__), LITERATURE=lit,
                             _psi4_one=energy, json=json, time=time)
            systems = [SimpleNamespace(name="good", tag="H2_diss", geom_xyz=lambda: "good"),
                       SimpleNamespace(name="bad", tag="H2_diss", geom_xyz=lambda: "bad")]
            result = ns["compute_references"](systems)
            self.assertEqual(set(result["refs"]), {"good", "bad"})
            self.assertFalse(result["complete"])
            self.assertEqual(result["computed_systems"], ["good"])
            self.assertEqual(result["fallback_systems"], ["bad"])


if __name__ == "__main__":
    unittest.main()
