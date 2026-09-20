"""Small deterministic checks; not full grid/time convergence tests."""
import unittest
import json
import tempfile
from pathlib import Path
from unittest.mock import patch
import numpy as np
import run_phase14_active_matter_condensate_phase_separation as m


class Phase14Tests(unittest.TestCase):
    def test_e1_e2_default_plan_and_reuse(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / 'round'
            store = m.ExperimentCheckpoints(root, m.checkpoint_config(1.))
            real_completed = m.completed_sim
            requested = []
            def cheap_completed(cp, key, kwargs, snaps):
                requested.append((key, kwargs.copy(), snaps))
                return real_completed(cp, key, dict(kwargs, n=8, dt=.003, t_end=.009), (.009,))
            with patch.object(m, 'completed_sim', side_effect=cheap_completed), patch.object(m, 'log'):
                original = m.experiment_main(1., store)
                scan = m.experiment_phase_diagram(1., store)
                resumed = m.ExperimentCheckpoints(root, m.checkpoint_config(1.), resume=True)
                with patch.object(m.ActiveCondensateSim, 'run', side_effect=AssertionError('recomputed')):
                    restored = m.experiment_main(1., resumed)
                    restored_scan = m.experiment_phase_diagram(1., resumed)
            self.assertEqual(len(scan), 28)
            self.assertEqual(scan, restored_scan)
            self.assertEqual(original['active'].frames, restored['active'].frames)
            for key, kwargs, snaps in requested:
                self.assertEqual(kwargs['t_end'], 1000.)
                self.assertEqual(kwargs['n'], 160 if key.startswith('E1') else 112)
                self.assertEqual(kwargs['dt'], .25 if key.startswith('E1') else .5)
                self.assertEqual(snaps, m.SNAP_TIMES if key.startswith('E1') else (1000.,))

    def test_completed_sim_roundtrip_all_frames_and_arrays(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / 'round'
            config = m.checkpoint_config(1.)
            store = m.ExperimentCheckpoints(root, config)
            kwargs = dict(n=8, dt=.003, t_end=.017, k_atp=.02)
            original = m.completed_sim(store, 'E1_active', kwargs, (0., .007, .017))
            before = {p.relative_to(root): p.read_bytes() for p in root.rglob('*') if p.is_file()}
            resumed = m.ExperimentCheckpoints(root, config, resume=True)
            with patch.object(m.ActiveCondensateSim, 'run', side_effect=AssertionError('recomputed')):
                restored = m.completed_sim(resumed, 'E1_active', kwargs, (0., .007, .017))
            self.assertEqual(original.frames, restored.frames)
            self.assertEqual(original.ness_stats(), restored.ness_stats())
            for key, val in vars(original).items():
                if isinstance(val, np.ndarray):
                    np.testing.assert_array_equal(val, getattr(restored, key))
                elif key == 'snapshots':
                    self.assertEqual(val.keys(), restored.snapshots.keys())
                    for t in val:
                        np.testing.assert_array_equal(val[t], restored.snapshots[t])
                else:
                    self.assertEqual(val, getattr(restored, key))
            self.assertEqual(before, {p.relative_to(root): p.read_bytes() for p in root.rglob('*') if p.is_file()})

    def test_identity_explicit_resume_and_immutable_commits(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / 'round'
            store = m.ExperimentCheckpoints(root, {'quick': False})
            store.save('E2_test', {'n': 8}, {'value': np.arange(3)})
            with self.assertRaises(FileExistsError):
                m.ExperimentCheckpoints(root, {'quick': False})
            with self.assertRaisesRegex(ValueError, 'explicit'):
                store.load('E2_test', {'n': 8})
            with self.assertRaises(FileExistsError):
                store.save('E2_test', {'n': 8}, {})
            with self.assertRaisesRegex(ValueError, 'mismatch'):
                m.ExperimentCheckpoints(root, {'quick': True}, resume=True)
            resumed = m.ExperimentCheckpoints(root, {'quick': False}, resume=True)
            with self.assertRaisesRegex(ValueError, 'identity'):
                resumed.load('E2_test', {'n': 9})
            resumed.round_id = 'a-different-round'
            with self.assertRaisesRegex(ValueError, 'identity'):
                resumed.load('E2_test', {'n': 8})
            with patch.object(m, '__file__', __file__):
                with self.assertRaisesRegex(ValueError, 'mismatch'):
                    m.ExperimentCheckpoints(root, {'quick': False}, resume=True)
                with self.assertRaisesRegex(ValueError, 'source changed'):
                    store.save('E2_other', {}, {})

    def test_interrupted_commit_can_retry_without_overwriting(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / 'round'
            store = m.ExperimentCheckpoints(root, {})
            with patch.object(m.os, 'rename', side_effect=OSError('interrupted')):
                with self.assertRaises(OSError):
                    store.save('E1_active', {}, np.ones((2, 2)))
            staged = list(root.glob('.E1_active-*'))
            self.assertEqual(len(staged), 1)
            resumed = m.ExperimentCheckpoints(root, {}, resume=True)
            self.assertIsNone(resumed.load('E1_active', {}))
            resumed.save('E1_active', {}, np.zeros((2, 2)))
            np.testing.assert_array_equal(resumed.load('E1_active', {}), np.zeros((2, 2)))
            self.assertTrue(staged[0].exists())

    def test_corruption_missing_commit_and_object_arrays_fail_closed(self):
        for filename in ('payload.json', 'arrays.npz', 'commit.json'):
            with self.subTest(filename=filename), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp) / 'round'
                store = m.ExperimentCheckpoints(root, {})
                store.save('E2_test', {}, np.arange(3))
                (root / 'E2_test' / filename).write_bytes(b'corrupt')
                resumed = m.ExperimentCheckpoints(root, {}, resume=True)
                with self.assertRaises(ValueError):
                    resumed.load('E2_test', {})
        with tempfile.TemporaryDirectory() as tmp:
            store = m.ExperimentCheckpoints(Path(tmp) / 'round', {})
            with self.assertRaisesRegex(ValueError, 'object arrays'):
                store.save('bad', {}, np.array([object()], dtype=object))
            store.save('safe', {}, np.arange(2))
            root = store.root
            # Even with matching checksum, an object NPZ must not unpickle.
            with open(root / 'safe' / 'arrays.npz', 'wb') as stream:
                np.savez(stream, a0=np.array([{}], dtype=object))
            manifest_path = root / 'safe' / 'commit.json'
            manifest = json.loads(manifest_path.read_text())
            manifest['npz_sha256'] = m._sha((root / 'safe' / 'arrays.npz').read_bytes())
            manifest_path.write_text(json.dumps(manifest))
            with self.assertRaises(ValueError):
                m.ExperimentCheckpoints(root, {}, resume=True).load('safe', {})
            manifest_path.unlink()
            with self.assertRaises(FileNotFoundError):
                m.ExperimentCheckpoints(root, {}, resume=True).load('safe', {})

    def test_e3_resume_skips_completed_fingerprints_after_interruption(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / 'round'
            store = m.ExperimentCheckpoints(root, {})
            real_completed = m.completed_sim
            def cheap_completed(cp, key, kwargs, snaps):
                self.assertEqual(kwargs['n'], 128)
                self.assertEqual(kwargs['t_end'], 1000.)
                return real_completed(cp, key, dict(kwargs, n=8, dt=.003, t_end=.009), (.009,))
            real_frap = m.simulate_frap
            def frap(sim, k_atp, **kwargs):
                if k_atp == m.K_ATP_NOMINAL:
                    raise RuntimeError('simulated interruption')
                return real_frap(sim, k_atp, t_sim=.5)
            with patch.object(m, 'completed_sim', side_effect=cheap_completed), patch.object(m, 'simulate_frap', side_effect=frap):
                with self.assertRaisesRegex(RuntimeError, 'simulated interruption'):
                    m.experiment_fingerprints({}, 1., store)
            resumed = m.ExperimentCheckpoints(root, {}, resume=True)
            saved = resumed.load('E3_k0.000', dict(k_atp=0., s_quick=1.))
            self.assertIsInstance(saved['_frap_curve'], tuple)
            sim_run = m.ActiveCondensateSim.run
            calls = []
            def counted_run(sim, *a, **kw):
                calls.append(sim.k_atp)
                return sim_run(sim, *a, **kw)
            with patch.object(m, 'completed_sim', side_effect=cheap_completed), patch.object(m.ActiveCondensateSim, 'run', counted_run), patch.object(m, 'simulate_frap', side_effect=lambda sim, k, **kw: real_frap(sim, k, t_sim=.5)) as compute:
                result = m.experiment_fingerprints({}, 1., resumed)
            self.assertEqual(calls, [.08])
            self.assertEqual(compute.call_count, 2)
            for key in ('_frap_curve', '_saxs_curve'):
                for a, b in zip(saved[key], result['k0.000'][key]):
                    np.testing.assert_array_equal(a, b)

    def test_invalid_state_and_step_fail(self):
        sim = m.ActiveCondensateSim(k_atp=.02, n=8)
        for dt in (0, -1, np.nan, np.inf):
            with self.assertRaises(ValueError):
                sim._attempt(dt)
        sim.phi[0, 0] = np.nan
        with self.assertRaises(ValueError):
            sim._attempt(.001)

    def test_reaction_and_step_conserve_total(self):
        sim = m.ActiveCondensateSim(k_atp=.02, n=16)
        before = np.mean(sim.phi + sim.psi)
        self.assertTrue(sim._attempt(.001))
        self.assertAlmostEqual(np.mean(sim.phi + sim.psi), before, places=13)
        self.assertTrue(np.isfinite(sim.phi).all())
        self.assertTrue(np.isfinite(sim.psi).all())

    def test_curvature_matches_derivative(self):
        sim = m.ActiveCondensateSim(k_atp=.02, n=8)
        for p in (.05, .3, .7, .95):
            h = 1e-6
            fd = (sim.fprime(np.array(p+h))-sim.fprime(np.array(p-h)))/(2*h)
            self.assertAlmostEqual(float(fd), sim.fpp(p), places=6)

    def test_frames_describe_reached_concentrations(self):
        sim = m.ActiveCondensateSim(k_atp=.02, n=8)
        self.assertTrue(sim._attempt(.001))
        self.assertAlmostEqual(sim.frames[-1]['phi_mean'], sim.phi.mean(), places=14)
        self.assertAlmostEqual(sim.frames[-1]['psi_mean'], sim.psi.mean(), places=14)
        self.assertEqual(sim.frames[-1]['dt'], .001)

    def test_ness_weights_physical_time_not_frame_count(self):
        sim = m.ActiveCondensateSim(k_atp=.02, n=8)
        keys = ('S_diff', 'S_chem', 'S_total', 'cycle_flux', 'R_mean_um',
                'n_droplets', 'area_fraction')
        sim.frames = [dict(t=t, **dict.fromkeys(keys, val))
                      for t, val in ((.1, 10.), (.2, 10.), (1., 2.))]
        self.assertAlmostEqual(sim.ness_stats(1.)['R_mean_um'], 3.6)
        self.assertAlmostEqual(sim.ness_stats(.5)['R_mean_um'], 2.)
        for fraction in (0., 2., np.nan):
            with self.assertRaises(ValueError): sim.ness_stats(fraction)

    def test_frap_failure_never_drops_diffusion(self):
        sim = m.ActiveCondensateSim(k_atp=.02, n=8)
        with patch.object(m, 'cg', return_value=(np.full(64, np.nan), 1)) as solver:
            with self.assertRaisesRegex(RuntimeError, 'FRAP diffusion solve failed'):
                m.simulate_frap(sim, .02, t_sim=.01)
            self.assertEqual(solver.call_count, 2)

    def test_time_boundary_and_post_step_snapshots(self):
        sim = m.ActiveCondensateSim(k_atp=.02, n=8, dt=.25, t_end=.031)
        # Deterministic surrogate step makes a pre-step timestamp detectable.
        def step():
            sim.t += sim.dt
            sim.phi.fill(sim.t)
        with patch.object(sim, 'step', side_effect=step):
            sim.run(snapshot_times=(0., .007, .019, .031, 2.))
        self.assertAlmostEqual(sim.t, .031, places=13)
        for stamp in (.007, .019, .031):
            np.testing.assert_allclose(sim.snapshots[stamp], stamp, atol=1e-13)
        self.assertNotIn(2., sim.snapshots)


if __name__ == '__main__':
    unittest.main()
