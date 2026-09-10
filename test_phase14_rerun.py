"""Small deterministic checks; not full grid/time convergence tests."""
import unittest
from unittest.mock import patch
import numpy as np
import run_phase14_active_matter_condensate_phase_separation as m


class Phase14Tests(unittest.TestCase):
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
