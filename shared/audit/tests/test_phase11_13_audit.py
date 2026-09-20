"""Run with python -m unittest -v test_phase11_13_audit (NumPy/SciPy)."""
import ast
from pathlib import Path
import unittest
import numpy as np
from scipy import linalg

ROOT = Path(__file__).resolve().parent


def load(phase, names, ns):
    source = next(ROOT.glob(f'run_phase{phase}_*.py')).read_text(encoding='utf-8')
    t = ast.parse(source)
    nodes = [n for n in ast.walk(t) if isinstance(n, ast.FunctionDef) and n.name in names]
    assert set(n.name for n in nodes) == set(names)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), '<source>', 'exec'), ns)
    return ns


class Tests(unittest.TestCase):
    def test_compile(self):
        for phase in range(11, 14):
            p = next(ROOT.glob(f'run_phase{phase}_*.py'))
            compile(p.read_text(encoding='utf-8'), str(p), 'exec')

    def test_logabs_has_no_floor(self):
        ns = load(11, ['logabs'], {})
        class Net:
            def _fermi_layer(self, r): return np.array([-100., -1000.]), 1., np.array([1., 2.])
        np.testing.assert_array_equal(ns['logabs'](Net(), None), [-99., -998.])

    def test_scalar_audit(self):
        f = load(12, ['scalar_audit_metrics'], {'np': np})['scalar_audit_metrics']
        out = f([1., 1.000001], [0.0005, 0.0005])
        self.assertEqual(out['V_nonincrease_fraction'], 0.)
        self.assertAlmostEqual(out['H_absolute_drift'], 1e-6)
        self.assertAlmostEqual(out['H_drift_fixed_scale'], 4e-6)
        self.assertFalse(out['global_certificate'])
        self.assertEqual(f([1., 1.], [0., -1.])['V_nonincrease_fraction'], 1.)
        for H, dV in (([], [0.]), ([np.nan], [0.]), ([0.], [np.inf])):
            with self.assertRaises(ValueError): f(H, dV)

    def test_hamiltonian_transverse_constraint(self):
        p = next(ROOT.glob('run_phase12_*.py'))
        s = p.read_text(encoding='utf-8')
        self.assertNotIn('0.25 - Hb.std()', s)
        self.assertIn('0.25 - Hr.std(unbiased=False)', s)
        self.assertIn('nn_value_and_grad(netV, audit_off)', s)
        self.assertIn('pairs=pairs, transverse_samples=off_samples', s)

    def pcet(self):
        return load(13, ['proton_pes', 'fd_schrodinger', 'masked_well_states',
                         'vibronic_rate_hs', 'pcet_rate_geometry'],
                    dict(np=np, linalg=linalg, BOHR_A=.529177210903,
                         HARTREE_EV=27.211386245988, K_B_EV=8.617333262e-5))

    def test_identical_state_overlap(self):
        ns = self.pcet()
        for dx in (.02, .01):
            x = np.arange(-1., 1., dx)
            e, v = ns['fd_schrodinger'](x*x, x, 1.008)
            _, _, S = ns['vibronic_rate_hs'](x, e, v, e, v, 0., .8, 100.)
            np.testing.assert_allclose(S, np.eye(4), atol=1e-12)

    def test_pes_direction(self):
        f = self.pcet()['proton_pes']
        for dg in (-.12, 0., .12):
            v = f(np.array([1.02, 1.72]), dG_eV=dg)
            self.assertAlmostEqual(v[0]-v[1], dg)

    def test_pcet_rate_grid_convergence(self):
        f = self.pcet()['pcet_rate_geometry']
        coarse = np.array(f(.7, .4, -.12, .8, 100., 300., grid=.0012)[:2])
        fine = np.array(f(.7, .4, -.12, .8, 100., 300., grid=.0006)[:2])
        self.assertTrue(np.all(np.isfinite(fine)) and np.all(fine > 0))
        np.testing.assert_allclose(coarse/fine, np.ones(2), rtol=.1)


if __name__ == '__main__':
    unittest.main()
