"""Isolated regressions of real AST bodies; no legacy import/job side effects.

Run: python -m unittest -v test_phase6_10_audit
Only NumPy is required. These tests do not validate full scientific workflows.
"""
import ast
import datetime
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parent
FILES = {i: next(ROOT.glob(f'run_phase{i}_*.py')) for i in range(6, 11)}


def tree(i):
    return ast.parse(FILES[i].read_text(encoding='utf-8'))


def execute(nodes, ns):
    exec(compile(ast.fix_missing_locations(ast.Module(body=nodes, type_ignores=[])),
                 '<real-source-extract>', 'exec'), ns)
    return ns


def method(i, name, ns):
    node = next(n for n in ast.walk(tree(i))
                if isinstance(n, ast.FunctionDef) and n.name == name)
    return execute([node], ns)[name]


class AuditTests(unittest.TestCase):
    def test_sources_compile(self):
        for i in FILES:
            compile(tree(i), str(FILES[i]), 'exec')

    def test_phase6_bias_conversion(self):
        f = next(n for n in tree(6).body if isinstance(n, ast.FunctionDef)
                 and n.name == 'analyze_fes')
        scale = next(n for n in f.body if isinstance(n, ast.Assign)
                     and any(isinstance(t, ast.Name) and t.id == 'scale' for t in n.targets))
        for gamma in (2., 10., 50.):
            ns = execute([scale], {'GAMMA': gamma})
            free_energy = np.array([0., 4., 10.])
            bias = -(gamma-1)/gamma * free_energy
            np.testing.assert_allclose(ns['scale'] * bias, free_energy)

    def test_phase7_uhf_uses_unrestricted_data(self):
        assignment = next(n for n in ast.walk(tree(7)) if isinstance(n, ast.Assign)
                          and any(isinstance(t, ast.Name) and t.id == 'rel' for t in n.targets)
                          and isinstance(n.value, ast.Dict))
        ns = execute([assignment], dict(to_rel=lambda x: x, e_cas=[0.],
                                       e_rhf=[100.], e_uhf=[2.], e_uks=[3.]))
        self.assertEqual(ns['rel']['UHF'], [2.])

    def test_phase8_force_matches_energy_gradient(self):
        force = method(8, 'forces_and_nac', {'np': np})
        eigs = method(8, 'eigs', {'np': np})
        m = SimpleNamespace(sqrt_It=1., wR=0., wB=0., kR=0., kB=0.,
                            lR=0., lB=0., rel_q_lo=-10., rel_q_hi=10.)
        m.phi_of = lambda q: q[:, 0]
        m.D0 = m.DREL = lambda x: np.zeros_like(x)
        m.D1 = lambda x: np.ones_like(x)
        # Nonzero coupling slope checks both terms of the square-root derivative.
        m.DV = lambda x: np.full_like(x, .07)
        m.diabats = lambda q: (np.zeros(len(q)), 2.+q[:, 0], .3+.07*q[:, 0])
        m.eigs = lambda q: eigs(m, q)
        q = np.array([[0., 0., 0.], [.5, 0., 0.], [-1., 0., 0.]])
        for active in (0, 1):
            analytic = force(m, q, np.full(len(q), active))[0]
            for axis in range(3):
                qp, qm = q.copy(), q.copy()
                qp[:, axis] += 1e-6
                qm[:, axis] -= 1e-6
                def energy(x):
                    e = m.eigs(x)
                    return e[3] + (2*active-1)*e[5]
                numeric = -(energy(qp)-energy(qm))/2e-6
                np.testing.assert_allclose(analytic[:, axis], numeric, atol=1e-9)

    def test_phase9_generated_tip_isolation(self):
        t = tree(9)
        ns = dict(Path=Path, json=json, _dt=datetime)
        names = {'PHASE5_ANCHOR', 'REAGENT_WELLS', 'RESERVOIR_WELLS',
                 'SURROGATE', 'LIQUID_CLASSES'}
        for n in t.body:
            targets = n.targets if isinstance(n, ast.Assign) else ([n.target] if isinstance(n, ast.AnnAssign) else [])
            if any(isinstance(x, ast.Name) and x.id in names for x in targets):
                execute([n], ns)
        execute([n for n in t.body if isinstance(n, (ast.ClassDef, ast.FunctionDef))
                 and n.name in {'OT2Emitter', '_j', 'compile_ot2_protocol'}], ns)
        condition = dict(T_c=20., cat_molpct=2., t_h=.1, phi_tol=.5,
                         label='test', well='A1', vials=['A1', 'A2'], v_mix_uL=80.)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'protocol.py'
            ns['compile_ot2_protocol']([condition], 'audit', 'audit', path)
            generated = ast.parse(path.read_text(encoding='utf-8'))
        run = next(n for n in generated.body if isinstance(n, ast.FunctionDef) and n.name == 'run')
        # Every subsequent stock _dose in a straight-line block needs a new tip.
        dirty = False
        for n in run.body:
            if not isinstance(n, ast.Expr) or not isinstance(n.value, ast.Call):
                continue
            call = n.value
            if isinstance(call.func, ast.Attribute) and ast.unparse(call.func) == 'p20.pick_up_tip':
                dirty = False
            if isinstance(call.func, ast.Name) and call.func.id == '_dose' and ast.unparse(call.args[0]) == 'p20':
                self.assertFalse(dirty, 'stock accessed with destination-exposed tip')
                dirty = True
        dose = next(n for n in generated.body if isinstance(n, ast.FunctionDef) and n.name == '_dose')
        fn = execute([dose], {'AIR_GAP_UL': 10})['_dose']
        class Well:
            def __init__(self, name): self.name = name
            def bottom(self, *a): return self
            def top(self, *a): return self
        class Pipette:
            max_volume = 20.
            flow_rate = SimpleNamespace()
            dirty = False
            tips = 1
            delivered = 0.
            def aspirate(self, v, w):
                if w.name == 'stock':
                    assert not self.dirty, 'split transfer contaminated stock'
            def dispense(self, v, w):
                if w.name == 'dest':
                    self.dirty = True
                    self.delivered += v-10
            def air_gap(self, v): pass
            def blow_out(self, w): pass
            def touch_tip(self, w): self.dirty = True
            def drop_tip(self): pass
            def pick_up_tip(self):
                self.dirty = False
                self.tips += 1
        p = Pipette()
        fn(p, Well('stock'), Well('dest'), 25., 1., 1., prewet=2)
        self.assertEqual(p.tips, 3)
        self.assertEqual(p.delivered, 25.)

    def test_phase10_carbon_units(self):
        reward = method(10, '_reward', dict(np=np, C_A_FEED=1., M_P=1000.,
                         C_CAT_FEED=1., D_TUBE=1., V_REACTOR=1.,
                         W1=1., W2=1., W3=1., W4=1., W5=1.))
        q = 1e-6
        plant = SimpleNamespace(q_tot=q, Q={'A': 1., 'cat': 0.},
                                pbpr=1000*.45/q/1e5, T=np.array([0.]), Tc=0.)
        plant._rates = lambda: [np.array([1.])]*5
        plant._thermal_fields = lambda: (0, 0, np.array([0.]))
        obs = dict(truth_out={'P': (1/3600)/(q*1e3)}, dP=0., max_dT=0.)
        _, info = reward(SimpleNamespace(plant=plant), obs, None)
        self.assertAlmostEqual(info['carbon'], .5542)  # 1 kW / 1 kg per hour
        plant.Q['cat'] = 1.
        _, more = reward(SimpleNamespace(plant=plant), obs, None)
        expected = (1/6e7)*1e3*.05*640/1e3*85/(1/3600)
        self.assertAlmostEqual(more['carbon']-info['carbon'], expected)

    def test_phase10_observation_preserves_sensor_sample(self):
        fn = method(10, '_get_obs', {'np': np})
        def forbidden(**kw): raise AssertionError('unexpected second sensor sample')
        env = SimpleNamespace(plant=SimpleNamespace(observe=forbidden),
                              _dp_prev=0., _act_prev=np.zeros(5))
        obs = dict(est=dict(A=0., I=0., P=0.), ports_T=[293.15]*3,
                   dP=0., actions=np.zeros(5), max_dT=1., dT_rate=.4, vapor_margin=10.)
        vector = fn(env, obs)
        self.assertEqual(vector[14], 2.)
        self.assertIs(env._last_obs_full, obs)
        step = next(n for n in ast.walk(tree(10)) if isinstance(n, ast.FunctionDef)
                    and n.name == 'step' and 'obs_raw' in ast.unparse(n))
        self.assertIn('self._get_obs(obs_raw)', ast.unparse(step))


if __name__ == '__main__':
    unittest.main()
