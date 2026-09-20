"""Bounded math/ASE serialization tests; no QM, ML or optimization."""
import ast
import json
import math
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
from diagnose_phase4_neb import (atom_fmax, band_geometry, improved_tangent,
                                norm, project_force, replay_force_diagnostics)


class DiagnosticMathTests(unittest.TestCase):
    def test_monotonic_tangents(self):
        self.assertEqual(improved_tangent([2, 0, 0], [0, 3, 0], [0, 1, 2]), [0, 1, 0])
        self.assertEqual(improved_tangent([2, 0, 0], [0, 3, 0], [2, 1, 0]), [1, 0, 0])

    def test_extremum_weights(self):
        t = improved_tangent([1, 0, 0], [0, 1, 0], [0, 3, 2])
        self.assertAlmostEqual(t[1] / t[0], 3)
        self.assertAlmostEqual(norm(t), 1)

    def test_projection_and_spring(self):
        result = project_force([3, 4, 0], [1, 0, 0], spring_parallel=2)
        self.assertEqual(result['projected'], [2, 4, 0])
        self.assertEqual(result['perpendicular'], [0, 4, 0])

    def test_climbing_reflects_and_ignores_spring(self):
        r = project_force([3, 4, 0], [1, 0, 0], spring_parallel=99, climb=True)
        self.assertEqual(r['projected'], [-3, 4, 0])
        self.assertEqual(r['fmax'], 5)

    def test_global_tangent_not_per_atom_normalized(self):
        t = [1 / math.sqrt(2), 0, 0] * 2
        r = project_force([1, 0, 0, 0, 1, 0], t, climb=True)
        self.assertAlmostEqual(norm(r['projected']), math.sqrt(2))
        self.assertAlmostEqual(r['fmax'], math.sqrt(2))

    def test_fmax_is_max_atom_norm_not_component_or_band_norm(self):
        self.assertEqual(atom_fmax([3, 4, 0, 0, 0, 4]), 5)

    def test_path_and_spring_bound(self):
        g = band_geometry([[0, 0, 0], [1, 0, 0], [3, 0, 0], [4, 0, 0]], [0, 1, 3, 0])
        self.assertEqual(g['climbing_image'], 2)
        self.assertEqual(g['segment_lengths_A'], [1, 2, 1])
        self.assertAlmostEqual(g['images'][0]['projected_fmax_lower_bound_eV_A'], .1)
        self.assertEqual(g['images'][1]['projected_fmax_lower_bound_eV_A'], 0)

    def test_flat_tangent_rejected(self):
        with self.assertRaises(ValueError):
            improved_tangent([1, 0, 0], [1, 0, 0], [1, 1, 1])

    def test_nonfinite_and_shapes_rejected(self):
        for force in ([math.nan, 0, 0], [1, 2], []):
            with self.assertRaises(ValueError):
                atom_fmax(force)
        with self.assertRaises(ValueError):
            project_force([1, 0, 0], [2, 0, 0])
        with self.assertRaises(ValueError):
            band_geometry([[0, 0, 0]] * 3, [0, 1])

    def test_tied_climber_not_silently_selected(self):
        g = band_geometry([[i, 0, 0] for i in range(4)], [0, 1, 1, 0])
        self.assertIsNone(g['climbing_image'])
        self.assertEqual(g['climbing_candidates'], [1, 2])


class ProductionEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import numpy as np
        cls.np = np
        path = Path(__file__).with_name('run_phase4_reaction_mechanism.py')
        tree = ast.parse(path.read_text(encoding='utf-8'))
        names = {'_export_neb_candidate', '_save_neb_force_diagnostics'}
        selected = [node for node in tree.body if isinstance(node, ast.FunctionDef)
                    and node.name in names]
        assert len(selected) == 2
        cls.namespace = dict(np=np, json=json, Path=Path, __file__=str(path))
        exec(compile(ast.Module(body=selected, type_ignores=[]), str(path), 'exec'), cls.namespace)
        cls.tree = tree

    def test_export_without_calculator_results(self):
        from ase import Atoms
        from ase.io import read
        image = Atoms('H2', positions=[[0, 0, 0], [0, 0, .74]])
        image.calc = SimpleNamespace(name='no-results-calculator')
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'candidate.xyz'
            self.namespace['_export_neb_candidate'](path, image)
            loaded = read(path)
            self.assertGreater(path.stat().st_size, 0)
            self.np.testing.assert_array_equal(loaded.positions, image.positions)
            self.assertFalse(hasattr(image.calc, 'results'))

    def test_same_state_force_snapshot_and_replay_no_evaluations(self):
        from ase import Atoms
        from ase.calculators.singlepoint import SinglePointCalculator
        from ase.constraints import FixAtoms
        from ase.mep.neb import NEB
        class HalfForce:
            # Non-idempotent: a second adjustment would halve forces again.
            def adjust_forces(self, atoms, forces):
                forces *= .5

            def todict(self):
                return {'name': 'synthetic_HalfForce', 'kwargs': {}}

        images = []
        for x, energy in zip((0, 1, 3, 4), (0, 1, 3, 0)):
            image = Atoms('H2', positions=[[x, 0, 0], [x, 0, 1]])
            image.set_constraint([FixAtoms(indices=[1]), HalfForce()])
            image.calc = SinglePointCalculator(image, energy=energy,
                                               forces=[[.2, .3, 0], [.1, 0, .4]])
            images.append(image)
        neb = NEB(images, method='improvedtangent', climb=True, k=.1)
        projected = neb.get_forces()
        raw_before = neb.real_forces.copy()
        with tempfile.TemporaryDirectory() as directory:
            # Capturing evidence may not request even a cached calculator evaluation.
            with (patch.object(Atoms, 'get_forces', side_effect=AssertionError('extra force call')),
                  patch.object(Atoms, 'get_potential_energy', side_effect=AssertionError('extra energy call'))):
                self.namespace['_save_neb_force_diagnostics'](
                    Path(directory), neb, projected,
                    SimpleNamespace(name='synthetic fixed data', charge=0, mult=1), .05)
            saved = json.loads((Path(directory) / 'neb_force_diagnostics.json').read_text())
        self.assertEqual(saved['energies_eV'], [0, 1, 3, 0])
        self.assertEqual(saved['force_image_indices'], [1, 2])
        self.assertEqual(len(saved['raw_forces_eV_A']), 2)
        self.np.testing.assert_allclose(saved['raw_forces_eV_A'],
                                       [[[.2, .3, 0], [.1, 0, .4]]] * 2)
        self.np.testing.assert_allclose(saved['constraint_adjusted_forces_eV_A'],
                                       [[[.1, .15, 0], [0, 0, 0]]] * 2)
        self.assertIn('apply_constraint=False', saved['force_provenance']['raw'])
        self.assertEqual(saved['force_provenance']['extra_calculator_evaluations'], 0)
        self.np.testing.assert_array_equal(neb.real_forces, raw_before)
        self.np.testing.assert_array_equal(saved['projected_forces_eV_A'], projected.reshape(2, 2, 3))
        replay = replay_force_diagnostics(saved)
        self.assertLess(replay['projection_max_atom_error_eV_A'], 1e-12)
        self.assertEqual(replay['fmax_bookkeeping_error_eV_A'], 0)
        self.assertGreater(replay['constraint_adjustment_fmax_eV_A'], .4)
        self.assertFalse(replay['force_gate_passes'])
        saved['projected_forces_eV_A'][0][0][1] += .1
        self.assertGreater(replay_force_diagnostics(saved)['projection_max_atom_error_eV_A'], .09)

    def test_both_production_paths_save_final_force_and_safe_export(self):
        for name in ('stage2_neb', 'stage_bounded_refinement'):
            fn = next(n for n in self.tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
            calls = [n for n in ast.walk(fn) if isinstance(n, ast.Call)]
            saves = [n for n in calls if isinstance(n.func, ast.Name)
                     and n.func.id == '_save_neb_force_diagnostics']
            exports = [n for n in calls if isinstance(n.func, ast.Name)
                       and n.func.id == '_export_neb_candidate']
            self.assertEqual(len(saves), 1)
            self.assertEqual(len(exports), 1)
            self.assertLess(saves[0].lineno, exports[0].lineno)


if __name__ == '__main__':
    unittest.main()
