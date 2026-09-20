import ast
from pathlib import Path
import unittest
import numpy as np


class FieldTests(unittest.TestCase):
    def field(self):
        source = Path(__file__).with_name('run_phase16_megamachine_cryoem_transport.py')
        tree = ast.parse(source.read_text(encoding='utf-8'))
        cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'TrilinearField')
        ns = dict(np=np)
        exec(compile(ast.Module(body=[cls], type_ignores=[]), str(source), 'exec'), ns)
        grid = (np.sin(np.linspace(0, 3, 30))[:, None, None]
                * np.cos(np.linspace(0, 2, 24))[None, :, None]
                * np.linspace(1, 2, 20)[None, None, :])
        return ns['TrilinearField'](grid, np.zeros(3), np.ones(3)*0.1)

    def test_interior_gradient_matches_finite_difference(self):
        field = self.field()
        p = np.random.default_rng(7).uniform(.12, (field.n-1)*field.voxel-.12, (40, 3))
        gradient, inside = field.gradients(p)
        self.assertTrue(inside.all())
        eps = 1e-5
        fd = np.stack([(field.values(p+np.eye(3)[a]*eps)[0]
                        - field.values(p-np.eye(3)[a]*eps)[0])/(2*eps)
                       for a in range(3)], axis=1)
        np.testing.assert_allclose(gradient, fd, atol=1e-8, rtol=1e-7)

    def test_outside_support_force_is_zero(self):
        field = self.field()
        gradient, inside = field.gradients(np.array([[.5, 2.8, .5], [.5, .5, 2.8], [-.1, .5, .5]]))
        self.assertFalse(inside.any())
        np.testing.assert_array_equal(gradient, np.zeros((3, 3)))


if __name__ == '__main__':
    unittest.main()
