import unittest
import ast
from pathlib import Path
import re
import subprocess
import tempfile
from unittest.mock import patch
import numpy as np
from refine_phase4_geometric import atomic_units, EH_EV, BOHR_A


class AdapterTests(unittest.TestCase):
    def test_actual_hessian_parser_handles_utf8_headers(self):
        source = Path(__file__).with_name('run_phase4_reaction_mechanism.py')
        tree = ast.parse(source.read_text(encoding='utf-8'))
        nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef)
                 and n.name in {'xtb_hessian', '_parse_g98_modes'}]
        ns = dict(Path=Path, tempfile=tempfile, subprocess=subprocess, re=re,
                  np=np, XTB_EXE='xtb')
        exec(compile(ast.Module(body=nodes, type_ignores=[]), str(source), 'exec'), ns)
        def run(cmd, **kwargs):
            self.assertEqual(cmd[-2:], ['--uhf', '0'])
            self.assertEqual(kwargs['encoding'], 'utf-8')
            folder = Path(kwargs['cwd'])
            (folder/'vibspectrum').write_text('# frequency / cm⁻¹\n1 a -123.45 0 YES\n2 a 234.56 0 YES\n', encoding='utf-8')
            (folder/'g98.out').write_text('cm⁻¹\nFrequencies -- -123.45 234.56\nAtom AN\n1 1 1 0 0 0 1 0\n', encoding='utf-8')
            return subprocess.CompletedProcess(cmd, 0, stdout='Thermochemistry', stderr='')
        with patch.object(subprocess, 'run', run):
            result = ns['xtb_hessian']([1], [[0, 0, 0]])
        self.assertEqual(result['frequencies'], [-123.45, 234.56])
        self.assertEqual(len(result['modes']), 2)
        np.testing.assert_array_equal(result['modes'][1]['disp'], [[0, 1, 0]])

    def test_energy_and_negative_force_units(self):
        result = atomic_units(EH_EV, [[EH_EV / BOHR_A, 0, 0]])
        self.assertEqual(result['energy'], 1.)
        np.testing.assert_allclose(result['gradient'], [-1., 0, 0])

    def test_gradient_matches_analytic_energy_in_bohr(self):
        q = np.array([.1, .2, -.3])
        def energy(x):
            return float(((x * BOHR_A) ** 2).sum() / 2 / EH_EV)
        fd = [(energy(q + d) - energy(q - d)) / 2e-6 for d in np.eye(3) * 1e-6]
        np.testing.assert_allclose(atomic_units(0, -q * BOHR_A)['gradient'], fd, atol=1e-12)


if __name__ == '__main__':
    unittest.main()
