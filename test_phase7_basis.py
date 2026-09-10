import ast
from pathlib import Path
import unittest


class BasisTests(unittest.TestCase):
    def filter(self, records, tags, ref):
        path = Path(__file__).with_name('run_phase7_strong_correlation_wall.py')
        node = next(n for n in ast.parse(path.read_text(encoding='utf-8')).body
                    if isinstance(n, ast.FunctionDef) and n.name == 'consistent_basis_records')
        ns = {}
        exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), 'exec'), ns)
        return ns[node.name](records, tags, ref)

    def test_mixed_basis_is_not_a_curve(self):
        records = {'a': {'basis': 'def2-svp', 'energy': -400},
                   'b': {'basis': '6-31g', 'energy': -399},
                   'c': {'basis': 'DEF2-SVP', 'energy': -400.1}}
        kept, audit = self.filter(records, ['a', 'b', 'c'], 'a')
        self.assertIsNone(kept['b'])
        self.assertEqual(audit['excluded_tags'], ['b'])
        self.assertEqual(records['b']['energy'], -399)
        self.assertEqual(kept['c']['energy'], -400.1)

    def test_unknown_basis_not_assumed(self):
        kept, audit = self.filter({'a': {'energy': 0}}, ['a'], 'a')
        self.assertIsNone(kept['a'])
        self.assertIsNone(audit['reference_basis'])


if __name__ == '__main__':
    unittest.main()
