"""Exercise the real CLI finalization without running any chemistry."""
import argparse
import ast
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch


class DeliveryTests(unittest.TestCase):
    def run_main(self, plot_error=False, save_error=False):
        source = Path(__file__).with_name('run_phase3_complex_dynamics.py')
        node = next(n for n in ast.parse(source.read_text(encoding='utf-8')).body
                    if isinstance(n, ast.FunctionDef) and n.name == 'main')
        result = dict(timestamp='test', all_stages_ok=False, meta=dict(warnings=[]))
        ns = dict(argparse=argparse, Path=Path, sys=sys, RESULTS=result,
                  TARGET_PDB='7RPZ', NATIVE_LIG_CODE='test',
                  traceback=Mock(), _hr=Mock(), _log=Mock(), _warn=Mock(),
                  _shutdown=Mock(), stage6_figures=Mock(), write_json_atomic=Mock())
        for name in ('stage1a_ingest', 'stage2_docking', 'stage3_ff_duality',
                     'stage4_complex_md', 'stage5_mmgbsa'):
            ns[name] = Mock(return_value=True)
        if plot_error:
            ns['stage6_figures'].side_effect = RuntimeError('render failed')
        if save_error:
            ns['write_json_atomic'].side_effect = OSError('disk full')
        exec(compile(ast.Module([node], type_ignores=[]), str(source), 'exec'), ns)
        with tempfile.TemporaryDirectory() as folder, patch.object(sys, 'argv',
                ['test', '--out_dir', folder, '--fig_dir', folder]):
            code = ns['main']()
        ns['_shutdown'].assert_not_called()
        return code, result

    def test_success(self):
        code, result = self.run_main()
        self.assertEqual(code, 0)
        self.assertTrue(result['all_stages_ok'])

    def test_plot_failure_is_not_success(self):
        code, result = self.run_main(plot_error=True)
        self.assertEqual(code, 1)
        self.assertFalse(result['all_stages_ok'])

    def test_serialization_failure_is_not_success(self):
        code, result = self.run_main(save_error=True)
        self.assertEqual(code, 1)
        self.assertFalse(result['all_stages_ok'])


if __name__ == '__main__':
    unittest.main()
