"""Compatibility tests independent of scientific dependencies."""
import base64
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import run_phase


class LayoutTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / 'repo'
        self.root.mkdir()
        self.workspace = Path(self.temp.name) / 'run'
        self.files = {
            'driver.py': ('projects/phase01/code/driver.py',
                b'from pathlib import Path\nfrom helper import value\nfrom pkg import util\n'
                b'assert value == util.answer == 42\n'
                b'assert (Path(__file__).parent / "results/data.txt").read_text() == "original"\n'
                b'Path("results/data.txt").write_text("changed")\n'),
            'helper.py': ('projects/phase02/code/helper.py', b'value = 42\n'),
            'pkg/__init__.py': ('shared/code/pkg/__init__.py', b''),
            'pkg/util.py': ('projects/phase03/code/util.py', b'from helper import value\nanswer = value\n'),
            'results/data.txt': ('projects/phase01/results/data.txt', b'original'),
            'REPORT.md': ('projects/phase01/reports/REPORT.md', b'[data](../results/data.txt)\n'),
        }
        self.manifest = {'version': 1, 'files': []}
        for old, (new, data) in self.files.items():
            path = self.root/new; path.parent.mkdir(parents=True, exist_ok=True); path.write_bytes(data)
            self.manifest['files'].append({'old': old, 'new': new, 'organized_sha256': run_phase.digest(data)})
        original = {'REPORT.md': base64.b64encode(b'[data](results/data.txt)\n').decode()}
        (self.root/'shared/original_document_bytes.json.gz').write_bytes(gzip.compress(json.dumps(original).encode()))

    def test_cross_phase_imports_relative_package_and_output_isolation(self):
        run_phase.materialize(self.root, self.workspace, self.manifest)
        result = subprocess.run([sys.executable, 'driver.py'], cwd=self.workspace, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.workspace/'results/data.txt').read_text(), 'changed')
        self.assertEqual((self.root/'projects/phase01/results/data.txt').read_text(), 'original')
        self.assertEqual((self.workspace/'REPORT.md').read_bytes(), b'[data](results/data.txt)\n')

    def test_existing_directory_is_never_silently_overwritten(self):
        self.workspace.mkdir()
        (self.workspace/'keep').write_text('keep')
        with self.assertRaises(ValueError):
            run_phase.materialize(self.root, self.workspace, self.manifest)
        self.assertEqual((self.workspace/'keep').read_text(), 'keep')

    def test_reuse_retains_outputs_and_rejects_source_changes(self):
        run_phase.materialize(self.root, self.workspace, self.manifest)
        (self.workspace/'results/data.txt').write_text('new result')
        run_phase.materialize(self.root, self.workspace, self.manifest, reuse=True)
        self.assertEqual((self.workspace/'results/data.txt').read_text(), 'new result')
        (self.root/'projects/phase02/code/helper.py').write_text('value = 99')
        with self.assertRaises(ValueError):
            run_phase.materialize(self.root, self.workspace, self.manifest, reuse=True)

    def test_reuse_rejects_modified_workspace_code(self):
        run_phase.materialize(self.root, self.workspace, self.manifest)
        (self.workspace/'driver.py').write_text('raise SystemExit(0)')
        with self.assertRaises(ValueError):
            run_phase.materialize(self.root, self.workspace, self.manifest, reuse=True)

    def test_rejects_traversal_and_repository_destination(self):
        for name in ('../outside', '/absolute', 'C:/absolute', 'C:relative', 'a\\b'):
            with self.subTest(name=name), self.assertRaises(ValueError):
                run_phase.safe_path(self.root, name)
        for dest in (self.root, self.root/'projects'/'new', self.root.parent):
            with self.subTest(dest=dest), self.assertRaises(ValueError):
                run_phase.materialize(self.root, dest, self.manifest)

    def test_later_source_edits_are_used_in_new_workspaces(self):
        (self.root/'projects/phase02/code/helper.py').write_text('value = 99')
        (self.root/'projects/phase01/reports/REPORT.md').write_text('Updated report')
        run_phase.materialize(self.root, self.workspace, self.manifest)
        self.assertEqual((self.workspace/'helper.py').read_text(), 'value = 99')
        self.assertEqual((self.workspace/'REPORT.md').read_text(), 'Updated report')


if __name__ == '__main__':
    unittest.main()
