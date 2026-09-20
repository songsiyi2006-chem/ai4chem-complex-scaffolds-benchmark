"""Test native phase dispatch without importing scientific packages or launching jobs."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import run_phase


class NativePhaseTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root/'docs').mkdir()
        entry = self.root/'projects/phase28/run.py'
        entry.parent.mkdir(parents=True)
        entry.write_text('raise SystemExit(0)\n')
        self.manifest = {'phases': {'1': {'title': 'legacy', 'folder': 'projects/phase01', 'entrypoint': 'legacy.py'}}, 'files': []}
        self.native = {'28': {'title': 'native', 'folder': 'projects/phase28', 'entrypoint': 'projects/phase28/run.py', 'execution': 'native'}}
        self.manifest_path = self.root/'docs/layout_manifest.json'
        self.native_path = self.root/'docs/native_phases.json'
        self.manifest_path.write_text(json.dumps(self.manifest))
        self.native_path.write_text(json.dumps(self.native))
        self.addCleanup(patch.stopall)
        patch.object(run_phase, 'ROOT', self.root).start()
        patch.object(run_phase, 'MANIFEST', self.manifest_path).start()
        patch.object(run_phase, 'NATIVE_PHASES', self.native_path).start()

    def test_dispatch_help_without_materializing_legacy_files(self):
        with patch.object(run_phase.subprocess, 'call', return_value=0) as call, patch.object(run_phase, 'materialize') as materialize:
            self.assertEqual(run_phase.main(['28','--','--help']), 0)
            materialize.assert_not_called()
            self.assertEqual(call.call_args.args[0][-1], '--help')
            self.assertEqual(call.call_args.kwargs['cwd'], self.root)

    def test_workspace_forwarded_as_native_output_and_return_code_propagated(self):
        out = self.root/'work/new-results'
        with patch.object(run_phase.subprocess, 'call', return_value=7) as call:
            self.assertEqual(run_phase.main(['28','--workspace',str(out),'--python','custom-python','--','--self-test']), 7)
            args = call.call_args.args[0]
            self.assertEqual(args[0], 'custom-python')
            self.assertEqual(args[2:], ['--out',str(out.resolve()),'--self-test'])

    def test_rejects_ambiguous_output_and_legacy_flags(self):
        for args in (['28','--workspace','x','--','--out=y'], ['28','--workspace','x','--','--output-dir=y'], ['28','--prepare-only'], ['28','--module','phase25_27.analysis']):
            with self.subTest(args=args), self.assertRaises(SystemExit) as cm, patch.object(run_phase.subprocess, 'call') as call:
                run_phase.main(args)
            self.assertEqual(cm.exception.code, 2)
            call.assert_not_called()

    def test_cannot_shadow_old_phase_or_escape_repository(self):
        for native in ({'1':self.native['28']}, {'28':dict(self.native['28'],entrypoint='../outside.py')}):
            self.native_path.write_text(json.dumps(native))
            with self.assertRaises(ValueError):run_phase.phase_registry(self.manifest)


if __name__ == '__main__':
    unittest.main()
