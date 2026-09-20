"""No real launch: validate one-shot ordering, provenance and exit recording."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch
import run_prepared_phase3 as runner


class PreparedTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name) / 'source'
        self.root.mkdir()
        self.attempt = self.root.parent / 'phase1-19-rerun-20260910/phase03/prepared_test'
        self.attempt.mkdir(parents=True)
        (self.attempt / runner.SCRIPT).write_text('print("test")')
        (self.attempt / 'results_phase3').mkdir()
        (self.attempt / 'results_phase3/stage1a.json').write_text('{}')
        record = dict(status='prepared_not_started',
                      restart_source_sha256=runner.digest(self.attempt / runner.SCRIPT),
                      reused_file_sha256={'stage1a.json': runner.digest(self.attempt / 'results_phase3/stage1a.json')})
        (self.attempt / 'restart_provenance.json').write_text(json.dumps(record))
        self.queue = self.root / 'results_rerun_audit/queue_status.json'
        self.queue.parent.mkdir()
        self.queue.write_text('{"status":"queue_finished"}')
        patcher = patch.object(runner, 'ROOT', self.root)
        patcher.start()
        self.addCleanup(patcher.stop)

    def execute(self, code=0):
        process = Mock(pid=1234)
        process.wait.return_value = code
        with patch.object(sys, 'argv', ['test', str(self.attempt), '--after-queue', str(self.queue)]), \
                patch.object(runner, 'available_gib', return_value=4), \
                patch.object(runner, 'runtime_metadata', return_value={}), \
                patch.object(runner.subprocess, 'Popen', return_value=process) as spawn:
            return runner.main(), spawn

    def test_default_full_command_and_completion(self):
        code, spawn = self.execute()
        self.assertEqual(code, 0)
        self.assertEqual(spawn.call_args.args[0][-1], runner.SCRIPT)
        self.assertNotIn('--skip_md', spawn.call_args.args[0])
        result = json.loads((self.attempt / 'rerun_status.json').read_text())
        self.assertEqual(result['status'], 'process_completed')
        self.assertEqual(result['scientific_acceptance'], 'pending output review')
        with self.assertRaises(FileExistsError):
            self.execute()

    def test_failed_child_stays_failed(self):
        self.assertEqual(self.execute(7)[0], 7)
        result = json.loads((self.attempt / 'rerun_status.json').read_text())
        self.assertEqual(result['status'], 'process_failed')

    def test_mutated_source_or_input_refused(self):
        runner.validate(self.attempt)
        (self.attempt / 'results_phase3/stage1a.json').write_text('{"changed":true}')
        with self.assertRaisesRegex(ValueError, 'identity mismatch'):
            runner.validate(self.attempt)
        self.assertFalse((self.attempt / 'rerun_status.json').exists())

    def test_failed_queue_does_not_launch(self):
        self.queue.write_text('{"status":"queue_failed"}')
        with patch.object(runner.subprocess, 'Popen') as spawn:
            with self.assertRaisesRegex(RuntimeError, 'did not finish'):
                self.execute()
            spawn.assert_not_called()
        self.assertFalse((self.attempt / 'rerun_status.json').exists())
        state = json.loads((self.attempt / 'restart_launch_status.json').read_text())
        self.assertEqual(state['status'], 'launcher_failed')


if __name__ == '__main__':
    unittest.main()
