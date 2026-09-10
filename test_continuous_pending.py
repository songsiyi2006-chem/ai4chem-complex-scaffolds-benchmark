import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import run_continuous_pending as q


class QueueTests(unittest.TestCase):
    def test_finite_queue_preserves_failures_and_continues(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)/'repo'
            root.mkdir()
            results = [dict(phase=16, exit_code=1), dict(phase=12, exit_code=0)]
            with patch.object(q, 'ROOT', root), patch.object(q, 'available_gib', return_value=5), \
                 patch.object(q, 'run_one', side_effect=results) as run, \
                 patch('sys.argv', ['queue', '--phases', '16', '12']):
                q.main()
            status = json.loads(next(root.glob('results_rerun_audit/*/queue_status.json')).read_text())
            self.assertEqual(status['status'], 'queue_finished')
            self.assertEqual(status['completed'], results)
            self.assertEqual(run.call_count, 2)
            self.assertIn('pending', status['scientific_acceptance'])

    def test_unresolved_attempt_blocks_duplicate_run(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)/'repo'
            root.mkdir()
            prior = root.parent/'phase1-19-rerun-20260910/phase16/old'
            prior.mkdir(parents=True)
            (prior/'rerun_status.json').write_text('{"status":"running"}')
            with patch.object(q, 'ROOT', root), patch.object(q, 'run_one') as run, \
                 patch('sys.argv', ['queue', '--phases', '16']):
                with self.assertRaises(RuntimeError):
                    q.main()
                run.assert_not_called()


if __name__ == '__main__':
    unittest.main()
