"""Fail-closed tests for campaign checkpoint reuse; no expensive QC jobs."""
import copy
import tempfile
import unittest
from pathlib import Path
import resume_phase7_rerun as r


class ResumeTests(unittest.TestCase):
    def point(self):
        return dict(R=1.4, tag='R1.40', basis='def2-svp', theory={
            name: dict(theory=name, basis='def2-svp', converged=True,
                       energy_eh=-10., s2=2. if name.endswith('triplet') else 0., ds2=0.)
            for name in r.METHODS})

    def test_all_six_converged_methods(self):
        r.valid_point(self.point(), '7A', 1.4, 'def2-svp')
        for name in r.METHODS:
            data = self.point()
            del data['theory'][name]
            with self.assertRaises(ValueError):
                r.valid_point(data, '7A', 1.4, 'def2-svp')

    def test_nonfinite_and_unconverged_rejected(self):
        for change in ({'energy_eh': float('nan')}, {'converged': False},
                       {'basis': 'sto-3g'}, {'s2': 1., 'ds2': 0.}):
            data = self.point()
            data['theory']['uhf'].update(change)
            with self.assertRaises(ValueError):
                r.valid_point(data, '7A', 1.4, 'def2-svp')

    def test_path_escape_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                r.checked_path(Path(tmp), '../outside')

    def test_changed_cache_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root/'checkpoint.json'
            path.write_text('{}', encoding='utf-8')
            hashes = {path.name: r.digest(path)}
            r.verify_hashes(root, hashes)
            path.write_text('{"changed":true}', encoding='utf-8')
            with self.assertRaises(ValueError):
                r.verify_hashes(root, hashes)


if __name__ == '__main__':
    unittest.main()
