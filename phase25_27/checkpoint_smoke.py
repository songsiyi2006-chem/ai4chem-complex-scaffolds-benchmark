"""Tiny real Psi4 cache/replay/integrity regression, not a chemistry benchmark."""
import json
from pathlib import Path
import tempfile
import numpy as np
import psi4
from psi4.driver.task_base import AtomicComputer
from .atomic_checkpoints import atomic_checkpoints


def main():
    psi4.set_memory('300 MB');psi4.set_num_threads(1)
    with tempfile.TemporaryDirectory() as td:
        out=Path(td);psi4.core.IOManager.shared_object().set_default_path(td)
        psi4.core.set_output_file(str(out/'psi4.out'),False)
        mol=psi4.geometry('0 1\nO 0 0 0\nH 0 0 1\nH 1 0 0\nsymmetry c1\nno_reorient\nno_com')
        def task():return AtomicComputer(molecule=mol,driver='gradient',method='hf',basis='sto-3g',keywords={'scf_type':'df'})
        with atomic_checkpoints(out/'jobs'):
            a=task();a.compute();b=task();b.compute()
            np.testing.assert_allclose(a.result.return_result,b.result.return_result,atol=0,rtol=0)
            events=json.loads((out/'jobs/progress.json').read_text())
            assert len(events)==2 and not events[0]['cached'] and events[1]['cached']
            binding=next((out/'jobs').glob('*.binding.json'));saved=binding.read_bytes()
            meta=json.loads(saved);meta['result_sha256']='invalid';binding.write_text(json.dumps(meta))
            try:task().compute()
            except ValueError as exc:assert 'integrity mismatch' in str(exc)
            else:raise AssertionError('Corrupt checkpoint was accepted')
            binding.write_bytes(saved)
        psi4.core.close_outfile();psi4.core.clean()
    print('PASS native gradient computation, identical replay, corrupt-checkpoint rejection')


if __name__=='__main__':main()
