"""Persist native Psi4/QCEngine atomic jobs without changing the potential surface.

The exact QCSchema input is hashed. Only successful matching AtomicResults can
be reused. This wrapper is local to this process and does not edit installed code.
"""
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import time
import numpy as np


@contextmanager
def atomic_checkpoints(directory):
    from psi4.driver.task_base import AtomicComputer
    from qcelemental.models.v2 import AtomicResult
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    original=AtomicComputer.compute
    progress=directory/'progress.json';events=json.loads(progress.read_text(encoding='utf8')) if progress.exists() else []

    def wrapped(self,client=None):
        if client is not None:raise ValueError('This pilot checkpoints local jobs only')
        if self.computed:return
        plan=json.loads(self.plan().model_dump_json())
        canonical=json.dumps(plan,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
        key=hashlib.sha256(canonical).hexdigest();resultfile=directory/(key+'.result.json')
        inputfile=directory/(key+'.input.json');cached=False;t=time.monotonic()
        inputfile.write_bytes(canonical)
        if resultfile.exists():
            result=AtomicResult.model_validate_json(resultfile.read_text(encoding='utf8'))
            if not result.success:raise ValueError('Failed checkpoint cannot be reused')
            # The filename binds the result to its exact submitted QCSchema input.
            binding=directory/(key+'.binding.json')
            expected=json.loads(binding.read_text(encoding='utf8'))
            if expected['input_sha256']!=key or hashlib.sha256(resultfile.read_bytes()).hexdigest()!=expected['result_sha256']:
                raise ValueError('Checkpoint integrity mismatch')
            actual=json.loads(result.input_data.model_dump_json())
            for field in ['model','driver','keywords']:
                if actual['specification'][field]!=plan['specification'][field]:raise ValueError('Checkpoint model/options mismatch')
            for field in ['symbols','molecular_charge','molecular_multiplicity']:
                if actual['molecule'][field]!=plan['molecule'][field]:raise ValueError('Checkpoint molecular identity mismatch')
            if not np.array_equal(actual['molecule']['geometry'],plan['molecule']['geometry']):raise ValueError('Checkpoint geometry mismatch')
            self.result=result;self.computed=True;cached=True
        else:
            original(self,client=None)
            if not self.result.success:raise ValueError('Atomic job failed')
            data=self.result.model_dump_json().encode();resultfile.write_bytes(data)
            (directory/(key+'.binding.json')).write_bytes(json.dumps(dict(input_sha256=key,result_sha256=hashlib.sha256(data).hexdigest()),indent=2).encode())
        events.append(dict(input_sha256=key,cached=cached,wall_seconds=time.monotonic()-t,successful=True))
        progress.write_bytes((json.dumps(events,indent=2)+'\n').encode())
        print('Native Hessian atomic job',len(events),'cached' if cached else 'computed',round(events[-1]['wall_seconds'],1),flush=True)

    AtomicComputer.compute=wrapped
    try:yield
    finally:AtomicComputer.compute=original
