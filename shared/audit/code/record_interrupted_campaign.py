"""Record observed lost processes without inventing their exit codes."""
import argparse
import ctypes
import json
import os
from pathlib import Path
import time


def process_state(pid):
    if os.name != 'nt':
        try:
            os.kill(pid, 0)
            return True, None
        except ProcessLookupError:
            return False, None
    kernel = ctypes.WinDLL('kernel32', use_last_error=True)
    kernel.OpenProcess.restype = ctypes.c_void_p
    kernel.OpenProcess.argtypes = [ctypes.c_ulong, ctypes.c_int, ctypes.c_ulong]
    kernel.CloseHandle.argtypes = [ctypes.c_void_p]
    kernel.GetExitCodeProcess.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_ulong)]
    handle = kernel.OpenProcess(0x1000, False, pid)
    if handle:
        code = ctypes.c_ulong()
        try:
            if not kernel.GetExitCodeProcess(handle, ctypes.byref(code)):
                raise OSError(ctypes.get_last_error(), 'Cannot read process exit status')
            return code.value == 259, None if code.value == 259 else code.value
        finally:
            kernel.CloseHandle(handle)
    if ctypes.get_last_error() == 87:  # ERROR_INVALID_PARAMETER: no such PID
        return False, None
    raise OSError(ctypes.get_last_error(), 'Cannot establish process absence')


def record(attempt, expected_pid, evidence):
    path = Path(attempt) / 'rerun_status.json'
    rec = json.loads(path.read_text(encoding='utf-8'))
    if rec.get('pid') != expected_pid or rec.get('status') != 'running':
        raise ValueError('Unexpected PID/status; refusing to change record')
    active, exit_code = process_state(expected_pid)
    if active:
        raise RuntimeError('PID exists; do not mark an active or reused PID interrupted')
    rec.update(status='process_interrupted', exit_code=exit_code,
               interruption_observed_at=time.strftime('%Y-%m-%d %H:%M:%S'),
               interruption_evidence=evidence,
               exit_evidence='GetExitCodeProcess on recorded PID; null means unavailable; an exit code does not establish full calculation completion',
               scientific_acceptance='not accepted: computation interrupted before required completion artifacts')
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(rec, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    temp.replace(path)
    print(f'Recorded interruption, not completion: {attempt}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('attempt', type=Path)
    ap.add_argument('--expected-pid', type=int, required=True)
    ap.add_argument('--evidence', required=True)
    args = ap.parse_args()
    record(args.attempt, args.expected_pid, args.evidence)
