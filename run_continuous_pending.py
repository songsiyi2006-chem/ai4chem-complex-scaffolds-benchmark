"""One continuously running, finite campaign queue; no recurring scheduler.

Preserves all previous attempts. Failure records are retained for human/agent
review; this runner does not claim scientific acceptance or edit code itself.
"""
import argparse
import json
import os
from pathlib import Path
import time
from rerun_phase1_19_campaign import ROOT, run_one, save
from resume_phase7_rerun import available_gib

MIN_FREE = {16: 1.6, 12: 2.5, 13: 3.5, 15: 2.0, 8: 3.0, 3: 3.5}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--phases', nargs='+', type=int, choices=sorted(MIN_FREE), required=True)
    args = ap.parse_args()
    if len(set(args.phases)) != len(args.phases):
        ap.error('Duplicate phase')
    folder = ROOT/'results_rerun_audit'/('continuous_'+time.strftime('%Y%m%dT%H%M%S'))
    folder.mkdir(parents=True, exist_ok=False)
    state = dict(status='running', pid=os.getpid(), queue=args.phases, completed=[],
                 scientific_acceptance='pending separate scientific review',
                 mode='single continuous process; finite queue; no recurring trigger')
    path = folder/'queue_status.json'
    campaign = ROOT.parent/'phase1-19-rerun-20260910'
    save(path, state)
    try:
        for phase in args.phases:
            # Fail closed on any unresolved attempt, not just a stale PID.
            for prior in (campaign/f'phase{phase:02d}').glob('*/rerun_status.json'):
                if json.loads(prior.read_text(encoding='utf-8')).get('status') == 'running':
                    raise RuntimeError(f'Unresolved existing Phase{phase} attempt: {prior}')
            state.update(current_phase=phase, activity='waiting_for_memory')
            save(path, state)
            while available_gib() < MIN_FREE[phase]:
                time.sleep(5)
            state['activity'] = 'calculating'
            save(path, state)
            # Default scientific parameters. No fast mode, reduced samples or
            # relaxed criteria. Each run retains source/runtime fingerprints.
            result = run_one(phase, campaign, profile='molecular')
            state['completed'].append(result)
            save(path, state)
        state['status'] = 'queue_finished'
        state['activity'] = 'awaiting_scientific_review'
    except BaseException as exc:
        state.update(status='queue_failed', error=repr(exc))
        raise
    finally:
        save(path, state)


if __name__ == '__main__':
    main()
