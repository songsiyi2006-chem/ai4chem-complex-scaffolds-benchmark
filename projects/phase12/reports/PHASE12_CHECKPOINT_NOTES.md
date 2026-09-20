# Phase12 explicit training checkpoints

New full-size run (choose a NEW attempt directory):

```powershell
& 'C:/Users/HUIWEI/.codex/tools/chem-ai4s/venv/Scripts/python.exe' run_phase12_hamiltonian_law_discovery.py --checkpoint-dir results_rerun_audit/phase12-checkpoint-attempt-01
```

Resume that attempt with the same source, environment and options:

```powershell
& 'C:/Users/HUIWEI/.codex/tools/chem-ai4s/venv/Scripts/python.exe' run_phase12_hamiltonian_law_discovery.py --checkpoint-dir results_rerun_audit/phase12-checkpoint-attempt-01 --resume
```

- Coordinator integration: fresh full runs now checkpoint by default to `results_phase12/training_checkpoints`, which is unique inside each campaign snapshot. `--no-checkpoint` explicitly disables this; resume still requires an explicit existing directory. Existing checkpoint directories are never overwritten. Training, sample sizes, steps, loss, selection and held-out audit are unchanged. No full training was run for this change.
- Only H/V scalar training is resumed. Upstream simulation/discovery is recomputed, and its exact training arrays (shape, dtype, SHA256) and symbolic RHS including coefficients must match before a network resumes. Validation data never enters checkpoint training metadata or fitting.
- Manifest rejects changed source, quick/full mode and runtime before upstream work. Per-network metadata also rejects changed total steps, learning rate, seed, training arrays, pairs, transverse samples and RHS. A changed script, even comments, requires a new attempt. The total cosine schedule cannot be extended during resume.
- A checkpoint contains model weights, scaling buffer, Adam state, cosine scheduler, completed step count, Python/NumPy/Torch CPU RNG states. Both networks are CPU models. An already complete H or V performs no further training. V starts fresh if interruption occurred before its initial checkpoint.
- Initial step zero, every 100 completed steps (configurable with `--checkpoint-every N`) and the final step are saved. An interrupted partial optimizer step is never published; resume replays from the last successful checkpoint. At most N steps of training are lost during a failed checkpoint write. There is no termination handler trying to save potentially inconsistent state.
- Atomic same-directory replace follows flush/fsync. There are two published network files, one manifest, and fixed `.tmp` staging names reused after hard interruption; storage does not grow with training steps. Atomic replacement protects the last successful file against normal process/save interruption, not all hardware/power/filesystem failures. Use only one process per checkpoint directory.
- A fresh run refuses an existing checkpoint directory; resume requires its manifest and Hamiltonian weights. No old attempt is moved, deleted or converted. Old interrupted runs lacking weights **cannot resume**. Existing final H/V export files alone lack optimizer/scheduler/RNG state and are not resumable training checkpoints.
- Normal scientific outputs still use the existing `results_phase12` and `figures_phase12` locations. The checkpoint directory isolates training state only; retain the coordinator's existing output snapshot/archive procedure before starting a new run to preserve prior scientific artifacts.
- Microtests: `python -m unittest test_phase12_checkpoint -v`. Eight-row, five-step CPU tests cover both losses, interruption during an optimizer step, exact final weights/optimizer/scheduler/RNG equality, completed-state no-op, disabled-checkpoint equivalence, mismatch/missing-state rejection, old-attempt preservation and failed atomic writes. These verify software behavior, not scientific validity.
