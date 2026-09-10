# Phase 3 production checkpoint contract

Scope: `run_phase3_complex_dynamics.py`, `test_phase3_checkpoint.py`, and this note.

Coordinator integration additionally fixes nonzero exit status when the final
renderer or result serialization fails. `test_phase3_delivery.py` exercises
the actual main function with mocked stages (three tests, no molecular engine).
The fresh prepared restart uses the integrated source, not the earlier prepared
copy named at the end of this worker note. See PHASE1_19_CLOSEOUT.md for the
finite queue dependency; no DCD resume is enabled by that launcher.
Implementation uses the explicitly permitted checkpoint-writer + fail-closed
resume option. **Production DCD append/resume is not implemented or enabled.**
`--resume_md` exits before creating directories or running any stage. Calling
Stage4 directly with `resume_md=True` also fails before writing/engine work.
An integrity validator is available for inspection, not permission to resume.

## Persistence and ordering

Fresh Stage4 exclusively claims `stage4_checkpoints/` before system construction.
Existing partial directories/start structures are preserved; `--force_rerun`
cannot overwrite Stage4 artifacts (CLI checks before upstream stages too).
An error during preparation requires a new attempt directory.

After unchanged minimization, 2000-step 0.5-fs heating, 2000-step 1-fs heating,
and requested equilibration (default 5000 steps), production advances to each
absolute reporter boundary and finally to the requested final step. The default
remains **100000 production steps at 1 fs, 200 frames every 500 steps, 40 MMGBSA
frames**. Scientific gates, target structures, analysis formulas and sampling
are unchanged. No checkpoint is promised before the first production report.

For each completed chunk:

1. DCD and analyzer reports must have returned successfully. Check absolute
   progress, DCD writer frame count, analyzer frame count and five-column rows.
2. Flush and fsync DCD; record its full streaming SHA256 and byte size.
3. Create an exact `Context.createCheckpoint()` binary, including velocities,
   time and engine RNG state. XML alone is never treated as restart state.
4. Exclusively write and fsync immutable generation files `state.chk`,
   `progress.json`, `identity.json`; hash all three.
5. Write/fsync a unique temporary manifest, then atomically replace
   `manifest.json` LAST. No existing generation is deleted or overwritten.

The manifest declares `resume_supported: false`. Failed data writes or manifest
replacement leave no new committed generation. If DCD advanced past an older
manifest, the validator rejects its byte/hash mismatch; it never truncates,
appends, silently replays steps, or selects an orphan generation. Generation
retention uses extra disk space, not a growing in-memory binary cache. Files
are fsynced; filesystem/directory metadata and power-loss durability on Windows
are not claimed to be a fully transactional multi-file guarantee.

Identity includes source/input hashes (including missing optional inputs), exact
System and production Integrator XML, configured RNG seed, OpenMM full version,
platform/properties, host identity, sampling/target config, analysis index maps,
full-precision reference positions and topology PDB hash. Progress includes
absolute and production steps, physical context time, exact report steps, chunk
boundaries, all cumulative five-column metric rows and PLIF counts. Hash/identity
validation does not itself prove scientific correctness or portable restart.

## Verification and limits

Existing chem-ai4s runtime reused; live research_status reports OpenMM 8.6.0.
Reference-platform toy contexts only, one atom, fixed seeds and a handful of
steps. Tests use temporary directories and AST-load only the relevant functions.
No ML imports, package installs, full MD, campaign launches or snapshot changes.
No git commit/push. Main owner handles integration and runtime. No callable
subagent interface was available to this owner; independent checks ran in parallel.

Each check uses `OPENMM_CPU_THREADS=1`, `OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `PYTHONDONTWRITEBYTECODE=1` and
`C:/Users/HUIWEI/.codex/tools/chem-ai4s/venv/Scripts/python.exe`.
Exact bounded test commands (inside that interpreter):

```python
import subprocess, sys
subprocess.run([sys.executable, '-B', '-m', 'unittest', '-v',
                'test_phase3_checkpoint'], timeout=30, check=True)
# Separate independent process, same thread limits:
subprocess.run([sys.executable, '-B', '-m', 'unittest', '-v',
                'test_phase3_completion'], timeout=30, check=True)
```

Tests cover actual binary RNG continuation on Reference, identity mismatches,
DCD advance versus stale manifest, corrupt checkpoint data, frame mismatch,
disk/flush/manifest faults, nonfinite metrics, resume/force rejection and final
nonreport chunks. The existing strict default 200-frame completion suite is
retained without changes. Reference toy success does not establish CPU-platform
restart or scientific validity of the molecular run. DCDReporter private fields
`_out`, `_dcd`, `_modelCount` are tested against the installed OpenMM version;
unsupported layouts must fail rather than imply continuity.

Final verification: checkpoint suite **9/9 passed in 1.691 s**; unchanged
completion suite **13/13 passed in 1.204 s**, run in independent parallel
processes with the limits above. Scoped `git diff --check` passed. This is
checkpoint/guard verification, not a Phase3 scientific completion result.

Per current assignment: legacy `old3MD` has only 27/200 frames and no binary
checkpoint, so cannot resume. The prepared `phase03/prepared_20260910T233241`
attempt was not started or modified here; it needs integration of this source
by its owner before a fresh run. These campaign statuses are supplied context,
not a new runtime audit by this owner.
