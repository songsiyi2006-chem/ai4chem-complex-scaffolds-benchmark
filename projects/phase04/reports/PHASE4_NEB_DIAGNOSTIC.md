# Phase4 saved-band diagnostic

## Coordinator execution: fixed-band probe completed

The main coordinator executed `audit_phase4_forces.py` on the saved bounded
band, independently of the worker's unexecuted recipe below. Nine GFN2-xTB
energy/gradient evaluations completed, with no optimizer or coordinate updates.
The original input hashes were unchanged. Evidence is under
`force_postaudit_20260910T235842/` within the original bounded attempt.

- Recomputed maximum projected atom force: **0.15066260951037777 eV/A**.
- Independent projection replay maximum atom error: **4.163336342344337e-17 eV/A**.
- fmax bookkeeping error: **0.0 eV/A**; no applied constraints.
- Original force target remains **0.05 eV/A**, and still fails.

This closes the missing-force diagnostic gap for this saved band: the current
force projection/bookkeeping agrees with independent replay. Nonconvergence is
real under this model and geometry, not an identified projection bookkeeping bug.
It does not establish an optimized TS, IRC, or valid reaction rate. Older raw
records and the worker's original findings below are preserved as historical context.

Scope now includes the actual production script `run_phase4_reaction_mechanism.py`
(there is no `run_phase4_reaction_pathway.py`), `diagnose_phase4_neb.py`,
`test_phase4_neb_diagnostic.py`, and this report. Read chem-ai4s and its quantum reference, PHASE3_5_RERUN_NOTES.md,
the actual Phase4 attempts, production source, and installed ASE NEB source.
No calculator, QM, NEB rerun, environment installation, snapshot modification,
commit or push. Main owns integration/runtime. No callable parallel-agent tool
was available in this worker; independent source/evidence reads were parallel.

## Authorized bounded implementation

Integration provenance check: the installed phase2ff ASE 3.29.0
`Lib/site-packages/ase/mep/neb.py:474-476` assigns adjusted `forces[i-1]` from
`images[i].get_forces()`, but separately assigns `real_forces[i]` from
`images[i].get_forces(apply_constraint=False)`. The parallel branches use the
same distinction. Thus the current raw label is correct for this installation.
The JSON now records exact provenance: raw arrays are copied from that final
evaluation, adjusted arrays are reconstructed by applying constraints once to
a raw copy, and projected arrays are the actual final NEB return value.
No additional calculator call obtains raw forces. The regression includes a
non-idempotent half-force constraint as well as FixAtoms, explicitly checks raw
versus once-adjusted values, and would detect a second adjustment. This verifies
the installed implementation, not an assumption about all ASE versions.

The previously recommended export repair is now implemented through
`_export_neb_candidate`, called by both stage2 and bounded refinement. It builds
a fresh calculator-free Atoms object from numbers and positions. The actual
ASE extxyz roundtrip passes with a calculator stub lacking `results`.

Both production paths now call `_save_neb_force_diagnostics` immediately after
their existing final `neb.get_forces()`, following any best-band restoration.
It writes `neb_force_diagnostics.json` in the requested output directory only
when a future authorized run reaches this point. No existing campaign outputs
were modified. Captured values include full-precision positions/cells/PBC,
atom numbers, constraints, same-evaluation `neb.energies`, actual interior raw
forces, constraint-adjusted forces, actual NEB projected forces, per-image and
overall max atom force, the unchanged target, climb/index/method/springs, units,
engine name/charge/multiplicity where available, ASE version, source hash and
float64 geometry hash. Unknown engine version is explicitly null; no engine
process is invoked to obtain it. Endpoint force placeholders are omitted.

The helper uses saved raw arrays and applies the recorded ASE force constraints
to copies; it never asks a calculator for energy/force. Nonfinite captured
numerical evidence raises instead of serializing invalid JSON. This is evidence
collection only, and does not change optimizer steps, force thresholds or TS/IRC
acceptance. Capture is final-state only, not an iteration trajectory.

`replay_force_diagnostics` now checks same-state projection and fmax from this
JSON using nonperiodic improved-tangent equations and the recorded per-spring
values. The CLI compares the replay fmax/energies with stage2 and full-precision
positions with rounded XYZ. Older attempts without force records retain the
historical evidence-gap diagnosis below. Replay assumes the current identity
preconditioning; unsupported method, units or periodicity are rejected.

Exact current regression command (existing phase2ff environment, all four thread
limits below set to 1; no QM, optimizer, ML import or environment installation):

```powershell
$env:OMP_NUM_THREADS='1'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:NUMEXPR_NUM_THREADS='1'
& 'C:/Users/HUIWEI/miniconda3/envs/phase2ff/python.exe' -B -m unittest -v test_phase4_neb_diagnostic
```

**13 tests passed in 2.062 seconds**, shell duration 2.921 seconds (<30 s).
Tests extract only the two production helpers via AST to avoid script startup.
ASE SinglePointCalculator supplies synthetic arrays, not quantum calculations.
A fixed-atom constrained four-image band verifies actual ASE projected forces
against independent replay; force/energy accessors are patched to raise during
capture, proving zero additional evaluations. Raw arrays remain unchanged,
endpoint rows are omitted, the synthetic force gate remains failed, and a
deliberately perturbed projected component is detected. Temporary test files
are cleaned automatically; historical snapshots remain untouched.

## Findings

### Next allocated probe: one fixed-band evaluation, NOT EXECUTED

Use this exact PowerShell invocation from the assigned repository directory.
It imports the current production module without calling main, builds the nine
saved images with current `_PerAtomCalc`/GFN2-xTB, calls `NEB.get_forces()` once,
and invokes the production diagnostics writer. No optimizer is constructed.
The unique output directory is newly created under this repository; an existing
directory is rejected. Input files are hashed before and checked after the
probe. The 0.05 eV/A target is read from the failed record and must remain 0.05.

A 60-second shared wall-clock budget applies to xTB subprocess invocations,
including their normal wait/cleanup; Python import and final serialization add
overhead. This is an allocation ceiling, not a completion-time promise. A timeout
preserves failure provenance in the new folder and does not schedule a retry.
`subprocess.run` handles termination of its timed-out xTB child. No historical
record is overwritten and no TS/IRC acceptance is produced.

```powershell
Set-Location 'C:/Users/HUIWEI/Documents/Codex/2026-09-08/https-github-com-songsiyi2006-chem-ai4chem-2/work/phase1-5-fixes'
$env:OMP_NUM_THREADS='1'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:NUMEXPR_NUM_THREADS='1'
$env:PATH='C:/Users/HUIWEI/miniconda3/envs/phase2ff/Library/bin;C:/Users/HUIWEI/miniconda3/envs/phase2ff;' + $env:PATH
@'
import hashlib, json, time
from datetime import datetime
from pathlib import Path
from unittest.mock import patch
from ase.io import read
from ase.mep.neb import NEB
import run_phase4_reaction_mechanism as p
from diagnose_phase4_neb import replay_force_diagnostics

source = Path('../phase1-19-rerun-20260910/phase04/20260910T181601/results_phase4').resolve()
band = source / 'neb_final_path.xyz'
stage = source / 'stage2.json'
inputs = [band, stage]
digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
hashes = {str(path): digest(path) for path in inputs}
previous = json.loads(stage.read_text())
assert previous['fmax_target'] == 0.05
assert p.XTB_EXE, 'Existing xTB executable required; no fallback permitted'
images = read(str(band), index=':')
assert len(images) == previous['n_images'] == 9
assert all(len(im) == 19 and not im.constraints and not im.pbc.any() for im in images)
engine = p.XTBWrap(charge=0, mult=1)
for image in images:
    image.calc = p._PerAtomCalc(engine)
neb = NEB(images, climb=True, k=0.1, method='improvedtangent',
          parallel=False, remove_rotation_and_translation=False)
out = Path('phase4_force_probe_' + datetime.now().strftime('%Y%m%dT%H%M%S%f')).resolve()
out.mkdir(exist_ok=False)
meta = dict(scope='one_fixed_band_force_evaluation_no_optimizer',
            acceptance='no_TS_no_IRC', input_sha256=hashes,
            production_source_sha256=digest(Path(p.__file__)),
            xtb_executable=str(p.XTB_EXE), xtb_executable_sha256=digest(Path(p.XTB_EXE)),
            force_target_eV_A=0.05, spring_eV_A2=0.1, status='started')
provenance = out / 'probe_provenance.json'
provenance.write_text(json.dumps(meta, indent=2), encoding='utf-8')
original_run = p.subprocess.run
deadline = time.monotonic() + 60.0
def bounded_run(*args, **kwargs):
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise TimeoutError('Fixed-band probe exhausted 60-second allocation')
    kwargs['timeout'] = min(kwargs.get('timeout', remaining), remaining)
    return original_run(*args, **kwargs)
try:
    with patch.object(p.subprocess, 'run', bounded_run):
        projected = neb.get_forces()
    p._save_neb_force_diagnostics(out, neb, projected, engine, 0.05)
    path = out / 'neb_force_diagnostics.json'
    evidence = json.loads(path.read_text(encoding='utf-8'))
    # Each endpoint energy evaluation also populated this wrapper's raw-force
    # cache. Read only exact-state caches; never request another evaluation.
    all_raw = []
    for image in images:
        state = (image.get_atomic_numbers().tobytes(), image.get_positions().tobytes())
        assert image.calc._state == state and image.calc._values is not None
        all_raw.append(image.calc._values[1].copy().tolist())
    evidence['all_image_raw_forces_eV_A'] = all_raw
    evidence['all_image_raw_force_provenance'] = (
        'Exact-state _PerAtomCalc._values[1] caches from this single band evaluation; '
        'includes endpoints evaluated for energies; no extra engine calls')
    evidence['state'] = 'one_evaluation_of_saved_rounded_XYZ_no_geometry_updates'
    evidence['probe_input_sha256'] = hashes
    path.write_text(json.dumps(evidence, indent=2, allow_nan=False), encoding='utf-8')
    replay = replay_force_diagnostics(evidence)
    (out / 'force_replay.json').write_text(json.dumps(replay, indent=2), encoding='utf-8')
    assert all(digest(path) == hashes[str(path)] for path in inputs), 'Input changed during probe'
    meta.update(status='evaluated_not_scientifically_accepted', engine_calls=engine.n_calls,
                evidence_sha256=digest(out / 'neb_force_diagnostics.json'),
                input_hashes_rechecked=True)
    print(json.dumps(dict(output=str(out), engine_calls=engine.n_calls, replay=replay), indent=2))
except Exception as exc:
    meta.update(status='failed', error=repr(exc), engine_calls=engine.n_calls)
    raise
finally:
    provenance.write_text(json.dumps(meta, indent=2), encoding='utf-8')
'@ | & 'C:/Users/HUIWEI/miniconda3/envs/phase2ff/python.exe' -B -
```

Expected work is nine xTB energy/gradient subprocesses, one per image, with cache
reuse for ASE's repeated accessor calls. The new JSON contains all nine raw
force arrays plus seven interior adjusted/projected arrays; NEB has no projected
endpoint forces. A successful probe supplies missing force evidence and a
projection/fmax comparison. Any change from the historical fmax reflects a
fresh evaluation of the saved eight-decimal XYZ and current source, not a
retroactive correction or acceptance of the historical attempt. No probe was
executed while preparing this recipe. Scope frozen after recipe handoff.

The saved evidence confirms reported NEB nonconvergence, but does **not** identify
whether the residual is a physical gradient or a projection implementation error.
Neither completed attempt saved raw or projected force arrays. The only force
evidence is scalar fmax in stage2.json and sparsely printed iteration logs.
Plain final-path XYZ and candidate XYZ contain positions only; xTB gradients
were created in temporary directories and were not retained. Do not manufacture
force reconstruction from energies sampled along a path: they do not determine
the transverse gradient.

| Check | Full 20260910T141256 | Bounded 20260910T181601 |
| --- | --- | --- |
| Final reported max atom force, eV/A | 0.15217299888877053 | 0.15066262106165354 |
| Original force threshold, eV/A | 0.05 | 0.05 |
| Saved converged flag | false | false |
| Images / atoms | 9 / 19 | 9 / 19 |
| Interior highest-energy image (zero-based) | 4 | 4 |
| Endpoint displacement vs saved IDPP, A | 0 | 0 |
| Segment lengths in full 57D space, A | 2.79388 to 4.29124 | 2.79446 to 4.28707 |
| Smallest pair distance over all images, A | 1.002286 | 1.002226 |
| Largest spring-only fmax lower bound, eV/A | 0.0142346 | 0.0143558 |
| Candidate XYZ bytes | 1093 | 0 |

No zero/undefined improved tangents occur in either saved band. Reconstructed
tangents have unit norm over all atoms, as ASE requires. Segment lengths are not
single-atom displacements or bond lengths. Minimum contacts alone do not validate
chemical topology. Spring bounds assume the inspected source default k=0.1
eV/A^2; neither recorded launch command overrides it. They are below 0.05, so
spacing alone does not prove failure. The diagnostic CLI explicitly labels its
spring input as an assumption rather than silently recovering it from metadata.

### Tangent, projection, constraint and bookkeeping assessment

Installed ASE 3.29.0 `mep/neb.py` ImprovedTangentMethod selects forward/backward
directions on monotonic energies and energy-weighted directions at extrema.
It normalizes once over all Cartesian components. Ordinary-image force is
`F - dot(F,t)*t + (k_right*L_right-k_left*L_left)*t`. Climbing-image force is
`F - 2*dot(F,t)*t`, without a spring term. The diagnostic reproduces these
equations for nonperiodic molecular images and identity preconditioning.
Degenerate tangents are flagged (ASE's zero-vector fallback is not interpreted
as a valid direction). Tied highest interior energies are explicitly ambiguous.

Climbing reflection preserves the full vector norm, but may change the maximum
per-atom norm. Thus a raw CI-image max atom force cannot simply be substituted
for the NEB convergence metric. Tests cover this distinction and global rather
than per-atom tangent normalization. The lower bound used in the table follows
from orthogonality: `||F_perp+s*t|| >= |s|`, and maximum atom norm is at least
that norm divided by sqrt(number of atoms). This is conditional on correct
projection, not independent proof of its implementation.

Inspected source attaches no atom constraints and uses nonperiodic molecular
XYZ images. Fixed NEB endpoints are intentional and remain unchanged. The cache
is keyed by exact atom numbers/positions and returns copied forces, avoiding
in-place projection corrupting a cached force. Restored best coordinates are
followed by a fresh `neb_fmax` call; cache state changes trigger recalculation.
No demonstrated tangent, constraint, force-copy or final-fmax bookkeeping bug
was found. Missing force arrays prevent excluding such errors experimentally.

Both saved flags correctly reject their reported force values; image counts,
highest interior image, and forward energy-difference arithmetic agree.
Full-run logging shows tail force 4.6001 and restoration of best step41. Bounded
logging shows tail force 0.269966 and restoration before the final 0.150663
record. Logging only every ten iterations explains why the bounded saved best
value need not appear among printed values; no best-step index was persisted.
Neither discrepancy establishes faulty convergence accounting.

## Definite independent source fix for main

The bounded run's actual traceback is not solely the expected figure gate.
At `run_phase4_reaction_mechanism.py:914`, this export fails:

```python
write(str(out / 'ts_candidate.xyz'), images[i_ts])
```

ASE extxyz tries `_PerAtomCalc.results`, which the custom wrapper does not
implement. The file is opened first and remains zero bytes. The failure occurs
after stage2.json, energies and final path were saved, so it does not explain the
NEB force or justify another optimization.

The minimal production change identified during diagnosis is now applied via
the shared export helper: pass a fresh calculator-free Atoms object:

```python
from ase import Atoms
write(str(out / 'ts_candidate.xyz'),
      Atoms(numbers=images[i_ts].get_atomic_numbers(),
            positions=images[i_ts].get_positions()))
```

This writes only a geometric candidate; retain `return passed`, the 0.05 force
threshold, failed status and downstream stage3 stop. It is an export repair,
not a TS or IRC validation. Main can verify it using an in-memory XYZ roundtrip
with a calculator stub lacking `results`, requiring no QM or NEB execution.

The following was the original instrumentation recommendation. Final-state
capture is now implemented as described above; per-step/best-update capture is
not implemented. For a future separately authorized calculation, persist diagnostic evidence at
the already-performed force evaluation, especially when updating the best band
and after restoration. Record full-precision positions, energies, atom order,
charge0/singlet convention, GFN2-xTB/version, method, climb flag/index, springs,
constraints, stage/step and units. Copy `neb.real_forces[1:-1]` immediately after
`f = neb.get_forces()` and copy `f` reshaped to interior-image/atom/xyz dimensions.
ASE's endpoint `real_forces` entries are zero placeholders, not evaluated
physical endpoint forces. With actual constraints, separately retain adjusted
forces or their transform before projection. Save records with geometry/source
hashes; compare reconstructed projection and fmax to the same-state saved arrays.
This would separate real residuals, constraint effects, projection differences,
and stale convergence records without adding calculator calls. Do not change
tangents, springs, thresholds or optimization limits on the present evidence.

## Original read-only reproduction and validation (before implementation)

Run from the assigned `work/phase1-5-fixes` directory. PRIMARY below is the existing
`C:/Users/HUIWEI/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe`.
Both modules use the standard library only. `-B` avoids pycache writes.

```powershell
$env:OMP_NUM_THREADS='1'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:NUMEXPR_NUM_THREADS='1'
& 'C:/Users/HUIWEI/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -B -m unittest -v test_phase4_neb_diagnostic
& 'C:/Users/HUIWEI/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -B diagnose_phase4_neb.py '../phase1-19-rerun-20260910/phase04/20260910T141256' '../phase1-19-rerun-20260910/phase04/20260910T181601'
```

10 mathematical tests passed in 0.001 s (shell call 0.542 s). Actual two-attempt
read-only audit exited0 in a 0.368 s shell call. Each was below the 30 s budget;
no new environment, native engine, ML import, or heavy runtime was used. These
tests verify diagnostic math, not production chemistry. The CLI prints JSON to
stdout, includes evidence SHA256s, and writes nothing. Incomplete earlier attempts
are reported as lacking saved final-band metrics. It does not parse force files
or regenerate missing force evidence; `project_force` is an in-memory math API.

Key read-only evidence hashes:

| Artifact | SHA256 |
| --- | --- |
| Full final path | 48dff5d033981ec75278aa5bade2d151376c39060b59a507e8f56ed07f473943 |
| Bounded final path | 17530240607656bf62c8be2d1c6e04a5e50b4491c2912d75fd9d4382cd9134df |
| Bounded stage2.json | 0ca708bd45bbd6747362885cbcc5c51afd0fd0b30d80b8d057cba9aa4cfd114f |
| Bounded run.log | 05c8558359b2f1e3684088be96011806068ddb4f5a9c8a156c434075d21c193a |
| Bounded source snapshot | 5df76dc703f16232b9a12ab32bccf7e4bcf6bfac07ff536fe1bf1a32d48ed463 |

Scientific disposition remains **no accepted TS, no IRC, no validated barrier or
rate**. A scalar force gate and software test success cannot confer acceptance.
