# Phase25–27 implementation

Actual local continuation: [Chinese progress report](../PHASE25_PROGRESS_ZH.md) / [English progress report](../PHASE25_PROGRESS_EN.md). Do not regenerate prospective inputs over an active campaign.

## Continuing the real local pilot

Run Psi4 through its existing environment. On Windows, **use `-X utf8`**: the native finite-difference driver reads UTF-8 output through QCSchema; a GBK default caused the retained first frequency failure. `dftd3`, `cffi` and `pycparser` must be discoverable by this interpreter. In this workspace they are in the isolated `work/qm_dependencies` extension, selected through the process `PYTHONPATH`.

```sh
python -X utf8 -m phase25_27.stationary_pilot --input results_phase25_27/phase25/starting_structures/I01_E.xyz --tag YOUR_NEW_OPT_TAG --stage opt --memory-mb 500 --threads 2
python -X utf8 -m phase25_27.stationary_pilot --input results_phase25_27/phase25/stationary_pilot/YOUR_NEW_OPT_TAG/optimized.xyz --tag YOUR_NEW_FREQ_TAG --stage freq --memory-mb 500 --threads 3
```

Optimization defaults to PBE-D3BJ/def2-SVP, 75x302 grid and GAU_TIGHT geometry criteria. Frequency uses the same method/basis/grid and native finite differences of gradients. Each successful native QCSchema job is hash-bound and checkpointed. Reusing the same frequency tag requires `--resume`, identical input bytes and method/stage; the exact native QCSchema input hash also guards basis/options. A partial Hessian never yields accepted frequencies. To run the real replay/corruption check, use `python -X utf8 -m phase25_27.checkpoint_smoke` in the Psi4 environment.

The finite local queue is `python -m phase25_27.local_queue --config PATH_TO_CONFIG`. Its config records repository root, Psi4 interpreter, dependency path, thread/memory limits and explicit jobs. It adopts matching live module/tag/cwd processes before submitting anything. Its status and logs are stored beside that config (`work/local_campaign` in this workspace), and it does not publish automatically. Keep the computer/processes running. Inspect status before starting another queue; do not terminate other projects' Python processes.

Run the inexpensive pilot in a Python environment with the analysis dependencies, passing the existing xTB executable:

```sh
python -m phase25_27.ternary_starts
python -m phase25_27.xtb_pilot --exe PATH_TO_XTB --tag YOUR_NEW_XTB_TAG --hessian --threads 1 --inputs PATH_TO_XYZ
python -m phase25_27.pilot_validation
python -m phase25_27.thermo_sensitivity --exe PATH_TO_XTB
python -m phase25_27.pilot_summary
python -m phase25_27.ternary_analysis
python -m unittest phase25_27.test_analysis phase25_27.test_pilot -v
```

`ternary_starts` requires a fresh destination. `xtb_pilot` uses GFN2-xTB/ALPB(toluene), defaults to vtight optimization, then the same-surface Hessian; `--solvent ch2cl2` selects the second solvent and `--single-point` runs coordinate-fixed controls. The earliest execution ledgers concatenated repeated frequency printouts. After the finite campaigns finish, `normalize_xtb_ledgers` corrects that metadata and archives each original ledger without changing raw outputs. The independent `pilot_validation` parser reads raw output, verifies 3N, and is authoritative. It rejects imaginary/near-zero internal modes even when the engine thermochemistry cutoff reports zero imaginary frequencies. `thermo_sensitivity` reuses saved Hessians with native `xtb thermo` at cutoffs 50/100 and validates the 50 replay against the original. No fitted inertia model substitutes for that replay.

`snapshot_live` is workspace-specific: it freezes completed, checked E-Hessian atomic records and snapshots the local queue. Live/queued outputs listed in [publication omissions](../results_phase25_27/publication_omissions.json) are excluded from the hash manifest and Git until completion review. After updating audited data and snapshot, regenerate `progress_reports`, `pilot_figures`, `evidence_audit`, and `release_audit`. Check that reports and manifest refer to a single frozen snapshot before pushing.

## Initial release workflow

The package prepares prospective inputs and validates analysis primitives. It does **not** implement a production transition-state search, grand-canonical electronic-structure engine, solvent free-energy sampler, SOC/NAC trajectory engine or an external experimental oracle. Missing chemical results stay null.

From the repository root:

```sh
python -m pip install -r requirements_phase25_27.txt
python -m phase25_27.build_inputs
python -m phase25_27.network_spec
python -m phase25_27.build_surface_seeds
python -m phase25_27.evidence_audit
python -m unittest phase25_27.test_analysis -v
```

The first command installs software; no quantum jobs are submitted. `build_inputs` regenerates deterministic prospective structures and resets their unmeasured labels. Preserve an existing campaign directory before extending it with production data.

Optional original-source recovery (pypdf needed only for extraction):

```sh
python -m pip install pypdf
python -m phase25_27.fetch_sources --download-dir ../source_cache --extract
```

Pinned source hashes are checked before extraction. Coordinate redistribution must retain its source attribution and CC BY-NC 4.0 conditions reported by the Figshare deposits. Full manuscripts/SI PDFs are not redistributed in this release.

Optional **real isolated-imine single-point** pilot in a separately installed Psi4 1.11 environment:

```sh
python -m pip install dftd3==1.6.0 cffi pycparser
python -m phase25_27.local_qm_pilot --threads 4 --run-id YOUR_UNIQUE_RUN_ID --methods pbe-d3bj wb97x-d
```

The command performs gas-phase single points on full isolated E/Z imines at def2-SVP, with 3 GB requested memory and a 75×302 DFT grid. It does not optimize, compute frequencies, solvent free energies or catalytic complexes. Existing result records append, raw log filenames must be unique, and missing backends/interruptions are not converted to numerical estimates. Current successful outputs are included for inspection.

After source extraction and the pilot ledger exist:

```sh
python -m phase25_27.release_audit
```

This audits hashes, structure and split invariants, and writes `software_validation.json` and `file_manifest.json`. The scientific acceptance state remains NOT_PASSED regardless of software test success.

Interfaces for future production data are described in the three phase protocols and bilingual reports. `Network` supports closed homogeneous solution kinetics only; it is **not** a completed surface site/transport model. `competing_risk_incidence` assumes terminal event labels and non-informative censoring; it is not a trajectory event detector. The TS metadata check is a prefilter, not an independent quantum-output parser.
