# Phase25–27 implementation

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
