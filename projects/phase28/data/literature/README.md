# Public-source numerical data

The three `*.csv.gz` files preserve every numeric time/ordinate pair from the published Figure 2 source workbooks, with the original Excel row. They contain 196,496 source rows, **not** independent electrode replicates. Original workbooks, full articles and supplementary PDFs are not redistributed here. Download URLs, access times and SHA256 are recorded in [source_manifest.json](source_manifest.json); source-specific terms remain applicable. Repository code licensing does not assert ownership of these third-party measurements.

`trace_metadata.json` records axis/caption checks. The Nature 2023 constant-potential current axis has an unresolved normalization; its raw ordinate is preserved and no absolute current/charge is inferred. The other two traces are **anode potentials versus NHE**, not full-cell voltages. None enters the synthetic degradation fit.

`extraction_manifest.json` links each extracted file to source SHA256, sheet and columns. Duplicate timestamps are retained; descriptive statistics include both raw-record-weighted and timestamp-balanced endpoint windows. Hourly percentile bands show within-trace spread and are not confidence intervals for an electrode population.

`nature2023_energy_source_audit.csv` retains columns B–E and the stored mean/s.d. separately. Figure 2c reports an ODC-based extrapolated energy metric, so failure to reproduce its stored mean by direct arithmetic does not establish a paper error. The mapping/transformation is unresolved in this bounded audit.

Re-download the two numerical workbooks into a new folder using `python projects/phase28/code/fetch_sources.py --out work/source_recheck`; then run `literature_analysis.py --source work/source_recheck --destination work/reextracted`. The fetcher verifies hashes and never overwrites an existing directory.
