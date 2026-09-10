# Phase 24 validation record — 2026-09-10

The complete fixed-seed pipeline completed with exit code 0. Final run: seed 24,
98.19 seconds in the available Windows primary Python environment. No real
instrument or chemical experiment was run. Original PPTX was not published.

## Checks performed

- 10 Phase 24 self-tests passed: 48 valid chiral structures, plate geometry and
  volumes, HPLC recovery, malformed traces, Pareto/joint targets, database
  uniqueness and provenance, qEI exclusion of observed points, protocol syntax,
  JSON serialization and rejection of misaligned source plates.
- All 41 existing Phase 1–19 audit regression tests passed. Two pre-existing
  Phase 17 escape-sequence warnings remain unrelated to this addition.
- Independent artifact validation passed: 342 manifest hashes, 312 raw-trace
  hashes, 312 unique conditions, 9,216-row condition catalog, SQLite integrity
  and foreign keys, matching 24/96/96/96 plate maps, and acquisition training
  IDs restricted to previous batches.
- Official Opentrons **8.8.2**, API **2.15**: all three batch protocols completed,
  each producing 676 simulator commands, 36 P20 and 12 P300 tip-column pickups
  (288 and 96 individual tips), 25 C commissioning temperature and a final
  temperature-module deactivation. The exported root protocol matches round 3
  byte-for-byte. Simulation used default calibration, not physical deck offsets.
- All four generated 300-dpi figures were visually inspected for readable axes,
  labels, plate positions and explicit simulation disclosure.

Machine-readable records: [artifact audit](results_phase24/artifact_validation.json),
[official simulator](results_phase24/opentrons_validation.json),
[metrics](results_phase24/metrics.json), [versions and hashes](results_phase24/manifest.json).

## Numerical results and interpretation

312 synthetic measurements passed QC. 126 met both yield >=90% and |ee| >=95%.
The initial 24 already contained one target hit. Best joint condition: ID 5632,
round 1, C1, yield 96.1502%, |ee| 98.4340%, conversion 98.0056%.
The chromatography RMSEs against hidden synthetic truth were 0.00137, 0.00424
and 0.00214 percentage points for conversion, product yield and signed ee.
Minimum estimated R/S resolution was 2.0248. These low errors reflect the
shared known synthetic peak model and do not establish real analytical precision.

The toy landscape has strongly aligned yield/ee objectives and an easy feasible
region. The result demonstrates plumbing and acquisition execution, not a hard
chemistry benchmark, newly discovered asymmetric reaction, statistically proven
algorithm advantage, calibrated uncertainty, or guaranteed three-round success.

## Environment and execution notes

Primary Python 3.12.14, NumPy 2.3.5, SciPy 1.18.1, Matplotlib 3.11.1,
scikit-learn 1.9.0, RDKit 2026.03.6; OMP/BLAS thread count limited to two.
The pre-existing `phase2ff` environment crashed in native numerical code, so it
was not used for accepted calculations. RDKit EHT printed short-distance
warnings for approximately 0.97–0.98 A O–H bonds in MMFF amino-alcohol structures;
no descriptor or conformer failure was silently replaced with arbitrary numbers.

Opentrons was installed separately because its NumPy constraints differ.
Version 9.1.2 rejected OT-2 protocols; 8.8.2 was then installed and all official
simulations completed. Default-calibration warnings are retained as a limitation.

The supplied PPTX changed from 16 to 15 slides during work. The current file was
re-read and its SHA-256 verified; the deleted visualization slide did not change
the platform requirements or budget text. Printed page labels 14/15 correspond
to physical file slides 13/14 in the current version.

Existing measurement directories are protected against overwrite. For a new
run, choose a new `--out` path, for example `phase24_runs/seed24`. The generated
database, raw CSVs and provenance hashes should be kept together. No live HTTP,
MQTT, vendor hardware, pressure/inert handling, sealing or LC integration has
been tested. The chemistry-enabled protocol branch is not commissioned.
