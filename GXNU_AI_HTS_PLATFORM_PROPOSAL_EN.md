# GXNU AI-assisted high-throughput catalysis pilot: executive brief

## Page 1: Scope and demonstrated software capability

The current 15-slide GXNU platform deck dated 4 September 2026 calls for inert handling, parallel synthesis and structured research data (slides 3–4). References to 14/15 use printed labels, corresponding to file slides 13/14. This pilot supplies an executable software/data template for commissioning. It is not a commissioned synthesis platform. Slide 12 specifies performance-based, brand-neutral procurement. The OT-2 implementation is one replaceable adapter, not a purchase recommendation.

The virtual library contains 48 stereochemical structures (24 connectivities and their enantiomers), six metal precursors, eight bases and four solvents: 9,216 conditions. RDKit supplies TPSA, Gasteiger charges, geometric steric-envelope proxies, a point-charge dipole proxy and an extended-Huckel frontier gap. These are not DFT or measured ligand properties; the molecules are not certified stock items.

Twenty-four diverse seeds precede three 96-well recommendations. A GP uses correlated posterior samples for greedy joint qEI on min(yield/90, |ee|/95). Chromatographic fits, not oracle labels, enter the learning database. This run produced 312 synthetic records, 312 passing QC, and 126 simultaneous target hits. The best joint condition yielded 96.15% with 98.43% |ee| in the simulation. Success is reported as observed and is not guaranteed in real chemistry.

Deliverables include a standalone script, SQLite database, raw detector traces, plate/source maps, an API 2.15 protocol, four 300-dpi figures and an auditable metrics manifest. Existing supplier controllers and safety interlocks remain responsible for real equipment.

<div style="page-break-after: always;"></div>

## Page 2: Integration and decision gates

Slide 11 describes JSON task specifications, HTTP requests and MQTT state returns. Local schemas and example events demonstrate this contract; no live server, broker or vendor link has been commissioned. Preserve task IDs, well positions, provenance, units and raw-file hashes across adapters. Reject duplicates and quarantine QC failures before model updates.

Eight-channel pipettes actuate all eight tips together. Three pre-arrayed source plates enable different chemistry in every destination well while a 12-channel reservoir supplies common substrate. Preparation of those source plates is a separate, reviewed workflow. Dispensing totals 80 uL/well, 288 P20 tips and 96 P300 tips. Default execution uses water at 25 C. A 60 C chemistry branch requires independent approval of materials, atmosphere, sealing and thermal/pressure risks. No pressurized hydrogen handling is implemented.

Proposed gates: weeks 0–2 for schemas and software tests; weeks 3–4 for gravimetric water/dye commissioning; weeks 5–8 for a supervisor-selected benchmark and chiral HPLC validation. These are planning assumptions. Students can assist with labels, database checks, simulations and supervised calibration records. Qualified staff retain authority over utilities, gas systems, chemical handling and experiment approval.

Slide 14 lists a 3.50-million-CNY total, 3.00-million-CNY discounted price and 198,094 CNY of renovation. The inclusion relationships require written reconciliation. No price gap or incremental OT-2 budget is inferred here. The next approval requires a defined substrate/product, calibrated chiral method, stock-solubility evidence and documented vendor interfaces.
