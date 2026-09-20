# New public-data inputs

`production_fe_summary_inputs.csv`: Figure2e/f, columnsA:C, 25timepoints each; means and SD from four independent experiments per caption. Raw replicate values and temporal pairing are unavailable. Units were checked visually: molCl2, FEpercent, hours.

`polarization_inputs.csv`: Figure2aA:B(iR-corrected1a),D:E(1a), andFigure2hD:E(after100h),G:H(after200h).777pairs retain Excel rows/columns. The source has no numerical before-durability trace within2h;2a is only an unmatched cross-panel proxy. Post-test iR status is unresolved.

`source_manifest.json` records publisher URL, exact workbook hash, axis/caption review and extraction hashes. Raw XLSX is not redistributed. Re-extract with the existing `additional_public_analysis.py --extract-workbook SOURCE.xlsx --extract-to NEW_DIRECTORY` option; an unexpected source hash is rejected.

These are published numerical facts, not new measurements by this repository. Original2023/2026 long traces remain in `../literature/` and are referenced by hash by the extension. No synthetic values enter these analyses.
