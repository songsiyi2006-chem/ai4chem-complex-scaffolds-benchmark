# Basis-consistent scan analysis

The current same-campaign continuation has a 6-31g fallback at R=1.75 A;
the other 7A points use def2-svp. Raw absolute energies from different basis
sets must not be joined into a single potential-energy curve. The merger now
selects the basis at the reference point (or the first documented basis if
unavailable), masks other/missing bases, and records exclusions in
`basis_consistency`. Raw checkpoint files remain unchanged. Missing points
remain missing, not zero or interpolated successes.

Two actual-source tests cover mixed and unknown bases. This corrects analysis,
not the failed SCF point. The already-running old source snapshot must remain
untouched for its hash-checked continuation; once finished, a separately
provenanced merge with the corrected script is required. Do not certify its
original mixed-basis merged curve as the corrected result.
