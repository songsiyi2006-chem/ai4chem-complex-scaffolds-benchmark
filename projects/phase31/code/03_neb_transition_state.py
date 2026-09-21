"""ASE CI-NEB interface and separate strict TS/IRC evidence gate."""
import math
from common import EvidenceBlocked


def run_neb(reactant, product, calculator_factory, out, images=7, fmax=0.03, steps=150):
    from ase.mep import NEB
    from ase.optimize import FIRE
    from ase.io import write
    from pathlib import Path
    if reactant.get_chemical_symbols() != product.get_chemical_symbols():
        raise ValueError("atom order/identity must match")
    if not 5 <= images <= 31 or fmax <= 0 or steps < 1:
        raise ValueError("invalid bounded NEB settings")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    chain = [reactant.copy()] + [reactant.copy() for _ in range(images - 2)] + [product.copy()]
    for i, atoms in enumerate(chain):
        atoms.calc = calculator_factory(i)  # Distinct calculators, no shared mutable checkpoint.
    if len({id(a.calc) for a in chain}) != images:
        raise ValueError("each image needs a distinct calculator")
    neb = NEB(chain, climb=True)
    neb.interpolate(method="idpp")
    opt = FIRE(neb, logfile=str(out / "neb.log"), trajectory=str(out / "neb.traj"))
    ok = bool(opt.run(fmax=fmax, steps=steps))
    energies = [float(a.get_potential_energy()) for a in chain]
    if not all(math.isfinite(e) for e in energies):
        raise ValueError("nonfinite image energy")
    write(out / "images.xyz", chain)
    return dict(neb_converged=ok, energies_eV=energies, ts_validated=False,
                remaining=["endpoint_minimum_audit", "saddle_hessian", "mode_assignment", "bidirectional_IRC"])


def audit_ts(evidence):
    reasons = []
    if not evidence.get("neb_converged"):
        reasons.append("NEB unconverged")
    freqs = evidence.get("frequencies_cm1", [])
    if not freqs or not all(math.isfinite(x) for x in freqs) or sum(x < -20 for x in freqs) != 1:
        reasons.append("requires one significant imaginary Hessian mode")
    for key in ("endpoint_minima_verified", "imaginary_mode_reaction_verified",
                "irc_forward_connected", "irc_reverse_connected", "same_hamiltonian"):
        if evidence.get(key) is not True:
            reasons.append(key)
    paths = evidence.get("verified_artifact_hashes", {})
    for key in ("hessian", "irc_forward", "irc_reverse", "reactant", "product"):
        value = paths.get(key, "")
        if len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
            reasons.append("missing verified hash: " + key)
    # This contract does not independently authenticate a caller's scientific assertions.
    return dict(contract_passed=not reasons, reasons=reasons, requires_independent_artifact_review=True)


def irc(*args, **kwargs):
    raise EvidenceBlocked("No validated mass-weighted IRC driver is wired; NEB/downhill paths are not renamed IRC")
