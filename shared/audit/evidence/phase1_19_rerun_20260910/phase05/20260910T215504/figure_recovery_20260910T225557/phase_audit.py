"""Shared numerical safeguards for the Phase 1-5 audit (no import side effects)."""
import math

AUDIT_VERSION = 'phase1-5-validation-v1'


def converged_ensemble(ids, results):
    """Keep finite, converged structures only; fail closed if none survive."""
    if len(ids) != len(results):
        raise ValueError('Conformer/result length mismatch')
    accepted = [(cid, float(energy)) for cid, (flag, energy) in zip(ids, results)
                if flag == 0 and math.isfinite(energy)]
    if not accepted:
        raise ValueError('No finite converged conformer; 3D analysis is unavailable')
    return [x[0] for x in accepted], [x[1] for x in accepted]


def periodic_angle_error(actual, target):
    return abs((actual-target+180.) % 360.-180.)


def ts_validation(freqs_ts, freqs_r, force_ts, force_r, neb_ok,
                  mode_overlap, force_tolerance=.05):
    """Numerical saddle screening, NOT IRC/chemical-connectivity validation."""
    finite = bool(freqs_ts and freqs_r) and all(math.isfinite(x) for x in [
        *freqs_ts, *freqs_r, force_ts, force_r])
    checks = dict(finite=finite, neb_converged=bool(neb_ok),
                  one_significant_imaginary=sum(x < -20. for x in freqs_ts) == 1,
                  reactant_minimum=bool(freqs_r) and min(freqs_r) >= -20.,
                  stationary=math.isfinite(force_ts) and math.isfinite(force_r)
                  and max(force_ts, force_r) <= force_tolerance,
                  reactive_mode=mode_overlap is not None and math.isfinite(mode_overlap)
                  and mode_overlap >= .3)
    return dict(version=AUDIT_VERSION, passed=all(checks.values()), checks=checks,
                irc_verified=False, force_tolerance_ev_a=force_tolerance,
                imaginary_threshold_cm=-20., mode_overlap_threshold=.3,
                scope='Numerical saddle candidate; IRC connectivity not established')


def exploratory_corrections(drop, stereo):
    """Explicit scenario caps preserve sign and never imply a measured correction."""
    if not all(math.isfinite(x) for x in (drop, stereo)):
        raise ValueError('Non-finite exploratory correction')
    return dict(d1_drop=max(-8., min(8., drop)),
                d2a_drop=max(-8., min(8., drop)),
                ddG_stereo=max(-1.5, min(1.5, stereo)),
                evidence='assigned_sensitivity_scenario_not_prediction')


def cross_nonbonded_system(nb, ligand_indices):
    """Independent OpenMM kernel for unscreened interfragment LJ+Coulomb.

    This checks the explicit pair sum, not cutoff/GB/PME or binding free energy.
    Cross-fragment exceptions are rejected: this check is for noncovalent binding.
    """
    from openmm import System, CustomNonbondedForce
    n = nb.getNumParticles(); ligand = set(map(int, ligand_indices))
    receptor = set(range(n))-ligand
    if not ligand or not receptor or not ligand <= set(range(n)):
        raise ValueError('Invalid receptor/ligand partition')
    if nb.getNumParticleParameterOffsets() or nb.getNumExceptionParameterOffsets():
        raise ValueError('Parameter offsets require a separate alchemical validation')
    for k in range(nb.getNumExceptions()):
        i,j,*_ = nb.getExceptionParameters(k)
        if (i in ligand) != (j in ligand):
            raise ValueError('Cross-fragment exception: cannot use noncovalent pair check')
    force = CustomNonbondedForce(
        '4*sqrt(epsilon1*epsilon2)*((sigma/r)^12-(sigma/r)^6)'
        '+138.935457644382*q1*q2/r; sigma=(sigma1+sigma2)/2')
    for name in ('q','sigma','epsilon'): force.addPerParticleParameter(name)
    system = System()
    for k in range(n):
        system.addParticle(1.)
        force.addParticle(nb.getParticleParameters(k))
    force.addInteractionGroup(ligand,receptor)
    force.setNonbondedMethod(CustomNonbondedForce.NoCutoff)
    system.addForce(force)
    return system
