"""Synthetic effective barriers and an explicit external-DFT handoff contract.

No DFT, AIMD, ligand search or constant-potential calculation is performed here.
The arbitrary parameters exercise coupling and are not literature fits.
"""
import math
from common import KB_EV, KB_H, LABEL, R, F

BASE = {
    "Cu": dict(ads=.84, CO=.68, C2H4=.73, HCOOH=.86, H2=.88),
    "Ag": dict(ads=.82, CO=.61, C2H4=.94, HCOOH=.88, H2=.91),
}
ALPHA = dict(ads=.35, CO=.03, C2H4=.32, HCOOH=.4, H2=.38)
LIGAND = dict(ads=.3, CO=-.6, C2H4=.8, HCOOH=-.3, H2=.2)


def rhe_to_she(potential_RHE_V, bulk_pH, temperature_K):
    return potential_RHE_V - math.log(10) * R * temperature_K / F * bulk_pH


def rates(config, diffuse_potential_V=0., local_pH=None):
    config.validate()
    local_pH = config.bulk_pH if local_pH is None else local_pH
    # Nominal CO reference used only as a declared driving-coordinate convention.
    drive = -.11 - config.potential_RHE_V
    barriers = {key: base - ALPHA[key] * drive + LIGAND[key] * config.ligand_shift_eV
                + .04 * diffuse_potential_V for key, base in BASE[config.material].items()}
    if min(barriers.values()) <= .05:
        raise ValueError("Barrier outside activated-rate envelope; no silent clipping")
    constants = {key: KB_H * config.temperature_K * math.exp(-barrier / (KB_EV * config.temperature_K))
                 for key, barrier in barriers.items()}
    constants["H2"] *= 10 ** (-.15 * (local_pH - 7.))  # declared empirical proton sensitivity
    return constants, barriers


def packet(config):
    constants, barriers = rates(config)
    return dict(schema_version="phase30.demo.v1", evidence_type=LABEL,
                dft_executed=False, aimd_executed=False, material=config.material,
                potential_RHE_V=config.potential_RHE_V,
                potential_SHE_V=rhe_to_she(config.potential_RHE_V, config.bulk_pH, config.temperature_K),
                barriers_eV=barriers, effective_constants_s=constants,
                adsorption_energies_eV=None,
                missing="Adsorption free energies and constant-potential PCET barriers require external DFT",
                standard_state="CO2 activity = c/33 mol m^-3; unit free-site fraction",
                network=[
                    {"id": "ads", "lumped_step": "CO2 + * + 2(H+ + e-) -> CO* + H2O"},
                    {"id": "CO", "lumped_step": "CO* -> CO + *"},
                    {"id": "C2H4", "lumped_step": "2CO* + 8(H+ + e-) -> C2H4 + 2H2O + 2*"},
                    {"id": "HCOOH", "lumped_step": "CO2 + 2(H+ + e-) -> HCOOH"},
                    {"id": "H2", "lumped_step": "2(H+ + e-) -> H2"}],
                limitations=["irreversible lumped network, not detailed-balance mechanism",
                             "ligand_shift_eV is a hypothetical perturbation, not a ligand identity",
                             "HCOOH reports acid-equivalent moles; local formate speciation not resolved"])
