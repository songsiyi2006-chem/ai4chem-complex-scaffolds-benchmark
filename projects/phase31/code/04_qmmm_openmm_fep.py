"""QM/MM contracts and independently testable free-energy estimators.

No automatic protein repair, force-field guessing, softcore or REST2 emulation.
"""
import math
import numpy as np
from scipy.special import logsumexp
from common import EvidenceBlocked


def validate_partition(n_atoms, qm, mm, boundary_bonds, link_atoms_reviewed=False):
    if type(n_atoms) is not int or n_atoms <= 0:
        raise ValueError("invalid atom count")
    if any(type(i) is not int for i in list(qm) + list(mm)):
        raise ValueError("atom indices must be integers")
    if len(set(qm)) != len(qm) or len(set(mm)) != len(mm) or set(qm) & set(mm):
        raise ValueError("overlap or duplicate QM/MM atoms")
    if not qm or not mm or set(qm) | set(mm) != set(range(n_atoms)):
        raise ValueError("partition must cover each atom exactly once")
    for a, b in boundary_bonds:
        if not ((a in qm and b in mm) or (b in qm and a in mm)):
            raise ValueError("boundary bond is not across QM/MM")
    if boundary_bonds and not link_atoms_reviewed:
        raise EvidenceBlocked("link atoms and charge redistribution require review")
    return True


def lambda_schedule(values):
    x = np.asarray(values, dtype=float)
    if (x.ndim != 1 or len(x) < 3 or not np.all(np.isfinite(x)) or
            x[0] != 0 or x[-1] != 1 or np.any(np.diff(x) <= 0)):
        raise ValueError("lambda must increase strictly from 0 to 1, at least 3 windows")
    return x


def plan(values=None):
    x = lambda_schedule(np.linspace(0, 1, 17) if values is None else values)
    return dict(status="SPECIFICATION_ONLY", lambda_values=x.tolist(), independent_replica_seeds=[3101, 3102, 3103],
                legs=["ligand_in_solvent", "ligand_in_binary_VHL", "ligand_in_ternary_BRD4_VHL"],
                temperature_K=298.15, units="kJ/mol", softcore="REQUIRES_VALIDATED_IMPLEMENTATION",
                rest2="distinct replica dimension, NOT lambda windows",
                required=["atom_mapping", "forcefield_parameters", "charge_and_protonation",
                          "equilibration", "restraint_and_standard_state_corrections",
                          "finite_size_charge_corrections", "overlap_and_effective_samples",
                          "replicate_and_time_convergence", "closed_thermodynamic_cycle"])


def ti(values, mean_dhdl, covariance_of_means):
    x = lambda_schedule(values)
    y = np.asarray(mean_dhdl, dtype=float)
    cov = np.asarray(covariance_of_means, dtype=float)
    if y.shape != x.shape or cov.shape != (len(x), len(x)) or not np.all(np.isfinite(y)) or not np.all(np.isfinite(cov)):
        raise ValueError("shape/nonfinite TI inputs")
    if not np.allclose(cov, cov.T) or np.linalg.eigvalsh(cov).min() < -1e-10:
        raise ValueError("covariance must be symmetric positive semidefinite")
    weights = np.zeros(len(x))
    weights[:-1] += np.diff(x) / 2
    weights[1:] += np.diff(x) / 2
    return dict(delta_g_kj_mol=float(weights @ y),
                standard_error_kj_mol=float(np.sqrt(max(0, weights @ cov @ weights))),
                quadrature_error_included=False, convergence_proven=False)


def forward_fep(delta_u_kj_mol, temperature_K=298.15):
    du = np.asarray(delta_u_kj_mol, dtype=float)
    if du.ndim != 1 or len(du) < 2 or not np.all(np.isfinite(du)) or not math.isfinite(temperature_K) or temperature_K <= 0:
        raise ValueError("finite energy differences and positive temperature required")
    rt = 0.00831446261815324 * temperature_K
    logw = -du / rt
    dg = -rt * (logsumexp(logw) - math.log(len(du)))
    weights = np.exp(logw - logsumexp(logw))
    return dict(delta_g_kj_mol=float(dg), weight_ess=float(1 / (weights @ weights)),
                autocorrelation_corrected=False, convergence_proven=False)


def binding_cycle(complex_leg, solvent_leg, restraint, standard_state, finite_size):
    terms = (complex_leg, solvent_leg, restraint, standard_state, finite_size)
    if any(v is None or not math.isfinite(v) for v in terms):
        raise EvidenceBlocked("all cycle terms required; no missing correction becomes zero")
    return complex_leg - solvent_leg + restraint + standard_state + finite_size


def execute_sampling(*args, **kwargs):
    raise EvidenceBlocked("QM/MM forces, validated alchemical topology and PLUMED/REST2 drivers not wired; no trajectories fabricated")
