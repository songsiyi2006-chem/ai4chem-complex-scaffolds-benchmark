"""Conditional equilibrium thermodynamics; no inferred cellular degradation."""
import math
import numpy as np


def alpha_from_kd(binary_Kd_M, conditional_Kd_M, temperature_K=298.15, reversible=True):
    if not reversible:
        raise ValueError("irreversible covalent capture needs kinetics, not equilibrium Kd")
    if any(not math.isfinite(x) or x <= 0 for x in (binary_Kd_M, conditional_Kd_M, temperature_K)):
        raise ValueError("Kd in mol/L and T in kelvin must be finite and positive")
    log_alpha = math.log(binary_Kd_M) - math.log(conditional_Kd_M)
    return dict(alpha=math.exp(log_alpha), delta_delta_g_kj_mol=-0.00831446261815324 * temperature_K * log_alpha,
                definition="same binding leg: binary vs conditional in presence of other partner",
                proves_ubiquitination=False)


def alpha_from_free_energies(binary_dg, conditional_dg, covariance, temperature_K=298.15):
    cov = np.asarray(covariance, dtype=float)
    if (cov.shape != (2, 2) or not np.all(np.isfinite(cov)) or not np.allclose(cov, cov.T)
            or np.linalg.eigvalsh(cov).min() < -1e-10):
        raise ValueError("PSD 2x2 covariance required, including shared-reference covariance")
    if any(not math.isfinite(x) for x in (binary_dg, conditional_dg, temperature_K)) or temperature_K <= 0:
        raise ValueError("invalid thermodynamics")
    rt = 0.00831446261815324 * temperature_K
    contrast = np.array([-1., 1.])
    ddg = conditional_dg - binary_dg
    se = float(np.sqrt(max(0, contrast @ cov @ contrast)))
    loga = -ddg / rt
    return dict(log_alpha=loga, alpha=math.exp(loga), delta_delta_g_kj_mol=ddg,
                standard_error_ddg=se,
                alpha_interval_normal_approx_95=[math.exp(loga - 1.96 * se / rt), math.exp(loga + 1.96 * se / rt)],
                requires_equilibrated_same_standard_state=True)
