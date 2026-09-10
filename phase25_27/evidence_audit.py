"""Create explicit acceptance ledger and analytical precision diagnostics.

The software run has no authority to promote missing physical evidence to PASS.
"""
from .build_inputs import OUT,dump
from .analysis import R_KCAL,KB_EV,wilson_interval,standard_state_correction
import math


def main():
    records=[]
    gates={25:[
      ('25.1','PARTIAL','12 prospective catalyst-axis entries, 8 E/Z imines, 2 solvents; absolute catalyst CIP and full reactive atom mapping pending'),
      ('25.2','NOT_RUN','19-species/26-edge proposed topology; no computed complete R/S pathway network or validated exchange timescales'),
      ('25.3','PARTIAL','144 MMFF starting conformers; no complex/TS DFT ensemble or correlated-wavefunction reference'),
      ('25.4','NOT_RUN','No accepted TS, frequencies, IRC, standard-state solution free energies, or explicit-solvent PMFs'),
      ('25.5','NOT_RUN','Kinetic and uncertainty kernels tested; no physical rates, conversion or signed ee'),
      ('25.6','PARTIAL','96 conditions, family-grouped split, 10 fixed seeds; no external labels or cost-matched model/AL evaluation'),
      ('25.7','NOT_PASSED','No ensemble ddG convergence, free-energy replication or external interval-coverage evidence')],
     26:[
      ('26.1','PARTIAL','48 condition entries and 3 dry slab seeds; 0 complete solvated interfaces'),
      ('26.2','NOT_RUN','Legendre/charge/potential checks implemented; no grand-canonical DFT path points'),
      ('26.3','NOT_RUN','No explicit water/cation sampling or thickness/area convergence'),
      ('26.4','NOT_RUN','Competing reactions specified; no activation PMFs or committor tests'),
      ('26.5','NOT_RUN','No calibrated site-conserving surface model or mass-transfer solution'),
      ('26.6','NOT_RUN','Controls and degree-of-rate-control protocol specified, no evaluated mechanisms'),
      ('26.7','NOT_PASSED','0/3 independent interfaces; 0.05 eV and 0.05 V targets untested')],
     27:[
      ('27.1','PARTIAL','6-member design; literature parent coordinates recovered; derivative/adduct/isomer structures not complete'),
      ('27.2','NOT_RUN','Published active-space details audited; no new converged multireference benchmark'),
      ('27.3','NOT_RUN','No new SOC/NAC/gradients/crossing searches'),
      ('27.4','NOT_RUN','3 excitation windows specified; 0 nonadiabatic trajectories'),
      ('27.5','PARTIAL','Censor-aware competing-risk analysis tested; no measured trajectories'),
      ('27.6','NOT_RUN','No generated absorption/TA spectra, state lifetimes or quantum yields'),
      ('27.7','NOT_PASSED','0/3 independent trajectory batches; branch CI precision not evaluated')]}
    for phase,items in gates.items():
        for id,status,evidence in items:records.append(dict(phase=phase,criterion=id,status=status,evidence=evidence))
    dump(OUT/'acceptance.json',dict(overall='NOT_SCIENTIFICALLY_COMPLETE',criteria=records,
      phase25_selectivity='选择性未判定 / selectivity undetermined',phase26_mechanism='UNDETERMINED',phase27_branch_quantum_yields=None))
    T=298.15
    dump(OUT/'analytical_precision.json',dict(evidence_kind='analytical_design_diagnostics_not_chemical_predictions',T_K=T,
       RT_kcal_mol=R_KCAL*T,kBT_eV=KB_EV*T,standard_1atm_to_1M_kcal_mol=standard_state_correction(T),
       ee_magnitude_for_hypothetical_0_5kcal_single_ddG=math.tanh(.5/(2*R_KCAL*T)),
       rate_ratio_for_hypothetical_0_05eV_barrier_change=math.exp(.05/(KB_EV*T)),
       wilson_n385_near_half=wilson_interval(192,385),
       warnings=['single-ddG ee identity requires validated kinetic assumptions','binomial Wilson precision assumes independent equal-horizon fully observed absorbing outcomes; censored trajectories require competing-risk treatment','none of these diagnostics estimates an actual system response']))


if __name__=='__main__':main()
