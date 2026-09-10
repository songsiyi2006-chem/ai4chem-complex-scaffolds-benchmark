# 第五阶段：审计后的探索性化学模型

以下是条件模型输出，不是已验证势垒降低、分离收率或实验 ee。

The RRHO fallback units are corrected but its translation/rotation estimates remain approximate. Designed kinetic caps require --exploratory-kinetics and retain the sign of the stereo difference. The reaction graph and several rates/floors are assigned, not autonomously discovered. Electronic calculations and full dynamics must be rerun before comparing with historical outputs.

```json
{
  "audit_version": "phase1-5-validation-v1",
  "scope": "assigned reaction network and approximate constrained-energy proxies",
  "barrier_reduction_proven": false,
  "reason": "Unconstrained saddle stationarity and IRC connectivity not established",
  "baseline_conditional_kinetics": {
    "yield_target": 4.646846863903552e-16,
    "ee_pct": 0.0,
    "selectivity_Pt_Pside": 0.06839426792485295,
    "conversion_R": 7.216449660063518e-15,
    "poly_decomp": 4.124694330090995e-45,
    "rate_constants": {
      "r1": 1000000000.0,
      "r2": 3.3437646310056584e+83,
      "r3": 212447379253.63406,
      "r4": 477209091.610657,
      "r5": 8496117.44153304,
      "r6": 8496117.44153304,
      "r7": 248445306.87478092,
      "r8": 9935.733320983056,
      "r9": 1000000000.0,
      "r10": 500.0,
      "r11": 4.753287900222205e-26,
      "r12": 2.0309071779013788e-19
    },
    "rate_constant_units": {
      "r1": "M^-1 s^-1",
      "r2": "s^-1",
      "r3": "s^-1",
      "r4": "s^-1",
      "r5": "s^-1",
      "r6": "s^-1",
      "r7": "s^-1",
      "r8": "s^-1",
      "r9": "M^-1 s^-1",
      "r10": "M^-1 s^-1",
      "r11": "s^-1",
      "r12": "s^-1"
    },
    "internal_time_origin_shift_s": 1e-09,
    "rate_kinds": [
      "assigned",
      "derived_with_assigned_floors",
      "derived_with_assigned_floors",
      "derived_with_assigned_floors",
      "derived_with_assigned_floors",
      "derived_with_assigned_floors",
      "derived_with_assigned_floors",
      "assigned",
      "assigned",
      "assigned",
      "assigned",
      "derived_with_assigned_floors"
    ],
    "nfev": 362,
    "njev": 17,
    "nlu": 111
  },
  "exploratory_scenario": {},
  "constrained_energy_proxy": {
    "winner_smiles": "O=P1(O)Oc2c(c3cc(C(C)C)ccc3)cc3ccccc3c2-c2c(C(C)(C)C)cc3ccccc3c2O1",
    "winner_motifs": [
      "iPr-Ph",
      "tBu"
    ],
    "ddG_bind_winner_TS_minus_R_kcal": 42865.37053668965,
    "barrier_drop_vs_baseline_kcal": -25.176990856271004,
    "baseline": {
      "smiles": "O=P1(O)Oc2ccc3ccccc3c2-c2ccc3ccccc3c2O1",
      "ddG_bind_kcal": 42840.193545833376
    },
    "target_kcal": 4.0,
    "claim_proven": false,
    "threshold_exceeded": false,
    "evidence_status": "unvalidated_constrained_complex_energy_proxy",
    "reason": "No unconstrained saddle/gradient/IRC verification for both catalysts"
  },
  "assigned_parameters": [
    {
      "name": "dGTS2a_trap",
      "value": "8.0",
      "rationale": "post-RDS cation capture; computed raw = 59.73 kcal (relaxed scan collapses into product well)"
    },
    {
      "name": "dGTS2b_trap",
      "value": "6.0",
      "rationale": "C3 deprotonation capture; computed raw = 9.05 kcal"
    },
    {
      "name": "k_on(R+Cat→RC)",
      "value": "1000000000.0",
      "rationale": "Smoluchowski diffusion-limited cap"
    },
    {
      "name": "k_dimer(2I1→P_poly)",
      "value": "500.0",
      "rationale": "cationic oligomerization surrogate"
    },
    {
      "name": "ΔG‡(P_elim→Q+H2)",
      "value": "52.0",
      "rationale": "dehydrogenative aromatization reference barrier"
    }
  ]
}
```
