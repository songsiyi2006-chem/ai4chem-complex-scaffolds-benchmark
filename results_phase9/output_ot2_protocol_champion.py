# ==========================================================================
# Phase-9 SDL Champion batch — optimized Pareto conditions
# Compiled by run_phase9_self_driving_lab_compiler.py  (Phase 9 / Module 9A)
# Generated: 2026-09-05T08:18:51
# Reaction  : 2-methyl-azirino[1,2-a]indole -> (R)-3-methyl-2,3-dihydroquinoline-imine + (S)-enantiomer (C10H11N asymmetric ring expansion)
# Catalyst  : 3,3'-bis(4-CF3-phenyl)/3-[(iPr)phenyl]-BINOL phosphoric acid (Phase-5C designed winner)
# Safety    : guardrails G1 (dT_ad < 30 K), G2 (T < T_boil - 15 degC),
#             G3 (microfluidic dP < 15 bar) verified at compile time.
# Liquid classes: volatile DCM / toluene / viscous DMSO / MeOH / EtOH
# ==========================================================================

from opentrons import protocol_api

metadata = {
    'protocolName': 'Phase-9 SDL Champion batch — optimized Pareto conditions',
    'author': 'AI4Chem Phase-9 Self-Driving Lab Compiler',
    'description': 'Top-5 feasible Pareto conditions from the closed-loop campaign',
    'apiLevel': '2.15',
}

CONDITIONS = [
    {
        "label": "CHAMP-1",
        "well": "A1",
        "T_c": 30.849,
        "cat_molpct": 8.06,
        "t_h": 20.804,
        "phi_tol": 0.291
    },
    {
        "label": "CHAMP-2",
        "well": "B1",
        "T_c": 45.202,
        "cat_molpct": 2.14,
        "t_h": 23.824,
        "phi_tol": 0.751
    },
    {
        "label": "CHAMP-3",
        "well": "C1",
        "T_c": 60.806,
        "cat_molpct": 1.038,
        "t_h": 15.932,
        "phi_tol": 0.978
    },
    {
        "label": "CHAMP-4",
        "well": "D1",
        "T_c": 40.303,
        "cat_molpct": 9.516,
        "t_h": 22.027,
        "phi_tol": 0.987
    }
]

REAGENT_WELLS = {
    "substrate": "A1",
    "catalyst": "A2",
    "is_stock": "A3"
}
SOLVENT_WELLS = {
    "dcm": "A1",
    "toluene": "A2",
    "hplc_diluent": "A3",
    "meoh": "A4",
    "etoh": "A5"
}

AIR_GAP_UL = 10  # standard air gap for all non-aqueous classes


def _dose(pipette, src, dst, volume_ul, flow_asp, flow_disp, prewet=0):
    """Liquid-class-aware dose: split by pipette capacity, air gap,
    blow-out and touch-tip. Valid for all non-aqueous organic classes."""
    pipette.flow_rate.aspirate = flow_asp
    pipette.flow_rate.dispense = flow_disp
    for _ in range(prewet):  # condition the tip with the actual liquid
        pipette.aspirate(min(volume_ul, pipette.max_volume), src.bottom(4))
        pipette.dispense(min(volume_ul, pipette.max_volume), src.bottom(2))
    remaining = float(volume_ul)
    while remaining > 1e-9:
        v = min(remaining, pipette.max_volume - AIR_GAP_UL)
        pipette.aspirate(v, src.bottom(4))
        pipette.air_gap(AIR_GAP_UL)
        pipette.dispense(v + AIR_GAP_UL, dst.bottom(3))
        pipette.blow_out(dst.top(-2))
        pipette.touch_tip(dst)
        remaining -= v


def run(ctx: protocol_api.ProtocolContext) -> None:
    # ---- deck layout -------------------------------------------------
    tips300 = ctx.load_labware('opentrons_96_tiprack_300ul', '1',
                                   'P300 filter tip rack')
    tips20 = ctx.load_labware('opentrons_96_tiprack_20ul', '2',
                                  'P20 filter tip rack')
    solv = ctx.load_labware('usascientific_12_reservoir_22ml', '3',
                                'Solvent reservoir')
    reag = ctx.load_labware('usascientific_12_reservoir_22ml', '4',
                                'Reagent reservoir')
    temp_mod = ctx.load_module('temperature module gen2', '7')
    plate = temp_mod.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt',
                                      '96-well reaction plate')
    stocks = ctx.load_labware('opentrons_24_aluminumblock_nest_1.5ml_snapcap', '9',
                                  'Catalyst master stock (aluminum block)')
    vials = ctx.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap',
                                 '11', 'HPLC vial rack')

    # ---- instruments --------------------------------------------------
    p300 = ctx.load_instrument('p300_single_gen2', 'left', tip_racks=[tips300])
    p20 = ctx.load_instrument('p20_single_gen2', 'right', tip_racks=[tips20])

    ctx.comment('PHASE-9 self-driving lab batch: ' + str(len(CONDITIONS)) + ' conditions + blank')
    p300.home()

    # Liquid-class flow rates (uL/s), calibrated for non-aqueous organics
    FLOWS = {
        'dcm': (35, 50),
        'toluene': (70, 90),
        'dmso': (12, 18),
        'meoh': (55, 75),
        'etoh': (55, 75),
        'hplc_diluent': (80, 100),
    }


    # ---- experiment 1: CHAMP-1 -> plate A1 ----
    T_SET = 30.8
    assert 4.0 <= T_SET <= 95.0, 'temperature outside Temperature Module range'
    temp_mod.set_temperature(T_SET)
    temp_mod.await_temperature(T_SET)
    v_cat_ul = 32.24  # catalyst dose at 25 mM
    v_mix_ul = 107.76          # blend makeup volume
    phi_tol  = 0.2908           # toluene volume fraction
    v_tol = v_mix_ul * phi_tol
    v_dcm = v_mix_ul - v_tol

    # (1) solvent blend into the reaction well (instrument chosen by volume)
    if v_tol >= 30.0:
        p300.pick_up_tip()
        _dose(p300, solv.wells_by_name()[SOLVENT_WELLS['toluene']],
              plate.wells_by_name()[CONDITIONS[0]['well']], v_tol, FLOWS['toluene'][0], FLOWS['toluene'][1], prewet=2)
        p300.drop_tip()
    elif v_tol >= 5.0:
        p20.pick_up_tip()
        _dose(p20, solv.wells_by_name()[SOLVENT_WELLS['toluene']],
              plate.wells_by_name()[CONDITIONS[0]['well']], v_tol, FLOWS['toluene'][0], FLOWS['toluene'][1], prewet=2)
        p20.drop_tip()
    else:
        ctx.comment('toluene fraction below P20 minimum -> documented dosing deviation')
    if v_dcm >= 30.0:
        p300.pick_up_tip()
        _dose(p300, solv.wells_by_name()[SOLVENT_WELLS['dcm']],
              plate.wells_by_name()[CONDITIONS[0]['well']], v_dcm, FLOWS['dcm'][0], FLOWS['dcm'][1], prewet=2)
        p300.drop_tip()
    elif v_dcm >= 5.0:
        p20.pick_up_tip()
        _dose(p20, solv.wells_by_name()[SOLVENT_WELLS['dcm']],
              plate.wells_by_name()[CONDITIONS[0]['well']], v_dcm, FLOWS['dcm'][0], FLOWS['dcm'][1], prewet=2)
        p20.drop_tip()
    else:
        ctx.comment('dcm fraction below P20 minimum -> documented dosing deviation')
    p20.pick_up_tip()
    _dose(p20, reag.wells_by_name()[REAGENT_WELLS['substrate']],
          plate.wells_by_name()[CONDITIONS[0]['well']], 40.0,
          FLOWS['toluene'][0], FLOWS['toluene'][1], prewet=1)
    # (2) catalyst stock (viscous-free toluene class, volumetric accuracy)
    _dose(p20, reag.wells_by_name()[REAGENT_WELLS['catalyst']],
          plate.wells_by_name()[CONDITIONS[0]['well']], v_cat_ul,
          FLOWS['toluene'][0] * 0.8, FLOWS['toluene'][1] * 0.8, prewet=2)
    p20.mix(8, 18, plate.wells_by_name()[CONDITIONS[0]['well']])
    p20.drop_tip()

    # (3) reaction hold  t = 20.80 h  (robot idles; temp module holds)
    ctx.delay(minutes=1248.3)

    # (4) quench: MeOH + internal standard (viscous DMSO IS class)
    p20.pick_up_tip()
    q_well = plate.wells_by_name()[CONDITIONS[0]['well']]
    _dose(p20, solv.wells_by_name()[SOLVENT_WELLS['meoh']], q_well, 15.0,
          FLOWS['meoh'][0], FLOWS['meoh'][1], prewet=1)
    _dose(p20, reag.wells_by_name()[REAGENT_WELLS['is_stock']], q_well, 5.0,
          FLOWS['dmso'][0], FLOWS['dmso'][1], prewet=3)
    p20.mix(5, 18, q_well)
    p20.drop_tip()

    # (5) HPLC serial dilution: 5 uL reaction + 195 uL diluent (1:40),
    #     then 10 uL + 190 uL diluent (1:20) -> 1:800 analytical dilution
    v1 = vials.wells_by_name()['A1']
    v2 = vials.wells_by_name()['B1']
    p300.pick_up_tip()
    _dose(p300, solv.wells_by_name()[SOLVENT_WELLS['hplc_diluent']], v1, 195.0,
          FLOWS['hplc_diluent'][0], FLOWS['hplc_diluent'][1])
    p300.drop_tip()
    p20.pick_up_tip()
    p20.aspirate(5.0, q_well.bottom(2))
    p20.air_gap(AIR_GAP_UL)
    p20.dispense(5.0 + AIR_GAP_UL, v1.bottom(3))
    p20.blow_out(v1.top(-2))
    p20.touch_tip(v1)
    p20.mix(8, 18, v1)
    p20.drop_tip()
    p300.pick_up_tip()
    _dose(p300, solv.wells_by_name()[SOLVENT_WELLS['hplc_diluent']], v2, 190.0,
          FLOWS['hplc_diluent'][0], FLOWS['hplc_diluent'][1])
    p300.drop_tip()
    p20.pick_up_tip()
    p20.aspirate(10.0, v1.bottom(3))
    p20.air_gap(AIR_GAP_UL)
    p20.dispense(10.0 + AIR_GAP_UL, v2.bottom(3))
    p20.blow_out(v2.top(-2))
    p20.touch_tip(v2)
    p20.mix(8, 18, v2)
    p20.drop_tip()

    # (6) EtOH tip-conditioning rinse between organic experiments
    p20.pick_up_tip()
    p20.flow_rate.aspirate = FLOWS['etoh'][0]
    p20.flow_rate.dispense = FLOWS['etoh'][1]
    p20.aspirate(20.0, solv.wells_by_name()[SOLVENT_WELLS['etoh']].bottom(2))
    p20.dispense(20.0, solv.wells_by_name()[SOLVENT_WELLS['etoh']].bottom(2))
    p20.drop_tip()

    # ---- experiment 2: CHAMP-2 -> plate B1 ----
    T_SET = 45.2
    assert 4.0 <= T_SET <= 95.0, 'temperature outside Temperature Module range'
    temp_mod.set_temperature(T_SET)
    temp_mod.await_temperature(T_SET)
    v_cat_ul = 8.56  # catalyst dose at 25 mM
    v_mix_ul = 131.44          # blend makeup volume
    phi_tol  = 0.7507           # toluene volume fraction
    v_tol = v_mix_ul * phi_tol
    v_dcm = v_mix_ul - v_tol

    # (1) solvent blend into the reaction well (instrument chosen by volume)
    if v_tol >= 30.0:
        p300.pick_up_tip()
        _dose(p300, solv.wells_by_name()[SOLVENT_WELLS['toluene']],
              plate.wells_by_name()[CONDITIONS[1]['well']], v_tol, FLOWS['toluene'][0], FLOWS['toluene'][1], prewet=2)
        p300.drop_tip()
    elif v_tol >= 5.0:
        p20.pick_up_tip()
        _dose(p20, solv.wells_by_name()[SOLVENT_WELLS['toluene']],
              plate.wells_by_name()[CONDITIONS[1]['well']], v_tol, FLOWS['toluene'][0], FLOWS['toluene'][1], prewet=2)
        p20.drop_tip()
    else:
        ctx.comment('toluene fraction below P20 minimum -> documented dosing deviation')
    if v_dcm >= 30.0:
        p300.pick_up_tip()
        _dose(p300, solv.wells_by_name()[SOLVENT_WELLS['dcm']],
              plate.wells_by_name()[CONDITIONS[1]['well']], v_dcm, FLOWS['dcm'][0], FLOWS['dcm'][1], prewet=2)
        p300.drop_tip()
    elif v_dcm >= 5.0:
        p20.pick_up_tip()
        _dose(p20, solv.wells_by_name()[SOLVENT_WELLS['dcm']],
              plate.wells_by_name()[CONDITIONS[1]['well']], v_dcm, FLOWS['dcm'][0], FLOWS['dcm'][1], prewet=2)
        p20.drop_tip()
    else:
        ctx.comment('dcm fraction below P20 minimum -> documented dosing deviation')
    p20.pick_up_tip()
    _dose(p20, reag.wells_by_name()[REAGENT_WELLS['substrate']],
          plate.wells_by_name()[CONDITIONS[1]['well']], 40.0,
          FLOWS['toluene'][0], FLOWS['toluene'][1], prewet=1)
    # (2) catalyst stock (viscous-free toluene class, volumetric accuracy)
    _dose(p20, reag.wells_by_name()[REAGENT_WELLS['catalyst']],
          plate.wells_by_name()[CONDITIONS[1]['well']], v_cat_ul,
          FLOWS['toluene'][0] * 0.8, FLOWS['toluene'][1] * 0.8, prewet=2)
    p20.mix(8, 18, plate.wells_by_name()[CONDITIONS[1]['well']])
    p20.drop_tip()

    # (3) reaction hold  t = 23.82 h  (robot idles; temp module holds)
    ctx.delay(minutes=1429.4)

    # (4) quench: MeOH + internal standard (viscous DMSO IS class)
    p20.pick_up_tip()
    q_well = plate.wells_by_name()[CONDITIONS[1]['well']]
    _dose(p20, solv.wells_by_name()[SOLVENT_WELLS['meoh']], q_well, 15.0,
          FLOWS['meoh'][0], FLOWS['meoh'][1], prewet=1)
    _dose(p20, reag.wells_by_name()[REAGENT_WELLS['is_stock']], q_well, 5.0,
          FLOWS['dmso'][0], FLOWS['dmso'][1], prewet=3)
    p20.mix(5, 18, q_well)
    p20.drop_tip()

    # (5) HPLC serial dilution: 5 uL reaction + 195 uL diluent (1:40),
    #     then 10 uL + 190 uL diluent (1:20) -> 1:800 analytical dilution
    v1 = vials.wells_by_name()['C1']
    v2 = vials.wells_by_name()['D1']
    p300.pick_up_tip()
    _dose(p300, solv.wells_by_name()[SOLVENT_WELLS['hplc_diluent']], v1, 195.0,
          FLOWS['hplc_diluent'][0], FLOWS['hplc_diluent'][1])
    p300.drop_tip()
    p20.pick_up_tip()
    p20.aspirate(5.0, q_well.bottom(2))
    p20.air_gap(AIR_GAP_UL)
    p20.dispense(5.0 + AIR_GAP_UL, v1.bottom(3))
    p20.blow_out(v1.top(-2))
    p20.touch_tip(v1)
    p20.mix(8, 18, v1)
    p20.drop_tip()
    p300.pick_up_tip()
    _dose(p300, solv.wells_by_name()[SOLVENT_WELLS['hplc_diluent']], v2, 190.0,
          FLOWS['hplc_diluent'][0], FLOWS['hplc_diluent'][1])
    p300.drop_tip()
    p20.pick_up_tip()
    p20.aspirate(10.0, v1.bottom(3))
    p20.air_gap(AIR_GAP_UL)
    p20.dispense(10.0 + AIR_GAP_UL, v2.bottom(3))
    p20.blow_out(v2.top(-2))
    p20.touch_tip(v2)
    p20.mix(8, 18, v2)
    p20.drop_tip()

    # (6) EtOH tip-conditioning rinse between organic experiments
    p20.pick_up_tip()
    p20.flow_rate.aspirate = FLOWS['etoh'][0]
    p20.flow_rate.dispense = FLOWS['etoh'][1]
    p20.aspirate(20.0, solv.wells_by_name()[SOLVENT_WELLS['etoh']].bottom(2))
    p20.dispense(20.0, solv.wells_by_name()[SOLVENT_WELLS['etoh']].bottom(2))
    p20.drop_tip()

    # ---- experiment 3: CHAMP-3 -> plate C1 ----
    T_SET = 60.8
    assert 4.0 <= T_SET <= 95.0, 'temperature outside Temperature Module range'
    temp_mod.set_temperature(T_SET)
    temp_mod.await_temperature(T_SET)
    v_cat_ul = 4.15  # catalyst dose at 25 mM
    v_mix_ul = 135.85          # blend makeup volume
    phi_tol  = 0.9784           # toluene volume fraction
    v_tol = v_mix_ul * phi_tol
    v_dcm = v_mix_ul - v_tol

    # (1) solvent blend into the reaction well (instrument chosen by volume)
    if v_tol >= 30.0:
        p300.pick_up_tip()
        _dose(p300, solv.wells_by_name()[SOLVENT_WELLS['toluene']],
              plate.wells_by_name()[CONDITIONS[2]['well']], v_tol, FLOWS['toluene'][0], FLOWS['toluene'][1], prewet=2)
        p300.drop_tip()
    elif v_tol >= 5.0:
        p20.pick_up_tip()
        _dose(p20, solv.wells_by_name()[SOLVENT_WELLS['toluene']],
              plate.wells_by_name()[CONDITIONS[2]['well']], v_tol, FLOWS['toluene'][0], FLOWS['toluene'][1], prewet=2)
        p20.drop_tip()
    else:
        ctx.comment('toluene fraction below P20 minimum -> documented dosing deviation')
    if v_dcm >= 30.0:
        p300.pick_up_tip()
        _dose(p300, solv.wells_by_name()[SOLVENT_WELLS['dcm']],
              plate.wells_by_name()[CONDITIONS[2]['well']], v_dcm, FLOWS['dcm'][0], FLOWS['dcm'][1], prewet=2)
        p300.drop_tip()
    elif v_dcm >= 5.0:
        p20.pick_up_tip()
        _dose(p20, solv.wells_by_name()[SOLVENT_WELLS['dcm']],
              plate.wells_by_name()[CONDITIONS[2]['well']], v_dcm, FLOWS['dcm'][0], FLOWS['dcm'][1], prewet=2)
        p20.drop_tip()
    else:
        ctx.comment('dcm fraction below P20 minimum -> documented dosing deviation')
    p20.pick_up_tip()
    _dose(p20, reag.wells_by_name()[REAGENT_WELLS['substrate']],
          plate.wells_by_name()[CONDITIONS[2]['well']], 40.0,
          FLOWS['toluene'][0], FLOWS['toluene'][1], prewet=1)
    # (2) catalyst stock (viscous-free toluene class, volumetric accuracy)
    _dose(p20, reag.wells_by_name()[REAGENT_WELLS['catalyst']],
          plate.wells_by_name()[CONDITIONS[2]['well']], v_cat_ul,
          FLOWS['toluene'][0] * 0.8, FLOWS['toluene'][1] * 0.8, prewet=2)
    p20.mix(8, 18, plate.wells_by_name()[CONDITIONS[2]['well']])
    p20.drop_tip()

    # (3) reaction hold  t = 15.93 h  (robot idles; temp module holds)
    ctx.delay(minutes=955.9)

    # (4) quench: MeOH + internal standard (viscous DMSO IS class)
    p20.pick_up_tip()
    q_well = plate.wells_by_name()[CONDITIONS[2]['well']]
    _dose(p20, solv.wells_by_name()[SOLVENT_WELLS['meoh']], q_well, 15.0,
          FLOWS['meoh'][0], FLOWS['meoh'][1], prewet=1)
    _dose(p20, reag.wells_by_name()[REAGENT_WELLS['is_stock']], q_well, 5.0,
          FLOWS['dmso'][0], FLOWS['dmso'][1], prewet=3)
    p20.mix(5, 18, q_well)
    p20.drop_tip()

    # (5) HPLC serial dilution: 5 uL reaction + 195 uL diluent (1:40),
    #     then 10 uL + 190 uL diluent (1:20) -> 1:800 analytical dilution
    v1 = vials.wells_by_name()['A2']
    v2 = vials.wells_by_name()['B2']
    p300.pick_up_tip()
    _dose(p300, solv.wells_by_name()[SOLVENT_WELLS['hplc_diluent']], v1, 195.0,
          FLOWS['hplc_diluent'][0], FLOWS['hplc_diluent'][1])
    p300.drop_tip()
    p20.pick_up_tip()
    p20.aspirate(5.0, q_well.bottom(2))
    p20.air_gap(AIR_GAP_UL)
    p20.dispense(5.0 + AIR_GAP_UL, v1.bottom(3))
    p20.blow_out(v1.top(-2))
    p20.touch_tip(v1)
    p20.mix(8, 18, v1)
    p20.drop_tip()
    p300.pick_up_tip()
    _dose(p300, solv.wells_by_name()[SOLVENT_WELLS['hplc_diluent']], v2, 190.0,
          FLOWS['hplc_diluent'][0], FLOWS['hplc_diluent'][1])
    p300.drop_tip()
    p20.pick_up_tip()
    p20.aspirate(10.0, v1.bottom(3))
    p20.air_gap(AIR_GAP_UL)
    p20.dispense(10.0 + AIR_GAP_UL, v2.bottom(3))
    p20.blow_out(v2.top(-2))
    p20.touch_tip(v2)
    p20.mix(8, 18, v2)
    p20.drop_tip()

    # (6) EtOH tip-conditioning rinse between organic experiments
    p20.pick_up_tip()
    p20.flow_rate.aspirate = FLOWS['etoh'][0]
    p20.flow_rate.dispense = FLOWS['etoh'][1]
    p20.aspirate(20.0, solv.wells_by_name()[SOLVENT_WELLS['etoh']].bottom(2))
    p20.dispense(20.0, solv.wells_by_name()[SOLVENT_WELLS['etoh']].bottom(2))
    p20.drop_tip()

    # ---- experiment 4: CHAMP-4 -> plate D1 ----
    T_SET = 40.3
    assert 4.0 <= T_SET <= 95.0, 'temperature outside Temperature Module range'
    temp_mod.set_temperature(T_SET)
    temp_mod.await_temperature(T_SET)
    v_cat_ul = 38.07  # catalyst dose at 25 mM
    v_mix_ul = 101.93          # blend makeup volume
    phi_tol  = 0.9870           # toluene volume fraction
    v_tol = v_mix_ul * phi_tol
    v_dcm = v_mix_ul - v_tol

    # (1) solvent blend into the reaction well (instrument chosen by volume)
    if v_tol >= 30.0:
        p300.pick_up_tip()
        _dose(p300, solv.wells_by_name()[SOLVENT_WELLS['toluene']],
              plate.wells_by_name()[CONDITIONS[3]['well']], v_tol, FLOWS['toluene'][0], FLOWS['toluene'][1], prewet=2)
        p300.drop_tip()
    elif v_tol >= 5.0:
        p20.pick_up_tip()
        _dose(p20, solv.wells_by_name()[SOLVENT_WELLS['toluene']],
              plate.wells_by_name()[CONDITIONS[3]['well']], v_tol, FLOWS['toluene'][0], FLOWS['toluene'][1], prewet=2)
        p20.drop_tip()
    else:
        ctx.comment('toluene fraction below P20 minimum -> documented dosing deviation')
    if v_dcm >= 30.0:
        p300.pick_up_tip()
        _dose(p300, solv.wells_by_name()[SOLVENT_WELLS['dcm']],
              plate.wells_by_name()[CONDITIONS[3]['well']], v_dcm, FLOWS['dcm'][0], FLOWS['dcm'][1], prewet=2)
        p300.drop_tip()
    elif v_dcm >= 5.0:
        p20.pick_up_tip()
        _dose(p20, solv.wells_by_name()[SOLVENT_WELLS['dcm']],
              plate.wells_by_name()[CONDITIONS[3]['well']], v_dcm, FLOWS['dcm'][0], FLOWS['dcm'][1], prewet=2)
        p20.drop_tip()
    else:
        ctx.comment('dcm fraction below P20 minimum -> documented dosing deviation')
    p20.pick_up_tip()
    _dose(p20, reag.wells_by_name()[REAGENT_WELLS['substrate']],
          plate.wells_by_name()[CONDITIONS[3]['well']], 40.0,
          FLOWS['toluene'][0], FLOWS['toluene'][1], prewet=1)
    # (2) catalyst stock (viscous-free toluene class, volumetric accuracy)
    _dose(p20, reag.wells_by_name()[REAGENT_WELLS['catalyst']],
          plate.wells_by_name()[CONDITIONS[3]['well']], v_cat_ul,
          FLOWS['toluene'][0] * 0.8, FLOWS['toluene'][1] * 0.8, prewet=2)
    p20.mix(8, 18, plate.wells_by_name()[CONDITIONS[3]['well']])
    p20.drop_tip()

    # (3) reaction hold  t = 22.03 h  (robot idles; temp module holds)
    ctx.delay(minutes=1321.6)

    # (4) quench: MeOH + internal standard (viscous DMSO IS class)
    p20.pick_up_tip()
    q_well = plate.wells_by_name()[CONDITIONS[3]['well']]
    _dose(p20, solv.wells_by_name()[SOLVENT_WELLS['meoh']], q_well, 15.0,
          FLOWS['meoh'][0], FLOWS['meoh'][1], prewet=1)
    _dose(p20, reag.wells_by_name()[REAGENT_WELLS['is_stock']], q_well, 5.0,
          FLOWS['dmso'][0], FLOWS['dmso'][1], prewet=3)
    p20.mix(5, 18, q_well)
    p20.drop_tip()

    # (5) HPLC serial dilution: 5 uL reaction + 195 uL diluent (1:40),
    #     then 10 uL + 190 uL diluent (1:20) -> 1:800 analytical dilution
    v1 = vials.wells_by_name()['C2']
    v2 = vials.wells_by_name()['D2']
    p300.pick_up_tip()
    _dose(p300, solv.wells_by_name()[SOLVENT_WELLS['hplc_diluent']], v1, 195.0,
          FLOWS['hplc_diluent'][0], FLOWS['hplc_diluent'][1])
    p300.drop_tip()
    p20.pick_up_tip()
    p20.aspirate(5.0, q_well.bottom(2))
    p20.air_gap(AIR_GAP_UL)
    p20.dispense(5.0 + AIR_GAP_UL, v1.bottom(3))
    p20.blow_out(v1.top(-2))
    p20.touch_tip(v1)
    p20.mix(8, 18, v1)
    p20.drop_tip()
    p300.pick_up_tip()
    _dose(p300, solv.wells_by_name()[SOLVENT_WELLS['hplc_diluent']], v2, 190.0,
          FLOWS['hplc_diluent'][0], FLOWS['hplc_diluent'][1])
    p300.drop_tip()
    p20.pick_up_tip()
    p20.aspirate(10.0, v1.bottom(3))
    p20.air_gap(AIR_GAP_UL)
    p20.dispense(10.0 + AIR_GAP_UL, v2.bottom(3))
    p20.blow_out(v2.top(-2))
    p20.touch_tip(v2)
    p20.mix(8, 18, v2)
    p20.drop_tip()

    # (6) EtOH tip-conditioning rinse between organic experiments
    p20.pick_up_tip()
    p20.flow_rate.aspirate = FLOWS['etoh'][0]
    p20.flow_rate.dispense = FLOWS['etoh'][1]
    p20.aspirate(20.0, solv.wells_by_name()[SOLVENT_WELLS['etoh']].bottom(2))
    p20.dispense(20.0, solv.wells_by_name()[SOLVENT_WELLS['etoh']].bottom(2))
    p20.drop_tip()

    # ---- shutdown ------------------------------------------------------
    temp_mod.deactivate()
    p300.home()
    ctx.comment('Batch complete: HPLC-ready vials on slot 11.')
