import unittest
import numpy as np
from physics import simulate,observations,product_kg,kaplan_meier,no_batch_leakage,FARADAY_C_MOL
from control import schedule,choose_action,account

class PhysicsTests(unittest.TestCase):
    def test_mass_conservation_nonnegative(self):
        for kir,kdes,kc in [(0,0,.12),(.003,.0015,.12),(.003,.0015,.0045),(.4,.6,0)]:
            states=simulate(np.tile([0,.75,1.25],20),kir=kir,kdes=kdes,clearance=kc)
            np.testing.assert_allclose(states.sum(axis=1),1,atol=1e-13)
            self.assertGreaterEqual(states.min(),-1e-12)
    def test_voltage_route_nonidentifiability(self):
        load=np.linspace(.75,1.25,30)
        np.testing.assert_allclose(observations(load,kir=.004,kdes=.001)[:,0],observations(load,kir=.001,kdes=.004)[:,0],atol=1e-14)
        self.assertGreater(np.max(abs(observations(load,kir=.004,kdes=.001)[:,2]-observations(load,kir=.001,kdes=.004)[:,2])),.01)
    def test_faraday_charge(self):
        # 2 F coulombs produces one mole Cl2 at FE=1.
        mol=product_kg([1],dt_h=2*FARADAY_C_MOL/3600,area_cm2=1,j_ref=1,fe=1)[0]/.07090
        self.assertAlmostEqual(mol,1)
    def test_censoring(self):
        rows=kaplan_meier([1,2,2,3],[True,True,False,False])
        self.assertAlmostEqual(rows[-1][-1],.5)
        self.assertEqual(rows[1][1:4],(3,1,1))
    def test_batch_leakage_rejected(self):
        with self.assertRaises(ValueError): no_batch_leakage([{"batch_id":"a","split":"train"},{"batch_id":"a","split":"blind"}])
    def test_fair_schedules_and_constraints(self):
        t=np.arange(48); prices=.1+.05*np.sin(t/3); avail=np.where(t%17<3,1,1.25)
        for strategy in ("constant","rule","mpc"):
            load=schedule(strategy,prices,avail)
            self.assertAlmostEqual(sum(load),48)
            self.assertLessEqual(max(abs(np.diff(np.r_[1,load]))),.25+1e-12)
            self.assertLessEqual(max(abs(np.cumsum(load-1))),2+1e-12)
            self.assertTrue(np.all(load<=avail))
            self.assertAlmostEqual(account(load,prices)["qualified_Cl2_kg"],account(np.ones(48),prices)["qualified_Cl2_kg"])
    def test_no_future_price_access(self):
        prices=np.array([.1,.08,.07,.06,.2,.3]); avail=np.full(6,1.25)
        changed=prices.copy(); changed[3:]=[99,99,99]
        np.testing.assert_array_equal(schedule("mpc",prices,avail)[:3],schedule("mpc",changed,avail)[:3])
    def test_jensen_negative_result(self):
        fixed=simulate(np.ones(20))[-1,0]
        dynamic=simulate(np.tile([.75,1.25],10))[-1,0]
        self.assertLess(dynamic,fixed)
    def test_zero_current_no_product(self):
        self.assertEqual(float(product_kg([0,0]).sum()),0)
    def test_unavailable_baseload_rejected(self):
        with self.assertRaises(ValueError): schedule("constant",[.1,.1],[1,.8])
    def test_unsupported_config_rejected(self):
        import json
        from pathlib import Path
        from analysis import validate_supported_config
        cfg=json.loads((Path(__file__).resolve().parents[1]/"configs/local.json").read_text())
        validate_supported_config(cfg)
        cfg["faradaic_efficiency"]=.5
        with self.assertRaises(ValueError): validate_supported_config(cfg)

if __name__=="__main__": unittest.main()
