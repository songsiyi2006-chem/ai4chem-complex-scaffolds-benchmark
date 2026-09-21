"""Small independent invariants; no claims of DFT/CFD or industrial accuracy."""
import contextlib
from dataclasses import asdict
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from common import Config, F, R, MW_KG, dump, module, sha, validate_micro_artifact
import pipeline

MICRO = module("01_micro_surface.py")
MESO = module("02_meso_microkinetics.py")
MACRO = module("03_macro_flow_cfd.py")
TEA = module("04_tea_lca_model.py")
AGENT = module("05_agent_closed_loop.py")


class PhysicsTests(unittest.TestCase):
    def test_reference_potential_sign_and_temperature(self):
        self.assertAlmostEqual(MICRO.rhe_to_she(0,7,298.15),-.4141,places=3)
        self.assertEqual(MICRO.rhe_to_she(-.7,0,298.15),-.7)

    def test_dimension_and_booleans_fail_closed(self):
        for kw in ({"cells":True},{"cells":3},{"material":"Ni"},{"temperature_K":float("nan")},
                   {"recovery_fraction":0},{"velocity_m_s":0},{"discount_rate":-1},
                   {"co2_inlet_mol_m3":100},{"height_m":1},{"batch_size":4}):
            with self.subTest(kw=kw), self.assertRaises(ValueError):
                Config().with_changes(**kw)
        with self.assertRaises(TypeError):
            Config(unrecognized_unit=1)

    def test_rates_are_finite_and_explicit_synthetic(self):
        p = MICRO.packet(Config())
        self.assertFalse(p["dft_executed"])
        self.assertIsNone(p["adsorption_energies_eV"])
        self.assertTrue(all(0 < k < 1e14 for k in p["effective_constants_s"].values()))

    def test_zero_carbon_does_not_create_carbon_products(self):
        s = MESO.solve(Config(),0,True)
        for p in ("CO","C2H4","HCOOH"):
            self.assertEqual(s["flux_mol_m2_s"][p],0)
        self.assertEqual(s["FE"]["H2"],1.)

    def test_ode_matches_quadratic_root_and_charge(self):
        c = Config()
        for material in ("Cu","Ag"):
            s = MESO.solve(c.with_changes(material=material),15.,True)
            self.assertLess(s["ode_validation"]["analytic_error"],1e-7)
            self.assertAlmostEqual(sum(s["FE"].values()),1.)
            r = s["flux_mol_m2_s"]
            independently = F*(2*r["CO"]+2*r["HCOOH"]+12*r["C2H4"]+2*r["H2"])
            self.assertAlmostEqual(independently,s["current_A_m2"],places=8)

    def test_macro_carbon_heat_and_grid_refinement(self):
        c = Config()
        a,b = MACRO.channel(c),MACRO.channel(c.with_changes(cells=24))
        self.assertLess(abs(a["carbon_residual_mol_s"]),1e-12)
        self.assertLess(abs(a["heat_residual_W"]),1e-9)
        self.assertAlmostEqual(a["power_W"]*c.heat_fraction,
                               a["heat_removed_W"]+a["sensible_heat_W"],places=10)
        self.assertLess(abs(a["current_A"]/b["current_A"]-1),.1)
        self.assertTrue(all(r["dissolved_CO2_mol_m3"]>=0 and r["gas_CO2_mol_s"]>=-1e-15 for r in a["cells"]))
        self.assertLessEqual(a["carbon_utilization"],1.)

    def test_pressure_scales_with_velocity_in_laminar_closure(self):
        c=Config()
        a,b=MACRO.channel(c),MACRO.channel(c.with_changes(velocity_m_s=.1))
        self.assertAlmostEqual(b["pressure_drop_Pa"]/a["pressure_drop_Pa"],2.)

    def test_unsupported_envelope_rejected(self):
        with self.assertRaises(ValueError):
            MACRO.channel(Config(max_current_A_m2=1.))
        with self.assertRaises(ValueError):
            MESO.solve(Config(),-1.)


class EngineeringTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.flow=MACRO.channel(Config())

    def test_zero_degradation_zero_discount_analytic_product(self):
        c=Config(degradation_per_operating_h=0,discount_rate=0)
        t=TEA.assess(c,self.flow)
        expected=(self.flow["products_mol_s"][c.target]*c.stack_kW*c.stacks*1000/
                  self.flow["power_W"]*MW_KG[c.target]*3600*c.recovery_fraction*
                  8760*c.capacity_factor*c.years)
        self.assertAlmostEqual(t["product_lifetime_kg"],expected,places=6)
        self.assertEqual(t["total_replacements"],0)
        self.assertAlmostEqual(t["discounted_product_kg"],expected,places=6)
        self.assertIsNone(t["LCOE_USD_kWh"])

    def test_discounted_cost_components_reconcile(self):
        t=TEA.assess(Config(),self.flow)
        self.assertAlmostEqual(sum(t["discounted_cost_USD"].values())/t["discounted_product_kg"],
                               t["LCO_product_USD_kg"])
        self.assertTrue(all(x >=0 for x in t["GHG_inventory_kgCO2e"].values()))

    def test_no_false_economies_of_scale(self):
        c=Config()
        a,b=TEA.assess(c,self.flow,1),TEA.assess(c,self.flow,10)
        self.assertAlmostEqual(b["product_lifetime_kg"]/a["product_lifetime_kg"],10.)
        self.assertAlmostEqual(a["LCO_product_USD_kg"],b["LCO_product_USD_kg"])
        self.assertAlmostEqual(a["screening_GHG_kgCO2e_kg"],b["screening_GHG_kgCO2e_kg"])

    def test_events_reconcile_scheduled_time_and_exact_cycle_integral(self):
        c=Config()
        rows,trace=TEA.lifetime(c)
        scheduled=8760*c.capacity_factor*c.years
        self.assertAlmostEqual(sum(r["operating_h"]+r["downtime_h"] for r in rows),scheduled)
        self.assertGreater(sum(r["replacements"] for r in rows),0)
        self.assertTrue(all(r["effective_product_h"]<=r["operating_h"]+1e-8 for r in rows))
        self.assertTrue(all(c.replacement_threshold-1e-8<=r["activity"]<=1 for r in trace))
        nochange=Config(degradation_per_operating_h=0)
        self.assertGreater(TEA.assess(c,self.flow)["LCO_product_USD_kg"],
                           TEA.assess(nochange,self.flow)["LCO_product_USD_kg"])

    def test_zero_electricity_price_and_emissions_not_negative_credit(self):
        t=TEA.assess(Config(electricity_USD_kWh=0,grid_kgCO2e_kWh=0),self.flow)
        self.assertEqual(t["discounted_cost_USD"]["electricity"],0)
        self.assertEqual(t["GHG_inventory_kgCO2e"]["electricity"],0)
        self.assertGreater(t["screening_GHG_kgCO2e_kg"],0)


class OptimizerTests(unittest.TestCase):
    def test_pareto_direction_ties_and_nonfinite(self):
        np.testing.assert_array_equal(AGENT.pareto_mask([[1,3],[2,2],[3,3],[1,3]]),
                                      [True,True,False,True])
        with self.assertRaises(ValueError):
            AGENT.pareto_mask([[float("nan"),1]])

    def test_gp_and_ei_are_finite_nonnegative(self):
        x=np.array([[0.],[1.]])
        mean,sd=AGENT.gp_predict(x,np.array([2.,3.]),np.array([[.5],[.7]]))
        self.assertTrue(np.isfinite(mean).all())
        self.assertTrue((sd>0).all())
        self.assertTrue((AGENT.expected_improvement(2,mean,sd)>=0).all())

    def test_closed_loop_causality_duplicates_and_observed_front(self):
        # Cheap deterministic oracle: assesses controller plumbing, not chemistry.
        def assess(c,flow):
            return {"LCO_product_USD_kg":2+(c.potential_RHE_V+.6)**2,
                    "screening_GHG_kgCO2e_kg":1+c.velocity_m_s,
                    "initial_recovered_kg_h":3+c.ligand_shift_eV}
        def wrapped(c):
            return dict(c=c,FE={c.target:.7},carbon_residual_mol_s=0)
        c=Config(initial_samples=3,rounds=2,batch_size=3)
        with patch.object(AGENT.macro,"channel",wrapped),patch.object(AGENT.tea,"assess",assess):
            result=AGENT.optimize(c)
            repeat=AGENT.optimize(c)
        self.assertEqual(result,repeat)
        rows=result["observations"]
        self.assertEqual(len(rows),9)
        self.assertEqual(len({r["candidate_id"] for r in rows}),9)
        rounds={r["candidate_id"]:r["round"] for r in rows}
        for p in result["proposal_ledger"]:
            if p.get("event")=="proposal":
                self.assertTrue(all(rounds[i]<p["round"] for i in p["training_candidate_ids"]))
        self.assertTrue({r["candidate_id"] for r in result["pareto"]}<={r["candidate_id"] for r in rows})

    def test_zero_burden_scenario_does_not_break_acquisition(self):
        c=Config(initial_samples=2,rounds=1,batch_size=3)
        with patch.object(AGENT.macro,"channel",return_value={"FE":{"C2H4":.5},"carbon_residual_mol_s":0}), \
             patch.object(AGENT.tea,"assess",return_value={"LCO_product_USD_kg":1.,"screening_GHG_kgCO2e_kg":0.,
                                                         "initial_recovered_kg_h":1.}):
            result=AGENT.optimize(c)
        self.assertEqual(result["evaluations"],5)
        self.assertTrue(all(np.isfinite(r["acquisition"]) for r in result["proposal_ledger"] if r.get("event")=="proposal"))


class EvidenceTests(unittest.TestCase):
    def test_atomic_json_refuses_nonfinite(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/"x.json"
            with self.assertRaises(ValueError): dump(p,{"x":float("nan")})
            self.assertFalse(p.exists())

    def test_workspace_no_clobber_source_and_mutated_config_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp)/"run"
            pipeline.prepare(out,Config())
            with self.assertRaises(FileExistsError): pipeline.prepare(out,Config())
            dump(out/"resolved_config.json",asdict(Config(potential_RHE_V=-.6)))
            with self.assertRaisesRegex(ValueError,"changed"): pipeline.execute(out)
        with self.assertRaises(ValueError): pipeline.prepare(pipeline.PROJECT/"results/new",Config())

    def test_failure_is_recorded_and_cannot_be_reexecuted(self):
        with tempfile.TemporaryDirectory() as temp:
            out=pipeline.prepare(Path(temp)/"run",Config())
            with patch.object(pipeline.MICRO,"packet",side_effect=ArithmeticError("injected")):
                with self.assertRaises(ArithmeticError): pipeline.execute(out)
            self.assertEqual(json.loads((out/"run_status.json").read_text())["status"],"FAILED")
            self.assertTrue((out/"failure.txt").is_file())
            self.assertTrue((out/"artifact_manifest.json").is_file())
            with self.assertRaises(ValueError): pipeline.execute(out)

    def test_micro_evidence_hash_units_and_path_guard(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); (root/"raw.txt").write_text("computed output")
            data=dict(schema_version="phase30.micro.v1",evidence_type="computed",
                      method="grand canonical",functional="PBE",solvation="external",
                      potential_reference="SHE",temperature_K=298.15,surface_id="Cu(111)",
                      geometry_sha256="0"*64,convergence={"scf":True},
                      records=[dict(reaction_id="ads",barrier_eV=.8,potential_SHE_V=-1.,unit="eV")],
                      raw_files=[dict(path="raw.txt",sha256=sha(root/"raw.txt"))])
            p=root/"input.json"; dump(p,data)
            self.assertFalse(validate_micro_artifact(p)["chemistry_validated"])
            data["raw_files"][0]["path"]="../escape.txt"; dump(p,data)
            with self.assertRaises(ValueError): validate_micro_artifact(p)
            data["records"][0]["unit"]="kcal/mol"; dump(p,data)
            with self.assertRaises(ValueError): validate_micro_artifact(p)

    def test_production_mode_never_silently_runs_demo(self):
        spec=importlib.util.spec_from_file_location("phase30_cli",pipeline.PROJECT/"run.py")
        cli=importlib.util.module_from_spec(spec); spec.loader.exec_module(cli)
        with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit) as error:
            cli.main(["--mode","production"])
        self.assertEqual(error.exception.code,2)


if __name__=="__main__": unittest.main()
