"""Analytic/unit fixtures are NOT molecular simulation results."""
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from common import EvidenceBlocked, load_module, source_hashes, ROOT
from swarm.memory_ledger import Ledger
from swarm.coordinator import prepare, validate_resume, work
from swarm.sandbox_executor import QuantumExecutionSandbox
from swarm.provider_contract import Response, invoke_reserved

qm = load_module("02_pyscf_qm_engine")
ts = load_module("03_neb_transition_state")
fep = load_module("04_qmmm_openmm_fep")
retro = load_module("05_retrosynthesis_check")
coop = load_module("06_cooperativity_eval")


class LedgerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "queue.sqlite"
        self.ledger = Ledger(self.path, 1000)

    def tearDown(self):
        self.ledger.close()
        self.tmp.cleanup()

    def test_fencing_expired(self):
        self.ledger.add("a", "Alpha")
        a = self.ledger.claim("first", 1, now=0)
        b = self.ledger.claim("second", 10, now=2)
        with self.assertRaises(RuntimeError):
            self.ledger.finish(a, "SUCCEEDED", {}, now=3)
        self.ledger.finish(b, "SUCCEEDED", {}, now=3)

    def test_no_completion_after_expiry(self):
        self.ledger.add("a", "Alpha")
        a = self.ledger.claim("x", 1, now=0)
        with self.assertRaises(RuntimeError):
            self.ledger.finish(a, "SUCCEEDED", {}, now=1)

    def test_atomic_claim(self):
        self.ledger.add("a", "Alpha")
        def worker(n):
            ledger = Ledger(self.path, 1000)
            try:
                return ledger.claim(str(n))
            finally:
                ledger.close()
        with ThreadPoolExecutor(max_workers=2) as pool:
            claims = list(pool.map(worker, range(2)))
        self.assertEqual(sum(c is not None for c in claims), 1)

    def test_dependency_audit_after_block(self):
        self.ledger.add("a", "Alpha")
        self.ledger.add("b", "Omega", ["a"])
        a = self.ledger.claim("worker")
        self.assertIsNone(self.ledger.claim("other"))
        self.ledger.finish(a, "BLOCKED", {"reason":"missing input"})
        self.assertEqual(self.ledger.claim("worker")["id"], "b")

    def test_forward_or_cyclic_dependencies_forbidden(self):
        with self.assertRaises(ValueError):
            self.ledger.add("a", "Alpha", ["b"])
        with self.assertRaises(ValueError):
            self.ledger.add("a", "Alpha", ["a"])

    def test_idempotent_task(self):
        self.ledger.add("a", "Alpha")
        self.ledger.add("a", "Alpha")
        with self.assertRaises(ValueError):
            self.ledger.add("a", "Beta")
        self.assertEqual(len(self.ledger.snapshot()["tasks"]), 1)

    def test_role_budget(self):
        self.ledger.reserve("r", "Alpha", 100)
        with self.assertRaises(RuntimeError):
            self.ledger.reserve("r2", "Alpha", 1)

    def test_budget_reservation_atomic(self):
        def reserve(n):
            ledger = Ledger(self.path, 1000)
            try:
                try:
                    ledger.reserve(str(n), "Alpha", 60)
                    return True
                except RuntimeError:
                    return False
            finally:
                ledger.close()
        with ThreadPoolExecutor(max_workers=2) as pool:
            self.assertEqual(sum(pool.map(reserve, range(2))), 1)

    def test_cached_not_double_counted(self):
        self.ledger.reserve("r", "Alpha", 100)
        self.ledger.settle("r", "receipt", 50, 20, 40)
        self.ledger.settle("r", "receipt", 50, 20, 40)
        self.assertEqual(self.ledger.snapshot()["actual_provider_tokens"], 70)

    def test_duplicate_receipt_rejected(self):
        self.ledger.reserve("a", "Alpha", 40)
        self.ledger.reserve("b", "Alpha", 40)
        self.ledger.settle("a", "same", 10, 10)
        with self.assertRaises(ValueError):
            self.ledger.settle("b", "same", 10, 10)

    def test_unknown_usage_keeps_reservation(self):
        self.ledger.reserve("r", "Alpha", 90)
        self.assertEqual(self.ledger.snapshot()["unresolved_reserved_tokens"], 90)

    def test_overrun_fully_charged(self):
        self.ledger.reserve("r", "Alpha", 90)
        self.ledger.settle("r", "receipt", 1000, 50)
        self.assertEqual(self.ledger.snapshot()["actual_provider_tokens"], 1050)
        with self.assertRaises(RuntimeError):
            self.ledger.reserve("other", "Beta", 1)

    def test_receipt_conflict(self):
        self.ledger.reserve("r", "Alpha", 40)
        self.ledger.settle("r", "id", 10, 10)
        with self.assertRaises(ValueError):
            self.ledger.settle("r", "id", 11, 10)

    def test_cached_subset(self):
        self.ledger.reserve("r", "Alpha", 40)
        with self.assertRaises(ValueError):
            self.ledger.settle("r", "id", 10, 10, 11)

    def test_event_chain_and_immutability(self):
        self.ledger.add("a", "Alpha")
        self.assertEqual(len(self.ledger.verify_events()), 64)
        with self.assertRaises(sqlite3.IntegrityError):
            self.ledger.db.execute("UPDATE events SET hash='wrong'")

    def test_frozen_cap(self):
        with self.assertRaises(ValueError):
            Ledger(self.path, 2000)

    def test_numerical_failure_not_blacklist(self):
        with self.assertRaises(ValueError):
            self.ledger.reject_molecule("a"*64, "SCF_FAILURE", "Omega")

    def test_reviewed_blacklist(self):
        self.ledger.reject_molecule("a"*64, "INVALID_VALENCE", "reviewer")
        self.assertEqual(self.ledger.db.execute("SELECT count(*) FROM blacklist").fetchone()[0], 1)

    def test_provider_once_and_response_recovered(self):
        class FakeProvider:
            calls = 0
            def invoke(self, **kwargs):
                self.calls += 1
                return Response("test_fixture_receipt", 10, 5, 3, {"fixture": True})
        provider = FakeProvider()
        a = invoke_reserved(self.ledger, provider, "r", "Alpha", ["hash"], 50)
        b = invoke_reserved(self.ledger, provider, "r", "Alpha", ["hash"], 50)
        self.assertEqual(a, b)
        self.assertEqual(provider.calls, 1)
        self.assertEqual(self.ledger.snapshot()["actual_provider_tokens"], 15)

    def test_provider_unknown_usage_retained(self):
        class FakeProvider:
            def invoke(self, **kwargs):
                raise TimeoutError("test fixture")
        with self.assertRaises(TimeoutError):
            invoke_reserved(self.ledger, FakeProvider(), "r", "Alpha", ["hash"], 50)
        self.assertEqual(self.ledger.snapshot()["unresolved_reserved_tokens"], 50)

    def test_provider_payload_immutable(self):
        class FakeProvider:
            def invoke(self, **kwargs):
                return Response("test_fixture_receipt", 10, 5, 3, {})
        invoke_reserved(self.ledger, FakeProvider(), "r", "Alpha", ["hash"], 50)
        with self.assertRaises(ValueError):
            invoke_reserved(self.ledger, FakeProvider(), "r", "Alpha", ["different"], 50)


class ScienceTests(unittest.TestCase):
    def test_original_method_blocked(self):
        with self.assertRaises(EvidenceBlocked):
            qm.validate({"xc":"wb97x-d"})

    def test_no_unauthorized_method_fallback(self):
        with self.assertRaises(EvidenceBlocked):
            qm.validate({"xc":"pbe0"})

    def test_physics_fingerprint(self):
        spec = dict(atoms=[["H", [0,0,0]]], charge=0, spin=1, xc="pbe0", basis="def2-tzvp")
        altered = {**spec, "spin":0}
        self.assertNotEqual(qm.fingerprint(spec), qm.fingerprint(altered))
        self.assertEqual(qm.fingerprint(spec), qm.fingerprint({**spec, "max_cycle":999}))

    def test_neb_alone_not_ts(self):
        result = ts.audit_ts(dict(neb_converged=True, frequencies_cm1=[-400, 50, 100]))
        self.assertFalse(result["contract_passed"])

    def test_no_fake_irc(self):
        with self.assertRaises(EvidenceBlocked):
            ts.irc()

    def test_partition(self):
        self.assertTrue(fep.validate_partition(4,[0,1],[2,3],[]))
        with self.assertRaises(ValueError):
            fep.validate_partition(4,[0,1],[1,2,3],[])

    def test_link_atom_review(self):
        with self.assertRaises(EvidenceBlocked):
            fep.validate_partition(4,[0,1],[2,3],[(1,2)])

    def test_ti_analytic_linear(self):
        result = fep.ti([0,.5,1], [0,1,2], np.eye(3) * .04)
        self.assertAlmostEqual(result["delta_g_kj_mol"], 1)
        self.assertAlmostEqual(result["standard_error_kj_mol"], np.sqrt(.015))

    def test_ti_full_covariance(self):
        result = fep.ti([0,.5,1], [1,1,1], np.ones((3,3)) * .04)
        self.assertAlmostEqual(result["standard_error_kj_mol"], .2)

    def test_ti_bad_covariance(self):
        with self.assertRaises(ValueError):
            fep.ti([0,.5,1], [0,1,2], -np.eye(3))

    def test_lambda_endpoints(self):
        for x in ([.1,.5,1], [0,.5,.5,1], [0,float("nan"),1]):
            with self.assertRaises(ValueError):
                fep.lambda_schedule(x)

    def test_fep_constant(self):
        result = fep.forward_fep([5.,5.,5.])
        self.assertAlmostEqual(result["delta_g_kj_mol"], 5)
        self.assertAlmostEqual(result["weight_ess"], 3)
        self.assertFalse(result["convergence_proven"])

    def test_fep_stability(self):
        result = fep.forward_fep([-10000,10000,0])
        self.assertTrue(np.isfinite(result["delta_g_kj_mol"]))

    def test_missing_correction(self):
        with self.assertRaises(EvidenceBlocked):
            fep.binding_cycle(1,2,0,None,0)

    def test_cycle_sign(self):
        self.assertEqual(fep.binding_cycle(5,2,1,2,-1), 5)

    def test_cooperativity_sign(self):
        result = coop.alpha_from_kd(20e-9, 1e-9)
        self.assertAlmostEqual(result["alpha"],20)
        self.assertLess(result["delta_delta_g_kj_mol"], 0)

    def test_covalent_not_kd(self):
        with self.assertRaises(ValueError):
            coop.alpha_from_kd(1,1,reversible=False)

    def test_shared_reference_cancels(self):
        result = coop.alpha_from_free_energies(-1,-2,[[1,1],[1,1]])
        self.assertEqual(result["standard_error_ddg"],0)
        self.assertGreater(result["alpha"],1)

    def test_unknown_route(self):
        self.assertIsNone(retro.audit_route(None)["synthesis_feasible"])

    def test_route_cycle(self):
        with self.assertRaises(ValueError):
            retro.audit_route({"steps":[{"id":"a","parents":["b"]},{"id":"b","parents":["a"]}]})

    def test_route_longest_not_total(self):
        nodes = [{"id":"root"}] + [{"id":str(i),"parents":["root"]} for i in range(20)]
        self.assertEqual(retro.audit_route({"steps":nodes})["longest_linear_steps"],2)

    def test_route_policy(self):
        nodes = [{"id":str(i),"parents":[str(i-1)] if i else []} for i in range(13)]
        self.assertEqual(retro.audit_route({"steps":nodes})["status"],"REJECT_POLICY")

    def test_sandbox_timeout_bound(self):
        with self.assertRaises(ValueError):
            QuantumExecutionSandbox(0)


class PipelineTests(unittest.TestCase):
    def test_two_cli_worker_processes(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = prepare(Path(tmp)/"workers", "plan", 100_000_000)
            workers = [subprocess.Popen([sys.executable, str(ROOT/"run.py"), "--out", str(out), "--worker"],
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE) for _ in range(2)]
            for worker in workers:
                stdout, stderr = worker.communicate(timeout=30)
                self.assertEqual(worker.returncode, 0, stderr.decode(errors="replace"))
            ledger = Ledger(out/"ledger.sqlite")
            try:
                tasks = ledger.snapshot()["tasks"]
                self.assertEqual(len(tasks), 7)
                self.assertTrue(all(t["fence"] == 1 for t in tasks))
                self.assertTrue(all(t["status"] in {"SUCCEEDED", "BLOCKED"} for t in tasks))
            finally:
                ledger.close()

    def test_plan_resume_and_no_reexecution(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = prepare(Path(tmp)/"run", "plan", 100_000_000)
            first = work(out, "unitworker")
            second = work(out, "again")
            self.assertEqual(first, second)
            self.assertEqual(first["actual_provider_tokens"], 0)
            self.assertEqual(len(first["tasks"]), 7)
            self.assertFalse(any(t["status"] == "FAILED" for t in first["tasks"]))

    def test_source_tamper_blocks_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = prepare(Path(tmp)/"run", "plan", 100)
            p = out/"manifest.json"
            spec = json.loads(p.read_text(encoding="utf-8"))
            spec["source_hashes"] = {}
            p.write_text(json.dumps(spec), encoding="utf-8")
            with self.assertRaises(ValueError):
                validate_resume(out)

    def test_no_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(FileExistsError):
                prepare(Path(tmp), "plan", 100)

    def test_actual_macrocycle_minipilot(self):
        from rdkit import Chem
        gen = load_module("01_macrocycle_generator")
        with tempfile.TemporaryDirectory() as tmp:
            result = gen.generate(Path(tmp)/"conformers", smiles=[gen.benchmark_smiles()[0]], conformers=1)
            record = result["records"][0]
            self.assertGreaterEqual(record["ring_size"],12)
            self.assertFalse(result["absolute_ring_strain_computed"])
            self.assertEqual(result["proposed_leads"],[])
            mol = next(iter(Chem.SDMolSupplier(str(Path(tmp)/"conformers"/record["sdf"]))))
            self.assertIsNotNone(mol)
            self.assertIn("@", Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True))
