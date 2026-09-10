"""AST-only geometry fixtures: no entrypoint, torch, training, MD or probe."""
import ast
import copy
import math
from pathlib import Path
import builtins
from unittest.mock import patch
from types import SimpleNamespace
import unittest

import numpy as np


SOURCE = Path(__file__).with_name("run_phase19_active_inference_denovo_enzyme.py")
TREE = ast.parse(SOURCE.read_text(encoding="utf-8"))


def fixture():
    names = {"place_atom", "step_forward", "build_canonical_helix",
             "close_peptide_loop", "require_constructed_backbone",
             "heavy_clash_diagnostic", "realize_backbone", "backbone_torsions",
             "dihedral", "ramachandran_allowed"}
    nodes = [n for n in TREE.body if isinstance(n, ast.FunctionDef) and n.name in names]
    assert {n.name for n in nodes} == names
    ns = dict(np=np, math=math, HELIX_PHI=-57., HELIX_PSI=-47.)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SOURCE), "exec"), ns)
    return ns


class LoopClosureTests(unittest.TestCase):
    def setUp(self):
        self.ns = fixture()
        self.helix = self.ns["build_canonical_helix"](12)

    def solve(self, left, right, **kwargs):
        return self.ns["close_peptide_loop"](left, right, **kwargs)

    def test_exact_solution_and_determinism_leave_fixed_atoms_exact(self):
        before = copy.deepcopy(self.helix)
        obstacle = self.helix[5:]
        a = self.solve(self.helix[0], self.helix[4], obstacles=obstacle,
                       min_res=3, max_res=3)
        b = self.solve(self.helix[0], self.helix[4], obstacles=obstacle,
                       min_res=3, max_res=3)
        self.ns["require_constructed_backbone"]([self.helix[0]] + a + self.helix[4:])
        for x, y in zip(a, b):
            for n in x:
                np.testing.assert_array_equal(x[n], y[n])
        for x, y in zip(before, self.helix):
            for n in x:
                np.testing.assert_array_equal(x[n], y[n])

    def test_noncanonical_endpoint_requires_torsion_optimization(self):
        chain = [self.helix[0]]
        for phi, psi in [(-62., -42.), (-52., -53.), (-64., -39.), (-57., -47.)]:
            chain.append(self.ns["step_forward"](chain, phi, psi))
        loop = self.solve(chain[0], chain[-1], min_res=3, max_res=3)
        actual = [chain[0]] + loop + [chain[-1]]
        self.ns["require_constructed_backbone"](actual)
        projected = self.ns["step_forward"](loop, -57., -47.)
        for n in ("N", "CA"):
            np.testing.assert_allclose(projected[n], chain[-1][n], atol=1e-6, rtol=0)
        self.assertGreater(np.linalg.norm(loop[1]["CA"] - self.helix[2]["CA"]), .01)
        for phi, psi, omega in self.ns["backbone_torsions"](actual)[1:]:
            self.assertLess(abs(abs(omega) - 180.), 1e-4)

    def test_both_right_endpoint_atoms_are_constraints(self):
        for name in ("N", "CA"):
            right = copy.deepcopy(self.helix[4])
            right[name] += np.array([25., 0., 0.])
            with self.assertRaisesRegex(AssertionError, "BACKBONE_REJECT"):
                self.solve(self.helix[0], right, min_res=3, max_res=3,
                           max_evaluations=100)

    def test_future_obstacle_and_exact_overlap_cannot_be_ignored(self):
        # First new N is fixed by the left carbonyl for every torsion solution.
        # An unrelated future rod atom there makes this closure impossible.
        obstacle = {n: self.helix[1]["N"].copy() for n in ("N", "CA", "C", "O")}
        with self.assertRaisesRegex(AssertionError, "BACKBONE_REJECT"):
            self.solve(self.helix[0], self.helix[4], obstacles=[obstacle],
                       min_res=3, max_res=3, max_evaluations=100)

    def test_budget_unreachable_invalid_and_length_reject(self):
        right = {n: p + 100. for n, p in self.helix[4].items()}
        with self.assertRaisesRegex(AssertionError, "evaluations=0/"):
            self.solve(self.helix[0], right)
        with self.assertRaisesRegex(AssertionError, "evaluations=1/1"):
            self.solve(self.helix[0], self.helix[4], max_evaluations=1)
        with self.assertRaisesRegex(AssertionError, "200 residues"):
            self.solve(self.helix[0], self.helix[4], obstacles=[self.helix[5]] * 197)
        with self.assertRaises(ValueError):
            self.solve(self.helix[0], self.helix[4], max_res=9)
        right = copy.deepcopy(self.helix[4])
        right["N"][0] = np.nan
        with self.assertRaisesRegex(AssertionError, "invalid loop endpoint"):
            self.solve(self.helix[0], right)

    def test_realization_integrates_all_fixed_rods_and_keeps_pins(self):
        ns = self.ns
        chain = ns["build_canonical_helix"](62)
        rods = [chain[i:i+6] for i in range(0, 64, 8)]
        saved = copy.deepcopy(rods)
        tz = SimpleNamespace(ca_glu=rods[1][1]["CA"].copy(),
                             ca_trp=rods[4][1]["CA"].copy(),
                             n_don1=rods[6][0]["N"].copy())
        def check(ok, message):
            if not ok:
                raise AssertionError(message)
        solver = ns["close_peptide_loop"]
        calls = []
        def checked(left, right, obstacles, **kwargs):
            k = len(calls) + 1
            # Every fixed residue other than the two endpoints, plus prior loops.
            self.assertEqual(len(obstacles), 46 + 2 * (k - 1))
            for rod in rods[k+1:]:
                for residue in rod:
                    self.assertTrue(any(r is residue for r in obstacles))
            calls.append(kwargs)
            return solver(left, right, obstacles, **kwargs)
        ns.update(N_HELIX=8, N_HELIX_RES=48, HELIX_STARTS=list(range(0, 48, 6)),
                  HELIX_LENS=[6]*8, MOTIF_POS=dict(GLU=7, TRP=25, ASN=36, SER=42),
                  tassert=check, place_free_bundle=lambda *a, **k: (rods, 0.),
                  _best_chain_order=lambda r: (list(range(8)), 5.),
                  kabsch=lambda p, x: (np.eye(3), np.zeros(3)), log=lambda m: None,
                  close_peptide_loop=checked)
        sample = dict(x=np.array([r["CA"] for rod in rods for r in rod]),
                      R=np.tile(np.eye(3), (48, 1, 1)), mask=np.zeros(48))
        atoms, slots, report = ns["realize_backbone"](sample, tz, np.random.default_rng(7))
        self.assertEqual(report["loop_sizes"], [2]*7)
        self.assertEqual(report["construction"]["backbone_clashes"], 0)
        self.assertLessEqual(len(atoms), 200)
        for old, new in zip(saved, rods):
            for a, b in zip(old, new):
                for n in a:
                    np.testing.assert_array_equal(a[n], b[n])
        np.testing.assert_array_equal(atoms[slots["GLU"]]["CA"], tz.ca_glu)
        np.testing.assert_array_equal(atoms[slots["TRP"]]["CA"], tz.ca_trp)
        np.testing.assert_array_equal(atoms[slots["ASN"]]["N"], tz.n_don1)

    def test_no_torch_import_or_spline_fallback(self):
        # Other suites legitimately load torch. Guard this operation's imports,
        # not global interpreter history (and do not dump sys.modules on failure).
        original_import = builtins.__import__
        def guarded_import(name, *args, **kwargs):
            if name == 'torch' or name.startswith('torch.'):
                raise AssertionError('Geometry-only fixture attempted a torch import')
            return original_import(name, *args, **kwargs)
        with patch.object(builtins, '__import__', guarded_import):
            ns = fixture()
            chain = ns['build_canonical_helix'](5)
            loop = ns['close_peptide_loop'](chain[0], chain[4], min_res=3, max_res=3)
            self.assertEqual(len(loop), 3)
        self.assertFalse(any(isinstance(n, ast.FunctionDef) and n.name == "build_loop_spline"
                             for n in ast.walk(TREE)))

    def test_actual_assembly_reserves_remaining_loops_at_200_boundary(self):
        # Execute only the actual assembly block; no large geometry computation.
        fn = next(n for n in TREE.body if isinstance(n, ast.FunctionDef)
                  and n.name == "realize_backbone")
        def assigns(node, name):
            return isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == name for t in node.targets)
        start = next(i for i, n in enumerate(fn.body) if assigns(n, "atoms"))
        end = next(i for i, n in enumerate(fn.body) if assigns(n, "slot_idx"))
        wrapper = ast.parse("def assemble():\n    pass\n").body[0]
        wrapper.body = copy.deepcopy(fn.body[start:end]) + ast.parse(
            "return atoms, loop_sizes").body
        limits = []
        def fake_solver(left, right, obstacles, max_res):
            limits.append(max_res)
            return [{} for _ in range(max_res)]
        ns = dict(rods=[[{} for _ in range(23)] for _ in range(8)],
                  order=list(range(8)), N_HELIX=8, close_peptide_loop=fake_solver)
        exec(compile(ast.fix_missing_locations(ast.Module(body=[wrapper], type_ignores=[])),
                     str(SOURCE), "exec"), ns)
        atoms, sizes = ns["assemble"]()
        self.assertEqual(len(atoms), 200)
        self.assertEqual(limits, [4] + [2]*6)
        ns["rods"] = [[{} for _ in range(24)] for _ in range(8)]
        limits.clear()
        with self.assertRaisesRegex(AssertionError, "BACKBONE_REJECT.*200 residues"):
            ns["assemble"]()
        self.assertEqual(limits, [])


if __name__ == "__main__":
    unittest.main()
