"""Single-thread, bounded actual-source tests; never import the Phase19 entrypoint."""
import ast
import copy
import math
from pathlib import Path
import unittest

import numpy as np


SOURCE = Path(__file__).with_name("run_phase19_active_inference_denovo_enzyme.py")
TREE = ast.parse(SOURCE.read_text(encoding="utf-8"))


def production_functions(names, namespace=None):
    ns = dict(np=np, math=math, HELIX_PHI=-57., HELIX_PSI=-47.)
    if namespace:
        ns.update(namespace)
    nodes = [node for node in TREE.body if isinstance(node, ast.FunctionDef)
             and node.name in names]
    assert {n.name for n in nodes} == set(names)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SOURCE), "exec"), ns)
    return ns


class ClashTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ns = production_functions({"heavy_clash_diagnostic", "place_atom",
                                       "build_canonical_helix", "require_constructed_backbone"})

    def audit(self, atoms, seq):
        return self.ns["heavy_clash_diagnostic"](atoms, seq)

    def test_peptide_topology_matches_independent_shortest_path_oracle(self):
        # All 8 atoms coincide: every nonexcluded pair MUST remain a clash.
        atoms = [{n: np.zeros(3) for n in ("N", "CA", "C", "O")} for _ in range(2)]
        adjacency = np.eye(8, dtype=bool)
        for a, b in [(0, 1), (1, 2), (2, 3), (2, 4), (4, 5), (5, 6), (6, 7)]:
            adjacency[a, b] = adjacency[b, a] = True
        excluded = (adjacency.astype(int) @ adjacency.astype(int)) > 0
        expected = int(np.triu(~excluded, 1).sum())
        result = self.audit(atoms, ["GLY"] * 2)
        self.assertEqual(result["n_nonbonded_clashes"], expected)
        self.assertEqual(result["n_previously_hidden_local_clashes"], expected)
        self.assertEqual(result["n_legacy_nonlocal_clashes"], 0)
        self.assertEqual(result["pair_classes"]["backbone_backbone"], expected)
        self.assertEqual(result["n_close_bonded_or_1_3_pairs_excluded"], 28 - expected)

    def test_ring_closure_missing_intermediates_and_one_four_retained(self):
        # CD1-NE1-CE2 is 1-3 even when NE1 isn't present in coordinates.
        result = self.audit([{n: np.zeros(3) for n in ("CD1", "CE2", "CH2")}], ["TRP"])
        self.assertEqual(result["n_nonbonded_clashes"], 1)  # CD1--CH2 is 1-5
        self.assertEqual(result["n_close_bonded_or_1_3_pairs_excluded"], 2)
        # CA-CB-CG-CD: 1-4 remains counted, with threshold unchanged.
        self.assertEqual(self.audit([dict(CA=np.zeros(3), CD=np.zeros(3))],
                                    ["GLU"])["n_nonbonded_clashes"], 1)
        self.assertEqual(self.audit([dict(CA=np.zeros(3), CD=np.array([2.65, 0., 0.]))],
                                    ["GLU"])["n_nonbonded_clashes"], 0)

    def test_actual_source_helices_keep_real_backbone_collisions(self):
        rod = self.ns["build_canonical_helix"](6)
        atoms = copy.deepcopy(rod) + [{n: p + .1 for n, p in r.items()} for r in rod]
        before = copy.deepcopy(atoms)
        result = self.audit(atoms, ["GLY"] * len(atoms))
        pts = np.array([p for r in atoms for p in r.values()])
        owners = np.repeat(np.arange(len(atoms)), 4)
        distances = np.linalg.norm(pts[:, None] - pts[None, :], axis=-1)
        legacy = int(np.triu((distances < 2.65) &
                            (abs(owners[:, None] - owners[None, :]) >= 2), 1).sum())
        self.assertEqual(result["n_legacy_nonlocal_clashes"], legacy)
        self.assertGreaterEqual(legacy, 24)
        self.assertEqual(result["pair_classes"]["backbone_backbone"],
                         result["n_nonbonded_clashes"])
        for old, new in zip(before, atoms):
            for name in old:
                np.testing.assert_array_equal(old[name], new[name])
        print("actual_source_overlapping_6_residue_rods:", result, flush=True)

    def test_hydrogens_invalid_input_and_bounded_examples(self):
        self.assertEqual(self.audit([dict(H=np.zeros(3), CA=np.zeros(3),
                                         **{"1HB": np.zeros(3)})], ["ALA"])
                         ["n_nonbonded_clashes"], 0)
        with self.assertRaises(ValueError):
            self.audit([dict(CA=np.zeros(3))], [])
        with self.assertRaises(ValueError):
            self.audit([dict(CA=np.zeros(3))], ["UNK"])
        with self.assertRaises(ValueError):
            self.audit([dict(CA=np.full(3, np.nan))], ["ALA"])
        result = self.audit([dict(CA=np.zeros(3)) for _ in range(30)], ["ALA"] * 30)
        self.assertEqual(result["n_nonbonded_clashes"], 435)
        self.assertEqual(len(result["worst_pairs"]), 12)

    def test_static_audit_preserves_historical_count_and_adds_diagnostic(self):
        ns = production_functions({"static_fold_audit", "heavy_clash_diagnostic",
                                   "backbone_torsions", "dihedral", "ramachandran_allowed",
                                   "_aa_class"},
                                  dict(HYDRO=set(), POLAR=set(), CHARGE_P1=set(), CHARGE_M1=set()))
        atoms = self.ns["build_canonical_helix"](6)
        result = ns["static_fold_audit"](atoms, ["GLY"] * 6)
        self.assertEqual(result["n_clashes"], 0)
        self.assertEqual(result["clash_diagnostic"]["n_nonbonded_clashes"], 0)
        self.assertEqual(result["n_heavy"], 24)

    def test_actual_loop_spline_ignores_next_n_and_can_be_collinear(self):
        # Diagnose without changing closure, targets, sampling or running realization.
        node = next(n for n in ast.walk(TREE) if isinstance(n, ast.FunctionDef)
                    and n.name == "build_loop_spline")
        ns = dict(np=np, math=math, cen=np.zeros(3))
        exec(compile(ast.Module(body=[node], type_ignores=[]), str(SOURCE), "exec"), ns)
        a, b = np.array([10., 0., 0.]), np.array([27.2, 0., 0.])
        loop = ns["build_loop_spline"](a, b - 1., b)
        changed_n = ns["build_loop_spline"](a, b + 20., b)
        angles = []
        for r, other in zip(loop, changed_n):
            for name in r:
                np.testing.assert_array_equal(r[name], other[name])
            n, c = r["N"] - r["CA"], r["C"] - r["CA"]
            angles.append(math.degrees(math.acos(np.clip(n @ c /
                          (np.linalg.norm(n) * np.linalg.norm(c)), -1., 1.))))
        np.testing.assert_allclose(angles, 180.)
        print("actual_source_collinear_spline_N_CA_C_angles_deg:", angles, flush=True)
        with self.assertRaisesRegex(AssertionError, r"BACKBONE_REJECT.*N-CA-C=180"):
            self.ns["require_constructed_backbone"](loop)

    def test_constructor_accepts_canonical_geometry_without_movement(self):
        atoms = self.ns["build_canonical_helix"](12)
        before = copy.deepcopy(atoms)
        self.assertEqual(self.ns["require_constructed_backbone"](atoms),
                         dict(constructor_geometry_valid=True, backbone_clashes=0))
        for old, new in zip(before, atoms):
            for name in old:
                np.testing.assert_array_equal(old[name], new[name])

    def test_constructor_rejects_open_junctions_and_backbone_intersections(self):
        rod = self.ns["build_canonical_helix"](6)
        separated = copy.deepcopy(rod)
        for residue in separated[3:]:
            for name in residue:
                residue[name] += np.array([30., 0., 0.])
        with self.assertRaisesRegex(AssertionError, r"junction 2->3 C-N="):
            self.ns["require_constructed_backbone"](separated)
        overlap = copy.deepcopy(rod) + [{n: p + .1 for n, p in r.items()} for r in rod]
        with self.assertRaisesRegex(AssertionError, r"backbone_clashes=124"):
            self.ns["require_constructed_backbone"](overlap)
        malformed = copy.deepcopy(rod)
        malformed[0]["N"][:] = np.nan
        with self.assertRaisesRegex(AssertionError, "invalid coordinates"):
            self.ns["require_constructed_backbone"](malformed)

    def test_realization_rejects_actual_spline_before_downstream_gates(self):
        # Supply small deterministic rods at the assembly boundary, then execute
        # the actual realization body and its spline; no flow, search or packing.
        from types import SimpleNamespace
        ns = production_functions({"realize_backbone", "require_constructed_backbone",
                                   "heavy_clash_diagnostic", "geometry_probe",
                                   "residue_frames_from_atoms"})
        chain = self.ns["build_canonical_helix"](48)
        rods = [chain[i:i + 6] for i in range(0, 48, 6)]
        tz = SimpleNamespace(ca_glu=rods[1][1]["CA"].copy(),
                             ca_trp=rods[4][1]["CA"].copy(),
                             n_don1=rods[6][0]["N"].copy())
        saved = copy.deepcopy(rods)
        def check(condition, message):
            if not condition:
                raise AssertionError(message)
        def forbidden(*args):
            raise RuntimeError("downstream torsion gate must not receive malformed backbone")
        ns.update(N_HELIX=8, N_HELIX_RES=48, HELIX_STARTS=list(range(0, 48, 6)),
                  HELIX_LENS=[6] * 8, MOTIF_POS=dict(GLU=7, TRP=25, ASN=36, SER=42),
                  tassert=check, place_free_bundle=lambda *a, **k: (rods, 0.),
                  _best_chain_order=lambda r: (list(range(8)), 1.335),
                  kabsch=lambda p, x: (np.eye(3), np.zeros(3)), backbone_torsions=forbidden)
        sample = dict(x=np.array([r["CA"] for r in chain]),
                      R=np.tile(np.eye(3), (48, 1, 1)), mask=np.zeros(48))
        with self.assertRaisesRegex(AssertionError, "BACKBONE_REJECT: implement peptide closure") as caught:
            ns["realize_backbone"](sample, tz, np.random.default_rng(7))
        print("actual_realization_fixture_rejection:", caught.exception, flush=True)
        # The actual probe caller also stops, without producing a success JSON.
        # All expensive upstream generation remains replaced by the tiny fixture.
        ns.update(CONFIG=dict(SEED=7), S_PRIOR_MEAN=np.zeros(6), Theozyme=lambda s: tz,
                  design_sequence=forbidden, pack_sidechains=forbidden, static_fold_audit=forbidden)
        with self.assertRaisesRegex(AssertionError, "BACKBONE_REJECT"):
            ns["geometry_probe"]()
        for old_rod, new_rod in zip(saved, rods):
            for old, new in zip(old_rod, new_rod):
                for name in old:
                    np.testing.assert_array_equal(old[name], new[name])
        np.testing.assert_array_equal(tz.ca_glu, rods[1][1]["CA"])
        np.testing.assert_array_equal(tz.ca_trp, rods[4][1]["CA"])

    def test_actual_scaffold_packing_registers_sidechains_and_penalizes_coincidence(self):
        # Compile the unchanged production scaffold block inside a small wrapper.
        # Skip only catalytic inverse fitting, which needs a full probe allocation.
        from test_phase19_orientation import fixture
        ns, _, _ = fixture()
        fn = next(n for n in TREE.body if isinstance(n, ast.FunctionDef)
                  and n.name == "pack_sidechains")
        def assignment(node, name):
            return isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == name for t in node.targets)
        start = next(i for i, n in enumerate(fn.body) if assignment(n, "all_atoms"))
        end = next(i for i, n in enumerate(fn.body) if assignment(n, "signs"))
        wrapper = ast.parse("def scaffold(atoms, seq, slot_idx, outward):\n    pass\n").body[0]
        wrapper.body = copy.deepcopy(fn.body[start:end]) + ast.parse(
            "return ref, ref_res, clash_score").body
        for node in TREE.body:
            if assignment(node, "ROTAMER_CHI"):
                exec(compile(ast.Module(body=[node], type_ignores=[]), str(SOURCE), "exec"), ns)
        exec(compile(ast.fix_missing_locations(ast.Module(body=[wrapper], type_ignores=[])),
                     str(SOURCE), "exec"), ns)
        atoms = ns["build_canonical_helix"](3)
        saved = copy.deepcopy(atoms)
        ref, owners, score = ns["scaffold"](atoms, ["ALA", "SER", "SER"],
                                            {"GLU": -1, "TRP": -2}, [1.] * 3)
        self.assertEqual(len(ref), sum(map(len, atoms)))
        for i, residue in enumerate(atoms):
            np.testing.assert_array_equal(ref[owners == i], np.array(list(residue.values())))
            for name in saved[i]:
                np.testing.assert_array_equal(residue[name], saved[i][name])
        candidate = atoms[0]["CB"]
        distances = np.linalg.norm(ref - candidate, axis=1)
        expected = float(np.clip(2.9 - distances[owners != 2], 0., None).sum())
        self.assertAlmostEqual(score([candidate], 2), expected)
        self.assertGreaterEqual(score([candidate], 2), 2.9)
        # Demonstrate why refreshing matters: the new SER OG sees ALA CB.
        self.assertTrue(any(np.array_equal(p, candidate) for p in ref[owners == 0]))


if __name__ == "__main__":
    unittest.main()
