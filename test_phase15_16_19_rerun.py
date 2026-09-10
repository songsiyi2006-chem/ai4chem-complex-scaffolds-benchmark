"""Cheap real-source regression checks; these are not production reruns."""
import ast
import json
import math
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parent


def tree(phase):
    return ast.parse(next(ROOT.glob(f"run_phase{phase}_*.py")).read_text(encoding="utf-8"))


def functions(phase, names, **extra):
    ns = dict(np=np, math=math, **extra)
    nodes = [n for n in ast.walk(tree(phase))
             if isinstance(n, ast.FunctionDef) and n.name in names]
    assert set(names) == {n.name for n in nodes}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), "<actual-source>", "exec"), ns)
    return ns


class RerunTests(unittest.TestCase):
    def test_phase19_branch_angles_planarity_and_bond_lengths(self):
        f = functions(19, ["_trigonal_third"])["_trigonal_third"]
        c = np.zeros(3)
        a = np.array([1., 0., 0.])
        b = np.array([-.5, np.sqrt(3)/2, 0.])
        p = f(c, a, b, 1.25, 120.)
        self.assertAlmostEqual(np.linalg.norm(p), 1.25)
        self.assertAlmostEqual(p[2], 0.)
        self.assertAlmostEqual(p @ a / 1.25, -.5)
        self.assertAlmostEqual(p @ b / 1.25, -.5)
        b = np.array([-1/3, np.sqrt(8)/3, 0.])
        p = f(c, a, b, 1.53, math.degrees(math.acos(-1/3)))
        self.assertAlmostEqual(np.linalg.norm(p), 1.53)
        self.assertAlmostEqual(p @ a / 1.53, -1/3)
        self.assertAlmostEqual(p @ b / 1.53, -1/3)
        with self.assertRaises(ValueError):
            f(c, a, c, 1.25, 120.)

    def test_phase19_trp_named_ring_bonds_close_without_duplicate_atoms(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        m = Chem.AddHs(Chem.MolFromSmiles("c1ccc2[nH]ccc2c1"))
        self.assertEqual(AllChem.EmbedMolecule(m, randomSeed=485), 0)
        AllChem.MMFFOptimizeMolecule(m)
        xyz = np.asarray(m.GetConformer().GetPositions())
        symbols = [a.GetSymbol() for a in m.GetAtoms()]
        geom = next(n.value for n in tree(19).body if isinstance(n, ast.Assign)
                    and any(isinstance(t, ast.Name) and t.id == "SC_GEOM" for t in n.targets))
        ns = functions(19, ["build_sidechain", "place_atom", "place_cb", "_trigonal_third",
                            "_indole_attach", "_pca_frame", "_rot_about_axis"],
                       SC_GEOM=ast.literal_eval(geom), _l_cb_reference_sign=lambda: 1.,
                       _cached_geom=lambda *a: (xyz, symbols))
        r = dict(N=np.array([0.,0.,0.]), CA=np.array([1.458,0.,0.]), C=np.array([2.,1.42,0.]))
        out = ns["build_sidechain"]("TRP", r, [60., 90.], 1.)
        bonds = [("CG","CD1"),("CD1","NE1"),("NE1","CE2"),("CE2","CD2"),
                 ("CD2","CG"),("CD2","CE3"),("CE3","CZ3"),("CZ3","CH2"),
                 ("CH2","CZ2"),("CZ2","CE2")]
        for a,b in bonds:
            self.assertTrue(1.25 < np.linalg.norm(out[a]-out[b]) < 1.55, (a,b))
        ring = np.array([out[n] for n in out if n != "CB"])
        self.assertLess(np.linalg.svd(ring-ring.mean(0), compute_uv=False)[-1], .02)
        self.assertEqual(len(np.unique(np.round(ring, 5), axis=0)), 9)
        asp = ns["build_sidechain"]("ASP", r, [60., 180.], 1.)
        self.assertAlmostEqual(np.linalg.norm(asp["OD2"]-asp["CG"]), 1.25)
        plane = np.cross(asp["OD1"]-asp["CG"], asp["CB"]-asp["CG"])
        self.assertAlmostEqual((asp["OD2"]-asp["CG"]) @ plane, 0.)

    def test_phase15_sampling_is_fixed_and_accounted(self):
        fn = next(n for n in tree(15).body if isinstance(n, ast.FunctionDef)
                  and n.name == "run_openmm_allostery")
        source = ast.unparse(fn)
        self.assertNotIn("np.clip", source)
        self.assertNotIn("planned_s", source)
        self.assertIn("executed_steps != requested_steps", source)
        calls = [n for n in ast.walk(fn) if isinstance(n, ast.Call)
                 and isinstance(n.func, ast.Attribute) and n.func.attr == "step"]
        self.assertEqual(len(calls), 1)  # all integration uses the accounting wrapper
        config = next(n.value for n in tree(15).body if isinstance(n, ast.Assign)
                      and any(isinstance(t, ast.Name) and t.id == "CONFIG" for t in n.targets))
        values = {k.arg: ast.literal_eval(k.value) for k in config.keywords}
        self.assertEqual(values["MD_PRODUCTION_PS"], 300.)
        self.assertEqual(values["UMB_SAMPLE_PS"], 150.)
        self.assertEqual((70 + 300 + 8 * (12 + 150)) / .002, 833000)

    def test_phase15_allostery_summary_uses_current_schema(self):
        fn = next(n for n in tree(15).body if isinstance(n, ast.FunctionDef)
                  and n.name == "run_openmm_allostery")
        keys = {n.slice.value for n in ast.walk(fn) if isinstance(n, ast.Subscript)
                and isinstance(n.value, ast.Name) and n.value.id == "allo"
                and isinstance(n.slice, ast.Constant)}
        self.assertEqual(keys, {"dG_latch_shift_kcal"})

    def test_phase15_terminal_carboxylate_oxygens_are_separate(self):
        ns = functions(15, ["place_atom", "build_helix_backbone"])
        terminal = ns["build_helix_backbone"](3)[-1]
        self.assertAlmostEqual(np.linalg.norm(terminal["OXT"] - terminal["C"]), 1.25)
        self.assertGreater(np.linalg.norm(terminal["OXT"] - terminal["O"]), 2.)

    def test_fresh_stage_assembly_keeps_input_hashes_and_does_not_overwrite(self):
        ns = functions(15, ["integrate_spin_yields", "assemble_fresh_stages"],
                       Path=Path, json=json, log=lambda _: None,
                       __file__="run_phase15_quantum_biology_spin_allostery.py")
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            spin_path, allo_path = folder / "fresh_spin.json", folder / "allostery_results.json"
            t = np.linspace(0., 1., 101)
            curve = dict(t_s=t.tolist(), PS=np.exp(-t).tolist(), PT=np.zeros_like(t).tolist(),
                         survival=np.exp(-t).tolist(), phi_S=.6, phi_T=0.)
            record = dict(config=dict(K_S=1., K_T=1.), spin=dict(full_dynamics={
                          name:curve for name in ("theta0","theta30","theta60","theta90","Bzero")}))
            raw = json.dumps(record)
            spin_path.write_text(raw, encoding="utf-8")
            allo_path.write_text(json.dumps(dict(states={"test":{}}, wham={"test":{}})), encoding="utf-8")
            combined = ns["assemble_fresh_stages"](spin_path, allo_path)
            self.assertFalse(combined["stage_assembly"]["single_default_run"])
            self.assertFalse(combined["stage_assembly"]["spin_dynamics_recomputed"])
            self.assertEqual(len(combined["stage_assembly"]["inputs"]), 2)
            self.assertEqual(spin_path.read_text(encoding="utf-8"), raw)
            with self.assertRaises(FileExistsError):
                ns["assemble_fresh_stages"](spin_path, allo_path)

    def test_phase15_sidechains_bond_to_correct_centers_and_keep_chirality(self):
        spec = next(n.value for n in tree(15).body if isinstance(n, ast.Assign)
                    and any(isinstance(t, ast.Name) and t.id == "SIDECHAIN" for t in n.targets))
        ns = functions(15, ["sidechain_geometry", "place_atom", "build_helix_backbone"],
                       SIDECHAIN=ast.literal_eval(spec))
        frame = ns["build_helix_backbone"](3)[1]
        for residue in ("ALA", "ASN", "ASP", "GLN", "LYS"):
            atoms = ns["sidechain_geometry"](frame, residue, -1.)
            self.assertAlmostEqual(np.linalg.norm(atoms["CB"] - frame["CA"]), 1.53)
            self.assertLess(np.linalg.det(np.column_stack([frame["N"]-frame["CA"],
                                  frame["C"]-frame["CA"], atoms["CB"]-frame["CA"]])), 0.)
            if residue in ("ASN", "ASP", "GLN"):
                center, oxygen, terminal, length = {
                    "ASN": ("CG", "OD1", "ND2", 1.33),
                    "ASP": ("CG", "OD1", "OD2", 1.25),
                    "GLN": ("CD", "OE1", "NE2", 1.33)}[residue]
                self.assertAlmostEqual(np.linalg.norm(atoms[terminal]-atoms[center]), length)
                self.assertGreater(np.linalg.norm(atoms[terminal]-atoms[oxygen]), 2.)

    def test_phase15_bias_matches_unweighted_sampled_centroids(self):
        seen = []
        class Force:
            def __init__(self, *args): pass
            def addGlobalParameter(self, *args): pass
            def addGroup(self, group, weights):
                seen.append((group, weights))
                return len(seen)-1
            def addBond(self, *args): pass
        ns = functions(15, ["add_bias"], CustomCentroidBondForce=Force,
                       unit=SimpleNamespace(kilojoule_per_mole=1., nanometer=1.))
        ns["add_bias"](SimpleNamespace(addForce=lambda _: None), [[1,2,3], [4,5]])
        self.assertEqual(seen, [([1,2,3], [1.,1.,1.]), ([4,5], [1.,1.])])

    def test_so3_log_float32_near_pi_roundtrip_and_random_batches(self):
        from scipy.spatial.transform import Rotation
        ns = functions(19, ["_so3_log_batch", "so3_log"])
        axes = np.array([[1.,0.,0.], [1.,-2.,3.], [-1.,1.,0.]])
        axes /= np.linalg.norm(axes, axis=1)[:, None]
        vectors = np.concatenate([axes * angle for angle in (0., 1e-10, math.pi-1e-4,
                                                                            math.pi-1e-8, math.pi)])
        matrices = Rotation.from_rotvec(vectors).as_matrix().astype(np.float32)
        identity = np.tile(np.eye(3), (len(matrices),1,1))
        logs = ns["_so3_log_batch"](identity, matrices)
        self.assertLessEqual(np.linalg.norm(logs, axis=1).max(), math.pi + 1e-12)
        np.testing.assert_allclose(Rotation.from_rotvec(logs).as_matrix(), matrices, atol=2e-7)
        rng = np.random.default_rng(875)
        a = Rotation.random(5000, random_state=rng).as_matrix().astype(np.float32)
        b = Rotation.random(5000, random_state=rng).as_matrix().astype(np.float32)
        logs = ns["_so3_log_batch"](a, b)
        np.testing.assert_allclose(Rotation.from_rotvec(logs).as_matrix(),
                                   a.transpose(0,2,1).astype(float) @ b.astype(float), atol=3e-7)
        self.assertLessEqual(np.linalg.norm(logs, axis=1).max(), math.pi + 1e-12)
        np.testing.assert_allclose(Rotation.from_rotvec(ns["so3_log"](matrices[-1])).as_matrix(),
                                   matrices[-1], atol=2e-7)

    def test_so3_log_rejects_bad_rotations_instead_of_hiding_them(self):
        ns = functions(19, ["_so3_log_batch"])
        for bad in (np.diag([1.,1.,-1.]), np.eye(3)*2, np.full((3,3),np.nan)):
            with self.assertRaises(ValueError):
                ns["_so3_log_batch"](np.eye(3), bad)

    def test_simpson_yields_against_equal_rate_analytic_trace(self):
        ns = functions(15, ["integrate_spin_yields", "lindblad_grid"], CONFIG=dict(T_MAX_US=4.))
        t, _ = ns["lindblad_grid"](None)
        for rate in (1e6, 1e7):
            survival = np.exp(-rate * t)
            result = ns["integrate_spin_yields"](t, survival * .7, survival * .3, survival, rate, rate)
            self.assertLess(abs(result["closure_error"]), 2e-4)
            self.assertLess(abs(result["closure_error"]), abs(result["trapezoid_closure_error"]) / 10.)
            self.assertAlmostEqual(result["equal_rate_analytic_total"], 1 - math.exp(-rate*t[-1]))
            self.assertEqual(result["equal_rate_trace_max_error"], 0.)
            self.assertFalse(result["individual_yields_grid_converged"])

    def test_simpson_does_not_renormalize_bad_data(self):
        ns = functions(15, ["integrate_spin_yields"])
        t = np.linspace(0, 1, 101)
        result = ns["integrate_spin_yields"](t, np.exp(-t)*2, np.zeros_like(t), np.exp(-t), 1., 1.)
        self.assertGreater(result["closure_error"], .6)
        self.assertGreater(result["trace_population_mismatch_max"], .9)

    def test_strict_candidate_acceptance_rejects_missing_failed_and_bad_measurements(self):
        ns = functions(19, ["candidate_gate_failures"], CONFIG=dict(RMSF_GATE=.8, BARRIER_GATE=12.0, MD_PS_CAND=40.))
        gate = ns["candidate_gate_failures"]
        valid = dict(feas=True, const_rmsd=.1, md=dict(n_frames=10, simulated_ns=.04,
                     rmsf_constellation_A=.7, ca_rmsd_to_design_A=1.), qmmm=dict(barrier_kcal=10.))
        self.assertEqual(gate(valid), [])
        for rmsf in (.80001, 2., 25., float("nan")):
            self.assertTrue(gate(valid | dict(md=valid["md"] | dict(rmsf_constellation_A=rmsf))))
        self.assertTrue(gate(valid | dict(md=valid["md"] | dict(error="simulation failed"))))
        self.assertTrue(gate(valid | dict(md=valid["md"] | dict(n_frames=1))))
        for duration in (0.001, 0.0399, float("nan")):
            self.assertIn("md_below_full_declared_production_duration",
                          gate(valid | dict(md=valid["md"] | dict(simulated_ns=duration))))
        self.assertTrue(gate(valid | dict(qmmm={})))
        self.assertTrue(gate(valid | dict(qmmm=dict(barrier_kcal=12.6))))
        self.assertTrue(gate(valid | dict(qmmm=dict(barrier_kcal=12.4))))
        self.assertTrue(gate(valid | dict(const_rmsd=.31)))
        self.assertEqual(gate(valid | dict(qmmm={}), require_qm=False), [])
        config = next(n for n in tree(19).body if isinstance(n, ast.Assign)
                      and any(isinstance(t, ast.Name) and t.id == "CONFIG" for t in n.targets))
        self.assertEqual(next(k.value.value for k in config.value.keywords if k.arg == "RMSF_GATE"), .8)
        self.assertEqual(next(k.value.value for k in config.value.keywords if k.arg == "BARRIER_GATE"), 12.)

    def test_saved_acceptance_audit_withdraws_12_point_4_without_mutating_source(self):
        ns = functions(19, ["candidate_gate_failures", "audit_saved_acceptance"],
                       CONFIG=dict(RMSF_GATE=.8, BARRIER_GATE=12., MD_PS_CAND=40.), Path=Path, json=json, log=lambda _: None)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "phase19_results.json"
            result = dict(config=dict(BARRIER_GATE=12.5), generations=[dict(generation=1,
                          candidates=[dict(candidate=0, feasible=True, constellation_rmsd_A=.1,
                          md=dict(n_frames=10, simulated_ns=.04, rmsf_constellation_A=.7,
                                  ca_rmsd_to_design_A=1.), qmmm=dict(barrier_kcal=12.4))])])
            raw = json.dumps(result)
            path.write_text(raw, encoding="utf-8")
            audit = ns["audit_saved_acceptance"](path)
            self.assertEqual(audit["acceptance_status"], "no_accepted_design")
            self.assertEqual(audit["snapshot_barrier_gate"], 12.5)
            self.assertEqual(audit["audit_barrier_gate"], 12.)
            self.assertEqual(path.read_text(encoding="utf-8"), raw)

    def test_compile_owned_phases(self):
        for phase in (15, 16, 19):
            compile(tree(phase), str(phase), "exec")

    def test_rot_between_nonunit_and_antiparallel_vectors(self):
        ns = functions(19, ["_rot_between"])
        for a, b in (([2.,0.,0.], [-3.,0.,0.]), ([1.,2.,3.], [-1.,-2.,-3.]),
                     ([2.,1.,0.], [0.,0.,4.])):
            a, b = np.array(a), np.array(b)
            rotation = ns["_rot_between"](a, b)
            np.testing.assert_allclose(rotation @ (a / np.linalg.norm(a)),
                                       b / np.linalg.norm(b), atol=1e-12)
            np.testing.assert_allclose(rotation.T @ rotation, np.eye(3), atol=1e-12)
            self.assertAlmostEqual(np.linalg.det(rotation), 1.)

    def test_pyscf_bohr_input_is_explicit(self):
        ns = functions(15, ["ub3lyp_single"],
                       gto=None, dft=None, log=lambda _: None)
        class MF:
            converged = True
            grids = SimpleNamespace()
            def kernel(self): return -1.0
            def spin_square(self): return (0.75, 2.0)
        class Fake:
            def M(self, **kw):
                self.kw = kw
                return kw
            def UKS(self, m): return MF()
        fake = Fake()
        ns.update(gto=fake, dft=fake)
        xyz = [("H", (0., 0., 1.889726))]
        ns["ub3lyp_single"](xyz, 0, 1, "test")
        self.assertEqual(fake.kw["unit"], "Bohr")
        self.assertEqual(fake.kw["atom"], xyz)

    def test_hfcc_cli_dispatch(self):
        main = next(n for n in tree(15).body if isinstance(n, ast.FunctionDef) and n.name == "main")
        branch = next(n for n in main.body if isinstance(n, ast.If) and "hfcc" in ast.unparse(n.test))
        self.assertIn("args.stage == 'hfcc'", ast.unparse(branch.test))
        self.assertEqual(ast.unparse(branch.body[0]), "hfcc_stage()")

    def test_glu_rigid_transform_matches_direct_nerf(self):
        config = next(n.value for n in tree(19).body if isinstance(n, ast.Assign)
                      and any(isinstance(t, ast.Name) and t.id == "SC_GEOM" for t in n.targets))
        ns = functions(19, ["place_atom", "place_cb", "build_sidechain", "_trigonal_third"],
                       SC_GEOM=ast.literal_eval(config), _l_cb_reference_sign=lambda: -1.)
        frame = dict(N=np.array([-1., .8, 0.]), CA=np.zeros(3), C=np.array([1.3, .4, 0.]))
        rng = np.random.default_rng(87)
        for chis in ((-180., -180., -180.), (12., 60., 108.), (156., -60., 36.)):
            rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
            rotation[:, 0] *= np.linalg.det(rotation)
            sc = ns["build_sidechain"]("GLU", frame, chis, 1.)
            self.assertAlmostEqual(np.linalg.norm(sc["CB"] - frame["CA"]), 1.53, places=12)
            direct = ns["build_sidechain"]("GLU", {k: rotation @ v for k, v in frame.items()}, chis, 1.)
            for key in sc:
                np.testing.assert_allclose(rotation @ sc[key], direct[key], atol=2e-14)

    def test_dihedral_ignores_bond_parallel_components(self):
        ns = functions(19, ["dihedral"])
        points = [np.array(p, float) for p in ((0,1,0), (0,0,0), (1,0,0), (1,0,1))]
        angle = ns["dihedral"](*points)
        self.assertAlmostEqual(angle, 90.)
        points[0] += np.array([2.,0.,0.])
        points[3] += np.array([3.,0.,0.])
        self.assertAlmostEqual(ns["dihedral"](*points), angle)

    def test_md_reference_uses_original_pdb_atom_indices(self):
        fun = next(n for n in tree(19).body if isinstance(n, ast.FunctionDef) and n.name == "openmm_fold_check")
        assignment = next(n for n in fun.body if isinstance(n, ast.Assign)
                          and ast.unparse(n.targets[0]) == "ca_design")
        self.assertIn("pdb.topology.atoms()", ast.unparse(assignment))

    def test_rmsf_is_root_mean_square_not_mean_distance(self):
        fun = next(n for n in tree(19).body if isinstance(n, ast.FunctionDef) and n.name == "openmm_fold_check")
        expr = next(n.value for n in fun.body if isinstance(n, ast.Assign)
                    and ast.unparse(n.targets[0]) == "rmsf")
        aligned = np.zeros((3, 1, 3))
        aligned[2, 0, 0] = 3.
        actual = eval(compile(ast.Expression(expr), "<actual-source>", "eval"),
                      dict(np=np, aligned=aligned, mu=aligned.mean(0)))
        np.testing.assert_allclose(actual, [10 * np.sqrt(2.)])

    def test_glu_clearance_uses_placed_stems_and_rejects_all_clashing(self):
        names = ("CB", "CG", "CD", "OE1", "OE2")
        xyz = np.array([[[0.,0.,0.], [0.,0.,0.], [0.,1.,0.], [-4.3,0.,-1.3], [-4.,0.,-1.]]])
        ns = functions(19, ["_select_glu_anchor"],
                       _glu_rotamer_grid=lambda: (names, [(0.,0.,0.)], xyz))
        o = np.zeros(3)
        # Origin coincides with unplaced stem, but placed stems clear it.
        result = ns["_select_glu_anchor"](o, np.array([0.,0.,1.]), np.zeros((1,3)), [np.eye(3)])
        self.assertIsNotNone(result)
        np.testing.assert_allclose(result[1] + result[2]["OE1"], o)
        clash = (result[1] + result[2]["CG"])[None, :]
        self.assertIsNone(ns["_select_glu_anchor"](o, np.array([0.,0.,1.]), clash, [np.eye(3)]))

    def test_acetate_alignment_preserves_distances_and_faces_away(self):
        fun = next(n for n in tree(19).body if isinstance(n, ast.FunctionDef) and n.name == "qmmm_kemp_scan")
        body = next(n.body for n in fun.body if isinstance(n, ast.If))
        start = next(i for i, n in enumerate(body) if isinstance(n, ast.Assign) and ast.unparse(n.targets[0]) == "src_x")
        end = next(i for i, n in enumerate(body) if isinstance(n, ast.Assign) and ast.unparse(n.targets[0]) == "base")
        ace = np.array([[0.,0.,0.], [.2,-.7,.5], [1.,.2,.3], [0.,0.,0.]])
        ns = dict(np=np, ace_pos=ace, attack=np.array([1.,0.,0.]),
                  tz=SimpleNamespace(zhat=np.array([0.,0.,1.])), a=np.array([3.5,0.,0.]))
        exec(compile(ast.Module(body=body[start:end+1], type_ignores=[]), "<actual-source>", "exec"), ns)
        rotation = ns["R_ace"]
        np.testing.assert_allclose(ns["base"][1] - ns["base"][3],
                                   np.array([np.linalg.norm(ace[1]-ace[3]),0.,0.]), atol=1e-12)
        np.testing.assert_allclose(rotation.T @ rotation, np.eye(3), atol=1e-12)
        np.testing.assert_allclose(np.linalg.det(rotation), 1., atol=1e-12)

    def test_transport_rejects_short_unsorted_and_nonfinite_inputs(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "umbrella").mkdir()
            ns = functions(16, ["analyze_transport"], RES=root)
            baseline = dict(z0=np.array([-25., 0., 25.]), force=np.zeros((3,4)),
                            contacts=np.zeros((3,4)), eq_contact_max=np.zeros(3),
                            gamma=np.ones(3), M_cx=1.)
            for changes in (dict(force=np.zeros((3,1))), dict(force=np.full((3,4), np.nan)),
                            dict(z0=np.array([25., 0., -25.])), dict(gamma=np.zeros(3))):
                np.savez(root / "umbrella/windows_receptor.npz", **(baseline | changes))
                with self.assertRaises(ValueError):
                    ns["analyze_transport"]({}, False)


if __name__ == "__main__":
    unittest.main()
