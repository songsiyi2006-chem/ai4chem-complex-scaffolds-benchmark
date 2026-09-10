"""Bounded actual-source geometry fixture; no model, MD, or QM imports."""
import ast
import math
from pathlib import Path
from types import SimpleNamespace
import unittest

import numpy as np


def fixture():
    # One small in-memory indole, matching the prior diagnostic's seed/method.
    from rdkit import Chem
    from rdkit.Chem import AllChem
    molecule = Chem.AddHs(Chem.MolFromSmiles("c1ccc2[nH]ccc2c1"))
    assert AllChem.EmbedMolecule(molecule, randomSeed=485) == 0
    assert AllChem.MMFFOptimizeMolecule(molecule) == 0
    xyz = np.asarray(molecule.GetConformer().GetPositions())
    symbols = [a.GetSymbol() for a in molecule.GetAtoms()]
    source = Path(__file__).with_name("run_phase19_active_inference_denovo_enzyme.py")
    tree = ast.parse(source.read_text(encoding="utf-8"))
    names = {"place_atom", "build_canonical_helix", "canon_helix",
             "residue_frames_from_atoms", "place_helix_anchor", "_rot_between",
             "_rot_about_axis", "rotate_about", "_canon_frame_constants",
             "_canon_cbhat", "_canon_cghat", "place_cb", "_l_cb_reference_sign",
             "build_sidechain", "_indole_attach", "_pca_frame", "_trigonal_third",
             "trp_target_reachability", "trp_optimal_anchor_roll", "_trp_ca_target",
             "_aimed_rod"}
    nodes = [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name in names]
    assert {n.name for n in nodes} == names
    geom = next(n.value for n in tree.body if isinstance(n, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id == "SC_GEOM" for t in n.targets))
    ns = dict(np=np, math=math, SC_GEOM=ast.literal_eval(geom),
              _cached_geom=lambda *args: (xyz, symbols), _CANON_HELIX=None,
              _L_CB_SIGN=None, HELIX_PHI=-57., HELIX_PSI=-47.,
              sub_heavy_s=np.array([[1000., 1000., 1000.]]))
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(source), "exec"), ns)
    ns["_NHAT"], ns["_CHAT"], ns["_KHAT"], ns["_BHAT"] = ns["_canon_frame_constants"]()
    ns["_CBHAT"], ns["_CGHAT"] = ns["_canon_cbhat"](), ns["_canon_cghat"]()
    tz = SimpleNamespace(stack_centroid=np.array([0., 0., 4.]),
                         stack_inplane=ns["rotate_about"](np.array([0., 1., 0.]),
                                                          np.array([0., 0., 1.]), math.radians(25)),
                         zhat=np.array([0., 0., 1.]), yhat=np.array([0., 1., 0.]))
    tz.ca_trp = ns["_trp_ca_target"](tz)
    return ns, tz, tree


class OrientationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ns, cls.tz, cls.tree = fixture()

    def test_actual_aimed_rod_reaches_bound_without_moving_target_or_geometry(self):
        from scipy.optimize import least_squares
        ns, tz = self.ns, self.tz
        target, ca = tz.ne1.copy(), tz.ca_trp.copy()
        direction = target - ca
        axis = np.cross(direction, [0., 0., 1.])
        axis /= np.linalg.norm(axis)
        old = ns["_aimed_rod"](ca, axis, 18, target, htype="CB")
        new = ns["_aimed_rod"](ca, axis, 18, target, htype="TRP")
        before = ns["trp_target_reachability"](old[1], target)
        after = ns["trp_target_reachability"](new[1], target)
        bound = after["target_CA_distance_A"] - after["maximum_CA_reach_A"]
        self.assertGreater(before["minimum_NE1_error_A"], bound + 1.)
        self.assertAlmostEqual(after["minimum_NE1_error_A"], bound, places=9)
        self.assertGreater(bound, .3)  # NE1 still cannot individually pass .30 A.
        np.testing.assert_array_equal(target, tz.ne1)
        np.testing.assert_array_equal(ca, tz.ca_trp)
        np.testing.assert_array_equal(new[1]["CA"], ca)
        p = np.array([v for r in old for v in r.values()])
        q = np.array([v for r in new for v in r.values()])
        np.testing.assert_allclose(np.linalg.norm(p[:, None]-p, axis=-1),
                                   np.linalg.norm(q[:, None]-q, axis=-1), atol=1e-12)
        # Chi-only numerical inverse placement independently attains the bound.
        best = min((least_squares(lambda chi: ns["build_sidechain"](
            "TRP", new[1], chi, 1.)["NE1"]-target, seed, max_nfev=70,
            ftol=1e-10, xtol=1e-10, gtol=1e-10)
                    for seed in ([0., 0.], [90., 90.], [-90., 180.])),
                   key=lambda fit: np.linalg.norm(fit.fun))
        self.assertAlmostEqual(np.linalg.norm(best.fun), bound, places=6)
        # Same torsions before/after: every sidechain pair distance and chirality preserved.
        sc0 = ns["build_sidechain"]("TRP", old[1], best.x, 1.)
        sc1 = ns["build_sidechain"]("TRP", new[1], best.x, 1.)
        x, y = np.array(list(sc0.values())), np.array(list(sc1.values()))
        np.testing.assert_allclose(np.linalg.norm(x[:, None]-x, axis=-1),
                                   np.linalg.norm(y[:, None]-y, axis=-1), atol=1e-12)
        def chirality(r, sc):
            return np.linalg.det(np.stack([r["N"]-r["CA"], r["C"]-r["CA"], sc["CB"]-r["CA"]]))
        self.assertAlmostEqual(chirality(old[1], sc0), chirality(new[1], sc1), places=12)
        print(f"fixture old_band_min={before['minimum_NE1_error_A']:.9f} "
              f"optimized_band_min={after['minimum_NE1_error_A']:.9f} "
              f"chi_fit={np.linalg.norm(best.fun):.9f} "
              f"conditional_six_anchor_RMSD={bound/math.sqrt(6):.9f}", flush=True)

    def test_interior_and_unattainable_roll_keep_original(self):
        ns, tz = self.ns, self.tz
        rod = ns["place_helix_anchor"](tz.ca_trp, [0., 0., 1.], 18, 1, 0.)
        self.assertIs(ns["trp_optimal_anchor_roll"](rod, [0., 0., 1.], tz.ca_trp), rod)
        axis = ns["place_cb"](rod[1])-tz.ca_trp
        self.assertIs(ns["trp_optimal_anchor_roll"](rod, axis, tz.ca_trp + 20*axis), rod)

    def test_twelve_perpendicular_axes_attain_same_radial_bound(self):
        ns, tz = self.ns, self.tz
        direction = tz.ne1-tz.ca_trp
        direction /= np.linalg.norm(direction)
        axis0 = np.cross(direction, [0., 0., 1.])
        axis0 /= np.linalg.norm(axis0)
        for angle in np.linspace(-math.pi, math.pi, 12, endpoint=False):
            axis = ns["rotate_about"](axis0, direction, angle)
            rod = ns["_aimed_rod"](tz.ca_trp, axis, 18, tz.ne1, htype="TRP")
            diagnostic = ns["trp_target_reachability"](rod[1], tz.ne1)
            self.assertAlmostEqual(diagnostic["minimum_NE1_error_A"],
                                   diagnostic["target_CA_distance_A"]-
                                   diagnostic["maximum_CA_reach_A"], places=9)

    def test_production_wiring_and_collision_filter_retained(self):
        function = next(n for n in self.tree.body if isinstance(n, ast.FunctionDef)
                        and n.name == "place_free_bundle")
        trp_calls = [n for n in ast.walk(function) if isinstance(n, ast.Call)
                     and isinstance(n.func, ast.Name) and n.func.id == "_aimed_rod"
                     and "tz.ca_trp" in ast.unparse(n)]
        self.assertEqual(len(trp_calls), 1)
        self.assertIn("htype='TRP'", ast.unparse(trp_calls[0]))
        ns, tz = self.ns, self.tz
        axis = np.cross(tz.ne1-tz.ca_trp, [0., 0., 1.])
        rod = ns["_aimed_rod"](tz.ca_trp, axis, 18, tz.ne1, htype="TRP")
        saved = ns["sub_heavy_s"]
        try:
            ns["sub_heavy_s"] = rod[10]["CA"][None]
            self.assertIsNone(ns["_aimed_rod"](tz.ca_trp, axis, 18, tz.ne1, htype="TRP"))
        finally:
            ns["sub_heavy_s"] = saved


if __name__ == "__main__":
    unittest.main()
