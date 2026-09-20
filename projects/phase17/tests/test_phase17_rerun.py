"""Cheap SCF control-flow tests; no physical heavy-atom computation."""
import ast
import math
from pathlib import Path
import unittest
import numpy as np

SOURCE = Path(__file__).with_name("run_phase17_relativistic_actinide_quantum.py")


def scf_namespace(states=None):
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
    names = {"SCFConvergenceError", "scf_residuals", "atomic_scf", "orbital_density",
             "normalize_scf_orbitals", "scf_density_from_orbitals",
             "_alpha_mats", "angular_pair_contraction", "breit_pair_energy",
             "hydrogenic_radial_grid", "coulomb_pair_energy_ref"}
    nodes = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))
             and n.name in names]
    x = np.array([1., 2.])
    p = np.sqrt([.5, .5])
    ns = dict(np=np, math=math, C_LIGHT=137., NR_DR=1.,
              nr_grid=lambda: x, nucleus_potential=lambda *a, **k: (np.zeros(2), 0.),
              config_blocks=lambda c: ({(1, -1): 2.}, [-1]),
              guess_density=lambda *a: (2*p*p/(4*np.pi*x*x), 2.),
              nr_bound_states=states or (lambda *a, **k: [(-1., p)]),
              hartree_potential=lambda *a: np.zeros(2), xalpha_vx=lambda r: np.zeros(2),
              total_energy_from_eigenvalues=lambda *a: -2., log=lambda s: None,
              _orbital_record=lambda *a: {})
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SOURCE), "exec"), ns)
    return ns


class SCFTests(unittest.TestCase):
    def test_hydrogenic_grid_scaling_and_refinement(self):
        ns = scf_namespace()
        errors = []
        for count in [256, 512]:
            scaled_energies = []
            for z in [1.,80.]:
                x,w = ns['hydrogenic_radial_grid'](z,count)
                p=2*z**1.5*x*np.exp(-z*x)
                energy=ns['coulomb_pair_energy_ref'](p,x,w)
                scaled_energies.append(energy/z)
            np.testing.assert_allclose(scaled_energies[0],scaled_energies[1],rtol=1e-13)
            errors.append(abs(scaled_energies[0]-.625))
        self.assertLess(errors[1], errors[0]/3)

    def test_sparse_angular_matches_dense_kernels(self):
        ns = scf_namespace()
        mats = ns["_alpha_mats"]()
        K9 = np.array([np.kron(a,b) for a in mats for b in mats])
        rng = np.random.default_rng(17)
        psi = rng.normal(size=(23,16)) + 1j*rng.normal(size=(23,16))
        direction = rng.normal(size=(23,3))
        direction /= np.linalg.norm(direction,axis=1)[:,None]
        wr = rng.random(23)
        aa = K9[0]+K9[4]+K9[8]
        T = np.einsum('pi,pj,ijab->pab',direction,direction,K9.reshape(3,3,16,16))
        for kernel, op in [('coulomb',np.broadcast_to(np.eye(16),(23,16,16))),
                           ('gaunt',np.broadcast_to(-aa,(23,16,16))),
                           ('breit',-.5*(aa+T))]:
            expected = np.einsum('p,pa,pab,pb->ab',wr,psi.conj(),op,psi)
            actual = ns['angular_pair_contraction'](psi,direction,wr,kernel,K9)
            np.testing.assert_allclose(actual,expected,rtol=2e-13,atol=2e-13)

    def test_radial_ratio_and_chunk_independence(self):
        ns = scf_namespace()
        A = np.zeros((4,4,4,4,2),complex)
        A[0,0,0,0] = [0.,1.]
        tab = {'u':np.array([0.,1.]),'A_dir':A}
        for chunk in [1,2,128]:
            energy=ns['breit_pair_energy'](tab,np.ones(2),None,np.array([1.,4.]),np.ones(2),chunk=chunk)
            self.assertAlmostEqual(energy.real,1.375)
            self.assertAlmostEqual(energy.imag,0.)

    def test_grouped_radial_matches_dense_reference(self):
        ns=scf_namespace()
        rng=np.random.default_rng(170)
        x=np.array([.1,.3,.8,1.7])
        w=np.array([.1,.2,.4,.5])
        P=rng.normal(size=4); Q=rng.normal(size=4)
        u=np.linspace(0,1,5)
        A=rng.normal(size=(4,4,4,4,5))+1j*rng.normal(size=(4,4,4,4,5))
        ratio=np.minimum.outer(x,x)/np.maximum.outer(x,x)
        fr=[P,P,Q,Q]; expected=0j
        for a,b,g,d in np.ndindex(4,4,4,4):
            vals=np.interp(ratio,u,A[a,b,g,d].real)+1j*np.interp(ratio,u,A[a,b,g,d].imag)
            expected+=np.sum((fr[a]*fr[g]*w)[:,None]*(fr[b]*fr[d]*w)[None,:]*vals/np.maximum.outer(x,x))
        actual=ns['breit_pair_energy']({'u':u,'A_dir':A},P,Q,x,w,chunk=2)
        np.testing.assert_allclose(actual,expected,rtol=2e-13,atol=2e-13)

    def test_orbital_quadrature_normalization(self):
        ns = scf_namespace()
        w = np.array([.0003, .0003])
        original = {(1,-1): (-1., np.sqrt([.5,.5]), None),
                    (2,-1): (-.5, np.array([2.,3.]), np.array([.4,.2]))}
        normalized = ns["normalize_scf_orbitals"](original, w)
        for key, (energy, p, q) in normalized.items():
            self.assertAlmostEqual(np.sum((p*p+(0 if q is None else q*q))*w), 1.)
            self.assertEqual(energy, original[key][0])
        with self.assertRaises(FloatingPointError):
            ns["normalize_scf_orbitals"]({(1,-1): (-1.,np.zeros(2),None)}, w)

    def test_scf_final_density_independent_of_eigenvector_scale(self):
        ns = scf_namespace(lambda *a, **k: [(-1., np.sqrt([.5,.5]) * .03)])
        atom = ns["atomic_scf"](2, [(1,0,2.,0.)], relativistic=False)
        self.assertTrue(atom["scf"]["converged"])
        self.assertLess(atom["scf"]["final"]["charge_error"], 1e-12)

    def test_unnormalized_density_is_rejected_not_rescaled(self):
        ns = scf_namespace()
        with self.assertRaises(FloatingPointError):
            ns["scf_density_from_orbitals"]({(1,-1): (-1.,np.ones(2),None)},
                                             {(1,-1): 2.},np.array([1.,2.]),np.ones(2))

    def test_failed_candidate_can_continue_to_convergence(self):
        calls = [0]
        def states(*a, **kw):
            calls[0] += 1
            return [(-1., np.sqrt([.9,.1] if calls[0] == 10 else [.5,.5]))]
        ns = scf_namespace(states)
        atom = ns["atomic_scf"](2,[(1,0,2.,0.)],relativistic=False,max_iter=12)
        self.assertTrue(atom["scf"]["converged"])
        self.assertEqual(atom["iterations"],10)
        self.assertEqual(len(atom["scf"]["verification_history"]),2)

    def test_constant_energy_does_not_hide_density_cycle(self):
        calls = [0]
        def states(*a, **kw):
            calls[0] += 1
            return [(-1., np.sqrt([.9, .1] if calls[0] % 2 else [.1, .9]))]
        ns = scf_namespace(states)
        with self.assertRaises(ns["SCFConvergenceError"]) as caught:
            ns["atomic_scf"](2, [(1,0,2.,0.)], relativistic=False, max_iter=12)
        d = caught.exception.diagnostics
        self.assertFalse(d["converged"])
        self.assertEqual(d["iterations"], 12)
        self.assertEqual(d["history"][-1]["energy_change_eh"], 0.)
        self.assertGreater(d["final"]["density_l1_per_electron"], 3e-6)

    def test_fixed_point_passes_and_initial_delta_is_unknown(self):
        ns = scf_namespace()
        atom = ns["atomic_scf"](2, [(1,0,2.,0.)], relativistic=False)
        self.assertTrue(atom["scf"]["converged"])
        self.assertIsNone(atom["scf"]["history"][0]["energy_change_eh"])
        self.assertLess(atom["scf"]["final"]["density_l1_per_electron"], 3e-6)

    def test_final_solve_must_also_converge(self):
        calls = [0]
        def states(*a, **kw):
            calls[0] += 1
            return [(-1., np.sqrt([.5,.5] if calls[0] <= 9 else
                                   ([.9,.1] if calls[0] % 2 else [.1,.9])))]
        ns = scf_namespace(states)
        with self.assertRaises(ns["SCFConvergenceError"]):
            ns["atomic_scf"](2, [(1,0,2.,0.)], relativistic=False, max_iter=12)

    def test_orbital_failure_is_not_accepted(self):
        ns = scf_namespace()
        ns.update(log_grid=lambda *a: (np.array([1.,2.]), None,None,np.ones(2)),
                  l_of_kappa=lambda k: 0,
                  dirac_orbital_shoot=lambda *a, **k: (-1., np.ones(2),np.zeros(2),False))
        with self.assertRaisesRegex(RuntimeError, "shooting failed"):
            ns["atomic_scf"](2, [(1,0,2.,0.)])

    def test_invalid_iteration_budget(self):
        with self.assertRaises(ValueError):
            scf_namespace()["atomic_scf"](2, [], max_iter=0)


class PostauditSafetyTests(unittest.TestCase):
    def test_kb_positive_branch_and_basis_refinement(self):
        import scipy.linalg as sla
        tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
        ns = {"np": np, "math": math, "sla": sla, "C_LIGHT": 137.035999084,
              "log": lambda *a: None}
        selected = [n for n in tree.body if isinstance(n, ast.FunctionDef)
                    and n.name in ("log_grid", "sommerfeld_1s", "gaussian_kb_demo")]
        exec(compile(ast.Module(body=selected, type_ignores=[]), str(SOURCE), "exec"), ns)
        errors = []
        for n in (10, 24, 32):
            _, exact, diagnostics = ns["gaussian_kb_demo"](n_basis=n, return_diagnostics=True)
            d = diagnostics["balanced"]
            self.assertLess(d["lowest_spectrum_eh"], -ns["C_LIGHT"]**2)
            self.assertGreater(d["bound_1s_eh"], -ns["C_LIGHT"]**2)
            self.assertLess(d["bound_1s_eh"], 0)
            self.assertLess(d["eigen_residual_relative"], 1e-8)
            errors.append(abs(d["bound_1s_eh"]-exact)/abs(exact))
        self.assertLess(errors[-1], 2e-4)
        self.assertLess(errors[-1], errors[1])
        self.assertLess(errors[1], errors[0]/50)

    def test_axisymmetric_laplacian_and_axis_descriptor(self):
        tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
        ns = {"np": np, "math": math}
        selected = [n for n in tree.body if isinstance(n, ast.FunctionDef)
                    and n.name in ("_laplacian_plane", "promolecular_axis_descriptor")]
        exec(compile(ast.Module(body=selected, type_ignores=[]), str(SOURCE), "exec"), ns)
        x = np.linspace(-2, 4, 121)
        z = np.linspace(0, 3, 61)
        rho = 1 + (x[:, None]-1)**2 + 2*z[None, :]**2
        lap = ns["_laplacian_plane"](rho, x[1]-x[0], z[1]-z[0])
        np.testing.assert_allclose(lap[1:-1, :-1], 10, atol=2e-11)
        self.assertTrue(np.isnan(lap[0]).all())
        q = ns["promolecular_axis_descriptor"](rho, x, z, 3)
        self.assertAlmostEqual(q["rho_bcp"], 1)
        self.assertAlmostEqual(q["lap_bcp"], 10)
        self.assertAlmostEqual(q["axis_gradient_residual"], 0)

    def test_merge_rejects_unrelated_change(self):
        from audit_phase17_qtaim_figure3 import validate_merge
        prior = dict(atoms={"unchanged": 1}, breit={"Gaunt": 4, "hydrogenic_breit_cm": [1]})
        candidate = dict(atoms={"unchanged": 1}, breit={"Gaunt": 4, "hydrogenic_breit_cm": [2]})
        validate_merge(prior, candidate)
        candidate["breit"]["Gaunt"] = 5
        with self.assertRaises(RuntimeError):
            validate_merge(prior, candidate)

    def test_completion_gate(self):
        from audit_phase17_scaled_hydrogenic import require_completed
        valid = dict(status="process_completed", stages=[
            dict(stage=s, exit_code=0) for s in ("breit", "multiplet", "figures")])
        require_completed(valid)
        with self.assertRaises(RuntimeError):
            require_completed(dict(valid, status="running"))
        valid["stages"][-1]["exit_code"] = 1
        with self.assertRaises(RuntimeError):
            require_completed(valid)

    def test_scaled_coulomb_retains_measured_error(self):
        from audit_phase17_scaled_hydrogenic import derive_coulomb
        rows = derive_coulomb([dict(Z=z, J_quadrature=.622) for z in (1,20,80)])
        for row in rows:
            self.assertAlmostEqual(row["J_quadrature"], .622*row["Z"])
            self.assertAlmostEqual(row["rel_err"], .003/.625)
            self.assertIn("derived", row["provenance"])

    def test_immutable_prior(self):
        from audit_phase17_scaled_hydrogenic import preserve
        import tempfile
        import stat
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"prior.json"
            try:
                preserve(path, b"original")
                preserve(path, b"original")
                with self.assertRaises(RuntimeError):
                    preserve(path, b"changed")
                self.assertEqual(path.read_bytes(), b"original")
            finally:
                if path.exists():
                    path.chmod(stat.S_IWRITE | stat.S_IREAD)

    def test_anchor_implementation_verification(self):
        from audit_phase17_scaled_hydrogenic import verify_anchor_source
        code = b"def breit_pair_energy(): pass\ndef breit_pair_tabulation(): pass\ndef angular_pair_contraction(): pass\n"
        verify_anchor_source(code, code)
        with self.assertRaises(RuntimeError):
            verify_anchor_source(code, code.replace(b"pass", b"return 1", 1))


if __name__ == "__main__":
    unittest.main()
