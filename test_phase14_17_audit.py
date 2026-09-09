"""Isolated actual-source regressions; NumPy/SciPy only, no long jobs."""
import ast
import json
import math
from pathlib import Path
import tempfile
from types import SimpleNamespace as NS
import unittest
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply

ROOT = Path(__file__).resolve().parent


def tree(i):
    return ast.parse(next(ROOT.glob(f'run_phase{i}_*.py')).read_text(encoding='utf-8'))


def load(i, names, ns):
    nodes = [n for n in ast.walk(tree(i)) if isinstance(n, ast.FunctionDef) and n.name in names]
    assert set(n.name for n in nodes) == set(names)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), '<actual-source>', 'exec'), ns)
    return ns


class Tests(unittest.TestCase):
    def test_compile(self):
        for i in range(14, 18): compile(tree(i), f'phase{i}', 'exec')

    def test_reduced_entropy_units(self):
        t = tree(14)
        d = next(n for n in ast.walk(t) if isinstance(n, ast.Dict)
                 and any(isinstance(k, ast.Constant) and k.value == 'S_total' for k in n.keys))
        expr = next(v for k, v in zip(d.keys, d.values)
                    if isinstance(k, ast.Constant) and k.value == 'S_total')
        val = eval(compile(ast.Expression(expr), '<actual-source>', 'eval'),
                   dict(s_diff=0., s_chem=19.4, T_K=310.))
        self.assertEqual(val, 19.4)
        self.assertNotIn('S_total_kB_um2_s', ast.unparse(t))

    def test_singlet_trace_and_nonnegative_populations(self):
        ns = dict(np=np, sp=sp, expm_multiply=expm_multiply,
                  CONFIG=dict(K_S=1e6, K_T=1e6, XI_E=0., T_MAX_US=1., N_T=20))
        nodes = [n for n in tree(15).body if
                 (isinstance(n, (ast.FunctionDef, ast.ClassDef)) and
                  n.name in ('spin_ops', 'RadicalPairSpinSystem')) or
                 (isinstance(n, ast.Assign) and any(isinstance(a, ast.Name)
                  and a.id in ('SX', 'SY', 'SZ') for a in n.targets))]
        exec(compile(ast.Module(body=nodes, type_ignores=[]), '<actual-source>', 'exec'), ns)
        s = ns['RadicalPairSpinSystem']([], [])
        t, ps, pt, surv, ys, yt = s.propagate(sp.csr_matrix((4,4), dtype=complex))
        np.testing.assert_allclose(surv, np.exp(-1e6*t), atol=1e-12)
        np.testing.assert_allclose(ps, surv, atol=1e-12)
        np.testing.assert_allclose(pt, 0., atol=1e-12)
        self.assertGreater(ys, 0.)
        self.assertAlmostEqual(yt, 0.)

    def test_wham_known_harmonic_pmf_and_latch_units(self):
        ns = dict(np=np, math=math, KB=1.380649e-23, N_AVOGADRO=6.02214076e23,
                  KCAL_MOL_J=4184., CONFIG=dict(MD_TEMP_K=310., UMB_R_RELEASE_NM=2.4))
        f = load(15, ['wham_pmfs'], ns)['wham_pmfs']
        beta = 4184/(ns['KB']*310*ns['N_AVOGADRO'])
        centers = .5*(np.linspace(.5,3.,29)[:-1]+np.linspace(.5,3.,29)[1:])
        true = .4*(centers-1.6)**2
        windows = []
        for r0 in (1., 1.7, 2.4):
            p = np.exp(-beta*(true+.5*(4./4.184)*(centers-r0)**2))
            counts = np.maximum(1, np.rint(100000*p/p.sum())).astype(int)
            windows.append(dict(k=4., r0_nm=r0, samples_nm=np.repeat(centers, counts)))
        out = f(dict(FAD_oxid=dict(windows=windows,latch_occupied_frac=.2),
                     FAD_radan=dict(windows=windows,latch_occupied_frac=.4)))
        got = np.array(out['FAD_oxid']['pmf_raw'])
        np.testing.assert_allclose(got, true-true.min(), atol=.003)
        self.assertAlmostEqual(out['allostery']['dG_latch_shift_kcal'], -math.log(2)/beta)
        with self.assertRaises(ValueError):
            f(dict(FAD_oxid=dict(windows=[dict(k=4.,r0_nm=1.,samples_nm=[])],latch_occupied_frac=.2)))

    def test_three_lj_force_paths(self):
        ns = load(16, ['pair_forces_local','pair_forces_split','pair_forces'],
                  dict(np=np, EPS_MAT=np.ones((1,1)), SIGMA=np.ones(1)))
        eng = NS(box=100., pl_lj=NS(i=np.array([0]),j=np.array([1]),rc2=100.),
                 sub_code=np.zeros(2,dtype=int), fcap=1e9, pl_idx=np.arange(2),
                 p_eps=np.ones(1),p_sig=np.ones(1),q_pairs=([],[],[]))
        eng._scatter = lambda frc,idx,vec: np.add.at(frc,idx,vec)
        for r in (1.1, 1.5, 2.):
            pos = np.array([[r,0.,0.],[0.,0.,0.]])
            v = lambda z: 4*(z**-12-z**-6)
            fd = -(v(r+1e-6)-v(r-1e-6))/2e-6
            local = ns['pair_forces_local'](eng,pos,np.zeros_like(pos),[0],np.array([0,-1]))
            frc = np.zeros_like(pos)
            split = ns['pair_forces_split'](eng,pos,frc,[0])
            self.assertAlmostEqual(local[0,0], fd, places=6)
            self.assertAlmostEqual(split[0,0], fd, places=6)
            self.assertAlmostEqual(frc[1,0], -fd, places=6)
            frc[:] = 0
            ns['pair_forces'](eng,pos,frc)
            self.assertAlmostEqual(frc[0,0], fd, places=6)
            np.testing.assert_allclose(frc.sum(0),0.)

    def test_transport_sign_bulk_reference_and_si_units(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td); (root/'umbrella').mkdir()
            ns = dict(np=np,math=math,json=json,RES=root,KT=1.,BETA=1.,
                      TRAPZ=np.trapezoid,ETA_WATER=.001,NM2PS_TO_UM2S=1e6,KJ_TO_KCAL=1/4.184)
            f = load(16,['analyze_transport'],ns)['analyze_transport']
            z = np.linspace(-25,25,21)
            for kind in ('receptor','inert'):
                force = np.repeat((.032*z)[:,None],20,axis=1)
                np.savez(root/'umbrella'/f'windows_{kind}.npz',z0=z,force=force,
                         contacts=np.zeros_like(force),eq_contact_max=np.zeros(len(z)),
                         gamma=np.ones(len(z)),M_cx=1.,rh=2.)
            out = f({},True)
            with np.load(root/'umbrella'/'pmf_receptor.npz') as d:
                expected_g = 10*(1-(z/25)**2)
                expected_g -= expected_g[np.abs(z)>22].mean()
                np.testing.assert_allclose(d['G'],expected_g,atol=1e-12)
                self.assertLess(d['G'].min(),0.)  # negative wells are retained
            expected_d = 1.6605e-21/(6*math.pi*.001*2e-9)*1e12
            self.assertAlmostEqual(out['receptor']['D_cal_factor'],expected_d/1e6)
            self.assertAlmostEqual(out['P_ratio_receptor_over_inert'],1.)

    def test_exchange_derivative_and_double_counting(self):
        ns = load(17,['xalpha_vx','slater_x_energy','hartree_potential',
                      'total_energy_from_eigenvalues'],dict(np=np,XALPHA=.7))
        x = np.array([1.,2.,3.]); w = np.ones(3)*.1; rho = np.array([.2,.1,.05])
        dv = 4*np.pi*x*x*w
        vx = ns['xalpha_vx'](rho)
        for i in range(3):
            rp,rm = rho.copy(),rho.copy(); rp[i]+=1e-6; rm[i]-=1e-6
            numeric = (ns['slater_x_energy'](rp,x,w)-ns['slater_x_energy'](rm,x,w))/2e-6/dv[i]
            self.assertAlmostEqual(numeric,vx[i],places=8)
        u = ns['hartree_potential'](rho,x,w)
        ex = ns['slater_x_energy'](rho,x,w)
        expected = -5.-.5*np.sum(rho*u*dv)-ex/3
        actual = ns['total_energy_from_eigenvalues'](-5.,rho,u,vx,x,w)
        self.assertAlmostEqual(actual,expected)


if __name__ == '__main__': unittest.main()
