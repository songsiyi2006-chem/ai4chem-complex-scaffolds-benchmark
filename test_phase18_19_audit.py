"""Targeted source regressions, NumPy/SciPy only. No production imports."""
import ast
import math
from pathlib import Path
from types import SimpleNamespace as NS
import unittest
import numpy as np
from scipy.optimize import Bounds, LinearConstraint
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parent


def tree(i):
    return ast.parse(next(ROOT.glob(f'run_phase{i}_*.py')).read_text(encoding='utf-8'))


def load(i,names,ns):
    nodes = [n for n in ast.walk(tree(i)) if isinstance(n,ast.FunctionDef) and n.name in names]
    assert set(n.name for n in nodes)==set(names)
    exec(compile(ast.Module(body=nodes,type_ignores=[]),'<actual-source>','exec'),ns)
    return ns


class Tests(unittest.TestCase):
    def test_compile(self):
        for i in (18,19): compile(tree(i),str(i),'exec')

    def test_lp_only_relaxes_binary_bounds(self):
        seen = {}
        def milp(**kw): seen.update(kw); return NS(success=True)
        ns = load(18,['solve_lp'],dict(np=np,milp=milp,Bounds=Bounds,LinearConstraint=LinearConstraint))
        model = dict(net=NS(ridx={'BIOMASS':0}),nv=4,nR=1,nD=1,nZ=1,zoff=2,doff=3,
                     c=np.array([-1.,0,0,0]),lbv=np.array([0.,-8.,0.,-2000.]),
                     ubv=np.array([10.,-1.,1.,2000.]),A_eq=sp.csr_matrix((1,4)),
                     A_ub=sp.csr_matrix((1,4)),b_ub=np.ones(1))
        ns['solve_lp'](model)
        np.testing.assert_equal(seen['bounds'].lb,model['lbv'])
        np.testing.assert_equal(seen['bounds'].ub,model['ubv'])

    def test_gated_thermodynamics_excludes_boundary_and_uses_effective_eps(self):
        ns = load(18,['extract'],dict(np=np,RT=1.,EPS_T=.5,
                   zero_loop_certificate=lambda *a: {'pass':True}))
        net = NS(rxns=[dict(kind='enz',stoich={'a':1},parent='E',dirn='F'),
                       dict(kind='boundary',stoich={'b':1},parent='B',dirn='F')],
                 midx={'a':0,'b':1},crowd=np.ones(2))
        model = dict(nR=2,nD=0,doff=2,didx={},dmet=[],unknown_mets=[],nZ=0,
                     n_blocked=0,effective_eps=.02)
        out = ns['extract'](net,np.array([-.025,10.]),np.array([1.,1.]),model,1.,1.)
        self.assertEqual(out['n_violations'],0)
        self.assertEqual(out['worst_dG_active'],-.025)

    def test_qm_complete_fragment_embedding_and_peak_geometry(self):
        seen = []
        class Calc:
            def __init__(self,numbers,q,mp,mq,charge,solvent):
                self.numbers=numbers; self.qm_charges=q; self.charge=charge
                self.solvent=solvent; self.n_calls=0
            def _embedding(self,pos):
                seen.append(self.qm_charges.copy())
                return 0.,np.zeros_like(pos)
        def optimize(numbers,pos,i,j,target,charge,solvent,maxcyc,return_charges=False):
            p=pos.copy()
            if not return_charges: return p,0.
            k=len(seen); p[j]=[2.+k,0.,0.]
            return p,[0.,3.,1.][k],np.full(len(numbers),.1)
        ace_sym=['C','C','O','O','H','H','H']
        ace=np.array([[0,0,0],[1,0,0],[1,1,0],[2,0,0],[0,0,1],[0,1,0],[0,0,-1]],float)
        ns = load(19,['qmmm_kemp_scan'],dict(np=np,math=math,
             _cached_geom=lambda *a:(ace,ace_sym),_atomic_numbers=lambda syms:[1]*len(syms),
             XtbMMCalculator=Calc,_xtb_constrained_opt=optimize,
             _amber_charges=lambda enz,syms,pos,skip:(np.ones((1,3)),np.ones(1),np.zeros(len(syms))),
             tassert=lambda cond,msg: self.assertTrue(cond,msg)))
        tz = NS(pos_sub=np.array([[0,1,0],[0,2,0],[0,0,0],[1,0,0]],float),
                sym=['O','N','C','H'],i_H3=3,i_C3=2,i_N8=1,i_O7=0,
                o_base=np.array([3.5,0.,0.]),zhat=np.array([0.,0.,1.]))
        out=ns['qmmm_kemp_scan'](tz,{},[],3)
        self.assertEqual(out['qm_atoms'],11)  # 4 substrate + complete 7-atom acetate
        self.assertEqual(out['barrier_kcal'],3.)
        self.assertEqual(out['ts_d_CH'],3.)  # peak, not final point at 4 A
        self.assertTrue(all(np.all(q==.1) for q in seen))
        self.assertFalse(out['transition_state_validated'])

    def test_embedding_gradient(self):
        f=load(19,['_embedding'],{'np':np})['_embedding']
        engine=NS(mm_pos=np.array([[3.,0.,0.]]),mm_charges=np.array([-.2]),
                  qm_charges=np.array([.3]),cutoff=12.)
        p=np.zeros((1,3)); e,g=f(engine,p)
        h=1e-5; pp=p.copy(); pm=p.copy(); pp[0,0]+=h; pm[0,0]-=h
        numeric=(f(engine,pp)[0]-f(engine,pm)[0])/(2*h)
        self.assertAlmostEqual(g[0,0],numeric,places=7)
        self.assertNotEqual(e,0.)

    def test_dynamic_conserved_moieties(self):
        ns=dict(np=np)
        names={'POOLS','PIDX','NPOOL','ESETS','EIDX','NY','O2_SAT','Y_AA'}
        nodes=[n for n in tree(18).body if isinstance(n,ast.Assign)
               and any(isinstance(t,ast.Name) and t.id in names for t in n.targets)]
        exec(compile(ast.Module(body=nodes,type_ignores=[]),'<constants>','exec'),ns)
        ns['mu_of']=lambda *args: .01
        load(18,['sat','dyn_rhs'],ns)
        P=dict(ngam=1.,gam=2.,h2o2_rate=lambda t:0.,krela=.1,kspot=.1,
               kcyc=.1,kpde=.1,W=.1,kdeg=.1)
        rng=np.random.default_rng(18)
        for _ in range(8):
            y=rng.uniform(.1,2.,ns['NY']); dy=ns['dyn_rhs'](0.,y,P)
            for group in [('atp','adp','amp'),('nadh','nad'),('nadph','nadp'),
                          ('q','qh2'),('coa','accoa')]:
                self.assertAlmostEqual(sum(dy[ns['PIDX'][m]] for m in group),0.,places=11)


if __name__=='__main__': unittest.main()
