"""Fast regressions: python -m unittest -v test_phase1_5_audit.

RDKit and OpenMM are required for integration tests; no xTB/MD jobs or shutdown.
AST extraction executes the real function bodies without legacy import side effects.
"""
import ast
import json
import math
from pathlib import Path
import tempfile
import traceback
import types
import sys
import unittest
from unittest.mock import patch

import numpy as np
from phase_audit import (AUDIT_VERSION, converged_ensemble, periodic_angle_error,
                         ts_validation, exploratory_corrections, cross_nonbonded_system)

ROOT = Path(__file__).resolve().parent


def load_functions(file, names, namespace):
    tree = ast.parse((ROOT/file).read_text(encoding='utf-8'))
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    assert {n.name for n in nodes} == set(names)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), file, 'exec'), namespace)
    return namespace


class AuditTests(unittest.TestCase):
    def test_all_sources_parse(self):
        for name in ('molecule_benchmark.py','run_heavy_dynamics_benchmark.py',
                     'run_phase3_complex_dynamics.py','run_phase4_reaction_mechanism.py',
                     'run_phase5_chemical_world_model.py'):
            ast.parse((ROOT/name).read_text(encoding='utf-8'))

    def test_conformer_filter(self):
        ids,energies=converged_ensemble([0,1,2],[(1,-100),(0,5),(0,float('nan'))])
        self.assertEqual(ids,[1]); self.assertEqual(energies,[5])
        with self.assertRaises(ValueError): converged_ensemble([0],[(1,0)])

    def test_real_phase1_function_rejects_nonconverged(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        ns=load_functions('molecule_benchmark.py',['embed_and_optimize'],
                          dict(Any=object,Dict=dict,SEED=2026,AllChem=AllChem))
        with patch.object(AllChem,'MMFFOptimizeMoleculeConfs',return_value=[(1,-100.),(0,1.)]):
            result=ns['embed_and_optimize'](Chem.MolFromSmiles('CCCC'),2)
        self.assertEqual(result['conf_ids'],[1]);self.assertEqual(result['e_min'],1.)
        with patch.object(AllChem,'MMFFOptimizeMoleculeConfs',return_value=[(1,-100.)]):
            result=ns['embed_and_optimize'](Chem.MolFromSmiles('CCCC'),1)
        self.assertFalse(result['ok'])

    def test_circular_angles(self):
        self.assertEqual(periodic_angle_error(-170,190),0.)
        self.assertEqual(periodic_angle_error(-179,179),2.)

    def test_phase2_minimum_saved_and_unbiased(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        ns=dict(Path=Path,List=list,Any=object,Dict=dict,math=math,traceback=traceback,
                REGISTRY=[dict(id='TEST',smiles='CCCC')],SEED=2026,
                RESULTS={'stage2_torsion_scans':{}},_hr=lambda *a:None,_log=lambda *a:None,
                _pick_scan_torsion=lambda m:dict(torsion=(0,1,2,3),mode='test'))
        load_functions('run_heavy_dynamics_benchmark.py',['stage2_torsion_scans'],ns)
        with tempfile.TemporaryDirectory() as td:
            ns['stage2_torsion_scans'](Path(td),12)
            rec=ns['RESULTS']['stage2_torsion_scans']['TEST']
            self.assertEqual(rec['status'],'ok',rec)
            mol=Chem.SDMolSupplier(str(Path(td)/'TEST_scan_min.sdf'),removeHs=False)[0]
            angle=AllChem.GetDihedralDeg(mol.GetConformer(),0,1,2,3)
            self.assertLess(periodic_angle_error(angle,rec['angle_of_minimum_deg']),5.2)
            self.assertGreater(periodic_angle_error(angle,330),20)
            props=AllChem.MMFFGetMoleculeProperties(mol)
            energy=AllChem.MMFFGetMoleculeForceField(mol,props).CalcEnergy()
            self.assertAlmostEqual(energy,float(mol.GetProp('energy_kcal_mol')),places=3)
            self.assertFalse(rec['energy_includes_restraint'])

    def test_openmm_cross_lj_and_coulomb(self):
        from openmm import NonbondedForce, Context, VerletIntegrator, Platform, unit
        nb=NonbondedForce()
        params=[(.4,.25,.8),(-.3,.31,.5),(.2,.27,.4)]
        for values in params: nb.addParticle(*values)
        system=cross_nonbonded_system(nb,[0,1])
        integrator=VerletIntegrator(.001)
        ctx=Context(system,integrator,Platform.getPlatformByName('Reference'))
        positions=np.array([[0,0,0],[.4,0,0],[0,.6,0]])
        ctx.setPositions(positions)
        actual=ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        expected=0.
        for i in (0,1):
            r=np.linalg.norm(positions[i]-positions[2]);q,s,e=params[i];q2,s2,e2=params[2]
            a=((s+s2)/2/r)**6
            expected+=4*math.sqrt(e*e2)*(a*a-a)+138.935457644382*q*q2/r
        self.assertAlmostEqual(actual,expected,places=8)
        del ctx,integrator
        nb.addException(0,2,0,.25,0)
        with self.assertRaises(ValueError):cross_nonbonded_system(nb,[0,1])

    def test_ts_gate(self):
        valid=ts_validation([-200,100,200],[100,200],.01,.01,True,.8)
        self.assertTrue(valid['passed']);self.assertFalse(valid['irc_verified'])
        for ft,fr,force,neb,overlap in [([100], [100],.01,True,.8),
            ([-200,-100],[100],.01,True,.8),([-200],[100],.2,True,.8),
            ([-200],[-100],.01,True,.8),([-200],[100],.01,False,.8),
            ([-200],[100],.01,True,None)]:
            self.assertFalse(ts_validation(ft,fr,force,.01,neb,overlap)['passed'])

    def test_phase4_invalid_ts_blocks_figures(self):
        ns=load_functions('run_phase4_reaction_mechanism.py',['stage4_figures'],
                          dict(Path=Path,RESULTS={'stage3_ts':{}}))
        with self.assertRaises(RuntimeError):ns['stage4_figures'](Path('unused'))

    def test_phase4_failed_stationarity_withholds_rate(self):
        # Exercise the actual stage with realistic list-valued normal modes.
        class Atoms:
            def get_atomic_numbers(self): return np.array([1,1])
            def get_positions(self): return np.array([[0.,0,0],[1.,0,0]])
        class Calculator:
            def get_forces(self,atoms): return np.array([[1.,0,0],[-1.,0,0]])
        io=types.ModuleType('ase.io');io.read=lambda *a,**k:Atoms()
        hs=dict(frequencies=[-200.,100.,200.],modes=[dict(freq=-200.,disp=[[1.,0,0],[-1.,0,0]])])
        hr=dict(frequencies=[100.,200.,300.],modes=[])
        with tempfile.TemporaryDirectory() as td:
            ns=dict(Path=Path,json=json,np=np,XTB_EXE='mock',XTBWrap=Calculator,
                    xtb_hessian=unittest.mock.Mock(side_effect=[hs,hr]),
                    RESULTS={'stage1_geometry':{'broken_bonds_R_idx':[[0,1]],'formed_bonds_R_idx':[]},
                             'stage2_neb':{'converged':True}},
                    _log=lambda *a:None,_warn=lambda *a:None)
            load_functions('run_phase4_reaction_mechanism.py',['stage3_ts_verify'],ns)
            with patch.dict(sys.modules,{'ase':types.ModuleType('ase'),'ase.io':io}):
                result=ns['stage3_ts_verify'](Path(td),types.SimpleNamespace(fmax=.05))
            self.assertFalse(result)
            record=json.loads((Path(td)/'stage3.json').read_text())
            self.assertIsNone(record['thermochemistry_298K'])
            self.assertFalse(record['validation']['passed'])

    def test_phase5_report_no_proven_claim(self):
        with tempfile.TemporaryDirectory() as td:
            ns=dict(ROOT=Path(td),RESULTS={},AUDIT_VERSION=AUDIT_VERSION,json=json)
            load_functions('run_phase5_chemical_world_model.py',['_write_audited_report'],ns)
            ns['_write_audited_report']('EN')
            report=(Path(td)/'WORLD_MODEL_REPORT_EN.md').read_text(encoding='utf-8')
            self.assertIn('"barrier_reduction_proven": false',report)
            self.assertNotIn('CLAIM PROVEN',report)

    def test_rrho_entropy_units(self):
        ns=load_functions('run_phase5_chemical_world_model.py',['_rrho_gibbs_fallback'],
                          dict(math=math,T_REF=298.15))
        actual=ns['_rrho_gibbs_fallback'](0,[],3)
        expected=-298.15*(26.34+21.34)/1000/627.5095
        self.assertAlmostEqual(actual,expected,places=14)

    def test_scenario_signs_and_invalid(self):
        a=exploratory_corrections(-20,-3)
        self.assertEqual(a['d1_drop'],-8);self.assertEqual(a['ddG_stereo'],-1.5)
        with self.assertRaises(ValueError):exploratory_corrections(float('nan'),1)

    def test_module_d_requires_optin(self):
        ns=load_functions('run_phase5_chemical_world_model.py',['module_D'],
                          dict(EXPLORATORY_KINETICS=False))
        with self.assertRaisesRegex(RuntimeError,'opt-in'):ns['module_D']()

    def test_old_phase5_cache_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            file=Path(td)/'old.json';file.write_text(json.dumps({'module_A':{'G_kcal':123}}))
            ns=dict(RESULTS_PATH=file,RESULTS={},XA={},AUDIT_VERSION=AUDIT_VERSION,
                    json=json,_warn=lambda *a:None)
            load_functions('run_phase5_chemical_world_model.py',['_hydrate_from_results'],ns)
            ns['_hydrate_from_results']()
            self.assertEqual(ns['XA'],{})


if __name__ == '__main__':unittest.main()
