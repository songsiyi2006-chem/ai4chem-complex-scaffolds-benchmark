"""Analytic and adversarial software tests; no chemical reference labels."""
import json
import math
from pathlib import Path
import unittest
import numpy as np
from .analysis import *


class AnalysisTests(unittest.TestCase):
    def test_ensemble_degeneracy_and_energy_zero(self):
        g,w=ensemble_free_energy([0,0],[1,2]); self.assertAlmostEqual(g,-R_KCAL*298.15*math.log(3));np.testing.assert_allclose(w,[1/3,2/3])
        gp,wp=ensemble_free_energy([10000,10000],[1,2]); self.assertAlmostEqual(gp-g,10000);np.testing.assert_allclose(w,wp)

    def test_missing_energies_rejected(self):
        with self.assertRaises(ValueError): ensemble_free_energy([],[])
        with self.assertRaises(ValueError): ensemble_free_energy([float('nan')],[1])
        with self.assertRaises(ValueError): ensemble_free_energy([0],[0])

    def test_mirror_and_zero_conversion(self):
        a=signed_selectivity(9,1);b=signed_selectivity(1,9)
        self.assertAlmostEqual(a['ee'],-b['ee']);self.assertAlmostEqual(a['effective_ddG_kcal'],-b['effective_ddG_kcal'])
        self.assertIsNone(signed_selectivity(0,0)['ee']);self.assertEqual(signed_selectivity(1,1)['major'],'racemic')

    def test_sign_uncertainty_abstains(self):
        self.assertEqual(selectivity_interval(np.linspace(-.1,.2,1000))['major'],'selectivity_undetermined')
        self.assertEqual(selectivity_interval(np.linspace(.1,.2,1000))['major'],'R')

    def test_common_mode_uncertainty_cancels(self):
        samples=correlated_energy_draws([10,11],[[1,1],[1,1]],104729)
        np.testing.assert_allclose(samples[:,1]-samples[:,0],1,atol=1e-7)
        with self.assertRaises(ValueError): correlated_energy_draws([0,0],[[1,2],[2,1]],1)

    def test_reversible_network_matches_exact_solution(self):
        net=Network(['A','B'],[[1,0]],[[0,1]],[0,1],[18])
        t=np.linspace(0,10,101); c0=[.1,0]; sol=net.integrate(c0,t,[[1,1]])
        eq=.1*net.kf[0]/(net.kf[0]+net.kr[0])
        exact=eq*(1-np.exp(-(net.kf[0]+net.kr[0])*t))
        np.testing.assert_allclose(sol.y[1],exact,rtol=2e-6,atol=1e-10)
        self.assertAlmostEqual(net.kf[0]/net.kr[0],math.exp(-1/(R_KCAL*298.15)))

    def test_second_order_standard_state_and_balance(self):
        net=Network(['A','A2'],[[2,0]],[[0,1]],[0,-1],[18])
        self.assertAlmostEqual(net.flux([.2,0])[0]/net.flux([.1,0])[0],4)
        with self.assertRaises(ValueError): net.integrate([.1,0],[0,1],[[1,1]])
        self.assertGreater(standard_state_correction(),1.89);self.assertLess(standard_state_correction(),1.90)

    def test_ch_slow_or_disconnected_rejected(self):
        self.assertFalse(curtin_hammett_diagnostic([[-1,1],[1,-1]],[1,1])['supported'])
        self.assertTrue(curtin_hammett_diagnostic([[-1000,1000],[1000,-1000]],[1,1])['supported'])
        self.assertFalse(curtin_hammett_diagnostic(np.zeros((2,2)),[1,1])['supported'])

    def test_grand_energy_charge_and_potential(self):
        self.assertAlmostEqual(electrode_grand_energy(-100,101,100,-4),-96)
        p=dict(G_eV=-100,N_e=101,N_e_neutral=100,mu_e_eV=-4,actual_U_V_SHE=-1,total_cell_charge_e=-1,electrode_charge_partition_e=-.8)
        audit_potential_path([p,p,p],-1)
        with self.assertRaises(ValueError): audit_potential_path([p,p,p],-.8)
        with self.assertRaises(ValueError): audit_potential_path([dict(p,total_cell_charge_e=1)]*3,-1)

    def test_censoring_is_retained(self):
        result=competing_risk_incidence([1,1.5,2],['A','censored','B'])
        self.assertAlmostEqual(result[-1]['CIF']['A'],1/3)
        self.assertAlmostEqual(result[-1]['CIF']['B'],2/3)
        self.assertEqual(result[1]['censored'],1)
        all_censored=competing_risk_incidence([1,2],['censored','censored'])[-1]
        self.assertEqual(all_censored['survival'],1)

    def test_wilson_precision_does_not_imply_mechanism(self):
        self.assertLess(wilson_interval(192,385)['halfwidth'],.05)
        self.assertGreater(wilson_interval(50,100)['halfwidth'],.05)

    def test_path_rejects_missing_mixed_surface_wrong_endpoint(self):
        record=dict(opt_surface='A',freq_surface='A',irc_surface='A',imaginary_frequencies_cm1=[-500],imaginary_mode_is_reaction=True,irc_forward_endpoint='R',irc_reverse_endpoint='P',expected_endpoints=['R','P'],raw_log_sha256='0'*64,geometry_sha256='1'*64)
        self.assertTrue(validate_stationary_path(record))
        for change in [dict(freq_surface='B'),dict(irc_forward_endpoint='wrong'),dict(imaginary_frequencies_cm1=[]),dict(raw_log_sha256=None)]:
            with self.assertRaises(ValueError):validate_stationary_path(dict(record,**change))


class InputTests(unittest.TestCase):
    root=Path(__file__).resolve().parents[1]/'results_phase25_27'
    def test_no_enantiomer_or_family_leak(self):
        import csv
        with (self.root/'phase25'/'conditions.csv').open(encoding='utf8') as f: rows=list(csv.DictReader(f))
        self.assertEqual(len(rows),96)
        for pair in {r['pair_id'] for r in rows}:
            for sid in {r['substrate_id'] for r in rows}:
                self.assertEqual(len({r['split'] for r in rows if r['pair_id']==pair and r['substrate_id']==sid}),1)
        train=[r for r in rows if r['split']=='train']
        for r in rows:
            if r['split'] in ['test_catalyst','test_both']: self.assertNotIn(r['catalyst_family'],{x['catalyst_family'] for x in train})
            if r['split'] in ['test_substrate','test_both']: self.assertNotIn(r['substrate_family'],{x['substrate_family'] for x in train})

    def test_production_values_remain_unfilled(self):
        net=json.loads((self.root/'phase25'/'network.json').read_text())
        self.assertTrue(all(x['Gts_kcal_mol'] is None for x in net['reactions']))
        self.assertEqual(json.loads((self.root/'phase27'/'protocol.json').read_text())['trajectories_run'],0)

    def test_catalyst_connectivity_matches_primary_structure(self):
        from rdkit import Chem
        cats=json.loads((self.root/'phase25'/'catalysts.json').read_text())
        m=Chem.MolFromSmiles(next(x['mapped_connectivity_smiles'] for x in cats if x['pair_id']=='C03'))
        for a in m.GetAtoms(): a.SetAtomMapNum(0)
        ref=Chem.MolFromSmiles('OP1(=O)Oc2c(cc3ccccc3c2-c4c(O1)c(cc5ccccc45)-c6cc(cc(c6)C(F)(F)F)C(F)(F)F)-c7cc(cc(c7)C(F)(F)F)C(F)(F)F')
        self.assertEqual(Chem.MolToSmiles(m),Chem.MolToSmiles(ref))

    def test_actual_imine_geometry_preserves_EZ(self):
        from rdkit import Chem
        from rdkit.Chem import rdMolTransforms
        for sid in ['I01','I02','I03','I04']:
            for ez in ['E','Z']:
                for mol in Chem.SDMolSupplier(str(self.root/'phase25'/'starting_structures'/f'{sid}_{ez}.sdf'),removeHs=False):
                    self.assertIsNotNone(mol)
                    angle=abs(rdMolTransforms.GetDihedralDeg(mol.GetConformer(),9,1,2,3))
                    self.assertEqual(angle>90,ez=='E')

    def test_mirror_preserves_distances_and_inverts_axis(self):
        from rdkit import Chem
        from rdkit.Chem import rdMolTransforms
        audit=json.loads((self.root/'phase25'/'search_audit.json').read_text())
        for cid in ['C01','C02','C03','C04','C05','C06']:
            records=audit[cid]['records'];best=min((x for x in records if x['converged']),key=lambda x:x['MMFF_energy_kcal_mol'])['start']
            orig=list(Chem.SDMolSupplier(str(self.root/'phase25'/'starting_structures'/f'{cid}.sdf'),removeHs=False))[best]
            mirror=Chem.SDMolSupplier(str(self.root/'phase25'/'starting_structures'/f'{cid}_mirror.sdf'),removeHs=False)[0]
            p=orig.GetConformer().GetPositions();q=mirror.GetConformer().GetPositions()
            np.testing.assert_allclose(np.linalg.norm(p[:,None]-p,axis=2),np.linalg.norm(q[:,None]-q,axis=2),atol=1e-7)
            self.assertAlmostEqual(rdMolTransforms.GetDihedralDeg(orig.GetConformer(),4,13,14,15),-rdMolTransforms.GetDihedralDeg(mirror.GetConformer(),4,13,14,15),places=5)

    def test_component_maps_are_disjoint(self):
        from rdkit import Chem
        sets=[]
        for name in ['C04','I01_E','HE_dimethyl']:
            mol=Chem.SDMolSupplier(str(self.root/'phase25'/'starting_structures'/f'{name}.sdf'),removeHs=False)[0]
            maps=[a.GetAtomMapNum() for a in mol.GetAtoms()];self.assertEqual(len(maps),len(set(maps)));self.assertNotIn(0,maps);sets.append(set(maps))
        for i in range(3):
            for j in range(i):self.assertFalse(sets[i]&sets[j])


if __name__=='__main__':unittest.main()
