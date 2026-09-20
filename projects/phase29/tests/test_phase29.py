import json
import unittest
import numpy as np
from rdkit import Chem
from phase29 import ROOT, panel_mol, carbonyl_partner, proposed_product, canonical, simulate, kinetic_rhs, selectivity, group_split, conformal_radius, manufacturing_metrics, acquisition_score

class Phase29Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.panel=json.loads((ROOT/'configs/preregistered_panel.json').read_text(encoding='utf8'))
        cls.base=json.loads((ROOT/'configs/compute.json').read_text(encoding='utf8'))['kinetics']
    def test_panel_unlabelled_and_size(self):
        self.assertEqual(len(self.panel['substrates']),12)
        self.assertTrue(all(p['yield'] is None for p in self.panel['substrates']))
    def test_unique_heavy_atom_maps(self):
        for p in self.panel['substrates']:
            mol=panel_mol(p)
            maps=[a.GetAtomMapNum() for m in [mol,carbonyl_partner()] for a in m.GetAtoms()]
            self.assertEqual(len(maps),len(set(maps)))
            for site in [2,4,6]:self.assertEqual(sorted(maps),sorted(a.GetAtomMapNum() for a in proposed_product(mol,site).GetAtoms()))
    def test_element_formula_and_charge_conserved(self):
        from rdkit.Chem.rdMolDescriptors import CalcMolFormula
        for p in self.panel['substrates']:
            mol=panel_mol(p);reactants=Chem.CombineMols(mol,carbonyl_partner())
            for site in [2,4,6]:
                prod=proposed_product(mol,site)
                self.assertEqual(CalcMolFormula(prod),CalcMolFormula(reactants))
                self.assertEqual(Chem.GetFormalCharge(prod),0)
    def test_site_graph_distance(self):
        mol=panel_mol(self.panel['substrates'][1])
        for site,distance in [(2,1),(4,3),(6,1)]:
            p=proposed_product(mol,site);idx={a.GetAtomMapNum():a.GetIdx() for a in p.GetAtoms()}
            self.assertIsNotNone(p.GetBondBetweenAtoms(idx[site],idx[100]))
            self.assertEqual(len(Chem.GetShortestPath(p,idx[1],idx[site]))-1,distance)
    def test_parent_symmetry_not_counted_as_distinct(self):
        mol=panel_mol(self.panel['substrates'][0])
        self.assertEqual(canonical(proposed_product(mol,2)),canonical(proposed_product(mol,6)))
        self.assertNotEqual(canonical(proposed_product(mol,2)),canonical(proposed_product(mol,4)))
    def test_occupied_C2_rejected(self):
        mol=Chem.MolFromSmiles('[n:1]1[c:2](-c2ccccc2)[cH:3][cH:4][cH:5][cH:6]1')
        with self.assertRaises(ValueError):proposed_product(mol,2)
    def test_family_holdout_not_leaked(self):
        development={p['family'] for p in self.panel['substrates'] if p['split']=='development'}
        prospective={p['family'] for p in self.panel['substrates'] if p['split']=='prospective_holdout'}
        self.assertFalse(development & prospective)
    def test_ode_conserves_substrate_equivalents(self):
        t,y=simulate(self.base)
        np.testing.assert_allclose(y[:,:5].sum(1),1,atol=1e-8)
        self.assertGreaterEqual(y[:,:5].min(),-1e-8)
    def test_relay_removed_no_product(self):
        _,y=simulate(dict(self.base,relay_setpoint=0))
        np.testing.assert_allclose(y[-1,1:4],0,atol=1e-12)
    def test_pulse_charge_and_mass_balance(self):
        grid=np.arange(.005,9,.01)
        self.assertAlmostEqual(np.sum(np.where(grid%1<.5,2.,0.))*.01,9)
        _,y=simulate(dict(self.base,pulse=True))
        np.testing.assert_allclose(y[:,:5].sum(1),1,atol=1e-7)
    def test_zero_denominator_abstains(self):
        self.assertIsNone(selectivity([1,0,0,0,0,0,0])['s2_C2_over_C2_plus_C4'])
    def test_duplicate_families_stay_together(self):
        groups=np.repeat(np.arange(40),3);splits=group_split(groups)
        sets=[set(groups[s]) for s in splits]
        self.assertFalse(sets[0]&sets[1] or sets[0]&sets[2] or sets[1]&sets[2])
    def test_finite_sample_conformal_quantile(self):
        self.assertEqual(conformal_radius(np.arange(1,10),.1),9)
        self.assertTrue(np.isinf(conformal_radius([.1,.2,.3,.4],.1)))
        with self.assertRaises(ValueError):conformal_radius([])
    def test_manufacturing_units(self):
        m=manufacturing_metrics(1,1,1,1,.5,100,100,50,1)
        self.assertAlmostEqual(m['product_g'],50)
        self.assertAlmostEqual(m['energy_kWh_per_kg'],.02)
        self.assertAlmostEqual(m['PMI_without_water'],2)
        self.assertAlmostEqual(m['PMI_including_water'],3)
        self.assertAlmostEqual(m['STY_g_L_h'],50)
    def test_zero_yield_not_reported_as_infinite_advantage(self):
        with self.assertRaises(ValueError):manufacturing_metrics(1,1,1,1,0,100,100,0,1)
    def test_acquisition_rejects_impurity_and_values_yield(self):
        self.assertIsNone(acquisition_score(.9,.7,.5,1,.1))
        self.assertGreater(acquisition_score(.9,.7,.1,1,.1),acquisition_score(.9,.1,.1,1,.1))
        self.assertGreater(acquisition_score(.9,.7,.1,1,.1),acquisition_score(.9,.7,.1,3,.1))

if __name__=='__main__':unittest.main()
