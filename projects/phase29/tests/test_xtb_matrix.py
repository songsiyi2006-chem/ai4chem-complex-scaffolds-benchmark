import json
import unittest
from xtb_matrix import ROOT, HARTREE_EV, sha, parse_output, state_difference, ordering_reversals

class XTBMatrixTests(unittest.TestCase):
    def test_frozen_scope_144_single_points(self):
        c=json.loads((ROOT/'configs/xtb_matrix.json').read_text(encoding='utf8'))
        self.assertEqual(len(c['species'])*len(c['methods'])*len(c['environments'])*len(c['states']),144)
        self.assertEqual(sha(ROOT/'configs/preregistered_panel.json'),c['source_panel_sha256'])
    def test_original_three_geometries_exact_hash(self):
        c=json.loads((ROOT/'configs/xtb_matrix.json').read_text(encoding='utf8'))
        original=[s for s in c['species'] if s['geometry_policy']=='exact_original_preflight']
        self.assertEqual({s['id'] for s in original},{'P01','P04','P09'})
        for s in original:self.assertEqual(sha(ROOT/s['geometry_source']),s['geometry_sha256'])
    def test_warning_does_not_disappear_on_convergence(self):
        x=parse_output('convergence criteria satisfied\nTOTAL ENERGY -1.0000 Eh','normal termination of xtb\nIEEE_DIVIDE_BY_ZERO',0)
        self.assertEqual(x['status'],'success_with_warnings');self.assertEqual(x['warning_count'],1)
    def test_missing_scc_rejected_despite_printed_energy(self):
        x=parse_output('TOTAL ENERGY -1.0000 Eh','normal termination of xtb',0)
        self.assertEqual(x['status'],'failed_or_unverified');self.assertIsNone(x['energy_hartree'])
    def test_state_difference_sign_and_units(self):
        a=dict(molecule_id='P01',geometry_sha256='a',gfn=1,environment='gas',charge=1,uhf=0,status='success',energy_hartree=-10)
        b=dict(a,charge=0,uhf=1,energy_hartree=-10.1)
        self.assertAlmostEqual(state_difference(a,b),-.1*HARTREE_EV)
    def test_comparison_rejects_each_mismatch(self):
        a=dict(molecule_id='P01',geometry_sha256='a',gfn=1,environment='gas',charge=1,uhf=0,status='success',energy_hartree=-10)
        b=dict(a,charge=0,uhf=1,energy_hartree=-10.1)
        for k,v in [('molecule_id','P02'),('geometry_sha256','other'),('gfn',2),('environment','dmf'),('uhf',0),('status','failed')]:
            with self.subTest(field=k),self.assertRaises(ValueError):state_difference(a,dict(b,**{k:v}))
    def test_order_reversal_counts_pairs_not_molecules(self):
        x=ordering_reversals({'a':1,'b':2,'c':3},{'a':3,'b':2,'c':1})
        self.assertEqual(x['rank_reversals'],3);self.assertEqual(x['n_pairwise'],3)
    def test_ties_and_missing_are_not_forced_into_reversals(self):
        x=ordering_reversals({'a':1,'b':1.00000001,'c':3},{'a':2,'b':1})
        self.assertEqual(x['rank_reversals'],0);self.assertEqual(x['n_common'],2);self.assertEqual(len(x['ties_within_1e-6_eV']),1)

if __name__=='__main__':unittest.main()
