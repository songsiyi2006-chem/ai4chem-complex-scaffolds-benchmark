import unittest
from dft_crosscheck import pair_difference, HARTREE_EV
from compare_quantum_models import matching_gas_pairs


class DFTPairTests(unittest.TestCase):
    def setUp(self):
        self.cation=dict(molecule_id='Q01',geometry_sha256='fixed',method='pbe0',basis='def2-svp',environment='gas',
                         charge=1,multiplicity=1,status='converged',energy_hartree=-200.)
        self.radical=dict(self.cation,charge=0,multiplicity=2,energy_hartree=-200.2)

    def test_same_state_reference_and_units(self):
        self.assertAlmostEqual(pair_difference(self.cation,self.radical),-.2*HARTREE_EV)

    def test_reject_mismatched_model_environment_or_geometry(self):
        for key,value in [('basis','def2-tzvp'),('environment','MeCN'),('geometry_sha256','different'),('molecule_id','Q03')]:
            with self.subTest(key=key),self.assertRaises(ValueError):
                pair_difference(self.cation,dict(self.radical,**{key:value}))

    def test_reject_failed_scf_and_wrong_spin(self):
        for change in [{'status':'failed'},{'multiplicity':1},{'charge':-1}]:
            with self.subTest(change=change),self.assertRaises(ValueError):
                pair_difference(self.cation,dict(self.radical,**change))

    def test_cross_method_comparison_requires_matching_gas_geometries(self):
        dft=dict(molecule_id='Q01',environment='gas',geometry_sha256='fixed',basis='def2-svp',
                 delta_E_radical_minus_cation_eV=-5.,radical_spin_squared=.75)
        xtb=[dict(molecule_id='P01',environment='gas',geometry_sha256='fixed',gfn=g,delta_model_eV=-6.) for g in (1,2)]
        self.assertEqual(len(matching_gas_pairs([dft],xtb)),2)
        with self.assertRaises(ValueError):matching_gas_pairs([dict(dft,environment='dmf')],xtb)
        with self.assertRaises(ValueError):matching_gas_pairs([dict(dft,geometry_sha256='other')],xtb)
        with self.assertRaises(ValueError):matching_gas_pairs([dft],xtb[:1])


if __name__=='__main__':unittest.main()
