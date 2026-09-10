"""Regression checks for actual-output failure modes and complete-system inputs."""
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from .pilot_validation import frequency_block,qrrho_entropy,audit_job,OUT
from .ternary_analysis import rmsd


class PilotTests(unittest.TestCase):
    def test_alignment_does_not_superimpose_a_mirror(self):
        a=np.array([[0.,0,0],[1,0,0],[0,2,0],[0,0,3]])
        rot=np.array([[0.,-1,0],[1,0,0],[0,0,1]])
        self.assertLess(rmsd(a,a@rot+4),1e-12)
        mirror=a.copy();mirror[:,0]*=-1
        self.assertGreater(rmsd(a,mirror),.1)

    def test_repeated_modes_not_double_counted(self):
        block='projected vibrational frequencies (cm-1)\neigval : 0 0 0 0 0 0\neigval : -18.26 200 300\nother output\n'
        modes=frequency_block(block+block,3)
        self.assertEqual(len(modes),3);self.assertAlmostEqual(modes[0],-18.26)

    def test_wrong_mode_count_and_disagreement_rejected(self):
        with self.assertRaises(ValueError):frequency_block('projected vibrational frequencies\neigval : 1 2\n',3)
        a='projected vibrational frequencies\neigval : 0 0 0 0 0 0 1 2 3\nother\n'
        with self.assertRaises(ValueError):frequency_block(a+a.replace('1 2 3','1 2 4'),3)

    def test_negative_mode_rejected_despite_engine_zero_imaginaries(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)
            for f in ['input.xyz','xtbopt.xyz']: (p/f).write_text('3\nwater\nO 0 0 0\nH 1 0 0\nH 0 1 0\n')
            (p/'hessian.out').write_text('projected vibrational frequencies\neigval : 0 0 0 0 0 0 -18.26 200 300\n# imaginary freq. 0\n')
            row=dict(name='water',status='HESSIAN_COMPUTED_REVIEW_REQUIRED',method='GFN2-xTB',solvation='ALPB(toluene)',optimization_converged=True,hessian_exit_code=0)
            result=audit_job(p,row)
            self.assertFalse(result['accepted_local_minimum'])
            self.assertEqual(result['frequencies_below_minus1_cm1'],[-18.26])

    def test_thermochemistry_rejects_imaginary_modes(self):
        with self.assertRaises(ValueError):qrrho_entropy([-1,200],50)
        self.assertTrue(np.isfinite(qrrho_entropy([20,50,1000],100)))

    def test_ternary_atom_maps_and_both_faces(self):
        rows=json.loads((OUT/'phase25/ternary_starts/manifest.json').read_text())
        self.assertEqual(len(rows),8)
        for ez in ['E','Z']:
            for face in [-1,1]:self.assertEqual(sum(r['imine_E_Z']==ez and r['approach_face']==face for r in rows),2)
        for r in rows:
            self.assertEqual(r['atom_count'],117);self.assertEqual(len(set(r['atom_maps'])),117)
            self.assertEqual(r['components'].count(0),58);self.assertEqual(r['components'].count(1),28);self.assertEqual(r['components'].count(2),31)
            self.assertEqual(r['product_CIP'],'UNASSIGNED');self.assertIsNone(r['statistical_weight'])
            self.assertGreaterEqual(r['min_intercomponent_distance_over_covalent_radius_sum'],.72)


if __name__=='__main__':unittest.main()
