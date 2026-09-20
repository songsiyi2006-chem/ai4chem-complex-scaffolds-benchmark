import unittest
import numpy as np
from additional_public_analysis import mean_interval,contrast,interpolate_current,unique_time_medians,window_change

class PublicAnalysisTests(unittest.TestCase):
    def test_t_not_z_small_n(self):
        lo,hi=mean_interval(10,2,4)
        self.assertAlmostEqual(hi-10,3.182446305284263,places=8)
        self.assertGreater(hi-10,1.96)
    def test_independent_equal_variance_welch_df(self):
        r=contrast(1,2,2,2,n=4)
        self.assertAlmostEqual(r['df'],6)
        self.assertAlmostEqual(r['standard_error'],np.sqrt(2))
    def test_covariance_sensitivity(self):
        low=contrast(99,.1,98,.1,rho=1);high=contrast(99,.1,98,.1,rho=-1)
        self.assertEqual(low['standard_error'],0)
        self.assertAlmostEqual(high['standard_error'],.1)
        with self.assertRaises(ValueError):contrast(1,1,1,1,rho=1.1)
    def test_interpolation_and_no_extrapolation(self):
        a,b,width=interpolate_current([1,2,3],[10,20,30],15)
        self.assertAlmostEqual(a,1.5);self.assertAlmostEqual(b,1.5);self.assertEqual(width,10)
        with self.assertRaises(ValueError):interpolate_current([1,2],[10,20],30)
    def test_no_silent_curve_sort(self):
        with self.assertRaises(ValueError):interpolate_current([1,2,3],[10,30,20],15)
    def test_duplicate_timestamp_balancing(self):
        x,y=unique_time_medians([0,0,0,1,2],[1,3,5,9,11])
        np.testing.assert_array_equal(x,[0,1,2]);np.testing.assert_array_equal(y,[3,9,11])
    def test_window_overlap_rejected(self):
        with self.assertRaises(ValueError):window_change(np.arange(10),np.arange(10)+1,6,3)

if __name__=='__main__':unittest.main()
