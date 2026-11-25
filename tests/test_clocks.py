import importlib.util,pathlib,unittest,math
import numpy as np
s=importlib.util.spec_from_file_location('clocks',pathlib.Path(__file__).resolve().parents[1]/'experiments'/'SL-E021'/'clocks.py')
a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
class Clocks(unittest.TestCase):
    def test_survival(self):
        self.assertEqual(a.survive(0),1)
        self.assertAlmostEqual(a.survive(a.HALF),.5)
        self.assertAlmostEqual(10*a.survive(0),20*a.survive(a.HALF))
        self.assertAlmostEqual(a.survive(2*a.HALF),.25)
    def test_filter(self):
        self.assertLess(a.decay_filter(1)['amplitude'],a.decay_filter(28)['amplitude'])
        self.assertAlmostEqual(a.decay_filter(1e10)['amplitude'],1)
        for p in (1,7,28):
            h=a.LAM/complex(a.LAM,2*math.pi/p)
            self.assertLess(abs(a.quad_filter(p,200000)-h),2e-6)
    def test_identifiability(self):
        self.assertLess(a.collision()['error'],1e-14)
        self.assertEqual(a.matrix(1)['nullity'],7)
        self.assertEqual(a.matrix(2)['nullity'],6)
        self.assertEqual(a.matrix(8)['rank'],8)
        self.assertGreater(a.matrix(8)['condition_nonzero'],1e4)
    def test_fisher(self):
        x=a.dual(100,a.HALF);y=a.dual(10000,a.HALF)
        self.assertAlmostEqual(x['sd_days']/y['sd_days'],10)
        self.assertAlmostEqual(a.pulse_age(100,a.HALF)['sd_days'],1/a.LAM/10)
        self.assertLess(a.dual(100,a.HALF,math.inf)['sd_days'],x['sd_days'])
    def test_pair(self):
        self.assertAlmostEqual(a.pairing(a.HALF,0)['true_accept'],.5)
        self.assertEqual(a.pairing(a.HALF,0)['accidental'],0)
