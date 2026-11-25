import importlib.util,pathlib,unittest,math
import numpy as np
s=importlib.util.spec_from_file_location('records',pathlib.Path(__file__).resolve().parents[1]/'experiments'/'SL-E020'/'records.py')
a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
class Records(unittest.TestCase):
    def test_memory(self):
        self.assertEqual(a.effective(0),0)
        self.assertEqual(a.effective(10,math.inf),10)
        self.assertAlmostEqual(a.effective(1e6),a.HALF/a.LN2)
        self.assertAlmostEqual(a.rate(1000)/a.rate(1),1000)
        self.assertLess(a.effective(28),28)
    def test_information(self):
        for mu in (.01,1,100):
            self.assertAlmostEqual(a.binary_poisson_mi(0,mu),a.exclusive(mu),delta=1e-12)
        self.assertEqual(a.binary_poisson_mi(1,1),0)
        self.assertGreater(a.binary_poisson_mi(100,100.01),0)
        self.assertLess(a.binary_poisson_mi(100,110),a.chi(100,.1))
        self.assertAlmostEqual(a.marked(1e-8)['marked_kl'],200e-16/a.LN2,delta=1e-25)
    def test_gate(self):
        E=np.array([861000.,862000.,863000.]);g=np.array([0.,.001,0.])
        self.assertAlmostEqual(a.tail(E,g,862000),.5)
        for x in a.threshold(E,g)['rows']:
            self.assertLessEqual(abs(x['delta']),x['absolute_bound']+1e-15)
            self.assertLess(x['oracle_error'],1e-10)
    def test_resonance(self):
        r=a.resonance(.0005)
        self.assertLess(r['integral_rel_error'],1e-14)
        self.assertLess(r['average_upper_cm2'],r['peak_cm2'])
        self.assertGreater(r['width_eV'],0)
