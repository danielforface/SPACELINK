import importlib.util,pathlib,unittest,math
import numpy as np
s=importlib.util.spec_from_file_location('electrical',pathlib.Path(__file__).resolve().parents[1]/'experiments'/'SL-E023'/'electrical.py')
a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
class Electrical(unittest.TestCase):
    def test_mass_and_cells(self):
        r=a.material()
        self.assertAlmostEqual(r['captures_100t_day'],a.RATE_GA*100000*a.FRAC)
        for s in (1e-6,.1,1):
            c=a.cells(s)
            self.assertAlmostEqual(c['aggregate_true_year'],r['captures_100t_day']*365)
            self.assertAlmostEqual(c['mean_false_hazard_per_cell_day']*c['cells'],.1*r['captures_100t_day'])
        self.assertLess(a.latch(100)['chi_bits_day'],a.latch(0)['chi_bits_day'])
    def test_coulomb(self):
        self.assertAlmostEqual(a.coulomb(10,4)['energy_eV']/a.coulomb(100,4)['energy_eV'],10)
        self.assertGreater(a.coulomb(10,4)['over_kBT'],a.coulomb(10,300)['over_kBT'])
    def test_gate(self):
        self.assertAlmostEqual(float(a.weight(np.array([0.]))[0]),2*math.pi*a.ALPHA*4)
        E=np.array([a.Q-100,a.Q,a.Q+100]);g=np.array([0,.01,0])
        self.assertEqual(a.interval(E,g,a.Q,a.Q,.1),0)
        den=a.interval(E,g,a.Q,E[-1],.1)
        self.assertGreater(den,0)
        self.assertAlmostEqual(den,a.interval(E,g,a.Q,a.Q+20,.1)+a.interval(E,g,a.Q+20,E[-1],.1))
        self.assertLess(abs(den-a.interval(E,g,a.Q,E[-1],.1,'mid')),1e-8)
    def test_direction(self):
        x=a.momentum(100);y=a.momentum(1000)
        self.assertGreater(x['within_patch'],y['within_patch'])
        self.assertAlmostEqual(x['within_patch']/y['within_patch'],100,delta=1e-5)
