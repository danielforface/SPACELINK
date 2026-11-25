import importlib.util,pathlib,unittest,math
s=importlib.util.spec_from_file_location('scars',pathlib.Path(__file__).resolve().parents[1]/'experiments'/'SL-E022'/'scars.py')
a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
class Scars(unittest.TestCase):
    def test_recoil(self):
        for E in (.3,.42,.862,1.442,5,10):
            lo,hi=a.range_recoil(E)
            self.assertLessEqual(lo,hi)
            self.assertLess(abs(a.exact(E,-1)-hi),max(.01,.001*hi))
            self.assertEqual(a.fraction(E,hi+1,0),0)
            self.assertEqual(a.fraction(E,0,-1/3),1)
        self.assertLess(a.recoil(a.Q*a.Q)*1e6,1)
    def test_tag(self):
        for x in (.9,.99):
            v=a.tag(x,114e-9,0)
            self.assertAlmostEqual(-math.expm1(-v['window_s']/114e-9),x)
            self.assertEqual(v['accidental'],0)
    def test_harvest(self):
        self.assertAlmostEqual(a.geom(1)/(1/(4*math.pi*a.D*a.D)),1)
        r=a.harvest(365,100)
        self.assertLess(r['stored_normal'],r['stored_ideal'])
        self.assertAlmostEqual(r['MI_ideal']/r['arrivals_ideal'],.5)
        self.assertEqual(a.onoff(0),0)
    def test_gate(self):
        g=a.gate(28)
        self.assertAlmostEqual(g['inventory_after_off']-g['initial_inventory'],g['missed_emissions'])
        self.assertEqual(g['missed_emissions'],g['excess_after_resume'])
