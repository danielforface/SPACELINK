import importlib.util, unittest
from pathlib import Path
spec=importlib.util.spec_from_file_location("audit",Path(__file__).resolve().parents[1]/"experiments/SL-E002/audit.py")
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
class Audit(unittest.TestCase):
    def test_coefficient(self): self.assertLess(abs(a.coefficient()/a.VFACTOR-1),.01)
    def test_zchannel(self):
        z=a.zchannel(); self.assertAlmostEqual(z["capacity"],.37,delta=.01); self.assertAlmostEqual(z["capacity"],z["grid_capacity"],delta=1e-7);self.assertAlmostEqual(z["reported_rate"],.1,delta=.005)
    def test_vacuum(self):
        S=a.propagate([(100000,0,.5)]);self.assertEqual(a.contrast(S),0);self.assertLess(a.diagnostics(S),1e-9)
    def test_analytic(self):
        for L in (1,100,1000,100000):
            ls=[(L,19.3,74/183.84)]; expected=a.analytic2(*ls[0]); got=a.leave(a.propagate(ls,s13=0),0)
            self.assertAlmostEqual(got,expected,delta=1e-10)
            if expected>1e-10:self.assertLess(abs(got/expected-1),1e-6)
    def test_layers(self):
        for ls in a.sample_profiles():
            S=a.propagate(ls);self.assertLess(a.diagnostics(S),1e-9)
            for i in range(3):self.assertLessEqual(a.leave(S,i),a.bound(ls)+1e-12)
            self.assertLess(abs(a.contrast(S,(1/3,1/3,1/3))),1e-10)
            self.assertLessEqual(abs(a.contrast(S)),a.bound(ls)+1e-12)
