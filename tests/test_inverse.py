import importlib.util,pathlib,unittest,numpy as np,math
s=importlib.util.spec_from_file_location("e17",pathlib.Path(__file__).resolve().parents[1]/"experiments"/"SL-E017"/"inverse.py")
a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
class Inverse(unittest.TestCase):
    def test_echo(self):
        for E in (859630,862000,868730):
            for M in (a.X,a.slab(E,a.VS,a.LS),np.eye(2)):
                for off in (0,a.Q,100000):
                    self.assertAlmostEqual(a.echo(E,M,off)[0],a.direct_phase(E,M,off),delta=1e-15)
        self.assertAlmostEqual(abs(a.slab(a.E0,a.VS,a.LS)[0,1])**2,1)
    def test_null(self):
        self.assertAlmostEqual(a.echo(a.E0,a.X,a.Q,0)[0],0,delta=1e-15)
        for t in (0,.123,10,1e6):
            D=np.diag(np.exp(1j*np.array([1,-1])*t))
            np.testing.assert_allclose(D@a.X@D,a.X,atol=1e-15)
    def test_maxswap(self):
        for v in (a.V,a.VS):
            for L in (1,10,1000,10000):
                self.assertLessEqual(abs(a.slab(a.E0,v,L)[0,1])**2,a.maxswap(a.E0,v)+1e-14)
