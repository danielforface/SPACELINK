import importlib.util,pathlib,unittest,math
import numpy as np
s=importlib.util.spec_from_file_location("natural",pathlib.Path(__file__).resolve().parents[1]/"experiments"/"SL-E019"/"natural.py")
a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
class Natural(unittest.TestCase):
    def test_null_and_conservation(self):
        E=np.array([a.E0]);A=np.eye(3)[None,:,:]
        z=a.measure(E,A,0)
        self.assertLess(abs(z[0][0]),1e-24)
        self.assertEqual(a.bound(0),0)
        z=a.measure(E,A,1)
        self.assertEqual(z[2][0],0)
        self.assertLess(z[-1],1e-20)
        self.assertGreaterEqual(z[3][0],0)
        self.assertLessEqual(z[3][0],a.bound(1))
    def test_natural_preparation(self):
        E=np.array([a.E0]);tab=np.array([[1,1.67,0,0]])
        A=a.vac(E,a.Q)[:,:,None]*a.earth(E,100000,5000,tab)
        z=a.measure(E,A,1)
        self.assertGreater(abs(z[2][0]),abs(z[1][0]))
        self.assertGreaterEqual(z[3][0],0)
        self.assertLessEqual(z[3][0],a.bound(1))
        self.assertLess(a.two_state_oracle(),1e-14)
    def test_arbitrary_preparation(self):
        rng=np.random.default_rng(113)
        E=np.array([a.E0])
        for i in range(100):
            A,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
            for L in (.1,1,10):
                z=a.measure(E,A[None,:,:],L)
                self.assertGreaterEqual(z[3][0],-1e-19)
                self.assertLessEqual(z[3][0],a.bound(L)*(1+1e-8))
                self.assertLess(z[-1],1e-12)
    def test_sidebands(self):
        E=np.array([0.,1.,2.]);g=np.array([0.,1.,0.])
        r=a.sidebands(E,g)
        self.assertAlmostEqual(r[-1]["TV_bound"]/r[0]["TV_bound"],1000)
        self.assertAlmostEqual(r[0]["shift_eV"],4.135667696e-9)
