import unittest,importlib.util,pathlib,math
import numpy as np
def load(name,folder,file):
    s=importlib.util.spec_from_file_location(name,pathlib.Path(__file__).resolve().parents[1]/"experiments"/folder/file)
    a=importlib.util.module_from_spec(s);s.loader.exec_module(a);return a
a=load("e008","SL-E008","coherence.py");b=load("e009","SL-E009","optical.py")
class NewTests(unittest.TestCase):
    def test_bounds(self):
        E=np.linspace(860000,864000,301);g=np.ones(len(E))/(E[-1]-E[0])
        self.assertEqual(a.amplitudes(0,E[0]),(0,0))
        for L in (1,1000):
            a.run(E,g,L)
            D=np.array([a.delta(e,L)[0,1] for e in E])
            dm=7.5e-5;dist=1e6;K=dm*dist/(2*a.HB)
            actual=abs(np.trapz(g*D*np.exp(-1j*K/E),E))
            self.assertLessEqual(actual,a.entry_bound(E,g,L,dist,dm)+1e-14)
    def test_scales(self):
        for dm in (7.5e-5,2.5e-3):
            s=a.scales(a.DIST,dm);w=s["energy_full_width_eV"];K=dm*a.DIST/(2*a.HB)
            self.assertAlmostEqual(K*w/(a.E0**2-w*w/4),1)
            self.assertAlmostEqual(math.exp(-s["packet_separation_m"]**2/(8*s["packet_sigma_half_overlap_m"]**2)),.5)
    def test_band(self):
        x=np.array([1.,2.,3.]);y=2*x
        self.assertAlmostEqual(b.band(x,y,1.2,2.8),2.8**2-1.2**2)
        self.assertAlmostEqual(b.band(x,y,1.2,2)+b.band(x,y,2,2.8),b.band(x,y,1.2,2.8))
    def test_ppm(self):
        for M in (2,4):
            for q0,q1 in ((0,.3),(.01,.5),(.1,.100001),(.8,.99),(.2,1),(.1,.1)):
                z=b.ppm(M,q0,q1)
                self.assertAlmostEqual(z,b.brute(M,q0,q1),delta=1e-12)
                self.assertGreaterEqual(z,0);self.assertLessEqual(z,math.log2(M)+1e-12)
        self.assertEqual(b.ppm(16,0,1),4)
    def test_guard(self):
        z=b.rate(100,10,.001,0,4);y=b.rate(100,10,.001,.001,4)
        self.assertAlmostEqual(z["ppm"]/y["ppm"],2);self.assertAlmostEqual(z["ook"]/y["ook"],2)
        self.assertAlmostEqual(y["mean_fraction"],.125)
        self.assertEqual(b.rate(0,10,.001,1e-6,4)["ppm"],0)
