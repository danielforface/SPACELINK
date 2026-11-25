import unittest,importlib.util,pathlib,math
import numpy as np
def load(name,folder,file):
    s=importlib.util.spec_from_file_location(name,pathlib.Path(__file__).resolve().parents[1]/"experiments"/folder/file)
    a=importlib.util.module_from_spec(s);s.loader.exec_module(a);return a
a=load("solar6","SL-E006","solar.py");b=load("opt7","SL-E007","optical.py")
class SolarOptical(unittest.TestCase):
    def test_mixture(self):
        np.testing.assert_allclose(a.mixture(-20),a.W,atol=1e-12)
        for d in (0,1,2,2.2):
            p=a.mixture(d);self.assertAlmostEqual(sum(p),1);self.assertLess(abs((p-a.reduced(d))@a.W),.003)
    def test_trace(self):
        for L in (0,1,1000,100000):
            S=a.e3.oracle([(L,19.3,74/183.84)])
            r=a.trace_metrics(S,a.e5.P);self.assertGreaterEqual(r["coherent_trace_distance"]+1e-15,r["dephased_tv"])
        self.assertAlmostEqual(a.trace_metrics(np.eye(3),a.e5.P)["coherent_trace_distance"],0)
    def test_phase(self):
        E=np.linspace(860000,864000,10001);g=np.ones(len(E))/(E[-1]-E[0])
        for K in (0,1e8,1e10,1e11):
            actual=abs(np.trapz(g*np.exp(-1j*K/E),E))
            self.assertLessEqual(actual,a.envelope(E,g,K)+1e-6)
        self.assertAlmostEqual(a.envelope(E,g,1e13)/a.envelope(E,g,2e13),2)
    def test_airy(self):
        self.assertEqual(b.airy(0),1);self.assertLess(b.airy(3.83170597),1e-14)
        for x in (1,3,7):self.assertAlmostEqual(b.airy(x),b.airy(x,1024),places=10)
    def test_optics(self):
        s,g=b.budget("mirror",1,0,1e-7,1)
        s2,g2=b.budget("mirror",10,0,1e-7,1)
        self.assertAlmostEqual(s2/s,10);self.assertAlmostEqual((g2-100)/(g-100),10)
        _,g5=b.budget("mirror",1,0,1e-7,5);self.assertAlmostEqual((g5-100)/(g-100),25)
    def test_click(self):
        self.assertEqual(b.channel(0,100,1e-3)["bits_slot"],0)
        for p in (.1,.5,.9):
            self.assertAlmostEqual(b.mutual(p,.01,.7),b.joint(p,.01,.7))
        self.assertLess(b.channel(1000,10000,1e-4)["bits_slot"],b.channel(1000,0,1e-4)["bits_slot"])
        self.assertLessEqual(b.channel(10000,0,.1)["bits_slot"],1)
        q=.6;p=1/(q*(1+2**(b.h2(q)/q)))
        self.assertGreaterEqual(b.mutual(p,0,q)+1e-12,b.channel(-math.log1p(-q),0,1)["bits_slot"])
