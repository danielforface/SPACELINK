import unittest,importlib.util,pathlib,math
import numpy as np
def load(name,folder,file):
    s=importlib.util.spec_from_file_location(name,pathlib.Path(__file__).resolve().parents[1]/"experiments"/folder/file)
    a=importlib.util.module_from_spec(s);s.loader.exec_module(a);return a
a=load("e14","SL-E014","receiver.py");b=load("e15","SL-E015","scatter.py");c=load("e16","SL-E016","relay.py")
class NewTests(unittest.TestCase):
    def test_receiver(self):
        y,x=a.sample(10,100,12,.6,31);yy,xx=a.sample(10,100,12,.6,31)
        np.testing.assert_equal(y,yy);np.testing.assert_equal(x,xx);self.assertLessEqual(x.max(),31)
        m=a.calibrate(np.array([[0,10],[10,0]]),np.array([0,1]),31)
        np.testing.assert_equal(a.decode(np.array([[0,10],[10,0]]),m)["differential"],[0,1])
    def test_code(self):
        x=np.array([[10,0]]*3+[[0,1]]*2+[[0,10]]*3+[[1,0]]*2)
        for v in a.repeat_decode(x).values():np.testing.assert_equal(v,[1,0])
    def test_null(self):
        for v in a.run(0,.6,31)["errors"].values():self.assertLess(abs(v["mean"]-.5),.025)
    def test_form(self):
        self.assertEqual(float(b.ff(0)),1)
        self.assertAlmostEqual(b.integ(.862,184,0,2,point_only=True),1)
        for h in (1e-16,1e-8,.1):
            self.assertAlmostEqual(b.integ(.862,184,0,h,point_only=True)/(h-h*h/4),1)
            self.assertAlmostEqual(b.integ(.862,184,0,h,point_only=True,backward=True)/(h*h/4),1)
    def test_scatter(self):
        z=b.run(.862,184,74);self.assertLessEqual(z["form_factor_integral"],1)
        self.assertLess(z["convergence"],1e-10)
        self.assertAlmostEqual(b.point(1.724,184,74)/b.point(.862,184,74),4)
    def test_relay(self):
        z=c.relay(225e9,.5,1)
        self.assertLessEqual(z["passive_rx_W"],.5*z["first_rx_W"])
        self.assertAlmostEqual(c.optical(1e9,.22,5)/c.optical(2e9,.22,5),4)
        self.assertAlmostEqual(z["passive_ratio"],c.relay(225e9,.5,1)["passive_rx_W"]/z["direct_rx_W"])
        self.assertAlmostEqual(c.radio(1e9)["transmitter_W"][2],122.49)
