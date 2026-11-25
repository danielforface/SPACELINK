import unittest,importlib.util,pathlib,math,random
import numpy as np
def load(name,folder,file):
    s=importlib.util.spec_from_file_location(name,pathlib.Path(__file__).resolve().parents[1]/"experiments"/folder/file)
    a=importlib.util.module_from_spec(s);s.loader.exec_module(a);return a
a=load("e11","SL-E011","count.py");b=load("e12","SL-E012","solar_diffraction.py");c=load("e13","SL-E013","ambient.py")
class NewTests(unittest.TestCase):
    def test_count(self):
        rng=random.Random(20261125)
        for _ in range(100):
            x=10**rng.uniform(-4,4);y=x+10**rng.uniform(-4,4);p=rng.uniform(.1,.9)
            I,ber,err=a.count_mi(x,y,p);q0,q1=a.e9.probs(y-x,x,1)
            self.assertLess(err,1e-9);self.assertGreaterEqual(I+1e-9,a.e9.o7.mutual(p,q0,q1))
            self.assertLessEqual(I,a.e9.o7.h2(p)+1e-9)
        self.assertEqual(a.count_mi(1,1,.5)[0],0)
        self.assertAlmostEqual(a.count_mi(0,.4,.5)[0],a.e9.o7.mutual(.5,0,-math.expm1(-.4)))
    def test_reference(self):
        for x,y in ((.01,.02),(1,3),(10,12),(0,10),(.001,.002)):
            self.assertAlmostEqual(a.count_mi(x,y,.25)[0],a.reference(x,y,.25),delta=1e-12)
    def test_mixture(self):
        den=np.array([-20.,2.]);w=np.array([.5,.5]);f=np.array([.5,.1])
        p,q,e=b.mixture(862000.,den,w,f);self.assertAlmostEqual(sum(p),1);self.assertAlmostEqual(sum(q),1);self.assertLess(e,.003)
    def test_diffraction(self):
        x=b.diffraction(.001);y=b.diffraction(.01)
        self.assertAlmostEqual(y["first_zero_angle_rad"]/x["first_zero_angle_rad"],.1)
        self.assertAlmostEqual(y["weak_born_average_bound"]/x["weak_born_average_bound"],10)
        self.assertAlmostEqual(x["resonant_scalar"]/x["kappa_L"]**2,1,places=9)
        z=np.linspace(-100,100,20001)
        self.assertAlmostEqual(float(np.trapz(np.sinc(z/math.pi)**2,z)),math.pi,delta=.02)
    def test_radio(self):
        x=c.radio(70,1,1e9,100);y=c.radio(70,1,2e9,100)
        self.assertAlmostEqual(x["rx_W"]/y["rx_W"],4)
        G=4*math.pi/c.LAM**2;Gr=4*math.pi*c.RX/c.LAM**2
        self.assertAlmostEqual(c.transfer(1,1e9)/(G*Gr*(c.LAM/(4*math.pi*1e9))**2),1)
        self.assertLessEqual(x["awgn_ref_bits_s"],x["rx_W"]/(c.K*100*math.log(2))*(1+1e-12))
    def test_thermal(self):
        x=c.thermal(.09,1,280,5,10);y=c.thermal(.09,1,360,5,10)
        self.assertAlmostEqual(x["radiator_m2"]/y["radiator_m2"],(360/280)**4)
        self.assertAlmostEqual(x["collector_m2"]*1361,x["heat_W"]+4)
