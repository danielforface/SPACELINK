import importlib.util, pathlib, unittest, random, math
spec=importlib.util.spec_from_file_location("screen",pathlib.Path(__file__).resolve().parents[1]/"experiments"/"SL-E004"/"screen.py")
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
class Tagged(unittest.TestCase):
    def test_perfect(self):
        for f in (1e-15,1e-9,.5,1):
            self.assertAlmostEqual(a.tagged(20,f,.01,1,0)/a.count(20,f,.01)*f,1)
    def test_no_tag(self):
        for t in (0,.1,.5,1):
            self.assertAlmostEqual(a.tagged(20,.001,.01,t,t)/a.count(20,.001,.01),1)
    def test_data_processing(self):
        rng=random.Random(20261124)
        for _ in range(1000):
            f=10**rng.uniform(-15,0);t=rng.random();l=rng.random()
            c=a.count(20,f,.01);p=a.tagged(20,f,.01,1,0);v=a.tagged(20,f,.01,t,l)
            self.assertGreaterEqual(v/c,1-1e-12);self.assertLessEqual(v/p,1+1e-12)
    def test_background(self):
        self.assertLess(a.tagged(20,.01,.01,1,0,100),a.tagged(20,.01,.01,1,0))
        self.assertEqual(a.tagged(20,0,.01,1,0),0)
    def test_scatter(self):
        self.assertEqual(a.cevns(L=0)["scatter_probability"],0)
        self.assertAlmostEqual(a.cevns(L=1000)["tau"]/a.cevns(L=1)["tau"],1000)
        self.assertAlmostEqual(a.cevns(E=1.724)["sigma_cm2"]/a.cevns()["sigma_cm2"],4)
    def test_spin_gravity(self):
        self.assertEqual(a.spin(0,1),0);self.assertAlmostEqual(a.spin(1e-15,10)/a.spin(1e-15,1),100)
        self.assertAlmostEqual(a.gravity(1000,1)/a.gravity(1000,2),2)
