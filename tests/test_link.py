import importlib.util, pathlib, unittest
spec=importlib.util.spec_from_file_location("link",pathlib.Path(__file__).resolve().parents[1]/"experiments"/"SL-E003"/"link.py")
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
class Link(unittest.TestCase):
    def test_oracle(self):
        r=a.compare();self.assertEqual(r["cases"],112);self.assertLess(r["max_absolute_probability_error"],1e-9);self.assertLess(r["max_relative_offdiagonal_error"],1e-5)
    def test_geometry(self):
        self.assertEqual(a.coverage(0,.1),0);self.assertEqual(a.coverage(1e9,.1),1);self.assertEqual(a.coverage(1e9,.1,zrx=.9),0)
        self.assertAlmostEqual(a.coverage(2,.1)/a.coverage(1,.1),4)
    def test_shell(self):
        self.assertEqual(a.shell_fraction(0),0);self.assertEqual(a.shell_fraction(1),1);self.assertAlmostEqual(a.shell_fraction(1e-8)/1e-16,.5)
    def test_cross_section(self):
        for e in (True,False):
            for t in (0,.25):
                self.assertAlmostEqual(a.xs(T0=t,electron=e),a.xs_quad(T0=t,electron=e),delta=1e-8)
    def test_information(self):
        self.assertEqual(a.info_bound(1,0,0),0);self.assertLess(a.info_bound(1,10,.1),a.info_bound(1,0,.1))
        self.assertAlmostEqual(a.info_bound(1,0,.2)/a.info_bound(1,0,.1),4)
        self.assertEqual(a.link(0,1)["info_upper_bits_day"],0)
        self.assertAlmostEqual(a.link(1,2)["info_upper_bits_day"]/a.link(1,1)["info_upper_bits_day"],16)
