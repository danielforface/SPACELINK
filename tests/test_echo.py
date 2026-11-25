import importlib.util,pathlib,unittest,math
import numpy as np
s=importlib.util.spec_from_file_location("e18",pathlib.Path(__file__).resolve().parents[1]/"experiments"/"SL-E018"/"echo.py")
a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
class Echo(unittest.TestCase):
    def test_harmonics(self):
        n=a.choose()
        for M in a.modes(a.E0,n).values():
            self.assertAlmostEqual(a.echo(a.E0,M,a.imprint(a.E0))[0],a.direct(a.E0,M,a.imprint(a.E0)),delta=1e-14)
            self.assertLess(np.max(abs(M.conj().T@M-np.eye(3))),1e-10)
    def test_zero(self):
        for M in a.modes(a.E0,a.choose()).values():
            self.assertEqual(a.echo(a.E0,M,np.zeros((3,3))),(0.,0.))
        self.assertAlmostEqual(float(np.trace(a.imprint(a.E0)).real),0.,delta=1e-20)
    def test_geometry_delay(self):
        rows=[(.1,1.)]
        self.assertGreater(a.coverage(rows,1,1e6),a.coverage(rows,1,1e9))
        self.assertEqual(a.coverage(rows,1,1),1)
        self.assertEqual(a.tv(0,10),0)
        self.assertAlmostEqual(a.delay(1e6,2.5e-3)/a.delay(1e6,7.5e-5),100/3)
        self.assertLessEqual(a.tv(1e5,1e-12),2e-7)
    def test_independent_two_state(self):
        n=a.choose()
        self.assertLess(abs(a.twostate(a.E0,n)-abs(np.linalg.matrix_power(a.cell(a.E0),n)[1,0])**2),.02)
