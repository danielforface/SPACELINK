import importlib.util,pathlib,unittest,math,hashlib
spec=importlib.util.spec_from_file_location("screen5",pathlib.Path(__file__).resolve().parents[1]/"experiments"/"SL-E005"/"screen.py")
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
class MassAmbient(unittest.TestCase):
    def test_mass(self):
        for ls in a.e3.old.sample_profiles():
            r=a.mass_metrics(ls)
            self.assertLessEqual(r["tv"],r["B"]+1e-12)
            self.assertLessEqual(r["kl_bits"],r["chi_bits"]+1e-18)
            self.assertLessEqual(r["chi_bits"],r["column_ceiling_bits"]+1e-18)
        self.assertEqual(a.mass_metrics([(100,0,.5)])["kl_bits"],0)
    def test_kl(self):
        self.assertAlmostEqual(a.kl_term(1e-9)/1e-18,.5)
        self.assertAlmostEqual(a.kl_term(.1),1.1*math.log(1.1)-.1)
    def test_shell(self):
        rows=[(.1,.5),(.2,.5)]
        self.assertEqual(a.shell_coverage(rows,1,0),0)
        self.assertEqual(a.shell_coverage(rows,1,1e10),1)
        self.assertEqual(a.shell_coverage(rows,1,1e10,zrx=.9),0)
    def test_parse(self):
        with self.assertRaises(ValueError):a.parse(b"bad","bad",13,10)
        data="".join(f"{i/1000} .01\n" for i in range(1,101)).encode()
        sha=hashlib.sha1(f"blob {len(data)}\0".encode()+data).hexdigest()
        rows,n=a.parse(data,sha,2,1);self.assertEqual(len(rows),100);self.assertAlmostEqual(n,1)
    def test_b16_density(self):
        data=''.join(f'{i*.0005:.4f} '+ ' '.join(['2']*11)+'\n' for i in range(1,1001)).encode()
        sha=hashlib.sha1(f'blob {len(data)}\0'.encode()+data).hexdigest()
        rows,n=a.parse(data,sha,12,7);self.assertAlmostEqual(n,1);self.assertEqual(len(rows),1000)
    def test_optics(self):
        f=lambda D,dt,dr:a.laser(1,D,dt,dr,1550e-9)
        self.assertAlmostEqual(f(1e11,.22,5)/f(2e11,.22,5),4)
        self.assertAlmostEqual(f(1e11,.44,5)/f(1e11,.22,5),4)
        self.assertAlmostEqual(a.mirror(1,1e11,5)/a.mirror(1,2e11,5),4)
