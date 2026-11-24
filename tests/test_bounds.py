import importlib.util
import math
from pathlib import Path
import unittest
p=Path(__file__).resolve().parents[1]/"experiments"/"SL-E001"/"bounds.py"
spec=importlib.util.spec_from_file_location("bounds",p)
b=importlib.util.module_from_spec(spec); spec.loader.exec_module(b)
class TestBounds(unittest.TestCase):
    def test_zero(self):
        self.assertEqual(b.phase(1,.5,0),0)
        self.assertEqual(b.transition_bound(0),0)
    def test_scale(self):
        self.assertAlmostEqual(b.phase(1,.5,1000),3.8e-11/b.HBARC)
        self.assertAlmostEqual(b.transition_bound(1e-5)/b.transition_bound(1e-6),100,places=6)
    def test_saturation(self):
        self.assertEqual(b.transition_bound(4*math.pi),1)
    def test_distance(self):
        self.assertAlmostEqual(b.integration_days(.01,au=2)/b.integration_days(.01),4)
    def test_invalid(self):
        with self.assertRaises(ValueError): b.integration_days(0)
        with self.assertRaises(ValueError): b.phase(-1,.5,1)
