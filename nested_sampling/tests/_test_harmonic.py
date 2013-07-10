import unittest
import numpy as np

from nested_sampling.harmonic import Harmonic

class TestHarmonic(unittest.TestCase):
    def setUp(self):
        self.ndim = 3
        self.harmonic = Harmonic(3)
    
    def test_energy(self):
        x = np.zeros(self.ndim)
        e = self.harmonic.getEnergy(x)
        self.assertAlmostEqual(e, 0., 7)
    
    def test_e2(self):
        x = np.ones(self.ndim)
        x *= 1. / np.linalg.norm(x)
        e = self.harmonic.getEnergy(x)
        self.assertAlmostEqual(e, 0.5, 7)
    
    def test_random_config(self):
        r = 1.
        x = self.harmonic.get_random_configuration(r)
        xnorm = np.linalg.norm(x)
        self.assertLessEqual(xnorm, r)

if __name__ == "__main__":
    unittest.main()