import unittest

from _test_mc_walker import TestMCWalker
from nested_sampling import Harmonic
from nested_sampling.models.harmonic_nowalk import HarmonicSampler

class TestHarmonicSampler(TestMCWalker):
    def setUp(self):
        self.prepare()
        self.mcwalker = HarmonicSampler(self.harmonic, self.ndim)
        self.runwalker()
    
    def test_large_stepsize(self):
        # skipt this test
        pass

    
if __name__ == "__main__":
    unittest.main()