import unittest
import numpy as np

from nested_sampling import MonteCarloWalker, Harmonic

class TestMCWalker(unittest.TestCase):
    def setUp(self):
        self.prepare()
        self.runwalker()
    def prepare(self):
        self.ndim = 6
        self.harmonic = Harmonic(self.ndim)
        self.x0 = np.zeros(self.ndim)
        self.energy = self.harmonic.get_energy(self.x0)
        
    
        self.mciter = 100
        self.mcwalker = MonteCarloWalker(self.harmonic, mciter=self.mciter)
        self.stepsize = 0.1
        self.Emax = 10.
        self.seed = None
    
    def runwalker(self):    
        self.res = self.mcwalker(self.x0, self.stepsize, 
                                 self.Emax, self.energy, seed=self.seed)
        
    def test1(self):
        self.assert_(np.all(self.x0 != self.res.x))
    def test2(self):
        self.assert_(self.energy != self.res.energy)
        self.assert_(self.res.energy < self.Emax)
    def test3(self):
        self.assert_(self.res.naccept > 0)
    def test4(self):
        self.assert_(self.res.nsteps == self.mciter)
    
    def test_large_stepsize(self):
        self.stepsize = 1e4
        self.runwalker()
        self.assert_(self.res.naccept == 0)
        self.assert_(self.energy == self.res.energy)
        self.assert_(np.all(self.x0 == self.res.x))

    
    
if __name__ == "__main__":
    unittest.main()  
