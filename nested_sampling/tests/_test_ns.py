import unittest
import numpy as np

from nested_sampling import NestedSampling, MonteCarloWalker, Harmonic, Replica

class TestNS(unittest.TestCase):
    def setUp(self):
        self.setUp1()

    def setUp1(self, nproc=1):
        self.ndim = 3
        self.harmonic = Harmonic(self.ndim)
        self.nreplicas = 10
        self.stepsize = 0.1
        self.nproc = nproc
        
        self.mc_runner = MonteCarloWalker(self.harmonic, mciter=40)

        replicas = []
        for i in xrange(self.nreplicas):
            x = self.harmonic.get_random_configuration()
            replicas.append(Replica(x, self.harmonic.get_energy(x)))
        self.ns = NestedSampling(replicas, self.mc_runner, 
                                 stepsize=0.1, nproc=nproc, verbose=False)
        
        self.Emax0 = self.ns.replicas[-1].energy
        
        self.niter = 100
        for i in xrange(self.niter):
            self.ns.one_iteration()
        self.Emax = self.ns.replicas[-1].energy
        self.Emin = self.ns.replicas[0].energy
    
    def test1(self):
        self.assert_(len(self.ns.replicas) == self.nreplicas)
        self.assert_(self.Emax < self.Emax0)
        self.assert_(self.Emin < self.Emax)
        self.assert_(self.Emin >= 0)
        self.assert_(self.ns.stepsize != self.stepsize)
        self.assertEqual(len(self.ns.max_energies), self.niter * self.nproc)

class testNSPar(TestNS):
    def setUp(self):
        self.setUp1(nproc=3)
    
    
if __name__ == "__main__":
    unittest.main()  