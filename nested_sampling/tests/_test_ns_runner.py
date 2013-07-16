import unittest
from nested_sampling import NestedSampling, MonteCarloWalker, Harmonic, run_nested_sampling


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

        self.ns = NestedSampling(self.harmonic, self.nreplicas, self.mc_runner, 
                                 stepsize=0.1, nproc=nproc, verbose=False)


        self.etol = 0.01
        run_nested_sampling(self.ns, label="test", etol=self.etol)
        self.Emax = self.ns.replicas[-1].energy
        self.Emin = self.ns.replicas[0].energy

    
    def test1(self):
        self.assert_(len(self.ns.replicas) == self.nreplicas)
        self.assert_(self.ns.stepsize != self.stepsize)
        self.assertGreater(len(self.ns.max_energies), 0)
        self.assertLess(self.Emin, self.Emax)
        self.assertLess(self.Emax - self.Emax, self.etol)
        self.assertLess(self.Emax, 2. * self.etol)
        

if __name__ == "__main__":
    unittest.main()  