import unittest

from lj_run import run_nested_sampling_lj, LJClusterNew

class TestLJNS(unittest.TestCase):
    
    def setUp(self):
        natoms = 6
        self.system = LJClusterNew(natoms)
        self.db = self.system.create_database("lj31_small.db")
        self.minima = self.db.minima()[:10]
        self.nreplicas = 10
        self.mciter = 100
        self.maxiter=10000
        
        self.kwargs = {"verbose":False}
        
        self.e1 = -12.712062
        self.e2 = -12.302928

    def assertNS(self, ns):
        self.assertLess(ns.replicas[0].energy, self.e2 + 0.5)

    def test_comp(self):
        ns = run_nested_sampling_lj(self.system, nreplicas=self.nreplicas, mciter=self.mciter, 
                            maxiter=self.maxiter, label="test", minima=self.minima, 
                            use_compiled=True, nproc=1, **self.kwargs)
        self.assertNS(ns)
        
    def test_ser(self):    
        ns = run_nested_sampling_lj(self.system, nreplicas=self.nreplicas, mciter=self.mciter, 
                            maxiter=self.maxiter, label="test", minima=self.minima, 
                            use_compiled=False, nproc=1, **self.kwargs)
        self.assertNS(ns)

    def test_comp_par(self):    
        ns = run_nested_sampling_lj(self.system, nreplicas=self.nreplicas, mciter=self.mciter, 
                            maxiter=self.maxiter, label="test", minima=self.minima, 
                            use_compiled=True, nproc=2, **self.kwargs)
        self.assertNS(ns)

    def test4_par(self):    
        ns = run_nested_sampling_lj(self.system, nreplicas=self.nreplicas, mciter=self.mciter, 
                            maxiter=self.maxiter, label="test", minima=self.minima, 
                            use_compiled=False, nproc=2, **self.kwargs)
        self.assertNS(ns)


if __name__ == "__main__":
    unittest.main()