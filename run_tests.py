import unittest


from lj_run import run_nested_sampling_lj, LJClusterNew



class TestLJNS(unittest.TestCase):
    
    def setUp(self):
        natoms = 31
        self.system = LJClusterNew(natoms)
        self.db = self.system.create_database("lj31_small.db")
        self.minima = self.db.minima()[:10]
        self.nreplicas = 10
        self.mciter = 10
        self.maxiter=100

    def test_comp(self):    
        run_nested_sampling_lj(self.system, nreplicas=self.nreplicas, mciter=self.mciter, 
                            maxiter=self.maxiter, label="test", minima=self.minima, use_compiled=True, nproc=1)
    def test_ser(self):    
        run_nested_sampling_lj(self.system, nreplicas=self.nreplicas, mciter=self.mciter, 
                            maxiter=self.maxiter, label="test", minima=self.minima, use_compiled=False, nproc=1)
    def test_comp_par(self):    
        run_nested_sampling_lj(self.system, nreplicas=self.nreplicas, mciter=self.mciter, 
                            maxiter=self.maxiter, label="test", minima=self.minima, use_compiled=True, nproc=2)
    def test4_par(self):    
        run_nested_sampling_lj(self.system, nreplicas=self.nreplicas, mciter=self.mciter, 
                            iterscale=1, label="test", minima=self.minima, use_compiled=False, nproc=2)


if __name__ == "__main__":
    unittest.main()