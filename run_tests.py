import unittest


from lj_run import run_nested_sampling, LJClusterNew



class TestLJNS(unittest.TestCase):
    
    def setUp(self):
        natoms = 31
        self.system = LJClusterNew(natoms)
        self.db = self.system.create_database("lj31_small.db")
        self.minima = self.db.minima()[:10]
        self.nreplicas = 10
        self.mciter = 10
        self.iterscale = 1

    def test1(self):    
        run_nested_sampling(self.system, nreplicas=self.nreplicas, mciter=self.mciter, 
                            iterscale=1, label="test", minima=self.minima, use_compiled=True, nproc=1)
    def test2(self):    
        run_nested_sampling(self.system, nreplicas=self.nreplicas, mciter=self.mciter, 
                            iterscale=1, label="test", minima=self.minima, use_compiled=False, nproc=1)
    def test3(self):    
        run_nested_sampling(self.system, nreplicas=self.nreplicas, mciter=self.mciter, 
                            iterscale=1, label="test", minima=self.minima, use_compiled=True, nproc=2)
    def test4(self):    
        run_nested_sampling(self.system, nreplicas=self.nreplicas, mciter=self.mciter, 
                            iterscale=1, label="test", minima=self.minima, use_compiled=False, nproc=2)


if __name__ == "__main__":
    unittest.main()