from nested_sampling import Replica

class Workitem(Replica):
    def __init__(self, x, energy, Emax, stepsize, seed, mc_runner, itemId, niter=0, from_random=True,):
        super(Workitem,self).__init__(x, energy, niter, from_random) 
        
        self.Emax = Emax
        self.seed = seed
        self.mc_runner = mc_runner
        self.stepsize = stepsize
        self.itemId=itemId
        self.processedBy=None
    
    def __str__(self):
        return "<Workitem id=%s>" % str(self.itemId)