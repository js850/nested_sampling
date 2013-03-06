"""classes and routines for running nested sampling"""
import random
import numpy as np

class MonteCarloChain(object):
    """from an initial configuration do a monte carlo chain with niter iterations
    
    the step will be accepted subject only to a maximum energy criterion
    """
    def __init__(self, potential, x0, takestep, Emax):
        self.potential = potential
        self.x = x0
        self.takestep = takestep
        self.Emax = Emax
        self.energy = self.potential.getEnergy(self.x)
        
        self.nsteps = 0
        self.naccept = 0
    
    def step(self):
        xnew = self.x.copy()
        
        # displace configuration
        self.takestep(xnew)
        
        # get energy
        e = self.potential.getEnergy(xnew)
        
        accept = e < self.Emax
        
        if accept:
            self.x = xnew
            self.energy
            self.naccept += 1
        
        self.nsteps += 1

       
class Replica(object):
    def __init__(self, x, energy):
        self.x = x.copy()
        self.energy = float(energy)

class NestedSampling(object):
    """the main class for implementing nested sampling"""
    def __init__(self, system, nreplicas, mciter=100):
        self.system = system
        self.mciter=mciter
        
        self.max_energies = []
        
        self.setup_replicas(nreplicas)
    
    def create_replica(self):
        x = self.system.get_random_configuration()
        pot = self.system.get_potential()
        e = pot.getEnergy(x)
        return Replica(x, e)
    
    def sort_replicas(self):
        self.replicas.sort(key=lambda r:r.energy)
    
    def setup_replicas(self, nreplicas):
        self.replicas = [self.create_replica() for i in range(nreplicas)]
        self.sort_replicas()
        print "min replica energy", self.replicas[0].energy
        print "max replica energy", self.replicas[-1].energy
    
    def sample_replica(self, Emax):
        # choose a replica randomly
        rstart = random.choice(self.replicas)
        
        # do a monte carlo iteration
        mc = MonteCarloChain(self.system.get_potential(), rstart.x, self.system.get_takestep(), Emax)
        for i in range(self.mciter):
            mc.step()
        
        return Replica(mc.x, mc.energy)
            
    
    def one_iteration(self):
        # pull out the replica with the largest energy
        rmax = self.replicas.pop()
        
        # add store it for later analysis
        self.max_energies.append(rmax.energy)
        
        # get a new replica to replace it
        rnew = self.sample_replica(rmax.energy)
        
        self.replicas.append(rnew)
        self.sort_replicas()
        
if __name__ == "__main__":
    from bh_sampling import LJClusterNew
    natoms = 13
    nreplicas = 100
    system = LJClusterNew(natoms)
    
    ns = NestedSampling(system, nreplicas)
    for i in range(10):
        ns.one_iteration()
    print ns.max_energies
    