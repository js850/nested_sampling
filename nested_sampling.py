"""classes and routines for running nested sampling"""
import random
import numpy as np
import sys

class MonteCarloChain(object):
    """from an initial configuration do a monte carlo chain with niter iterations
    
    the step will be accepted subject only to a maximum energy criterion
    """
    def __init__(self, potential, x0, takestep, Emax, accept_tests = None):
        self.potential = potential
        self.x = x0
        self.takestep = takestep
        self.Emax = Emax
        self.energy = self.potential.getEnergy(self.x)
        self.accept_tests = accept_tests
        
        self.nsteps = 0
        self.naccept = 0
    
    def step(self):
        xnew = self.x.copy()
        
        # displace configuration
        self.takestep(xnew)
        
        # get energy
        e = self.potential.getEnergy(xnew)
        
        accept = e < self.Emax

        if accept and self.accept_tests is not None:
            for test in self.accept_tests:
                a = test(energy=e, coords=xnew)
                if not a:
                    accept = False
                    break

        if accept:
            self.x = xnew
            self.energy = e
            self.naccept += 1
        
        self.nsteps += 1

       
class Replica(object):
    def __init__(self, x, energy):
        self.x = x.copy()
        self.energy = float(energy)

class NestedSampling(object):
    """the main class for implementing nested sampling"""
    def __init__(self, system, nreplicas, takestep, mciter=100, accept_tests=None):
        self.system = system
        self.takestep = takestep
        self.mciter=mciter
        self.accept_tests = accept_tests
        
        self.max_energies = []
        
        self.setup_replicas(nreplicas)
    
        self.iter_number = 0
        
        print "nreplicas", len(self.replicas)
        print "mciter", self.mciter
        sys.stdout.flush()
    
    
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
    
    def adjust_step_size(self, mc):
        f = 0.8
        target_ratio = 0.7
        max_stepsize = 0.5
        ratio = float(mc.naccept) / mc.nsteps
        if ratio < target_ratio:
            # reduce step size
            self.takestep.stepsize *= f
        else:
            self.takestep.stepsize /= f
        if self.takestep.stepsize > max_stepsize:
            self.takestep.stepsize = max_stepsize
    
    def do_monte_carlo_chain(self, x0, Emax, energy=np.nan):
        mc = MonteCarloChain(self.system.get_potential(), x0, self.takestep, Emax, accept_tests=self.accept_tests)
        for i in xrange(self.mciter):
            mc.step()

        verbose = True
        if verbose:
            # print some data
            dist = np.linalg.norm(mc.x - x0)
            print "step:", self.iter_number, "%accept", float(mc.naccept) / mc.nsteps, \
                "Enew", mc.energy, "Eold", energy, "Emax", Emax, "Emin", self.replicas[0].energy, \
                "stepsize", self.takestep.stepsize, "distance", dist
        
        self.adjust_step_size(mc)
        return mc
    
    def sample_replica(self, Emax):
        # choose a replica randomly
        rstart = random.choice(self.replicas)
        
        # do a monte carlo iteration
        mc = self.do_monte_carlo_chain(rstart.x, Emax, rstart.energy)
        
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
        
        self.iter_number += 1
        
if __name__ == "__main__":
    from lj_run import LJClusterNew
    from pygmin.takestep import RandomDisplacement
    natoms = 13
    nreplicas = 1000
    mciter = 100
    system = LJClusterNew(natoms)
    
    
    ns = NestedSampling(system, nreplicas, RandomDisplacement(stepsize=0.5), mciter=mciter)
    for i in range(nreplicas * 300):
        ediff = ns.replicas[-1].energy - ns.replicas[0].energy
        if ediff < .01: break  
        ns.one_iteration()
    print ns.max_energies
    print "min replica energy", ns.replicas[0].energy
    print "max replica energy", ns.replicas[-1].energy
    
    with open("max_energies", "w") as fout:
        fout.write( "\n".join([ str(e) for e in ns.max_energies]) )
        fout.write("\n")
        
    
    if True:
        import pylab as pl
        e = np.array(ns.max_energies)
        pl.plot(e)
        pl.show()

    
