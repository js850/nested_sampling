"""this is a simplified version of nested sampling without
any parallel components.  For bug testing"""

import sys
from itertools import izip
import numpy as np
from nested_sampling import Replica


class NestedSamplingSerial(object):
    """the main class for implementing nested sampling

    Parameters
    ----------
    system : pele System
        is the particular system of interest, say LJCluster
    nreplicas : integer
        number of replicas
    mc_runner : callable
        does Markov Chain walk to create a new configuration from an old configuration.
        It should return an object with attributes x, energy, nsteps, naccept, etc.
    mciter : integer
        number of steps in Markov Chain (sampling)
    accept_test : list of callables
        it's an array of pointers to functions. The dereferenced functions
        operate a set of tests on the energy/configuration.
    nproc : int
        number of processors to use for parallel nested sampling
    
    Attributes
    ----------
    max_energies : list
        array of stored energies (at each step the highest energy configuration
        is stored and replaced by a valid configuration)
    replicas : list
        list of objects of type Replica
        """
    def __init__(self, system, nreplicas, mc_runner, mciter=100, 
                 stepsize=None, verbose=True):
        self.system = system
        self.mciter = mciter
        self.verbose = verbose
        self.nreplicas = nreplicas
        self.mc_runner = mc_runner
        self.stepsize = stepsize

        self.max_energies = []
        
        self.setup_replicas(nreplicas)
    
        self.iter_number = 0
        
        print "nreplicas", len(self.replicas)
        print "mciter", self.mciter
        sys.stdout.flush()    
    
    def do_monte_carlo_chain(self, rold, Emax, **kwargs):
        """
        from an initial configuration do a monte carlo chain with niter iterations, 
        
        the step will be accepted subject only to a maximum energy criterion.
        Re-iterates for mciter steps.  After each step update stepize to meet
        target_ratio. 
        """
        
        x0, energy = rold.x.copy(), rold.energy
        x0save = x0.copy()
        seed = np.random.randint(0, sys.maxint)
        res = self.mc_runner(x0, self.mciter, self.stepsize, Emax, energy, seed)
        rnew = rold
        rnew.x = res.x
        rnew.energy = res.energy
        rnew.niter += res.nsteps
    
        if self.verbose:
            # print some data
            dist = np.linalg.norm(res.x - x0save)
            print "step:", self.iter_number, "%accept", float(res.naccept) / res.nsteps, \
                "Enew", res.energy, "Eold", energy, "Emax", Emax, "Emin", self.replicas[0].energy, \
                "stepsize", self.stepsize, "distance", dist #, "%reject_config", float(res.nreject_config) / res.nsteps, res.nsteps - res.naccept 
        
#        for res, r in izip(mclist, configs):
#            if res.naccept == 0:
#                sys.stderr.write("WARNING: zero steps accepted in the Monte Carlo chain %d\n")
#                sys.stdout.write("WARNING: zero steps accepted in the Monte Carlo chain %d\n")
#                print >> sys.stderr, "WARNING: step:", self.iter_number, "%accept", float(res.naccept) / res.nsteps, \
#                    "Enew", res.energy, "Eold", r.energy, "Emax", Emax, "Emin", self.replicas[0].energy, \
#                    "stepsize", self.stepsize #, "distance", dist removed because dist is undefined for multiple processors

#        self.adjust_step_size(mclist)
        return rnew
    
    def create_replica(self):
        """
        creates a random configuration, evaluates its energy and creates the corresponding Replica object
        """
        x = self.system.get_random_configuration()
        pot = self.system.get_potential()
        e = pot.getEnergy(x)
#        if self.verbose: print "pot=", e
        return Replica(x, e)
    
    def sort_replicas(self):
        """
        sorts the replicas in decreasing order of energy
        """
        self.replicas.sort(key=lambda r:r.energy)
    
    def setup_replicas(self, nreplicas):
        """
        creates nreplicas replicas and sorts them in decreasing order of energy
        """
        self.replicas = [self.create_replica() for i in range(nreplicas)]
        self.sort_replicas()
        print "min replica energy", self.replicas[0].energy
        print "max replica energy", self.replicas[-1].energy
        
    def pop_replica(self):
        """
        remove the replica with the largest energy (last item in the list) and store it in the max_energies array. 
        
        This is then used for integration.
        """
        # pull out the replica with the largest energy
        rmax = self.replicas.pop()
        # print "Epop:", rmax.energy, "n_replicas", len(self.replicas)
        
        # add store it for later analysis
        self.max_energies.append(rmax.energy)
        return rmax

    def get_starting_configuration_from_replicas(self):
        # choose a replica randomly
        i = np.random.randint(0, len(self.replicas))
        r = self.replicas[i].copy()
        return r

    def get_starting_configuration(self, Emax):
        return self.get_starting_configuration_from_replicas()
    
    def add_new_replica(self, rnew):
#        rnew = Replica(x, energy)
        self.replicas.append(rnew)
        return rnew
    
    def one_iteration(self):
        Emax = self.replicas[-1].energy
       
        rstart = self.get_starting_configuration(Emax)
        self.pop_replica()
        
        rnew = self.do_monte_carlo_chain(rstart, Emax)

        self.add_new_replica(rnew)
            
        self.iter_number += 1        
        self.sort_replicas()
        
        if self.nreplicas != len(self.replicas):
            raise Exception("number of replicas (%d) is not correct (%d)" % (len(self.replicas), self.nreplicas))

    def finish(self):
        pass