"""classes and routines for running nested sampling"""
import random
import numpy as np
import sys
import multiprocessing as mp
import copy
from itertools import izip

#this import fixes some bugs in how multiprocessing deals with exceptions
import pele.utils.fix_multiprocessing
from mc_walker import MCRunner


class Replica(object):
    """ Replica is simply a pair of types coordinates (x) and energy
    
    also attached is some additional information
    
    Parameters
    ----------
    x : array
    energy : float
    niter : int
        the number of MC iterations this structure has already been through
    from_random : bool
        if true, this replica started life as a completely random configuration 
    """
    def __init__(self, x, energy, niter=0, from_random=True):
        self.x = x.copy()
        self.energy = float(energy)
        self.niter = niter
        self.from_random = True

    def copy(self):
        return copy.deepcopy(self)


class NestedSampling(object):
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
                  stepsize=None, nproc=1, verbose=True,
                  max_stepsize=0.5):
        self.system = system
        self.mciter = mciter
        self.nproc = nproc
        self.verbose = verbose
        self.nreplicas = nreplicas
        self.mc_runner = mc_runner
        self.stepsize = stepsize
        self.max_stepsize = max_stepsize

        self.max_energies = []
        self.store_all_energies = True
        
        self.setup_replicas(nreplicas)
    
        self.iter_number = 0
        
        print "nreplicas", len(self.replicas)
        print "mciter", self.mciter
        sys.stdout.flush()

        #initialize the parallel workers to do the Monte Carlo Walk
        if self.nproc > 1:
            self.connlist = []
            self.workerlist = []
            for i in range(self.nproc):
                parent_conn, child_conn = mp.Pipe()
                worker = MCRunner(child_conn, self.mc_runner)
                worker.daemon = True
                self.connlist.append(parent_conn)
                self.workerlist.append(worker)
                worker.start()

    def _do_monte_carlo_chain_parallel(self, configs, Emax):
        # pass the workers the starting configurations for the MC walk
        for conn, r in izip(self.connlist, configs):
            seed = np.random.randint(0, sys.maxint)
            message = ("do mc", r.x, self.mciter, self.stepsize, Emax, r.energy, seed)
            conn.send(message)

        # receive the results back from the workers
        results = [conn.recv() for conn in self.connlist]
        
        # update the replicas
        for r, mc in izip(configs, results):
            r.x = mc.x
            r.energy = mc.energy
            r.niter += mc.nsteps
        configs

        # print some data
        if self.verbose:
            accrat = float(sum(mc.naccept for mc in results))/ sum(mc.nsteps for mc in results)
            print "step:", self.iter_number, "%accept", accrat, "Emax", Emax, "Emin", self.replicas[0].energy, \
                "stepsize", self.stepsize

        return configs, results

    def _do_monte_carlo_chain(self, r, Emax):
        rsave = r.copy()
        seed = np.random.randint(0, sys.maxint)
        
        # do the monte carlo walk
        result = self.mc_runner(r.x, self.mciter, self.stepsize, Emax, r.energy, seed)

        # update the replica
        r.x = result.x
        r.energy = result.energy
        r.niter += result.nsteps
    
        # print some data
        if self.verbose:
            # print some data
            dist = np.linalg.norm(result.x - rsave.x0)
            print "step:", self.iter_number, "%accept", float(result.naccept) / result.nsteps, \
                "Enew", result.energy, "Eold", rsave.energy, "Emax", Emax, "Emin", self.replicas[0].energy, \
                "stepsize", self.stepsize, "distance", dist #, "%reject_config", float(result.nreject_config) / result.nsteps, result.nsteps - result.naccept 

        return r, result

    def do_monte_carlo_chain(self, configs, Emax):
        """
        from an initial configuration do a monte carlo chain with niter iterations, 
        
        the step will be accepted subject only to a maximum energy criterion.
        Re-iterates for mciter steps.  After each step update stepize to meet
        target_ratio. 
        """
        assert len(configs) == self.nproc

        if self.nproc > 1:
            rnewlist, results = self._do_monte_carlo_chain_parallel(configs, Emax)

        else:
            rnew, result = self._do_monte_carlo_chain(configs[0], Emax)
            rnewlist = [rnew]
            results = [result]

        for result, r in izip(results, configs):
            if result.naccept == 0:
                sys.stderr.write("WARNING: zero steps accepted in the Monte Carlo chain %d\n")
                sys.stdout.write("WARNING: zero steps accepted in the Monte Carlo chain %d\n")
                print >> sys.stderr, "WARNING: step:", self.iter_number, "%accept", float(result.naccept) / result.nsteps, \
                    "Enew", result.energy, "Eold", r.energy, "Emax", Emax, "Emin", self.replicas[0].energy, \
                    "stepsize", self.stepsize #, "distance", dist removed because dist is undefined for multiple processors

        self.adjust_step_size(results)
        return rnewlist
    
    def create_replica(self):
        """
        creates a random configuration, evaluates its energy and creates the corresponding Replica object
        """
        x = self.system.get_random_configuration()
        e = self.system.getEnergy(x)
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
    
    def adjust_step_size(self, mc):
        """
        Adjust the stepsize to keep the acceptance ratio at a good value

        Notes
        -----
        when ratio naccept/nsteps < target_ratio then decrease step size by a
        correction factor f (new-step-size = old-step-size*f, where 0<f<1),
        else if naccept/nsteps > target_ratio, then increases the step size
        (new-step-size = old-step-size/f, where 0<f<1), although this is
        increased with an upper bound: max_stepsize.
        """
        if self.stepsize is None: return
        f = 0.8
        target_ratio = 0.5
        max_stepsize = self.max_stepsize # these should to be passed
        ratio = float(sum(m.naccept for m in mc)) / sum(m.nsteps for m in mc)  
        
        if ratio < target_ratio:
            # reduce step size
            self.stepsize *= f
        else:
            # increase stepsize
            self.stepsize /= f
        if self.stepsize > max_stepsize:
            self.stepsize = max_stepsize
        
    def pop_replicas(self):
        """
        remove the replicas with the largest energy and store them in the max_energies array 
        """
        # pull out the replicas with the largest energy
        for i in xrange(self.nproc):
            r = self.replicas.pop()
            if i == 0 or self.store_all_energies:
                self.max_energies.append(r.energy)
        
    def get_starting_configurations_from_replicas(self):
        """use existing replicas as starting configurations"""
        # choose a replica randomly
        assert len(self.replicas) == (self.nreplicas - self.nproc)
        rlist = random.sample(self.replicas, self.nproc)
        
        # make a copy of the replicas so we don't modify the old ones
        rlist = [r.copy() for r in rlist]
        return rlist

    def get_starting_configurations(self, Emax):
        """return nproc replicas to be used as starting configurations"""
        return self.get_starting_configurations_from_replicas()
    
    def add_new_replicas(self, rlist):
        """add new replicas the the list and keep the list sorted"""
        for r in rlist:
            self.replicas.append(r)
        self.sort_replicas()
    
    def _get_new_Emax(self):
        return self.replicas[-self.nproc].energy
    
    def one_iteration(self):
        """do one iteration of the nested sampling algorithm"""
        Emax = self._get_new_Emax()

        self.pop_replicas()
        
        rlist = self.get_starting_configurations(Emax)
                    
        rlist = self.do_monte_carlo_chain(rlist, Emax)

        self.add_new_replicas(rlist)
            
        self.iter_number += 1    
        
        if self.nreplicas != len(self.replicas):
            raise Exception("number of replicas (%d) is not correct (%d)" % (len(self.replicas), self.nreplicas))

    def finish(self):
        """do any clean up that needs to be done at the end of the simulation"""
        if self.nproc > 1:
            for conn, worker in izip(self.connlist, self.workerlist):
                conn.send("kill")
                worker.join()
                worker.terminate()
                worker.join()
        print "copying live replicas"
        for replica in reversed(self.replicas):
            self.max_energies.append(replica.energy)


if __name__ == "__main__":
    from lj_run import LJClusterNew
    from pele.takestep import RandomDisplacement
    natoms = 6
    nreplicas = 10
    mciter = 1000
    system = LJClusterNew(natoms)
    
    takestep = RandomDisplacement(stepsize=0.5)
    mcrunner = MonteCarloChain(system.get_potential(), takestep, 
                               system.get_config_tests())
    
    ns = NestedSampling(system, nreplicas, mcrunner, mciter=mciter)
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

    
