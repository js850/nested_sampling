"""classes and routines for running nested sampling"""
import random
import numpy as np
import sys
import multiprocessing as mp
import copy
from itertools import izip

#this import fixes some bugs in how multiprocessing deals with exceptions
import pygmin.utils.fix_multiprocessing
from parallel_ns import MCRunner

class MonteCarloChain(object):
    """Class for doing a Monte Carlo chain
    

    Parameters
    -----------
    potential : pygmin Potentials
        attribute of system with member function getEnergy (in essence a
        particular potential energy function)
    x : array
        are the coordinates
    takestep : callable takestep object
        take a random montecarlo step, imported from pygmin: takestep(x) makes
        a move from x
    Emax : float
        energy upperbound
    energy : pygmin LJCluster
        energy of a given configuration
    accept_test : list of callables
        it's an array of pointers to functions. The dereferenced functions
        operate a set of tests on the energy/configuration.
    events : list fo callables
        it's an array of pointers to functions. This is general and not
        compulsury e.g. can use if you want to do something with the new
        configuration for the guy.

    Attributes
    ----------
    nsteps : integer
        tot. number of steps
    naccepts : integer
        tot. number of accepted configurations
    xnew : array
        new proposed configuration
    accept : boolean
        true or false if energy constraint is satisfied

    Notes
    -----
    defines parameters for a monte carlo chain from an initial configuration x
    with niter iterations. The step will be accepted subject only to a maximum
    energy criterion and geometric constraints on the configuration.  

    See Also
    --------
    NestedSampling
    """
    def __init__(self, potential, takestep, accept_tests=None, events=None):
        self.potential = potential
        self.takestep = takestep
        self.accept_tests = accept_tests
        self.events = events
            
    def __call__(self, x0, mciter, stepsize, Emax, energy, seed=None):
        return self.run(x0, mciter, stepsize, Emax, energy, seed=seed)
    
    def run(self, x0, mciter, stepsize, Emax, energy, seed=None):        
        self.x = x0
        self.Emax = Emax
        self.energy = self.potential.getEnergy(self.x)
        self.mciter = mciter
        self.takestep.stepsize = stepsize
        
        self.nsteps = 0
        self.naccept = 0
        self.nreject_config = 0
        
        if not self.test_configuration(self.x, 0.):
            print "ERROR: initial configuration for monte carlo chain failed configuration test"
            from pygmin.utils.xyz import write_xyz
            with open("error.xyz", "w") as fout:
                write_xyz(fout, self.x)
        
        for i in xrange(self.mciter):
            self.step()
#        print self.nsteps, self.naccept, self.nreject_config, self.takestep.stepsize, stepsize
        
        return self


        
    def test_configuration(self, x, e):
        if self.accept_tests is not None:
            for test in self.accept_tests:
                if not test(energy=e, coords=x):
                    self.nreject_config += 1
#                    print "rejecting config"    
                    return False
        return True

    def step(self):
        """
        Do one iteration in the Monte Carlo chain

        Notes
        -----
        copy current configuration x to xnew and then perform a trial move on
        xnew. Accept if the configuration of the new energy is less than Emax
        and if accept_tests are satisfied. If true update x to xnew. Finally,
        optionally applies all events functions to x.  
        """
        xnew = self.x.copy()
        
        # displace configuration
        self.takestep(xnew)
        
        # get energy
        e = self.potential.getEnergy(xnew)
        
#        print e, self.Emax
        accept = e < self.Emax

        if accept:
            accept = self.test_configuration(xnew, e)
#            print "accepting", self.naccept

        if accept:
            self.x = xnew
            self.energy = e
            self.naccept += 1
        
        if self.events is not None:
            for event in self.events:
                event(coords=self.x, x=self.x, energy=self.energy, accept=accept)

        self.nsteps += 1

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

#def mc_runner_wrapper(x_tuple):
#    mc_runner = x_tuple[0]
#    x0, mciter, stepsize, Emax, seed = x_tuple[1:]
#    return mc_runner(x0, mciter, stepsize, Emax, seed) 

class NestedSampling(object):
    """the main class for implementing nested sampling

    Parameters
    ----------
    system : pygmin System
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
                 stepsize=None, nproc=1, triv_paral=True, verbose=True):
        self.system = system
        self.mciter = mciter
        self.nproc = nproc
        self.verbose = verbose
        self.triv_paral = triv_paral
        self.nreplicas = nreplicas
        self.mc_runner = mc_runner
        self.stepsize = stepsize

        self.max_energies = []
        
        self.setup_replicas(nreplicas)
    
        self.iter_number = 0
        
        print "nreplicas", len(self.replicas)
        print "mciter", self.mciter
        sys.stdout.flush()

        #initialize the pool
        if self.nproc > 1:
#            self.pool = mp.Pool(processes=self.nproc)
            self.connlist = []
            self.workerlist = []
            for i in range(self.nproc):
                parent_conn, child_conn = mp.Pipe()
                worker = MCRunner(child_conn, self.mc_runner)
                self.connlist.append(parent_conn)
                self.workerlist.append(worker)
                worker.start()

    
    
    def do_monte_carlo_chain(self, configs, Emax, **kwargs):
        """
        from an initial configuration do a monte carlo chain with niter iterations, 
        
        the step will be accepted subject only to a maximum energy criterion.
        Re-iterates for mciter steps.  After each step update stepize to meet
        target_ratio. 
        """
        if self.nproc > 1:

            try:
                # pass the workers the starting configurations for the MC walk
                for conn, r in izip(self.connlist, configs):
                    seed = np.random.randint(0, sys.maxint)
                    message = ("do mc", r.x, self.mciter, self.stepsize, Emax, r.energy, seed)
                    conn.send(message)
                # recieve the results back from the workers
                mclist = [conn.recv() for conn in self.connlist]
            except:
                #this is bad, it shouldn't be done here
                print "exception caught during parallel MC iteration.  Terminating child processes"
                for worker in self.workerlist:
                    worker.terminate()
                    worker.join()
                print "done terminating child processes"
                raise
            
            for r, mc in izip(configs, mclist):
                r.x = mc.x
                r.energy = mc.energy
                r.niter += mc.nsteps
            rnewlist = configs

            if self.verbose:
                # print some data
                accrat = float(sum(mc.naccept for mc in mclist))/ sum(mc.nsteps for mc in mclist)
                print "step:", self.iter_number, "%accept", accrat, "Emax", Emax, "Emin", self.replicas[0].energy, \
                    "stepsize", self.stepsize


        else:
            rold = configs[0]
            x0, energy = rold.x.copy(), rold.energy
            x0save = x0.copy()
            attempts = 0
            success = 0
            while (success < 1) and (attempts < self.nreplicas):
                if attempts > 0:
                    print "warning: redoing monte carlo walk"
                seed = np.random.randint(0, sys.maxint)
                mc = self.mc_runner(x0, self.mciter, self.stepsize, Emax, energy, seed)
                attempts += 1
                success = mc.naccept
            rnew = rold
            rnew.x = mc.x
            rnew.energy = mc.energy
            #print "do_montecarlo_chain rold.energy", rold.energy
            #print "do_montecarlo_chain mc.energy", mc.energy
            rnew.niter += mc.nsteps
            rnewlist = [rnew]
            mclist = [mc]
        
            if self.verbose:
                # print some data
                dist = np.linalg.norm(mc.x - x0save)
                print "step:", self.iter_number, "%accept", float(mc.naccept) / mc.nsteps, \
                    "Enew", mc.energy, "Eold", energy, "Emax", Emax, "Emin", self.replicas[0].energy, \
                    "stepsize", self.stepsize, "distance", dist #, "%reject_config", float(mc.nreject_config) / mc.nsteps, mc.nsteps - mc.naccept 
        
        for mc, r in izip(mclist, configs):
            if mc.naccept == 0:
                sys.stderr.write("WARNING: zero steps accepted in the Monte Carlo chain %d\n")
                sys.stdout.write("WARNING: zero steps accepted in the Monte Carlo chain %d\n")
                print >> sys.stderr, "WARNING: step:", self.iter_number, "%accept", float(mc.naccept) / mc.nsteps, \
                    "Enew", mc.energy, "Eold", r.energy, "Emax", Emax, "Emin", self.replicas[0].energy, \
                    "stepsize", self.stepsize #, "distance", dist removed because dist is undefined for multiple processors

        self.adjust_step_size(mclist)
        return rnewlist
    
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
        target_ratio = 0.7
        max_stepsize = 0.5 # these should to be passed
        ratio = float(sum(m.naccept for m in mc))/ sum(m.nsteps for m in mc)  
        
        if ratio < target_ratio:
            # reduce step size
            self.stepsize *= f
        else:
            self.stepsize /= f
        if self.stepsize > max_stepsize:
            self.stepsize = max_stepsize
        
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
    
    def pop_replica_par(self, rtuple):
        """
        removes the other processes for parallel implementation
        
        remove nproc-1 replicas.  the one with the highest energy has already been
        removed.  
        
        non trivial
        -----------
        remove the nproc-1 highest energy replicas
        
        Trivial
        -------
        randomly remove nproc - 1 replicas stored in rtuple[1]
        for trivial parallelization we m
        
        """
        #remove other nproc-1 replica
        if self.triv_paral is True:
            indices = copy.copy(rtuple[1])
            # we need to remove one replica randomly from indices.  This will be the replica
            # that is duplicated before doing the monte carlo walk.    
            try:
                # if one of the indices is the highest energy one, then remove this one
                indices.remove(self.nreplicas-1) # this one has already been removed
            except ValueError:
                # indices is already randomly arranged, so removing the first one is fine
                indices = indices[1:]
            indices.sort(reverse=True)
            # sort the indices so the correct indices are removed
            # even though the list is being modified
            for i in indices:
                if i == self.nreplicas - 1: 
                    continue # this is already popped
                self.replicas.pop(i)  
        else:
            length = self.nproc-1
            for i in xrange(length):
                self.replicas.pop()
#        print len(self.replicas), len(indices)
#        if len(self.replicas) != self.nreplicas - self.nproc:
#            raise Exception("wrong number of replicas %d" % (len(self.replicas)))

    def get_starting_configuration_from_replicas(self):
        # choose a replica randomly
        rlist_int = random.sample(xrange(len(self.replicas)), self.nproc)
        configlist = [copy.deepcopy(self.replicas[i]) for i in rlist_int]
        return configlist,rlist_int

    def get_starting_configuration(self, Emax):
        return self.get_starting_configuration_from_replicas()
    
    def add_new_replica(self, rnew):
#        rnew = Replica(x, energy)
        self.replicas.append(rnew)
        return rnew
    
    def one_iteration(self):
#        rmax = self.pop_replica()
#        Emax = rmax.energy
        if self.nproc > 1 and self.triv_paral is False:
            Emax = self.replicas[-self.nproc].energy
        else:
            Emax = self.replicas[-1].energy
       
        rtuple = self.get_starting_configuration(Emax)
        configs = rtuple[0]

        self.pop_replica()
        if self.nproc > 1:
            self.pop_replica_par(rtuple)
        
#        if self.verbose: print "one_iteration Emax", Emax
        rlist = self.do_monte_carlo_chain(configs, Emax)

        for r in rlist:
            self.add_new_replica(r)
            
        self.iter_number += 1        
        self.sort_replicas()
        
        if self.nreplicas != len(self.replicas):
            raise Exception("number of replicas (%d) is not correct (%d)" % (len(self.replicas), self.nreplicas))

    def finish(self):
        if self.nproc > 1:
            for conn, worker in izip(self.connlist, self.workerlist):
                conn.send("kill")
                worker.join()
                worker.terminate()
                worker.join()
#            self.pool.close()
#            self.pool.join()


if __name__ == "__main__":
    from lj_run import LJClusterNew
    from pygmin.takestep import RandomDisplacement
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

    
