"""classes and routines for running nested sampling"""
import random
import numpy as np
import sys
import copy
from itertools import izip
import Pyro4
import Pyro4.util
import multiprocessing as mp
#this import fixes some bugs in how multiprocessing deals with exceptions
import nested_sampling.utils.fix_multiprocessing
from nested_sampling._mc_walker import MCWalkerParallelWrapper
try:
    import queue
except ImportError:
    import Queue as queue

class Replica(object):
    """object to represent the state of a system
    
    also attached is some additional information
    
    Parameters
    ----------
    x : array
        the structural coordinates
    energy : float
        the energy of the structure
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
        """return a complete copy of self"""
        return copy.deepcopy(self)

class Forwarditem(Replica):
    """wrapper to Replica for distributed parallel calculations"""
    
    def __init__(self, replica, Emax, stepsize, seed, itemId):
        super(Forwarditem,self).__init__(replica.x, replica.energy, replica.niter, replica.from_random) 
        
        self.Emax = Emax
        self.seed = seed
        self.stepsize = stepsize
        self.Id=itemId
    
    def __str__(self):
        return "<Workitem id=%s>" % str(self.Id)

class NestedSampling(object):
    """the main class for implementing nested sampling

    Parameters
    ----------
    replicas : list
        list of objects of type Replica
    mc_walker: callable
        class of type MonteCarloWalker
        performs a Monte Carlo walk to sample phase space.
        It should return an object with attributes x, energy, nsteps, naccept, etc.
    nproc : int
        number of processors to use for parallel nested sampling
    verbose : bool
        print status messages
    iprint : int
        if verbose is true, then status messages will be printed every iprint iterations
    cpfile: str
        checkpoint file name
    cpfreq: int
        checkpointing frequency in number of steps
    cpstart: bool
        start calculation from checkpoint binary file
    dispatcher_URI: str
        address (URI) of dispatcher (required for distributed parallelisation)
    serializer: str
        choice of serializer
    Attributes
    ----------
    nreplicas : integer
        number of replicas
    stepsize : float
        starting stepsize. It is then adapted to meet some target MC acceptance
        by default 0.5
    max_energies : list
        array of stored energies (at each step the highest energy configuration
        is stored and replaced by a valid configuration)
    store_all_energies: bool
        store the energy of all the nproc replicas replaced when running NS
        in parallel
    """
    def __init__(self, replicas, mc_walker, stepsize=0.1, nproc=1, verbose=True,
                  max_stepsize=0.5, iprint=1, cpfile=None, cpfreq=10000, cpstart = False, dispatcher_URI=None, serializer='pickle'):
        self.nproc = nproc
        self.verbose = verbose
        self.iprint = iprint
        self.replicas = replicas
        self.nreplicas = len(self.replicas)
        self.sort_replicas()
        self.mc_walker = mc_walker
        self.stepsize = stepsize
        self.max_stepsize = max_stepsize
        self.cpfreq = cpfreq
        self.cpfile = cpfile
        self.cpstart = cpstart
        self.max_energies = []
        self.store_all_energies = True
        
        self.iter_number = 0
        self.failed_mc_walks = 0
        self._mc_niter = 0 # total number of monte carlo iterations
        
        #pyro
        self.serializer= serializer
        self.dispatcher_URI = dispatcher_URI
                
        if self.verbose:
            print "nreplicas", len(self.replicas)
            sys.stdout.flush()

        
        if self.nproc > 1:
            if dispatcher_URI != None:
                self._set_up_serializer()
            else:
                self._set_up_multiproc_parallelization()
    
    #===========================================================================
    # multiprocessing functions for cpu locked parallelisation (on single node)
    #===========================================================================
    
    def _set_up_multiproc_parallelization(self):
        #initialize the parallel workers to do the Monte Carlo Walk
        self.connlist = []
        self.workerlist = []
        for _ in xrange(self.nproc):
            parent_conn, child_conn = mp.Pipe()
            worker = MCWalkerParallelWrapper(child_conn, self.mc_walker)
            worker.daemon = True
            self.connlist.append(parent_conn)
            self.workerlist.append(worker)
            worker.start()
            
    def _do_monte_carlo_chain_parallel_multiproc(self, configs, Emax):
        """run all the monte carlo walkers in parallel"""
        # pass the workers the starting configurations for the MC walk
        for conn, r in izip(self.connlist, configs):
            seed = np.random.randint(0, sys.maxint)
            message = ("do mc", r.x, self.stepsize, Emax, r.energy, seed)
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
        if self.verbose and self.iter_number % self.iprint == 0:
            accrat = float(sum(mc.naccept for mc in results))/ sum(mc.nsteps for mc in results)
            print "step:", self.iter_number, "%accept", accrat, "Emax", Emax, "Emin", self.replicas[0].energy, \
                "stepsize", self.stepsize

        return configs, results
    
    #===========================================================================
    # end of multiprocessing based functions
    #===========================================================================
    
    #===========================================================================
    # Pyro4 functions for distributed computing parallelisation
    #===========================================================================

    def _set_up_serializer(self):
        sys.excepthook = Pyro4.util.excepthook
        Pyro4.config.SERIALIZER = self.serializer
        Pyro4.config.SERIALIZERS_ACCEPTED.add(self.serializer)       
        self.dispatcher = Pyro4.Proxy(self.dispatcher_URI)
        
    def collectresults(self):
        results={}
        results_list = []
        
        while len(results) < self.nproc:
            try:
                item = self.dispatcher.getResult()
                #print("Got result: %s (from %s)" % (item, item.processedBy))
                results[item.Id] = item
            except queue.Empty:
                pass
                #print("Not all results available yet (got %d out of %d). Work queue size: %d" % \
                #        (len(results), self.nproc, self.dispatcher.workQueueSize()))
    
        if self.dispatcher.resultQueueSize() > 0:
            print("there's still stuff in the dispatcher result queue, that is odd...")
            print "queue size", self.dispatcher.resultQueueSize()
        
        for key in sorted(results.iterkeys()):
            results_list.append(results[key])
            
        return results_list
       
    def _do_monte_carlo_chain_parallel_distributed(self, configs, Emax):
        """run all the monte carlo walkers in parallel"""
        
        for i in xrange(len(configs)):
            seed = np.random.randint(0, sys.maxint)
            replica = configs[i]
            #crete item to be put in the queue
            item = Forwarditem(replica, Emax, self.stepsize, seed, i)
            #put item in the queue
            self.dispatcher.putWork(item)
        
        #collect results from queue    
        results = self.collectresults()
                
        # update the replicas
        for r, mc in izip(configs, results):
            r.x = mc.x
            r.energy = mc.energy
            r.niter += mc.nsteps
        
        # print some data
        if self.verbose and self.iter_number % self.iprint == 0:
            accrat = float(sum(mc.naccept for mc in results))/ sum(mc.nsteps for mc in results)
            print "step:", self.iter_number, "%accept", accrat, "Emax", Emax, "Emin", self.replicas[0].energy, \
                "stepsize", self.stepsize

        return configs, results

    #===========================================================================
    # end of Pyro4 based functions
    #===========================================================================

    def _do_monte_carlo_chain_parallel(self, configs, Emax):
        
        if self.dispatcher_URI != None:
            rnewlist, result = self._do_monte_carlo_chain_parallel_distributed(configs, Emax)
        else:
            rnewlist, result = self._do_monte_carlo_chain_parallel_multiproc(configs, Emax)
        return rnewlist, result

    def _do_monte_carlo_chain(self, r, Emax):
        """do the monte carlo walk"""
        assert self.nproc == 1
        rsave = r.copy()
        seed = np.random.randint(0, sys.maxint)
        
        # do the monte carlo walk
        result = self.mc_walker(r.x, self.stepsize, Emax, r.energy, seed)

        # update the replica
        r.x = result.x
        r.energy = result.energy
        r.niter += result.nsteps
    
        # print some data
        if self.verbose and self.iter_number % self.iprint == 0:
            # print some data
            dist = np.linalg.norm(result.x - rsave.x)
            print "step:", self.iter_number, "%accept", float(result.naccept) / result.nsteps, \
                "Enew", result.energy, "Eold", rsave.energy, "Emax", Emax, "Emin", self.replicas[0].energy, \
                "stepsize", self.stepsize, "distance", dist #, "%reject_config", float(result.nreject_config) / result.nsteps, result.nsteps - result.naccept 

        return r, result

    def do_monte_carlo_chain(self, configs, Emax):
        """
        from an initial configuration do a monte carlo walk 
        
        the steps will be accepted subject only to a maximum energy criterion.
        At the end of the walk, update stepize to meet the target_ratio. 
        """
        assert len(configs) == self.nproc

        if self.nproc > 1:
            rnewlist, results = self._do_monte_carlo_chain_parallel(configs, Emax)
            
        else:
            rnew, result = self._do_monte_carlo_chain(configs[0], Emax)
            rnewlist = [rnew]
            results = [result]

        self._mc_niter += sum((result.nsteps for result in results))

        for result, r in izip(results, configs):
            if result.naccept == 0:
                self.failed_mc_walks += 1
#                sys.stderr.write("WARNING: zero steps accepted in the Monte Carlo chain\n")
#                sys.stdout.write("WARNING: zero steps accepted in the Monte Carlo chain\n")
                sys.stdout.write("WARNING: step: %d accept %g Enew %g Eold %g Emax %g Emin %g stepsize %g\n" % 
                                 (self.iter_number, float(result.naccept) / result.nsteps,
                                  result.energy, r.energy, Emax, self.replicas[0].energy,
                                  self.stepsize))
#                if True:
#                    from pele.utils.xyz import write_xyz
#                    write_xyz(open("error.xyz", "w"), r.x, title="energy %g %d" % (r.energy, self.iter_number))
#                if True:
##                    print "fail stdout"
##                    sys.stdout.write("fail stderr\n")
#                    sys.stdout.write("testing configuration\n")
#                    tests = self.system.get_config_tests()
#                    for test in tests:
#                        if not test(coords=r.x):
#                            sys.stdout.write("    test failed\n")
#                            raise Exception("configuration failed configuration tests")
#                        else:
#                            sys.stdout.write("    test passed\n")
#                if True:
#                    print "testing energy of replica"
#                    testenergy = self.system.get_energy(r.x)
#                    assert np.abs(testenergy - r.energy) < 1e-5
                    
        self.adjust_step_size(results)
        return rnewlist
    
    def adjust_step_size(self, results):
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
        naccept = sum(m.naccept for m in results)
        nsteps = sum(m.nsteps for m in results)
        ratio = float(naccept) / nsteps  
        
        if ratio < target_ratio:
            # reduce step size
            self.stepsize *= f
        else:
            # increase stepsize
            self.stepsize /= f
        if self.stepsize > max_stepsize:
            self.stepsize = max_stepsize
    
    def sort_replicas(self):
        """
        sorts the replicas in decreasing order of energy
        """
        self.replicas.sort(key=lambda r:r.energy)
        
    def pop_replicas(self):
        """
        remove the replicas with the largest energies and store them in the max_energies array 
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
        self.starting_replicas = rlist # save in case the user wants to see it
        
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
        # get the new Emax
        Emax = self._get_new_Emax()

        # remove the npar replicas with the highest energy
        self.pop_replicas()
        
        # get starting configurations for the monte carlo walk
        rlist = self.get_starting_configurations(Emax)
        
        # do the monte carlo walk
        rlist = self.do_monte_carlo_chain(rlist, Emax)
        self.new_replicas = rlist # for optional later access

        # add the new replicas and keep the list sorted
        self.add_new_replicas(rlist)
            
        self.iter_number += 1    
        
        if self.nreplicas != len(self.replicas):
            raise Exception("number of replicas (%d) is not correct (%d)" % (len(self.replicas), self.nreplicas))
        
#        print "copying live replicas"
#        for replica in reversed(self.replicas):
#            self.max_energies.append(replica.energy)
        
                