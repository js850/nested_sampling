import multiprocessing as mp

import numpy as np

from nested_sampling import Result

def random_displace(x, stepsize):
    x += np.random.uniform(low=-stepsize, high=stepsize, size=x.shape)
    return x

class MonteCarloChain(object):
    """Class for doing a Monte Carlo chain
    

    Parameters
    -----------
    potential : pele Potentials
        attribute of system with member function getEnergy (in essence a
        particular potential energy function)
    x : array
        are the coordinates
    takestep : callable takestep object
        take a random montecarlo step, imported from pele: takestep(x) makes
        a move from x
    Emax : float
        energy upperbound
    energy : pele LJCluster
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
    def __init__(self, potential, takestep=random_displace, accept_tests=None, events=None):
        self.potential = potential
        self.takestep = takestep
        self.accept_tests = accept_tests
        self.events = events
            
    def __call__(self, x0, mciter, stepsize, Emax, energy, seed=None):
        return self.run(x0, mciter, stepsize, Emax, energy, seed=seed)
    
    def run(self, x0, mciter, stepsize, Emax, energy, seed=None):        

        # make sure we're not starting with a bad configuration
        assert energy <= Emax
        if self.accept_tests is not None:
            if not self.accept_tests(self.x, 0.):
                raise Exception("initial configuration for monte carlo chain failed configuration test")
        
        naccept = 0
        x = x0.copy()
        for i in xrange(mciter):
            xnew = x.copy()
            
            # displace configuration
            self.takestep(xnew, stepsize)
            
            # get energy
            e = self.potential.getEnergy(xnew)
            
            accept = e < Emax
            
            # perform configuration tests if they exist
            if accept and self.accept_tests is not None:
                accept = self.accept_tests(xnew, e)
    
            if accept:
                x = xnew
                energy = e
                naccept += 1
            
            # process callback functions if they exist
            if self.events is not None:
                for event in self.events:
                    event(coords=x, x=x, energy=energy, accept=accept)
    
        
        res = Result()
        res.x = x
        res.Emax = Emax
        res.energy = energy
        res.nsteps = mciter
        res.naccept = naccept
        return res


        
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


class MCRunner(mp.Process):
    def __init__(self, conn, mc_runner):
        mp.Process.__init__(self)
        self.conn = conn
        self.mc_runner = mc_runner

    def do_MC(self, x0, mciter, stepsize, Emax, energy, seed):
        return self.mc_runner(x0, mciter, stepsize, Emax, energy, seed) 
     
    def run(self):
        while 1:
            message = self.conn.recv()
            #print "message", message
            if message == "kill":
                #print "terminating", self.name
                return
            elif message[0] == "do mc":
                #print "received message: calculating gradient"
#                print "recieved message", message[0]
                x0, mciter, stepsize, Emax, energy, seed = message[1:]
                res = self.do_MC(x0, mciter, stepsize, Emax, energy, seed)
#                print "sending results back"
                self.conn.send(res)
            else:
                print "error: unknown message: %s\n%s" % (self.name, message)


