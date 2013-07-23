import multiprocessing as mp

import numpy as np

from nested_sampling import Result

def random_displace(x, stepsize):
    x += np.random.uniform(low=-stepsize, high=stepsize, size=x.shape)
    return x

class MonteCarloWalker(object):
    """Class for doing a Monte Carlo chain walk
    

    Parameters
    -----------
    potential :
        attribute of system with member function get_energy (in essence a
        particular potential energy function)
    x : array
        are the coordinates
    takestep : callable takestep object
        take a random montecarlo step, imported from pele: takestep(x) makes
        a move from x
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
    def __init__(self, potential, takestep=random_displace, accept_test=None, 
                  events=None, mciter=100):
        self.potential = potential
        self.takestep = takestep
        self.accept_test = accept_test
        self.events = events
        self.mciter = mciter
            
    def __call__(self, x0, stepsize, Emax, energy, seed=None):
        return self.run(x0, stepsize, Emax, energy, seed=seed)
    
    def run(self, x0, stepsize, Emax, energy, seed=None):        

        # make sure we're not starting with a bad configuration
        assert energy <= Emax
        if self.accept_test is not None:
            if not self.accept_test(self.x, 0.):
                raise Exception("initial configuration for monte carlo chain failed configuration test")
        
        naccept = 0
        x = x0.copy()
        for i in xrange(self.mciter):
            xnew = x.copy()
            
            # displace configuration
            self.takestep(xnew, stepsize)
            
            # get energy
            e = self.potential.get_energy(xnew)
            
            accept = e < Emax
            
            # perform configuration tests if they exist
            if accept and self.accept_test is not None:
                accept = self.accept_test(xnew, e)
    
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
        res.nsteps = self.mciter
        res.naccept = naccept
        return res


        
class MCWalkerParallelWrapper(mp.Process):
    """define a worker to do the monte carlo chain in a separate process
    """
    def __init__(self, conn, mc_runner):
        mp.Process.__init__(self)
        self.conn = conn
        self.mc_runner = mc_runner

    def do_MC(self, x0, stepsize, Emax, energy, seed):
        return self.mc_runner(x0, stepsize, Emax, energy, seed) 
     
    def run(self):
        while 1:
            message = self.conn.recv()
            if message == "kill":
                return
            elif message[0] == "do mc":
                x0, stepsize, Emax, energy, seed = message[1:]
                res = self.do_MC(x0, stepsize, Emax, energy, seed)
                self.conn.send(res)
            else:
                print "error: unknown message: %s\n%s" % (self.name, message)


