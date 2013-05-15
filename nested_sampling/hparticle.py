"""classes and functions related to a particle in a n dimensional harmonic potential"""

import numpy as np
import random
from scipy.special import gamma, gammaln
from scipy.misc import factorial
from pygmin.utils.rotations import vec_random_ndim
from pygmin.potentials import BasePotential
from pygmin.systems import BaseSystem
from nested_sampling import MonteCarloChain
import copy
from pygmin.optimize import Result

__all__ =["HarPotential","HarParticle","HarRunner"]

class HarPotential(BasePotential):
    """
    Defines an harmonic potential, returns energy depending on coordinates
    
    Parameters
    ----------
    centre: list of int
        centre of harmonic well
    kappa: list of float
        spring constants
    
    """
    def __init__(self, centre, kappa):
        self.centre = np.asfarray(centre)
        self.kappa = np.asfarray(kappa)
                
    def getEnergy(self, coords):
        coords = np.asfarray(coords)
        dist_vec = coords - self.centre
        return 0.5 * dist_vec.dot(self.kappa * dist_vec)
#        E_vec = 0.5 * self.kappa * dist_vec * dist_vec      
#        E = 0.
#        for e in E_vec:
#            E += e
##        print "getEnergy ", E
#        return E
   
class HarParticle(BaseSystem):
    """ Defines a particle in an harmonic potential of n dimensions 
    
    Notes
    -----
    with n degrees of freedom and spring constant k:
     
        Z = (2*pi*T/k)^(n/2)
        <E> = -n*T/2
        <E^2> = T^2*n*(n+2)/4
        Cv = n/2
    """
    
    def __init__(self, ndim, centre=None, kappa=None, Eground=0., Emax_init=10.):
        self.ndim = ndim
        #kappa is an ndimensional array containing n spring constants
        if kappa is None:
            self.kappa = np.ones(self.ndim)
        else:
            self.kappa = np.asfarray(kappa)
        assert(self.kappa.size == ndim)
        # use numpy vectors for coordinates, centre must be a numpy vector
        if centre is None:
            self.centre = np.zeros(self.ndim)
        else:
            self.centre = np.array(centre)
        self.Eground = Eground
        self.Emax_init = Emax_init
        
        print "Emax_init", self.Emax_init
            
    def get_potential(self):
        return HarPotential(self.centre, self.kappa) 
    
    def vector_random_uniform_hypersphere(self):
        """
        return a normalised vector sampled uniformly in a hypersphere of dimension ndim
        must then multiply each component by sqrt(norm_max)
        """
        u = vec_random_ndim(self.ndim)
        #draw the magnitude of the vector from a power law density:
        #draws samples in [0, 1] from a power distribution with positive exponent k/2 - 1.
        p = np.random.power(self.ndim)
        return p * u
    
    def get_config_tests(self):
        """
        needed to return true at all tests in run_nestedsampling
        """ 
        return True
    
    def get_config_tests_in(self, coords, radius):
        coords_norm = np.linalg.norm(coords-self.centre)
#        print "coords_norm", coords_norm
        if coords is not None and radius > coords_norm:
            return True
        else:
            return False
     
    def get_random_configuration_Emax(self, Emax):
        """make sure they're all inside the radius, get_config_test is not strictly necessary, consider removing it"""
        #radius is a scalar corresponding to the max distance from the centre
#        print "Emax", float(Emax)
        radius = np.sqrt(2. * (float(Emax) - self.Eground))
#        print "Emax - Eground =", (float(Emax) - self.Eground)
#        print "radius =", radius
#        coords = np.zeros(self.ndim)
        coords = self.vector_random_uniform_hypersphere() * radius
#        print "sqrt.radius ", np.sqrt(radius)
        assert(self.get_config_tests_in(coords, radius))
        return coords
    
    def get_random_configuration(self):
        """make sure they're all inside the radius, get_config_test is not strictly necessary, consider removing it"""
        #radius is a scalar corresponding to the max distance from the centre
        Emax = self.Emax_init
        return self.get_random_configuration_Emax(Emax)

class HarRunner(object):
    
    def __init__(self, system):
        self.system = system
        self.pot = system.get_potential()
    
    def __call__(self, x0, mciter, stepsize, Emax, energy, seed=None):
        return self.run(Emax)  
            
    def run(self, Emax):
        self.Emax = Emax
        res = Result()
        res.mciter = 100
        res.nsteps = 100
        res.naccept = 70
        res.x = self.system.get_random_configuration_Emax(self.Emax)
        res.energy = self.pot.getEnergy(res.x)
        return res

def test1():
    import matplotlib.pyplot as plt
    v = np.random.uniform(-4,4,size=10000)
    energies = v**2 / 2.
    plt.hist(energies, bins=50)
    plt.show()

def test():
    import matplotlib.pyplot as plt
    n = 1000
    hp = HarParticle(n, centre=[0.]*n, kappa=[1.]*n)
    pot = hp.get_potential()
    
    coords = np.array([hp.get_random_configuration_Emax(10.) for i in range(10000)])
    energies = [pot.getEnergy(x) for x in coords]
    print coords
#    plt.plot(vals, 'o')
    plt.hist(coords[:,0], bins=50)
    plt.show()
    plt.hist(energies, bins=50)
    plt.show()
    
if __name__ == "__main__":
    test()
#    test1()