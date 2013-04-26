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
        E_vec = 0.5 * self.kappa * dist_vec * dist_vec      
        E = 0.
        for e in E_vec:
            E += e
        print "getEnergy ", E
        return E
   
class HarParticle(BaseSystem):
    """ Defines a particle in an harmonic potential of n dimensions """
    
    def __init__(self, ndim, centre, kappa, Eground, Emax_init):
        self.ndim = ndim
        #kappa is an ndimensional array containing n spring constants
        self.kappa = np.asfarray(kappa)
        assert(self.kappa.size == ndim)
        # use numpy vectors for coordinates, centre must be a numpy vector
        self.centre = centre
        self.Eground = Eground
        self.Emax_init = Emax_init
            
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
        print "coords_norm", coords_norm
        if coords is not None and radius > coords_norm:
            return True
        else:
            return False
     
    def get_random_configuration_Emax(self, Emax):
        """make sure they're all inside the radius, get_config_test is not strictly necessary, consider removing it"""
        #radius is a scalar corresponding to the max distance from the centre
        print "Emax", float(Emax)
        radius = np.sqrt(2. * (float(Emax) - self.Eground))
        print "Emax - Eground =", (float(Emax) - self.Eground)
        print "radius =", radius
        coords = np.zeros(self.ndim)
        coords = self.vector_random_uniform_hypersphere() * radius
        print "sqrt.radius ", np.sqrt(radius)
        assert(self.get_config_tests_in(coords, radius) == True)
        return coords
    
    def get_random_configuration(self):
        """make sure they're all inside the radius, get_config_test is not strictly necessary, consider removing it"""
        #radius is a scalar corresponding to the max distance from the centre
        Emax = self.Emax_init
        radius = np.sqrt(2. * (Emax - self.Eground))
        coords = np.zeros(self.ndim)
        coords = self.vector_random_uniform_hypersphere() * radius
        assert(self.get_config_tests_in(coords, radius) == True)
        return coords

class HarRunner(object):
    
    def __init__(self, system):
        self.system = system
        self.pot = system.get_potential()
    
    def __call__(self, x0, mciter, stepsize, Emax, seed=None):
        return self.run(Emax)  
            
    def run(self, Emax):        
        self.Emax = Emax
        self.mciter = 1
        self.nsteps = 1
        self.naccept = 0.7
        self.x = self.system.get_random_configuration_Emax(self.Emax)
        self.energy = self.pot.getEnergy(self.x)
        return self

        
        