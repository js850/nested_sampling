"""classes and functions related to a particle in a n dimensional harmonic potential"""

import numpy as np
import random
from scipy.special import gamma, gammaln
from scipy.misc import factorial
from pygmin.accept_tests import SphericalContainer
from pygmin.utils.rotations import vec_random_ndim

__all__ =["HParticle"]

class HParticle():
    """
    Defines a particle in an harmonic potential of n dimensions
    """
    
    def __init__(self, ndim, Eground, radius, centre):
     # use numpy vectors for coordinates, centre must be a numpy vector   
        self.ndim = ndim
        self.Eground = Eground
        self.radius = radius
        self.centre = np.asarray(centre)
            
    def get_potential(self):
        
        potential = Eground + 0.5 *     
        return 
    
    def vector_random_uniform_hypersphere(self):
    """return a vector sampled uniformly in a hypersphere of dimension ndim"""
    u = vec_random_ndim(self.ndim)
    #draw the magnitude of the vector from a power law density:
    #draws samples in [0, 1] from a power distribution with positive exponent k/2 - 1.
    p = np.random.power(0.5 * self.ndim)
    return p * u
    
    def get_config_tests(self, coords):
        coords_norm np.linalg.norm(coords)
        if self.radius < coords_norm:
            return False
        else:
            return coords
     
    def get_random_configuration(self):
        """make sure they're all inside the radius"""
        test = self.get_config_tests()[0]
        coords = np.zeros([1,self.ndim])
        coords[i,:] = vector_random_uniform_hypersphere(self.ndim) * self.radius
        assert(np.linalg.norm(coords[i,:]) <= self.radius)
        assert(test.accept(coords.flatten()))
        return coords.flatten()
    