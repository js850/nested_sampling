import numpy as np
from nested_sampling.utils.rotations import vec_random_ndim

def vector_random_uniform_hypersphere(k):
    """return a vector sampled uniformly in a hypersphere of dimension k"""
    u = vec_random_ndim(k)
    #draw the magnitude of the vector from a power law density:
    #draws samples in [0, 1] from a power distribution with positive exponent k - 1.
    p = np.random.power(k)
    return p * u



class Harmonic(object):
    def __init__(self, ndim):
        self.ndim = ndim
    
    def getEnergy(self, x):
        assert len(x) == self.ndim
        return 0.5 * x.dot(x)
    
    def get_random_configuration(self, radius=10.):
        """ return a random vector sampled uniformly from within a hypersphere of dimensions self.ndim"""
        x = vector_random_uniform_hypersphere(self.ndim) * radius
        return x

