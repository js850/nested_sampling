import numpy as np

from nested_sampling.utils.rotations import vector_random_uniform_hypersphere

class Harmonic(object):
    def __init__(self, ndim):
        self.ndim = ndim
    
    def get_energy(self, x):
        assert len(x) == self.ndim
        return 0.5 * x.dot(x)
    
    def get_random_configuration(self, radius=10.):
        """ return a random vector sampled uniformly from within a hypersphere of dimensions self.ndim"""
        x = vector_random_uniform_hypersphere(self.ndim) * radius
        return x

