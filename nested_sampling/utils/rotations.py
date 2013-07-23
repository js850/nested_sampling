"""
Functions related to random vectors

IMPORTED FROM PELE (POTENTIAL ENERY LANDSCAPE EXPLORATOR)
source can be found at https://github.com/pele-python/pele
"""
import numpy as np

def vec_random_ndim(n):
    """n-dimensional uniform random unit vector"""
    v = np.random.normal(size=n)
    v /= np.linalg.norm(v)
    return v

def vector_random_uniform_hypersphere(k):
    """return a vector sampled uniformly in a hypersphere of dimension k"""
    u = vec_random_ndim(k)
    # draw the magnitude of the vector from a power law density:
    # draws samples in [0, 1] from a power distribution with positive exponent k - 1.
    p = np.random.power(k)
    return p * u
