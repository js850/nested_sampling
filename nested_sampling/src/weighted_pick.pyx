import sys
import numpy as np
cimport numpy as np
cimport cython

@cython.boundscheck(False) # turn of bounds-checking for entire function
def weighted_pick_cython(np.ndarray[double, ndim=1, mode="c"] weights):
    cdef unsigned int nweights = len(weights)
    cdef unsigned int i
    cdef double r, s
    
    if nweights == 0:
        raise ValueError("weights must not have zero length")
    r = np.random.uniform(0., np.sum(weights))
    s = 0.0
  
    for i in xrange(nweights):
        s += weights[i]
        if r < s: break
    
    return i