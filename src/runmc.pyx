import sys
import numpy as np
cimport numpy as np

#cimport _random_displace
cdef extern:
    int mc(double *x0, double *xreturn, int natoms, long int mciter, double stepsize, double Emax, double radius,
           long int seed)

def mc_cython(np.ndarray[double, ndim=1, mode="c"] x0,
               int mciter, double stepsize, double Emax, double radius,
               ):
    cdef np.ndarray[double, ndim=1, mode="c"] x
    x = np.zeros(x0.size);
    seed = np.random.randint(0, sys.maxint)

    mc(<double*>x0.data, <double*>x.data, len(x0)/3, mciter, stepsize, Emax, radius, seed)
    return x
