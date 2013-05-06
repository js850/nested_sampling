import sys
import numpy as np
cimport numpy as np

#cimport _random_displace
cdef extern:
    int mcising(long int *spins, double *Energy, long int nspins, 
                long int niter, double Emax, long int seed,
                long int *neighbor_list, long int *nbegin, long int *nend, double macheps)

def mc_ising_c(np.ndarray[long int, ndim=1, mode="c"] spins,
                long int mciter, double Emax, long int seed,
                np.ndarray[long int, ndim=1, mode="c"] neighbor_list,
                np.ndarray[long int, ndim=1, mode="c"] nbegin,
                np.ndarray[long int, ndim=1, mode="c"] nend,
                double energy
                ):
#    cdef np.ndarray[double, ndim=1, mode="c"] x
#    x = np.zeros(x0.size);
    cdef np.ndarray[double, ndim=1, mode="c"] Eretarray 
    Eretarray = np.zeros(1) + energy

    naccept = mcising(<long int*>spins.data, <double*> Eretarray.data, len(spins), 
                      mciter, Emax, seed,
                      <long int*>neighbor_list.data,
                      <long int*>nbegin.data,
                      <long int*>nend.data
                      )
    Ereturn = Eretarray[0]
    return spins, Ereturn, naccept
