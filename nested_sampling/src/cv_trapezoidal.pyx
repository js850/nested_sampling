from cpython cimport bool
import sys
import numpy as np
cimport numpy as np

cdef extern:
    void compute_dos_log(double* gl, int N, double P, double K, int live)

cdef extern:
    void compute_dos_alpha_log(double* gl, double* rn, int N, double K, int live)
    
cdef extern: 
    void heat_capacity_loop(double* El, double* gl, double* wl, double* Cvl, 
                            double* U, double* U2, double * Tlist,
                            int N, int nT, double ndof, 
                            )
    
def compute_cv_c(np.ndarray[double, ndim=1, mode="c"] E_list,
                 double P, double K, double Tmin, double Tmax,
                 int nT, double ndof, bool live):
    cdef int N
    cdef np.ndarray[double, ndim=1, mode="c"] dos_list
    cdef np.ndarray[double, ndim=1, mode="c"] logw_list
    cdef np.ndarray[double, ndim=1, mode="c"] cv_list = np.zeros(nT)
    cdef np.ndarray[double, ndim=1, mode="c"] U = np.zeros(nT)
    cdef np.ndarray[double, ndim=1, mode="c"] U2 = np.zeros(nT)
    cdef np.ndarray[double, ndim=1, mode="c"] T = np.zeros(nT)
    
    N = np.size(E_list)
    dos_list = np.zeros(N)
    logw_list = np.zeros(N)
#    cv_list = np.zeros(nT)
#    U = np.zeros(nT)
#    U2 = np.zeros(nT)
      
    compute_dos_log(<double*>dos_list.data, N, P, K, <bint>live)
    
#    print "dos list",dos_list

    T = np.linspace(Tmin, Tmax, nT)
    
    print "len U", len(U), len(U2)
    heat_capacity_loop(<double*>E_list.data, <double*>dos_list.data, 
                       <double*>logw_list.data, <double*>cv_list.data, 
                       <double*>U.data, <double*>U2.data, <double*>T.data,
                       N, nT, ndof,
                       )
    
    print "logw list",logw_list
    
    
    return T, cv_list, U, U2

def compute_alpha_cv_c(np.ndarray[double, ndim=1, mode="c"] E_list,
                       np.ndarray[double, ndim=1, mode="c"] rn_list,
                       double P, double K, double Tmin, double Tmax,
                       int nT, double ndof, bool live):
    
    cdef int N
    cdef int logscale
    cdef np.ndarray[double, ndim=1, mode="c"] dos_list
    cdef np.ndarray[double, ndim=1, mode="c"] logw_list
    cdef np.ndarray[double, ndim=1, mode="c"] cv_list = np.zeros(nT)
    cdef np.ndarray[double, ndim=1, mode="c"] U = np.zeros(nT)
    cdef np.ndarray[double, ndim=1, mode="c"] U2 = np.zeros(nT)
    cdef np.ndarray[double, ndim=1, mode="c"] T = np.zeros(nT)
    
    N = np.size(E_list)
    dos_list = np.zeros(N)
    logw_list = np.zeros(N)
    
    compute_dos_alpha_log(<double*>dos_list.data, <double*>rn_list.data, N, K, live)
    
#    print "dos list", dos_list
    
    T = np.linspace(Tmin, Tmax, nT)
    
    print "len U", len(U), len(U2)
    heat_capacity_loop(<double*>E_list.data, <double*>dos_list.data, 
                       <double*>logw_list.data, <double*>cv_list.data, 
                       <double*>U.data, <double*>U2.data, <double*>T.data,
                       N, nT, ndof,
                       )
    
    print "logw list",logw_list
    
    
    return T, cv_list, U, U2
