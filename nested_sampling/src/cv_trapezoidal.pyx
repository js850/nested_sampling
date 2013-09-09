from cpython cimport bool
import sys
import numpy as np
cimport numpy as np

cdef extern:
    void compute_dos(double* gl, int N, double P, double K, int live)

cdef extern:
    void compute_dos_imp(double* gl, int N, double P, double K, int live)

cdef extern:
    void compute_dos_log(double* gl, int N, double P, double K, int live)

cdef extern:
    void renorm_energies(double* El, int N, double Emin)

cdef extern:
    void compute_dos_alpha_log(double* gl, double* rn, int N, double K, int live)
    
cdef extern: 
    void heat_capacity_loop(double* El, double* gl, double* wl, double* Cvl, int N, double Tmin, double Tmax, int nT, double ndof, int logscale)
    
def compute_cv_c(np.ndarray[double, ndim=1, mode="c"] E_list,
                 double P, double K, double Tmin, double Tmax,
                 int nT, double ndof, int imp, bool live):
    cdef double Emin
    cdef int N
    cdef int logscale
    cdef np.ndarray[double, ndim=1, mode="c"] dos_list
    cdef np.ndarray[double, ndim=1, mode="c"] logw_list
    cdef np.ndarray[double, ndim=1, mode="c"] cv_list
    #Emin = E_list[-1]
    N = np.size(E_list)
    dos_list = np.zeros(N)
    logw_list = np.zeros(N)
    cv_list = np.zeros(nT)
    logscale = 1
    
    #renorm_energies(<double*>E_list.data, N, Emin)
    
    if imp == 0 and logscale == 0:
        compute_dos(<double*>dos_list.data, N, P, K, <bint>live)
    elif logscale == 0:
        compute_dos_imp(<double*>dos_list.data, N, P, K, <bint>live)
    else:
        compute_dos_log(<double*>dos_list.data, N, P, K, <bint>live)
    
    print "dos list",dos_list
    
    heat_capacity_loop(<double*>E_list.data,<double*>dos_list.data,<double*>logw_list.data,<double*>cv_list.data, N, Tmin, Tmax, nT, ndof, logscale)
    
    print "logw list",logw_list
    
    return cv_list

def compute_alpha_cv_c(np.ndarray[double, ndim=1, mode="c"] E_list,
                       np.ndarray[double, ndim=1, mode="c"] rn_list,
                       double P, double K, double Tmin, double Tmax,
                       int nT, double ndof, int imp, bool live):
    cdef double Emin
    cdef int N
    cdef int logscale
    cdef np.ndarray[double, ndim=1, mode="c"] dos_list
    cdef np.ndarray[double, ndim=1, mode="c"] logw_list
    cdef np.ndarray[double, ndim=1, mode="c"] cv_list
    
    N = np.size(E_list)
    dos_list = np.zeros(N)
    logw_list = np.zeros(N)
    cv_list = np.zeros(nT)
    
    logscale = 1
    ##sample the compression ratios
        
    if imp == 0 or logscale == 0:
        raise ValueError('Method only implemented for improved brkf and log-space')
    else:
        compute_dos_alpha_log(<double*>dos_list.data, <double*>rn_list.data, N, K, live)
    
    print "dos list",dos_list
    
    heat_capacity_loop(<double*>E_list.data,<double*>dos_list.data,<double*>logw_list.data,<double*>cv_list.data, N, Tmin, Tmax, nT, ndof, logscale)
    
    print "logw list",logw_list
    
    return cv_list
