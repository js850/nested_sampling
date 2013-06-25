import sys
import numpy as np
cimport numpy as np

cdef extern:
    void compute_dos(double* gl, int N, double P, double K)

cdef extern:
    void compute_dos_imp(double* gl, int N, double P, double K)

cdef extern:
    void renorm_energies(double* El, int N, double Emin)
    
cdef extern: 
    void heat_capacity_loop(double* El, double* gl, double* wl, double* Cvl, int N, double Tmin, double Tmax, int nT, double ndof)
    
def compute_cv_c(np.ndarray[double, ndim=1, mode="c"] E_list,
                 double P, double K, double Tmin, double Tmax,
                 int nT, double ndof, int imp):
    cdef double Emin
    cdef int N
    cdef np.ndarray[double, ndim=1, mode="c"] dos_list
    cdef np.ndarray[double, ndim=1, mode="c"] logw_list
    cdef np.ndarray[double, ndim=1, mode="c"] cv_list
    #Emin = E_list[-1]
    N = np.size(E_list)
    dos_list = np.zeros(N)
    logw_list = np.zeros(N)
    cv_list = np.zeros(nT)
        
    #renorm_energies(<double*>E_list.data, N, Emin)
    
    if imp == 0:
        compute_dos(<double*>dos_list.data, N, P, K)
    else:
        compute_dos_imp(<double*>dos_list.data, N, P, K)
    
    heat_capacity_loop(<double*>E_list.data,<double*>dos_list.data,<double*>logw_list.data,<double*>cv_list.data, N, Tmin, Tmax, nT, ndof)
    
    return cv_list

#cdef extern:
#    double heat_capacity(double* El, double* wl, int N, double T, double ndof)
    
#def compute_cv_c_single(np.ndarray[double, ndim=1, mode="c"] E_list, 
#                      double K, double P, double T, double ndof, int imp):
#    cdef double Emin
#    cdef int N
#    cdef np.ndarray[double, ndim=1, mode="c"] dos_list
#    Emin = E_list[-1]
#    N = np.size(E_list)
#    dos_list = np.zeros(N)
#    
#    renorm_energies(<double*>E_list.data, N, Emin)
#    if (P == 1) or (imp == 0): 
#        compute_dos(<double*>dos_list.data, N, P, K)
#    else:
#        compute_dos_imp(<double*>dos_list.data, N, P, K)
#    return heat_capacity(<double*>E_list.data,<double*>dos_list.data, N, T, ndof)
