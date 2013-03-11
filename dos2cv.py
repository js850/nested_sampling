from itertools import izip
import argparse
import numpy as np


def compute_cv_old(energies, ldos, T, K=1.):
    """compute the heat capacity and other thermodynamic quantities from density of states
    
    Parameters
    ----------
    energies : array
        list of energies
    ldos : array
        log of density of states at energy given by energies
    T : float or array of floats
        temperatures at which to do the calculations.  should be normalized by the Bolzmann constant (you should pass k_B*T)
    K : int
        number of degrees of vibrational freedom (3*N-6 for clusters) 
    """
    beta = 1./T
    Z = 0.
    U = 0.
    U2 = 0.
    Emin = min(energies)
    for E, lomega in izip(energies, ldos):
        Zpref = np.exp(lomega - beta * (E - Emin))        
        
        Z += Zpref
        U += Zpref * E
        U2 += Zpref * E**2
    
    U /= Z
    U2 /= Z
    Cv = float(K)/2 + (U2 - U**2) * beta**2 # this is not quite right
    return Z, Cv, U, U2

def compute_cv(energies, ldos, T, K=1.):
    """compute the heat capacity and other thermodynamic quantities from density of states
    
    this version is faster and should be more accurate
    
    Parameters
    ----------
    energies : array
        list of energies
    ldos : array
        log of density of states at energy given by energies
    T : float or array of floats
        temperatures at which to do the calculations.  should be normalized by the Bolzmann constant (you should pass k_B*T)
    K : int
        number of degrees of vibrational freedom (3*N-6 for clusters) 
    """
    beta = np.array(1./T)
    # lZ[T, E] = log(  Omega(E) * exp(-E / T) 
    lZ =  ldos[np.newaxis,:] - beta[:,np.newaxis] * energies[np.newaxis,:]
    # subtract out the smallest value to avoid overflow issues when lZ is exponentiated
    lZmax = np.max(lZ,axis=1) #  maximum lZ for each temperature
    lZ -= lZmax[:,np.newaxis]

    # compute Z, <E> and <E**2> 
    Zpref = np.exp(lZ)
    Z = np.sum(Zpref, axis=1 )
    U = np.sum(Zpref * energies[np.newaxis,:], axis=1 )
    U2 = np.sum(Zpref * energies[np.newaxis,:]**2, axis=1 )
    
    U /= Z
    U2 /= Z
    # compute Cv from the energy fluctuations
    Cv = float(K)/2 + (U2 - U**2) * beta**2 # this is not quite right
    Z *= np.exp(lZmax)
    return Z, Cv, U, U2


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="load density of states and compute cv.  energies in first column, log(dos) in second column")
#    parser.add_argument("--db", type=str, nargs=1, help="database filename",
#                        default="otp.db")
    parser.add_argument("K", type=int, help="number of vibrational degrees of freedom")
    parser.add_argument("fname", type=str, help="filenames with energies and dos")
    args = parser.parse_args()
    print args.fname

    data = np.genfromtxt(args.fname)
    energies = data[:,0]
    ldos = data[:,0]
    
    Tmin = .01
    Tmax = .6
    nT = 1000
    dT = (Tmax-Tmin) / nT
    
    T = np.array([Tmin + dT*i for i in range(nT)])
    Z, Cv, U, U2 = compute_cv(energies, ldos, T, args.K)
    
    print "saving data to cv"
    with open("cv", "w") as fout:
        fout.write("#T Cv <E> <E**2>\n")
        for vals in zip(T, Cv, U, U2):
            fout.write("%g %g %g %g\n" % vals)
    
    import pylab as pl
    pl.plot(T, Cv, '-')
    pl.xlabel("T")
    pl.ylabel("Cv")
    print "saving plot to cv.pdf"
    pl.savefig("cv.pdf")
        
    
