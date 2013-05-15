import argparse
import numpy as np

def compute_Z(energies, T, K, P=1, ndof=0):
    """
    compute the heat capacity and other quantities from nested sampling history
    
    Parameters
    ----------
    energies : list of floats
        list of Emax energy limits
    T : list of floats
        temperatures for which to perform the calculation
    K : int
        number of nested sampling replicas
    P : int
        number of cores used in parallel
    ndof : int
        number of degrees of freedom
        
    Notes
    -----
    from the formula::
    
        alpha = K / (K+1)
        Z = sum_n( (alpha**(n-1) - alpha**(n+1) * np.exp(-beta * energies[n]) )
    
    this can be rewritten::
    
        exp(x) = lim (n->inf) ( (1 + x/n)**n )
        Z = sum_n( (np.exp(-float(n-1) / K) - np.exp(-float(n+1) / K)) * np.exp(-beta * energies[n]) )
    
    or, even better, to try to avoid overflow errors, pull out a factor `np.exp(-float(n-1) / K)`
    
    for parallel, the situation is the same except with a different alpha::
    
        alpha = 1. - float(P) / (K + 1.)

    """
    beta = np.array(1./T)    
    K = float(K)
    E = energies[1:-1]
    n = np.array(xrange(1,len(E)+1), np.float)
    
    # compute the term in the partition function sum for each temperature and each energy
    # lZ[iT, iE] is the contribution to the partition function at temperature 1./beta[iT] from the data at energy E[iE]
    if P == 1:
#        lZ = (-(n[np.newaxis,:]-1.) / K - beta[:,np.newaxis] * E[np.newaxis,:])  + np.log(1. - np.exp(-2. / K))
        lZ = (-(n[np.newaxis,:]) / K - beta[:,np.newaxis] * E[np.newaxis,:])  - np.log(K)
#        lZ = (- beta[:,np.newaxis] * E[np.newaxis,:])  - np.log(K) + n[np.newaxis,:] * np.log(K/(K+1))
#        lZ = (-(n[np.newaxis,:]+1) / K - beta[:,np.newaxis] * E[np.newaxis,:])  - np.log((1.+2*K)/K**2)
#        lZ = (- beta[:,np.newaxis] * E[np.newaxis,:])  - np.log((1.+2*K)/K**2) + (n[np.newaxis,:]+1) * np.log(K/(K+1))
    else:
        a = 1. - float(P) / (K + 1.)
        lZ = n[np.newaxis,:] * np.log(a) + (-beta[:,np.newaxis] * E[np.newaxis,:]) + np.log(1 - a)
        
    # subtract out the smallest value to avoid overflow issues when lZ is exponentiated
    lZmax = np.max(lZ,axis=1) #  maximum lZ for each temperature
    lZ -= lZmax[:,np.newaxis]

    # compute Z, <E> and <E**2> 
    Zpref = np.exp(lZ)
    Z = np.sum(Zpref, axis=1 )
    U = np.sum(Zpref * E[np.newaxis,:], axis=1 )
    U2 = np.sum(Zpref * E[np.newaxis,:]**2, axis=1 )
    
    U /= Z
    U2 /= Z
    # compute Cv from the energy fluctuations
    Cv = (U2 - U**2) * beta**2 + float(ndof)/2 # this is not quite right
    
    lZfinal = np.log(Z) + lZmax
#    Z *= np.exp(lZmax)
        
    return lZfinal, Cv, U, U2


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="load energy intervals and compute cv")
#    parser.add_argument("--db", type=str, nargs=1, help="database filename",
#                        default="otp.db")
    parser.add_argument("K", type=int, help="number of replicas")
    parser.add_argument("fname", type=str, help="filenames with energies")
    parser.add_argument("-P", type=int, help="number of cores for parallel run", default=1)
    parser.add_argument("--ndof", type=int, help="number of degrees of freedom", default=0)
    args = parser.parse_args()
    print args.fname

    energies = np.genfromtxt(args.fname)
    
    P = args.P
    print "parallel nprocessors", P
    
    Tmin = .02
    Tmax = 8
    nT = 1000
    dT = (Tmax-Tmin) / nT
    
    T = np.array([Tmin + dT*i for i in range(nT)])
    lZ, Cv, U, U2 = compute_Z(energies, T, args.K, P=P, ndof=args.ndof)
    
    with open("cv", "w") as fout:
        fout.write("#T Cv <E> <E**2> logZ\n")
        for vals in zip(T, Cv, U, U2, lZ):
            fout.write("%g %g %g %g %g\n" % vals)
    
    import pylab as pl
    pl.plot(T, Cv, 'o-')
    pl.xlabel("T")
    pl.ylabel("Cv")
    pl.savefig("cv.pdf")
        
    