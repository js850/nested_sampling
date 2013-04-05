"""routines to calculate Cv from the harmonic approximation"""

import argparse
import numpy as np


def compute_cv(minima, T, k):
    """compute the heat capacity and other thermodynamic quantities from a database using the harmonic approximation
    
    Parameters
    ----------
    minima : list of Minimum objects
        mimima from which to compute the thermodynamic computations
    T : float or array of floats
        temperatures at which to do the calculations.  should be normalized by the Bolzmann constant (you should pass k_B*T)
    k : int
        number of degrees of vibrational freedom (3*N-6 for clusters) 
    
    Notes
    -----
    See DJW Energy Landscapes book page 371
    Za = P * exp(-beta Ea) / (beta * h * nu_bar)**k / O_a
    """
    beta = np.array(1./T)
    k = float(k)
    Emin = min([m.energy for m in minima])
    Z = 0.
    U = 0.
    U2 = 0.
    for m in minima:
        E = m.energy
        lZpref = -beta * (E - Emin) - m.fvib/2. - np.log(m.pgorder)
        Zpref = np.exp(lZpref)
        
        Z += Zpref
        U += Zpref * E
        U2 += Zpref * E**2
    
    U /= Z
    U2 /= Z
    Cv = k + (U2 - U**2) * beta**2   # this is not completely correct
    Z *= np.exp(beta * Emin) / beta**k
    return Z, U, U2, Cv
        
if __name__ == "__main__":
    from pygmin.storage import Database
    parser = argparse.ArgumentParser(description="load energy intervals and compute cv")
#    parser.add_argument("--db", type=str, nargs=1, help="database filename",
#                        default="otp.db")
    parser.add_argument("k", type=int, help="number of degrees of vibrational freedom")
    parser.add_argument("fname", type=str, help="database file name")
    args = parser.parse_args()
    print args.fname
    k = args.k
    
    dbfname = args.fname
    db = Database(dbfname)

    Tmin = .001
    Tmax = .5
    nT = 300
    dT = (Tmax-Tmin) / nT
    
    T = np.array([Tmin + dT*i for i in range(nT)])
    Z, U, U2, Cv = compute_cv(db.minima(), T, k)
    print Z, Cv
    
    with open("cv_ha", "w") as fout:
        fout.write("#T Cv <E> <E**2>\n")
        for vals in zip(T, Cv, U, U2):
            fout.write("%g %g %g %g\n" % vals)
    
    import pylab as pl
    pl.plot(T, Cv, 'o-')
    pl.xlabel("T")
    pl.ylabel("Cv")
    pl.savefig("cv_ha.pdf")
        
