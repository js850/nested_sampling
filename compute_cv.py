import argparse
import numpy as np


def compute_Z(energies, T, K):
    beta = 1./T
    N = len(energies)
    Z = 0.
    U = 0.
    U2 = 0.
    Cv = 0.
    Emin = energies[-1]
    Ediff = energies - Emin
    for n in xrange(1, len(energies)-2):
#        Z += (np.exp(-float(n-1) / K) - np.exp(-float(n+1) / K)) * np.exp(-beta * energies[n])
        E = Ediff[n]
        Zpref = np.exp(-float(n-1) / K - beta * E) * (1. - np.exp(-2. / K))
        
        
        Z += Zpref
        U += Zpref * (E + Emin)
        U2 += Zpref * (E + Emin)**2
    
    U /= Z
    U2 /= Z
    Cv = (U2 - U**2) * beta**2
        
    return Z, Cv, U, U2


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="load energy intervals and compute cv")
#    parser.add_argument("--db", type=str, nargs=1, help="database filename",
#                        default="otp.db")
    parser.add_argument("K", type=int, help="number of replicas")
    parser.add_argument("fname", type=str, help="filenames with energies")
    args = parser.parse_args()
    print args.fname

    energies = np.genfromtxt(args.fname)
    
    Tmin = .02
    Tmax = .5
    nT = 300
    dT = (Tmax-Tmin) / nT
    
    T = np.array([Tmin + dT*i for i in range(nT)])
    Z, Cv, U, U2 = compute_Z(energies, T, args.K)
    print Z, Cv
    
    with open("cv", "w") as fout:
        fout.write("#T Cv <E> <E**2>\n")
        for vals in zip(T, Cv, U, U2):
            fout.write("%g %g %g %g\n" % vals)
    
    import pylab as pl
    pl.plot(T, Cv, 'o-')
    pl.xlabel("T")
    pl.ylabel("Cv")
    pl.savefig("cv.pdf")
        
    
