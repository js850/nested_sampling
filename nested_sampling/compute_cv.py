from __future__ import division
import argparse
import numpy as np
import copy
from src.cv_trapezoidal import compute_cv_c

def get_energies(fnames,block=False):
    if len(fnames) == 1:
        return np.genfromtxt(fnames[0])
    if block == False:
        eall = []
        for fname in fnames:
            e = np.genfromtxt(fname)
            eall += e.tolist()
        eall.sort(key=lambda x: -x)
        return np.array(eall).flatten()
    else:
        eall = [[] for i in xrange(len(fnames))]
        for fname,i in zip(fnames,xrange(len(fnames))):
            e = np.genfromtxt(fname)
            eall[i] = copy.deepcopy(e.tolist())
        return np.array(eall)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="load energy intervals and compute cv", 
                                     epilog="if more than one file name is given the energies from all runs will be combined and sorted."
                                     "  the number of replicas will be the sum of the replicas used from all runs (automated!!!)")
    parser.add_argument("K", type=int, help="number of replicas")
    parser.add_argument("fname", nargs="+", type=str, help="filenames with energies")
    parser.add_argument("-P", type=int, help="number of cores for parallel run", default=1)
    parser.add_argument("--Tmin", type=float,help="set minimum temperature for Cv evaluation (default=0.01)",default=0.01)
    parser.add_argument("--Tmax", type=float,help="set maximum temperature for Cv evaluation (default=0.5)",default=0.5)
    parser.add_argument("--nT", type=int,help="set number of temperature in the interval Tmin-Tmax at which Cv is evaluated (default=500)",default=500)
    parser.add_argument("--ndof", type=int, help="number of degrees of freedom (default=0)", default=0)
    parser.add_argument("--live", action="store_true", help="use live replica energies (default=False), numerically unstable for K>2.5k.",default=False)
    parser.add_argument("--live_not_stored", action="store_true", help="turn this flag on if you're using a set of data that does not contain the live replica.",default=False)
    args = parser.parse_args()
    print args.fname

    energies = get_energies(args.fname,block=False)
    print "energies size", np.size(energies)
    
    P = args.P
    print "parallel nprocessors", P
    
    Tmin = args.Tmin
    Tmax = args.Tmax
    nT = args.nT
    dT = (Tmax-Tmin) / nT
    T = np.array([Tmin + dT*i for i in range(nT)])
    
    #in the improved brkf we save the energies of the replicas at the live replica but the ln(dos) underflows for these, hence this:
    if args.live_not_stored == False:
            energies = energies[:-args.K]
    else:
        assert args.live == False,"cannot use live replica under any circumstances if they have not been saved" 
    
    #make nd-arrays C contiguous 
    energies = np.array(energies, order='C')
        
    print "trapezoidal"
    Cv = compute_cv_c(energies, float(P), float(args.K*len(args.fname)), float(Tmin), float(Tmax), nT, float(args.ndof), args.live)
    
    with open("cv.dat", "w") as fout:
        fout.write("#T Cv\n")
        for vals in zip(T, Cv):
            fout.write("%.30f %.30f\n" % vals)
                
    import matplotlib
    matplotlib.use('PDF')
    import pylab as pl
    pl.plot(T, Cv)
    pl.xlabel("T")
    pl.ylabel("Cv")
    pl.savefig("cv.pdf")
        
    
