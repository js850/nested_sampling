import argparse
import numpy as np
from itertools import izip, chain
from nested_sampling import compute_heat_capacity,run_jackknife_variance, run_alpha_variance, get_energies

def main():   
    parser = argparse.ArgumentParser(description="load energy intervals and compute cv", 
                                     epilog="if more than one file name is given the energies from all runs will be combined and sorted."
                                     "  the number of replicas will be the sum of the replicas used from all runs (you need to input this number!)")
    parser.add_argument("K", type=int, help="number of replicas")
    parser.add_argument("fname", nargs="+", type=str, help="filenames with energies")
    parser.add_argument("-P", type=int, help="number of cores for parallel run", default=1)
    parser.add_argument("--Tmin", type=float,help="set minimum temperature for Cv evaluation (default=0.01)",default=0.01)
    parser.add_argument("--Tmax", type=float,help="set maximum temperature for Cv evaluation (default=0.5)",default=0.5)
    parser.add_argument("--nT", type=int,help="set number of temperature in the interval Tmin-Tmax at which Cv is evaluated (default=500)",default=500)
    parser.add_argument("--ndof", type=int, help="number of degrees of freedom (default=0)", default=0)
    parser.add_argument("--live", action="store_true", help="use live replica energies (default=False)",default=False)
    parser.add_argument("-o", type=str, default="cv", help="change the prefix of the output files")
    parser.add_argument("--stdev", type=str, default="none", help="choose what kind of error analysis want to perform: jackknife, alpha, none. Default: none")
    parser.add_argument("--B", action="store_true", help="keep data as they are in blocks (set true)," 
                                                "by default is randomised for a single set of data,"
                                                "while multiple sets are used as they are. To set B true, all sets must have the same K", default=False)
    parser.add_argument("-n", type=int, default=0, help="if --B=0 determine the number of subsets in which you want to split the data ")
    args = parser.parse_args()
    print args.fname
    print args
    
    #===========================================================================
    # determine the number of subsets in which to split the data
    #===========================================================================
    
    if args.B is True:
        nsubsets = len(args.fname)
        if args.live is True:
            nsubsets /= 2
            print "nsubsets=",nsubsets
    else:
        nsubsets = args.n        
        
    #===========================================================================
    # read in the energies
    #===========================================================================
    energies = get_energies(args.fname, args.B, args.live)
    print "energies size", np.size(energies)
    
    #===========================================================================
    # generate a merged set of energies to compute an unbiased estimate of the heat capacity
    #===========================================================================
    if args.B is True:
        energies_merged = [l for l in chain.from_iterable(energies)]
        energies_merged = np.array(np.sort(energies_merged)[::-1])
    else:
        energies_merged = energies
    
    print "parallel nprocessors", args.P
    print "replicas", args.K
    print "error analysis method: ",args.stdev
    
    #make nd-arrays C contiguous # js850> this will already be the case 
    #energies = np.array(energies, order='C')
    
    #===========================================================================
    # do the computation
    #===========================================================================
    T, Cv, U, U2 = compute_heat_capacity(energies_merged, args.K, npar=args.P, 
                                         ndof=args.ndof, Tmin=args.Tmin, Tmax=args.Tmax, 
                                         nT=args.nT, live_replicas=args.live)
    if args.stdev == "jackknife":
        T, Cv_stdev, Cv_singles, CvMom1 =  run_jackknife_variance(energies, nsubsets, args.K, 
                                                                  args.Tmin, args.Tmax, args.nT, args.P, 
                                                                  args.ndof, args.B, args.live)
        #expression for the elimination of the leading piece of bias in the mean
        Cv_best = nsubsets * Cv - (nsubsets -1 ) * CvMom1
    
    elif args.stdev == "alpha":
        T, Cv_stdev, Cv_singles, CvMom1 =  run_alpha_variance(energies, nsubsets, args.K, 
                                                                  args.Tmin, args.Tmax, args.nT, args.P, 
                                                                  args.ndof, args.B, args.live)
    
    #===========================================================================
    # print data
    #===========================================================================
    
    #print just the average estimate of thermodynamic quantities
    with open(args.o+".dat", "w") as fout:
        fout.write("#T Cv <E> <E^2>\n")
        for vals in izip(T, Cv, U, U2):
            fout.write("%.16g %.16g %.16g %.16g\n" % vals)

    if args.stdev is not "none":
        # print Cv with associated standard deviation
        with open('{o}_std_K{K}_Nsub{n}_d{ndof}_B{B}.dat'.format(o=args.o,K = args.K,n=nsubsets,ndof=args.ndof,B=args.B), "w") as fout:
            fout.write("#T Cv stdev\n")
            for vals in zip(T, Cv, Cv_stdev):
                fout.write("%g %g %g\n" % vals)
            
        #CvMom1 is also the jackknife estimate (average of jackknife averages)
        with open('{o}_{m}_variance_K{K}_Nsub{n}_d{ndof}_B{B}.dat'.format(o=args.o,m=args.stdev,K = args.K,n=nsubsets,ndof=args.ndof,B=args.B), "w") as fout:
            fout.write("#T CvMom1 stdev\n")
            for vals in zip(T, CvMom1, Cv_stdev):
                fout.write("%g %g %g\n" % vals)
        
        #single jackknife heat capacity curves
        with open('{o}_{m}_singles_K{K}_Nsub{n}_d{ndof}_B{B}.dat'.format(o=args.o,m=args.stdev,K = args.K,n=nsubsets,ndof=args.ndof,B=args.B), "w") as fout:
            for i in xrange(nsubsets):
                fout.write("#T Cv\n")
                for vals in zip(T, Cv_singles[i]):
                    fout.write("%g %g \n" % vals)
        
        if args.stdev is "jackknife":    
            with open('{o}_jack_best_K{K}_Nsub{n}_d{ndof}_B{B}.dat'.format(o=args.o,K = args.K,n=nsubsets,ndof=args.ndof,B=args.B), "w") as fout:
                fout.write("#T CvBest stdev\n")
                for vals in zip(T, Cv_best, Cv_stdev):
                    fout.write("%g %g %g\n" % vals)

    #===========================================================================
    # make plots and save them
    #===========================================================================
    
    # make a plot and save it
    import matplotlib
    matplotlib.use('PDF')
    import pylab as plt
    
    if args.stdev is "none":
        plt.plot(T, Cv)
        plt.xlabel("T")
        plt.ylabel("Cv")
        plt.savefig(args.o+".pdf")
    else:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.errorbar(T, Cv, yerr=Cv_stdev,ecolor='g', capsize=None)
        ax.set_xlabel("T")
        ax.set_ylabel("Cv")
        fig.savefig('{o}_stdev_est_K{K}_Nsub{n}_d{ndof}.pdf'.format(o=args.o,K = args.K,n=nsubsets,ndof=args.ndof))
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.errorbar(T, CvMom1, yerr=Cv_stdev,ecolor='g', capsize=None)
        ax.set_xlabel("T")
        ax.set_ylabel("Cv")
        fig.savefig('{o}_{m}_variance_est_K{K}_Nsub{n}_d{ndof}.pdf'.format(o=args.o,m=args.stdev,K = args.K,n=nsubsets,ndof=args.ndof))
        
        plt.figure()
        for i in xrange(nsubsets):
            plt.plot(T, Cv_singles[i])
        plt.xlabel("T")
        plt.ylabel("Cv")
        plt.savefig('{o}_{m}_singles_K{K}_Nsub{n}_d{ndof}.pdf'.format(o=args.o,m=args.stdev,K = args.K,n=nsubsets,ndof=args.ndof))
        
    
if __name__ == "__main__":
    main()
