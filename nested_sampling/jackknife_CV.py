import argparse
import numpy as np
import copy
from compute_cv import compute_Z, get_energies

class Jackknife_CV(object):
    
    def __init__(self, energies, nsubsets, K, T, P, ndof):
        self.E = np.array(energies)
        self.n = np.floor(float(K)/float(nsubsets))
        self.K = K
        self.nsubsets = nsubsets
        self.T = np.array(T)
        self.P = P
        self.ndof = ndof
    
    def __call__(self):
        Esplit = self.split_energies()
        EJack = self.jack_E_averages(Esplit)
        CvJack = self.jack_Cv_averages(EJack)
        return self.jack_Cv_stdev(EJack)
    
    def split_energies(self):
        """
        split the array of energies into n subsets of size n/K
        """
        Esplit = [[] for i in xrange(self.nsubsets)]
        for x in self.E:
            i = np.random.randint(0,self.nsubsets-1)
            Esplit[i].append(x)
        Esplit = np.array(Esplit)
        return Esplit
    
    def jack_E_averages(self, Esplit): 
        """
        return array of Jacknife averages (more like combined subsets than averages):    
        """
        EJack = [[] for i in xrange(self.nsubsets)]
        for i in xrange(self.nsubsets):
            EJack_tmp = copy.deepcopy(Esplit)
            EJack_tmp = np.delete(EJack_tmp, i, 0) 
            EJack_tmp = np.ravel(Esplit,order='F') 
            EJack[i] = np.sort(EJack_tmp)[::-1]
        EJack = np.array(EJack)
        return EJack
    
    def jack_Cv_averages(self, EJack):
        """
        returns the M(=self.nsubsets) Cv Jackknife averages (from the combined subsets)
        """
        CvJack = np.zeros((self.nsubsets,self.T.size))
        for i in xrange(self.nsubsets):
            CvJack[i,:] = compute_Z(np.array(EJack[i]), self.T, self.K, P=self.P, ndof=self.ndof)[1]
        return CvJack
    
    def jack_Cv_moments(self, CvJack):
        """
        return Cv expectation value from the Jackknife averages of Cv
        """
        CvMom1 = (1/self.nsubsets) * np.sum(CvJack,axis=0)               #first moments
        CvMom2 = (1/self.nsubsets) * np.sum(np.square(CvJack),axis=0)    #second moments
        return CvMom1, CvMom2
    
    def jack_Cv_stdev(self, CvJack):
        """
        returns the stdev associated with the heat capacity, it calculates the variance of the Jackknife
        estimate and then from this finds the standard deviation of the heat capacity estimate obtained 
        from the sample average
        """
        CvMom1, CvMom2 = self.jack_Cv_moments(self, CvJack)
        sigmasquare_jack = CvMom2 - np.square(CvMom1)
        sigma = sqrt(self.nsubsets-1)*np.sqrt(sigmasquare_jack) 
        return sigma
        
def run_jackknife(energies, nsubsets, K, T, P, ndof):
    """
    returns the stdev calculated by jackknifing
    """
    Cv_sigma = Jackknife_CV(energies, nsubsets, K, T, P, ndof)
    return Cv_sigma()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="load energy intervals and compute Cv stdev", 
                                     epilog="if more than one file name is given the energies from all runs will be combined and sorted."
                                     "  the number of replicas must be the sum of the replicas used from all runs")
    parser.add_argument("K", type=int, help="number of replicas")
    parser.add_argument("N", type=int, help="number of subsets for jackknifing")
    parser.add_argument("fname", nargs="+", type=str, help="filenames with energies")
    parser.add_argument("-P", type=int, help="number of cores for parallel run", default=1)
    parser.add_argument("--ndof", type=int, help="number of degrees of freedom", default=0)
    args = parser.parse_args()
    print args.fname

    energies = get_energies(args.fname)
    P = args.P
    print "parallel nprocessors", P
    
    Tmin = .02
    Tmax = 1
    nT = 1000
    dT = (Tmax-Tmin) / nT
    
    T = np.array([Tmin + dT*i for i in range(nT)])
    lZ, Cv, U, U2 = compute_Z(energies, T, args.K, P=P, ndof=args.ndof)
    Cv_stdev = run_jackknife(energies, args.N, args.K, T, P=P, ndof=args.ndof)
    
    with open("cv", "w") as fout:
        fout.write("#T Cv stdev <E> <E**2> logZ\n")
        for vals in zip(T, Cv, Cv_stdev, U, U2, lZ):
            fout.write("%g %g %g %g %g %g\n" % vals)
    
    import pylab as pl
    pl.figure()
    pl.plot.errorbars(T, Cv, yerr=Cv_stdev)
    pl.xlabel("T")
    pl.ylabel("Cv")
    pl.savefig("cv_std.pdf")
        
