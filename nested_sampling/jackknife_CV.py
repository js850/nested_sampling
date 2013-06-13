import argparse
import numpy as np
import copy
from compute_cv import compute_Z, get_energies
from itertools import chain

class Jackknife_CV(object):
    
    def __init__(self, energies, nsubsets, K, T, P, ndof, block):
        self.E = np.array(energies)
        self.n = np.floor(float(K)/float(nsubsets))
        self.K = K
        self.nsubsets = nsubsets
        self.T = np.array(T)
        self.P = P
        self.ndof = ndof
        self.block = block
    
    def __call__(self):
        Esplit = self.split_energies()
        EJack = self.jack_E_averages(Esplit)
        CvJack = self.jack_Cv_averages(EJack)
        return self.jack_Cv_stdev(CvJack)
    
    def split_energies_randomly(self):
        """
        split the array of energies into n subsets of size n/K randomly
        """
        Esplit = [[] for i in xrange(self.nsubsets)]
        for x in self.E:
            i = np.random.randint(0,self.nsubsets)
            Esplit[i].append(x)
        Esplit = np.array(Esplit)
        for i in xrange(self.nsubsets):
            print 'Esplit size',i, 'is',np.size(Esplit[i])
        return Esplit
    
    def split_energies_block(self):
        """
        split the array of energies into n subsets of size n/K as provided
        """
        ESplit = self.E
        return Esplit
    
    def split_energies(self):
        """
        returns the correct type of energy subsets, as chosen by the user
        """
        if self.block == False:
            Esplit = self.split_energies_randomly()
        else:
            Esplit = self.split_energies_block()
        return Esplit
    
    def jack_E_averages(self, Esplit): 
        """
        return array of Jacknife averages (more like combined subsets than averages):    
        """
        EJack = [[] for i in xrange(self.nsubsets)]
        for i in xrange(self.nsubsets):
            EJack_tmp = copy.deepcopy(Esplit)
            print 'EJack_tmp shape',np.shape(EJack_tmp) 
            EJack_tmp = np.delete(EJack_tmp, i, 0) 
            #EJack_tmp = np.ravel(Esplit,order='F')
            EJack_tmp = [l for l in chain.from_iterable(EJack_tmp)]
            print np.shape(EJack_tmp)
            EJack[i] = np.sort(EJack_tmp)[::-1]
        print np.shape(EJack)
        EJack = np.array(EJack)
        return EJack
    
    def jack_Cv_averages(self, EJack):
        """
        returns the M(=self.nsubsets) Cv Jackknife averages (from the combined subsets)
        """
        CvJack = np.zeros((self.nsubsets,self.T.size))
        for i in xrange(self.nsubsets):
            CvJack[i][:] = compute_Z(np.array(EJack[i][:]), self.T, self.K, P=self.P, ndof=self.ndof)[1]
        print 'CvJack ',CvJack
        return np.array(CvJack)
    
    def jack_Cv_moments(self, CvJack):
        """
        return Cv expectation value from the Jackknife averages of Cv
        """
        CvMom1 = (float(1)/float(self.nsubsets))*np.sum(CvJack,axis=0)               #first moments (1/self.nsubsets)
        CvMom2 = (float(1)/float(self.nsubsets))*np.sum(np.square(CvJack),axis=0)    #second moments
        print 'CvMom1',CvMom1,'CvMom2',CvMom2
        return CvMom1, CvMom2
    
    def jack_Cv_stdev(self, CvJack):
        """
        returns the stdev associated with the heat capacity, it calculates the variance of the Jackknife
        estimate and then from this finds the standard deviation of the heat capacity estimate obtained 
        from the sample average
        """
        CvMom1, CvMom2 = self.jack_Cv_moments(CvJack)
        sigmasquare_jack = CvMom2 - np.square(CvMom1)
        sigma = np.sqrt(self.nsubsets-1)*np.sqrt(sigmasquare_jack) 
        return sigma
        
def run_jackknife(energies, nsubsets, K, T, P, ndof, block):
    """
    returns the stdev calculated by jackknifing
    """
    Cv_sigma = Jackknife_CV(energies, nsubsets, K, T, P, ndof, block)
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
    parser.add_argument("--B", type=int, help="randomise energies in blocks or keep data as they are," 
                                                "by default is randomised for a single set of data,"
                                                "while multiple sets are used as they are", default=2)
    args = parser.parse_args()
    print args.fname

    #by default define automatically weather to split things in block or not
    if args.B is 2:
        args.B = len(args.fname) > 1

    energies = get_energies(args.fname,args.B)
    P = args.P
    print "parallel nprocessors", P
    
    Tmin = .02
    Tmax = 1
    nT = 1000
    dT = (Tmax-Tmin) / nT
    
    T = np.array([Tmin + dT*i for i in range(nT)])
    lZ, Cv, U, U2 = compute_Z(energies, T, args.K, P=P, ndof=args.ndof)
    Cv_stdev = run_jackknife(energies, args.N, args.K, T, P=P, ndof=args.ndof, block=args.B)
    
    with open("cv", "w") as fout:
        fout.write("#T Cv stdev <E> <E**2> logZ\n")
        for vals in zip(T, Cv, Cv_stdev, U, U2, lZ):
            fout.write("%g %g %g %g %g %g\n" % vals)
    
    import matplotlib.pyplot as plt
    plt.figure()
    plt.errorbar(T, Cv, yerr=Cv_stdev, ecolor='g', capsize=None)
    plt.xlabel("T")
    plt.ylabel("Cv")
    plt.savefig("cv_std.pdf")
        
