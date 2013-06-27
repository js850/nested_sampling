from __future__ import division
import argparse
import numpy as np
import copy
from compute_cv import compute_Z, get_energies
from src.cv_trapezoidal import  compute_cv_c
from itertools import chain

class Jackknife_CV(object):
    
    def __init__(self, energies, nsubsets, K, Tmin, Tmax, nT, P, ndof, block, rect, imp, live):
        self.E = np.array(energies)
        self.n = np.floor(float(K)/float(nsubsets))
        self.K = K
        self.nsubsets = nsubsets
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.nT = nT
        self.dT = (Tmax-Tmin) / nT
        self.T = np.array([Tmin + dT*i for i in range(nT)])
        self.P = P
        self.ndof = ndof
        self.block = block
        self.rect = rect
        self.imp = imp
        self.live = live
    
    def __call__(self):
        print 'Splitting energies...'
        Esplit = self.split_energies()
        print 'Calculating Jacknife averages (combining sets)'
        EJack = self.jack_E_averages(Esplit)
        print 'Producing Cv curves for each subset'
        CvSingle = self.Cv_singles(Esplit)
        print 'Calculating Cv Jacknife average (from the combination of subsets)'
        CvJack = self.jack_Cv_averages(EJack)
        return self.jack_Cv_stdev(CvJack), CvSingle
    
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
        Esplit = copy.deepcopy(self.E)
        return Esplit
    
    def split_energies(self):
        """
        returns the correct type of energy subsets, as chosen by the user
        """
        if self.block == False:
            print 'splitting energies at random in ',self.nsubsets,' subsets'
            Esplit = self.split_energies_randomly()
        else:
            print 'keeping energies as from input'
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
        if self.rect is 1:
            for i in xrange(self.nsubsets):
                CvJack[i][:] = compute_Z(np.array(EJack[i][:]), self.T, (self.K - self.n), P=self.P, ndof=self.ndof, imp=self.imp, live=self.live)[1]
        else:
            for i in xrange(self.nsubsets):
                CvJack[i][:] = compute_cv_c(np.array(EJack[i][:]), float(P), (self.K - self.n), float(self.Tmin), float(self.Tmax), self.nT, float(self.ndof), self.imp, self.live)
        #print 'CvJack ',CvJack
        return np.array(CvJack)
    
    def Cv_singles(self,Esplit):
        """
        returns the single Cvs
        """
        CvSingle = np.zeros((self.nsubsets,self.T.size))
        if self.rect is 1:
            for i in xrange(self.nsubsets):
                CvSingle[i][:] = compute_Z(np.array(Esplit[i][:]), self.T, (self.n), P=self.P, ndof=self.ndof, imp=self.imp, live=self.live)[1]
        else:
            for i in xrange(self.nsubsets):
                CvSingle[i][:] = compute_cv_c(np.array(Esplit[i][:]), float(P), (self.n), float(self.Tmin), float(self.Tmax), self.nT, float(self.ndof), self.imp, self.live)
        #print 'CvSingle ',CvSingle
        return np.array(CvSingle)
    
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
            
def run_jackknife(energies, nsubsets, K, Tmin, Tmax, nT, P, ndof, block, rect, imp, live):
    """
    returns the stdev calculated by jackknifing
    """
    Cv_sigma = Jackknife_CV(energies, nsubsets, K, Tmin, Tmax, nT, P, ndof, block, rect, imp, live)
    return Cv_sigma()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="load energy intervals and compute Cv stdev", 
                                     epilog="if more than one file name is given and Nsubs is same as the number of inputs, the energies from all runs will be kept as they are by default"
                                     " (if want to combine and sort them at random add --B 0)."
                                     "  The number of replicas will be the sum of the replicas used from all runs (automated!!!"
                                        "does not support sets with different number of replicas yet)")
    parser.add_argument("K", type=int, help="number of replicas for each single run (they all must be equal in the present implementation)")
    parser.add_argument("N", type=int, help="number of subsets for jackknifing")
    parser.add_argument("fname", nargs="+", type=str, help="filenames with energies")
    parser.add_argument("-P", type=int, help="number of cores for parallel run", default=1)
    parser.add_argument("--Tmin", type=float,help="set minimum temperature for Cv evaluation (default=0.01)",default=0.01)
    parser.add_argument("--Tmax", type=float,help="set maximum temperature for Cv evaluation (default=0.5)",default=0.5)
    parser.add_argument("--nT", type=int,help="set number of temperature in the interval Tmin-Tmax at which Cv is evaluated (default=500)",default=500)
    parser.add_argument("--ndof", type=int, help="number of degrees of freedom (default=0)", default=0)
    parser.add_argument("--imp", type=int, help="define whether to use improved Burkoff (use all energies and live replica energies (default=1), otherwise set to 0)", default=1)
    parser.add_argument("--rect", type=int, help="0 for trapezoidal from arithmetic mean (default=0),1 for rectangular from geometric mean", default=0)
    parser.add_argument("--live", action="store_true", help="use live replica energies (default=False), numerically unstable for K>2.5k.",default=False)
    parser.add_argument("--live_not_stored", action="store_true", help="turn this flag on if you're using a set of data that does not contain the live replica.",default=False)
    parser.add_argument("--B", type=int, help="randomise energies (set 0) or keep data as they are in blocks (set 1)," 
                                                "by default is randomised for a single set of data,"
                                                "while multiple sets are used as they are", default=2)
    args = parser.parse_args()
    print args.fname
    P = args.P
    
    ####################################################################################################
    #deal with input
    #by default define automatically weather to split things in block or not
    if args.B is 2:
        if (len(args.fname) > 1) and (args.N == len(args.fname)):
            args.B = 1
        else:
            args.B = 0

    energies = get_energies(args.fname,1) #always copy in blocks, then depending on B flatten or not
            
    #in the improved brkf we save the energies of the replicas at the live replica but the ln(dos) underflows for these, hence this:
    if args.live_not_stored == False:
        if len(args.fname) > 1:
            for i in xrange(len(args.fname)):
                    energies[i] = energies[i][:-args.K]
        else:
            energies = energies[:-args.K]
    else:
        assert args.live == False,"cannot use live replica under any circumstances if they have not been saved" 
    
    if len(args.fname) > 1:
        energies_Cv =  [l for l in chain.from_iterable(energies)] #provide the sorted flatten list of energies to calculate the unbiased estimate for the Cv
        energies_Cv = np.sort(energies_Cv)[::-1]
        if args.B is 0:
            energies = energies_Cv #if want to rearrange energies randomly flatten energies
    else:
        energies_Cv = energies
    ##########################################################################################################
     
    print "parallel nprocessors", P
    
    Tmin = args.Tmin
    Tmax = args.Tmax
    nT = args.nT
    dT = (Tmax-Tmin) / nT
    T = np.array([Tmin + dT*i for i in range(nT)])
    
    if args.rect is 1:
        print "rectangular"
        Cv = compute_Z(energies_Cv, T, args.K*len(args.fname), P=P, ndof=args.ndof, imp=args.imp, live=args.live)[1]
    else:
        print "trapezoidal"
        Cv = compute_cv_c(energies_Cv, float(P), float(args.K*len(args.fname)), float(Tmin), float(Tmax), nT, float(args.ndof), args.imp, args.live)
    
    Cv_stdev, Cv_singles = run_jackknife(energies, args.N, (args.K*len(args.fname)), Tmin, Tmax, nT, P, args.ndof, args.B, args.rect, args.imp, args.live)
    
    with open('cv_std_K{K}_Nsub{N}_d{ndof}_B{B}.dat'.format(K = args.K,N=args.N,ndof=args.ndof,B=args.B), "w") as fout:
        fout.write("#T Cv stdev\n")
        for vals in zip(T, Cv, Cv_stdev):
            fout.write("%g %g %g\n" % vals)
    
    with open('cv_singles_K{K}_Nsub{N}_d{ndof}_B{B}.dat'.format(K = args.K,N=args.N,ndof=args.ndof,B=args.B), "w") as fout:
        for i in xrange(args.N):
            fout.write("#T Cv\n")
            for vals in zip(T, Cv_singles[i]):
                fout.write("%g %g \n" % vals)
    
    #####################
    ####### PLOT ########
    #####################
    
    import matplotlib
    matplotlib.use('PDF')
    import matplotlib.pyplot as plt
    fig, (ax1,ax2) = plt.subplots(2,1,sharey=True)
    ax = ax2
    ax.errorbar(T, Cv, yerr=Cv_stdev,ecolor='g', capsize=None)
    ax.set_xlabel("T")
    ax.set_ylabel("Cv")
    ax = ax1
    ax.set_xlim([0,0.5])
    ax.errorbar(T, Cv, yerr=Cv_stdev,ecolor='g', capsize=None)
    ax.set_xlabel("T")
    ax.set_ylabel("Cv")
    fig.savefig('cv_std_K{K}_Nsub{N}_d{ndof}_B{B}.pdf'.format(K = args.K,N=args.N,ndof=args.ndof,B=args.B))
    
    plt.figure()
    for i in xrange(args.N):
        plt.plot(T, Cv_singles[i])
    plt.xlabel("T")
    plt.ylabel("Cv")
    plt.savefig('cv_singles_K{K}_Nsub{N}_d{ndof}_B{B}.pdf'.format(K = args.K,N=args.N,ndof=args.ndof,B=args.B))
        
    fig, (ax1,ax2) = plt.subplots(2,1,sharey=True)
    ax = ax2
    ax.errorbar(T, Cv, yerr=Cv_stdev,ecolor='g', capsize=None )
    for i in xrange(args.N):
        ax.plot(T, Cv_singles[i],'k')
    ax.set_xlabel("T")
    ax.set_ylabel("Cv")
    ax = ax1
    ax.set_xlim([0,0.5])
    ax.errorbar(T, Cv, yerr=Cv_stdev,ecolor='g', capsize=None )
    for i in xrange(args.N):
        ax.plot(T, Cv_singles[i],'k')
    ax.set_xlabel("T")
    ax.set_ylabel("Cv")
    fig.savefig('cv_combined_K{K}_Nsub{N}_d{ndof}_B{B}.pdf'.format(K = args.K,N=args.N,ndof=args.ndof,B=args.B))
    
    

