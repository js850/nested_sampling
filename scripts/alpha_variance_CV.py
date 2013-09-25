from __future__ import division
import argparse
import numpy as np
import copy
from compute_cv import compute_Z, get_energies
from src.cv_trapezoidal import compute_cv_c, compute_alpha_cv_c
from itertools import chain

class variance_CV(object):
    
    def __init__(self, energies, nsubsets, K, Tmin, Tmax, nT, P, ndof, imp, live):
        self.E = np.array(energies)
        self.K = K
        self.N = len(self.E)
        self.nsubsets = nsubsets
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.nT = nT
        self.dT = (Tmax-Tmin) / nT
        self.T = np.array([Tmin + dT*i for i in range(nT)])
        self.P = P
        self.ndof = ndof
        self.imp = imp
        self.live = live
                
    def __call__(self):
        #print 'Sampling alphas...'
        #alpha_sets = self.sample_alphas()
        #print 'Calculating Jacknife averages (Jackknife averaging alphas)'
        #alphaJack = self.jack_alpha_averages(alpha_sets)
        #print 'Producing Cv curves for each subset of alphas'
        #CvSingle = self.Cv_singles(alpha_sets)
        #print 'Calculating Cv Jacknife average'
        #CvJack = self.jack_Cv_averages(alphaJack)
        #CvMom1 = self.jack_Cv_moments(CvJack)[0]
        #return self.jack_Cv_stdev(CvJack), CvSingle, CvMom1
        
        #print 'Sampling alphas...'
        #alpha_sets = self.sample_alphas()
        #print 'Calculating Cvs'
        #alphaJack = self.jack_Cv_averages(alpha_sets) 
        #print 'Calculating Jacknife averages'
        #CvJack = self.jack_alpha_averages(alphaJack)
        #print 'Producing Cv curves for each subset of alphas'
        #CvSingle = self.Cv_singles(alpha_sets)
        #print 'Calculating moments'
        #CvMom1 = self.jack_Cv_moments(CvJack)[0]
        #return self.jack_Cv_stdev(CvJack), CvSingle, CvMom1
        
        print 'Sampling alphas...'
        alpha_sets = self.sample_alphas()
        print 'Calculating single Cvs'
        CvSingle = self.Cv_singles(alpha_sets)
        print 'Calculating moments'
        CvMom1, CvMom2 = self.Cv_moments(CvSingle)
        print 'Calculating standard deviation'
        sigma = np.sqrt(CvMom2 - np.square(CvMom1))
        print 'Plotting...'
        return sigma, CvSingle, CvMom1
                
    def make_random_alpha_list(self):
        rn_list = np.zeros(self.N)
        if self.live == 0:
            for i in xrange(self.N):
                rn_list[i] = np.random.beta(self.K-(i%int(self.P)),1)
        else:
            for i in xrange(self.N-self.K):
                rn_list[i] = np.random.beta(self.K-(i%int(self.P)),1)
            for j in xrange(self.K):
                i +=1
                rn_list[i] = np.random.beta(self.K-j,1)
        return np.array(rn_list)
    
    def sample_alphas(self):
        """
        create nsubsets of sampled compression factors 
        """
        alpha_sets = [[] for i in xrange(self.nsubsets)]
        for i in xrange(self.nsubsets):
            alpha = self.make_random_alpha_list()
            alpha_sets[i] = alpha
        alpha_sets = np.array(alpha_sets)
        for i in xrange(self.nsubsets):
            print 'alpha_tot size',i, 'is',np.size(alpha_sets[i])
        return alpha_sets
    
    def jack_alpha_averages(self, alpha_sets): 
        """
        return array of Jacknife averages:    
        """
        alphaJack = [[] for i in xrange(self.nsubsets)]
        for i in xrange(self.nsubsets):
            alphaJack_tmp = copy.deepcopy(alpha_sets)
            print 'alphaJack_tmp shape',np.shape(alphaJack_tmp) 
            alphaJack_tmp = np.delete(alphaJack_tmp, i, 0) 
            alphaJack_tmp =  (1. / (self.nsubsets-1)) * np.sum(alphaJack_tmp, axis=0)
            print np.shape(alphaJack_tmp)
            alphaJack[i] = copy.deepcopy(alphaJack_tmp)
        print np.shape(alphaJack)
        alphaJack = np.array(alphaJack)
        return alphaJack
    
    def jack_Cv_averages(self, alphaJack):
        """
        returns the M(=self.nsubsets) Cv Jackknife averages
        """
        CvJack = np.zeros((self.nsubsets, self.nT))
        for i in xrange(self.nsubsets):
            CvJack[i][:] = compute_alpha_cv_c(self.E , np.array(alphaJack[i][:]), float(P), self.K, float(self.Tmin), float(self.Tmax), self.nT, float(self.ndof), self.imp, self.live)
        #print 'CvJack ',CvJack
        return np.array(CvJack)
    
    def Cv_singles(self,alpha_sets):
        """
        returns the single Cvs
        """
        CvSingle = np.zeros((self.nsubsets,self.T.size))
        for i in xrange(self.nsubsets):
            CvSingle[i][:] = compute_alpha_cv_c(self.E, np.array(alpha_sets[i][:]), float(P), self.K, float(self.Tmin), float(self.Tmax), self.nT, float(self.ndof), self.imp, self.live)
        #print 'CvSingle ',CvSingle
        return np.array(CvSingle)
    
    def Cv_moments(self, CvJack):
        """    
        return Cv expectation value from the Cvs
        """
        CvMom1 = (float(1)/float(self.nsubsets))*np.sum(CvJack,axis=0)               #first moments (1/self.nsubsets)
        CvMom2 = (float(1)/float(self.nsubsets))*np.sum(np.square(CvJack),axis=0)    #second moments
        print 'CvMom1',CvMom1,'CvMom2',CvMom2
        return CvMom1, CvMom2
    
    def Cv_stdev(self, CvJack):
        """
        returns the stdev associated with the heat capacity, it calculates the variance of the estimates
        and then from this finds the standard deviation of the heat capacity estimate obtained 
        from the heat capacity average
        """
        CvMom1, CvMom2 = self.Cv_moments(CvJack)
        sigmasquare_jack = CvMom2 - np.square(CvMom1)
        sigma = np.sqrt(self.nsubsets-1)*np.sqrt(sigmasquare_jack) 
        return sigma
               
def run_variance(energies, nsubsets, K, Tmin, Tmax, nT, P, ndof, imp, live):
    """
    returns the stdev calculated by jackknifing
    """
    Cv_sigma = variance_CV(energies, nsubsets, K, Tmin, Tmax, nT, P, ndof, imp, live)
    return Cv_sigma()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="load energy intervals and compute Cv stdev", 
                                     epilog="if more than one file name is given and Nsubs is same as the number of inputs, the energies from all runs will be kept as they are by default"
                                     " (if want to combine and sort them at random add --B 0)."
                                     "  The number of replicas will be the sum of the replicas used from all runs (automated!!!"
                                        "does not support sets with different number of replicas yet)")
    parser.add_argument("K", type=int, help="number of replicas for each single run (they all must be equal in the present implementation, unless using live replica)")
    parser.add_argument("N", type=int, help="number of subsets for jackknifing")
    parser.add_argument("fname", nargs="+", type=str, help="filenames with energies")
    parser.add_argument("-P", type=int, help="number of cores for parallel run", default=1)
    parser.add_argument("--Tmin", type=float,help="set minimum temperature for Cv evaluation (default=0.01)",default=0.01)
    parser.add_argument("--Tmax", type=float,help="set maximum temperature for Cv evaluation (default=0.5)",default=0.5)
    parser.add_argument("--nT", type=int,help="set number of temperature in the interval Tmin-Tmax at which Cv is evaluated (default=500)",default=500)
    parser.add_argument("--ndof", type=int, help="number of degrees of freedom (default=0)", default=0)
    parser.add_argument("--imp", type=int, help="define whether to use improved Burkoff (use all energies energies (default=1), otherwise set to 0)", default=1)
    parser.add_argument("--live", action="store_true", help="use live replica energies (default=False), numerically unstable for K>2.5k.",default=False)
    parser.add_argument("--live_not_stored", action="store_true", help="turn this flag on if you're using a set of data that does not contain the live replica.",default=False)
    
    args = parser.parse_args()
    print args.fname
    P = args.P
    
    ####################################################################################################
    #deal with input
    energies = get_energies(args.fname,0) #always flattens the input
            
    #if args.live_not_stored == False and args.live == False:
    #    if len(args.fname) > 1:
    #        for i in xrange(len(args.fname)):
    #                energies[i] = energies[i][:-args.K]
    #    else:
    #        energies = energies[:-args.K]
    #elif args.live_not_stored == True:
    #    assert args.live == False,"cannot use live replica under any circumstances if they have not been saved" 
    
    #make nd-arrays C contiguous 
    energies = np.array(energies, order='C')
    ##########################################################################################################
     
    print "parallel nprocessors", P
    
    Tmin = args.Tmin
    Tmax = args.Tmax
    nT = args.nT
    dT = (Tmax-Tmin) / nT
    T = np.array([Tmin + dT*i for i in range(nT)])
    
    Cv_stdev, Cv_singles, CvMom1 = run_variance(energies, args.N, (args.K*len(args.fname)), Tmin, Tmax, nT, P, args.ndof, args.imp, args.live)
    
    #CvMom1 is also the jackknife estimate
    with open('cv_alpha_variance_est_K{K}_Nsub{N}_d{ndof}.dat'.format(K = args.K,N=args.N,ndof=args.ndof), "w") as fout:
        fout.write("#T CvMom1 stdev\n")
        for vals in zip(T, CvMom1, Cv_stdev):
            fout.write("%g %g %g\n" % vals)
    
    with open('cv_alpha_singles_K{K}_Nsub{N}_d{ndof}.dat'.format(K = args.K,N=args.N,ndof=args.ndof), "w") as fout:
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
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(T, CvMom1, yerr=Cv_stdev,ecolor='g', capsize=None)
    ax.set_xlabel("T")
    ax.set_ylabel("Cv")
    fig.savefig('cv_alpha_variance_est_K{K}_Nsub{N}_d{ndof}.pdf'.format(K = args.K,N=args.N,ndof=args.ndof))
    
    plt.figure()
    for i in xrange(args.N):
        plt.plot(T, Cv_singles[i])
    plt.xlabel("T")
    plt.ylabel("Cv")
    plt.savefig('cv_alpha_singles_K{K}_Nsub{N}_d{ndof}.pdf'.format(K = args.K,N=args.N,ndof=args.ndof))
        
    fig, (ax1,ax2) = plt.subplots(2,1,sharey=True)
    ax = ax2
    ax.errorbar(T, CvMom1, yerr=Cv_stdev,ecolor='g', capsize=None )
    for i in xrange(args.N):
        ax.plot(T, Cv_singles[i],'k')
    ax.set_xlabel("T")
    ax.set_ylabel("Cv")
    ax = ax1
    ax.set_xlim([0,0.5])
    ax.errorbar(T, CvMom1, yerr=Cv_stdev,ecolor='g', capsize=None )
    for i in xrange(args.N):
        ax.plot(T, Cv_singles[i],'k')
    ax.set_xlabel("T")
    ax.set_ylabel("Cv")
    fig.savefig('cv_alpha_combined_K{K}_Nsub{N}_d{ndof}.pdf'.format(K = args.K,N=args.N,ndof=args.ndof))
    
    

