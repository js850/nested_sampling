from __future__ import division
import argparse
import numpy as np
import copy
from ..src.cv_trapezoidal import compute_alpha_cv_c
from itertools import chain

def run_alpha_variance(energies, nsubsets, K, Tmin, Tmax, nT, P, ndof, live):
    """
    returns the stdev calculated by jackknifing
    """
    Cv_sigma = _alpha_variance(energies, nsubsets, K, Tmin, Tmax, nT, P, ndof, live)
    return Cv_sigma()

class _alpha_variance(object):
    
    def __init__(self, energies, nsubsets, K, Tmin, Tmax, nT, P, ndof, live):
        self.E = np.array(energies)
        self.K = K
        self.N = len(self.E)
        self.nsubsets = nsubsets
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.nT = nT
        self.dT = (Tmax-Tmin) / nT
        self.T = np.array([Tmin + self.dT*i for i in range(nT)])
        self.P = P
        self.ndof = ndof
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
        return self.T, sigma, CvSingle, CvMom1
                
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
    
    def Cv_singles(self,alpha_sets):
        """
        returns the single Cvs
        """
        CvSingle = np.zeros((self.nsubsets,self.T.size))
        for i in xrange(self.nsubsets):
            T, CvSingle[i][:], U, U2 = compute_alpha_cv_c(self.E, np.array(alpha_sets[i][:]), float(self.P), self.K, float(self.Tmin), float(self.Tmax), self.nT, float(self.ndof), self.live)
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
    
    ####the following functions are currently unused, they were implemented in order to carry out jackknife resampling
    ####on sets of data with different alphas (the script currently performs alpha resampling, see )
    
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
            CvJack[i][:] = compute_alpha_cv_c(self.E , np.array(alphaJack[i][:]), float(P), self.K, float(self.Tmin), float(self.Tmax), self.nT, float(self.ndof), self.live)
        #print 'CvJack ',CvJack
        return np.array(CvJack)           