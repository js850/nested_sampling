from __future__ import division
import argparse
import numpy as np
import copy
from ..src.cv_trapezoidal import compute_cv_c
from itertools import chain

def run_jackknife_variance(energies, nsubsets, K, Tmin, Tmax, nT, P, ndof, block, live):
    """
    returns the stdev calculated by jackknifing
    """
    Cv_sigma = _jackknife_variance(energies, nsubsets, K, Tmin, Tmax, nT, P, ndof, block, live)
    return Cv_sigma()

class _jackknife_variance(object):
    
    def __init__(self, energies, nsubsets, K, Tmin, Tmax, nT, P, ndof, block, live):
        self.E = np.array(energies)
        self.n = np.floor(float(K)/float(nsubsets))
        self.K = K
        self.nsubsets = nsubsets
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.nT = nT
        self.dT = (Tmax-Tmin) / nT
        self.T = np.array([Tmin + self.dT*i for i in range(nT)])
        self.P = P
        self.ndof = ndof
        self.block = block
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
        CvMom1 = self.jack_Cv_moments(CvJack)[0]
        sigma = self.jack_Cv_stdev(CvJack) 
        return self.T, sigma, CvSingle, CvMom1
    
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
        for i in xrange(self.nsubsets):
            T, CvJack[i][:], U, U2 = compute_cv_c(np.array(EJack[i][:]), float(self.P), (self.K - self.n), float(self.Tmin), float(self.Tmax), self.nT, float(self.ndof), self.live)
        #print 'CvJack ',CvJack
        return np.array(CvJack)
    
    def Cv_singles(self,Esplit):
        """
        returns the single Cvs
        """
        CvSingle = np.zeros((self.nsubsets,self.T.size))
        for i in xrange(self.nsubsets):
            T, CvSingle[i][:], U, U2 = compute_cv_c(np.array(Esplit[i][:]), float(self.P), (self.n), float(self.Tmin), float(self.Tmax), self.nT, float(self.ndof), self.live)
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