"""classes and functions related to a particle in a n dimensional harmonic potential"""

import numpy as np

from pele.potentials import BasePotential
from pele.systems import BaseSystem
from pele.optimize import Result
from scipy import linalg
from bh_sampling import vector_random_uniform_hypersphere

__all__ =["HarPotential","HarParticle","HarRunner"]

class HarPotential(BasePotential):
    """
    Defines an harmonic potential, returns energy depending on coordinates
    
    Parameters
    ----------
    centre: list of int
        centre of harmonic well
    kappa: nxn matrix of floats
        spring constants
    
    """
    def __init__(self, centre, kappa):
        self.centre = np.asfarray(centre)
        self.kappa = kappa #np.asfarray(kappa) #kappa must be a nxn (diagonal) tensor where n is the number of dimension
                
    def getEnergy(self, coords):
        coords = np.asfarray(coords)
        dist_vec = coords - self.centre
        F = np.dot(self.kappa, dist_vec)
        return 0.5 * np.dot(F,F)
        #return 0.5 * dist_vec.dot(self.kappa * dist_vec) STEFANO 29 JUNE
#        E_vec = 0.5 * self.kappa * dist_vec * dist_vec      
#        E = 0.
#        for e in E_vec:
#            E += e
##        print "getEnergy ", E
#        return E
   
class HarParticle(BaseSystem):
    """ Defines a particle in an harmonic potential of n dimensions 
    
    Notes
    -----
    with n degrees of freedom and spring constant k:
     
        Z = (2*pi*T/k)^(n/2)
        <E> = -n*T/2
        <E^2> = T^2*n*(n+2)/4
        Cv = n/2
    """
    
    def __init__(self, ndim, centre=None, kappa=None, Eground=0., Emax_init=10.):
        self.ndim = ndim
        #kappa is an n x n dimensional tensor  containing n spring constants
        if kappa is None:
            self.kappa = np.diag(np.ones(self.ndim)) #np.ones(self.ndim) STEFANO 29 JUNE, create diagonal matrix of ones
        else:
            self.kappa = kappa #np.asfarray(kappa)
        assert(np.sqrt(self.kappa.size) == ndim) #STEFANO 29 JUNE
         # use numpy vectors for coordinates, centre must be a numpy vector
        if centre is None:
            self.centre = np.zeros(self.ndim)
        else:
            self.centre = np.array(centre)
        self.Eground = Eground
        self.Emax_init = Emax_init
        self.kappa_sqrt = linalg.sqrtm(self.kappa) #np.sqrt(self.kappa)
        
        print "Emax_init", self.Emax_init
            
    def get_potential(self):
        return HarPotential(self.centre, self.kappa) 
    
    def vector_random_uniform_hypersphere(self):
        """
        return a normalised vector sampled uniformly in a hypersphere of dimension ndim
        must then multiply each component by sqrt(norm_max)
        """
        return vector_random_uniform_hypersphere(self.ndim)
    
    def get_config_tests(self):
        """
        needed to return true at all tests in run_nestedsampling
        """ 
        return True
    
    def get_config_tests_in(self, coords, radius):
        coords_norm = np.linalg.norm(coords-self.centre)
#        print "coords_norm", coords_norm
        return coords is not None and radius > coords_norm
     
    def get_random_configuration_Emax(self, Emax):
        """make sure they're all inside the radius, get_config_test is not strictly necessary, consider removing it"""
        #radius is a scalar corresponding to the max distance from the centre
        radius = np.sqrt(2. * (float(Emax) - self.Eground))

        coords = self.vector_random_uniform_hypersphere() * radius
        coords = np.dot(coords,linalg.inv(self.kappa_sqrt)) #/= self.kappa_sqrt
        coords += self.centre
        assert(self.get_config_tests_in(coords, radius))
        return coords
    
    def get_random_configuration(self):
        """make sure they're all inside the radius, get_config_test is not strictly necessary, consider removing it"""
        #radius is a scalar corresponding to the max distance from the centre
        Emax = self.Emax_init
        return self.get_random_configuration_Emax(Emax)

class HarRunner(object):
    
    def __init__(self, system):
        self.system = system
        self.pot = system.get_potential()
    
    def __call__(self, x0, mciter, stepsize, Emax, energy, seed=None):
        return self.run(Emax)  
            
    def run(self, Emax):
        self.Emax = Emax
        res = Result()
        res.mciter = 100
        res.nsteps = 100
        res.naccept = 70
        res.x = self.system.get_random_configuration_Emax(self.Emax)
        res.energy = self.pot.getEnergy(res.x)
        return res

def test1():
    import matplotlib.pyplot as plt
    v = np.random.uniform(-4,4,size=10000)
    energies = v**2 / 2.
    plt.hist(energies, bins=50)
    plt.show()

def test():
    import matplotlib.pyplot as plt
    n = 100
    hp = HarParticle(n)#, centre=[0.]*n, kappa=[1.]*n)
    pot = hp.get_potential()
    
    coords = np.array([hp.get_random_configuration_Emax(10.) for i in range(10000)])
    energies = [pot.getEnergy(x) for x in coords]
    print coords
    print "max energy", np.max(energies)
#    plt.plot(vals, 'o')
    plt.hist(coords[:,0], bins=50)
    plt.show()
    plt.hist(energies, bins=50)
    plt.show()
    
if __name__ == "__main__":
    test()
#    test1()