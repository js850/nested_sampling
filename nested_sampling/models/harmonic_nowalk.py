"""classes and functions related to a particle in a n dimensional harmonic potential"""

import numpy as np


from pele.systems import BaseSystem
from scipy.optimize import Result
from scipy import linalg

from harmonic import vector_random_uniform_hypersphere

   
def get_random_configuration_Emax(ndim, Emax, Eground=None, kappa_sqrt=None, x0=None):
    """return a random configuration for the harmonic well with energy less than Emax"""
    #radius is a scalar corresponding to the max distance from the centre
    E = Emax
    if Eground is not None:
        E -= Eground
    
    radius = np.sqrt(2. * Emax)
    x = vector_random_uniform_hypersphere(ndim) * radius
    
    if kappa_sqrt is not None:
        x /= kappa_sqrt
    
    if x0 is not None:
        x += x0

    return x
    

class HarmonicSampler(object):
    
    def __init__(self, potential, ndim):
        self.potential = potential
        self.ndim = ndim
    
    def __call__(self, x0, mciter, stepsize, Emax, energy, seed=None):
        return self.get_random_configuration(Emax)  
            
    def get_random_configuration(self, Emax):
        res = Result()
        res.mciter = 100
        res.nsteps = 100
        res.naccept = 50
        res.x = get_random_configuration_Emax(self.ndim, Emax)
        res.energy = self.potential.getEnergy(res.x)
        return res

#
# testing only below here
#

def test1():
    import matplotlib.pyplot as plt
    v = np.random.uniform(-4,4,size=10000)
    energies = v**2 / 2.
    plt.hist(energies, bins=50)
    plt.show()

def test2():
    import matplotlib.pyplot as plt
    from harmonic import Harmonic
    n = 100
    hp = Harmonic(n)#, centre=[0.]*n, kappa=[1.]*n)
    pot = hp
    hp = HarmonicSampler(pot, n)
    
    coords = np.array([get_random_configuration_Emax(n, 10.) for i in range(5000)])
    energies = [pot.getEnergy(x) for x in coords]
    print coords
    print "max energy", np.max(energies)
#    plt.plot(vals, 'o')
    plt.hist(coords[:,0], bins=50)
    plt.show()
    plt.clf()
    plt.hist(energies, bins=50)
    plt.show()



if __name__ == "__main__":
#    test()
    test2()