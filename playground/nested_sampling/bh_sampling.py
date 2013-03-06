"""
routines for using low energy minima found using basinhopping to 
improve sampling in "nested sampling" at low energies
"""
import numpy as np
from scipy.special import gamma, gammaln

from pygmin.systems import LJCluster
from pygmin.utils.rotations import vec_random_ndim

def vector_random_uniform_hypersphere(k):
    """return a vector sampled uniformly in a hypersphere of dimension k"""
    u = vec_random_ndim(k)
    # draw the magnitude of the vector from a power law distribution with power k-1
    p = np.random.power(k)
    return p * u
    

def sample_uniformly_in_basin_harmonic(m, Emax):
    """assuming the harmonic approximation return a configuration with energy less than Emax sampled uniformly from the basin defined by m
    
    this is exact in the harmonic approximation 
    """
    dx = np.zeros(m.coords.shape)

    k = None #TODO number of vibrational degrees of freedom = 3*N - 6
    eigs = m.eval_vecs #TODO

    # get uniform random vector in k dimensional hypersphere
    f = vector_random_uniform_hypersphere(k)
    i = 0
    for eval, evec in eigs:
        if eval > 1e-4:
            dx += Emax * f[i] * evec / eval #TODO check prefactor
        i += 1
    
    return m.coords + dx
    

def sample_uniformly_in_basin(m, Emax):
    """return a configuration with energy less than Emax sampled uniformly from the basin defined by m
    """
    # displace randomly from the minimum according to the eigenvalues and eigenvectors
    
    # now do a short monte carlo sampling to improve unbiased sampling
    
    raise NotImplementedError

def compute_log_phase_space_volume(m, Emax):
    """return the log (base e) of the phase space volume of minima m up to energy Emax
    
    V = Integral from m.energy to Emax of the harmonic density of states DoS(E)
    
        DoS(E) = (E - m.energy)^k / (Gamma(k) * prod_freq * O_k)
        
    k = number of vibrational degrees of freedom
    Gamma = gamma function
    prod_freq = the product of the frequencies (from the eigenvalues of the Hessian)
    O_k the order of the symmety point group

        V = (Emax - m.energy)**(k+1) / ((k+1) * np.gamma(k) * prod_freq * O_k)
    """
    # I probably want to do this in log space
    from numpy import log
    k = None
    prod_freq = None
    O_k = None
    logV = (k+1) * log(Emax - m.energy) - log(k+1) - gammaln(k) - log(prod_freq) - log(O_k)
    raise NotImplementedError
    return logV

def select_uniformly_weighted(weights):
    """sample uniformly from the objects with given weights
    
    return the selected index
    """
    raise NotImplementedError

def sample_from_database(db, Emax):
    # calculate the harmonic phase space volume of each minima and store it in list `weights`
    weights = []
    minima = db.minima()
    for m in minima:
        lV = compute_log_phase_space_volume(m, Emax)
        weights.append(np.exp(lV))
    
    # select minima uniformly given `weights`
    index = select_uniformly_weighted(weights)
    m = minima[index]
    
    # sample configuration uniformly from the basin of minima m 
    sample_uniformly_in_basin(m, Emax)
    pass
    
def generate_database(natoms=31):
    """return a database with all important low energy minima and all information
    necessary to compute the density of states
    """
    # define the system
    system = LJCluster(natoms)
    
    # use basinhopping to find the low energy minima and store them in a database
    db = system.create_database("lj31.db")
    bh = system.get_basinhopping(database=db)
    nsteps = 1000
    bh.run(nsteps)
    
    # get the point group information

    # get the frequencies
    
    # combine the point group information and frequencies to compute the density of states
    
    return db
    
if __name__ == "__main__":
    generate_database()