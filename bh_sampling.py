"""
routines for using low energy minima found using basinhopping to 
improve sampling in "nested sampling" at low energies
"""
import numpy as np
from scipy.special import gamma, gammaln

from database_eigenvecs import HessianEigs

from pygmin.utils.rotations import vec_random_ndim
from pygmin.utils.hessian import sort_eigs, get_eig
from pygmin.thermodynamics import logproduct_freq2, normalmodes

def vector_random_uniform_hypersphere(k):
    """return a vector sampled uniformly in a hypersphere of dimension k"""
    u = vec_random_ndim(k)
    # draw the magnitude of the vector from a power law density function with power k-1
    p = np.random.power(k)
    return p * u

def sample_uniformly_in_basin_harmonic(m, Emax, k):
    """assuming the harmonic approximation return a configuration with energy less than Emax sampled uniformly from the basin defined by m
    
    this is exact in the harmonic approximation 
    """
    nm = m.hessian_eigs[0]
    evals = nm.eigenvalues
    vectors = nm.eigenvectors
    nzero = len(evals) - k

    # get uniform random k dimensional unit vector
    f = vec_random_ndim(k)

    # the target increase in energy is sampled from a power law distribution
    dE_target = np.random.power(k) * (Emax - m.energy)
    
    # scale f according to dE_target
    f *= np.sqrt(2. * dE_target) #TODO check prefactor

    # create the random displacement vector
    dx = np.zeros(m.coords.shape)
    for i in range(k):
        if evals[i+nzero] > 1e-4:
            dx += f[i] * vectors[:,i+nzero] / np.sqrt(evals[i+nzero])
    
    return m.coords + dx
    

def sample_uniformly_in_basin(m, Emax, potential, k):
    """return a configuration with energy less than Emax sampled uniformly from the basin defined by m
    """
    # displace randomly from the minimum according to the eigenvalues and eigenvectors
    E = Emax + 1.
    count = 0
    while E > Emax:
        coords = sample_uniformly_in_basin_harmonic(m, Emax, k)
        E = potential.getEnergy(coords)
        
        # print some stuff
        stepsize = np.linalg.norm(coords - m.coords)
        print "created structure with energy", E, "Emax", Emax, "Emin", m.energy, count, stepsize
        count += 1
        
    
    # now do a short monte carlo sampling to improve unbiased sampling
    
    return coords, E

def compute_log_phase_space_volume(m, Emax, k):
    """return the log (base e) of the phase space volume of minima m up to energy Emax
    
    V = Integral from m.energy to Emax of the harmonic density of states DoS(E)
    
        DoS(E) = (E - m.energy)^k / (Gamma(k) * prod_freq * O_k)
        
    k = number of vibrational degrees of freedom
    Gamma = gamma function
    prod_freq = the product of the frequencies (from the eigenvalues of the Hessian)
    O_k the order of the symmety point group

        V = (Emax - m.energy)**(k+1) / ((k+1) * np.gamma(k) * prod_freq * O_k)
    """
    from numpy import log
    if m.energy > Emax:
        raise ValueError("Emax (%g) must be greater than the energy of the minimum (%g)" % (Emax, m.energy))
    logV = (k+1) * log(Emax - m.energy) - log(k+1) - gammaln(k) - m.fvib - log(m.pgorder)
    return logV

def weighted_pick(weights):
    """sample uniformly from the objects with given weights
    
    return the selected index
    """
    if len(weights) == 0:
        raise ValueError("weights must not have zero length")
    r = np.random.uniform(0., sum(weights))
    s = 0.0
#    print r, len(weights)
    for i in range(len(weights)):
        s += weights[i]
        if r < s: return i
    return len(weights) - 1

def sample_from_database(system, db, Emax):
    # calculate the harmonic phase space volume of each minima and store it in list `weights`
    lweights = []
    minima = []
    for m in db.minima():
        if m.energy < Emax:
            lV = compute_log_phase_space_volume(m, Emax, system.k)
            lweights.append(lV)
            minima.append(m)
    lweights = np.array(lweights)
    weights = np.exp(lweights - np.max(lweights))
    
    # select a minimum uniformly given `weights`
    print "weights", weights[:10]
    index = weighted_pick(weights)
#    print index, len(weights), len(minima)
    m = minima[index]
    
    # sample configuration uniformly from the basin of minima m 
    coords, E = sample_uniformly_in_basin(m, Emax, system.get_potential(), system.k)

    return coords, E
 
def populate_database(system, db, niter=1000):
    """return a database with all important low energy minima 
    
    """
    # use basinhopping to find the low energy minima and store them in a database
    bh = system.get_basinhopping(database=db)
    bh.run(niter)

def get_thermodynamic_information(system, db):
    """
    for each minima in database, get all information necessary to compute the density of states
    """
    # get the point group information
    print "getting the point group information"
    determine_pgorder = system.get_pgorder()
    for m in db.minima():
        m.pgorder = determine_pgorder(m.coords)
#        print m.pgorder
    db.session.commit()

    # get the frequencies
    print "getting the normal mode frequencies"
    for m in db.minima():
        # calculate the Hessian
        pot = system.get_potential()
        e, g, hess = pot.getEnergyGradientHessian(m.coords)
        
        # calculate the normal modes from the hessian
        freq, evec = normalmodes(hess, metric=system.get_metric_tensor())
        # calculate the log product of the positive normal mode frequencies
        n, lnf = logproduct_freq2(freq, system.nzero_modes)
        m.fvib = lnf
        
        # calculate the eigenvalues and eigenvectors of the hessian and attach them to the database
        eval, evec = get_eig(hess)
        if len(m.hessian_eigs) == 0: 
            nm = HessianEigs(m, eval, evec)

    db.session.commit()
    
    
    # combine the point group information and frequencies to compute the density of states
    
    return db
    
if __name__ == "__main__":
    # define the system
    from lj_run import LJClusterNew
    natoms = 13
    system = LJClusterNew(natoms)

    db = system.create_database("lj%d.db" % (natoms))
    if False:
        db = populate_database(system, db, niter=500)
    
    Emin = db.minima()[0].energy
    Emax = Emin + 1.

    k = system.k
    for i in range(10):
        coords, E = sample_from_database(system, db, Emax)
    
