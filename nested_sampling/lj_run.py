#import pickle
import sys
import argparse
import numpy as np
#from types import *
from hparticle import HarParticle, HarRunner
from ising_model import IsingSystem, IsingRunner, IsingRunnerC
from nested_sampling_runner import run_nested_sampling

#import database_eigenvecs
from nested_sampling import NestedSampling, MonteCarloChain
from bh_sampling import get_thermodynamic_information, vector_random_uniform_hypersphere,\
    NestedSamplingBS


from pygmin.takestep import RandomDisplacement
#from pygmin.utils.xyz import write_xyz
from pygmin.accept_tests import SphericalContainer
from pygmin.systems import LJCluster
from pygmin.mindist import PointGroupOrderCluster
from pygmin.optimize import Result

class SphericalContainerNew(SphericalContainer):
    def __call__(self, energy=None, coords=None):
        return self.accept(coords)

class LJClusterNew(LJCluster):
    """same as LJCluster, but attach some additional information"""
    def __init__(self, natoms):
        super(LJClusterNew, self).__init__(natoms)
        self.nzero_modes = 6
        self.k = 3 * natoms - self.nzero_modes
        self.radius = 2.5 # specific to lj31
    
#    def get_metric_tensor(self):
#        return None
    
    def get_pgorder(self):
        return PointGroupOrderCluster(self.get_compare_exact())
    
    def get_config_tests(self):
        return [SphericalContainerNew(self.radius, nocenter=True)]
    
    def get_random_configuration(self):
        """make sure they're all inside the radius"""
        from pygmin.accept_tests import SphericalContainer
        test = self.get_config_tests()[0]
        coords = np.zeros([self.natoms,3])
        for i in range(self.natoms):
            coords[i,:] = vector_random_uniform_hypersphere(3) * self.radius
            assert(np.linalg.norm(coords[i,:]) <= self.radius)
        assert(test.accept(coords.flatten()))
        return coords.flatten()

    def center_coords(self, x):
        x = x.reshape(-1,3)
        natoms = x.shape[0] 
        com = np.sum(x, 0) / natoms
        x -= com[np.newaxis, :]

try:
    from src.runmc import mc_cython
except ImportError:
    print "warning, can't import compiled mc"

class MonteCarloCompiled(object):
    def __init__(self, radius):
        self.radius = radius
    
    def __call__(self, x0, mciter, stepsize, Emax, energy, seed=None):
        if seed is None:
            seed = np.random.randint(0, sys.maxint)
        x, energy, naccept = mc_cython(x0, mciter, stepsize, Emax, self.radius, seed)
#        print ret
        res = Result()
        res.x0 = x0
        res.x = x
        res.nsteps = mciter
        res.naccept = naccept
        res.energy = energy
        return res

def run_nested_sampling_lj(system, nreplicas=300, mciter=1000, label="test", 
                           minima=None, use_compiled=True, nproc=1,
                           triv_paral=True, minprob=1, maxiter=1e100, **kwargs):
    takestep = RandomDisplacement(stepsize=0.07)
    accept_tests = system.get_config_tests()

    if type(system) is  LJClusterNew:
        if use_compiled:
            mc_runner = MonteCarloCompiled(system.radius)
 #       import pickle
 #       with open("testpickle", "w") as pout:
 #           pickle.dump(mc_runner, pout)
        else:
            mc_runner = MonteCarloChain(system.get_potential(), takestep, accept_tests=accept_tests)
        print "using the compiled MC = ", use_compiled
    elif type(system) is HarParticle:
        mc_runner = HarRunner(system)
        print "using HarRunner"
    elif type(system) is IsingSystem:
        mc_runner = IsingRunnerC(system)
        print "using IsingRunner"
    else:
        raise TypeError('system type is not known')
       
        
    print "using", nproc, "processors"
    
    if minima is not None:
        assert(len(minima) > 0)
        print "using", len(minima), "minima"
        ns = NestedSamplingBS(system, nreplicas, mc_runner, minima, 
                              mciter=mciter, stepsize=0.07,
                              nproc=nproc, triv_paral=triv_paral, minprob=minprob, **kwargs)
    else:
        ns = NestedSampling(system, nreplicas, mc_runner, 
                            mciter=mciter, stepsize=0.07,
                            nproc=nproc, triv_paral=triv_paral, **kwargs)
    etol = 0.01

    ns = run_nested_sampling(ns, label, etol, maxiter=maxiter)
    return ns

def main():
    parser = argparse.ArgumentParser(description="do nested sampling with basin sampling for lennard jones clusters")
#    parser.add_argument("--db", type=str, nargs=1, help="database filename",
#                        default="otp.db")
    parser.add_argument("-d", "--db", type=str, help="database file name")
    parser.add_argument("-K", "--nreplicas", type=int, help="number of replicas", default=300)
    parser.add_argument("-n", "--mciter", type=int, help="number of iterations in the monte carlo chain", default=10000)
    parser.add_argument("-m", "--nminima", type=int, help="number of minima to use from the database", default=100)
    parser.add_argument("-A", "--nAtoms", type=int, help="number of atoms", default=31)
    parser.add_argument("-C", "--not-compiled-mc", action="store_false", help="option to use the Markov chain routine from C source (unique to LJ systems)", 
                        default=True)
    parser.add_argument("-P", "--nproc", type=int, help="number of precessors", default=1)
    parser.add_argument("-p", "--trivparal", action="store_true", help="set whether to do trivial parallelisation, by default True",default=False)
    parser.add_argument("-T", "--get-thermodynamic-properties", action="store_true", help="recalculates the eigenvectors of the hessian and writes them to the database",default=False)
    parser.add_argument("-a", "--minprob", type=float, help="probability of sampling from minima as a/K, default a=1",default=1)
    parser.add_argument("-S", "--system", type=int, help="define system type: 1 is LJ \n2 is HarParticle \n3 is Ising",default=1)
    args = parser.parse_args()

    natoms = args.nAtoms
    nreplicas = args.nreplicas
    mciter = args.mciter
    nminima = args.nminima
    
    if args.system == 1:
        system = LJClusterNew(natoms)
        label = "lj%d" % (natoms)
    elif args.system ==  2:
        ndim = 3
        centre = [0. for i in xrange(ndim)]
        kappa = [1. for i in xrange(ndim)]
        Eground = 0.
        Emax_init = 10.
        system = HarParticle(ndim, centre, kappa, Eground, Emax_init)
        label = "HarParticle_%dD" % (ndim)
    elif args.system == 3:
        system = IsingSystem()
        label = "Ising"
    else:
        raise TypeError('system type not known')
    nproc = args.nproc
    
    triv_paral = args.trivparal
    minprob = args.minprob
        
    #create minima database if needed
    if args.db is not None:
        dbname = args.db
        db = system.create_database(dbname, createdb=False)
        print dbname, "nminima", len(db.minima())
        minima = db.minima()
        if len(minima) > nminima:
            minima = minima[:nminima]
    else:
        minima = None

    # get thermodynamic information from database
    if (minima is not None) and (args.get_thermodynamic_properties is True):
        get_thermodynamic_information(system, db)
    # exit(1)

    # run nested sampling
    ns = run_nested_sampling_lj(system, nreplicas=nreplicas, 
                             label=label, minima=minima, mciter=mciter, 
                             use_compiled=args.compiled_mc, nproc=nproc, triv_paral = triv_paral, minprob = minprob)

    with open(label + ".energies", "w") as fout:
        fout.write( "\n".join([ str(e) for e in ns.max_energies]) )
        fout.write("\n")
            
#    if True:
#        import pylab as pl
#        e = np.array(ns.max_energies)
#        pl.plot(e)
#        pl.show()
    
    
    
if __name__ == "__main__":
    main()
