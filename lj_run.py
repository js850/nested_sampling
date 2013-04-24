import pickle
import sys
import argparse
import numpy as np
from types import *
from hparticle import *

import database_eigenvecs
from nested_sampling import NestedSampling
from bh_sampling import sample_uniformly_in_basin,\
    get_thermodynamic_information, vector_random_uniform_hypersphere,\
    NestedSamplingBS


from pygmin.takestep import RandomDisplacement
from pygmin.utils.xyz import write_xyz
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


from src.runmc import mc_cython
class MonteCarloCompiled(object):
    def __init__(self, radius):
        self.radius = radius
    
    def __call__(self, x0, mciter, stepsize, Emax, seed=None):
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

def run_nested_sampling(system, nreplicas=300, mciter=1000, iterscale=300, label="test", minima=None, use_compiled=True, use_bs=False, nproc = 1, triv_paral = True, minprob = 1):
    takestep = RandomDisplacement(stepsize=0.07)
    accept_tests = system.get_config_tests()

    if type(system) is  LJClusterNew:
        if use_compiled:
            mc_runner = MonteCarloCompiled(system.radius)
 #       import pickle
 #       with open("testpickle", "w") as pout:
 #           pickle.dump(mc_runner, pout)
        else:
            mc_runner = None
        print "using the compiled MC = ", use_compiled
    elif type(system) is HarParticle:
        mc_runner = HarRunner(system)
        print "using HarRunner"
    else:
        raise TypeError('system type is not known')
       
        
    print "using", nproc, "processors"
    
    if use_bs:
        assert minima is not None
        assert(len(minima) > 0)
        print "using", len(minima), "minima"
        ns = NestedSamplingBS(system, nreplicas, takestep, minima, mciter=mciter, accept_tests=accept_tests, mc_runner=mc_runner, nproc = nproc, triv_paral = triv_paral, minprob = minprob)
    else:
        ns = NestedSampling(system, nreplicas, takestep, mciter=mciter, accept_tests=accept_tests, mc_runner=mc_runner, nproc = nproc, triv_paral = triv_paral)
    etol = 0.01
    isave = 0
    maxiter = nreplicas * iterscale
    i = 0
#    pkl = label + ".ns.pickle"
    with open(label+".energies", "w") as fout:
        while True:
            ediff = ns.replicas[-1].energy - ns.replicas[0].energy
            if i != 0 and i % 100 == 0:
                fout.write( "\n".join([ str(e) for e in ns.max_energies[isave:i]]) )
                fout.write("\n")
                fout.flush()
                isave = i

            if False and i % 10000 == 0:
                with open(label+".x", "w") as xout:
                    for r in ns.replicas:
                        write_xyz(xout, r.x)
                    
            if i % 1000 == 0:
                if i == 0: openas = "w"
                else: openas = "a"
                with open(label+".replicas", openas) as xout:
                    xout.write("#energy niter from_random\n")
                    for r in ns.replicas:
                        xout.write("%g %d %d\n" % (r.energy, r.niter, int(r.from_random))) 
                    xout.write("\n\n")

            if ediff < etol: break
            if i >= maxiter: break  
            ns.one_iteration()
            i += 1
        fout.write( "\n".join([ str(e) for e in ns.max_energies[isave:i]]) )
        fout.write("\n")
        fout.flush()
#    print ns.max_energies
    print "min replica energy", ns.replicas[0].energy
    print "max replica energy", ns.replicas[-1].energy
    ns.finish()
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
    parser.add_argument("-C", "--compiled-mc", type=bool, help="option to use the Markov chain routine from C source (unique to LJ systems)", 
                        default=True)
    parser.add_argument("-P", "--nproc", type=int, help="number of precessors", default=1)
    parser.add_argument("-p", "--trivparal", type=bool, help="set whether to do trivial parallelisation, by default True",default=True)
    parser.add_argument("-T", "--get-thermodynamic-properties", type=bool, help="recalculates the eigenvectors of the hessian and writes them to the database",default=False)
    parser.add_argument("-a", "--minprob", type=bool, help="probability of sampling from minima as a/K, default a=1",default=1)
    parser.add_argument("-S", "--system", type=int, help="define system type: 1 is LJ \n2 is HarParticle",default=1)
    args = parser.parse_args()

    natoms = args.nAtoms
    nreplicas = args.nreplicas
    mciter = args.mciter
    nminima = args.nminima
    
    if args.system == 1:
        system = LJClusterNew(natoms)
        label = "lj%d" % (natoms)
    elif args.system ==  2:
        ndim = 1
        centre = [0.]
        kappa = [1.]
        Eground = 0.
        Emax_init = 100.
        system = HarParticle(ndim, centre, kappa, Eground, Emax_init)
        label = "HarParticle_%dD" % (ndim)
    else:
        raise TypeError('system type not known')
    nproc = args.nproc
    
    triv_paral = args.trivparal
    minprob = args.minprob
        
    #create minima database if needed
    if args.db is not None:
        dbname = args.db
        db = system.create_database(dbname, createdb=False)
        use_bs = True
        print dbname, "nminima", len(db.minima())
        minima = db.minima()
        if len(minima) > nminima:
            minima = minima[:nminima]
    else:
        use_bs = False #can avoid writing this and do it more concisely?
        minima = None

    # get thermodynamic information from database
    if (use_bs is True) and (args.get_thermodynamic_properties is True):
        get_thermodynamic_information(system, db)
    # exit(1)

    # run nested sampling
    ns = run_nested_sampling(system, nreplicas=nreplicas, iterscale=1000000, 
                             label=label, minima=minima, mciter=mciter, 
                             use_compiled=args.compiled_mc, use_bs=use_bs, nproc=nproc, triv_paral = triv_paral, minprob = minprob)

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
