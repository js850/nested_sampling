import pickle
import sys
import argparse
import numpy as np

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


class MonteCarloCompiled(object):
    def __init__(self, system, radius):
        self.radius = radius
        self.system = system
    
    def __call__(self, x0, mciter, stepsize, Emax):
        from src.runmc import mc_cython
        x, naccept = mc_cython(x0, mciter, stepsize, Emax, self.radius)
#        print ret
        self.x0 = x0
        self.x = x
        self.nsteps = mciter
        self.naccept = naccept
        pot = self.system.get_potential()
        self.energy = pot.getEnergy(x)
        return self



def run_nested_sampling(system, nreplicas=300, mciter=1000, iterscale=300, label="test", minima=None, use_compiled=True):
    takestep = RandomDisplacement(stepsize=0.07)
    accept_tests = system.get_config_tests()

    if use_compiled:
        mc_runner = MonteCarloCompiled(system, system.radius)
    else:
        mc_runner = None
    print "using the compiled MC = ", use_compiled
    
    use_bs = True
    
    if use_bs:
        assert minima is not None
        assert(len(minima) > 0)
        print "using", len(minima), "minima"
        ns = NestedSamplingBS(system, nreplicas, takestep, minima, mciter=mciter, accept_tests=accept_tests, mc_runner=mc_runner)
    else:
        ns = NestedSampling(system, nreplicas, takestep, mciter=mciter, accept_tests=accept_tests, mc_runner=mc_runner)
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
            if False and i % 10000:
                with open(label+".x", "w") as xout:
                    for r in ns.replicas:
                        write_xyz(xout, r.x)
                    
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
    args = parser.parse_args()

    natoms = args.nAtoms
    nreplicas = args.nreplicas
    mciter = args.mciter
    nminima = args.nminima
    system = LJClusterNew(natoms)
    
    label = "lj%d" % (natoms)
    if args.db is None:
        dbname = label + ".db"
        dbname = "lj31_small.db"
        print >> sys.stderr, "warning, using debug database"
    else:
        dbname = args.db
    db = system.create_database(dbname, createdb=False)
    print dbname, "nminima", len(db.minima())
    
    # populate database
    # if False:
    # populate_database(system, db, niter=10000)
    
    # get thermodynamic information from database
    get_thermodynamic_information(system, db)
    # exit(1)

    # run nested sampling
    minima = db.minima()
    if len(minima) > nminima:
        minima = minima[:nminima]
    ns = run_nested_sampling(system, nreplicas=nreplicas, iterscale=1000000, 
                             label=label, minima=minima, mciter=mciter, 
                             use_compiled=args.compiled_mc)

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
