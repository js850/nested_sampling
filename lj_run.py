import pickle
import numpy as np

import database_eigenvecs
from nested_sampling import NestedSampling
from bh_sampling import sample_uniformly_in_basin,\
    populate_database, get_thermodynamic_information, vector_random_uniform_hypersphere,\
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
        self.radius = 2.5 # for lj31
    
    def get_metric_tensor(self):
        return None
    
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



def run_nested_sampling(system, nreplicas=300, mciter=1000, iterscale=300, label="test", database=None):
    takestep = RandomDisplacement(stepsize=0.5)
    accept_tests = system.get_config_tests()
    
    use_bs = True
    if use_bs:
        assert(len(database.minima()) > 0)
        ns = NestedSamplingBS(system, nreplicas, takestep, database, mciter=mciter, accept_tests=accept_tests)
    else:
        ns = NestedSampling(system, nreplicas, takestep, mciter=mciter, accept_tests=accept_tests)
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
            if True and i % 1000:
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

if __name__ == "__main__":
    natoms = 31
    nreplicas = 300
    mciter = 10000
    system = LJClusterNew(natoms)
    label = "lj%d" % (natoms)
    dbname = label + ".db"
    db = system.create_database(dbname)
    print dbname, "nminima", len(db.minima())
    
    # populate database
#    if False:
#        populate_database(system, db, niter=10000)
    
    # get thermodynamic information from database
    get_thermodynamic_information(system, db)
#    exit(1)

    # run nested sampling
    ns = run_nested_sampling(system, nreplicas=nreplicas, iterscale=1000000, label=label, database=db, mciter=mciter)

#    with open(label+".energies", "w") as fout:
#        fout.write( "\n".join([ str(e) for e in ns.max_energies]) )
#        fout.write("\n")
        
    
#    if True:
#        import pylab as pl
#        e = np.array(ns.max_energies)
#        pl.plot(e)
#        pl.show()
    
    
    
