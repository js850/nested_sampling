import pickle

import database_eigenvecs
from nested_sampling import NestedSampling
from bh_sampling import LJClusterNew, sample_uniformly_in_basin,\
    populate_database, get_thermodynamic_information

from bh_sampling import LJClusterNew
from pygmin.takestep import RandomDisplacement
from pygmin.utils.xyz import write_xyz

def run_nested_sampling(system, nreplicas=300, mciter=1000, iterscale=300, label="test"):
    takestep = RandomDisplacement(stepsize=0.5)
    accept_tests = system.get_config_tests()
    
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
    nreplicas = 500
    mciter = 100000
    system = LJClusterNew(natoms)
    label = "lj%d" % (natoms)
    dbname = label + ".db"
    db = system.create_database()
    
    # populate database
    if False:
        populate_database(system, db, niter=10000, mciter=mciter)
    
    # get thermodynamic information from database
    get_thermodynamic_information(system, db)
#    exit(1)

    # run nested sampling
    ns = run_nested_sampling(system, nreplicas=nreplicas, iterscale=1000000, label=label)

#    with open(label+".energies", "w") as fout:
#        fout.write( "\n".join([ str(e) for e in ns.max_energies]) )
#        fout.write("\n")
        
    
#    if True:
#        import pylab as pl
#        e = np.array(ns.max_energies)
#        pl.plot(e)
#        pl.show()
    
    
    
