import numpy as np

from lj_run import LJClusterNew
from pygmin.takestep import RandomDisplacement, AdaptiveStepsize, Reseeding

if __name__ == "__main__":
    natoms = 75
    system = LJClusterNew(natoms)
    
    ts = RandomDisplacement(stepsize=0.42)
    takestep = AdaptiveStepsize(ts)
    takestep = Reseeding(takestep,system.get_random_configuration(),maxnoimprove=80)

    label = "lj%d" % (natoms)
    dbname = label + ".db"
    db = system.create_database(dbname)
    bh = system.get_basinhopping(database=db, takestep=ts, temperature=0.45, insert_rejected=True)
    bh.run(1000000)
        

    
    
