import numpy as np

from lj_run import LJClusterNew
from pele.takestep import RandomDisplacement, AdaptiveStepsize

if __name__ == "__main__":
    natoms = 38
    system = LJClusterNew(natoms)
    
    ts = RandomDisplacement(stepsize=0.5)
    takestep = AdaptiveStepsize(ts)
    takestep = Reseeding(takestep,system.get_random_configuration(),maxnoimprove=80)

    label = "lj%d" % (natoms)
    dbname = label + ".db"
    db = system.create_database(dbname)
    bh = system.get_basinhopping(database=db, takestep=ts, temperature=0.5, insert_rejected=True)
    bh.run(1000000)
        

    
    
