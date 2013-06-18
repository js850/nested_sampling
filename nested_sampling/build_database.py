import numpy as np

from lj_run import LJClusterNew
from pygmin.takestep import RandomDisplacement, AdaptiveStepsize

if __name__ == "__main__":
    natoms = 38
    system = LJClusterNew(natoms)
    
    ts = RandomDisplacement(stepsize=0.5)
    takestep = AdaptiveStepsize(ts)

    label = "lj%d" % (natoms)
    dbname = label + ".db"
    db = system.create_database(dbname)
    bh = system.get_basinhopping(database=db, takestep=ts, temperature=3., insert_rejected=True)
    bh.run(1000000)
        

    
    
