from pele.gui import run_gui
from lj_run import LJClusterNew
import sys

natoms = int(sys.argv[1])
dbname = sys.argv[2]

#natoms = 31
system = LJClusterNew(natoms)
#label = "lj%d" % (natoms)
#dbname = label + ".db"
run_gui(system, dbname)
#db = system.create_database(dbname)
