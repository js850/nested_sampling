from pygmin.gui import run_gui
from lj_run import LJClusterNew

natoms = 31
system = LJClusterNew(natoms)
label = "lj%d" % (natoms)
dbname = label + ".db"
run_gui(system, dbname)
#db = system.create_database(dbname)
