from lj_run import LJClusterNew
import sys
from pele.landscape import Graph

natoms = int(sys.argv[1])
dbname = sys.argv[2]

system = LJClusterNew(natoms)
db = system.create_database(dbname)

while True:
    min1 = db.minima()[0]

    graph = Graph(db)

    all_connected = True
    for m2 in db.minima()[1:]:
        if not graph.areConnected(min1, m2):
            all_connected = False
            break
    if all_connected:
        print "minima are all connected, ending"
        exit(1) 

    connect = system.get_double_ended_connect(min1, m2, db, fresh_connect=True, load_no_distances=True)
    connect.connect()


