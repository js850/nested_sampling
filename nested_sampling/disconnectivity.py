from pele.landscape import Graph
from pele.utils.disconnectivity_graph import DisconnectivityGraph
import matplotlib.pyplot as plt
from lj_run import LJClusterNew
import sys

natoms = int(sys.argv[1])
dbname = sys.argv[2]
system = LJClusterNew(natoms)
database = system.create_database(dbname)

graph=Graph(database)
dg=DisconnectivityGraph(graph.graph)
dg.calculate()
dg.plot()
plt.savefig('discgraph_LJ{n}.pdf'.format(n = natoms))
plt.show()
