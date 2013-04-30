import numpy as np
import networkx as nx
import sys
from pygmin.systems import BaseSystem
from pygmin.potentials import BasePotential

def make_grid_graph(dim, periodic=True):
    """
    this is a wrapper for nx.grid_graph()
    
    grid_graph creates a graph where the nodes are tuples (ix, iy, ...) 
    where ix and iy are the x and y positions of the site.  It would be more useful to
    have the nodes be simple integers that could act as indices for lists and arrays.
    The spatial positions will be returned in a separate dict
    """
    G = nx.grid_graph(dim, periodic=periodic)
    Gnew = nx.Graph()
    spatial = dict()
    node2i = dict()
    

    for i, node in enumerate(G.nodes()):
        Gnew.add_node(i)
        spatial[i] = node
        node2i[node] = i
    
    for edge in G.edges(data=False):
        u = node2i[edge[0]]
        v = node2i[edge[1]]
        Gnew.add_edge(u, v)
    
    return Gnew, spatial
        


class IsingPotential(BasePotential):
    """
    the potential object for an Ising system
    
    the networkx graph G defines the connectivity of the lattice
    
    This is for a 2d square lattice. Generalising to other lattices only involves
    defining a new G
    """
    def __init__(self, G):
        self.G = G
        self.eps = sys.float_info.epsilon
        
    def getEnergy(self, spins):
        E = 0.
        for edge in self.G.edges():
            i, j = edge
            E -= float(spins[i] * spins[j])
        return E

    def getEnergyGradient(self, spins):
        raise Exception("gradient not defined for Ising spin system")

    def getEnergyChange(self, spins, i):
        """
        return the change in energy if spin[i] were to flip
        """
        E = 0.
        Sold = float(spins[i])
        Snew = - Sold
        deltaS = Snew - Sold
        for j in self.G[i]:
            E -= deltaS * spins[j]
        return E
    
    def randomiseEnergy(self, E):
        u = np.random.ranf()
        l = 1. + (u - 0.5) * self.eps
        Erand = E * l 
        return Erand

class IsingSystem(BaseSystem):
    """
    The system class for an Ising System
     
    the networkx graph G defines the connectivity of the lattice
    
    This is for a 2d square lattice.  Generalising to other lattices only involves
    defining a new G
    
    Parameters
    ----------
    L : int
        number of spins along one axis of the square lattice (this is specific
        to a square lattice, but it could be generalised easily)
    """
    def __init__(self, L=10):
        self.G, self.spatial = make_grid_graph([L,L], periodic=True)

    def get_nspins(self):
        return self.G.number_of_nodes()

    def get_potential(self):
        return IsingPotential(self.G)
    
    def get_random_configuration(self):
        spins = np.random.randint(2, size=self.get_nspins())
        # spins are 0 or 1 now.  make them -1 or 1
        spins = spins * 2 - 1
        return spins
    
    def get_spatial_arrangement(self):
        xylist = self.spatial.items()
        # sort by index
        xylist.sort(key = lambda x: x[0])
        xylist = [xi[1] for xi in xylist]
        xylist = np.array(xylist).flatten()
        return xylist
    
    def get_config_tests(self):
        """
        needed to return true at all tests in run_nestedsampling
        """ 
        return True
    
#        spatial_to_spin = dict()
#        for xy, node in self.spatial.iteritems():
#            spatial_to_spin[xy] = spins[node]
#        return spatial_to_spin

class IsingRunner(object):
    def __init__(self, system):
        self.system = system
        self.pot = system.get_potential()
    
    def __call__(self, x0, mciter, stepsize, Emax, seed=None):
        return self.run(x0, mciter, Emax)
    
    def run(self, x0, mciter, Emax):        
        self.x = x0
        self.Emax = Emax
        self.mciter = mciter
        self.energy = self.pot.getEnergy(self.x)
        self.nsteps = 0
        self.naccept = 0
        for i in xrange(self.mciter):
            self.step()
        if self.naccept > 0:
            self.energy = self.pot.randomiseEnergy(self.energy)
        return self
    
    def step(self):
        """
        Do one iteration in the Monte Carlo chain
        """
        #flip one spin
        i = np.random.randint(self.system.get_nspins())
        dE = self.pot.getEnergyChange(self.x, i)
        new_energy = self.energy + self.pot.getEnergyChange(self.x,i)
        accept = new_energy < self.Emax
        
        if accept: 
            self.x[i] = - self.x[i]
            self.energy = new_energy
            self.naccept += 1
        
        self.nsteps += 1
    

def test():
    system = IsingSystem(L=10)
    pot = system.get_potential()
    spins = system.get_random_configuration()
    print spins
    print pot.getEnergy(spins)
    print pot.getEnergyChange(spins, 0)
    
    spins = spins.astype(float)
    xylist = system.get_spatial_arrangement()
    xylist = xylist.reshape(-1,2)
    x = xylist[:,0]
    y = xylist[:,1]
    import matplotlib.pyplot as plt
    plt.quiver(x, y, 0, spins, pivot="middle")
    a = plt.gca()
    a.set_xlim([-1, max(x)+1])
    a.set_ylim([-1, max(y)+1])

    plt.show() 

if __name__ == "__main__":
    test()