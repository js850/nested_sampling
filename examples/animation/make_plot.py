import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib as mpl
from matplotlib import mlab
import numpy as np
import time
import bisect

from nested_sampling import NestedSampling, MonteCarloWalker, Replica
from nested_sampling.utils.rotations import vector_random_uniform_hypersphere
from scipy.optimize import Result
from _bisect import bisect_left

    

class Pot(object):
    def get_energy(self, coords):
        assert len(coords) == 2
        return f(coords[::2], coords[1::2])[0]

def do_nested_sampling(nreplicas=10, niter=200, mciter=1000, stepsize=.8, estop=-.9,
                       x0=[1,1], r0=2):
    path = []
    def mc_record_position_event(coords=None, **kwargs):
        if len(path) == 0 or not np.all(path[-1] == coords):
            path.append(coords)

    p = Pot()
    print p.get_energy(np.array([1,2.]))
    mc_walker = MonteCarloWalker(p, mciter=mciter, events=[mc_record_position_event])
    
    # initialize the replicas with random positions
    replicas = []
    for i in xrange(nreplicas):
        # choose points uniformly in a circle 
        coords = vector_random_uniform_hypersphere(2) * r0 + x0
#         coords = np.random.uniform(-1,3,size=2)
        r = Replica(coords, p.get_energy(coords))
        replicas.append(r)
    
        
    ns = NestedSampling(replicas, mc_walker, stepsize=stepsize)
    results = [Result()]
    results[0].replicas = [r.copy() for r in replicas]
    for i in xrange(niter):
        ns.one_iteration()
        new_res = Result()
        new_res.replicas = [r.copy() for r in replicas]
        new_res.starting_replica = ns.starting_replicas[0].copy()
        new_res.new_replica = ns.new_replicas[0].copy()
        path.insert(0, new_res.starting_replica.x)
        new_res.mc_path = path
        results.append(new_res)
        path = []
        
        if ns.replicas[-1].energy < estop:
            break
        
        
        
#    plt.plot(ns.max_energies)
#    plt.show()
    
    return ns, results
    

def ig(x, y, x0, y0, sigx, sigy):
    """inverse gausian"""
    return np.exp(-( ((x-x0)/sigx)**2 + ((y-y0)/sigy)**2 ))

def f1(x, y):
    z = -1. *   ig(x, y, 0, 0, .3, .3)
    z += -0.6 * ig(x, y, 1.5, .5, 1, 1)
    z += -0.2 * ig(x, y, 1, 1, 3, 3)
    return z

def f2(x, y):
#    z =  -1.0 * ig(x, y, 0, -1, 1.3, 0.7)
#    z += -0.3 * ig(x, y, 1, 0.7, 0.45, 0.5)
#    z += -2.0 * ig(x, y, -1, 1, 2, 2)
    z =  -1.2 * ig(x, y, 0, -1, 1.3, 0.7)
    z += -1.5 * ig(x, y, 1, 0.7, 0.45, 0.5)
    z += -1.0 * ig(x, y, -1, 1, 2, 2)
    return z

def styblinski_tang(x, y):
    e1 = x**4 + 16*x**2 + 5*x
    e2 = y**4 + 16*y**2 + 5*y
    return (e1 + e2) / 2

f = f2


def plot3d():
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X = np.arange(-1, 3, 0.1)
    Y = np.arange(-1, 3, 0.1)
    X, Y = np.meshgrid(X, Y)
    Z = f(X, Y)
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=mpl.cm.coolwarm,
            linewidth=0, antialiased=False)
    ax.set_zlim(-1.01, 1.01)


    plt.show()

class NSViewer(object):
    def __init__(self, ns, results, xlim=[-1,3], ylim=[-1,3], better_dos=None):
        self.ns = ns
        self.results = results
        self.xmin = xlim[0]
        self.xmax = xlim[1]
        self.ymin = ylim[0]
        self.ymax = ylim[1]
        self.X = np.arange(self.xmin, self.xmax, 0.1)
        self.Y = np.arange(self.ymin, self.ymax, 0.1)
        self.X, self.Y = np.meshgrid(self.X, self.Y)
        self.Z = f(self.X, self.Y)
        self.zmin = np.min(self.Z)
        print "zmin", self.zmin
        self.zmax = 0
        self.vmax = self.Z.max() + .05
        
        self.Zlinear = np.array(sorted(self.Z.flatten(), reverse=True))
        self.Zlinear = self.Zlinear.reshape([-1,1])
        self.Zlinear_sorted = np.array(sorted(self.Zlinear.reshape(-1)))
        
        self.fig = plt.figure(figsize=(12,7))
        n = 10
        
        self.axes1 = plt.subplot2grid((n, 2*n+1), (0, 0), colspan=n, rowspan=n)
        self.axes2 = plt.subplot2grid((n, 2*n+1), (0, n), colspan=1, rowspan=n)
        self.axes3 = plt.subplot2grid((n, 2*n+1), (0, n+1), colspan=n, rowspan=n)
        
#         self.axes1 = self.fig.add_subplot(1,2,1)
#         self.axes2 = self.fig.add_subplot(1,2,2)
        
        self.cmap = mpl.cm.gist_heat
        self.cmap = mpl.cm.hot
        self.cmap = mpl.cm.gist_earth
#        self.cmap.set_clim(vmax=1.)

        
        max_energy_sidebar_indices = [self.sidebar_e_to_index(e)
                                           for e in self.ns.max_energies]
        self.sidebar_line_segments = [[(-0.5, y), (0.5, y)] 
                 for y in max_energy_sidebar_indices]

        # copute the density of states
        self.better_dos = better_dos
        self.dos = self.compute_dos(self.ns.max_energies)
        self.get_exact_dos()
        
    
    def plot_background(self):
        ax = self.axes1
        ax.clear()
        ax.set_xlim(self.xmin, self.xmax)
        ax.set_ylim(self.ymin, self.ymax)

        im = ax.imshow(self.Z, interpolation='bilinear', origin='lower',
                    cmap=self.cmap, extent=(self.xmin, self.xmax, self.ymin, self.ymax),
                    vmax=self.vmax, aspect="auto")
        
        self.plot_sidebar_background()
    
    def plot_contours_old(self, i):
        if i == 0: return
        ax = self.axes1
        ax.contour(self.X, self.Y, self.Z, levels=self.ns.max_energies[:i], 
                    vmin=self.zmin, vmax=self.vmax, 
                    colors="k", linestyles="solid")
        self.plot_sidebar_contours_old(i)

    def plot_contours_new(self, i):
        ax = self.axes1
        ax.contour(self.X, self.Y, self.Z, levels=self.ns.max_energies[i:(i+1)], 
                    vmin=self.zmin, vmax=self.vmax, 
                    colors="k", linestyles="solid")
        self.plot_sidebar_contours_new(i)


    def plot_replicas(self, i):
        # plot replicas
        ax = self.axes1
        results = self.results
        xrep = np.array([r.x[0] for r in results[i].replicas])
        yrep = np.array([r.x[1] for r in results[i].replicas])
        ax.scatter(xrep, yrep, c='k', facecolors="none")
        
        self.plot_sidebar_replicas(i)
        
    def plot_replica_to_delete(self, i):
        r = self.results[i].replicas[-1]
        self.axes1.scatter(r.x[0], r.x[1], marker='x', s=160, c='r')
        
            
    def plot_mc_path(self, i):
        ax = self.axes1
#            r = self.results[i].starting_replica
#            ax.scatter(r.x[0], r.x[1], c='b', 
#                       )
        
        r = self.results[i+1].new_replica
        ax.scatter(r.x[0], r.x[1], c='r')
        
        path = np.array(self.results[i+1].mc_path)
        ax.plot(path[:,0], path[:,1], '--k', lw=.5)
    
    def sidebar_e_to_index(self, energy):
        n = self.Zlinear.size
        return n-1-bisect.bisect_left(self.Zlinear_sorted, energy)
    
    def plot_sidebar_background(self):
        ax = self.axes2
        ax.clear()
        ax.set_yticks([])
        ax.set_xticks([])
        ax.imshow(self.Zlinear, aspect="auto", cmap=self.cmap, vmax=self.vmax)
        ax.set_ylim(self.Zlinear.size, 0)
        ax.set_xlim(-0.5,0.5)
    
    def plot_sidebar_replicas(self, i):
        ax = self.axes2
        ypos = [self.sidebar_e_to_index(r.energy)
                for r in self.results[i].replicas]
        xpos = [0] * len(ypos)
        ax.scatter(xpos, ypos, c='k', facecolors="none")
    
    def plot_sidebar_contours_old(self, i):
        # plot line segments
        lc = mpl.collections.LineCollection(self.sidebar_line_segments[:(i)],
                                            colors='k')
        self.axes2.add_collection(lc)

    def plot_sidebar_contours_new(self, i):
        # plot line segments
        lc = mpl.collections.LineCollection(self.sidebar_line_segments[i:(i+1)],
                                            colors='k')
        self.axes2.add_collection(lc)
    
    def get_exact_dos(self):
        V = len(self.Zlinear_sorted)
        K = float(len(self.ns.replicas))
        n = len(self.ns.max_energies)
        for k in xrange(n):
            vk = 1./(K+1) * 1./(K+1) * (K/(K+1))**k
            ek = self.Zlinear_sorted[int(vk*V)]
        self.better_dos = Result()
        self.better_dos.energies = [ self.Zlinear_sorted[int(V/(K+1) * 1./(K+1) * (K/(K+1))**i)] for i in xrange(n)]
        self.better_dos.dos = self.compute_dos(self.better_dos.energies)
        

    def compute_dos(self, energies):
        K = float(len(self.ns.replicas))
        dos = [1./(K+1) * (K/(K+1))**i 
                    for i, e in enumerate(energies)]
#        if self.better_dos is not None:
#            self.better_dos.dos = np.array(self.better_dos.dos)
#            self.better_dos.dos *= len(self.better_dos.energies) / len(self.dos)
        return dos
    
    def plot_dos(self, i):
        self.axes3.clear()
        if self.better_dos is not None:
            self.axes3.plot(self.better_dos.energies, np.log(self.better_dos.dos), '--r')
        self.axes3.plot(self.ns.max_energies[:(i+1)], np.log(self.dos[:(i+1)]))
        
    
    def pause(self):
        if not hasattr(self, "_pause_skip"):
            self._pause_skip = 0
        if self._pause_skip > 0:
            self._pause_skip -= 1
            plt.pause(.03)
            return
        plt.pause(.05)
        n = raw_input("press any key to continue. enter a number to skip future pauses")
        try:
            n = int(n)
            self._pause_skip = n
        except ValueError:
            return
    
    def title(self, str):
        if not hasattr(self, "_title"):
            self._title = self.fig.suptitle(str)
        self._title.set_text(str)
        
        
    
    def run(self):
        plt.ion()
        for i in xrange(len(self.results)-1):
            self.plot_background()
            self.plot_contours_old(i)
            if i == 0:
                self.title("the replicas are uniformly distributed in space")
            else:
                self.title("the replicas are uniformly distributed in the space with energy < Emax")
            self.plot_replicas(i)
            self.pause()
            self.plot_contours_new(i)
#            self.pause()
            self.title("the replica with the highest energy is removed.  It's energy becomes the next Emax")
            self.plot_replica_to_delete(i)
            self.pause()
            self.title("a new replica is generated via a MC walk")
            self.plot_mc_path(i)
            self.plot_dos(i)
            print "finishing iteration", i
            self.pause()
        plt.show(block=True)


def main1():
    ns, results = do_nested_sampling(nreplicas=15, niter=100, mciter=20)
    viewer = NSViewer(ns, results)
    viewer.run()
    
def main_st():
    r0 = 3
    xlim = ylim = [-r0,r0]
    x0 = [0,0]
    ns, results = do_nested_sampling(nreplicas=15, niter=100, mciter=20,
                                     x0=x0, r0=r0)
    viewer = NSViewer(ns, results, xlim=xlim, ylim=ylim)
    viewer.run()


    
    

def get_better_dos_2(nreplicas, niter, x0, r0, xlim, ylim, load_pickle=True):
    import pickle
    fname = "f2.better_dos.pkl"
    if load_pickle:
        try:
            better_dos = pickle.load(open(fname, "rb"))
            return better_dos
        except:
            pass
    ns, results = do_nested_sampling(nreplicas=nreplicas, niter=niter, mciter=20,
                                     x0=x0, r0=r0, estop=-1.7)
    viewer = NSViewer(ns, results, xlim=xlim, ylim=ylim)
    better_dos = Result()
    better_dos.dos = viewer.dos
    better_dos.energies = ns.max_energies
    pickle.dump(better_dos, open(fname, "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    return better_dos

def main2():
    xlim = [-4, 3]
    ylim = [-3, 4]
    x0 = [np.mean(xlim), np.mean(ylim)]
    r0 = np.abs(xlim[1] - xlim[0]) / 2
    nreplicas = 15
    niter = 100
    # do a longer NS run to get a better estimate
    better_dos = None
#    if True:
#        better_dos = get_better_dos_2(nreplicas*100, niter*100, x0, r0, xlim, ylim, load_pickle=False)
        
    ns, results = do_nested_sampling(nreplicas=nreplicas, niter=niter, mciter=20,
                                     x0=x0, r0=r0, estop=-1.7)
    viewer = NSViewer(ns, results, xlim=xlim, ylim=ylim, 
                      better_dos=better_dos)
    viewer.run()



if __name__ == "__main__":
    main2()
    
    