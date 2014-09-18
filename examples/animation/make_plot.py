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

    

class Pot(object):
    def get_energy(self, coords):
        assert len(coords) == 2
        return f(coords[::2], coords[1::2])[0]

def do_nested_sampling(nreplicas=10, niter=200, mciter=1000, stepsize=.8, estop=-.9,
                       x0=[1,1], r0=2,
                       xlim=None, ylim=None, circle=False
                       ):
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
        if circle: 
            coords = vector_random_uniform_hypersphere(2) * r0 + x0
        else:
            coords = np.zeros(2)
            coords[0] = np.random.uniform(xlim[0], xlim[1])
            coords[1] = np.random.uniform(ylim[0], ylim[1])
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

def prod(a):
    from operator import mul
    return reduce(mul, a)

class NSViewer(object):
    def __init__(self, ns, results, xlim=[-1,3], ylim=[-1,3], show_dos=False):
        self.ns = ns
        self.results = results
        self.show_dos = show_dos
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
        self._pause_count = 0
        self._save_figs = False
        
        self.Zlinear = np.array(sorted(self.Z.flatten(), reverse=True))
        self.Zlinear = self.Zlinear.reshape([-1,1])
        self.Zlinear_sorted = np.array(sorted(self.Zlinear.reshape(-1)))
        
        # set up the subplot geometry
        n = 10
        if self.show_dos:
            nx = 2*n
            self.fig = plt.figure(figsize=(12,7))
        else:
            nx = n  
            self.fig = plt.figure(figsize=(8,7))
        shape = (n+1, nx)
        self.axes1 = plt.subplot2grid(shape, (0, 0), colspan=n, rowspan=n)
        self.axes2 = plt.subplot2grid(shape, (n, 0), colspan=nx, rowspan=1)
        if self.show_dos:
            self.axes3 = plt.subplot2grid(shape, (0, n+1), colspan=n, rowspan=n)
            self.axes3.set_yticks([])

#        self.cmap = mpl.cm.gist_heat
#        self.cmap = mpl.cm.hot
        self.cmap = mpl.cm.gist_earth

        
        max_energy_sidebar_indices = [self.sidebar_e_to_index(e)
                                           for e in self.ns.max_energies]
        self.sidebar_line_segments = [[(y, -0.5), (y, 0.5)] 
                 for y in max_energy_sidebar_indices]

        # compute the density of states
        self.better_dos = None
        self.dos = self.compute_dos(len(self.ns.max_energies))
        self.get_exact_dos()
        
    
    def plot_background(self, sidebar=True):
        ax = self.axes1
        ax.clear()
        ax.set_xlim(self.xmin, self.xmax)
        ax.set_ylim(self.ymin, self.ymax)
        ax.set_xticks([])
        ax.set_yticks([])

        im = ax.imshow(self.Z, interpolation='bilinear', origin='lower',
                    cmap=self.cmap, extent=(self.xmin, self.xmax, self.ymin, self.ymax),
                    vmax=self.vmax, aspect="auto")
        
        if sidebar:
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


    def plot_replicas(self, i, replicas=None, sidebar=True):
        if replicas is None:
            replicas = self.results[i].replicas
        # plot replicas
        ax = self.axes1
        xrep = np.array([r.x[0] for r in replicas])
        yrep = np.array([r.x[1] for r in replicas])
        ax.scatter(xrep, yrep, c='k', facecolors="none")
        
        if sidebar:
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

        path = np.array(self.results[i+1].mc_path)
        ax.plot(path[:,0], path[:,1], '--k', lw=.5)
        
        self.plot_new_replica(i)

    def plot_new_replica(self, i):
        r = self.results[i+1].new_replica
        self.axes1.scatter(r.x[0], r.x[1], c='r', linewidths=0, s=30)
        
        # on the sidebar also
        xpos = 0
        ypos = self.sidebar_e_to_index(r.energy)
        self.axes2.scatter(xpos, ypos, c='r', linewidths=0, s=30)
    
    def sidebar_e_to_index(self, energy):
        n = self.Zlinear.size
        return n-1-bisect.bisect_left(self.Zlinear_sorted, energy)
    
    def plot_sidebar_background(self):
        ax = self.axes2
        ax.clear()
        ax.set_xticks([])
        ax.set_yticks([])
        ax.imshow(self.Zlinear.transpose(), aspect="auto", cmap=self.cmap, vmax=self.vmax)
        ax.set_xlim(self.Zlinear.size, 0)
        ax.set_ylim(-0.5,0.5)
    
    def plot_sidebar_replicas(self, i, replicas=None):
        if replicas is None:
            replicas = self.results[i].replicas
        ax = self.axes2
        ypos = [self.sidebar_e_to_index(r.energy)
                for r in replicas]
        xpos = [0] * len(ypos)
        ax.scatter(ypos, xpos, c='k', facecolors="none")
    
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
        N = len(self.Zlinear_sorted)
        V = N-1-self.sidebar_e_to_index(self.ns.max_energies[0])
        print self.ns.max_energies[0], self.Zlinear_sorted[V]
        K = float(len(self.ns.replicas))
        n = len(self.ns.max_energies)
        self.better_dos = Result()
        self.better_dos.energies = [ self.Zlinear_sorted[np.round(V * (K/(K+1))**i)] for i in xrange(n)]
        self.better_dos.dos = self.compute_dos(len(self.better_dos.energies))
        
        # also make some random dos versions
        self.better_dos.random_energies = []
        for i in xrange(20):
            alphas = np.random.beta(K,1, size=n-1)
            elist = [ self.Zlinear_sorted[np.round(V * prod(alphas[:i]))] for i in xrange(1,n)]
            elist.insert(0, self.Zlinear_sorted[np.round(V)])
            self.better_dos.random_energies.append(elist)
            
        

    def compute_dos(self, niter):
        K = float(len(self.ns.replicas))
        dos = [1./(K+1) * (K/(K+1))**i 
                    for i in xrange(niter)]
#        if self.better_dos is not None:
#            self.better_dos.dos = np.array(self.better_dos.dos)
#            self.better_dos.dos *= len(self.better_dos.energies) / len(self.dos)
        return dos
    
    def plot_dos(self, i):
        self.axes3.clear()
        self.axes3.set_yticks([])
        if self.better_dos is not None:
            self.axes3.plot(self.better_dos.energies, np.log(self.better_dos.dos), '-r', lw=.5)
            for elist in self.better_dos.random_energies:
                self.axes3.plot(elist, np.log(self.better_dos.dos), '-.b', lw=.5)
        self.axes3.plot(self.ns.max_energies[:(i+1)], np.log(self.dos[:(i+1)]), 'k')
        
    
    def pause(self):
        if self._save_figs or True:
#            plt.savefig("animation/animation_%i.pdf"%self._pause_count, type="pdf", bbox_inches="tight")
            plt.savefig("animation/animation_%i.jpg"%self._pause_count, type="jpg", bbox_inches="tight")
        self._pause_count += 1
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
            if self.show_dos:
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

#def make_first_plot(viewer):
#    p = Pot()
#    replicas = []
#    for i in xrange(viewer.ns.nreplicas):
#        # choose points uniformly in a circle 
##        coords = vector_random_uniform_hypersphere(2) * r0 + x0
#        coords = np.zeros(2)
#        coords[0] = np.random.uniform(viewer.xmin, viewer.xmax)
#        coords[1] = np.random.uniform(viewer.ymin, viewer.ymax)
#        r = Replica(coords, p.get_energy(coords))
#        replicas.append(r)
#
#    viewer.plot_background()
#    viewer.plot_sidebar_background()
#    viewer.plot_replicas(0, replicas=replicas)
#    viewer.plot_sidebar_replicas(0, replicas=replicas)
#    plt.savefig("animation_initial.pdf", type="pdf")
#    plt.show()

def make_first_plots(viewer):
    viewer.fig.set_figheight(5)
    viewer.fig.set_figwidth(6)
    i=0
    viewer.plot_background(sidebar=False)
    viewer.plot_replicas(i, sidebar=False)
#    viewer.title("the replicas are uniformly distributed in space")
    viewer.axes2.set_xticks([])
    viewer.axes2.set_yticks([])
    plt.savefig("initial_animation_1.pdf", type="pdf", bbox_inches="tight")
#    plt.show()
    viewer.plot_sidebar_background()
    viewer.plot_sidebar_replicas(i)
    plt.savefig("initial_animation_2.pdf", type="pdf", bbox_inches="tight")
    viewer.plot_contours_new(i)
    plt.savefig("initial_animation_3.pdf", type="pdf", bbox_inches="tight")
    viewer.plot_new_replica(i)
    plt.savefig("initial_animation_4.pdf", type="pdf", bbox_inches="tight")
    # add the monte carlo path
    viewer.plot_mc_path(i)
    plt.savefig("initial_animation_4_1.pdf", type="pdf", bbox_inches="tight")
    
    # refresh and redraw
    viewer.plot_background()
    viewer.plot_contours_old(i+1)
    viewer.plot_replicas(i+1)
    viewer.plot_contours_new(i+1)
    viewer.plot_new_replica(i)
    plt.savefig("initial_animation_5.pdf", type="pdf", bbox_inches="tight")
    
    # after n iterations
    i = 20
    viewer.plot_background()
    viewer.plot_contours_old(i)
    viewer.plot_replicas(i)
    viewer.plot_contours_new(i)
#    viewer.plot_new_replica(i)
    plt.savefig("initial_animation_6.pdf", type="pdf", bbox_inches="tight")
    
    # after n iterations with the mc path
    viewer.plot_mc_path(i)
    viewer.plot_replica_to_delete(i)
    plt.savefig("initial_animation_7.pdf", type="pdf", bbox_inches="tight")
    
    plt.show()
#    viewer.plot_contours_new(i)
#    viewer.title("the replica with the highest energy is removed.  It's energy becomes the next Emax")
#    viewer.plot_replica_to_delete(i)


def main2():
    xlim = [-4, 3]
    ylim = [-3, 4]
    x0 = [np.mean(xlim), np.mean(ylim)]
    r0 = np.abs(xlim[1] - xlim[0]) / 2
    nreplicas = 15
    niter = 100
        
    ns, results = do_nested_sampling(nreplicas=nreplicas, niter=niter, mciter=20,
                                     x0=x0, r0=r0, estop=-1.7, xlim=xlim, ylim=ylim,
                                     circle=False)
    viewer = NSViewer(ns, results, xlim=xlim, ylim=ylim, show_dos=False)
    
    if False:
        import pickle
        ns.mc_walker.events = [] # can't pickle it
        pickle.dump(viewer, open("viewer.pkl", "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    
    if False:
        make_first_plots(viewer)
        return
    
    viewer.run()



if __name__ == "__main__":
    main2()
    
    
