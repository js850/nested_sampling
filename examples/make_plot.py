import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib as mpl
import numpy as np
import time

from nested_sampling import NestedSampling, MonteCarloWalker, Replica
from scipy.optimize import Result

    

class Pot(object):
    def get_energy(self, coords):
        assert len(coords) == 2
        return f(coords[::2], coords[1::2])[0]

def do_nested_sampling(nreplicas=10, niter=200, mciter=1000, stepsize=.8, estop=-.7):
    path = []
    def mc_record_position_event(coords=None, **kwargs):
        path.append(coords) 

    p = Pot()
    print p.get_energy(np.array([1,2.]))
    mc_walker = MonteCarloWalker(p, mciter=mciter, events=[mc_record_position_event])
    replicas = []
    for i in xrange(nreplicas):
        coords = np.random.uniform(-1,3, size=2)
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
        new_res.mc_path = path
        path = []
        results.append(new_res)
        
        if ns.replicas[-1].energy < estop:
            break
        
        
        
#    plt.plot(ns.max_energies)
#    plt.show()
    
    return ns, results
    

def ig(x, y, x0, y0, a, b):
    """inverse gausian"""
    return - a * np.exp(-( (x-x0)**2 + (y-y0)**2 ) / b**2)

def f(x, y):
    x01 = np.array([0.,0])
    x02 = np.array([1.,1])
    z1 = ig(x, y, 0, 0, 1, .3)
    z2 = ig(x, y, 1, 1, .8, 1)
    return z1 + z2


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

def contour(ns, i):
    plt.clf()
#    fig = plt.figure()
#    ax = fig.gca(projection='3d')
    X = np.arange(-1, 3, 0.1)
    Y = np.arange(-1, 3, 0.1)
    X, Y = np.meshgrid(X, Y)
    Z = f(X, Y)
    zmin = np.min(Z)
    zmax = 0
    
    im = plt.imshow(Z, interpolation='bilinear', origin='lower',
                cmap=mpl.cm.summer, extent=(-1,3,-1,3))
    
    

    plt.contour(X,Y,Z, levels=ns.max_energies[:(i+1)], vmin=zmin, vmax=zmax, 
                colors="k", linestyles="solid")


def plot_replicas(results, i):
    # plot replicas
    xrep = np.array([r.x[0] for r in results[i].replicas])
    yrep = np.array([r.x[1] for r in results[i].replicas])
    plt.scatter(xrep, yrep, c='k')
    
    if i > 0:
        r = results[i].starting_replica
        plt.scatter(r.x[0], r.x[1], c='b')
    
        r = results[i].new_replica
        plt.scatter(r.x[0], r.x[1], c='r')
        
        path = np.array(results[i].mc_path)
        plt.plot(path[:,0], path[:,1], '--k', lw=.5)

def main():
    plt.ion()
    ns, results = do_nested_sampling(nreplicas=15, niter=100, mciter=50)
    for i in xrange(len(results)):
        contour(ns, i)
        plot_replicas(results, i)
        plt.pause(.2)
        raw_input("press any key to continue")
    plt.show(block=True)

if __name__ == "__main__":
    main()
    
    