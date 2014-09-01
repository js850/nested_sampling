import numpy as np
import matplotlib as mpl
from matplotlib import pylab as plt

def plots1d(nimages=30, nreplicas=4, with_hist=True, show=True,
            height=5, width=8):
    plt.clf()
    fig = plt.gcf()
    fig.set_figheight(height)
    fig.set_figwidth(width)
    if with_hist:
        ncol_hist = 2
    else:
        ncol_hist = 0
    ax_list = [plt.subplot2grid((1,nimages+ncol_hist), (0,i), 
                               colspan=1, rowspan=nimages)
               for i in xrange(nimages)]
    
    for ax in ax_list:
        ax.set_ylim(0,1)
        ax.set_xticks([])
        ax.set_yticks([])
        ypos = np.random.uniform(0,1,nreplicas)
        ymax= ypos.max()
        ax.axhspan(0, ymax, alpha=0.2)
        xpos = np.zeros(nreplicas, )
        ax.scatter(xpos, ypos, c='k', facecolors="none")
        ax.scatter(0, ymax, c='r', linewidths=0, s=40)
    
    if with_hist:
        ax = plt.subplot2grid((1,nimages+ncol_hist), (0,nimages), 
                                   colspan=ncol_hist, rowspan=nimages)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylim(0,1)

            
        n = 100000
        rmax = np.random.beta(nreplicas, 1, size=n)
        ax.hist(rmax, bins=np.sqrt(n)/10, orientation='horizontal', normed=True)
        
        if False:
            y = np.arange(1,0,-.01)
            x = y**(nreplicas-1)
            x /= x.max()
            ax.plot(x,y)
            ax.relim()
        
            
    
    
    if show:
        plt.show()
    

def make_1d_plots():
    show=False    
    plots1d(nimages=1, nreplicas=1, with_hist=False, show=show, width=1)
    plt.savefig("plot1d_1_1_nohist.pdf", type="pdf", bbox_inches="tight")
   
    plots1d(nimages=1, nreplicas=2, with_hist=False, show=show, width=1)
    plt.savefig("plot1d_1_2_nohist.pdf", type="pdf", bbox_inches="tight")

    plots1d(nimages=10, nreplicas=2, with_hist=True, show=show)
    plt.savefig("plot1d_10_2.pdf", type="pdf", bbox_inches="tight")
    
    plots1d(nimages=10, nreplicas=10, with_hist=True, show=show)
    plt.savefig("plot1d_10_10.pdf", type="pdf", bbox_inches="tight")
    plt.show()

if __name__ == "__main__":
    make_1d_plots()
#    plots1d(nimages=10, nreplicas=10)
