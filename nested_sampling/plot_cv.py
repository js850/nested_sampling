import argparse
import numpy as np
import copy
from compute_cv import compute_Z, get_energies
from src.cv_trapezoidal import  compute_cv_c
from itertools import chain, cycle
import matplotlib
from matplotlib import rc
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

#######################SET LATEX OPTIONS###################
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern']})
#rc('text.latex',preamble=r'\usepackage{times}')
plt.rcParams.update({'font.size': 12})
plt.rcParams['xtick.major.pad'] = 8
plt.rcParams['ytick.major.pad'] = 8
##########################################################
#####################LINE STYLE CYCLER####################
lines = ["-","--","-.",":"]
linecycler = cycle(lines)
##########################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="load energy intervals and compute Cv stdev", 
                                     epilog="must write file name followed by a label string, otherwise raise the flag --nolabels")
    parser.add_argument("fname", nargs="+", type=str, help="filenames with heat capacity followed by label")
    parser.add_argument("--nolabels", action="store_true", help="turn on if do not want to insert labels",default=False)
    parser.add_argument("--show", action="store_true", help="show explorable plot preview",default=False)
    parser.add_argument("--leg_loc", type=int, help="define location of the legend (default upper right 1):"
                        " \"upper right=1\" \"upper left=2\" \"lower left=3\" \"lower right=4\" \"right=5\" "
                        "\"center left=6\" \"center right=7\" \"lower center=8\" \"upper center=9\" \"center=10\"",default=1)
    parser.add_argument("--xlabel", type=str, help="set x-label",default="T")
    parser.add_argument("--ylabel", type=str, help="set y-label",default="Cv")
    parser.add_argument("--title", type=str, help="set title",default=None)
    parser.add_argument("--linewidth", type=float, help="set line width",default=2.5)
    parser.add_argument("--grid", action="store_true", help="plot grid",default=False)
    
    args = parser.parse_args()
    print args.fname
    fname = args.fname
    nolabels = args.nolabels
    xlabel = args.xlabel
    ylabel = args.ylabel
    show = args.show
    title = args.title
    grid = args.grid
    linew = args.linewidth
    ####################################################################################################
    #deal with input
        
    if nolabels is False:
        all_data = [[] for i in xrange(len(fname)/2)]
        all_labels = np.array(["no_name" for i in xrange(len(fname)/2)])
        j = 0    
        for name,i in zip(fname,xrange(len(fname))):
            if ((i+1) % 2) is not 0:
                data = np.genfromtxt(name)
                all_data[j] = copy.deepcopy(data.tolist())
            else:
                all_labels[j] = name
                j += 1
        print all_labels
    else:
        all_data = [[] for i in xrange(len(fname))]
        all_labels = np.array(["no_name" for i in xrange(len(fname))])
        for name,i in zip(fname,xrange(len(fname))):
            data = np.genfromtxt(name)
            all_data[i] = copy.deepcopy(data.tolist())
    
    all_data = np.array(all_data)
    ######################################################################################################
    
    #####################
    ####### PLOT ########
    #####################
    
    ####SET COLOUR MAP######
    cm = plt.get_cmap('jet')
    ########################
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_color_cycle([cm(1.*i/len(all_data)) for i in xrange(len(all_data))])
    
    for data, label, i in zip(all_data, all_labels, xrange(len(all_data))): 
            color = cm(1.*i/5)
            ax.plot(data[:,0], data[:,1], next(linecycler), label = label, linewidth=linew)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if grid is True:
        ax.grid(True,linestyle='--',color='0.75')
    if title is not None:
        ax.set_title(title)
    if nolabels is False:
        ax.legend(frameon=False)
    if show is True:
         plt.show()
    
    fig.savefig('cv_plot.eps')
    