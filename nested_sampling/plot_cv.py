#    plot_xye.py
#    
#    Created by Stefano Martiniani in June 2013, stefano.martiniani@gmail.com
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import numpy as np
import copy
from itertools import cycle
import matplotlib
from matplotlib import rc
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

#####################LINE STYLE CYCLER####################
lines = ["-","--","-.",":"]
linecycler = cycle(lines)
##########################################################

#######################SET LATEX OPTIONS###################
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern']})
#rc('text.latex',preamble=r'\usepackage{times}')
plt.rcParams.update({'font.size': 12})
plt.rcParams['xtick.major.pad'] = 8
plt.rcParams['ytick.major.pad'] = 8
##########################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="load energy intervals and compute Cv stdev", 
                                     epilog="must write file name followed by a label string, otherwise raise the flag --nolabels")
    parser.add_argument("fname", nargs="+", type=str, help="filenames with heat capacity followed by label")
    parser.add_argument("--legloc", type=int, help="define location of the legend (default upper right 1):"
                        " \"upper right=1\" \"upper left=2\" \"lower left=3\" \"lower right=4\" \"right=5\" "
                        "\"center left=6\" \"center right=7\" \"lower center=8\" \"upper center=9\" \"center=10\"",default=1)
    parser.add_argument("--xlabel", type=str, help="set x-label",default="T")
    parser.add_argument("--ylabel", type=str, help="set y-label",default="Cv")
    parser.add_argument("--xtop", type=float, help="set x-axis top",default=None)
    parser.add_argument("--xbot", type=float, help="set x-axis bottom",default=None)
    parser.add_argument("--ytop", type=float, help="set y-axis top",default=None)
    parser.add_argument("--ybot", type=float, help="set y-axis bottom",default=None)
    parser.add_argument("--title", type=str, help="set title",default=None)
    parser.add_argument("--linewidth", type=float, help="set line width (default 1.8)",default=1.8)
    parser.add_argument("--colormap", type=str, help="set colormap (default Dark2)",default='Dark2')
    parser.add_argument("--ecolor", type=str, help="set error bars color (default light gray)",default='g')
    parser.add_argument("--ecap", type=float, help="set error bars cap size (default is None)",default=None)
    parser.add_argument("--output", type=str, help="set output name (default cv_plot)",default='cv_plot')
    parser.add_argument("--filetype", type=str, help="set output file format (default eps)",default='eps')
    parser.add_argument("--nolabels", action="store_true", help="turn on if do not want to insert labels",default=False)
    parser.add_argument("--show", action="store_true", help="show explorable plot preview",default=False)
    parser.add_argument("--ebar", action="store_true", help="plot error bars if in data file",default=False)
    parser.add_argument("--grid", action="store_true", help="plot grid",default=False)
    args = parser.parse_args()
    
    fname = args.fname
    legloc=args.legloc
    nolabels = args.nolabels
    xlabel = args.xlabel
    ylabel = args.ylabel
    show = args.show
    title = args.title
    grid = args.grid
    linew = args.linewidth
    xtop = args.xtop
    xbot = args.xbot
    ytop = args.ytop
    ybot = args.ybot
    colormap = args.colormap
    errcolor = args.ecolor
    ebar = args.ebar
    filename = args.output
    filetype = args.filetype
    ecap = args.ecap
    
    ####################################################################################################
    #DEAL WITH INPUT
        
    if nolabels is False:
        all_data = [[] for i in xrange(len(fname)/2)]
        all_labels = np.array(["no_name" for i in xrange(len(fname)/2)])
        j = 0    
        for name,i in zip(fname,xrange(len(fname))):
            if ((i+1) % 2) is not 0:
                data = np.genfromtxt(name)
                dshape = np.shape(data)
                if dshape[1] > 3:
                    new_data = np.zeros([dshape[0],3])
                    for k in xrange(dshape[0]):
                        for l in xrange(3):
                            new_data[k][l] = data[k][l]
                    data = new_data  
                elif dshape[1] < 3:
                    data = np.hstack((data, np.zeros((data.shape[0], 1), dtype=data.dtype))) ##adds an error of 0 if there's no error associated
                all_data[j] = data
            else:
                all_labels[j] = name
                j += 1
    else:
        all_data = [[] for i in xrange(len(fname))]
        all_labels = np.array(["no_name" for i in xrange(len(fname))])
        for name,i in zip(fname,xrange(len(fname))):
            data = np.genfromtxt(name)
            dshape = np.shape(data)
            if dshape[1] > 3:
                new_data = np.zeros([dshape[0],3])
                for k in xrange(dshape[0]):
                    for l in xrange(3):
                       new_data[k][l] = data[k][l]
                data = new_data
            elif dshape[1] < 3:
                data = np.hstack((data, np.zeros((data.shape[0], 1), dtype=data.dtype)))
            data = np.array(data)
            all_data[i] = data
            
    all_data = np.array(all_data)
    ######################################################################################################
    
    #####################
    ####### PLOT ########
    #####################
    
    ####SET COLOUR MAP######
    cm = plt.get_cmap(colormap)
    ########################
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_color_cycle([cm(1.*i/len(all_data)) for i in xrange(np.shape(all_data)[0])])
    for data, label, i in zip(all_data, all_labels, xrange(np.shape(all_data)[0])): 
            ax.plot(data[:,0], data[:,1], next(linecycler), label = label, linewidth=linew)
            if np.shape(data)[1] is 3 and ebar is True:
                ax.errorbar(data[:,0], data[:,1], yerr=data[:,2], ecolor=errcolor, capsize=ecap )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim((xbot,xtop))
    ax.set_ylim((ybot,ytop))
    
    if grid is True:
        ax.grid(True,linestyle='--',color='0.75')
    if title is not None:
        ax.set_title(title)
    if nolabels is False:
        ax.legend(frameon=False,loc=legloc)
    if show is True:
         plt.show()
    
    fig.savefig('{F}.{T}'.format(F=filename,T=filetype))
    