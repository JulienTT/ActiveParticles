#!/usr/bin/python3
import subprocess
import matplotlib

matplotlib.use('Agg')

from pylab import *
from matplotlib.axes import Axes
from matplotlib.patches import Polygon
from matplotlib.ticker import *
from matplotlib.font_manager import fontManager, FontProperties 
from numpy import loadtxt
import subprocess
import math
from matplotlib import pyplot as plt
#from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = 1.5*fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean     # height in inches
fig_size =  [fig_width,fig_height]


params = {'text.usetex': True,
          'legend.fontsize': 15,
          'legend.borderpad' : 0.,
          'legend.labelspacing' : 0,
          'legend.frameon' : 0,
          'legend.handlelength': 0.1,
          'figure.figsize': fig_size}
matplotlib.rcParams.update(params)

if len(sys.argv) > 1:
    datafile = sys.argv[1]
else:
    print("No data file name given. Please enter")
    datafile = raw_input("-> ")

output = datafile+".pdf"

if len(sys.argv) > 2:
    xmax = float(sys.argv[2])
    if (xmax<1.5):
        xmax=1.5
else:
    xmax=float(2.5)


data=loadtxt(datafile)

X=data[:,0]
Y=data[:,1]

font = {'fontname'   : 'Times',
        'color'      : 'k',
        'fontsize'   : 10}

#I define a frame and give it the name ax [X,Y,DX,DY]
ax=axes([0.1,0.15,.98-0.15,0.93-0.15])

#ax.xaxis.set_major_locator(MultipleLocator(.1))
#ax.yaxis.set_major_locator(MultipleLocator(Ly/10))
ax.set_title('Histogram of local densities')
ax.set_xlabel(r'$P(\rho)$')
ax.set_ylabel(r'$\rho$')
ax.set_xlim([0,xmax])
ax.plot(X,Y,'bo-',linewidth=2,markersize=5)

savefig(output, dpi=300)
plt.close('all')
