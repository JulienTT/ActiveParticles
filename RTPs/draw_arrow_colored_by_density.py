#! /usr/bin/python

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

font= FontProperties(size='10');

fig_width_pt = 251.0  # Get this from LaTeX using \showthe\columnwidth
                      # For PRL columnwidth is 246pt.

inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (sqrt(5)+1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width/1.2 #(fig_width/golden_mean*1.1) * 0.7    # height in inches
fig_size =  [fig_width,fig_height]

params = {'backend': 'Agg',
          'axes.labelsize': 10,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'figure.figsize': fig_size}

rcParams.update(params)

if len(sys.argv) > 1:
    datafile = sys.argv[1]
else:
    print "No data file name given. Please enter"
    datafile = raw_input("-> ")

output = datafile+".pdf"

print output

if len(sys.argv) > 2:
    rhomax=float(sys.argv[2])
else:
    rhomax=25
print rhomax

if len(sys.argv) > 3:
    Lx=float(sys.argv[3])
else:
    Lx=10

if len(sys.argv) > 4:
    Ly=float(sys.argv[4])
else:
    Ly=10

# I load the data into the array data
data=loadtxt(datafile)

#X is the first column, Y the second, etc.
X=data[:,0]
Y=data[:,1]
rho=data[:,2]
Jx=data[:,3]
Jy=data[:,4]

font = {'fontname'   : 'Times',
        'color'      : 'k',
        'fontsize'   : 10}

#I define a frame and give it the name ax [X,Y,DX,DY]
ax=axes([0.08,0.08,1.0-0.08,0.92-0.08])

ax.xaxis.set_major_locator(MultipleLocator(Lx/10))
ax.yaxis.set_major_locator(MultipleLocator(Ly/10))

#cmap = matplotlib.cm.get_cmap('viridis')
#normalize = matplotlib.colors.Normalize(vmin=0, vmax=rhomax)
#colors = [cmap(normalize(value)) for value in rho]

ax.set_xlim(0,Lx)
ax.set_ylim(0,Ly)

q1=ax.quiver(X,Y,Jx,Jy,rho,scale=30.0,clim=(0,rhomax))#,headlength=5.0,headaxislength=8.0,headwidth=7.0,clim=)
c1=plt.colorbar(q1,orientation='vertical')
c1.ax.set_xlabel(r'$\rho$')

savefig(output)
