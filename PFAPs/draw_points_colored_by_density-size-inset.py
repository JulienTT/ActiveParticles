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
from matplotlib import pyplot as plt

font= FontProperties(size='10');

#fig_width_pt = 251.0  # Get this from LaTeX using \showthe\columnwidth
                      # For PRL columnwidth is 246pt.
fig_width_pt = 351.0

                      
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

output = datafile+".png"

print output

if len(sys.argv) > 2:
    rhomax=float(sys.argv[2])
else:
    rhomax=25
print rhomax

if len(sys.argv) > 3:
    time=float(sys.argv[3])
else:
    time='time'

if len(sys.argv) > 4:
    Lx=float(sys.argv[4])
else:
    Lx=10

if len(sys.argv) > 5: 
   Ly=float(sys.argv[5])
else:
    Ly=10

if len(sys.argv) > 6:
    _radius=float(sys.argv[6])
else:
    _radius=1

    
#print 'time: {0}' .format(time)

# I load the data into the array data
data=loadtxt(datafile)

#X is the first column, Y the second, etc.
X=data[:,2]
Y=data[:,3]
rho=data[:,5]+1/3.141596

font = {'fontname'   : 'Times',
        'color'      : 'k',
        'fontsize'   : 10}

#I define a frame and give it the name ax [X,Y,DX,DY]
ax=axes([0.08,0.08,1.0-0.04,0.92-0.08])

#ax.xaxis.set_major_locator(MultipleLocator(Lx/10))
#ax.yaxis.set_major_locator(MultipleLocator(Ly/10))

#Erasing the x and y coordinates
plt.xticks([])
plt.yticks([])

cmap = matplotlib.cm.get_cmap('plasma')
#normalize = matplotlib.colors.Normalize(vmin=0, vmax=rhomax)
normalize = matplotlib.colors.Normalize(vmin=0, vmax=12)#temporaire pour le PRL TAPS
colors = [cmap(normalize(value)) for value in rho]

for a, b, color in zip (X, Y, colors):
    circle = plt.Circle((a, b), _radius, color=color, fill=True, linewidth=.5,zorder=0)
    ax.add_artist(circle)

#q1=ax.scatter(X,Y,s=0,c=colors,vmin=0, vmax=rhomax,edgecolors='none')
q1=ax.scatter(X,Y,s=0,c=colors,vmin=0, vmax=12,edgecolors='none')#temporaire pour TAPS PRL
#ax.scatter(X,Y,s=1,c='k',edgecolors='none',zorder=10)#sert a mettre les points noirs au centre des particules
ax.set_xlim(0,Lx)
ax.set_ylim(0,Ly)
#ax.set_title('time = {0}' .format(time))

cax, _ = matplotlib.colorbar.make_axes(ax)
cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize)
# Modify the fontsize of the colorbar labels 
ticklabs = cbar.ax.get_yticklabels()
cbar.ax.set_yticklabels(ticklabs,fontsize=20)




#axes([debutagauche,debutenbas,largeur,hauteur)
bx=axes([0.53,0.6,.3,.3])

for a, b, color in zip (X, Y, colors):
    circle = plt.Circle((a, b), _radius, color=color, fill=True, linewidth=.5,zorder=0)
    bx.add_artist(circle)

q2=bx.scatter(X,Y,s=0,c=colors,vmin=0, vmax=rhomax,edgecolors='none')
#bx.scatter(X,Y,s=1,c='k',edgecolors='none',zorder=10)
bx.set_xlim(Lx/7,2*Lx/7)
bx.set_ylim(Ly/7,2*Ly/7)
plt.xticks([])
plt.yticks([])

#c1=plt.colorbar(q1,orientation='vertical')
#c1.ax.set_xlabel(r'$\rho$')

savefig(output, dpi=300)
