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

fig_height_pt = 120
inches_per_pt = 1.0/72.27               # Convert pt to inch
fig_height = fig_height_pt*inches_per_pt
fig_width = fig_height*1.2  # width in inches
fig_size =  [fig_width,fig_height]

params = {'backend': 'Agg',
          'axes.labelsize': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'text.usetex': True,
          'figure.figsize': fig_size}
rcParams.update(params)

output = "PhaseDiagram.png"
data1=loadtxt("DiagPhase-Dr-2-LG.txt")
data2=loadtxt("DiagPhase-Dr-2-MIPS.txt")

#X is the first column, Y the second, etc.
vLG=data1[:,0]
rhoLG_G=data1[:,1]
rhoLG_L=data1[:,3]

#for i in len(rhoLG_G):
#    rhoLG[i]=rhoLG_G[i]
#    vLGALL[i]=vLG[i]
#for i in len(rhoLG_L):
#    rhoLG[i+len(rhoLG_G)]=rhoLG_L[i]
#    vLGALL[i+len(rhoLG_G)]=vLG[i]


vMIPS=data2[:,0]
rhoMIPS_G=data2[:,1]
rhoMIPS_L=data2[:,3]

font = {'fontname'   : 'Times',
        'color'      : 'k',
        'fontsize'   : 8}

#I define a frame and give it the name ax [X,Y,DX,DY]
ax=axes([0.13,0.13,.99-0.13,0.99-0.13])

ax.xaxis.set_major_locator(MultipleLocator(0.3))
ax.yaxis.set_major_locator(MultipleLocator(10))

#plt.xticks([])
#plt.yticks([])

#ax.plot(rhoLG_G,vLG,c='b',linewidth=.5)

ax.fill_betweenx(vLG,rhoLG_G,rhoLG_L,color='cornflowerblue')
ax.scatter(rhoLG_G,vLG,marker='o',s=3,c='b')
ax.scatter(rhoLG_L,vLG,marker='o',s=3,c='b')

#ax.plot(rhoLG_L,vLG,c='b',linewidth=.5)

ax.fill_betweenx(vMIPS,rhoMIPS_G,rhoMIPS_L,color='peachpuff')
ax.scatter(rhoMIPS_G,vMIPS,marker='p',s=3,c='r')
#ax.plot(rhoMIPS_G,vMIPS,c='r',linewidth=.5)

ax.scatter(rhoMIPS_L,vMIPS,marker='p',s=3,c='r')
#ax.plot(rhoMIPS_L,vMIPS,c='r',linewidth=.5)

ax.set_xlim(0,1.4)
ax.set_ylim(0.05,25)

#ax.set_yscale('log')
#ax.set_title('time = {0}' .format(time))
ax.text(1.3, -3.5, r'$\rho_0$')

#ax.set_xlabel(r'$\rho_0$')
ax.text(-.2, 23, r'$v_0$')

savefig(output, dpi=300)
