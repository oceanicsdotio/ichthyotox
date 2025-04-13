#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
from pylab import *
import matplotlib.tri as tri
import math

nplots=1
fontsize  = 10
linewidth = 1
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.0
fheight = 6.0
default_alpha = 1.0

for_screen=True
if for_screen:
    bg_color = [0.0,0.0,0.0,default_alpha]
    overlay_color = [1.0,1.0,1.0,default_alpha]
    label_color = [1.0,1.0,1.0,default_alpha]
    default_dpi = 150
else:
    bg_color = [1.0,1.0,1.0,default_alpha]
    overlay_color = [0.0,0.0,0.0,default_alpha]
    label_color = [0.0,0.0,0.0,default_alpha]
    default_dpi = 300

rc('text', usetex=False)
rc('font', **{'family':'sans-serif', 'sans-serif':['Arial']})
rc('mathtext', default='sf')
rc('lines', markeredgewidth=1)
rc('lines', linewidth=linewidth)
rc('axes', labelsize=fontsize)
rc('axes', linewidth=(linewidth+1)//2)
rc('xtick', labelsize=fontsize)
rc('ytick', labelsize=fontsize)
rc('xtick.major', pad=5)
rc('ytick.major', pad=5)


# load data
ini = open('../100/fish_ini.dat', 'r')
nfish = int(str.strip(ini.readline()))
time = np.loadtxt('../100/fish_position.dat', usecols=[0], unpack=True)
time = time/24.
start = 0
end = len(time) - 1
dwidth = 4
fpsx = zeros((nfish,end-start+1))
fpsy = zeros((nfish,end-start+1))
data = zeros((nfish*dwidth+1,end-start+1))
step_alpha = 1./240.

# Figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.05, left=hpadding, right=1.0-hpadding/2., hspace=0.20)
ax = []
ax.append(fig.add_subplot(1,1,1))

for ii in range (0,nplots):
    ax[ii].patch.set_facecolor(bg_color)
    ax[ii].spines['top'].set_color(overlay_color)
    ax[ii].spines['bottom'].set_color(overlay_color)
    ax[ii].spines['left'].set_color(overlay_color)
    ax[ii].spines['right'].set_color(overlay_color)
    ax[ii].xaxis.label.set_color(label_color)
    ax[ii].yaxis.label.set_color(label_color)
    ax[ii].tick_params(axis='x', colors=label_color)
    ax[ii].tick_params(axis='y', colors=label_color)

# Experiment A, Particle 199
data = np.loadtxt('../100/fish_position.dat', unpack=True)
ii = 199
event_start=19.5
xcol = ii*dwidth + 2
ycol = ii*dwidth + 3
fpsx[ii,:] = data[xcol,:] 
fpsy[ii,:] = data[ycol,:] 
# plot path and end markers
print "Drawing paths..."
plt.scatter(fpsx[ii,int(240*(event_start+1))], fpsy[ii,int(240*(event_start+1))], s=80, color='white', zorder=10, edgecolors='black') # end markers
for jj in range (int(240*event_start), int(240*(event_start+1))):
    if (sqrt( (fpsx[ii,jj+1] - fpsx[ii,jj])**2 + (fpsy[ii,jj+1] - fpsy[ii,jj])**2 ) < 200.0):
        ax[0].plot( (fpsx[ii,jj], fpsx[ii,jj+1]), (fpsy[ii,jj], fpsy[ii,jj+1]), linewidth=1, linestyle='-', color='white', alpha=step_alpha*(jj-240*event_start), aa=True)

# Experiment B, Particle 199
data = np.loadtxt('../101/fish_position.dat', unpack=True)
ii = 199
event_start=28
xcol = ii*dwidth + 2
ycol = ii*dwidth + 3
fpsx[ii,:] = data[xcol,:] 
fpsy[ii,:] = data[ycol,:] 
# plot path and end markers
print "Drawing paths..."
plt.scatter(fpsx[ii,int(240*(event_start+1))], fpsy[ii,int(240*(event_start+1))], s=80, color='red', zorder=10, edgecolors='black') # end markers
for jj in range (int(240*event_start), int(240*(event_start+1))):
    if (sqrt( (fpsx[ii,jj+1] - fpsx[ii,jj])**2 + (fpsy[ii,jj+1] - fpsy[ii,jj])**2 ) < 200.0):
        ax[0].plot( (fpsx[ii,jj], fpsx[ii,jj+1]), (fpsy[ii,jj], fpsy[ii,jj+1]), linewidth=1, linestyle='-', color='red', alpha=step_alpha*(jj-240*event_start), aa=True)

# Experiment C, Particle 199
data = np.loadtxt('../102/fish_position.dat', unpack=True)
ii = 199
event_start=3
xcol = ii*dwidth + 2
ycol = ii*dwidth + 3
fpsx[ii,:] = data[xcol,:] 
fpsy[ii,:] = data[ycol,:] 
# plot path and end markers
print "Drawing paths..."
plt.scatter(fpsx[ii,int(240*(event_start+1))], fpsy[ii,int(240*(event_start+1))], s=80, color='green', zorder=10, edgecolors='black') # end markers
for jj in range (int(240*event_start), int(240*(event_start+1))):
    if (sqrt( (fpsx[ii,jj+1] - fpsx[ii,jj])**2 + (fpsy[ii,jj+1] - fpsy[ii,jj])**2 ) < 200.0):
        ax[0].plot( (fpsx[ii,jj], fpsx[ii,jj+1]), (fpsy[ii,jj], fpsy[ii,jj+1]), linewidth=1, linestyle='-', color='green', alpha=step_alpha*(jj-240*event_start), aa=True)


# Experiment D, Particle 199
data = np.loadtxt('../103/fish_position.dat', unpack=True)
ii = 199
event_start=10
xcol = ii*dwidth + 2
ycol = ii*dwidth + 3
fpsx[ii,:] = data[xcol,:] 
fpsy[ii,:] = data[ycol,:] 
# plot path and end markers
print "Drawing paths..."
plt.scatter(fpsx[ii,240*(event_start+1)], fpsy[ii,240*(event_start+1)], s=80, color='blue', zorder=10, edgecolors='black') # end markers
for jj in range (int(240*event_start), int(240*(event_start+1))):
    if (sqrt( (fpsx[ii,jj+1] - fpsx[ii,jj])**2 + (fpsy[ii,jj+1] - fpsy[ii,jj])**2 ) < 200.0):
        ax[0].plot( (fpsx[ii,jj], fpsx[ii,jj+1]), (fpsy[ii,jj], fpsy[ii,jj+1]), linewidth=1, linestyle='-', color='blue', alpha=step_alpha*(jj-240*event_start), aa=True)


ax[0].set_xlabel(r'X')
ax[0].set_xlim(0,500)
ax[0].xaxis.set_major_locator(MultipleLocator(125))
ax[0].xaxis.grid(True)
ax[0].set_ylabel(r'Y')
ax[0].set_ylim(0,500)
ax[0].yaxis.set_major_locator(MultipleLocator(125))
ax[0].yaxis.grid(False)


plt.savefig('./fish_singlepaths_screen.tiff', edgecolor='none', facecolor=bg_color, dpi=default_dpi)



