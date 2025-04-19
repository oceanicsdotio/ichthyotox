#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
from pylab import *
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi

nplots=1
fontsize  = 10
linewidth = 1
hpadding = 0.2
vpadding = 0.1
marginWidth = 3.5
fheight = 3.5
default_alpha = 1.0


rc('text', usetex=False)
rc('font', **{'family':'serif', 'serif':['Times New Roman']})
rc('mathtext', default='sf')
rc('lines', markeredgewidth=1)
rc('lines', linewidth=linewidth)
rc('axes', labelsize=fontsize)
rc('axes', linewidth=(linewidth+1)//2)
rc('xtick', labelsize=fontsize)
rc('ytick', labelsize=fontsize)
rc('legend', fontsize=fontsize)
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
fig.subplots_adjust(top=0.95, bottom=vpadding+0.05, left=hpadding, right=0.95, hspace=0.20)
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
plt.scatter(fpsx[ii,int(240*(event_start+1))], fpsy[ii,int(240*(event_start+1))], s=80, color='black', zorder=10, edgecolors='black', label='Control (A)') # end markers
for jj in range (int(240*event_start), int(240*(event_start+1))):
    if (sqrt( (fpsx[ii,jj+1] - fpsx[ii,jj])**2 + (fpsy[ii,jj+1] - fpsy[ii,jj])**2 ) < 200.0):
        ax[0].plot( (fpsx[ii,jj], fpsx[ii,jj+1]), (fpsy[ii,jj], fpsy[ii,jj+1]), linewidth=1, linestyle='-', color='black', alpha=step_alpha*(jj-240*event_start), aa=True)

# Experiment B, Particle 199
data = np.loadtxt('../101/fish_position.dat', unpack=True)
ii = 199
event_start=28
xcol = ii*dwidth + 2
ycol = ii*dwidth + 3
fpsx[ii,:] = data[xcol,:] 
fpsy[ii,:] = data[ycol,:] 
# plot path and end markers
plt.scatter(fpsx[ii,int(240*(event_start+1))], fpsy[ii,int(240*(event_start+1))], s=80, color='red', zorder=10, edgecolors='red', label='Formation (B)') # end markers
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
plt.scatter(fpsx[ii,int(240*(event_start+1))], fpsy[ii,int(240*(event_start+1))], s=80, color='green', zorder=10, edgecolors='green', label='Intensification (C)') # end markers
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
plt.scatter(fpsx[ii,240*(event_start+1)], fpsy[ii,240*(event_start+1)], s=80, color='blue', zorder=10, edgecolors='blue', label='Decline (D)') # end markers
for jj in range (int(240*event_start), int(240*(event_start+1))):
    if (sqrt( (fpsx[ii,jj+1] - fpsx[ii,jj])**2 + (fpsy[ii,jj+1] - fpsy[ii,jj])**2 ) < 200.0):
        ax[0].plot( (fpsx[ii,jj], fpsx[ii,jj+1]), (fpsy[ii,jj], fpsy[ii,jj+1]), linewidth=1, linestyle='-', color='blue', alpha=step_alpha*(jj-240*event_start), aa=True)


ax[0].set_xlabel(r'X (m)')
ax[0].set_xlim(0,500)
ax[0].xaxis.set_major_locator(MultipleLocator(125))
ax[0].xaxis.set_minor_locator(MultipleLocator(25))
ax[0].xaxis.grid(True)
ax[0].set_ylabel(r'Y (m)')
ax[0].set_ylim(0,500)
ax[0].yaxis.set_major_locator(MultipleLocator(125))
ax[0].yaxis.set_minor_locator(MultipleLocator(25))
ax[0].yaxis.grid(False)

# legend = ax[0].legend(loc='upper right')
# legend.get_frame().set_facecolor('none')
# legend.get_frame().set_edgecolor('none')
# text = legend.get_texts()
# text[0].set_color('black')
# text[1].set_color('red')
# text[2].set_color('green')
# text[3].set_color('blue')


plt.savefig('./fish_singlepaths_color.tiff', edgecolor='none', facecolor=bg_color, dpi=default_dpi)



