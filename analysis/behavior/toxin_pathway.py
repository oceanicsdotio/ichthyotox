#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from pylab import *
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi

nplots=1
fontsize  = 10
linewidth = 1
hpadding = 0.1
vpadding = 0.05
marginWidth = 8.5
fheight = 6
default_alpha = 1.0


# load data
ini = open('../100/fish_ini.dat', 'r')
nfish = int(str.strip(ini.readline()))
time = np.loadtxt('../100/fish_state.dat', usecols=[0], unpack=True)
start = 0
end = len(time) - 1
dwidth = 4

data = zeros((nfish*dwidth+1,end-start+1))
path = zeros((nfish,end-start+1))
path_avg1 = zeros(end-start+1)
path_avg2 = zeros(end-start+1)
time = time / 24.

# figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.05, left=hpadding, right=1.0-(hpadding))
ax = []
ax.append(fig.add_subplot(1,1,1))

ax[0].patch.set_facecolor(bg_color)
ax[0].spines['top'].set_color(overlay_color)
ax[0].spines['bottom'].set_color(overlay_color)
ax[0].spines['left'].set_color(overlay_color)
ax[0].spines['right'].set_color(overlay_color)
ax[0].xaxis.label.set_color(label_color)
ax[0].yaxis.label.set_color(label_color)
ax[0].tick_params(axis='x', colors=label_color)
ax[0].tick_params(axis='y', colors=label_color)


data = np.loadtxt('../100/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    pathcol = ii*dwidth + 4
    path[ii,:] = data[pathcol,:]
path_avg1 = np.mean(path[:,:], axis=0)
ax[0].plot(time, log10(path_avg1), linestyle='-', linewidth=1, color='white', aa=True, label='Control (A)')

data = np.loadtxt('../101/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    pathcol = ii*dwidth + 4
    path[ii,:] = data[pathcol,:]
path_avg1 = np.mean(path[:,:], axis=0)
ax[0].plot(time, log10(path_avg1), linestyle='-', linewidth=1, color='red', aa=True, label='Formation (B)')

data = np.loadtxt('../102/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    pathcol = ii*dwidth + 4
    path[ii,:] = data[pathcol,:]
path_avg1 = np.mean(path[0:nfish/2-1,:], axis=0)
path_avg2 = np.mean(path[nfish/2:nfish-1,:], axis=0)
ax[0].plot(time, log10(path_avg1), linestyle='--', linewidth=1, color='green', aa=True)
ax[0].plot(time, log10(path_avg2), linestyle='-', linewidth=1, color='green', aa=True, label='Intensification (C)')

data = np.loadtxt('../103/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    pathcol = ii*dwidth + 4
    path[ii,:] = data[pathcol,:]
path_avg1 = np.mean(path[0:nfish/2-1,:], axis=0)
path_avg2 = np.mean(path[nfish/2:nfish-1,:], axis=0)
ax[0].plot(time, log10(path_avg1), linestyle='--', linewidth=1, color='blue', aa=True)
ax[0].plot(time, log10(path_avg2), linestyle='-', linewidth=1, color='blue', aa=True, label='Decline (D)')


ax[0].set_xlabel(r'Days')
ax[0].xaxis.set_major_locator(MultipleLocator(5))
ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax[0].set_xlim(0,30)
ax[0].xaxis.grid(False); 

ax[0].set_ylabel('Mean pathway partitioning ratio', rotation=90)

#ax[0].yaxis.set_major_locator(MultipleLocator(1))
#ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax[0].yaxis.grid(False)
ax[0].set_frame_on(True)
ax[0].set_ylim(0.0,10.0)

legend = ax[0].legend(loc='upper right')
legend.get_frame().set_facecolor('none')
text = legend.get_texts()
text[0].set_color('white')
text[1].set_color('red')
text[2].set_color('green')
text[3].set_color('blue')


plt.savefig('./fish_toxin_pathway_screen.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
