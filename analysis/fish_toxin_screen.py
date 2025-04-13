#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
from pylab import *

nplots=1
fontsize  = 10
linewidth = 1
hpadding = 0.1
vpadding = 0.05
marginWidth = 8.5
fheight = 6
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
rc('legend', fontsize=fontsize)
rc('xtick.major', pad=5)
rc('ytick.major', pad=5)


# load data
ini = open('../100/fish_ini.dat', 'r')
nfish = int(str.strip(ini.readline()))
time = np.loadtxt('../100/fish_state.dat', usecols=[0], unpack=True)
start = 0
end = len(time) - 1
dwidth = 4

data = zeros((nfish*dwidth+1,end-start+1))
mass = zeros((nfish,end-start+1))
tox = zeros((nfish,end-start+1))
path = zeros((nfish,end-start+1))
ratio = zeros((nfish,end-start+1))
ratio_avg = zeros(end-start+1)
ratio_avg1 = zeros(end-start+1)
ratio_avg2 = zeros(end-start+1)
ratio_std = zeros(end-start+1)
time = time / 24.

# figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.05, left=hpadding, right=1.0-(hpadding))
ax = []
ax.append(fig.add_subplot(1,1,1))
#ax.append(ax[0].twinx())

ax[0].patch.set_facecolor(bg_color)
ax[0].spines['top'].set_color(overlay_color)
ax[0].spines['bottom'].set_color(overlay_color)
ax[0].spines['left'].set_color(overlay_color)
ax[0].spines['right'].set_color(overlay_color)
ax[0].xaxis.label.set_color(label_color)
ax[0].yaxis.label.set_color(label_color)
ax[0].tick_params(axis='x', colors=label_color)
ax[0].tick_params(axis='y', colors=label_color)

# ax[1].xaxis.label.set_color(label_color)
# ax[1].yaxis.label.set_color(label_color)
# ax[1].tick_params(axis='y', colors=label_color)

print "Reading State" # experiment A
data = np.loadtxt('../100/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    toxcol = ii*dwidth + 3
    pathcol = ii*dwidth + 4
    mass[ii,:] = data[masscol,:]
    tox[ii,:] = data[toxcol,:]
    path[ii,:] = data[pathcol,:]
print "Calculating Mean Tox Load and Pathway Partitioning"
ratio = tox/mass
ratio_avg = np.mean(ratio[0:nfish/2-1,:], axis=0)
ax[0].plot(time, 1000*ratio_avg, linestyle='-', linewidth=1, color='white', aa=True, label='Control (A)')
#######################################
print "Reading State" # experiment B
data = np.loadtxt('../101/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    toxcol = ii*dwidth + 3
    pathcol = ii*dwidth + 4
    mass[ii,:] = data[masscol,:]
    tox[ii,:] = data[toxcol,:]
    path[ii,:] = data[pathcol,:]
print "Calculating Mean Tox Load and Pathway Partitioning"
ratio = tox/mass
ratio_avg = np.mean(ratio[:,:], axis=0)
ratio_std = np.std(ratio[:,:], axis=0)
ax[0].plot(time, 1000*ratio_avg, linestyle='-', linewidth=1, color='red', aa=True, label='Formation (B)', zorder = 2)
ax[0].fill_between(time, 1000*(ratio_avg-ratio_std), 1000*(ratio_avg+ratio_std), facecolor=[1.0,0.0,0.0,0.25], edgecolor='red', zorder=2)
#######################################
print "Reading State" # experiment C
data = np.loadtxt('../102/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    toxcol = ii*dwidth + 3
    pathcol = ii*dwidth + 4
    mass[ii,:] = data[masscol,:]
    tox[ii,:] = data[toxcol,:]
    path[ii,:] = data[pathcol,:]
print "Calculating Mean Tox Load and Pathway Partitioning"
ratio = tox/mass
ratio_avg = np.mean(ratio[:,:], axis=0)
ratio_std = np.std(ratio[:,:], axis=0)
ratio_avg1 = np.mean(ratio[0:nfish/2-1,:], axis=0)
ratio_avg2 = np.mean(ratio[nfish/2:nfish-1,:], axis=0)
ax[0].fill_between(time, 1000*(ratio_avg-ratio_std), 1000*(ratio_avg+ratio_std), facecolor=[0.0,1.0,0.0,0.25], edgecolor='green', zorder=1)
ax[0].plot(time, 1000*ratio_avg1, linestyle=':', linewidth=1, color='green', aa=True, zorder = 1)
ax[0].plot(time, 1000*ratio_avg2, linestyle='-', linewidth=1, color='green', aa=True, zorder = 1, label='Intensification (C)')
#######################################
print "Reading State" # experiment D
data = np.loadtxt('../103/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    toxcol = ii*dwidth + 3
    pathcol = ii*dwidth + 4
    mass[ii,:] = data[masscol,:]
    tox[ii,:] = data[toxcol,:]
    path[ii,:] = data[pathcol,:]
print "Calculating Mean Tox Load and Pathway Partitioning"
ratio = tox/mass
ratio_avg = np.mean(ratio[:,:], axis=0)
ratio_std = np.std(ratio[:,:], axis=0)
ratio_avg1 = np.mean(ratio[0:nfish/2-1,:], axis=0)
ratio_avg2 = np.mean(ratio[nfish/2:nfish-1,:], axis=0)
ax[0].fill_between(time, 1000*(ratio_avg-ratio_std), 1000*(ratio_avg+ratio_std), facecolor=[0.0,0.0,1.0,0.25], edgecolor='blue', zorder=3)
ax[0].plot(time, 1000*ratio_avg1, linestyle=':', linewidth=1, color='blue', aa=True, zorder = 3)
ax[0].plot(time, 1000*ratio_avg2, linestyle='-', linewidth=1, color='blue', aa=True, zorder = 3, label='Decline (D)')
#######################################
# intoxication threshold
ax[0].plot((time[start], time[end]),(0.005, 0.005), linestyle='--', linewidth=1, color='white', aa=True, zorder = 5)
ax[0].plot((time[start], time[end]),(0.015, 0.015), linestyle='--', linewidth=1, color='white', aa=True, zorder = 5)


ax[0].set_xlabel(r'Days')
ax[0].xaxis.set_major_locator(MultipleLocator(5))
ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax[0].set_xlim(0,30)
ax[0].xaxis.grid(False); 

ax[0].set_ylabel('Mean toxin ratio, ug/g', rotation=90)

#ax[0].yaxis.set_major_locator(MultipleLocator(1))
#ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax[0].yaxis.grid(False);
ax[0].set_frame_on(True)
ax[0].set_ylim(0,)

legend = ax[0].legend(loc='upper left')
legend.get_frame().set_facecolor('none')
text = legend.get_texts()
text[0].set_color('white')
text[1].set_color('red')
text[2].set_color('green')
text[3].set_color('blue')


plt.savefig('./fish_toxin_screen.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
