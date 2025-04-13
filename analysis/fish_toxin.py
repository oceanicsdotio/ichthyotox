#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
from pylab import *

nplots=3
fontsize  = 10
linewidth = 1
uniformPadding = 0.15
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = marginWidth + 0.5
default_alpha = 1.0
lineRGBA = [1.0, 0.0, 0.0, default_alpha]
style = ['-',':']
div=150
show_grid=True

for_screen=False
if for_screen:
    bg_color = [0.0,0.0,0.0,default_alpha]
    overlay_color = [1.0,1.0,1.0,default_alpha]
    label_color = [0.5,0.5,0.5,default_alpha]
    default_dpi = 150
else:
    bg_color = [1.0,1.0,1.0,default_alpha]
    overlay_color = [0.0,0.0,0.0,default_alpha]
    label_color = [0.0,0.0,0.0,default_alpha]
    default_dpi = 300


rc('text', usetex=False)
rc('font', **{'family':'serif', 'serif':['Times New Roman']})
#rc('font', weight='bold')
rc('mathtext', default='sf')
rc('lines', markeredgewidth=1)
rc('lines', linewidth=linewidth)
#rc('axes', labelsize=fontsize)
rc('axes', labelsize=fontsize)
rc('axes', linewidth=(linewidth+1)//2)
rc('xtick', labelsize=fontsize)
rc('ytick', labelsize=fontsize)
#rc('legend', fontsize=2*fontsize/3)
rc('xtick.major', pad=5)
rc('ytick.major', pad=5)


# load data
ini = open('../100/fish_ini.dat', 'r')
#nfish = int(str.strip(ini.readline()))
nfish=4
time = np.loadtxt('../100/fish_state.dat', usecols=[0], unpack=True)
start = 0
end = len(time) - 1
dwidth = 4

mass = zeros((nfish,end-start+1))
tox = zeros((nfish,end-start+1))
path = zeros((nfish,end-start+1))
ratio = zeros((nfish,end-start+1))

ratio_avg1 = zeros(end-start+1)
ratio_avg2 = zeros(end-start+1)
path_avg1 = zeros(end-start+1)
path_avg2 = zeros(end-start+1)

time = time / 24.


# figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.05, left=hpadding, right=1.0-(hpadding), hspace=0.20)
ax = []
ax.append(fig.add_subplot(3,1,1))
ax.append(fig.add_subplot(3,1,2))
ax.append(fig.add_subplot(3,1,3))

ax[0].patch.set_facecolor(bg_color)
ax[0].spines['top'].set_color(overlay_color)
ax[0].spines['bottom'].set_color(overlay_color)
ax[0].spines['left'].set_color(overlay_color)
ax[0].spines['right'].set_color(overlay_color)
ax[0].xaxis.label.set_color(label_color)
ax[0].yaxis.label.set_color(label_color)
ax[0].tick_params(axis='x', colors=label_color)
ax[0].tick_params(axis='y', colors=label_color)


print "Reading State" # experiment B
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    toxcol = ii*dwidth + 3
    pathcol = ii*dwidth + 4
    #mass[ii,:], tox[ii,:], path[ii,:] = np.loadtxt('../100/fish_state.dat', usecols=(masscol, toxcol, pathcol), unpack=True)
print "Calculating Mean Tox Load and Pathway Partitioning"
# ratio = tox/mass
# ratio_avg1 = 1000*np.mean(ratio[0:nfish/2-1,:], axis=0)
# ratio_avg2 = 1000*np.mean(ratio[nfish/2:nfish-1,:], axis=0)
# path_avg1 = np.mean(path[0:nfish/2-1,:], axis=0)
# path_avg2 = np.mean(path[nfish/2:nfish-1,:], axis=0)
#ax[0].plot(time, ratio_avg1, linestyle='-', linewidth=1, color='red', aa=True)
#ax[0].plot(time, ratio_avg2, linestyle='-', linewidth=1, color='blue', aa=True)

ax.append(ax[0].twinx())
# ax[0+nplots].plot(time, path_avg1, linestyle='--', linewidth=1, color='red', aa=True)
# ax[0+nplots].plot(time, path_avg2, linestyle='--', linewidth=1, color='blue', aa=True)



print "Reading State" # experiment C
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    toxcol = ii*dwidth + 3
    pathcol = ii*dwidth + 4
    mass[ii,:], tox[ii,:], path[ii,:] = np.loadtxt('../102/fish_state.dat', usecols=(masscol, toxcol, pathcol), unpack=True)
print "Calculating Mean Tox Load and Pathway Partitioning"
ratio = tox/mass
ratio_avg1 = 1000*np.mean(ratio[0:nfish/2-1,:], axis=0)
ratio_avg2 = 1000*np.mean(ratio[nfish/2:nfish-1,:], axis=0)
path_avg1 = np.mean(path[0:nfish/2-1,:], axis=0)
path_avg2 = np.mean(path[nfish/2:nfish-1,:], axis=0)
ax[1].plot(time, ratio_avg1, linestyle='-', linewidth=1, color='red', aa=True)
ax[1].plot(time, ratio_avg2, linestyle='-', linewidth=1, color='blue', aa=True)

ax.append(ax[1].twinx())
ax[1+nplots].plot(time, path_avg1, linestyle='--', linewidth=1, color='red', aa=True)
ax[1+nplots].plot(time, path_avg2, linestyle='--', linewidth=1, color='blue', aa=True)

print "Reading State" # experiment D
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    toxcol = ii*dwidth + 3
    pathcol = ii*dwidth + 4
    mass[ii,:], tox[ii,:], path[ii,:] = np.loadtxt('../103/fish_state.dat', usecols=(masscol, toxcol, pathcol), unpack=True)
print "Calculating Mean Tox Load and Pathway Partitioning"
ratio = tox/mass
ratio_avg1 = 1000*np.mean(ratio[0:nfish/2-1,:], axis=0)
ratio_avg2 = 1000*np.mean(ratio[nfish/2:nfish-1,:], axis=0)
path_avg1 = np.mean(path[0:nfish/2-1,:], axis=0)
path_avg2 = np.mean(path[nfish/2:nfish-1,:], axis=0)
ax[2].plot(time, ratio_avg1, linestyle='-', linewidth=1, color='red', aa=True)
ax[2].plot(time, ratio_avg2, linestyle='-', linewidth=1, color='blue', aa=True)

ax.append(ax[2].twinx())
ax[2+nplots].plot(time, path_avg1, linestyle='--', linewidth=1, color='red', aa=True)
ax[2+nplots].plot(time, path_avg2, linestyle='--', linewidth=1, color='blue', aa=True)


sublabels = [[r'Days', r'Toxin (g/g)'],[r'Days', r'Toxin (g/g)'],[r'Days', r'Toxin (g/g)'],[r'Days', r'Toxin (g/g)']]
subtitles = [r'(A)', r'(B)',r'(C)', r'(D)']
for ii in range(0,nplots):
    # axis adjustment
    if (ii == nplots-1):
        ax[ii].set_xlabel(sublabels[ii][0])
        ax[ii].xaxis.set_major_locator(MultipleLocator(5))
        ax[ii].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    else:
        ax[ii].set_xlabel(' ')
        ax[ii].xaxis.set_major_locator(LinearLocator(0))

        
    ax[ii].set_xlim(0,30)
    ax[ii].xaxis.grid(False); 
    ax[ii].set_ylabel(sublabels[ii][1], rotation=90)
    
    #ax[ii].yaxis.set_major_locator(mticker.MultipleLocator(0.05))
    #ax[ii].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[ii].yaxis.grid(False);
    
    ax[ii].set_title(subtitles[ii], fontsize=fontsize, color=label_color)
    ax[ii].title.set_position([1.0, 1.0])
    ax[ii].set_frame_on(True)

ax[0+nplots].set_ylim(0.0,0.1)
ax[1+nplots].set_ylim(0.0,0.1)
ax[2+nplots].set_ylim(0.0,0.1)

#plt.show()
plt.savefig('./fish_toxin.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
