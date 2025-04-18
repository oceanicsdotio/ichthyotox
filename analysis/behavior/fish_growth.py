#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
from pylab import *

nplots=4
fontsize  = 10
linewidth = 1
uniformPadding = 0.15
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = marginWidth + 0.5
default_alpha = 1.0

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
time = np.loadtxt('../100/fish_state.dat', usecols=[0], unpack=True)
start = 0
end = len(time) - 1
dwidth = 4

data = zeros((nfish*dwidth+1,end-start+1))
mass = zeros((nfish,end-start+1))
suitability = zeros((nfish,end-start+1))
xpos = zeros((nfish,end-start+1))
mass_avg1 = zeros(end-start+1)
mass_avg2 = zeros(end-start+1)
suit_avg1 = zeros(end-start+1)
suit_avg2 = zeros(end-start+1)
mass_std = zeros(end-start+1)
suit_std = zeros(end-start+1)
time = time / 24.


# figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.05, left=hpadding, right=1.0-(hpadding/2.), hspace=0.20)
ax = []
ax.append(fig.add_subplot(4,1,1))
ax.append(fig.add_subplot(4,1,2))
ax.append(fig.add_subplot(4,1,3))
ax.append(fig.add_subplot(4,1,4))

ax[0].patch.set_facecolor(bg_color)
ax[0].spines['top'].set_color(overlay_color)
ax[0].spines['bottom'].set_color(overlay_color)
ax[0].spines['left'].set_color(overlay_color)
ax[0].spines['right'].set_color(overlay_color)
ax[0].xaxis.label.set_color(label_color)
ax[0].yaxis.label.set_color(label_color)
ax[0].tick_params(axis='x', colors=label_color)
ax[0].tick_params(axis='y', colors=label_color)

#######################################
print "Reading State" # experiment A
data = np.loadtxt('../100/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol, :]
print "Calculating Mean Mass"
mass_avg1 = np.mean(mass[0:99,:], axis=0)
mass_avg2 = np.mean(mass[100:199,:], axis=0)
ax[0].plot(time, mass_avg1, linestyle='-', linewidth=1, color='red', aa=True)
ax[0].plot(time, mass_avg2, linestyle='-', linewidth=1, color='blue', aa=True)
#######################################
ax.append(ax[0].twinx())
print "Reading Position" # experiment A
data = np.loadtxt('../100/fish_position.dat', unpack=True)
for ii in range(0, nfish):
    xcol = ii*dwidth + 2
    xpos[ii,:] = data[xcol,:]
print "Calculating Mean Suitability"
suitability = 0.5*(1.0 + sin( 2.0 *pi * (xpos[:,:]+125.0) / 500.0 ))
suit_avg1 = np.mean(suitability[0:99,:], axis=0)
suit_avg2 = np.mean(suitability[100:199,:], axis=0)
ax[0+nplots].plot(time, suit_avg1, linestyle='--', linewidth=1, color='red', aa=True)
ax[0+nplots].plot(time, suit_avg2, linestyle='--', linewidth=1, color='blue', aa=True)
#######################################

#######################################
ax.append(ax[1].twinx())
#######################################
print "Reading State" # experiment C
data = np.loadtxt('../102/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol,:]
print "Calculating Mean Mass"
mass_avg1 = np.mean(mass[0:99,:], axis=0)
mass_avg2 = np.mean(mass[100:199,:], axis=0)
ax[2].plot(time, mass_avg1, linestyle='-', linewidth=1, color='red', aa=True)
ax[2].plot(time, mass_avg2, linestyle='-', linewidth=1, color='blue', aa=True)
#######################################
ax.append(ax[2].twinx())
print "Reading Position" # experiment C
data = np.loadtxt('../102/fish_position.dat', unpack=True)
for ii in range(0, nfish):
    xcol = ii*dwidth + 2
    xpos[ii,:] = data[xcol,:]
print "Calculating Mean Suitability"
suitability = 0.5*(1.0 + sin( 2.0 *pi * (xpos[:,:]+125.0) / 500.0 ))
suit_avg1 = np.mean(suitability[0:99,:], axis=0)
suit_avg2 = np.mean(suitability[100:199,:], axis=0)
ax[2+nplots].plot(time, suit_avg1, linestyle='--', linewidth=1, color='red', aa=True)
ax[2+nplots].plot(time, suit_avg2, linestyle='--', linewidth=1, color='blue', aa=True)
#######################################

sublabels = [[r'Days', r'Mass (g)'],[r'Days', r'Mass (g)'],[r'Days', r'Mass (g)'],[r'Days', r'Mass (g)']]
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

        
    ax[ii].set_xlim(time[start],time[end])
    ax[ii].xaxis.grid(False); 
    
    ax[ii].set_ylabel(sublabels[ii][1], rotation=90)
    
    #ax[ii].yaxis.set_major_locator(mticker.MultipleLocator(0.05))
    #ax[ii].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[ii].yaxis.grid(False);
    
    ax[ii].set_title(subtitles[ii], fontsize=fontsize, color=label_color)
    ax[ii].title.set_position([1.0, 1.0])
    ax[ii].set_frame_on(True)

# ax[0].set_ylim(0.0,0.1)
# ax[1].set_ylim(0.0,0.35)
# ax[2].set_ylim(0.0,0.45)
# ax[3].set_ylim(0.0,0.4)

#plt.show()
plt.savefig('./fish_growth.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
