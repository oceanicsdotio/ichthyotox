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
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.05, left=hpadding, right=1.0-hpadding)
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

#######################################
print "Reading State" # experiment A
data = np.loadtxt('../100/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol, :]
print "Calculating Mean Mass"
mass_avg2 =  0.001*np.mean(mass[:,:], axis=0)
ax[0].plot(time, mass_avg2, linestyle='-', linewidth=1, color='white', aa=True, label='Control (A)')
#######################################
print "Reading State" # experiment B
data = np.loadtxt('../101/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol, :]
print "Calculating Mean Mass"
mass_avg2 =  0.001*np.mean(mass[:,:], axis=0)
ax[0].plot(time, mass_avg2, linestyle='-', linewidth=1, color='red', aa=True, label='Formation (B)')
#######################################
print "Reading State" # experiment C
data = np.loadtxt('../102/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol,:]
print "Calculating Mean Mass"
mass_avg1 =  0.001*np.mean(mass[0:nfish/2-1,:], axis=0)
mass_avg2 =  0.001*np.mean(mass[nfish/2:nfish-1,:], axis=0)
ax[0].plot(time, mass_avg1, linestyle='-', linewidth=2, color='green', aa=True)
ax[0].plot(time, mass_avg2, linestyle='-', linewidth=1, color='green', aa=True, label='Intensification (C)')
#######################################
print "Reading State" # experiment D
data = np.loadtxt('../103/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol,:]
print "Calculating Mean Mass"
mass_avg1 = 0.001*np.mean(mass[0:nfish/2-1,:], axis=0)
mass_avg2 = 0.001*np.mean(mass[nfish/2:nfish-1,:], axis=0)
ax[0].plot(time, mass_avg1, linestyle='-', linewidth=2, color='blue', aa=True)
ax[0].plot(time, mass_avg2, linestyle='-', linewidth=1, color='blue', aa=True, label='Decline (D)')
#######################################

ax[0].set_xlabel(r'Days')
ax[0].xaxis.set_major_locator(MultipleLocator(5))
ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax[0].set_xlim(0,30)
ax[0].xaxis.grid(False); 

ax[0].set_ylabel('Mass, g', rotation=90)

ax[0].yaxis.set_major_locator(MultipleLocator(1.0))
ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax[0].yaxis.grid(False);
ax[0].set_frame_on(True)
ax[0].set_ylim(11.0,17.0)

legend = ax[0].legend(loc='upper left')
legend.get_frame().set_facecolor(bg_color)
text = legend.get_texts()
text[0].set_color('white')
text[1].set_color('red')
text[2].set_color('green')
text[3].set_color('blue')

#plt.show()
plt.savefig('./fish_growth_screen.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
