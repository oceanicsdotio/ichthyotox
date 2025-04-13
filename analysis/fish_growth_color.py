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
marginWidth = 6.5
fheight = 3.0
default_alpha = 1.0

for_screen=False
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
time = np.loadtxt('../100/fish_state.dat', usecols=[0], unpack=True)
start = 0
end = len(time) - 1
dwidth = 4

data = zeros((nfish*dwidth+1,end-start+1))
mass = zeros((nfish,end-start+1))
suitability = zeros((nfish,end-start+1))
xpos = zeros((nfish,end-start+1))
mass_control = zeros(end-start+1)
mass_avg_all = zeros(end-start+1)
mass_avg_surface = zeros(end-start+1)
mass_avg_bottom = zeros(end-start+1)
time = time / 24.


# figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.1, left=hpadding, right=1.0-hpadding)
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
mass_control =  np.mean(mass[:,:], axis=0)
ax[0].plot(time, mass_control, linestyle='-', linewidth=1, color='black', aa=True, label='Control (A)')
#######################################
print "Reading State" # experiment B
data = np.loadtxt('../101/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol, :]
print "Calculating Mean Mass"
mass_avg =  np.mean(mass[:,:], axis=0)
ax[0].plot(time, mass_avg, linestyle='-', linewidth=1, color='red', aa=True, label='Formation (B)')
#######################################
print "Reading State" # experiment C
data = np.loadtxt('../102/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol,:]
print "Calculating Mean Mass"
mass_avg_surface = np.mean(mass[0:nfish/2-1,:], axis=0)
mass_avg_bottom = np.mean(mass[nfish/2:nfish-1,:], axis=0)
ax[0].plot(time, mass_avg_bottom, linestyle='-', linewidth=1, color='green', aa=True, label='Intensification (C)')
ax[0].plot(time, mass_avg_surface, linestyle=':', linewidth=1, color='green', aa=True)
#######################################
print "Reading State" # experiment D
data = np.loadtxt('../103/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol,:]
print "Calculating Mean Mass"
mass_avg_surface = np.mean(mass[0:nfish/2-1,:], axis=0)
mass_avg_bottom = np.mean(mass[nfish/2:nfish-1,:], axis=0)
ax[0].plot(time, mass_avg_bottom, linestyle='-', linewidth=1, color='blue', aa=True, label='Decline (D)')
ax[0].plot(time, mass_avg_surface, linestyle=':', linewidth=1, color='blue', aa=True)
#######################################
print "Reading State" # experiment D
data = np.loadtxt('../302/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol,:]
print "Calculating Mean Mass"
mass_avg_bottom = np.mean(mass[:,:], axis=0)
ax[0].plot(time, mass_avg_bottom, linestyle='-', linewidth=1, color='orange', aa=True, label='No Depuration (D)')
#######################################
print "Reading State" # min growth case, stationary
data = np.loadtxt('../200/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol, :]
print "Calculating Mean Mass"
min_growth =  np.mean(mass[:,:], axis=0)
#######################################
print "Reading State" # max growth case, stationary
data = np.loadtxt('../201/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol, :]
print "Calculating Mean Mass"
max_growth =  np.mean(mass[:,:], axis=0)
ax[0].fill_between(time, max_growth, min_growth, color='yellow', alpha=0.1)
#######################################

ax[0].set_xlabel(r'Days')
ax[0].xaxis.set_major_locator(MultipleLocator(5))
ax[0].xaxis.set_minor_locator(MultipleLocator(1))
ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax[0].set_xlim(0,30)
ax[0].xaxis.grid(False); 

ax[0].set_ylabel('Mass (g)', rotation=90)

ax[0].yaxis.set_major_locator(MultipleLocator(1))
ax[0].yaxis.set_minor_locator(MultipleLocator(0.5))
ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax[0].yaxis.grid(False);
ax[0].set_frame_on(True)
ax[0].set_ylim(11,23)

legend = ax[0].legend(loc='upper left')
legend.get_frame().set_facecolor('none')
legend.get_frame().set_edgecolor('none')
text = legend.get_texts()
text[0].set_color('black')
text[1].set_color('red')
text[2].set_color('green')
text[3].set_color('blue')

#plt.show()
plt.savefig('./fish_growth_color.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
