#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
from pylab import *

nplots=1
fontsize = 10
linewidth = 1
uniformPadding = 0.15
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
ini = open('../100/cyanobacteria_ini.dat', 'r')
ncolony = int(str.strip(ini.readline()))
time = np.loadtxt('../100/cyanobacteria_position.dat', usecols=[0], unpack=True)
start = 0
end = len(time) - 1
dwidth = 4


data = zeros((ncolony*dwidth+1, end-start+1))
cpsz = zeros((ncolony,end-start+1))
z_avg = zeros(end-start+1)
z_std = zeros(end-start+1)
time = time / 24.


# figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.1, left=hpadding, right=1.0-hpadding/2.)
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

#######################################
print "Reading File 100" # experiment A
data = np.loadtxt('../100/cyanobacteria_position.dat', unpack=True)
for ii in range(0, ncolony):
    zcol = ii*dwidth + 4
    cpsz[ii,:] = data[zcol,:]
print "Calculating Mean"
z_avg = np.mean(cpsz, axis=0)
z_std = np.std(cpsz, axis=0)
ax[0].fill_between(time, (z_avg-z_std), (z_avg+z_std), facecolor=[0.0,0.0,0.0,0.25], edgecolor='none', zorder = 1)
ax[0].plot(time, z_avg, linestyle='-', linewidth=1, color='black', aa=True, label='Control (A)', zorder = 5)
#######################################
print "Reading File 101" # experiment B
data = np.loadtxt('../101/cyanobacteria_position.dat', unpack=True)
for ii in range(0, ncolony):
    zcol = ii*dwidth + 4
    cpsz[ii,:] = data[zcol,:]
print "Calculating Mean"
z_avg = np.mean(cpsz, axis=0)
z_std = np.std(cpsz, axis=0)

ax[0].fill_between(time, (z_avg-z_std), (z_avg+z_std), facecolor=[1.0,0.0,0.0,0.25], edgecolor='none', zorder = 2)
ax[0].plot(time, z_avg, linestyle='-', linewidth=1, color='red', aa=True, label='Formation (B)', zorder = 6)
#######################################
print "Reading File 102" # experiment C
data = np.loadtxt('../102/cyanobacteria_position.dat', unpack=True)
for ii in range(0, ncolony):
    zcol = ii*dwidth + 4
    cpsz[ii,:] = data[zcol,:]
print "Calculating Mean"
z_avg = np.mean(cpsz, axis=0)
z_std = np.std(cpsz, axis=0)

ax[0].fill_between(time, (z_avg-z_std), (z_avg+z_std), facecolor=[0.0,1.0,0.0,0.25], edgecolor='none', zorder = 3)
ax[0].plot(time, z_avg, linestyle='-', linewidth=1, color='green', aa=True, label='Intensification (C)', zorder = 7)
#######################################
print "Reading File 103" # experiment D
data = np.loadtxt('../103/cyanobacteria_position.dat', unpack=True)
for ii in range(0, ncolony):
    zcol = ii*dwidth + 4
    cpsz[ii,:] = data[zcol,:]
print "Calculating Mean"
z_avg = np.mean(cpsz, axis=0)
z_std = np.std(cpsz, axis=0)

ax[0].fill_between(time, (z_avg-z_std), (z_avg+z_std), facecolor=[0.0,0.0,1.0,0.25], edgecolor='none', zorder = 4)
ax[0].plot(time, z_avg, linestyle='-', linewidth=1, color='blue', aa=True, label='Decline (D)', zorder = 8)
#######################################



ax[0].set_xlabel('Days')
ax[0].xaxis.set_major_locator(MultipleLocator(5))
ax[0].xaxis.set_minor_locator(MultipleLocator(1))
ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax[0].set_xlim(time[start],time[end])
ax[0].xaxis.grid(False); 

ax[0].set_ylabel('Depth (m)', rotation=90)
ax[0].yaxis.set_major_locator(MultipleLocator(1.0))
ax[0].yaxis.set_minor_locator(MultipleLocator(0.2))
ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax[0].yaxis.grid(False);
ax[0].set_ylim(-5.0,0.0)


ax[0].set_frame_on(True)

legend = ax[0].legend(loc='center right')
legend.get_frame().set_facecolor('none')
legend.get_frame().set_edgecolor('none')
text = legend.get_texts()
text[0].set_color('black')
text[1].set_color('red')
text[2].set_color('green')
text[3].set_color('blue')


plt.savefig('./cyano_paths_single.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
