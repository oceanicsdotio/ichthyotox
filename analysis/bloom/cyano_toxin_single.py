#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
import matplotlib.patches as mpatches
from pylab import *

nplots=1
fontsize  = 10
linewidth = 1
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = 3
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
ini = open('../100/cyanobacteria_ini.dat', 'r')
ncolony = int(str.strip(ini.readline()))
time = np.loadtxt('../100/cyanobacteria_state.dat', usecols=[0], unpack=True)
start = 0
end = len(time) - 1
dwidth = 4

microcystin = zeros((ncolony, end-start+1))
carbohydrate = zeros((ncolony, end-start+1))
protein = zeros((ncolony, end-start+1))
data = zeros((ncolony*dwidth+1, end-start+1))
mc_avg = zeros(end-start+1)
mc_std = zeros(end-start+1)
time = time / 24.

# figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.1, left=hpadding, right=1.0-hpadding)
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


print "Reading File 101" # experiment B
data = np.loadtxt('../101/cyanobacteria_state.dat', unpack=True)
for ii in range(0, ncolony):
    ccol = ii*dwidth + 2
    pcol = ii*dwidth + 3
    mcol = ii*dwidth + 4
    carbohydrate[ii,:] = data[ccol, :]
    protein[ii,:] = data[pcol, :]
    microcystin[ii,:] = data[mcol, :]
print "Plotting"
mc_avg = np.mean(microcystin/(protein+carbohydrate), axis=0)
mc_std = np.std(microcystin/(protein+carbohydrate), axis=0)
ax[0].plot(time, mc_avg, linestyle='-', linewidth=1, color='red', aa=True, label='Formation (B)')
ax[0].fill_between(time, (mc_avg-mc_std), (mc_avg+mc_std), facecolor=[1.0,0.0,0.0,0.25], edgecolor='none', zorder = 2)

print "Reading File 102" # experiment C
data = np.loadtxt('../102/cyanobacteria_state.dat', unpack=True)
for ii in range(0, ncolony):
    ccol = ii*dwidth + 2
    pcol = ii*dwidth + 3
    mcol = ii*dwidth + 4
    carbohydrate[ii,:] = data[ccol, :]
    protein[ii,:] = data[pcol, :]
    microcystin[ii,:] = data[mcol, :]
print "Plotting"
mc_avg = np.mean(microcystin/(protein+carbohydrate), axis=0)
mc_std = np.std(microcystin/(protein+carbohydrate), axis=0)
ax[0].plot(time, mc_avg, linestyle='-', linewidth=1, color='green', aa=True, label='Intensification (C)')
ax[0].fill_between(time, (mc_avg-mc_std), (mc_avg+mc_std), facecolor=[0.0,1.0,0.0,0.25], edgecolor='none', zorder = 2)

print "Reading File 103" # experiment D
data = np.loadtxt('../103/cyanobacteria_state.dat', unpack=True)
for ii in range(0, ncolony):
    ccol = ii*dwidth + 2
    pcol = ii*dwidth + 3
    mcol = ii*dwidth + 4
    carbohydrate[ii,:] = data[ccol, :]
    protein[ii,:] = data[pcol, :]
    microcystin[ii,:] = data[mcol, :]
print "Plotting"
mc_avg = np.mean(microcystin/(protein+carbohydrate), axis=0)
mc_std = np.std(microcystin/(protein+carbohydrate), axis=0)
ax[0].plot(time, mc_avg, linestyle='-', linewidth=1, color='blue', aa=True, label='Decline (D)')
ax[0].fill_between(time, (mc_avg-mc_std), (mc_avg+mc_std), facecolor=[0.0,0.0,1.0,0.25], edgecolor='none', zorder = 2)


ax[0].set_xlabel(r'Days')
ax[0].xaxis.set_major_locator(mticker.MultipleLocator(5))
ax[0].xaxis.set_minor_locator(mticker.MultipleLocator(1))
ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax[0].set_xlim(time[start],time[end])
ax[0].xaxis.grid(False); 
ax[0].set_ylabel(r'Toxicity (g/g)', rotation=90); 
ax[0].set_ylim(0,)
ax[0].yaxis.set_major_locator(MultipleLocator(0.1)); 
ax[0].yaxis.set_minor_locator(MultipleLocator(0.02)); 
ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax[0].yaxis.grid(False)
ax[0].set_frame_on(True)


legend = ax[0].legend(loc='upper left')
legend.get_frame().set_facecolor('none')
legend.get_frame().set_edgecolor('none')
text = legend.get_texts()
text[0].set_color('red')
text[1].set_color('green')
text[2].set_color('blue')

#plt.show()
plt.savefig('./cyano_toxin_single.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
