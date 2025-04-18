#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
from pylab import *

nplots=1
fontsize  = 10
linewidth = 1
uniformPadding = 0.15
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = 3.0
default_alpha = 1.0

for_screen=True
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
rc('font', **{'family':'sans-serif', 'sans-serif':['Arial']})
rc('mathtext', default='sf')
rc('lines', markeredgewidth=1)
rc('lines', linewidth=linewidth)
rc('axes', labelsize=fontsize)
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

data = zeros((ncolony*dwidth+1,end-start+1))
carbohydrate = zeros((ncolony,end-start+1))
protein = zeros((ncolony,end-start+1))
totalc = zeros((ncolony,end-start+1))
#tot_carb = zeros(end-start+1)
pro_avg = zeros(end-start+1)
carb_avg = zeros(end-start+1)
carb_std = zeros(end-start+1)
pro_std = zeros(end-start+1)
time = time / 24.
volume = 500.0*500.0*5.0


# figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.1, left=hpadding, right=1.0-(hpadding/2.), hspace=0.20)
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


print "Reading File 100" # experiment A
data = np.loadtxt('../100/cyanobacteria_state.dat', unpack=True)
for ii in range(0, ncolony):
    ccol = ii*dwidth + 2
    pcol = ii*dwidth + 3
    carbohydrate[ii,:] = data[ccol,:]
    protein[ii,:] = data[pcol,:]
print "Calculating Mean"
totalc = (protein)/volume
sum_avg = np.mean(totalc, axis=0)
sum_std = np.std(totalc, axis=0)
ax[0].plot(time, sum_avg, linestyle='-', linewidth=1, color=overlay_color, aa=True, label='Control (A)')
#ax[0].fill_between(time, (sum_avg+sum_std), (sum_avg-sum_std), facecolor=[1.0,1.0,1.0,0.25], edgecolor='none', zorder=3)
print "(A) Total Protein: ", sum(protein[:,end])
print "(A) Mean Biomass: ", sum_avg[end], "+/-", sum_std[end]
####################
print "Reading File 101" # experiment B
data = np.loadtxt('../101/cyanobacteria_state.dat', unpack=True)
for ii in range(0, ncolony):
    ccol = ii*dwidth + 2
    pcol = ii*dwidth + 3
    carbohydrate[ii,:] = data[ccol,:]
    protein[ii,:] = data[pcol,:]
print "Calculating Mean"
totalc = (protein)/volume
sum_avg = np.mean(totalc, axis=0)
sum_std = np.std(totalc, axis=0)
ax[0].plot(time, sum_avg, linestyle='-', linewidth=1, color='red', aa=True, label='Formation (B)')
#ax[0].fill_between(time, (sum_avg+sum_std), (sum_avg-sum_std), facecolor=[1.0,0.0,0.0,0.25], edgecolor='none', zorder=3)
print "(B) Total Protein: ", sum(protein[:,end])
print "(B) Mean Biomass: ", sum_avg[end], "+/-", sum_std[end]
####################
print "Reading File 102" # experiment C
data = np.loadtxt('../102/cyanobacteria_state.dat', unpack=True)
for ii in range(0, ncolony):
    ccol = ii*dwidth + 2
    pcol = ii*dwidth + 3
    carbohydrate[ii,:] = data[ccol,:]
    protein[ii,:] = data[pcol,:]
print "Calculating Mean"
totalc = (protein)/volume
sum_avg = np.mean(totalc, axis=0)
sum_std = np.std(totalc, axis=0)
ax[0].plot(time, sum_avg, linestyle='-', linewidth=1, color='green', aa=True, label='Intensification (C)')
#ax[0].fill_between(time, (sum_avg+sum_std), (sum_avg-sum_std), facecolor=[0.0,1.0,0.0,0.25], edgecolor='none', zorder=3)
print "(C) Total Protein: ", sum(protein[:,end])
print "(C) Mean Biomass: ", sum_avg[end], "+/-", sum_std[end]
####################
print "Reading File 103" # experiment D
data = np.loadtxt('../103/cyanobacteria_state.dat', unpack=True)
for ii in range(0, ncolony):
    ccol = ii*dwidth + 2
    pcol = ii*dwidth + 3
    carbohydrate[ii,:] = data[ccol,:]
    protein[ii,:] = data[pcol,:]
print "Calculating Mean"
totalc = (protein)/volume
sum_avg = np.mean(totalc, axis=0)
sum_std = np.std(totalc, axis=0)
ax[0].plot(time, sum_avg, linestyle='-', linewidth=1, color='blue', aa=True, label='Decline (D)')
#ax[0].fill_between(time, (sum_avg+sum_std), (sum_avg-sum_std), facecolor=[0.0,0.0,1.0,0.25], edgecolor='none', zorder=3)
print "(D) Total Protein: ", sum(protein[:,end])
print "(D) Mean Biomass: ", sum_avg[end], "+/-", sum_std[end]


ax[0].set_xlabel(r'Days')
ax[0].xaxis.set_major_locator(mticker.MultipleLocator(5))
ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax[0].set_xlim(time[start],time[end])
ax[0].xaxis.grid(False)
ax[0].set_ylabel(r'Mean Protein (g/m^3)', rotation=90)

ax[0].set_ylim(0,0.25)
ax[0].yaxis.set_major_locator(mticker.MultipleLocator(0.05))
ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax[0].yaxis.grid(False)
ax[0].set_frame_on(True)

legend = ax[0].legend(loc='center left')
legend.get_frame().set_facecolor('none')
legend.get_frame().set_edgecolor('none')
text = legend.get_texts()
text[0].set_color('white')
text[1].set_color('red')
text[2].set_color('green')
text[3].set_color('blue')
#plt.show()
plt.savefig('./cyano_carbon_protein_screen.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
