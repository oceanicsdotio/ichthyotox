#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, tri
from matplotlib.patches import Rectangle
from matplotlib.patches import Polygon


fontsize  = 10
linewidth = 1
uniformPadding = 0.15
hpadding = 0.1
vpadding = 0.1
marginWidth = 6.5
fheight = (marginWidth-0.5)/2
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
ini = open('../000/cyanobacteria_ini.dat', 'r')
ncolony = int(str.strip(ini.readline()))
time = np.loadtxt('../000/cyanobacteria_state.dat', usecols=[0], unpack=True)
start = 0
end = len(time) - 1
dwidth = 4

carbohydrate = zeros((ncolony,end-start+1))
protein = zeros((ncolony,end-start+1))
tot_carb = zeros(end-start+1)
pro_avg = zeros(end-start+1)
carb_avg = zeros(end-start+1)
time = time / 24.
volume = 500.0*500.0*5.0


# figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.05, left=hpadding, right=1.0-hpadding)
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


print "Reading File 000" # random / 10^-5
for ii in range(0, ncolony):
    ccol = ii*dwidth + 2
    pcol = ii*dwidth + 3
    carbohydrate[ii,:], protein[ii,:] = np.loadtxt('../000/cyanobacteria_state.dat', usecols=(ccol, pcol), unpack=True)
print "Calculating Mean"
pro_avg = np.mean(protein, axis=0)
carb_avg = np.mean(carbohydrate, axis=0)
tot_carb = float(ncolony)*(pro_avg + carb_avg)
ax[0].plot(time, tot_carb/volume, linestyle='-', color=overlay_color, aa=True)

print "Reading File 003" # random / 10^-3
for ii in range(0, ncolony):
    ccol = ii*dwidth + 2
    pcol = ii*dwidth + 3
    carbohydrate[ii,:], protein[ii,:] = np.loadtxt('../003/cyanobacteria_state.dat', usecols=(ccol, pcol), unpack=True)
print "Calculating Mean"
pro_avg = np.mean(protein, axis=0)
carb_avg = np.mean(carbohydrate, axis=0)
tot_carb = float(ncolony)*(pro_avg + carb_avg)
ax[0].plot(time, tot_carb/volume, linestyle='--', color=overlay_color, aa=True)

print "Reading File 004" # random / 10^-1
for ii in range(0, ncolony):
    ccol = ii*dwidth + 2
    pcol = ii*dwidth + 3
    carbohydrate[ii,:], protein[ii,:] = np.loadtxt('../004/cyanobacteria_state.dat', usecols=(ccol, pcol), unpack=True)
print "Calculating Mean"
pro_avg = np.mean(protein, axis=0)
carb_avg = np.mean(carbohydrate, axis=0)
tot_carb = float(ncolony)*(pro_avg + carb_avg)
ax[0].plot(time, tot_carb/volume, linestyle=':', color=overlay_color, aa=True)



sublabels = [[r'Days', r'Biomass (g/m^3)'], [r'Days', r'Carbon Ratio (g/g)']]
#subtitles = [r'(A)', r'(B)']
for ii in range(0,1):
    # axis adjustment
       

    ax[ii].set_xlabel(sublabels[ii][0])
    ax[ii].xaxis.set_major_locator(mticker.MultipleLocator(5)); 
    ax[ii].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))

        
    ax[ii].set_xlim(time[start],time[end])
    ax[ii].xaxis.grid(False); 
    
    ax[ii].set_ylabel(sublabels[ii][1], rotation=90); 
    ax[ii].set_ylim(0.0,15.0)
    ax[ii].yaxis.set_major_locator(mticker.MultipleLocator(5)); 
    ax[ii].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax[ii].yaxis.grid(False);
    
    #ax[ii].set_title(subtitles[ii], fontsize=fontsize, color=label_color)
    #ax[ii].title.set_position([0.5, 1.0])
    ax[ii].set_frame_on(True)

#plt.show()
plt.savefig('./cyano_carb_dependence.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
