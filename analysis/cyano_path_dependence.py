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
fheight = (marginWidth-0.5)
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
ndarks = int(time[end]/24)+1
dwidth = 4

cpsz = zeros((ncolony, end-start+1))
z_avg = zeros((ncolony, end-start+1))
time = time / 24.




# figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding, left=hpadding, right=1.0-hpadding)
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


print "Reading File 000" # random / 10^-5
for ii in range(0, ncolony):
    zcol = ii*dwidth + 4
    cpsz[ii,:] = np.loadtxt('../000/cyanobacteria_position.dat', usecols=[zcol], unpack=True)
print "Calculating Mean"
z_avg = np.mean(cpsz, axis=0)
ax[0].plot(time, z_avg, linestyle='-', color=overlay_color, aa=True)

print "Reading File 003" # random / 10^-3
for ii in range(0, ncolony):
    zcol = ii*dwidth + 4
    cpsz[ii,:] = np.loadtxt('../003/cyanobacteria_position.dat', usecols=[zcol], unpack=True)
print "Calculating Mean"
z_avg = np.mean(cpsz, axis=0)
ax[0].plot(time, z_avg, linestyle='--', color=overlay_color, aa=True)

print "Reading File 004" # random / 10^-1
for ii in range(0, ncolony):
    zcol = ii*dwidth + 4
    cpsz[ii,:] = np.loadtxt('../004/cyanobacteria_position.dat', usecols=[zcol], unpack=True)
print "Calculating Mean"
z_avg = np.mean(cpsz, axis=0)
ax[0].plot(time, z_avg, linestyle=':', color=overlay_color, aa=True)

print "Reading File 005" # sediment / 10^-5
for ii in range(0, ncolony):
    zcol = ii*dwidth + 4
    cpsz[ii,:] = np.loadtxt('../005/cyanobacteria_position.dat', usecols=[zcol], unpack=True)
print "Calculating Mean"
z_avg = np.mean(cpsz, axis=0)
ax[1].plot(time, z_avg, linestyle='-', color=overlay_color, aa=True)

print "Reading File 006" # sediment / 10^-3
for ii in range(0, ncolony):
    zcol = ii*dwidth + 4
    cpsz[ii,:] = np.loadtxt('../006/cyanobacteria_position.dat', usecols=[zcol], unpack=True)
print "Calculating Mean"
z_avg = np.mean(cpsz, axis=0)
ax[1].plot(time, z_avg, linestyle='--', color=overlay_color, aa=True)

print "Reading File 007" # sediment / 10^-1
for ii in range(0, ncolony):
    zcol = ii*dwidth + 4
    cpsz[ii,:] = np.loadtxt('../007/cyanobacteria_position.dat', usecols=[zcol], unpack=True)
print "Calculating Mean"
z_avg = np.mean(cpsz, axis=0)
ax[1].plot(time, z_avg, linestyle=':', color=overlay_color, aa=True)

print "Reading File 008" # surface / 10^-5
for ii in range(0, ncolony):
    zcol = ii*dwidth + 4
    cpsz[ii,:] = np.loadtxt('../008/cyanobacteria_position.dat', usecols=[zcol], unpack=True)
print "Calculating Mean"
z_avg = np.mean(cpsz, axis=0)
ax[2].plot(time, z_avg, linestyle='-', color=overlay_color, aa=True)

print "Reading File 009" # surface / 10^-3
for ii in range(0, ncolony):
    zcol = ii*dwidth + 4
    cpsz[ii,:] = np.loadtxt('../009/cyanobacteria_position.dat', usecols=[zcol], unpack=True)
print "Calculating Mean"
z_avg = np.mean(cpsz, axis=0)
ax[2].plot(time, z_avg, linestyle='--', color=overlay_color, aa=True)

print "Reading File 010" # surface / 10^-1
for ii in range(0, ncolony):
    zcol = ii*dwidth + 4
    cpsz[ii,:] = np.loadtxt('../010/cyanobacteria_position.dat', usecols=[zcol], unpack=True)
print "Calculating Mean"
z_avg = np.mean(cpsz, axis=0)
ax[2].plot(time, z_avg, linestyle=':', color=overlay_color, aa=True)



sublabels = [[r'Days', r'Depth (m)'], [r'Days', r'Depth (m)'], [r'Days', r'Depth (m)']]
subtitles = [r'(A)', r'(B)', r'(C)']
for ii in range(0,3):
    # axis adjustment
       
    if (ii == 2): 
        ax[ii].set_xlabel(sublabels[ii][0])
        ax[ii].xaxis.set_major_locator(mticker.MultipleLocator(5)); 
        ax[ii].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    else:
        ax[ii].set_xlabel(' ')
        ax[ii].xaxis.set_major_locator(LinearLocator(0));
        
    ax[ii].set_xlim(time[start],time[end])
    ax[ii].xaxis.grid(False); 
    
    ax[ii].set_ylabel(sublabels[ii][1], rotation=90); 
    ax[ii].set_ylim(-5.0, 0.0)
    ax[ii].yaxis.set_major_locator(LinearLocator(2)); 
    ax[ii].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax[ii].yaxis.grid(False);
    
    ax[ii].set_title(subtitles[ii], fontsize=fontsize, color=label_color)
    ax[ii].title.set_position([0.5, 1.0])
    ax[ii].set_frame_on(True)

#plt.show()
plt.savefig('./cyano_path_dependence.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
