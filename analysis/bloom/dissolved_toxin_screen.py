#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
from pylab import *

nplots=1
fontsize  = 10
hpadding = 0.15
vpadding = 0.1
marginWidth = 8.5
linewidth=1
fheight = 6.0
default_alpha = 1.0
show_grid=True

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
nlayers = 26
time = np.loadtxt('../102/dissolved_toxin.dat', usecols=[0], unpack=True)
start = 0
end = len(time) - 1

dissolved = zeros((nlayers,end-start+1))
avg_dissolved = zeros(end-start+1)
time = time / 24.
depth = linspace(0.0, -5.0, nlayers)

# figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.05, left=hpadding, right=1.0-hpadding/2.)
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


print "Reading File 012" # experiment C
data = np.loadtxt('../102/dissolved_toxin.dat', unpack=True)
for ii in range(0, nlayers):
    dissolved[ii,:] = data[ii+1,:]
print "Plotting"
dissolved = dissolved/(500.*500.)
avg_dissolved = np.mean(dissolved, axis=0)
ax[0].plot(log10(dissolved[:,240*3]), depth, linestyle='-', linewidth=1, color='green', aa=True, label='Intensification (C)')
ax[0].plot(log10(dissolved[:,end]), depth, linestyle='-', linewidth=2, color='green', aa=True)
ax[0].plot( (log10(avg_dissolved[end]), log10(avg_dissolved[end])), (-5.0, 0.0), linestyle='--', linewidth=1, color='green', aa=True)

print "Reading File 013" # experiment D
data = np.loadtxt('../103/dissolved_toxin.dat', unpack=True)
for ii in range(0, nlayers):
    dissolved[ii,:] = data[ii+1,:]
print "Plotting"
dissolved = dissolved/(500.*500.)
avg_dissolved = np.mean(dissolved, axis=0)
ax[0].plot(log10(dissolved[:,240*3]), depth, linestyle='-', linewidth=1, color='blue', aa=True, label='Decline (D)')
ax[0].plot( log10(dissolved[:,end]), depth, linestyle='-', linewidth=2, color='blue', aa=True)
ax[0].plot( (log10(avg_dissolved[end]), log10(avg_dissolved[end])), (-5.0, 0.0), linestyle='--', linewidth=1, color='blue', aa=True)



ax[0].set_xlabel('Log, base 10, of dissolved toxin concentration, g/m^3')
ax[0].xaxis.set_major_locator(MultipleLocator(0.5)); 
ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))


ax[0].xaxis.grid(False); 
ax[0].set_ylabel('Depth, meters', rotation=90); 
ax[0].set_ylim(-5.0,0.0)
ax[0].yaxis.set_major_locator(MultipleLocator(0.5)); 
ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax[0].yaxis.grid(False);
ax[0].set_frame_on(True)

legend = ax[0].legend(loc='upper left')
legend.get_frame().set_facecolor(bg_color)
text = legend.get_texts()
text[0].set_color('green')
text[1].set_color('blue')

plt.savefig('./dissolved_toxin_screen.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
