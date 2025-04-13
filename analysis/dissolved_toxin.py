#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
from pylab import *

nplots=2
fontsize  = 10
linewidth = 1
uniformPadding = 0.15
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = 2.*(marginWidth + 0.5)/4.0
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
nlayers = 26
time = np.loadtxt('../012/dissolved_toxin.dat', usecols=[0], unpack=True)
start = 0
end = len(time) - 1

dissolved = zeros((nlayers,end-start+1))
avg_dissolved = zeros(end-start+1)
time = time / 24.
volume = 500.0*500.0*5.0
depth = linspace(0.0, -5.0, nlayers)

# figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding-0.05, bottom=vpadding+0.1, left=hpadding, right=1.0-hpadding/2., hspace=0.2)
ax = []
ax.append(fig.add_subplot(2,1,1))
ax.append(fig.add_subplot(2,1,2))

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
print "Reading File 012" # experiment C
data = np.loadtxt('../012/dissolved_toxin.dat', unpack=True)
for ii in range(0, nlayers):
    dissolved[ii,:] = data[ii+1,:]
print "Plotting"
dissolved = dissolved/(500.*500.)
avg_dissolved = np.mean(dissolved, axis=0)
ax[0].plot((dissolved[:,240*3]), depth, linestyle='-', linewidth=1, color=overlay_color, aa=True)
ax[0].plot((dissolved[:,end]), depth, linestyle='-', linewidth=2, color=overlay_color, aa=True)
ax[0].plot( ((avg_dissolved[end]), (avg_dissolved[end])), (-5.0, 0.0), linestyle='--', linewidth=1, color=overlay_color, aa=True)

print "Mean Dissolved C: ", avg_dissolved[end]*500.0*500.0
#######################################
print "Reading File 013" # experiment D
data = np.loadtxt('../013/dissolved_toxin.dat', unpack=True)
for ii in range(0, nlayers):
    dissolved[ii,:] = data[ii+1,:]
print "Plotting"
dissolved = dissolved/(500.*500.)
avg_dissolved = np.mean(dissolved, axis=0)
ax[1].plot((dissolved[:,240*3]), depth, linestyle='-', linewidth=1, color=overlay_color, aa=True)
ax[1].plot( (dissolved[:,end]), depth, linestyle='-', linewidth=2, color=overlay_color, aa=True)
ax[1].plot( ((avg_dissolved[end]), (avg_dissolved[end])), (-5.0, 0.0), linestyle='--', linewidth=1, color=overlay_color, aa=True)

print "Mean Dissolved D: ", avg_dissolved[end]*500.0*500.0
#######################################

sublabels = [[r'Dissolved Toxin (g/m^3)', r'Depth (m)'],[r'Dissolved Toxin (g/m^3)', r'Depth (m)']]
subtitles = [r'(A)', r'(B)']
for ii in range(0,nplots):
    # axis adjustment
    if (ii == nplots-1):
      ax[ii].set_xlabel(sublabels[ii][0])
      ax[ii].xaxis.set_major_locator(MultipleLocator(0.05)); 
      ax[ii].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    else:
        ax[ii].set_xlabel(' ')
        ax[ii].xaxis.set_major_locator(LinearLocator(0)); 
    
    ax[ii].xaxis.grid(False); 
    ax[ii].set_ylabel(sublabels[ii][1], rotation=90); 
    ax[ii].set_ylim(-5.0,0.0)
    ax[ii].yaxis.set_major_locator(LinearLocator(6)); 
    ax[ii].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax[ii].yaxis.grid(False);
    
    ax[ii].set_title(subtitles[ii], fontsize=fontsize, color=label_color)
    ax[ii].title.set_position([1.0, 1.0])
    ax[ii].set_frame_on(True)
    
ax[0].set_xlim(0.0,0.25)
ax[1].set_xlim(0.0,0.25)

#plt.show()
plt.savefig('./dissolved_toxin.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
