#!/usr/bin/python
from pylab import *
import matplotlib.tri as tri
import math

nplots=1
fontsize  = 10
linewidth = 1
uniformPadding = 0.15
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = 3.0
default_alpha = 1.0
lineRGBA = [1.0, 0.0, 0.0, default_alpha]

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



# Calculate Domain Data
vert_x, vert_y = np.loadtxt('../100/mesh_node.dat', skiprows=1, usecols=(1,2), unpack=True)
ind1, ind2, ind3 = np.loadtxt('../100/mesh_elem.dat', dtype='i8', skiprows=1, usecols=(1,2,3), unpack=True)
indices = (ind1, ind2, ind3)
xmax = max(vert_x[:])
ymax = max(vert_y[:])
xmin = min(vert_x[:])
ymin = min(vert_y[:])

# Figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.1, left=hpadding, right=1.0-hpadding/2., hspace=0.20)
ax = []
ax.append(fig.add_subplot(1,2,1))
ax.append(fig.add_subplot(1,2,2))

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
print "(A) Reading Position" # experiment A
data = np.loadtxt('../100/fish_position.dat', unpack=True)
for ii in range(0, nfish):
    xcol = ii*dwidth + 2
    ycol = ii*dwidth + 3
    xpos[ii,:] = data[xcol,:]
    ypos[ii,:] = data[ycol,:]

ax[0].scatter(xpos[0:24,mid], ypos[0:24,mid], s=25, color='black', zorder=10, edgecolors='black') # mid markers
ax[0].scatter(xpos[100:124,mid], ypos[100:124,mid], s=25, color='black', zorder=10, edgecolors='black') # mid markers
ax[1].scatter(xpos[0:24,end], ypos[0:24,end], s=25, color='black', zorder=10, edgecolors='black') # final markers
ax[1].scatter(xpos[100:124,end], ypos[100:124,end], s=25, color='black', zorder=10, edgecolors='black') # final markers

#######################################
print "(B) Reading Position" # experiment B
data = np.loadtxt('../101/fish_position.dat', unpack=True)
for ii in range(0, nfish):
    xcol = ii*dwidth + 2
    ycol = ii*dwidth + 3
    xpos[ii,:] = data[xcol,:]
    ypos[ii,:] = data[ycol,:]
    
ax[0].scatter(xpos[0:24,mid], ypos[0:24,mid], s=25, color='red', zorder=10, edgecolors='red') # mid markers
ax[0].scatter(xpos[100:124,mid], ypos[100:124,mid], s=25, color='red', zorder=10, edgecolors='red') # mid markers
ax[1].scatter(xpos[0:24,end], ypos[0:24,end], s=25, color='red', zorder=10, edgecolors='red') # final markers
ax[1].scatter(xpos[100:124,end], ypos[100:124,end], s=25, color='red', zorder=10, edgecolors='red') # final markers
    
#######################################
print "(C) Reading Position" # experiment C
data = np.loadtxt('../102/fish_position.dat', unpack=True)
for ii in range(0, nfish):
    xcol = ii*dwidth + 2
    ycol = ii*dwidth + 3
    xpos[ii,:] = data[xcol,:]
    ypos[ii,:] = data[ycol,:]
print "(C) Plotting Mid-experiment Positions"
ax[0].scatter(xpos[0:24,mid], ypos[0:24,mid], s=25, color='green', zorder=10, edgecolors='green') # mid markers
ax[0].scatter(xpos[100:124,mid], ypos[100:124,mid], s=25, color='green', zorder=10, edgecolors='green') # mid markers
print "(C) Plotting Final Positions"
ax[1].scatter(xpos[0:24,end], ypos[0:24,end], s=25, color='green', zorder=10, edgecolors='green') # final markers
ax[1].scatter(xpos[100:124,end], ypos[100:124,end], s=25, color='green', zorder=10, edgecolors='green') # final markers

#######################################
print "(D) Reading Position" # experiment D
data = np.loadtxt('../103/fish_position.dat', unpack=True)
for ii in range(0, nfish):
    xcol = ii*dwidth + 2
    ycol = ii*dwidth + 3
    xpos[ii,:] = data[xcol,:]
    ypos[ii,:] = data[ycol,:]
print "(D) Plotting Mid-experiment Positions"
ax[0].scatter(xpos[0:24,mid], ypos[0:24,mid], s=25, color='blue', zorder=10, edgecolors='blue') # mid markers
ax[0].scatter(xpos[100:124,mid], ypos[100:124,mid], s=25, color='blue', zorder=10, edgecolors='blue') # mid markers
print "(D) Plotting Final Positions"
ax[1].scatter(xpos[0:24,end], ypos[0:24,end], s=25, color='blue', zorder=10, edgecolors='blue') # final markers
ax[1].scatter(xpos[100:124,end], ypos[100:124,end], s=25, color='blue', zorder=10, edgecolors='blue') # final markers
    
    
ax[0].set_xlabel(r'X (m)')
ax[0].set_xlim(xmin, xmax)
ax[0].xaxis.set_major_locator(MultipleLocator(125))
ax[0].xaxis.set_minor_locator(MultipleLocator(25))
ax[0].xaxis.grid(True)

ax[0].set_ylabel(r'Y (m)')
ax[0].set_ylim(ymin, ymax)
ax[0].yaxis.set_major_locator(MultipleLocator(125))
ax[0].yaxis.set_minor_locator(MultipleLocator(25))
ax[0].yaxis.grid(False)

ax[1].set_xlabel(r'X (m)')
ax[1].set_xlim(xmin, xmax)
ax[1].xaxis.set_major_locator(MultipleLocator(125))
ax[1].xaxis.set_minor_locator(MultipleLocator(25))
ax[1].xaxis.grid(True)

#ax[1].set_ylabel()
ax[1].set_ylim(ymin, ymax)
ax[1].yaxis.set_major_locator(MultipleLocator(125))
ax[1].yaxis.set_minor_locator(MultipleLocator(25))
ax[1].yaxis.grid(False)

plt.savefig('./fish_position_final_color.tiff', edgecolor='none', dpi=default_dpi)

