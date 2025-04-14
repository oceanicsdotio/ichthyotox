#!/usr/bin/python
from pylab import *
import matplotlib.tri as tri

nplots=1
fontsize  = 10
linewidth = 1
uniformPadding = 0.15
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = marginWidth + 0.5
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


# Load fish initial positions
ini = open('../100/fish_ini.dat', 'r')
nfish = int(str.strip(ini.readline()))
fpsx = zeros(nfish)
fpsy = zeros(nfish)
fpsz = zeros(nfish)
fpsx[:], fpsy[:], fpsz[:] = np.loadtxt('../100/fish_ini.dat', skiprows=1, usecols=(1,2,3), unpack=True)

# Calculate Domain Data
vert_x, vert_y = np.loadtxt('../100/mesh_node.dat', skiprows=1, usecols=(1,2), unpack=True)
ind1, ind2, ind3 = np.loadtxt('../100/mesh_elem.dat', dtype='i8', skiprows=1, usecols=(1,2,3), unpack=True)
indices = (ind1, ind2, ind3)
xmax = max(vert_x[:])
ymax = max(vert_y[:])
xmin = min(vert_x[:])
ymin = min(vert_y[:])
triang = tri.Triangulation(vert_x, vert_y) # Delauney triangulation
suitability = 0.5*(1.0 + sin(2.0*pi*(vert_x-125.)/500. )) # calculate suitability

# Figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.05, left=hpadding, right=1.0-hpadding/2., hspace=0.20)
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

# Plot suitability
#plt.tripcolor(triang, suitability, shading='flat', cmap=plt.cm.Blues)
#plt.colorbar()

# Plot grid
print("Drawing grid...")
for ii in range(0,len(ind1)):
    plt.plot( (vert_x[ind1[ii]-1], vert_x[ind2[ii]-1]), (vert_y[ind1[ii]-1], vert_y[ind2[ii]-1]), color='grey', linewidth=0.5, aa=True, zorder=2)
    plt.plot( (vert_x[ind2[ii]-1], vert_x[ind3[ii]-1]), (vert_y[ind2[ii]-1], vert_y[ind3[ii]-1]), color='grey', linewidth=0.5, aa=True, zorder=2)
    plt.plot( (vert_x[ind3[ii]-1], vert_x[ind1[ii]-1]), (vert_y[ind3[ii]-1], vert_y[ind1[ii]-1]), color='grey', linewidth=0.5, aa=True, zorder=2)

# Plot path and end markers
print("Drawing paths...")
for ii in range(0,nfish):
    if (fpsz[ii] < 2.5):
        plt.scatter(fpsx[ii], fpsy[ii], s=40, color='red', zorder=10, edgecolors='face') # ini markers
    else:
        plt.scatter(fpsx[ii], fpsy[ii], s=40, color='blue', zorder=10, edgecolors='face')

ax[0].set_xlabel(r'X (m)')
ax[0].set_xlim(xmin, xmax)
ax[0].xaxis.set_major_locator(LinearLocator(2))
ax[0].xaxis.grid(False)

ax[0].set_ylabel(r'Y (m)')
ax[0].set_ylim(ymin, ymax)
ax[0].yaxis.set_major_locator(LinearLocator(2))
ax[0].yaxis.grid(False)


#title(r'(A)', fontsize=fontsize, color=label_color)
ttl = ax[0].title
ttl.set_position([1.0, 1.0])
#plt.savefig('./fish_paths.tiff', edgecolor='none', dpi=default_dpi)
plt.show()
