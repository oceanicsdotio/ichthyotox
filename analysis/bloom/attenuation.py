#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
from pylab import *

fontsize  = 10
linewidth = 1
uniformPadding = 0.1
vpad = 0.1
hpad = 0.1
marginWidth = 7.0
default_alpha = 1.0
lineRGBA = [1.0, 0.0, 0.0, default_alpha]
style = ['-',':']

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
#rc('legend', fontsize=2*fontsize/3)
rc('xtick.major', pad=5)
rc('ytick.major', pad=5)

    
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, marginWidth/3.)) #Change this
fig.subplots_adjust(top=1.0-vpad, bottom=vpad, left=hpad/2., right=1.0-hpad/2.)
   
# set color scheme
ax = fig.add_subplot(1,2,1,projection='3d')

ax.patch.set_facecolor(bg_color)
ax.spines['top'].set_color(overlay_color)
ax.spines['bottom'].set_color(overlay_color)
ax.spines['left'].set_color(overlay_color)
ax.spines['right'].set_color(overlay_color)
ax.xaxis.label.set_color(label_color)
ax.yaxis.label.set_color(label_color)
ax.tick_params(axis='x', colors=label_color)
ax.tick_params(axis='y', colors=label_color)
     
   
def set_tick_sizes(ax, major, minor):
    for l in ax.get_xticklines() + ax.get_yticklines() + + ax.get_zticklines():
        l.set_markersize(major)
    for tick in ax.xaxis.get_minor_ticks() + ax.yaxis.get_minor_ticks() + ax.zaxis.get_minor_ticks():
        tick.tick1line.set_markersize(minor); tick.tick2line.set_markersize(minor); tick.tick3line.set_markersize(minor)
    ax.xaxis.LABELPAD      = 10.
    ax.xaxis.OFFSETTEXTPAD = 10.
    ax.yaxis.LABELPAD      = 10.
    ax.yaxis.OFFSETTEXTPAD = 10.
    ax.zaxis.LABELPAD      = 10.
    ax.zaxis.OFFSETTEXTPAD = 10.
    


# variable declarations
biomass = linspace(0.01,0.25,100)
time = linspace(0.0,24.0,100)
depth = linspace(-5.0,0.0,100)
depth2 = depth

time, depth = np.meshgrid(time, depth)
irradiance = (325. + 325.*cos(2.0*pi*time/12.)) * exp(0.15*depth)
irradiance[time[:]<6.0] = 0.0
irradiance[time[:]>18.0]=0.0
  

biomass, depth2 = np.meshgrid(biomass, depth2)
irradiance2 = (325. + 325.*cos(0.552*12.)) * exp(0.15*depth2 - 14.*biomass)

# surface and contours
ax.view_init(elev = 2.*45./3., azim=225.0)
surf = ax.plot_surface(time, depth, irradiance, rstride=5,cstride=5, color='white', edgecolors='black', linewidth=1.0, shade=False, antialiased=True)
cont1 = ax.contour(time, depth, irradiance, colors='black', extend3D=False, linewidth=2, zdir='x')


# axis adjustment
ax.set_xlabel(r't')
ax.set_xlim(0.0,24.0)
ax.xaxis.set_major_locator(LinearLocator(2))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.xaxis.grid(False)
ax.xaxis.set_rotate_label(False)
ax.w_xaxis.set_pane_color(bg_color)

ax.set_ylabel(r'Z')
ax.set_ylim(0.0,-5.0)
ax.yaxis.set_major_locator(LinearLocator(2))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.yaxis.grid(False)
ax.yaxis.set_rotate_label(False)
ax.w_yaxis.set_pane_color(bg_color)

ax.set_zlabel(r'I')
ax.set_zlim(0.0,650.0)
ax.zaxis.set_major_locator(LinearLocator(2))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.zaxis.grid(False)
ax.zaxis.set_rotate_label(False)
ax.w_zaxis.set_pane_color(bg_color)



title(r'(A)', fontsize=fontsize, color=label_color)
ttl = ax.title
ttl.set_position([1.0, 1.0])

ax2 = fig.add_subplot(1,2,2, projection='3d')
ax2.view_init(elev = 2.*45./3., azim=(225.+180.))
surf2 = ax2.plot_surface(log(biomass), depth2, irradiance2, rstride=5, cstride=5, color='white', edgecolors='black', linewidth=1.0, shade=False, antialiased=True)
cont2 = ax2.contour(log(biomass), depth2, irradiance2, colors=overlay_color, extend3D=False, linewidth=2, zdir='z')

# axis adjustment
ax2.set_xlabel(r'ln B')
ax2.set_xlim()
ax2.xaxis.set_major_locator(LinearLocator(2))
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax2.xaxis.grid(False)
ax2.xaxis.set_rotate_label(False)
ax2.w_xaxis.set_pane_color(bg_color)

ax2.set_ylabel('Z')
ax2.set_ylim()
ax2.yaxis.set_major_locator(LinearLocator(2))
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax2.yaxis.grid(False)
ax2.yaxis.set_rotate_label(False)
ax2.w_yaxis.set_pane_color(bg_color)

ax2.set_zlabel('I')
ax2.set_zlim(0.0,650.0)
ax2.zaxis.set_major_locator(LinearLocator(2))
ax2.zaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax2.zaxis.grid(False)
ax2.zaxis.set_rotate_label(False)
ax2.w_zaxis.set_pane_color(bg_color)

title(r'(B)', fontsize=fontsize, color=label_color)
ttl = ax2.title
ttl.set_position([1.0, 1.0])

#plt.show()
plt.savefig('./forcing_surfaces.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)