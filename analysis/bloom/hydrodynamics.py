#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi
from pylab import *

fontsize  = 10
linewidth = 1
uniformPadding = 0.1
marginWidth = 6.5
default_alpha = 1.0
lineRGBA = [1.0, 0.0, 0.0, default_alpha]
style = ['-',':']
div=100


fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, marginWidth/3.)) #Change this
fig.subplots_adjust(top=1.0-uniformPadding, bottom=uniformPadding, left=uniformPadding, right=1.0-uniformPadding)
   
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
     
  

radius = 10.0**-4.
temperature = linspace(20.0,30.0,div)
ratio = linspace(0.0,4.0,div)
salinity = 5.0

rho_pure = 999.842594 + 6.793952*10.0**(-2.0)*temperature - 9.095290*10.0**(-3.0)*temperature**(2.0) + 1.001685*10.0**(-4.0)*temperature**(3.0) - 1.120083*10.0**(-6.0)*temperature**(4.0) + 6.536332*10.0**(-9.0)*temperature**(5.0)


C = 0.824493 - 4.0899*10.0**(-3.0)*temperature + 7.6438*10.0**(-5.0)*temperature**(2.0) - 8.2467*10.0**(-7.0)*temperature**(3.0) + 5.3875*10.0**(-9.0)*temperature**(4.0)
D = -5.72466*10**(-3.0) + 1.0227*10.0**(-4.0)*temperature - 1.6546*10.0**(-6.0)*temperature**(2.0)
E = 4.8314*10.0**(-4.0)
density = rho_pure + salinity*C + salinity**(3.0/2.0)*D + E*salinity**2.0

A = 1.541 + 19.998*10.**(-2.)*temperature - 9.52*10.**(-5.)*temperature**(2.)
B = 7.974 - 7.561*10.**(-2.) + 4.724*10.**(-4.)*temperature**(2.)
visc_pure = 4.2844*10.**(-5.) + (0.157*(temperature+64.993)**(2.)-91.296)**(-1.)
visc = visc_pure*(1. + A*salinity + B*salinity**(2.)) 


temperature, rr1 = np.meshgrid(temperature, ratio, sparse=True)
density, rr2 = np.meshgrid(density, ratio, sparse=True)
visc, rr3 = np.meshgrid(visc, ratio, sparse=True)

velocity = zeros((div,div))
reynolds = zeros((div,div))

for ii in range(0,div): # X index value, temperature
    for jj in range(0,div): # Y index value, ratio
        velocity[ii,jj] = 1962. * radius**2. / visc[:,ii] * (density[:,ii] - 1072.1 - exp(4.644 - 0.7*rr1[jj]))
        reynolds[ii,jj] = 2. * radius * velocity[ii,jj] * density[:,ii] / visc[:,ii] / 3600.


# surface and contours
ax.view_init(elev = 45.0/2.0, azim=225.0-23.)
surf = ax.plot_surface(temperature, log(rr1), velocity, color=bg_color, edgecolor='gray', linewidth=1, antialiased=True, shade=False)

# axis adjustment
ax.set_xlabel(r'T')
ax.set_xlim(20.0,40.0)
ax.xaxis.set_major_locator(LinearLocator(2))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.xaxis.grid(False)
ax.xaxis.set_rotate_label(False)
ax.xaxis.set_pane_color(bg_color)

ax.set_ylabel(r'ln(C/M)')
ax.set_ylim()
ax.yaxis.set_major_locator(LinearLocator(2))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.yaxis.grid(False)
ax.yaxis.set_rotate_label(False)
ax.yaxis.set_pane_color(bg_color)

ax.set_zlabel(r'|V|')
ax.zaxis.set_major_locator(LinearLocator(2))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.3e'))
ax.zaxis.grid(False)
ax.zaxis.set_rotate_label(False)
ax.zaxis.set_pane_color(bg_color)

title(r'(A)', fontsize=fontsize, color=label_color)
ttl = ax.title
ttl.set_position([1.0, 1.0])


ax2 = fig.add_subplot(1,2,2,projection='3d')
ax2.view_init(elev = 45.0/2.0, azim=225.0-23.)
surf2 = ax2.plot_surface(temperature, log(rr1), reynolds, color=bg_color, edgecolor='gray', linewidth=1, antialiased=True, shade=False)

# axis adjustment
ax2.set_xlabel(r'T')
ax2.set_xlim()
ax2.xaxis.set_major_locator(LinearLocator(2))
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax2.xaxis.grid(False)
ax2.xaxis.set_rotate_label(False)
ax2.xaxis.set_pane_color(bg_color)

ax2.set_ylabel(r'ln(C/M)')
ax2.set_ylim()
ax2.yaxis.set_major_locator(LinearLocator(2))
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax2.yaxis.grid(False)
ax2.yaxis.set_rotate_label(False)
ax2.yaxis.set_pane_color(bg_color)

ax2.set_zlabel(r'Re')
ax2.zaxis.set_major_locator(LinearLocator(2))
ax2.zaxis.grid(False)
ax2.zaxis.set_rotate_label(False)
ax2.zaxis.set_pane_color(bg_color)

title(r"(B)", fontsize=fontsize, color=label_color)
ttl = ax2.title
ttl.set_position([1.0, 1.0])


plt.savefig('./figures/bloom/hydrodynamics.png', facecolor=bg_color, edgecolor='none', dpi=default_dpi)