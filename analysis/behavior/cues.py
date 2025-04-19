#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
from pylab import *
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi

nplots=1
fontsize  = 10
linewidth = 1
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = 7.0
default_alpha = 1.0


# load data
ini = open('../100/fish_ini.dat', 'r')
nfish = int(str.strip(ini.readline()))
time = np.loadtxt('../100/fish_state.dat', usecols=[0], unpack=True)
start = 0
end = len(time) - 1
dwidth = 4

data = zeros((nfish*dwidth+1,end-start+1))
mass = zeros((nfish,end-start+1))
tox = zeros((nfish,end-start+1))
xpos = zeros((nfish,end-start+1))
suitability = zeros((nfish,end-start+1))
ratio = zeros((nfish,end-start+1))
ratio_sum = zeros(end-start+1)
suit_sum = zeros(end-start+1)
time = time / 24.
nsamples = 4
pindex = [0,99,100,199]
epad = 0.1
lindex = 19



# figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.05, left=hpadding/2.0, right=1.0-(hpadding/2.0))
ax = []
ax.append(fig.add_subplot(1,1,1))
#ax.append(ax[0].twinx())

ax[0].patch.set_facecolor(bg_color)
ax[0].spines['top'].set_color(overlay_color)
ax[0].spines['bottom'].set_color(overlay_color)
ax[0].spines['left'].set_color(overlay_color)
ax[0].spines['right'].set_color(overlay_color)
ax[0].xaxis.label.set_color(label_color)
ax[0].yaxis.label.set_color(label_color)
ax[0].tick_params(axis='x', colors=label_color)
ax[0].tick_params(axis='y', colors=label_color)

# ax[1].xaxis.label.set_color(label_color)
# ax[1].yaxis.label.set_color(label_color)
# ax[1].tick_params(axis='y', colors=label_color)

#######################################
print("(A) Reading State") # experiment A
data = np.loadtxt('../100/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    toxcol = ii*dwidth + 3
    mass[ii,:] = data[masscol,:]
    tox[ii,:] = data[toxcol,:]
print("(A) Calculating Mean Tox Load")
ratio = 1.0*((1000*tox/mass)>0.005)
ratio_sum = np.sum(ratio, axis=0)/200.
ax[0].plot(time, ratio_sum[:]+lindex*(1+epad), linestyle='-', linewidth=1, color='black', aa=True, label='Control (A)', zorder=2)
lindex=lindex-1
for ii in range(0,nsamples):
    ax[0].plot(time, ratio[pindex[ii],:]+lindex*(1+epad), linestyle='-', linewidth=1, color='black', aa=True, zorder=2)
    lindex=lindex-1
#######################################
print("(A) Reading Position") # experiment A
data = np.loadtxt('../100/fish_position.dat', unpack=True)
for ii in range(0, nfish):
    xcol = ii*dwidth + 2
    xpos[ii,:] = data[xcol,:]
print ("(A) Calculating Mean Tox Load")
suitability = (0.5*(1.0 + sin(2.0*pi*(xpos-125.0)/500.0) ) > 0.5)
suit_sum = np.sum(suitability, axis=0)/200.
lindex = lindex+5
ax[0].fill_between(time, suit_sum[:]+lindex*(1+epad), lindex*(1+epad), facecolor=[0.0,0.0,0.0,0.5], edgecolor='none',  zorder=1)
lindex=lindex-1
for ii in range(0,nsamples):
    ax[0].fill_between(time, suitability[pindex[ii],:]+lindex*(1+epad), lindex*(1+epad), facecolor=[0.0,0.0,0.0,0.5], edgecolor='none',  zorder=1)
    lindex=lindex-1
#######################################
print ("(B) Reading State") # experiment B
data = np.loadtxt('../101/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    toxcol = ii*dwidth + 3
    mass[ii,:] = data[masscol,:]
    tox[ii,:] = data[toxcol,:]
print ("(B) Calculating Mean Tox Load")
ratio = 1.0*((1000*tox/mass)>0.005)
ratio_sum = np.sum(ratio, axis=0)/200.
ax[0].plot(time, ratio_sum[:]+lindex*(1+epad), linestyle='-', linewidth=1, color='red', aa=True, label='Formation (B)', zorder=2)
lindex=lindex-1
for ii in range(0,nsamples):
    ax[0].plot(time, ratio[pindex[ii],:]+lindex*(1+epad), linestyle='-', linewidth=1, color='red', aa=True, zorder=2)
    lindex=lindex-1
#######################################
print ("(B) Reading Position") # experiment B
data = np.loadtxt('../101/fish_position.dat', unpack=True)
for ii in range(0, nfish):
    xcol = ii*dwidth + 2
    xpos[ii,:] = data[xcol,:]
print ("(B) Calculating Mean Tox Load")
suitability = (0.5*(1.0 + sin(2.0*pi*(xpos-125.0)/500.0) ) > 0.5)
suit_sum = np.sum(suitability, axis=0)/200.
lindex = lindex+5
ax[0].fill_between(time, suit_sum[:]+lindex*(1+epad), lindex*(1+epad), facecolor=[1.0,0.0,0.0,0.5], edgecolor='none',  zorder=1)
lindex=lindex-1
for ii in range(0,nsamples):
    ax[0].fill_between(time, suitability[pindex[ii],:]+lindex*(1+epad), lindex*(1+epad), facecolor=[1.0,0.0,0.0,0.5], edgecolor='none',  zorder=1)
    lindex=lindex-1
#######################################
print ("(C) Reading State") # experiment C
data = np.loadtxt('../102/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    toxcol = ii*dwidth + 3
    mass[ii,:] = data[masscol,:]
    tox[ii,:] = data[toxcol,:]
print ("(C) Calculating Mean Tox Load")
ratio = 1.0*((1000*tox/mass)>0.005)
ratio_sum = np.sum(ratio, axis=0)/200.
ax[0].plot(time, ratio_sum[:]+lindex*(1+epad), linestyle='-', linewidth=1, color='green', aa=True, label='Intensification (C)', zorder=2)
lindex=lindex-1
for ii in range(0,nsamples):
    ax[0].plot(time, ratio[pindex[ii],:]+lindex*(1+epad), linestyle='-', linewidth=1, color='green', aa=True, zorder=2)
    lindex=lindex-1
#######################################
print ("(C) Reading Position") # experiment C
data = np.loadtxt('../102/fish_position.dat', unpack=True)
for ii in range(0, nfish):
    xcol = ii*dwidth + 2
    xpos[ii,:] = data[xcol,:]
print ("(C) Calculating Mean Tox Load")
suitability = (0.5*(1.0 + sin(2.0*pi*(xpos-125.0)/500.0) ) > 0.5)
suit_sum = np.sum(suitability, axis=0)/200.
lindex = lindex+5
ax[0].fill_between(time, suit_sum[:]+lindex*(1+epad), lindex*(1+epad), facecolor=[0.0,1.0,0.0,0.3], edgecolor='none',  zorder=1)
lindex=lindex-1
for ii in range(0,nsamples):
    ax[0].fill_between(time, suitability[pindex[ii],:]+lindex*(1+epad), lindex*(1+epad), facecolor=[0.0,1.0,0.0,0.3], edgecolor='none',  zorder=1)
    lindex=lindex-1
#######################################
print ("(D) Reading State") # experiment D
data = np.loadtxt('../103/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    toxcol = ii*dwidth + 3
    mass[ii,:] = data[masscol,:]
    tox[ii,:] = data[toxcol,:]
print ("(D) Calculating Mean Tox Load")
ratio = 1.0*((1000*tox/mass)>0.005)
ratio_sum = np.sum(ratio, axis=0)/200.
ax[0].plot(time, ratio_sum[:]+lindex*(1+epad), linestyle='-', linewidth=1, color='blue', aa=True, label='Decline (D)', zorder=2)
lindex=lindex-1
for ii in range(0,nsamples):
    ax[0].plot(time, ratio[pindex[ii],:]+lindex*(1+epad), linestyle='-', linewidth=1, color='blue', aa=True, zorder=2)
    lindex=lindex-1
#######################################
print ("(D) Reading Position") # experiment D
data = np.loadtxt('../103/fish_position.dat', unpack=True)
for ii in range(0, nfish):
    xcol = ii*dwidth + 2
    xpos[ii,:] = data[xcol,:]
print ("(D) Calculating Mean Tox Load")
suitability = (0.5*(1.0 + sin(2.0*pi*(xpos-125.0)/500.0) ) > 0.5)
suit_sum = np.sum(suitability, axis=0)/200.
lindex = lindex+5
ax[0].fill_between(time, suit_sum[:]+lindex*(1+epad), lindex*(1+epad), facecolor=[0.0,0.0,1.0,0.5], edgecolor='none',  zorder=1)
lindex=lindex-1
for ii in range(0,nsamples):
    ax[0].fill_between(time, suitability[pindex[ii],:]+lindex*(1+epad), lindex*(1+epad), facecolor=[0.0,0.0,1.0,0.5], edgecolor='none',  zorder=1)
    lindex=lindex-1
#######################################

ax[0].set_xlabel(r'Days')
ax[0].xaxis.set_major_locator(MultipleLocator(5))
ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax[0].xaxis.set_minor_locator(MultipleLocator(1))
ax[0].set_xlim(0,30)
ax[0].xaxis.grid(False)
ax[0].set_ylabel('Behavioral Cues', rotation=90)
ax[0].yaxis.set_major_locator(LinearLocator(0))
ax[0].yaxis.grid(False)
ax[0].set_frame_on(True)
ax[0].set_ylim(0,22)

ax[0].text(x=0.05*30, y=21.25, s='Control (A)', color='black', fontsize=fontsize)
ax[0].text(x=0.05*30, y=15.75, s='Formation (B)', color='red', fontsize=fontsize)
ax[0].text(x=0.05*30, y=10.25, s='Intensification (C)', color='green', fontsize=fontsize)
ax[0].text(x=0.05*30, y=4.75, s='Decline (D)', color='blue', fontsize=fontsize)


if __name__ == "__main__":
    plt.savefig('./fish_cues_color.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
