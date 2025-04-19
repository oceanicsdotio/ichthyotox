#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from pylab import *
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi

nplots=1
fontsize  = 10
linewidth = 1
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = 3.5
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
path = zeros((nfish,end-start+1))
xpos = zeros((nfish,end-start+1))
ypos = zeros((nfish,end-start+1))
suitability = zeros((nfish,end-start+1))
ratio = zeros((nfish,end-start+1))
ratio_avg = zeros(end-start+1)
ratio_avg1 = zeros(end-start+1)
ratio_avg2 = zeros(end-start+1)
ratio_std = zeros(end-start+1)
ratio_sum = zeros(end-start+1)
suit_sum = zeros(end-start+1)
residence = zeros((nfish,end-start+1))
time = time / 24.
nsamples = 4
pindex = [0,99,100,199]
epad = 0.05
lindex = 3
threshold = 20.0
ii=199
duration=1



# figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.1, left=hpadding/2.0, right=1.0-(hpadding/2.0))
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
print "(A) Reading State" # experiment A
data = np.loadtxt('../100/fish_state.dat', unpack=True)
event_start = 19.5
masscol = ii*dwidth + 2
toxcol = ii*dwidth + 3
mass[ii,:] = data[masscol,:]
tox[ii,:] = data[toxcol,:]
print "(A) Calculating Toxin Load"
ratio[ii,:] = 1.0*((1000*tox[ii,:]/mass[ii,:])>0.005)
ax[0].plot(time[240*event_start:240*(event_start+duration)]-event_start, ratio[ii,240*event_start:240*(event_start+duration)]+lindex*(1+epad), linestyle='--', linewidth=1, color='black', aa=True, zorder=2)
#######################################
print "(A) Reading Position" # experiment A
data = np.loadtxt('../100/fish_position.dat', unpack=True)
xcol = ii*dwidth + 2
ycol = ii*dwidth + 3
xpos[ii,:] = data[xcol,:]
print "(A) Calculating Suitability and Residence Time"
suitability[ii,:] = 1.0*(0.5*(1.0 + sin(2.0*pi*(xpos[ii,:]-125.0)/500.0) ) > 0.5)
ax[0].fill_between(time[240*event_start:240*(event_start+duration)]-event_start, suitability[ii,240*event_start:240*(event_start+duration)]+lindex*(1+epad), lindex*(1+epad), facecolor=[0.0,0.0,0.0,0.5], edgecolor='none',  zorder=1)
for jj in range(start,end+1):
    for kk in range(jj+1,end+1):
        distance = sqrt((ypos[ii,jj] - ypos[ii,kk])**2 + (xpos[ii,jj] - xpos[ii,kk])**2)
        if (distance > threshold):
            next_crossing = kk-1
            break
        if (kk==end):
            next_crossing = kk
    residence[ii,jj] = (time[next_crossing] - time[jj])
    for kk in range(jj-1,start-1,-1):
        distance = sqrt((ypos[ii,jj] - ypos[ii,kk])**2 + (xpos[ii,jj] - xpos[ii,kk])**2)
        if (distance > threshold):
            next_crossing = kk+1
            break
        if (kk==start):
            next_crossing = kk
    residence[ii,jj] = residence[ii,jj] + (time[jj] - time[next_crossing])
ax[0].plot(time[240*event_start:240*(event_start+duration)]-event_start, residence[ii,240*event_start:240*(event_start+duration)]+lindex*(1+epad), linestyle='-', linewidth=1, color='black', aa=True, zorder=2)
lindex=lindex-1
print
print "    Max Res. Time = ", np.max(residence[ii,240*event_start:240*(event_start+duration)])
print "    Min Res. Time = ", np.min(residence[ii,240*event_start:240*(event_start+duration)])
print
#######################################
print "(B) Reading State" # experiment B
data = np.loadtxt('../101/fish_state.dat', unpack=True)
event_start = 28
masscol = ii*dwidth + 2
toxcol = ii*dwidth + 3
mass[ii,:] = data[masscol,:]
tox[ii,:] = data[toxcol,:]
print "(B) Calculating Toxin Load"
ratio = 1.0*((1000*tox/mass[ii,:])>0.005)
ax[0].plot(time[240*event_start:240*(event_start+duration)]-event_start, ratio[ii,240*event_start:240*(event_start+duration)]+lindex*(1+epad), linestyle='--', linewidth=1, color='red', aa=True, zorder=2)
#######################################
print "(B) Reading Position" # experiment B
data = np.loadtxt('../101/fish_position.dat', unpack=True)
xcol = ii*dwidth + 2
xpos[ii,:] = data[xcol,:]
print "(B) Calculating Suitability and Residence Time"
suitability[ii,:] = 1.0*(0.5*(1.0 + sin(2.0*pi*(xpos[ii,:]-125.0)/500.0) ) > 0.5)
ax[0].fill_between(time[240*event_start:240*(event_start+duration)]-event_start, suitability[ii,240*event_start:240*(event_start+duration)]+lindex*(1+epad), lindex*(1+epad), facecolor=[1.0,0.0,0.0,0.5], edgecolor='none',  zorder=1)
for jj in range(start,end+1):
    for kk in range(jj+1,end+1):
        distance = sqrt((ypos[ii,jj] - ypos[ii,kk])**2 + (xpos[ii,jj] - xpos[ii,kk])**2)
        if (distance > threshold):
            next_crossing = kk-1
            break
        if (kk==end):
            next_crossing = kk
    residence[ii,jj] = (time[next_crossing] - time[jj])
    for kk in range(jj-1,start-1,-1):
        distance = sqrt((ypos[ii,jj] - ypos[ii,kk])**2 + (xpos[ii,jj] - xpos[ii,kk])**2)
        if (distance > threshold):
            next_crossing = kk+1
            break
        if (kk==start):
            next_crossing = kk
    residence[ii,jj] = residence[ii,jj] + (time[jj] - time[next_crossing])
ax[0].plot(time[240*event_start:240*(event_start+duration)]-event_start, residence[ii,240*event_start:240*(event_start+duration)]+lindex*(1+epad), linestyle='-', linewidth=1, color='red', aa=True, zorder=2)
lindex=lindex-1
print
print "    Max Res. Time = ", np.max(residence[ii,240*event_start:240*(event_start+duration)])
print "    Min Res. Time = ", np.min(residence[ii,240*event_start:240*(event_start+duration)])
print
#######################################
print "(C) Reading State" # experiment C
data = np.loadtxt('../102/fish_state.dat', unpack=True)
event_start = 3
masscol = ii*dwidth + 2
toxcol = ii*dwidth + 3
mass[ii,:] = data[masscol,:]
tox[ii,:] = data[toxcol,:]
print "(C) Calculating Toxin Load"
ratio[ii,:] = 1.0*((1000*tox[ii,:]/mass[ii,:])>0.005)
ax[0].plot(time[240*event_start:240*(event_start+duration)]-event_start, ratio[ii,240*event_start:240*(event_start+duration)]+lindex*(1+epad), linestyle='--', linewidth=1, color='green', aa=True, zorder=2)
#######################################
print "(C) Reading Position" # experiment C
data = np.loadtxt('../102/fish_position.dat', unpack=True)
xcol = ii*dwidth + 2
xpos[ii,:] = data[xcol,:]
print "(C) Calculating Suitability and Residence Time"
suitability[ii,:] = 1.0*(0.5*(1.0 + sin(2.0*pi*(xpos[ii,:]-125.0)/500.0) ) > 0.5)
ax[0].fill_between(time[240*event_start:240*(event_start+duration)]-event_start, suitability[ii,240*event_start:240*(event_start+duration)]+lindex*(1+epad), lindex*(1+epad), facecolor=[0.0,1.0,0.0,0.3], edgecolor='none',  zorder=1)
for jj in range(start,end+1):
    for kk in range(jj+1,end+1):
        distance = sqrt((ypos[ii,jj] - ypos[ii,kk])**2 + (xpos[ii,jj] - xpos[ii,kk])**2)
        if (distance > threshold):
            next_crossing = kk-1
            break
        if (kk==end):
            next_crossing = kk
    residence[ii,jj] = (time[next_crossing] - time[jj])
    for kk in range(jj-1,start-1,-1):
        distance = sqrt((ypos[ii,jj] - ypos[ii,kk])**2 + (xpos[ii,jj] - xpos[ii,kk])**2)
        if (distance > threshold):
            next_crossing = kk+1
            break
        if (kk==start):
            next_crossing = kk
    residence[ii,jj] = residence[ii,jj] + (time[jj] - time[next_crossing])
    
ax[0].plot(time[240*event_start:240*(event_start+duration)]-event_start, residence[ii,240*event_start:240*(event_start+duration)]+lindex*(1+epad), linestyle='-', linewidth=1, color='green', aa=True, zorder=2)
lindex=lindex-1
print
print "    Max Res. Time = ", np.max(residence[ii,240*event_start:240*(event_start+duration)])
print "    Min Res. Time = ", np.min(residence[ii,240*event_start:240*(event_start+duration)])
print
#######################################
print "(D) Reading State" # experiment D
data = np.loadtxt('../103/fish_state.dat', unpack=True)
event_start = 10
masscol = ii*dwidth + 2
toxcol = ii*dwidth + 3
mass[ii,:] = data[masscol,:]
tox[ii,:] = data[toxcol,:]
print "(D) Calculating Toxin Load"
ratio[ii,:] = 1.0*((1000*tox[ii,:]/mass[ii,:])>0.005)
ax[0].plot(time[240*event_start:240*(event_start+duration)]-event_start, ratio[ii,240*event_start:240*(event_start+duration)]+lindex*(1+epad), linestyle='--', linewidth=1, color='blue', aa=True, zorder=2)
#######################################
print "(D) Reading Position" # experiment D
data = np.loadtxt('../103/fish_position.dat', unpack=True)
xcol = ii*dwidth + 2
xpos[ii,:] = data[xcol,:]
print "(D) Calculating Suitability and Residence Time"
suitability[ii,:] = 1.0*(0.5*(1.0 + sin(2.0*pi*(xpos[ii,:]-125.0)/500.0) ) > 0.5)
ax[0].fill_between(time[240*event_start:240*(event_start+duration)]-event_start, suitability[ii,240*event_start:240*(event_start+duration)]+lindex*(1+epad), lindex*(1+epad), facecolor=[0.0,0.0,1.0,0.5], edgecolor='none',  zorder=1)
for jj in range(start,end+1):
    for kk in range(jj+1,end+1):
        distance = sqrt((ypos[ii,jj] - ypos[ii,kk])**2 + (xpos[ii,jj] - xpos[ii,kk])**2)
        if (distance > threshold):
            next_crossing = kk-1
            break
        if (kk==end):
            next_crossing = kk
    residence[ii,jj] = (time[next_crossing] - time[jj])
    for kk in range(jj-1,start-1,-1):
        distance = sqrt((ypos[ii,jj] - ypos[ii,kk])**2 + (xpos[ii,jj] - xpos[ii,kk])**2)
        if (distance > threshold):
            next_crossing = kk+1
            break
        if (kk==start):
            next_crossing = kk
    residence[ii,jj] = residence[ii,jj] + (time[jj] - time[next_crossing])
ax[0].plot(time[240*event_start:240*(event_start+duration)]-event_start, residence[ii,240*event_start:240*(event_start+duration)]+lindex*(1+epad), linestyle='-', linewidth=1, color='blue', aa=True, zorder=2)
lindex=lindex-1
print
print "    Max Res. Time = ", np.max(residence[ii,240*event_start:240*(event_start+duration)])
print "    Min Res. Time = ", np.min(residence[ii,240*event_start:240*(event_start+duration)])
print
#######################################

ax[0].set_xlabel(r'Days')
ax[0].xaxis.set_major_locator(LinearLocator(5))
ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax[0].xaxis.set_minor_locator(LinearLocator(21))
ax[0].set_xlim(0,)
ax[0].xaxis.grid(False); 
ax[0].set_ylabel('Behavioral Cues and Residence Time', rotation=90)
ax[0].yaxis.set_major_locator(LinearLocator(0))
ax[0].yaxis.grid(False);
ax[0].set_frame_on(True)
ax[0].set_ylim(-0.1,4.0+0.25)

# legend = ax[0].legend(loc='upper left')
# legend.get_frame().set_facecolor('none')
# text = legend.get_texts()
# text[0].set_color('white')
# text[1].set_color('red')
# text[2].set_color('green')
# text[3].set_color('blue')

ax[0].text(x=0.05, y=3.75, s='Control (A)', color='black', fontsize=fontsize)
ax[0].text(x=0.05, y=2.75, s='Formation (B)', color='red', fontsize=fontsize)
ax[0].text(x=0.05, y=1.75, s='Intensification (C)', color='green', fontsize=fontsize)
ax[0].text(x=0.05, y=0.75, s='Decline (D)', color='blue', fontsize=fontsize)


plt.savefig('./fish_singlecues_color.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
