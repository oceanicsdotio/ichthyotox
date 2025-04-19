#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
from pylab import *

nplots=1
fontsize  = 10
linewidth = 1
hpadding = 0.1
vpadding = 0.05
marginWidth = 8.5
fheight = 6
default_alpha = 1.0

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
residence = zeros((nfish,end-start+1))
ratio = zeros((nfish,end-start+1))
ratio_avg = zeros(end-start+1)
ratio_avg1 = zeros(end-start+1)
ratio_avg2 = zeros(end-start+1)
ratio_std = zeros(end-start+1)
ratio_sum = zeros(end-start+1)
suit_sum = zeros(end-start+1)
time = time / 24.
nsamples = 4
pindex = [0,99,100,199]
epad = 0.1
lindex = 0
threshold = 20.0
ii=199
duration=1


# figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.05, left=hpadding, right=1.0-(hpadding))
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
print "Reading State" # experiment A
event_start = 19.5
data = np.loadtxt('../100/fish_position.dat', unpack=True)
xcol = ii*dwidth + 2
ycol = ii*dwidth + 3
xpos[ii,:] = data[xcol,:]
ypos[ii,:] = data[ycol,:]
print "Calculating Mean Tox Load and Pathway Partitioning"


for jj in range(start,end):
    for kk in range(jj+1,end):
        distance = sqrt((ypos[ii,jj] - ypos[ii,kk])**2 + (xpos[ii,jj] - xpos[ii,kk])**2)
        if (distance > threshold):
            next_crossing = kk-1
            break
        if (kk==end):
            next_crossing = kk
    residence[ii,jj] = (time[next_crossing] - time[jj])
    for kk in range(jj-1,start,-1):
        distance = sqrt((ypos[ii,jj] - ypos[ii,kk])**2 + (xpos[ii,jj] - xpos[ii,kk])**2)
        if (distance > threshold):
            next_crossing = kk+1
            break
        if (kk==start):
            next_crossing = kk
    residence[ii,jj] = residence[ii,jj] + (time[jj] - time[next_crossing])
ax[0].plot(time[240*event_start:240*(event_start+duration)]-event_start, residence[ii,240*event_start:240*(event_start+duration)], color='white', aa=True, linewidth=1)
#######################################
print "Reading State" # experiment B
event_start = 28
data = np.loadtxt('../101/fish_position.dat', unpack=True)
xcol = ii*dwidth + 2
ycol = ii*dwidth + 3
xpos[ii,:] = data[xcol,:]
ypos[ii,:] = data[ycol,:]
print "Calculating Mean Tox Load and Pathway Partitioning"
for jj in range(start,end):
    for kk in range(jj+1,end):
        distance = sqrt((ypos[ii,jj] - ypos[ii,kk])**2 + (xpos[ii,jj] - xpos[ii,kk])**2)
        if (distance > threshold):
            next_crossing = kk-1
            break
        if (kk==end):
            next_crossing = kk
    residence[ii,jj] = (time[next_crossing] - time[jj])
    for kk in range(jj-1,start,-1):
        distance = sqrt((ypos[ii,jj] - ypos[ii,kk])**2 + (xpos[ii,jj] - xpos[ii,kk])**2)
        if (distance > threshold):
            next_crossing = kk+1
            break
        if (kk==start):
            next_crossing = kk
    residence[ii,jj] = residence[ii,jj] + (time[jj] - time[next_crossing])
ax[0].plot(time[240*event_start:240*(event_start+duration)]-event_start, residence[ii,240*event_start:240*(event_start+duration)], color='red')

#######################################

#######################################
print "Reading State" # experiment C
event_start = 3
data = np.loadtxt('../102/fish_position.dat', unpack=True)
xcol = ii*dwidth + 2
ycol = ii*dwidth + 3
xpos[ii,:] = data[xcol,:]
ypos[ii,:] = data[ycol,:]
print "Calculating Mean Tox Load and Pathway Partitioning"

for jj in range(start,end):
    for kk in range(jj+1,end):
        distance = sqrt((ypos[ii,jj] - ypos[ii,kk])**2 + (xpos[ii,jj] - xpos[ii,kk])**2)
        if (distance > threshold):
            next_crossing = kk-1
            break
        if (kk==end):
            next_crossing = kk
    residence[ii,jj] = (time[next_crossing] - time[jj])
    for kk in range(jj-1,start,-1):
        distance = sqrt((ypos[ii,jj] - ypos[ii,kk])**2 + (xpos[ii,jj] - xpos[ii,kk])**2)
        if (distance > threshold):
            next_crossing = kk+1
            break
        if (kk==start):
            next_crossing = kk
    residence[ii,jj] = residence[ii,jj] + (time[jj] - time[next_crossing])

ax[0].plot(time[240*event_start:240*(event_start+duration)]-event_start, residence[ii,240*event_start:240*(event_start+duration)], color='green')

#######################################

#######################################
print "Reading State" # experiment D
event_start = 10
data = np.loadtxt('../103/fish_position.dat', unpack=True)
xcol = ii*dwidth + 2
ycol = ii*dwidth + 3
xpos[ii,:] = data[xcol,:]
ypos[ii,:] = data[ycol,:]
print "Calculating Mean Tox Load and Pathway Partitioning"


for jj in range(start,end):
    for kk in range(jj+1,end):
        distance = sqrt((ypos[ii,jj] - ypos[ii,kk])**2 + (xpos[ii,jj] - xpos[ii,kk])**2)
        if (distance > threshold):
            next_crossing = kk-1
            break
        if (kk==end):
            next_crossing = kk
    residence[ii,jj] = (time[next_crossing] - time[jj])
    for kk in range(jj-1,start,-1):
        distance = sqrt((ypos[ii,jj] - ypos[ii,kk])**2 + (xpos[ii,jj] - xpos[ii,kk])**2)
        if (distance > threshold):
            next_crossing = kk+1
            break
        if (kk==start):
            next_crossing = kk
    residence[ii,jj] = residence[ii,jj] + (time[jj] - time[next_crossing])
ax[0].plot(time[240*event_start:240*(event_start+duration)]-event_start, residence[ii,240*event_start:240*(event_start+duration)], color='blue')

#######################################

ax[0].set_xlabel(r'Days')
ax[0].xaxis.set_major_locator(LinearLocator(5))
ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax[0].set_xlim(0,)
ax[0].xaxis.grid(False); 
ax[0].set_ylabel('Residence time, days (with 20 meter threshold)', rotation=90)
ax[0].yaxis.set_major_locator(LinearLocator(5))
ax[0].yaxis.grid(False);
ax[0].set_frame_on(True)
ax[0].set_ylim()

# legend = ax[0].legend(loc='upper left')
# legend.get_frame().set_facecolor('none')
# text = legend.get_texts()
# text[0].set_color('white')
# text[1].set_color('red')
# text[2].set_color('green')
# text[3].set_color('blue')


plt.savefig('./fish_singlesuit_screen.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
