#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
from pylab import *
from matplotlib import cm, tri

nplots=1
fontsize  = 10
linewidth = 1
hpadding = 0.1
vpadding = 0.05
marginWidth = 8.5
fheight = 6.0
default_alpha=1.0

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
rc('xtick.major', pad=5)
rc('ytick.major', pad=5)

# load data
ini = open('../100/fish_ini.dat', 'r')
nfish = int(str.strip(ini.readline()))
time = np.loadtxt('../100/fish_position.dat', usecols=[0], unpack=True)
time = time/24.0
start = 0
end = len(time) - 1
dwidth = 4
data = zeros((nfish*dwidth+1,end-start+1))
fpsx = zeros((nfish,end-start+1))
fpsy = zeros((nfish,end-start+1))
dx = zeros((nfish,end-start+1))
dy = zeros((nfish,end-start+1))
disp_x = zeros(nfish)
disp_y = zeros(nfish)
displacement = zeros((nfish,end-start+1))
msd_std_a = zeros((5,2))
msd_value_a = zeros((5,2))
msd_std_b = zeros((5,2))
msd_value_b = zeros((5,2))
msd_std_c = zeros((5,2))
msd_value_c = zeros((5,2))
msd_std_d = zeros((5,2))
msd_value_d = zeros((5,2))
residence = zeros((nfish,end-start+1))

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


# experiment phases (days)
phase_start = [0, 5, 10, 15, 20]
phase_end = [10, 15, 20, 25, 30]
window = 10 # steps = 1 hour
nphases = 5
width = 0.2

# load experiment data
data = np.loadtxt('../100/fish_position.dat', unpack=True) # Experiment A
for ii in range(0, nfish):
    xcol = ii*dwidth + 2
    ycol = ii*dwidth + 3
    fpsx[ii,:] = data[xcol,:]
    fpsy[ii,:] = data[ycol,:]
# calculate individual displacement over one hour time window
dx[:,1:end] = fpsx[:,1:end] - fpsx[:,0:end-1]
dy[:,1:end] = fpsy[:,1:end] - fpsy[:,0:end-1]
dx = dx - 500.0*(dx > 200)
dx = dx + 500.0*(dx < -200)
dy = dy - 500.0*(dy > 200)
dy = dy + 500.0*(dy < -200)
displacement[:,:] = 0.0
for tt in range(start+window, end): # loop of values to calculate
    disp_x[:] = 0.0
    disp_y[:] = 0.0
    for uu in range(tt-window, tt): # loop of steps in calc (10)
        disp_x[:] = disp_x[:] + dx[:,uu]
        disp_y[:] = disp_y[:] + dy[:,uu]
    displacement[:,tt] = displacement[:,tt] + disp_x**2 + disp_y**2 # displacement at end of one hour
    displacement[:,tt] = 0.1*displacement[:,tt] # average displace of one step, for each individual
# calculate summary stats
for ii in range(0, nphases):
    msd_std_a[ii,0] = np.std(displacement[:,phase_start[ii]*240:phase_end[ii]*240])
    msd_value_a[ii,0] = np.mean(displacement[:,phase_start[ii]*240:phase_end[ii]*240])
rects_a = ax[0].bar(np.arange(5), msd_value_a[:,0], width, color='white', yerr=msd_std_a[:,0], ecolor='white')


# load experiment data
data = np.loadtxt('../101/fish_position.dat', unpack=True) # experiment B
for ii in range(0, nfish):
    xcol = ii*dwidth + 2
    ycol = ii*dwidth + 3
    fpsx[ii,:] = data[xcol,:]
    fpsy[ii,:] = data[ycol,:]
# calculate individual displacement over one hour time window
dx[:,1:end] = fpsx[:,1:end] - fpsx[:,0:end-1]
dy[:,1:end] = fpsy[:,1:end] - fpsy[:,0:end-1]
dx = dx - 500.0*(dx > 200)
dx = dx + 500.0*(dx < -200)
dy = dy - 500.0*(dy > 200)
dy = dy + 500.0*(dy < -200)
displacement[:,:] = 0.0
for tt in range(start+window, end): # loop of values to calculate
    disp_x[:] = 0.0
    disp_y[:] = 0.0
    for uu in range(tt-window, tt): # loop of steps in calc (10)
        disp_x[:] = disp_x[:] + dx[:,uu]
        disp_y[:] = disp_y[:] + dy[:,uu]
    displacement[:,tt] = displacement[:,tt] + disp_x**2 + disp_y**2 # displacement at end of one hour
    displacement[:,tt] = 0.1*displacement[:,tt] # average displace of one step, for each individual
# calculate summary stats
for ii in range(0, nphases):
    msd_std_b[ii,0] = np.std(displacement[:,phase_start[ii]*240:phase_end[ii]*240])
    msd_value_b[ii,0] = np.mean(displacement[:,phase_start[ii]*240:phase_end[ii]*240])
rects_b = ax[0].bar(np.arange(5)+width, msd_value_b[:,0], width, color='red', yerr=msd_std_b[:,0], ecolor='red')


# load experiment data
data = np.loadtxt('../102/fish_position.dat', unpack=True) # experiment C
for ii in range(0, nfish):
    xcol = ii*dwidth + 2
    ycol = ii*dwidth + 3
    fpsx[ii,:] = data[xcol,:]
    fpsy[ii,:] = data[ycol,:]
# calculate individual displacement over one hour time window
dx[:,1:end] = fpsx[:,1:end] - fpsx[:,0:end-1]
dy[:,1:end] = fpsy[:,1:end] - fpsy[:,0:end-1]
dx = dx - 500.0*(dx > 200)
dx = dx + 500.0*(dx < -200)
dy = dy - 500.0*(dy > 200)
dy = dy + 500.0*(dy < -200)
displacement[:,:] = 0.0
for tt in range(start+window, end): # loop of values to calculate
    disp_x[:] = 0.0
    disp_y[:] = 0.0
    for uu in range(tt-window, tt): # loop of steps in calc (10)
        disp_x[:] = disp_x[:] + dx[:,uu]
        disp_y[:] = disp_y[:] + dy[:,uu]
    displacement[:,tt] = displacement[:,tt] + disp_x**2 + disp_y**2 # displacement at end of one hour
    displacement[:,tt] = 0.1*displacement[:,tt] # average displace of one step, for each individual
# calculate summary stats
for ii in range(0, nphases):
    msd_std_c[ii,0] = np.std(displacement[0:99,phase_start[ii]*240:phase_end[ii]*240])
    msd_value_c[ii,0] = np.mean(displacement[0:99,phase_start[ii]*240:phase_end[ii]*240])
    msd_std_c[ii,1] = np.std(displacement[100:199,phase_start[ii]*240:phase_end[ii]*240])
    msd_value_c[ii,1] = np.mean(displacement[100:199,phase_start[ii]*240:phase_end[ii]*240])
rects_c = ax[0].bar(np.arange(5)+2*width, msd_value_c[:,0], width/2, color='green', yerr=msd_std_c[:,0], ecolor='green')
rects_c = ax[0].bar(np.arange(5)+2.5*width, msd_value_c[:,1], width/2, color='green', yerr=msd_std_c[:,1], ecolor='green')



# load experiment data
data = np.loadtxt('../103/fish_position.dat', unpack=True) # experiment D
for ii in range(0, nfish):
    xcol = ii*dwidth + 2
    ycol = ii*dwidth + 3
    fpsx[ii,:] = data[xcol,:]
    fpsy[ii,:] = data[ycol,:]
# calculate individual displacement over one hour time window
dx[:,1:end] = fpsx[:,1:end] - fpsx[:,0:end-1]
dy[:,1:end] = fpsy[:,1:end] - fpsy[:,0:end-1]
dx = dx - 500.0*(dx > 200)
dx = dx + 500.0*(dx < -200)
dy = dy - 500.0*(dy > 200)
dy = dy + 500.0*(dy < -200)
displacement[:,:] = 0.0
for tt in range(start+window, end): # loop of values to calculate
    disp_x[:] = 0.0
    disp_y[:] = 0.0
    for uu in range(tt-window, tt): # loop of steps in calc (10)
        disp_x[:] = disp_x[:] + dx[:,uu]
        disp_y[:] = disp_y[:] + dy[:,uu]
    displacement[:,tt] = displacement[:,tt] + disp_x**2 + disp_y**2 # displacement at end of one hour
    displacement[:,tt] = 0.1*displacement[:,tt] # average displace of one step, for each individual
# calculate summary stats
for ii in range(0, nphases):
    msd_std_d[ii,0] = np.std(displacement[0:100,phase_start[ii]*240:phase_end[ii]*240])
    msd_value_d[ii,0] = np.mean(displacement[0:100,phase_start[ii]*240:phase_end[ii]*240])
    msd_std_d[ii,1] = np.std(displacement[100:199,phase_start[ii]*240:phase_end[ii]*240])
    msd_value_d[ii,1] = np.mean(displacement[100:199,phase_start[ii]*240:phase_end[ii]*240])
rects_d = ax[0].bar(np.arange(5)+3*width, msd_value_d[:,0], width/2, color='blue', yerr=msd_std_d[:,0], ecolor='blue')
rects_d = ax[0].bar(np.arange(5)+3.5*width, msd_value_d[:,1], width/2, color='blue', yerr=msd_std_d[:,1], ecolor='blue')

# plot path and end markers
ax[0].set_xlabel(r'Trajectory class')
ax[0].set_xlim()
#ax[0].xaxis.set_major_locator(MultipleLocator())
ax[0].xaxis.grid(False)

ax[0].set_ylabel(r'Mean square displacement')
ax[0].set_ylim(0,)
#ax[0].yaxis.set_major_locator(MultipleLocator())
ax[0].yaxis.grid(False)


plt.savefig('./fish_msd_screen.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
