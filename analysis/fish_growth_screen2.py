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
marginWidth = 6.5
fheight = 3.0
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
suitability = zeros((nfish,end-start+1))
xpos = zeros((nfish,end-start+1))
mass_control = zeros(end-start+1)
mass_avg_all = zeros(end-start+1)
mass_avg_surface = zeros(end-start+1)
mass_avg_bottom = zeros(end-start+1)
time = time / 24.


# figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.1, left=hpadding, right=1.0-hpadding)
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

#######################################
print "Reading State" # experiment A
data = np.loadtxt('../100/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol, :]
print "Calculating Mean Mass"
mass_control =  np.mean(mass[:,:], axis=0)
ax[0].plot(time, mass_control, linestyle='-', linewidth=1, color='white', aa=True, label='Control (A)')
#######################################
print "Reading State" # experiment B
data = np.loadtxt('../101/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol, :]
print "Calculating Mean Mass"
mass_avg =  np.mean(mass[:,:], axis=0)
ax[0].plot(time, mass_avg, linestyle='-', linewidth=1, color='red', aa=True, label='Formation (B)')
#######################################
print "Reading State" 
data = np.loadtxt('../300/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol, :]
print "Calculating Mean Mass"
mass_avg =  np.mean(mass[:,:], axis=0)
ax[0].plot(time, mass_avg, linestyle='-', linewidth=1, color='orange', aa=True, label='Low Suit.')
#######################################
print "Reading State" 
data = np.loadtxt('../301/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol, :]
print "Calculating Mean Mass"
mass_avg =  np.mean(mass[:,:], axis=0)
ax[0].plot(time, mass_avg, linestyle='-', linewidth=1, color='yellow', aa=True, label='Low Suit./Low Dep.')
#######################################
print "Reading State" 
data = np.loadtxt('../302/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol, :]
print "Calculating Mean Mass"
mass_avg =  np.mean(mass[:,:], axis=0)
ax[0].plot(time, mass_avg, linestyle='-', linewidth=1, color='purple', aa=True, label='Low Suit./Extra Impaired')
#######################################
print "Reading State" 
data = np.loadtxt('../303/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol, :]
print "Calculating Mean Mass"
mass_avg =  np.mean(mass[:,:], axis=0)
ax[0].plot(time, mass_avg, linestyle='-', linewidth=1, color='green', aa=True, label='Low Suit./Low Dep/Extra Impaired')
#######################################
print "Reading State" 
data = np.loadtxt('../304/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol, :]
print "Calculating Mean Mass"
mass_avg =  np.mean(mass[:,:], axis=0)
ax[0].plot(time, mass_avg, linestyle='-', linewidth=1, color='blue', aa=True, label='Extra Impaired')
#######################################
print "Reading State" 
data = np.loadtxt('../305/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol, :]
print "Calculating Mean Mass"
mass_avg =  np.mean(mass[:,:], axis=0)
ax[0].plot(time, mass_avg, linestyle='-', linewidth=1, color='grey', aa=True, label='Low Dep.')
#######################################
print "Reading State" 
data = np.loadtxt('../306/fish_state.dat', unpack=True)
for ii in range(0, nfish):
    masscol = ii*dwidth + 2
    mass[ii,:] = data[masscol, :]
print "Calculating Mean Mass"
mass_avg =  np.mean(mass[:,:], axis=0)
ax[0].plot(time, mass_avg, linestyle='--', linewidth=1, color='white', aa=True, label='Low Dep./Extra Impaired')
#######################################


ax[0].set_xlabel(r'Days')
ax[0].xaxis.set_major_locator(MultipleLocator(5))
ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax[0].set_xlim(0,30)
ax[0].xaxis.grid(False); 

ax[0].set_ylabel('Mass (g)', rotation=90)

ax[0].yaxis.set_major_locator(MultipleLocator(1))
ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax[0].yaxis.grid(False);
ax[0].set_frame_on(True)
ax[0].set_ylim(11,19)

legend = ax[0].legend(loc='upper left')
legend.get_frame().set_facecolor('none')
legend.get_frame().set_edgecolor('none')
text = legend.get_texts()
text[0].set_color('white')
text[1].set_color('red')
text[2].set_color('orange')
text[3].set_color('yellow')
text[4].set_color('purple')
text[5].set_color('green')
text[6].set_color('blue')
text[7].set_color('grey')
text[8].set_color('white')


#plt.show()
plt.savefig('./fish_growth_screen2.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
