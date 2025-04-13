#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc

fontsize  = 14
linewidth = 1
uniformPadding = 0.1
marginWidth = 7.0
for_screen=False


default_alpha = 0.25
default_dpi = 150
lineRGBA = [1.0, 0.0, 0.0, default_alpha]
style = ['-',':']


rc('text', usetex=False)
rc('font', **{'family':'sans-serif', 'sans-serif':['Avant Garde']})
#rc('font', weight='bold')
#rc('mathtext', default='sf')
rc('lines', markeredgewidth=1)
rc('lines', linewidth=linewidth)
rc('axes', labelsize=fontsize)
rc('axes', linewidth=(linewidth+1)//2)
rc('xtick', labelsize=2*fontsize/3)
rc('ytick', labelsize=2*fontsize/3)
rc('legend', fontsize=2*fontsize/3)
rc('xtick.major', pad=5)
rc('ytick.major', pad=5)



if for_screen:
    fig = plt.figure(facecolor='black', figsize=(marginWidth, marginWidth)) #Change this
    fig.subplots_adjust(top=1.0-uniformPadding, bottom=uniformPadding, left=uniformPadding, right=1.0-uniformPadding)

    bg_color = [0.0,0.0,0.0,1.0]
    overlay_color = [1.0,1.0,1.0,1.0]
    label_color = [0.5,0.5,0.5,1.0]
    
else:
    fig = plt.figure(figsize=(marginWidth, marginWidth)) #Change this
    fig.subplots_adjust(top=1.0-uniformPadding, bottom=uniformPadding, left=uniformPadding, right=1.0-uniformPadding)
    
    bg_color = [1.0,1.0,1.0,1.0]
    overlay_color = [0.0,0.0,0.0,1.0]
    label_color = [0.0,0.0,0.0,1.0]
   
   
# set color scheme
ax = fig.add_subplot(1,1,1)
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
    for l in ax.get_xticklines() + ax.get_yticklines():
        l.set_markersize(major)
    for tick in ax.xaxis.get_minor_ticks() + ax.yaxis.get_minor_ticks():
        tick.tick1line.set_markersize(minor); tick.tick2line.set_markersize(minor)
    ax.xaxis.LABELPAD      = 10.
    ax.xaxis.OFFSETTEXTPAD = 10.
