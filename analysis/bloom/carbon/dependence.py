#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
from pylab import *
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi
from analysis.bloom.carbon import State


fontsize  = 10
linewidth = 1
uniformPadding = 0.15
hpadding = 0.1
vpadding = 0.1
marginWidth = 6.5
fheight = (marginWidth-0.5)/2
default_alpha = 1.0
lineRGBA = [1.0, 0.0, 0.0, default_alpha]
style = ['-','--',':']
div=150
show_grid=True
dwidth = 4
start = 0


# figure and subplots
fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight))
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.05, left=hpadding, right=1.0-hpadding)
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

def read_data(ax, source, line_style):
    state = State(source)
    time = time / 24.
    volume = 500.0*500.0*5.0
    pro_avg = np.mean(state.protein, axis=0)
    carb_avg = np.mean(state.carbohydrate, axis=0)
    tot_carb = float(state.count)*(pro_avg + carb_avg)
    ax.plot(time, tot_carb/volume, linestyle=line_style, color=overlay_color, aa=True)


if __name__ == "__main__":
    data_dir = "/Users/keeney/Documents/Projects/2020_2017_Ichthyotox/data/"
    exp = ['100', '101', '103']
    for ii in range(0, len(exp)):
        source = data_dir + exp[ii]
        read_data(ax[0], source, style[ii])
        
    sublabels = [[r'Days', r'Biomass (g/m^3)'], [r'Days', r'Carbon Ratio (g/g)']]
    #subtitles = [r'(A)', r'(B)']
    for ii in range(0,1):
        # axis adjustment

        ax[ii].set_xlabel(sublabels[ii][0])
        ax[ii].xaxis.set_major_locator(mticker.MultipleLocator(5))
        ax[ii].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))

            
        # ax[ii].set_xlim(time[start],time[end])
        ax[ii].xaxis.grid(False)
        
        ax[ii].set_ylabel(sublabels[ii][1], rotation=90)
        ax[ii].set_ylim(0.0,15.0)
        ax[ii].yaxis.set_major_locator(mticker.MultipleLocator(5))
        ax[ii].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        ax[ii].yaxis.grid(False)
        
        #ax[ii].set_title(subtitles[ii], fontsize=fontsize, color=label_color)
        #ax[ii].title.set_position([0.5, 1.0])
        ax[ii].set_frame_on(True)
   
    plt.savefig('./figures/bloom/carb_dependence.png', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
