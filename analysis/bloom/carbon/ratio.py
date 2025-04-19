#!/usr/bin/python
from matplotlib.pyplot import figure, savefig
from matplotlib.ticker import MultipleLocator
import numpy as np
from pylab import *
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi, fontsize
from analysis.bloom.carbon import State

nplots=4
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = marginWidth + 0.5
dwidth = 4
start = 0


def multiple_plots():
    fig = figure(facecolor=bg_color, figsize=(marginWidth, fheight))
    fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.05, left=hpadding, right=1.0-hpadding/2., hspace=0.20)
    ax = []
    ax.append(fig.add_subplot(nplots,1,1))
    ax.append(fig.add_subplot(nplots,1,2))
    ax.append(fig.add_subplot(nplots,1,3))
    ax.append(fig.add_subplot(nplots,1,4))

    ax[0].patch.set_facecolor(bg_color)
    ax[0].spines['top'].set_color(overlay_color)
    ax[0].spines['bottom'].set_color(overlay_color)
    ax[0].spines['left'].set_color(overlay_color)
    ax[0].spines['right'].set_color(overlay_color)
    ax[0].xaxis.label.set_color(label_color)
    ax[0].yaxis.label.set_color(label_color)
    ax[0].tick_params(axis='x', colors=label_color)
    ax[0].tick_params(axis='y', colors=label_color)



def plot_experiment_data_subplot(ax, time, sum_avg, sum_std, overlay_color):
    ax.plot(time, (sum_avg+sum_std), linestyle='-', linewidth=1, color=overlay_color, aa=True)
    ax.plot(time, sum_avg, linestyle='-', linewidth=2,  color=overlay_color, aa=True)
    ax.plot(time, (sum_avg-sum_std), linestyle='-', linewidth=1,  color=overlay_color, aa=True)
    return ax




def format_single_plot():
# figure and subplots
    fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
    fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.1, left=hpadding, right=1.0-hpadding/2.0, hspace=0.20)
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

    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.set_xlabel('Days')
    # ax.set_xlim(time[start],time[end])
    ax.xaxis.grid(False)
    ax.set_ylabel('Mean Carbon Ratio (g/g) ', rotation=90)
    ax.yaxis.set_major_locator(MultipleLocator(1.0))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.grid(False)
    ax.set_frame_on(True)
    ax.set_ylim(0.0,4.0)


    legend = ax.legend(loc='upper right')
    legend.get_frame().set_facecolor('none')
    legend.get_frame().set_edgecolor('none')
    text = legend.get_texts()
    text[0].set_color('white')
    text[1].set_color('red')
    text[2].set_color('green')
    text[3].set_color('blue')
    return ax

if __name__ == "__main__":

    data_dir = "/Users/keeney/Documents/Projects/2020_2017_Ichthyotox/data/"
    iplot = 0
    for experiment in ['100', '101', '102', '103']:
        source = data_dir + experiment
        print(f"Reading data from {source}")
        state = State(source)
        ratio = state.carbon_ratio()
        ratio.plot(ax[iplot], state.days, ratio_mean, ratio_std_dev, overlay_color)
        iplot += 1

    sublabels = [
        [r'Days', r'Carbon Ratio (g/g)'],
        [r'Days', r'Carbon Ratio (g/g)'],
        [r'Days', r'Carbon Ratio (g/g)'],
        [r'Days', r'Carbon Ratio (g/g)']
    ]
    subtitles = [r'(A)', r'(B)',r'(C)', r'(D)']
    ylim = [4.0, 1.5, 0.5, 0.5]
    for ii in range(0, nplots):
        # axis adjustment
        if (ii == nplots-1):
            ax[ii].set_xlabel(sublabels[ii][0])
            ax[ii].xaxis.set_major_locator(mticker.MultipleLocator(5))
            ax[ii].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        else:
            ax[ii].set_xlabel(' ')
            ax[ii].xaxis.set_major_locator(LinearLocator(0))

            
        ax[ii].set_xlim(time[0],time[-1])
        ax[ii].xaxis.grid(False); 
        
        ax[ii].set_ylabel(sublabels[ii][1], rotation=90)

        ax[ii].yaxis.set_major_locator(mticker.MultipleLocator(0.5))
        ax[ii].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax[ii].yaxis.grid(False)
        
        ax[ii].set_title(subtitles[ii], fontsize=fontsize, color=label_color)
        ax[ii].title.set_position([1.0, 1.0])
        ax[ii].set_frame_on(True)
        
        ax[ii].set_ylim(0.0,ylim[ii])

    plt.savefig('./figures/bloom/carbon_ratio.png', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
