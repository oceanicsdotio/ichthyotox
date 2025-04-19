#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
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

fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight))
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



def summarize_experiment_data(carbohydrate, protein):
    ratio = carbohydrate / protein
    sum_avg = np.mean(ratio, axis=0)
    sum_std = np.std(ratio, axis=0)
    print ("Mean Ratio", sum_avg[-1], "+/-", sum_std[-1])
    print ("last day max", np.max(sum_avg[-241:-1]))
    return sum_avg, sum_std

def plot_experiment_data_subplot(ax, time, sum_avg, sum_std, overlay_color):
    ax.plot(time, (sum_avg+sum_std), linestyle='-', linewidth=1, color=overlay_color, aa=True)
    ax.plot(time, sum_avg, linestyle='-', linewidth=2,  color=overlay_color, aa=True)
    ax.plot(time, (sum_avg-sum_std), linestyle='-', linewidth=1,  color=overlay_color, aa=True)
    return ax


def read_data_and_plot_screen(experiment, color, label, zorder):
    # load data
    ini = open(f'{experiment}/cyanobacteria_ini.dat', 'r')
    ncolony = int(str.strip(ini.readline()))
    time = np.loadtxt(f'{experiment}/cyanobacteria_state.dat', usecols=[0], unpack=True)
    start = 0
    end = len(time) - 1
    dwidth = 4
    data = zeros((ncolony*dwidth+1, end-start+1))
    carbohydrate = zeros((ncolony,end-start+1))
    protein = zeros((ncolony,end-start+1))
    ratioc = zeros((ncolony,end-start+1))
    ratio_avg = zeros(end-start+1)
    time = time / 24.
    # volume = 500.0*500.0*5.0
    data = np.loadtxt(f'{experiment}/cyanobacteria_state.dat', unpack=True)
    for ii in range(0, ncolony):
        ccol = ii*dwidth + 2
        pcol = ii*dwidth + 3
        carbohydrate[ii,:] = data[ccol,:]
        protein[ii,:] = data[pcol,:]
    ratioc = carbohydrate / protein
    ratio_avg = np.mean(ratioc, axis=0)
    # ratio_std = np.std(ratioc, axis=0)
    ax[0].plot(time, ratio_avg, linestyle='-', linewidth=1, color=color, aa=True, label=label, zorder=zorder)
    #ax[0].fill_between(time, (ratio_avg-ratio_std), (ratio_avg+ratio_std), facecolor=[0.0,0.0,0.0,0.25], edgecolor='none', zorder=1)

def format_single_plot():
# figure and subplots
    fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
    fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.1, left=hpadding, right=1.0-hpadding/2.0, hspace=0.20)
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

    ax[0].xaxis.set_major_locator(MultipleLocator(5))
    ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax[0].set_xlabel('Days')
    # ax[0].set_xlim(time[start],time[end])
    ax[0].xaxis.grid(False)
    ax[0].set_ylabel('Mean Carbon Ratio (g/g) ', rotation=90)
    ax[0].yaxis.set_major_locator(MultipleLocator(1.0))
    ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax[0].yaxis.grid(False)
    ax[0].set_frame_on(True)
    ax[0].set_ylim(0.0,4.0)


    legend = ax[0].legend(loc='upper right')
    legend.get_frame().set_facecolor('none')
    legend.get_frame().set_edgecolor('none')
    text = legend.get_texts()
    text[0].set_color('white')
    text[1].set_color('red')
    text[2].set_color('green')
    text[3].set_color('blue')

if __name__ == "__main__":

    data_dir = "/Users/keeney/Documents/Projects/2020_2017_Ichthyotox/data/"
    iplot = 0
    for experiment in ['100', '101', '102', '103']:
        source = data_dir + experiment
        print(f"Reading data from {source}")
        state = State(source)
        ratio_mean, ratio_std_dev = summarize_experiment_data(state.carbohydrate, state.protein)
        plot_experiment_data_subplot(ax[iplot], state.days, ratio_mean, ratio_std_dev, overlay_color)
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
