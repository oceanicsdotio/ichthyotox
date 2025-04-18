#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from pylab import *
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi, fontsize

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

def read_experiment_data(experiment):
    ini = open(f'{experiment}/cyanobacteria_ini.dat', 'r')
    ncolony = int(str.strip(ini.readline()))
    time = np.loadtxt(f'{experiment}/cyanobacteria_state.dat', usecols=[0], unpack=True) / 24.0
    end = len(time) - 1
    data = np.loadtxt(f'{experiment}/cyanobacteria_state.dat', unpack=True)
    carbohydrate = zeros((ncolony,end-start+1))
    protein = zeros((ncolony,end-start+1))
    for ii in range(0, ncolony):
        ccol = ii*dwidth + 2
        pcol = ii*dwidth + 3
        carbohydrate[ii,:] = data[ccol,:]
        protein[ii,:] = data[pcol,:]
    return time, carbohydrate, protein

def summarize_experiment_data(carbohydrate, protein):
    ratio = carbohydrate / protein
    sum_avg = np.mean(ratio, axis=0)
    sum_std = np.std(ratio, axis=0)
    print ("Mean Ratio", sum_avg[-1], "+/-", sum_std[-1])
    print ("last day max", np.max(sum_avg[-241:-1]))
    return sum_avg, sum_std

def plot_experiment_data(ax, time, sum_avg, sum_std, overlay_color):
    ax.plot(time, (sum_avg+sum_std), linestyle='-', linewidth=1, color=overlay_color, aa=True)
    ax.plot(time, sum_avg, linestyle='-', linewidth=2,  color=overlay_color, aa=True)
    ax.plot(time, (sum_avg-sum_std), linestyle='-', linewidth=1,  color=overlay_color, aa=True)
    return ax



if __name__ == "__main__":

    data_dir = "/Users/keeney/Documents/Projects/2020_2017_Ichthyotox/data/"
    iplot = 0
    for experiment in ['100', '101', '102', '103']:
        source = data_dir + experiment
        print(f"Reading data from {source}")
        time, carbohydrate, protein = read_experiment_data(source)
        ratio_mean, ratio_std_dev = summarize_experiment_data(carbohydrate, protein)
        plot_experiment_data(ax[iplot], time, ratio_mean, ratio_std_dev, overlay_color)
        iplot += 1

    sublabels = [
        [r'Days', r'Carbon Ratio (g/g)'],
        [r'Days', r'Carbon Ratio (g/g)'],
        [r'Days', r'Carbon Ratio (g/g)'],
        [r'Days', r'Carbon Ratio (g/g)']
    ]
    subtitles = [r'(A)', r'(B)',r'(C)', r'(D)']
    ylim = [4.0, 1.5, 0.5, 0.5]
    for ii in range(0,nplots):
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
