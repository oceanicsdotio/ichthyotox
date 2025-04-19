#!/usr/bin/python
from matplotlib.pyplot import figure, savefig
import matplotlib.ticker as mticker
import numpy as np
from pylab import FormatStrFormatter, Axes
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi
from analysis.bloom.carbon import State

nplots=1
fontsize  = 10
linewidth = 1
uniformPadding = 0.15
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = 3.0
default_alpha = 1.0
start = 0
dwidth = 4

volume = 500.0*500.0*5.0


# figure and subplots
fig = figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.1, left=hpadding, right=1.0-(hpadding/2.), hspace=0.20)
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


def read_experiment_data(experiment, color, fill_bounds):
    state = State(experiment)
    
    volume = 500.0*500.0*5.0
    
    print ("Calculating Mean")
    totalc = state.carbon_concentration(volume)
    sum_avg = np.mean(totalc, axis=0)
    sum_std = np.std(totalc, axis=0)
    ax[0].plot(time, sum_avg, linestyle='-', linewidth=1, color=color, aa=True, label='Control (A)')
    ax[0].fill_between(time, (sum_avg+sum_std), (sum_avg-sum_std), facecolor=fill_bounds, edgecolor='none', zorder=3)
    print ("Total Protein: ", sum(state.protein[:,-1]))
    print ("Mean Biomass: ", sum_avg[-1], "+/-", sum_std[-1])



if __name__ == "__main__":
    data_dir = "/Users/keeney/Documents/Projects/2020_2017_Ichthyotox/data/"
    exp = ["100", "101", "102", "103"]
    overlay = [overlay_color, [1.0,0.0,0.0,default_alpha], [0.0,1.0,0.0,default_alpha], [0.0,0.0,1.0,default_alpha]]
    fillc = [[0.0,0.0,0.0,0.25], [1.0,0.0,0.0,0.25], [0.0,1.0,0.0,0.25], [0.0,0.0,1.0,0.25]]


    # figure and subplots
    fig = figure(facecolor=bg_color, figsize=(marginWidth, fheight))
    fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.1, left=hpadding, right=1.0-(hpadding/2.), hspace=0.20)
    ax: list[Axes] = []
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
    for ii in range(0, nplots):
        source = data_dir + exp[ii]
        state = State(source)
        carbon = state.carbon_concentration(volume)
        carbon.plot(ax[0], with_bounds=True)

    ax[0].set_xlabel(r'Days')
    ax[0].xaxis.set_major_locator(mticker.MultipleLocator(5))
    ax[0].xaxis.set_minor_locator(mticker.MultipleLocator(1))
    ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    # ax[0].set_xlim(time[start],time[end])
    ax[0].xaxis.grid(False)
    ax[0].set_ylabel(r'Total Carbon (g/m^3)', rotation=90)

    ax[0].set_ylim(0,0.55)
    ax[0].yaxis.set_major_locator(mticker.MultipleLocator(0.1))
    ax[0].yaxis.set_minor_locator(mticker.MultipleLocator(0.02))
    ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax[0].yaxis.grid(False)
    ax[0].set_frame_on(True)

    legend = ax[0].legend(loc='upper left')
    legend.get_frame().set_facecolor('none')
    legend.get_frame().set_edgecolor('none')
    text = legend.get_texts()
    # text[0].set_color('black')
    # text[1].set_color('red')
    # text[2].set_color('green')
    # text[3].set_color('blue')

    savefig('./figures/bloom/carbon.png', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
