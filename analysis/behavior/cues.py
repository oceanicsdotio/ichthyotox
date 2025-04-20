#!/usr/bin/python
"""Plot cue series for fish behavior triggers"""
from sys import argv
from matplotlib.pyplot import figure, savefig, LinearLocator, FormatStrFormatter, MultipleLocator
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi, fontsize
from analysis.behavior import State, Positions


if __name__ == "__main__":
    data_dir = argv[1]
    hpadding = 0.1
    vpadding = 0.05
    marginWidth = 6.5
    fheight = 7.0
    pindex = [0,99,100,199]
    lindex = 19

    # figure and subplots
    fig = figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
    fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.05, left=hpadding/2.0, right=1.0-(hpadding/2.0))
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

    def plot_experiment(experiment: str, lindex):
        state = State(experiment)
        state.plot_fraction_intoxicated(ax, lindex)
        state.plot_cues(ax, pindex, lindex)
        position = Positions(experiment)
        lindex += 5
        position.plot_fraction_suitable(ax, lindex)
        position.plot_suitability_cue(ax, pindex, lindex)

    plot_experiment(data_dir + "100/", lindex)

    ax.set_xlabel(r'Days')
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.set_xlim(0,30)
    ax.xaxis.grid(False)
    ax.set_ylabel('Behavioral Cues', rotation=90)
    ax.yaxis.set_major_locator(LinearLocator(0))
    ax.yaxis.grid(False)
    ax.set_frame_on(True)
    ax.set_ylim(0,22)

    ax.text(x=0.05*30, y=21.25, s='Control (A)', color='black', fontsize=fontsize)
    ax.text(x=0.05*30, y=15.75, s='Formation (B)', color='red', fontsize=fontsize)
    ax.text(x=0.05*30, y=10.25, s='Intensification (C)', color='green', fontsize=fontsize)
    ax.text(x=0.05*30, y=4.75, s='Decline (D)', color='blue', fontsize=fontsize)

    savefig('./figures/behavior/cues.png', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
