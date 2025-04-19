#!/usr/bin/python
"""Draw particle trajectories and suitability surface"""
from sys import argv
from matplotlib.pyplot import figure, savefig, title, MultipleLocator
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi
from analysis.behavior import Domain, Positions


fontsize  = 10
uniformPadding = 0.15
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = marginWidth + 0.5

if __name__ == "__main__":
    data_dir = argv[1]
    # load data
    fish = Positions(data_dir)

    # Figure and subplots
    fig = figure(facecolor=bg_color, figsize=(marginWidth, fheight))
    fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.05, left=hpadding, right=1.0-hpadding/2., hspace=0.20)
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
    # Calculate Domain Data
    mesh = Domain(data_dir)
    mesh.draw()

    fish.draw_final()
    fish.draw_trajectory(ax)

    ax.set_xlabel(r'X (m)')
    ax.set_xlim(*mesh.x_range())
    ax.xaxis.set_major_locator(MultipleLocator(125))
    ax.xaxis.grid(False)

    ax.set_ylabel(r'Y (m)')
    ax.set_ylim(*mesh.y_range())
    ax.yaxis.set_major_locator(MultipleLocator(125))
    ax.yaxis.grid(False)


    title(r'(A)', fontsize=fontsize, color=label_color)
    ttl = ax.title
    ttl.set_position([1.0, 1.0])
    savefig('./figures/behavior/paths.png', edgecolor='none', dpi=default_dpi)

