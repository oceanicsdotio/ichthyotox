#!/usr/bin/python
from sys import argv
from matplotlib.pyplot import figure, savefig, scatter, LinearLocator, get_cmap
from numpy import loadtxt, zeros
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi
from analysis.behavior import Domain

nplots=1
fontsize  = 10
linewidth = 1
uniformPadding = 0.15
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = marginWidth + 0.5


class InitialPosition:
    """Fish particle position data"""
    def __init__(self, experiment: str, fmt="dat"):
        ini_file = f'{experiment}/fish_ini.{fmt}'
        with open(ini_file, 'r', encoding="utf8") as ini:
            self.count = int(str.strip(ini.readline()))
        self.x = zeros(self.count)
        self.y = zeros(self.count)
        self.z = zeros(self.count)
        self.x[:], self.y[:], self.z[:] = loadtxt(ini_file, skiprows=1, usecols=(1,2,3), unpack=True)

    def draw(self):
        surface = self.z > -2.5
        scatter(self.x[surface], self.y[surface], s=40, color='black', zorder=10, edgecolors='black')
        scatter(self.x[~surface], self.y[~surface], s=40, color='none', zorder=10, edgecolors='black')



def single_figure(experiment: str):
    # Figure and subplots
    fig = figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
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

    mesh = Domain(experiment)
    mesh.draw()
    
    positions = InitialPosition(experiment)
    positions.draw()


    ax.set_xlabel(r'X (m)')
    ax.set_xlim(*mesh.x_range())
    ax.xaxis.set_major_locator(LinearLocator(2))
    ax.xaxis.grid(False)

    ax.set_ylabel(r'Y (m)')
    ax.set_ylim(*mesh.y_range())
    ax.yaxis.set_major_locator(LinearLocator(2))
    ax.yaxis.grid(False)

    savefig('./figures/behavior/initial_position.png', edgecolor='none', dpi=default_dpi)

if __name__ == "__main__":
    experiment = argv[1]
    single_figure(experiment)
