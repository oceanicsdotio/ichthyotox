#!/usr/bin/python
"""Default settings for matplotlib figures."""
from matplotlib.pyplot import Axes
from matplotlib import rc

IMAGE_FORMAT = "png"
SOURCE = "./data/"
fontsize = 12
linewidth = 1
uniform_padding = 0.1
margin_width = 7.0
for_screen = False
default_alpha = 0.25
default_dpi = 150
lineRGBA = [1.0, 0.0, 0.0, default_alpha]
style = ["-", ":"]


rc("text", usetex=False)
rc("font", **{"family": "sans-serif", "sans-serif": ["Arial"]})
rc("font", weight="normal")
rc("mathtext", default="sf")
rc("lines", markeredgewidth=1)
rc("lines", linewidth=linewidth)
rc("axes", labelsize=fontsize)
rc("axes", linewidth=(linewidth + 1) // 2)
rc("xtick", labelsize=2 * fontsize / 3)
rc("ytick", labelsize=2 * fontsize / 3)
rc("legend", fontsize=2 * fontsize / 3)
rc("xtick.major", pad=5)
rc("ytick.major", pad=5)


if for_screen:
    bg_color = [0.0, 0.0, 0.0, 1.0]
    overlay_color = [1.0, 1.0, 1.0, 1.0]
    label_color = [0.5, 0.5, 0.5, 1.0]

else:
    bg_color = [1.0, 1.0, 1.0, 1.0]
    overlay_color = [0.0, 0.0, 0.0, 1.0]
    label_color = [0.0, 0.0, 0.0, 1.0]

overlay = [
    overlay_color,
    [1.0, 0.0, 0.0, 1.0],
    [0.0, 1.0, 0.0, 1.0],
    [0.0, 0.0, 1.0, 1.0],
]
fill_color = [
    [0.0, 0.0, 0.0, 0.25],
    [1.0, 0.0, 0.0, 0.25],
    [0.0, 1.0, 0.0, 0.25],
    [0.0, 0.0, 1.0, 0.25],
]

def set_tick_sizes(ax: Axes, major, minor):
    """Animation tick size."""
    for l in ax.get_xticklines() + ax.get_yticklines():
        l.set_markersize(major)
    for tick in ax.xaxis.get_minor_ticks() + ax.yaxis.get_minor_ticks():
        tick.tick1line.set_markersize(minor)
        tick.tick2line.set_markersize(minor)
    ax.xaxis.LABELPAD = 10.0
    ax.xaxis.OFFSETTEXTPAD = 10.0
