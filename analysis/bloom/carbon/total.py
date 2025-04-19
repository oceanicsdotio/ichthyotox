#!/usr/bin/python
"""Combined Protein and Carbohydrate Analysis"""
from matplotlib.pyplot import figure, savefig
from matplotlib.ticker import MultipleLocator
from pylab import FormatStrFormatter
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi
from analysis.bloom.carbon import State

nplots = 1

hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = 3.0
default_alpha = 1.0
start = 0


def single_plot():

    # figure and subplots
    fig = figure(facecolor=bg_color, figsize=(marginWidth, fheight))  # Change this
    fig.subplots_adjust(
        top=1.0 - vpadding,
        bottom=vpadding + 0.1,
        left=hpadding,
        right=1.0 - (hpadding / 2.0),
        hspace=0.20,
    )
    ax = fig.add_subplot(1, 1, 1)

    ax.patch.set_facecolor(bg_color)
    ax.spines["top"].set_color(overlay_color)
    ax.spines["bottom"].set_color(overlay_color)
    ax.spines["left"].set_color(overlay_color)
    ax.spines["right"].set_color(overlay_color)
    ax.xaxis.label.set_color(label_color)
    ax.yaxis.label.set_color(label_color)
    ax.tick_params(axis="x", colors=label_color)
    ax.tick_params(axis="y", colors=label_color)
    ax.xaxis.grid(False)
    return ax


if __name__ == "__main__":
    data_dir = "/Users/keeney/Documents/Projects/2020_2017_Ichthyotox/data/"
    exp = ["100", "101", "102", "103"]
    overlay = [
        overlay_color,
        [1.0, 0.0, 0.0, default_alpha],
        [0.0, 1.0, 0.0, default_alpha],
        [0.0, 0.0, 1.0, default_alpha],
    ]
    fillc = [
        [0.0, 0.0, 0.0, 0.25],
        [1.0, 0.0, 0.0, 0.25],
        [0.0, 1.0, 0.0, 0.25],
        [0.0, 0.0, 1.0, 0.25],
    ]

    volume = 500.0 * 500.0 * 5.0
    axes = single_plot()
    for ii in range(0, nplots):
        source = data_dir + exp[ii]
        state = State(source)
        carbon = state.carbon_concentration(volume)
        carbon.plot(
            axes, time=state.days, color=overlay[ii], fill=fillc[ii], with_bounds=True
        )

    axes.set_xlabel(r"Days")
    axes.xaxis.set_major_locator(MultipleLocator(5))
    axes.xaxis.set_minor_locator(MultipleLocator(1))
    axes.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    # axes.set_xlim(time[start],time[end])

    axes.set_ylabel(r"Total Carbon (g/m^3)", rotation=90)

    axes.set_ylim(0, 0.55)
    axes.yaxis.set_major_locator(MultipleLocator(0.1))
    axes.yaxis.set_minor_locator(MultipleLocator(0.02))
    axes.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    axes.yaxis.grid(False)
    axes.set_frame_on(True)

    legend = axes.legend(loc="upper left")
    legend.get_frame().set_facecolor("none")
    legend.get_frame().set_edgecolor("none")
    text = legend.get_texts()
    # text[0].set_color('black')
    # text[1].set_color('red')
    # text[2].set_color('green')
    # text[3].set_color('blue')

    savefig(
        "./figures/bloom/carbon/total.png",
        facecolor=bg_color,
        edgecolor="none",
        dpi=default_dpi,
    )
