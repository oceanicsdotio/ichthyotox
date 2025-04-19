#!/usr/bin/python
"""Generate figures for bloom protein analysis."""
from matplotlib.pyplot import figure, savefig
from matplotlib.ticker import MultipleLocator
from pylab import FormatStrFormatter
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi
from analysis.bloom.carbon import State


def single_plot():
    """Multiple series on a single plot"""
    hpadding = 0.1
    vpadding = 0.05
    marginWidth = 6.5
    fheight = 3.0

    fig = figure(facecolor=bg_color, figsize=(marginWidth, fheight))
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
    ax.set_xlabel(r"Days")
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    # ax.set_xlim(time[0],time[end])
    ax.xaxis.grid(False)
    ax.set_ylabel(r"Mean Protein (g/m^3)", rotation=90)
    ax.set_ylim(0, 0.25)
    ax.yaxis.set_major_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.yaxis.grid(False)
    ax.set_frame_on(True)
    return ax


if __name__ == "__main__":

    VOLUME = 500.0 * 500.0 * 5.0
    data_dir = "/Users/keeney/Documents/Projects/2020_2017_Ichthyotox/data/"
    axes = single_plot()
    source = f"{data_dir}/100"
    state = State(source)
    protein = state.protein_concentration(VOLUME)
    protein.plot(axes, time=state.days, color="black", fill="grey", with_bounds=True)

    legend = axes.legend(loc="center left")
    legend.get_frame().set_facecolor("none")
    legend.get_frame().set_edgecolor("none")
    text = legend.get_texts()
    text[0].set_color("black")
    # text[1].set_color("red")
    # text[2].set_color("green")
    # text[3].set_color("blue")
    savefig(
        "./figures/bloom/carbon/protein.png",
        facecolor=bg_color,
        edgecolor="none",
        dpi=default_dpi,
    )
