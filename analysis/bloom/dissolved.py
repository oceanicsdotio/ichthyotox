#!/usr/bin/python
from matplotlib.pyplot import figure, savefig, MultipleLocator, LinearLocator, FormatStrFormatter
from numpy import loadtxt, mean, zeros, linspace
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi

nplots = 2
fontsize = 10
linewidth = 1
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = 2.0 * (marginWidth + 0.5) / 4.0
default_alpha = 1.0
lineRGBA = [1.0, 0.0, 0.0, default_alpha]
style = ["-", ":"]
div = 150


# figure and subplots
fig = figure(facecolor=bg_color, figsize=(marginWidth, fheight))
fig.subplots_adjust(
    top=1.0 - vpadding - 0.05,
    bottom=vpadding + 0.1,
    left=hpadding,
    right=1.0 - hpadding / 2.0,
    hspace=0.2,
)
ax = []
ax.append(fig.add_subplot(2, 1, 1))
ax.append(fig.add_subplot(2, 1, 2))

ax[0].patch.set_facecolor(bg_color)
ax[0].spines["top"].set_color(overlay_color)
ax[0].spines["bottom"].set_color(overlay_color)
ax[0].spines["left"].set_color(overlay_color)
ax[0].spines["right"].set_color(overlay_color)
ax[0].xaxis.label.set_color(label_color)
ax[0].yaxis.label.set_color(label_color)
ax[0].tick_params(axis="x", colors=label_color)
ax[0].tick_params(axis="y", colors=label_color)

class Profile:
    """Vertical profile of dissolved toxin concentration."""
    def __init__(self, experiment: str):
        # load data
        nlayers = 26
        time = loadtxt(f"{experiment}/dissolved_toxin.dat", usecols=[0], unpack=True)
        start = 0
        end = len(time) - 1
        shape = (nlayers, len(time))

        dissolved = zeros(shape)
        avg_dissolved = zeros(len(time))
        time = time / 24.0
        volume = 500.0 * 500.0 * 5.0
        self.depth = linspace(0.0, -5.0, nlayers)

        data = loadtxt(f"{experiment}/dissolved_toxin.dat", unpack=True)
        dissolved[:, :] = data[1:, :]
        dissolved = dissolved / (500.0 * 500.0)

        self.mean = mean(dissolved, axis=0)

    def plot(self):
        ax[0].plot(
            (dissolved[:, 240 * 3]),
            self.depth,
            linestyle="-",
            linewidth=1,
            color=overlay_color,
            aa=True,
        )
        ax[0].plot(
            (dissolved[:, -1]), self.depth, linestyle="-", linewidth=2, color=overlay_color, aa=True
        )
        ax[0].plot(
            ((avg_dissolved[-1]), (avg_dissolved[-1])),
            (-5.0, 0.0),
            linestyle="--",
            linewidth=1,
            color=overlay_color,
            aa=True,
        )

if __name__ == "__main__":

    #######################################
    print("Reading File 013")  # experiment D
    data = loadtxt("../013/dissolved_toxin.dat", unpack=True)
    for ii in range(0, nlayers):
        dissolved[ii, :] = data[ii + 1, :]
    dissolved = dissolved / (500.0 * 500.0)
    avg_dissolved = mean(dissolved, axis=0)
    ax[1].plot(
        (dissolved[:, 240 * 3]),
        depth,
        linestyle="-",
        linewidth=1,
        color=overlay_color,
        aa=True,
    )
    ax[1].plot(
        (dissolved[:, end]), depth, linestyle="-", linewidth=2, color=overlay_color, aa=True
    )
    ax[1].plot(
        ((avg_dissolved[end]), (avg_dissolved[end])),
        (-5.0, 0.0),
        linestyle="--",
        linewidth=1,
        color=overlay_color,
        aa=True,
    )

    print("Mean Dissolved D: ", avg_dissolved[end] * 500.0 * 500.0)
    #######################################

    sublabels = [
        [r"Dissolved Toxin (g/m^3)", r"Depth (m)"],
        [r"Dissolved Toxin (g/m^3)", r"Depth (m)"],
    ]
    subtitles = [r"(A)", r"(B)"]
    for ii in range(0, nplots):
        # axis adjustment
        if ii == nplots - 1:
            ax[ii].set_xlabel(sublabels[ii][0])
            ax[ii].xaxis.set_major_locator(MultipleLocator(0.05))
            ax[ii].xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        else:
            ax[ii].set_xlabel(" ")
            ax[ii].xaxis.set_major_locator(LinearLocator(0))

        ax[ii].xaxis.grid(False)
        ax[ii].set_ylabel(sublabels[ii][1], rotation=90)
        ax[ii].set_ylim(-5.0, 0.0)
        ax[ii].yaxis.set_major_locator(LinearLocator(6))
        ax[ii].yaxis.set_major_formatter(FormatStrFormatter("%.0f"))
        ax[ii].yaxis.grid(False)

        ax[ii].set_title(subtitles[ii], fontsize=fontsize, color=label_color)
        ax[ii].title.set_position([1.0, 1.0])
        ax[ii].set_frame_on(True)

    ax[0].set_xlim(0.0, 0.25)
    ax[1].set_xlim(0.0, 0.25)

    savefig(
        "./figures/bloom/toxin/dissolved.png", facecolor=bg_color, edgecolor="none", dpi=default_dpi
    )
