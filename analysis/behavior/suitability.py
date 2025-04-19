#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
from pylab import *
from sys import argv

nplots = 1
fontsize = 10
linewidth = 1
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = 3.0
default_alpha = 1.0

for_screen = False
if for_screen:
    bg_color = [0.0, 0.0, 0.0, default_alpha]
    overlay_color = [1.0, 1.0, 1.0, default_alpha]
    label_color = [1.0, 1.0, 1.0, default_alpha]
    default_dpi = 150
else:
    bg_color = [1.0, 1.0, 1.0, default_alpha]
    overlay_color = [0.0, 0.0, 0.0, default_alpha]
    label_color = [0.0, 0.0, 0.0, default_alpha]
    default_dpi = 300

if __name__ == "__main__":

    data_dir = argv[1]

    rc("text", usetex=False)
    rc("font", **{"family": "serif", "serif": ["Times New Roman"]})
    rc("mathtext", default="sf")
    rc("lines", markeredgewidth=1)
    rc("lines", linewidth=linewidth)
    rc("axes", labelsize=fontsize)
    rc("axes", linewidth=(linewidth + 1) // 2)
    rc("xtick", labelsize=fontsize)
    rc("ytick", labelsize=fontsize)
    rc("legend", fontsize=fontsize)
    rc("xtick.major", pad=5)
    rc("ytick.major", pad=5)

    # load data
    with open(f"{data_dir}/100/fish_ini.dat", "r", encoding="utf8") as ini:
        nfish = int(str.strip(ini.readline()))
    time = np.loadtxt(f"{data_dir}/100/fish_state.dat", usecols=[0], unpack=True)
    start = 0
    end = len(time) - 1
    dwidth = 4

    data = zeros((nfish * dwidth + 1, end - start + 1))
    mass = zeros((nfish, end - start + 1))
    suitability = zeros((nfish, end - start + 1))
    xpos = zeros((nfish, end - start + 1))
    mass_avg1 = zeros(end - start + 1)
    mass_avg2 = zeros(end - start + 1)

    suit_avg = zeros(end - start + 1)
    suit_smoothed = zeros(end - start + 1)
    suit_smoothed_s = zeros(end - start + 1)
    suit_smoothed_b = zeros(end - start + 1)
    suit_avg_surface = zeros(end - start + 1)
    suit_avg_bottom = zeros(end - start + 1)
    mass_std = zeros(end - start + 1)
    suit_std = zeros(end - start + 1)
    time = time / 24.0
    window = 10  # in both directions, 2 is therefore 5

    # figure and subplots
    fig = plt.figure(facecolor=bg_color, figsize=(marginWidth, fheight))  # Change this
    fig.subplots_adjust(
        top=1.0 - vpadding,
        bottom=vpadding + 0.1,
        left=hpadding,
        right=1.0 - hpadding / 2.0,
    )
    ax = []
    ax.append(fig.add_subplot(1, 1, 1))

    ax[0].patch.set_facecolor(bg_color)
    ax[0].spines["top"].set_color(overlay_color)
    ax[0].spines["bottom"].set_color(overlay_color)
    ax[0].spines["left"].set_color(overlay_color)
    ax[0].spines["right"].set_color(overlay_color)
    ax[0].xaxis.label.set_color(label_color)
    ax[0].yaxis.label.set_color(label_color)
    ax[0].tick_params(axis="x", colors=label_color)
    ax[0].tick_params(axis="y", colors=label_color)

    data = np.loadtxt(f"{data_dir}/100/fish_position.dat", unpack=True)
    for ii in range(0, nfish):
        xcol = ii * dwidth + 2
        xpos[ii, :] = data[xcol, :]
    suitability = 0.5 * (1.0 + sin(2.0 * pi * (xpos[:, :] - 125.0) / 500.0))
    suit_avg = np.mean(suitability[:, :], axis=0)
    for ii in range(0, 2 * window + 1):
        suit_smoothed[start + window : end - window] = suit_smoothed[
            start + window : end - window
        ] + suit_avg[start + ii : end + ii - 2 * window] / (2 * window + 1)
    ax[0].plot(
        time[start + window : end - window],
        suit_smoothed[start + window : end - window],
        linestyle="-",
        linewidth=1,
        color="black",
        aa=True,
        label="Control (A)",
        zorder=3,
    )

    data = np.loadtxt(f"{data_dir}/101/fish_position.dat", unpack=True)
    for ii in range(0, nfish):
        xcol = ii * dwidth + 2
        xpos[ii, :] = data[xcol, :]
    suitability = 0.5 * (1.0 + sin(2.0 * pi * (xpos[:, :] - 125.0) / 500.0))
    suit_avg = np.mean(suitability[:, :], axis=0)
    suit_smoothed = zeros(end - start + 1)
    for ii in range(0, 2 * window + 1):
        suit_smoothed[start + window : end - window] = suit_smoothed[
            start + window : end - window
        ] + suit_avg[start + ii : end + ii - 2 * window] / (2.0 * window + 1.0)
    ax[0].plot(
        time[start + window : end - window],
        suit_smoothed[start + window : end - window],
        linestyle="-",
        linewidth=1,
        color="red",
        aa=True,
        label="Formation (B)",
        zorder=4,
    )

    data = np.loadtxt(f"{data_dir}/102/fish_position.dat", unpack=True)
    for ii in range(0, nfish):
        xcol = ii * dwidth + 2
        xpos[ii, :] = data[xcol, :]
    suitability = 0.5 * (1.0 + sin(2.0 * pi * (xpos[:, :] - 125.0) / 500.0))
    suit_smoothed = zeros(end - start + 1)

    suit_avg_surface = np.mean(suitability[0 : nfish // 2 - 1, :], axis=0)
    for ii in range(0, 2 * window + 1):
        suit_smoothed_s[start + window : end - window] = suit_smoothed_s[
            start + window : end - window
        ] + suit_avg_surface[start + ii : end + ii - 2 * window] / (2.0 * window + 1.0)

    suit_avg_bottom = np.mean(suitability[nfish // 2 : nfish - 1, :], axis=0)
    for ii in range(0, 2 * window + 1):
        suit_smoothed_b[start + window : end - window] = suit_smoothed_b[
            start + window : end - window
        ] + suit_avg_bottom[start + ii : end + ii - 2 * window] / (2.0 * window + 1.0)

    ax[0].plot(
        time[start + window : end - window],
        suit_smoothed_s[start + window : end - window],
        linestyle=":",
        linewidth=1,
        color="green",
        aa=True,
    )
    ax[0].plot(
        time[start + window : end - window],
        suit_smoothed_b[start + window : end - window],
        linestyle="-",
        linewidth=1,
        color="green",
        aa=True,
        label="Intensification (C)",
    )

    data = np.loadtxt(f"{data_dir}/103/fish_position.dat", unpack=True)
    for ii in range(0, nfish):
        xcol = ii * dwidth + 2
        xpos[ii, :] = data[xcol, :]
    suitability = 0.5 * (1.0 + sin(2.0 * pi * (xpos[:, :] - 125.0) / 500.0))
    suit_smoothed = zeros(end - start + 1)
    suit_smoothed_s = zeros(end - start + 1)
    suit_smoothed_b = zeros(end - start + 1)

    suit_avg_surface = np.mean(suitability[0 : nfish // 2 - 1, :], axis=0)
    for ii in range(0, 2 * window + 1):
        suit_smoothed_s[start + window : end - window] = suit_smoothed_s[
            start + window : end - window
        ] + suit_avg_surface[start + ii : end + ii - 2 * window] / (2.0 * window + 1.0)

    suit_avg_bottom = np.mean(suitability[nfish // 2 : nfish - 1, :], axis=0)
    for ii in range(0, 2 * window + 1):
        suit_smoothed_b[start + window : end - window] = suit_smoothed_b[
            start + window : end - window
        ] + suit_avg_bottom[start + ii : end + ii - 2 * window] / (2.0 * window + 1.0)

    ax[0].plot(
        time[start + window : end - window],
        suit_smoothed_s[start + window : end - window],
        linestyle=":",
        linewidth=1,
        color="blue",
        aa=True,
    )
    ax[0].plot(
        time[start + window : end - window],
        suit_smoothed_b[start + window : end - window],
        linestyle="-",
        linewidth=1,
        color="blue",
        aa=True,
        label="Decline (D)",
    )

    ax[0].plot(
        (time[start], time[end]),
        (0.5, 0.5),
        linestyle="--",
        linewidth=1,
        color="black",
        aa=True,
        zorder=5,
    )

    ax[0].set_xlabel(r"Days")
    ax[0].xaxis.set_major_locator(MultipleLocator(5))
    ax[0].xaxis.set_minor_locator(MultipleLocator(1))
    ax[0].xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    ax[0].set_xlim(0, 30)
    ax[0].xaxis.grid(False)

    ax[0].set_ylabel("Suitability", rotation=90)

    ax[0].yaxis.set_major_locator(MultipleLocator(0.1))
    ax[0].yaxis.set_minor_locator(MultipleLocator(0.05))
    ax[0].yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax[0].yaxis.grid(False)
    ax[0].set_frame_on(True)
    ax[0].set_ylim(0.3, 0.6)

    legend = ax[0].legend(loc="lower left")
    legend.get_frame().set_facecolor("none")
    legend.get_frame().set_edgecolor("none")
    text = legend.get_texts()
    text[0].set_color("black")
    text[1].set_color("red")
    text[2].set_color("green")
    text[3].set_color("blue")

    # plt.show()
    plt.savefig(
        "./figures/behavior/suitabililty.png",
        facecolor=bg_color,
        edgecolor="none",
        dpi=default_dpi,
    )
