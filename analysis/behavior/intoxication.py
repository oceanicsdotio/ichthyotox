#!/usr/bin/python
from matplotlib.pyplot import MultipleLocator, FormatStrFormatter, figure, savefig
from numpy import zeros, loadtxt, mean, std
from sys import argv
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi

nplots = 1
fontsize = 10
linewidth = 1
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = 3
default_alpha = 1.0


if __name__ == "__main__":

    data_dir = argv[1]

    # load data
    with open(f"{data_dir}/100/fish_ini.dat", "r", encoding="utf8") as ini:
        nfish = int(str.strip(ini.readline()))
    time = loadtxt(f"{data_dir}/100/fish_state.dat", usecols=[0], unpack=True)
    start = 0
    end = len(time) - 1
    dwidth = 4

    data = zeros((nfish * dwidth + 1, end - start + 1))
    mass = zeros((nfish, end - start + 1))
    tox = zeros((nfish, end - start + 1))
    path = zeros((nfish, end - start + 1))
    ratio = zeros((nfish, end - start + 1))
    ratio_avg = zeros(end - start + 1)
    ratio_avg1 = zeros(end - start + 1)
    ratio_avg2 = zeros(end - start + 1)
    ratio_std = zeros(end - start + 1)
    time = time / 24.0

    # figure and subplots
    fig = figure(facecolor=bg_color, figsize=(marginWidth, fheight))  # Change this
    fig.subplots_adjust(
        top=1.0 - vpadding,
        bottom=vpadding + 0.1,
        left=hpadding,
        right=1.0 - (hpadding / 2.0),
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

    print("Reading State")  # experiment B
    data = loadtxt(f"{data_dir}/101/fish_state.dat", unpack=True)
    for ii in range(0, nfish):
        masscol = ii * dwidth + 2
        toxcol = ii * dwidth + 3
        pathcol = ii * dwidth + 4
        mass[ii, :] = data[masscol, :]
        tox[ii, :] = data[toxcol, :]
        path[ii, :] = data[pathcol, :]
    print("Calculating Mean Tox Load and Pathway Partitioning")
    ratio = tox / mass
    ratio_avg = mean(ratio[:, :], axis=0)
    ratio_std = std(ratio[:, :], axis=0)
    ax.plot(
        time,
        1000 * ratio_avg,
        linestyle="-",
        linewidth=1,
        color="red",
        aa=True,
        label="Formation (B)",
        zorder=2,
    )
    ax.fill_between(
        time,
        1000 * (ratio_avg - ratio_std),
        1000 * (ratio_avg + ratio_std),
        facecolor=[1.0, 0.0, 0.0, 0.25],
        edgecolor="none",
        zorder=2,
    )

    print("Reading State")  # experiment C
    data = loadtxt(f"{data_dir}/102/fish_state.dat", unpack=True)
    for ii in range(0, nfish):
        masscol = ii * dwidth + 2
        toxcol = ii * dwidth + 3
        pathcol = ii * dwidth + 4
        mass[ii, :] = data[masscol, :]
        tox[ii, :] = data[toxcol, :]
        path[ii, :] = data[pathcol, :]
    print("Calculating Mean Tox Load and Pathway Partitioning")
    ratio = tox / mass
    ratio_avg = mean(ratio[:, :], axis=0)
    ratio_std = std(ratio[:, :], axis=0)
    ratio_avg1 = mean(ratio[0 : nfish // 2 - 1, :], axis=0)
    ratio_avg2 = mean(ratio[nfish // 2 : nfish - 1, :], axis=0)
    ax.fill_between(
        time,
        1000 * (ratio_avg - ratio_std),
        1000 * (ratio_avg + ratio_std),
        facecolor=[0.0, 1.0, 0.0, 0.25],
        edgecolor="none",
        zorder=1,
    )
    ax.plot(
        time,
        1000 * ratio_avg1,
        linestyle=":",
        linewidth=1,
        color="green",
        aa=True,
        zorder=1,
    )
    ax.plot(
        time,
        1000 * ratio_avg2,
        linestyle="-",
        linewidth=1,
        color="green",
        aa=True,
        zorder=1,
        label="Intensification (C)",
    )

    print("Reading State")  # experiment D
    data = loadtxt(f"{data_dir}/103/fish_state.dat", unpack=True)
    for ii in range(0, nfish):
        masscol = ii * dwidth + 2
        toxcol = ii * dwidth + 3
        pathcol = ii * dwidth + 4
        mass[ii, :] = data[masscol, :]
        tox[ii, :] = data[toxcol, :]
        path[ii, :] = data[pathcol, :]
    print("Calculating Mean Tox Load and Pathway Partitioning")
    ratio = tox / mass
    ratio_avg = mean(ratio[:, :], axis=0)
    ratio_std = std(ratio[:, :], axis=0)
    ratio_avg1 = mean(ratio[0 : nfish // 2 - 1, :], axis=0)
    ratio_avg2 = mean(ratio[nfish // 2 : nfish - 1, :], axis=0)
    ax.fill_between(
        time,
        1000 * (ratio_avg - ratio_std),
        1000 * (ratio_avg + ratio_std),
        facecolor=[0.0, 0.0, 1.0, 0.25],
        edgecolor="none",
        zorder=3,
    )
    ax.plot(
        time,
        1000 * ratio_avg1,
        linestyle=":",
        linewidth=1,
        color="blue",
        aa=True,
        zorder=3,
    )
    ax.plot(
        time,
        1000 * ratio_avg2,
        linestyle="-",
        linewidth=1,
        color="blue",
        aa=True,
        zorder=3,
        label="Decline (D)",
    )

    # intoxication threshold
    ax.plot(
        (time[start], time[end]),
        (0.005, 0.005),
        linestyle="--",
        linewidth=1,
        color="black",
        aa=True,
        zorder=5,
    )
    ax.plot(
        (time[start], time[end]),
        (0.015, 0.015),
        linestyle="--",
        linewidth=1,
        color="black",
        aa=True,
        zorder=5,
    )

    ax.set_xlabel(r"Days")
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    ax.set_xlim(0, 30)
    ax.xaxis.grid(False)

    ax.set_ylabel("Intoxication (ug/g)", rotation=90)

    ax.yaxis.set_major_locator(MultipleLocator(0.01))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.yaxis.grid(False)
    ax.set_frame_on(True)
    ax.set_ylim(
        0,
    )

    legend = ax.legend(loc="upper left")
    legend.get_frame().set_facecolor("none")
    legend.get_frame().set_edgecolor("none")
    text = legend.get_texts()
    text[0].set_color("red")
    text[1].set_color("green")
    text[2].set_color("blue")

    savefig(
        "./figures/behavior/intoxication.png",
        facecolor=bg_color,
        edgecolor="none",
        dpi=default_dpi,
    )
