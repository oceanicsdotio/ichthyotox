#!/usr/bin/python
"""Plot attenuation surfaces"""
from matplotlib.pyplot import figure, savefig, title, LinearLocator, FormatStrFormatter
from numpy import linspace, pi, exp, cos, log, meshgrid
from mpl_toolkits.mplot3d import Axes3D
from analysis.defaults import (
    bg_color,
    overlay_color,
    label_color,
    default_dpi,
    fontsize,
    margin_width,
)

if __name__ == "__main__":

    vpad = 0.1
    hpad = 0.1
    fig = figure(
        facecolor=bg_color, figsize=(margin_width, margin_width / 3.0)
    )  # Change this
    fig.subplots_adjust(
        top=1.0 - vpad, bottom=vpad, left=hpad / 2.0, right=1.0 - hpad / 2.0
    )

    ax: Axes3D = fig.add_subplot(1, 2, 1, projection="3d")
    ax.patch.set_facecolor(bg_color)
    ax.spines["top"].set_color(overlay_color)
    ax.spines["bottom"].set_color(overlay_color)
    ax.spines["left"].set_color(overlay_color)
    ax.spines["right"].set_color(overlay_color)
    ax.xaxis.label.set_color(label_color)
    ax.yaxis.label.set_color(label_color)
    ax.tick_params(axis="x", colors=label_color)
    ax.tick_params(axis="y", colors=label_color)

    # variable declarations
    biomass = linspace(0.01, 0.25, 100)
    time = linspace(0.0, 24.0, 100)
    depth = linspace(-5.0, 0.0, 100)
    depth2 = depth

    time, depth = meshgrid(time, depth)
    irradiance = (325.0 + 325.0 * cos(2.0 * pi * time / 12.0)) * exp(0.15 * depth)
    irradiance[time[:] < 6.0] = 0.0
    irradiance[time[:] > 18.0] = 0.0
    biomass, depth2 = meshgrid(biomass, depth2)
    irradiance2 = (325.0 + 325.0 * cos(0.552 * 12.0)) * exp(
        0.15 * depth2 - 14.0 * biomass
    )

    # surface and contours
    ax.view_init(elev=30, azim=225)
    surf = ax.plot_surface(
        time,
        depth,
        irradiance,
        rstride=5,
        cstride=5,
        color="white",
        edgecolors="black",
        linewidth=1.0,
        shade=False,
        antialiased=True,
    )
    cont1 = ax.contour(
        time, depth, irradiance, colors="black", extend3D=False, linewidth=2, zdir="x"
    )

    # axis adjustment
    ax.set_xlabel(r"t")
    ax.set_xlim(0.0, 24.0)
    ax.xaxis.set_major_locator(LinearLocator(2))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    ax.xaxis.grid(False)
    ax.xaxis.set_rotate_label(False)
    ax.xaxis.set_pane_color(bg_color)

    ax.set_ylabel(r"Z")
    ax.set_ylim(0.0, -5.0)
    ax.yaxis.set_major_locator(LinearLocator(2))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    ax.yaxis.grid(False)
    ax.yaxis.set_rotate_label(False)
    ax.yaxis.set_pane_color(bg_color)

    ax.set_zlabel(r"I")
    ax.set_zlim(0.0, 650.0)
    ax.zaxis.set_major_locator(LinearLocator(2))
    ax.zaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    ax.zaxis.grid(False)
    ax.zaxis.set_rotate_label(False)
    ax.zaxis.set_pane_color(bg_color)

    title(r"(A)", fontsize=fontsize, color=label_color)
    ttl = ax.title
    ttl.set_position([1.0, 1.0])

    ax2: Axes3D = fig.add_subplot(1, 2, 2, projection="3d")
    ax2.view_init(elev=2.0 * 45.0 / 3.0, azim=(225.0 + 180.0))
    surf2 = ax2.plot_surface(
        log(biomass),
        depth2,
        irradiance2,
        rstride=5,
        cstride=5,
        color="white",
        edgecolors="black",
        linewidth=1.0,
        shade=False,
        antialiased=True,
    )
    cont2 = ax2.contour(
        log(biomass),
        depth2,
        irradiance2,
        colors=overlay_color,
        extend3D=False,
        linewidth=2,
        zdir="z",
    )

    ax2.set_xlabel(r"ln B")
    ax2.set_xlim()
    ax2.xaxis.set_major_locator(LinearLocator(2))
    ax2.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    ax2.xaxis.grid(False)
    ax2.xaxis.set_rotate_label(False)
    ax2.xaxis.set_pane_color(bg_color)

    ax2.set_ylabel("Z")
    ax2.set_ylim()
    ax2.yaxis.set_major_locator(LinearLocator(2))
    ax2.yaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    ax2.yaxis.grid(False)
    ax2.yaxis.set_rotate_label(False)
    ax2.yaxis.set_pane_color(bg_color)

    ax2.set_zlabel("I")
    ax2.set_zlim(0.0, 650.0)
    ax2.zaxis.set_major_locator(LinearLocator(2))
    ax2.zaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    ax2.zaxis.grid(False)
    ax2.zaxis.set_rotate_label(False)
    ax2.zaxis.set_pane_color(bg_color)

    title(r"(B)", fontsize=fontsize, color=label_color)
    ttl = ax2.title
    ttl.set_position([1.0, 1.0])

    savefig(
        "./figures/bloom/attenuation.png",
        facecolor=bg_color,
        edgecolor="none",
        dpi=default_dpi,
    )
