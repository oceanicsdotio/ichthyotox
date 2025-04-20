#!/usr/bin/python
"""Analyze vertical migration"""
from sys import argv
from analysis.bloom import Position, plot_and_save_experiments, single_time_series_plot
from analysis.defaults import IMAGE_FORMAT, SOURCE

if __name__ == "__main__":
    directory = argv[1] if len(argv) > 1 else SOURCE
    if len(argv) > 2:
        experiments = argv[2].split(",")
    else:
        experiments = ["100", "101", "102", "103"]
    axes = single_time_series_plot(
        size=(6.5, 3.0),
        padding=(0.1, 0.05),
        ylabel=r"depth (m)",
        yloc=1.0,
    )
    plot_and_save_experiments(
        axes=axes,
        directory=directory,
        experiments=experiments,
        callback=lambda exp: Position(exp).vertical(),
        save_to=f"./figures/bloom/depth.{IMAGE_FORMAT}",
    )
