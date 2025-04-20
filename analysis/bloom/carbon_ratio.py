#!/usr/bin/python
"""Plot Carbohydrate/Protein Ratio Analysis"""
from sys import argv
from analysis.bloom import State, plot_and_save_experiments, single_time_series_plot
from analysis.defaults import IMAGE_FORMAT, SOURCE

DEST = "./figures/bloom"

if __name__ == "__main__":
    DATA_DIR = argv[1] if len(argv) > 1 else SOURCE
    if len(argv) > 2:
        experiments = argv[2].split(",")
    else:
        experiments = ["100", "101", "102", "103"]
    axes = single_time_series_plot(
        size=(6.5, 3.0),
        padding=(0.15, 0.1),
        ylabel=r"carbon ratio",
        yloc=1.0,
    )
    plot_and_save_experiments(
        axes,
        directory=DATA_DIR,
        experiments=experiments,
        callback=lambda source: State(source).carbon_ratio(),
        save_to=f"{DEST}/carbon_ratio.{IMAGE_FORMAT}",
    )
