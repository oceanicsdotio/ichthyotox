#!/usr/bin/python
"""Combined Protein and Carbohydrate Analysis"""
from sys import argv
from analysis.defaults import IMAGE_FORMAT, SOURCE
from analysis.bloom import State, single_time_series_plot, plot_and_save_experiments

if __name__ == "__main__":
    directory = argv[1] if len(argv) > 1 else SOURCE
    if len(argv) > 2:
        experiments = argv[2].split(",")
    else:
        experiments = ["100", "101", "102", "103"]
    VOLUME = 500.0 * 500.0 * 5.0
    axes = single_time_series_plot(
        size=(6.5, 3.0),
        padding=(0.15, 0.1),
        ylabel=r"carbon (g/$m^3$)",
        yloc=0.1,
    )
    plot_and_save_experiments(
        axes=axes,
        directory=directory,
        experiments=experiments,
        callback=lambda exp: State(exp).carbon_concentration(VOLUME),
        save_to=f"./figures/bloom/carbon.{IMAGE_FORMAT}",
    )
