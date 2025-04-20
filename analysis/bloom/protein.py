#!/usr/bin/python
"""Generate figures for bloom protein analysis."""
from sys import argv
from analysis.bloom import State, single_time_series_plot, plot_and_save_experiments
from analysis.defaults import IMAGE_FORMAT, SOURCE

DEST = "./figures/bloom"
VOLUME = 500.0 * 500.0 * 5.0
EXPERIMENTS = ["100", "101", "102", "103"]

if __name__ == "__main__":
    DIR = argv[1] if len(argv) > 1 else SOURCE
    if len(argv) > 2:
        experiments = argv[2].split(",")
    else:
        experiments = EXPERIMENTS
    ax = single_time_series_plot(
        size=(6.5, 3.0),
        padding=(0.15, 0.1),
        ylabel=r"protein (g/$m^3$)",
        yloc=0.05,
    )
    plot_and_save_experiments(
        axes=ax,
        directory=DIR,
        experiments=experiments,
        callback=lambda exp: State(exp).protein_concentration(VOLUME),
        save_to=f"{DEST}/protein.{IMAGE_FORMAT}",
    )

