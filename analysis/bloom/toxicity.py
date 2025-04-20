#!/usr/bin/python
"""Particulate Toxin Analysis"""
from sys import argv
from analysis.bloom import State, single_time_series_plot, plot_and_save_experiments
from analysis.defaults import IMAGE_FORMAT

DEFAULT_SOURCE = "./data/"
DEST = "./figures/bloom"

if __name__ == "__main__":
    directory = argv[1] if len(argv) > 1 else DEFAULT_SOURCE
    sources = ["100", "101", "102", "103"]
    ax = single_time_series_plot(
        size=(6.5, 3.0),
        padding=(0.15, 0.1),
        ylabel=r"toxicity",
        yloc=0.1,
    )
    plot_and_save_experiments(
        axes=ax,
        directory=directory,
        experiments=sources,
        callback=lambda exp: State(exp).toxicity(),
        save_to=f"{DEST}/toxicity.png",
    )
