#!/usr/bin/python
"""Particulate Toxin Analysis"""
from analysis.bloom import State, single_time_series_plot, plot_and_save_experiments

if __name__ == "__main__":
    DATA_DIR = "/Users/keeney/Documents/Projects/2020_2017_Ichthyotox/data/"
    sources = ["100", "101", "102", "103"]
    ax = single_time_series_plot(
        size=(6.5, 3.0),
        padding=(0.1, 0.05),
        ylabel=r"Toxicity (g/g)",
        yloc=0.1,
    )
    plot_and_save_experiments(
        axes=ax,
        directory=DATA_DIR,
        experiments=sources,
        callback=lambda exp: State(exp).toxicity(),
        save_to="./figures/bloom/toxicity.png",
    )
