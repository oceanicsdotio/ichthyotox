#!/usr/bin/python
"""Generate figures for bloom protein analysis."""
from analysis.bloom import State, single_time_series_plot, plot_and_save_experiments

if __name__ == "__main__":
    DATA_DIR = "/Users/keeney/Documents/Projects/2020_2017_Ichthyotox/data/"
    experiments = ["100", "101", "102", "103"]
    VOLUME = 500.0 * 500.0 * 5.0

    ax = single_time_series_plot(
        size=(6.5, 3.0),
        padding=(0.1, 0.05),
        ylabel=r"Mean Protein (g/m^3)",
        yloc=0.05,
    )
    plot_and_save_experiments(
        axes=ax,
        directory=DATA_DIR,
        experiments=experiments,
        callback=lambda exp: State(exp).protein_concentration(VOLUME),
        save_to="./figures/bloom/protein.png",
    )

