#!/usr/bin/python
"""Combined Protein and Carbohydrate Analysis"""
from analysis.bloom import State, single_time_series_plot, plot_and_save_experiments

if __name__ == "__main__":
    DATA_DIR = "/Users/keeney/Documents/Projects/2020_2017_Ichthyotox/data/"
    VOLUME = 500.0 * 500.0 * 5.0
    experiments = ["100", "101", "102", "103"]
    axes = single_time_series_plot(
        size=(6.5, 3.0),
        padding=(0.1, 0.05),
        ylabel=r"Total Carbon (g/m^3)",
        yloc=0.1,
    )
    plot_and_save_experiments(
        axes=axes,
        directory=DATA_DIR,
        experiments=experiments,
        callback=lambda exp: State(exp).carbon_concentration(VOLUME),
        save_to="./figures/bloom/carbon.png",
    )
