#!/usr/bin/python
"""Plot Carbohydrate/Protein Ratio Analysis"""
from analysis.bloom import State, plot_and_save_experiments, single_time_series_plot

if __name__ == "__main__":
    DATA_DIR = "/Users/keeney/Documents/Projects/2020_2017_Ichthyotox/data/"
    experiments = ["100", "101", "102", "103"]
    axes = single_time_series_plot(
        size=(6.5, 3.0),
        padding=(0.1, 0.05),
        ylabel="Carbon Ratio (g/g)",
        yloc=1.0,
    )
    plot_and_save_experiments(
        axes,
        directory=DATA_DIR,
        experiments=experiments,
        callback=lambda source: State(source).carbon_ratio(),
        save_to="./figures/bloom/carbon_ratio.png",
    )
