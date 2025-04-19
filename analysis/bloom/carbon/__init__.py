#!/usr/bin/python
"""Bloom Carbon Analysis"""
import matplotlib.ticker as mticker
import numpy as np
from pylab import Axes
from analysis.defaults import linewidth


class Envelope:
    """
    Class to calculate and plot the mean and standard deviation of a series of data.
    """

    def __init__(self, series: np.ndarray[float]):
        self.mean = series.mean(axis=0)
        self.std_dev = series.std(axis=0)
        self.upper = self.mean + self.std_dev
        self.lower = self.mean - self.std_dev

    def plot(self, subplot: Axes, time, color, fill, with_bounds=True):
        """
        Plot the mean and standard deviation of the data series.
        """
        subplot.plot(
            time,
            self.mean,
            linestyle="-",
            linewidth=linewidth,
            color=color,
            aa=True,
            label="Control (A)",
        )
        if with_bounds:
            subplot.fill_between(
                time, self.upper, self.lower, facecolor=fill, edgecolor="none", zorder=3
            )


class State:
    """
    Carbon state of the system.
    """

    def __init__(self, experiment, file_format="dat", columns=4):
        ini_file = f"{experiment}/cyanobacteria_ini.{file_format}"
        filename = f"{experiment}/cyanobacteria_state.{file_format}"
        with open(ini_file, "r", encoding="utf-8") as initial_conditions:
            count = int(str.strip(initial_conditions.readline()))

        hours = np.loadtxt(
            filename,
            usecols=[0],
            unpack=True,
        )

        steps = len(hours)
        shape = (count, steps)
        data = np.zeros((count * columns + 1, steps))
        data = np.loadtxt(filename, unpack=True)

        self.count = count
        self.days = hours / 24.0
        self.carbohydrate = np.zeros(shape)
        self.protein = np.zeros(shape)

        for ii in range(0, count):
            carbohydrate_column = ii * columns + 2
            protein_column = ii * columns + 3
            self.carbohydrate[ii, :] = data[carbohydrate_column, :]
            self.protein[ii, :] = data[protein_column, :]

    def carbon_concentration(self, volume: float) -> Envelope:
        """
        Calculate the total carbon concentration in the system.
        """
        concentration = (self.protein + self.carbohydrate) / volume
        return Envelope(concentration)

    def carbon_ratio(self) -> Envelope:
        """
        Calculate the ratio of carbohydrate to protein.
        """
        ratio = self.carbohydrate / self.protein
        return Envelope(ratio)

    def protein_concentration(self, volume: float) -> Envelope:
        """
        Calculate the concentration of protein in the system.
        """
        concentration = self.protein / volume
        return Envelope(concentration)
