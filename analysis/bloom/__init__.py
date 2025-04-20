#!/usr/bin/python
"""Bloom Carbon Analysis"""
from sys import maxsize, float_info
from numpy import ndarray, zeros, loadtxt
from matplotlib.pyplot import figure, savefig, Line2D, MultipleLocator, FormatStrFormatter, Axes
from analysis.defaults import bg_color, overlay_color, label_color, linewidth, overlay, fill_color, default_dpi

def merge_range(previous: tuple, current: tuple) -> tuple:
    """
    Merge two ranges into a single range.
    """
    lower = min(previous[0], current[0])
    upper = max(previous[1], current[1])
    return lower, upper

def surface_plot():
    pass

def single_time_series_plot(
    size: tuple[float, float], padding: tuple[float, float], ylabel: str, yloc: float
) -> Axes:
    """Multiple series on a single plot"""
    fig = figure(facecolor=bg_color, figsize=size)
    hpadding, vpadding = padding
    fig.subplots_adjust(
        top=1.0 - vpadding,
        bottom=vpadding + 0.1,
        left=hpadding,
        right=1.0 - (hpadding / 2.0),
        hspace=0.20,
    )
    ax = fig.add_subplot(1, 1, 1)
    ax.patch.set_facecolor(bg_color)
    ax.spines["top"].set_color(overlay_color)
    ax.spines["bottom"].set_color(overlay_color)
    ax.spines["left"].set_color(overlay_color)
    ax.spines["right"].set_color(overlay_color)
    ax.xaxis.label.set_color(label_color)
    ax.yaxis.label.set_color(label_color)
    ax.tick_params(axis="x", colors=label_color)
    ax.tick_params(axis="y", colors=label_color)
    ax.xaxis.grid(False)
    ax.set_xlabel(r"Days")
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    ax.set_ylabel(ylabel, rotation=90)
    ax.yaxis.set_major_locator(MultipleLocator(yloc))
    # ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.yaxis.grid(False)
    ax.set_frame_on(True)
    return ax

def plot_and_save_experiments(
    axes: Axes, directory: str, experiments: list[str], callback, save_to: str
) -> None:
    """
    Plot and save the experiments."""
    drawn = []
    time_range = (maxsize, 0)
    value_range = (float_info.max, float_info.min)
    for experiment, color, fill in zip(experiments, overlay, fill_color):
        envelope: Envelope = callback(directory + experiment)
        time_range = merge_range(time_range, envelope.time_range())
        value_range = merge_range(value_range, envelope.value_range())
        line = envelope.plot(axes, color=color, fill=fill, with_bounds=True)
        drawn.append(line)

    axes.set_xlim(*time_range)
    axes.set_ylim(*value_range)
    legend = axes.legend(drawn, experiments, loc="best")
    frame = legend.get_frame()
    frame.set_facecolor("none")
    frame.set_edgecolor("none")
    savefig(
        save_to,
        facecolor=bg_color,
        edgecolor="none",
        dpi=default_dpi,
    )

class Envelope:
    """
    Class to calculate and plot the mean and standard deviation of a series of data.
    """

    def __init__(self, time: ndarray[float], series: ndarray[float], positive=True):
        self.time = time
        self.mean: ndarray[float] = series.mean(axis=0)
        self.std_dev: ndarray[float] = series.std(axis=0)
        self.upper: ndarray[float] = self.mean + self.std_dev
        self.lower: ndarray[float] = self.mean - self.std_dev
        if positive:
            self.lower = self.lower.clip(min=0.0)

    def time_range(self) -> tuple[float, float]:
        """
        Get the time range of the data.
        """
        return self.time[0], self.time[-1]

    def value_range(self) -> tuple[float, float]:
        """
        Get the range of values in the series.
        """
        upper: float = (self.mean + self.std_dev).max()
        lower: float = (self.mean - self.std_dev).min()
        return lower, upper

    def plot(self, subplot: Axes, color, fill, with_bounds=True) -> Line2D:
        """
        Plot the mean and standard deviation of the data series.
        """
        line = subplot.plot(
            self.time,
            self.mean,
            linestyle="-",
            linewidth=linewidth,
            color=color,
            aa=True,
        )
        if with_bounds:
            # subplot.fill_between(
            #     self.time,
            #     self.upper,
            #     self.lower,
            #     facecolor=fill,
            #     edgecolor="none",
            #     zorder=3,
            # )
            _ = subplot.plot(
                self.time,
                self.upper,
                linestyle="--",
                linewidth=linewidth,
                color=color,
                aa=True
            )
            _ = subplot.plot(
                self.time,
                self.lower,
                linestyle="--",
                linewidth=linewidth,
                color=color,
                aa=True
            )
        return line[0]

class Position:
    """
    Movement and position data of bloom particle system
    """

    def __init__(self, experiment, file_format="dat", columns=4):
        ini_file = f"{experiment}/cyanobacteria_ini.{file_format}"
        filename = f"{experiment}/cyanobacteria_position.{file_format}"
        with open(ini_file, "r", encoding="utf-8") as initial_conditions:
            count = int(str.strip(initial_conditions.readline()))
        hours = loadtxt(
            filename,
            usecols=[0],
            unpack=True,
        )
        steps = len(hours)
        shape = (count, steps)
        data = zeros((count * columns + 1, steps))
        data = loadtxt(filename, unpack=True)

        self.count = count
        self.days = hours / 24.0
        self.z = zeros(shape)

        for ii in range(0, count):
            depth_column = ii * columns + 4
            self.z[ii, :] = data[depth_column, :]

    def vertical(self):
        """
        Vertical position envelope
        """
        return Envelope(self.days, self.z)

class State:
    """
    Carbon and toxin state of the bloom system.
    """

    def __init__(self, experiment, file_format="dat", columns=4):
        ini_file = f"{experiment}/cyanobacteria_ini.{file_format}"
        filename = f"{experiment}/cyanobacteria_state.{file_format}"
        with open(ini_file, "r", encoding="utf-8") as initial_conditions:
            count = int(str.strip(initial_conditions.readline()))

        hours = loadtxt(
            filename,
            usecols=[0],
            unpack=True,
        )

        steps = len(hours)
        shape = (count, steps)
        data = zeros((count * columns + 1, steps))
        data = loadtxt(filename, unpack=True)

        self.count = count
        self.days = hours / 24.0
        self.carbohydrate = zeros(shape)
        self.protein = zeros(shape)
        self.toxin = zeros(shape)

        for ii in range(0, count):
            carbohydrate_column = ii * columns + 2
            protein_column = ii * columns + 3
            toxin_column = ii * columns + 4
            self.carbohydrate[ii, :] = data[carbohydrate_column, :]
            self.protein[ii, :] = data[protein_column, :]
            self.toxin[ii, :] = data[toxin_column, :]

    def carbon_concentration(self, volume: float) -> Envelope:
        """
        Calculate the total carbon concentration in the system.
        """
        concentration = (self.protein + self.carbohydrate) / volume
        return Envelope(self.days, concentration)

    def carbon_ratio(self) -> Envelope:
        """
        Calculate the ratio of carbohydrate to protein.
        """
        ratio = self.carbohydrate / self.protein
        return Envelope(self.days, ratio)

    def protein_concentration(self, volume: float) -> Envelope:
        """
        Calculate the concentration of protein in the system.
        """
        concentration = self.protein / volume
        return Envelope(self.days, concentration)
    
    def toxicity(self) -> Envelope:
        """
        Calculate the concentration of toxin in the system.
        """
        toxicity = self.toxin / (self.carbohydrate + self.protein)
        return Envelope(self.days, toxicity)
