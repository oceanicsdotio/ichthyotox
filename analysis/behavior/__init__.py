#!/usr/bin/python
from matplotlib.pyplot import tripcolor, get_cmap, colorbar, scatter
from matplotlib.tri import Triangulation
from numpy import sin, pi, loadtxt, zeros, sqrt, ndarray


def suitability(x):
    """Calculate suitability"""
    return 0.5 * (1.0 + sin(2.0 * pi * (x - 125.0) / 500.0))


class State:
    """Fish particle state data"""

    def __init__(self, experiment: str, fmt="dat", width=4):
        with open(f"{experiment}/fish_ini.{fmt}", "r", encoding="utf8") as ini:
            self.count = int(str.strip(ini.readline()))
                # load data
        time = loadtxt(f'{experiment}/fish_state.{fmt}', usecols=[0], unpack=True)
        start = 0
        end = len(time) - 1
        data = loadtxt(f'{experiment}/fish_state.{fmt}', unpack=True)
        self.days = time / 24.0
        self.mass = zeros((self.count, end-start+1))
        for ii in range(0, self.count):
            masscol = ii * width + 2
            self.mass[ii,:] = data[masscol, :]
    
    def mean_mass(self) -> ndarray[float]:
        """Return mean mass of fish"""
        return self.mass.mean(axis=0)

class Domain:
    """Horizontal domain"""

    def __init__(self, experiment: str, fmt="dat"):
        self.vert_x, self.vert_y = loadtxt(
            f"{experiment}/mesh_node.{fmt}", skiprows=1, usecols=(1, 2), unpack=True
        )
        self.ind1, self.ind2, self.ind3 = loadtxt(
            f"{experiment}/mesh_elem.{fmt}",
            dtype="i8",
            skiprows=1,
            usecols=(1, 2, 3),
            unpack=True,
        )
        self.triang = Triangulation(self.vert_x, self.vert_y)  # Delauney triangulation
        self.suitability = suitability(self.vert_x)  # calculate suitability

    def x_range(self):
        """Range of values for axis"""
        return self.vert_x.min(), self.vert_x.max()

    def y_range(self):
        """Range of values for axis"""
        return self.vert_y.min(), self.vert_y.max()

    def draw(self):
        """Plot suitability surface"""
        tripcolor(self.triang, self.suitability, shading="flat", cmap=get_cmap("Blues"))
        colorbar()

        # Plot grid
        # for ii in range(0,len(ind1)):
        #     plot( (vert_x[ind1[ii]-1], vert_x[ind2[ii]-1]), (vert_y[ind1[ii]-1], vert_y[ind2[ii]-1]), color='grey', linewidth=0.5, aa=True, zorder=2)
        #     plot( (vert_x[ind2[ii]-1], vert_x[ind3[ii]-1]), (vert_y[ind2[ii]-1], vert_y[ind3[ii]-1]), color='grey', linewidth=0.5, aa=True, zorder=2)
        #     plot( (vert_x[ind3[ii]-1], vert_x[ind1[ii]-1]), (vert_y[ind3[ii]-1], vert_y[ind1[ii]-1]), color='grey', linewidth=0.5, aa=True, zorder=2)


class Positions:
    """Fish particle position data"""
    def __init__(self, experiment: str, fmt="dat", columns=4):

        with open(f"{experiment}/fish_ini.{fmt}", "r", encoding="utf8") as ini:
            self.count = int(str.strip(ini.readline()))
        time = loadtxt(f"{experiment}/fish_position.{fmt}", usecols=[0], unpack=True)
        start = 0
        end = len(time) - 1
        shape = (self.count, end - start + 1)
        data = zeros((self.count * columns + 1, end - start + 1))
        self.x = zeros(shape)
        self.y = zeros(shape)
        data = loadtxt(f"{experiment}/fish_position.{fmt}", unpack=True)
        for ii in range(0, self.count):
            xcol = ii * columns + 2
            ycol = ii * columns + 3
            self.x[ii, :] = data[xcol, :]
            self.y[ii, :] = data[ycol, :]

    def draw_final(self):
        scatter(self.x[:,-1], self.y[:,-1], s=40, color='black', zorder=10, edgecolors='face') # end markers

    def draw_trajectory(self, ax):
        end = self.x.shape[1] - 1
        for ii in range(0, self.count): 
            for jj in range (0, end):
                if (sqrt( (self.x[ii,jj+1] - self.x[ii,jj])**2 + (self.y[ii,jj+1] - self.y[ii,jj])**2 ) < 250.0):
                    ax.plot( (self.x[ii,jj], self.x[ii,jj+1]), (self.y[ii,jj], self.y[ii,jj+1]), linewidth=1, linestyle='-', color='black', alpha=0.2, aa=True, zorder=7)

    def draw(self, ax, color):
        """Draw fish positions as scatter plot"""
        mid = 22 * 10 * 24
        ax[0].scatter(
            self.x[0:24, mid],
            self.y[0:24, mid],
            s=25,
            color=color,
            zorder=10,
            edgecolors="face",
        )  # mid markers
        ax[0].scatter(
            self.x[100:124, mid],
            self.y[100:124, mid],
            s=25,
            color=color,
            zorder=10,
            edgecolors="face",
        )  # mid markers
        ax[1].scatter(
            self.x[0:24, -1],
            self.y[0:24, -1],
            s=25,
            color="black",
            zorder=10,
            edgecolors="face",
        )  # final markers
        ax[1].scatter(
            self.x[100:124, -1],
            self.y[100:124, -1],
            s=25,
            color="black",
            zorder=10,
            edgecolors="face",
        )  # final markers
