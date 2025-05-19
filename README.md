# Ichthyotox

Drawing conclusions from drunk fish in dynamic environments.

## Overview

This repository contains code to simulate and analyze cyanobacteria blooms and the behavioral response of individual fish to toxins.

The simulation code written in Fortran, and can be be found in the `src/` directory. This is compiled using `gfortran`, with the resulting binaries output to `bin/`. Running the compiled programs will produce data in the `data/` directory.

Analysis and visualization code is written in Python, with dependencies managed with `pixi` for access to scientific computing libraries without directly using `conda`. The Python code can be found in `analysis/`. The visualizations are created in `figures/`.

For testing and post-processing, some of the Fortran modules are compiled with `f2py` and `meson` to create bindings so that numerical functions can be called from Python.

## Programs

### Random numbers

The program generates pseudorandom numbers with either uniform or normal distributions using the Box-Muller method. Randoms numbers are used for initial conditions and for movement algorithms for both blooms and fish.

### Forcing

You can `make bin/forcing` to compile a program that will let you input simulation parameters and create
a new experiment file. The command `make new` will create a test experiment.

### Blooms

Cyanobacteria blooms can be simulated without fish behavior.


