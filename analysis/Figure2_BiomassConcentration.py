#!/usr/bin/python
import sys, os
sys.path.append('Dropbox/UMaine/2017_Sandbox/Source')

from plotlib import DataView
from particles import ParticleSimulation
from numpy import sum

DATA_DIR = 'Dropbox/UMaine/2017_EcologicalComplexity/Data/'
volume = 500.0*500.0*5.0
params = ['protein', 'carbohydrate']
SIM_DIRS = ['100/', '101/', '102/', '103/']
toffsets = [0, 30, 60, 90 ]

# Setup figure for plotting objects
view = DataView(fwidth=3.25,fheight=1.5, verb=True)
for simulation, offset in zip(SIM_DIRS, toffsets):
    data = ParticleSimulation('cyanobacteria', \
        DATA_DIR+simulation, params, verb=True)
    protein = data._ParticleSimulation__state[0]
    carbs = data._ParticleSimulation__state[1]
    biomass = sum(protein + carbs, axis=0)/volume
    data.show(view, biomass, toffset=offset)

view.fshow(ylab='mgC/L', xloc=30, yloc=10, geo=False)