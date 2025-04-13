#!/usr/bin/python
import sys, os
sys.path.append('Dropbox/UMaine/2017_Sandbox/Source')

from plotlib import DataView
from particles import ParticleSimulation
from numpy import sum, linspace

DATA_DIR = 'Dropbox/UMaine/2017_EcologicalComplexity/Data/'
volume = 500.0*500.0*5.0
area = 500.*500.
nlayers = 26
window = 40 # steps = 4 hours

params = ['protein', 'carbohydrate', 'microcystin']
SIM_DIRS = ['100/', '101/', '102/', '103/']
toffsets = [0, 30, 60, 90 ]

# Setup figure for plotting objects
view = DataView(fwidth=3.25,fheight=1.5, verb=True)
for simulation, offset in zip(SIM_DIRS, toffsets):
    
    data = ParticleSimulation('fish', DATA_DIR+simulation, verb=True)
    
    sur, dem = data.calc_diffusivity(window)
    view.ax.plot(linspace(1.0,30.0,30)+offset, sur, color='white')

view.fshow(ylab='K (m^2/s)', xloc=30, yloc=2, geo=False)