#!/usr/bin/python
import sys, os
sys.path.append('Dropbox/UMaine/2017_Sandbox/Source')

from plotlib import DataView
from particles import ParticleSimulation
from numpy import sum

DATA_DIR = 'Dropbox/UMaine/2017_EcologicalComplexity/Data/'
volume = 500.0*500.0*5.0
area = 500.*500.
nlayers = 26

params = ['protein', 'carbohydrate', 'microcystin']
SIM_DIRS = ['100/', '101/', '102/', '103/']
toffsets = [0, 30, 60, 90 ]

# Setup figure for plotting objects
view = DataView(fwidth=3.25,fheight=1.5, verb=True)
for simulation, offset in zip(SIM_DIRS, toffsets):
    data = ParticleSimulation('cyanobacteria', \
        DATA_DIR+simulation, params, verb=True)
    protein = data._ParticleSimulation__state[0]
    carbs = data._ParticleSimulation__state[1]
    toxin = data._ParticleSimulation__state[2]
    toxicity = 1000*sum(toxin, axis=0)/sum(carbs+protein, axis=0)
    data.show(view, toxicity, toffset=offset)

view.fshow(ylab='Toxicity (ppt)', xloc=30, yloc=20, geo=False)




# def load_dissolved(sim_id):
#     dissolved = zeros((nlayers,len(time2)))
#     data = np.loadtxt('./'+sim_id+'/dissolved_toxin.dat', unpack=True)
#     for ii in range(0, nlayers):
#         dissolved[ii,:] = data[ii+1,:]/area*1000.
#     return np.mean(dissolved, axis=0)
# ax.append(ax[0].twinx())
# sur = load_dissolved('102')
# ax[1].plot(time2+60, sur, color='purple')
# sur = load_dissolved('103')
# ax[1].plot(time2+90, sur, color='purple')
# ax[1].set_ylabel('ExoMC (µg/L)', rotation=90)
# ax[1].yaxis.set_major_locator(MultipleLocator(20))
# 
# ax[1].yaxis.label.set_color('purple')
# ax[1].tick_params(axis='y', colors='purple')