#!/usr/bin/python
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import rc
from pylab import *

def ugrow(id):
	print("Reading State of "+id)
	data = np.loadtxt('../'+id+'/fish_state.dat', unpack=True)
	for ii in range(0, nfish):
		masscol = ii*dwidth + 2
		mass[ii,:] = data[masscol, :]
	return np.mean(mass[:,end])


# load data
ini = open('../100/fish_ini.dat', 'r')
nfish = int(str.strip(ini.readline()))
time = np.loadtxt('../100/fish_state.dat', usecols=[0], unpack=True)
start = 0
end = len(time) - 1
dwidth = 4

data = zeros((nfish*dwidth+1,end-start+1))
mass = zeros((nfish,end-start+1))
time = time / 24.


print(ugrow('100'))
print(ugrow('101'))
print(ugrow('102'))
print(ugrow('103'))
print(ugrow('300'))
print(ugrow('301'))
print(ugrow('302'))
print(ugrow('303'))
print(ugrow('304'))
print(ugrow('305'))
print(ugrow('306'))
print(ugrow('400'))
print(ugrow('401'))
print(ugrow('402'))
