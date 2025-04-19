#!/usr/bin/python
"""Plot residence time of fish in a given area"""
from matplotlib.pyplot import figure, savefig, LinearLocator, FormatStrFormatter
from numpy import zeros, loadtxt, sqrt
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi

nplots = 1
fontsize = 10
linewidth = 1
hpadding = 0.1
vpadding = 0.05
marginWidth = 8.5
fheight = 6
default_alpha = 1.0


# load data
ini = open("../100/fish_ini.dat", "r")
nfish = int(str.strip(ini.readline()))
time = loadtxt("../100/fish_state.dat", usecols=[0], unpack=True)
start = 0
end = len(time) - 1
dwidth = 4

data = zeros((nfish * dwidth + 1, end - start + 1))
mass = zeros((nfish, end - start + 1))
tox = zeros((nfish, end - start + 1))
path = zeros((nfish, end - start + 1))
xpos = zeros((nfish, end - start + 1))
ypos = zeros((nfish, end - start + 1))
residence = zeros((nfish, end - start + 1))
ratio = zeros((nfish, end - start + 1))
ratio_avg = zeros(end - start + 1)
ratio_avg1 = zeros(end - start + 1)
ratio_avg2 = zeros(end - start + 1)
ratio_std = zeros(end - start + 1)
ratio_sum = zeros(end - start + 1)
suit_sum = zeros(end - start + 1)
time = time / 24.0
nsamples = 4
pindex = [0, 99, 100, 199]
epad = 0.1
lindex = 0
threshold = 20.0
ii = 199
duration = 1


# figure and subplots
fig = figure(facecolor=bg_color, figsize=(marginWidth, fheight))  # Change this
fig.subplots_adjust(
    top=1.0 - vpadding, bottom=vpadding + 0.05, left=hpadding, right=1.0 - (hpadding)
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


print("Reading State")  # experiment A
event_start = 19.5
data = loadtxt("../100/fish_position.dat", unpack=True)
xcol = ii * dwidth + 2
ycol = ii * dwidth + 3
xpos[ii, :] = data[xcol, :]
ypos[ii, :] = data[ycol, :]
print("Calculating Mean Tox Load and Pathway Partitioning")


for jj in range(start, end):
    for kk in range(jj + 1, end):
        distance = sqrt(
            (ypos[ii, jj] - ypos[ii, kk]) ** 2 + (xpos[ii, jj] - xpos[ii, kk]) ** 2
        )
        if distance > threshold:
            next_crossing = kk - 1
            break
        if kk == end:
            next_crossing = kk
    residence[ii, jj] = time[next_crossing] - time[jj]
    for kk in range(jj - 1, start, -1):
        distance = sqrt(
            (ypos[ii, jj] - ypos[ii, kk]) ** 2 + (xpos[ii, jj] - xpos[ii, kk]) ** 2
        )
        if distance > threshold:
            next_crossing = kk + 1
            break
        if kk == start:
            next_crossing = kk
    residence[ii, jj] = residence[ii, jj] + (time[jj] - time[next_crossing])
ax.plot(
    time[240 * event_start : 240 * (event_start + duration)] - event_start,
    residence[ii, 240 * event_start : 240 * (event_start + duration)],
    color="white",
    aa=True,
    linewidth=1,
)
#######################################
print("Reading State")  # experiment B
event_start = 28
data = loadtxt("../101/fish_position.dat", unpack=True)
xcol = ii * dwidth + 2
ycol = ii * dwidth + 3
xpos[ii, :] = data[xcol, :]
ypos[ii, :] = data[ycol, :]
print("Calculating Mean Tox Load and Pathway Partitioning")
for jj in range(start, end):
    for kk in range(jj + 1, end):
        distance = sqrt(
            (ypos[ii, jj] - ypos[ii, kk]) ** 2 + (xpos[ii, jj] - xpos[ii, kk]) ** 2
        )
        if distance > threshold:
            next_crossing = kk - 1
            break
        if kk == end:
            next_crossing = kk
    residence[ii, jj] = time[next_crossing] - time[jj]
    for kk in range(jj - 1, start, -1):
        distance = sqrt(
            (ypos[ii, jj] - ypos[ii, kk]) ** 2 + (xpos[ii, jj] - xpos[ii, kk]) ** 2
        )
        if distance > threshold:
            next_crossing = kk + 1
            break
        if kk == start:
            next_crossing = kk
    residence[ii, jj] = residence[ii, jj] + (time[jj] - time[next_crossing])
ax.plot(
    time[240 * event_start : 240 * (event_start + duration)] - event_start,
    residence[ii, 240 * event_start : 240 * (event_start + duration)],
    color="red",
)

print("Reading State")  # experiment C
event_start = 3
data = loadtxt("../102/fish_position.dat", unpack=True)
xcol = ii * dwidth + 2
ycol = ii * dwidth + 3
xpos[ii, :] = data[xcol, :]
ypos[ii, :] = data[ycol, :]
print("Calculating Mean Tox Load and Pathway Partitioning")

for jj in range(start, end):
    for kk in range(jj + 1, end):
        distance = sqrt(
            (ypos[ii, jj] - ypos[ii, kk]) ** 2 + (xpos[ii, jj] - xpos[ii, kk]) ** 2
        )
        if distance > threshold:
            next_crossing = kk - 1
            break
        if kk == end:
            next_crossing = kk
    residence[ii, jj] = time[next_crossing] - time[jj]
    for kk in range(jj - 1, start, -1):
        distance = sqrt(
            (ypos[ii, jj] - ypos[ii, kk]) ** 2 + (xpos[ii, jj] - xpos[ii, kk]) ** 2
        )
        if distance > threshold:
            next_crossing = kk + 1
            break
        if kk == start:
            next_crossing = kk
    residence[ii, jj] = residence[ii, jj] + (time[jj] - time[next_crossing])

ax.plot(
    time[240 * event_start : 240 * (event_start + duration)] - event_start,
    residence[ii, 240 * event_start : 240 * (event_start + duration)],
    color="green",
)


print("Reading State")  # experiment D
event_start = 10
data = loadtxt("../103/fish_position.dat", unpack=True)
xcol = ii * dwidth + 2
ycol = ii * dwidth + 3
xpos[ii, :] = data[xcol, :]
ypos[ii, :] = data[ycol, :]
print("Calculating Mean Tox Load and Pathway Partitioning")


for jj in range(start, end):
    for kk in range(jj + 1, end):
        distance = sqrt(
            (ypos[ii, jj] - ypos[ii, kk]) ** 2 + (xpos[ii, jj] - xpos[ii, kk]) ** 2
        )
        if distance > threshold:
            next_crossing = kk - 1
            break
        if kk == end:
            next_crossing = kk
    residence[ii, jj] = time[next_crossing] - time[jj]
    for kk in range(jj - 1, start, -1):
        distance = sqrt(
            (ypos[ii, jj] - ypos[ii, kk]) ** 2 + (xpos[ii, jj] - xpos[ii, kk]) ** 2
        )
        if distance > threshold:
            next_crossing = kk + 1
            break
        if kk == start:
            next_crossing = kk
    residence[ii, jj] = residence[ii, jj] + (time[jj] - time[next_crossing])
ax.plot(
    time[240 * event_start : 240 * (event_start + duration)] - event_start,
    residence[ii, 240 * event_start : 240 * (event_start + duration)],
    color="blue",
)


ax.set_xlabel(r"Days")
ax.xaxis.set_major_locator(LinearLocator(5))
ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
ax.set_xlim(0, None)
ax.xaxis.grid(False)
ax.set_ylabel("Residence time, days (with 20 meter threshold)", rotation=90)
ax.yaxis.set_major_locator(LinearLocator(5))
ax.yaxis.grid(False)
ax.set_frame_on(True)
ax.set_ylim()


savefig(
    "./fish_singlesuit_screen.tiff",
    facecolor=bg_color,
    edgecolor="none",
    dpi=default_dpi,
)
