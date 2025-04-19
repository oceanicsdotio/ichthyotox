#!/usr/bin/python
from matplotlib.pyplot import MultipleLocator, FormatStrFormatter, figure, savefig
from numpy import zeros, loadtxt, mean
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi
from analysis.behavior import State

nplots=1
fontsize  = 10
linewidth = 1
hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = 3.0
default_alpha = 1.0

# load data
ini = open('../100/fish_ini.dat', 'r', encoding="utf8")
nfish = int(str.strip(ini.readline()))
time = loadtxt('../100/fish_state.dat', usecols=[0], unpack=True)
start = 0
end = len(time) - 1

time = time / 24.

fig = figure(facecolor=bg_color, figsize=(marginWidth, fheight))
fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.1, left=hpadding, right=1.0-hpadding)
ax = fig.add_subplot(1,1,1)

ax.patch.set_facecolor(bg_color)
ax.spines['top'].set_color(overlay_color)
ax.spines['bottom'].set_color(overlay_color)
ax.spines['left'].set_color(overlay_color)
ax.spines['right'].set_color(overlay_color)
ax.xaxis.label.set_color(label_color)
ax.yaxis.label.set_color(label_color)
ax.tick_params(axis='x', colors=label_color)
ax.tick_params(axis='y', colors=label_color)

mass_avg = State('../100/').mean_mass()
ax.plot(time, mass_avg, linestyle='-', linewidth=1, color='white', aa=True, label='Control (A)')

mass_avg = State('../101/').mean_mass()
ax.plot(time, mass_avg, linestyle='-', linewidth=1, color='red', aa=True, label='Formation (B)')

mass_avg = State('../300/').mean_mass()
ax.plot(time, mass_avg, linestyle='-', linewidth=1, color='orange', aa=True, label='Low Suit.')

mass_avg = State('../301/').mean_mass()
ax.plot(time, mass_avg, linestyle='-', linewidth=1, color='yellow', aa=True, label='Low Suit./Low Dep.')

mass_avg = State('../302/').mean_mass()
ax.plot(time, mass_avg, linestyle='-', linewidth=1, color='purple', aa=True, label='Low Suit./Extra Impaired')

mass_avg = State('../303/').mean_mass()
ax.plot(time, mass_avg, linestyle='-', linewidth=1, color='green', aa=True, label='Low Suit./Low Dep/Extra Impaired')

mass_avg = State('../304/').mean_mass()
ax.plot(time, mass_avg, linestyle='-', linewidth=1, color='blue', aa=True, label='Extra Impaired')

mass_avg = State('../305/').mean_mass()
ax.plot(time, mass_avg, linestyle='-', linewidth=1, color='grey', aa=True, label='Low Dep.')

mass_avg = State('../306/').mean_mass()
ax.plot(time, mass_avg, linestyle='--', linewidth=1, color='white', aa=True, label='Low Dep./Extra Impaired')


ax.set_xlabel(r'Days')
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.set_xlim(0,30)
ax.xaxis.grid(False)

ax.set_ylabel('Mass (g)', rotation=90)

ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.yaxis.grid(False)
ax.set_frame_on(True)
ax.set_ylim(11,19)

legend = ax.legend(loc='upper left')
legend.get_frame().set_facecolor('none')
legend.get_frame().set_edgecolor('none')
text = legend.get_texts()
text[0].set_color('white')
text[1].set_color('red')
text[2].set_color('orange')
text[3].set_color('yellow')
text[4].set_color('purple')
text[5].set_color('green')
text[6].set_color('blue')
text[7].set_color('grey')
text[8].set_color('white')


savefig('./fish_growth_screen2.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
