#!/usr/bin/python
from matplotlib.pyplot import hist, title, savefig, figure, setp, LinearLocator, FormatStrFormatter
from numpy import mean, var, linspace, loadtxt
from scipy.stats import norm
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi

fontsize  = 10
linewidth = 1
uniformPadding = 0.2
vpad = 0.2
hpad = 0.1
marginWidth = 7.0
default_alpha = 1.0
lineRGBA = [1.0, 0.0, 0.0, default_alpha]
style = ['-',':']

bg_color = [1.0,1.0,1.0,default_alpha]
overlay_color = [0.0,0.0,0.0,default_alpha]
label_color = [0.0,0.0,0.0,default_alpha]
default_dpi = 300

fig = figure(facecolor=bg_color, figsize=(marginWidth, marginWidth/3.)) #Change this
fig.subplots_adjust(top=1.0-vpad/2., bottom=vpad, left=hpad, right=1.0-hpad/2.)
   
# set color scheme
ax = fig.add_subplot(1,3,1)

ax.patch.set_facecolor(bg_color)
ax.spines['top'].set_color(overlay_color)
ax.spines['bottom'].set_color(overlay_color)
ax.spines['left'].set_color(overlay_color)
ax.spines['right'].set_color(overlay_color)
ax.xaxis.label.set_color(label_color)
ax.yaxis.label.set_color(label_color)
ax.tick_params(axis='x', colors=label_color)
ax.tick_params(axis='y', colors=label_color)


# variable declarations
gaussian = loadtxt('../random.dat')



print( "N=50: ", mean(gaussian[0:49]))
print( "N=500: ", mean(gaussian[0:499]))
print( "N=5000: ", mean(gaussian[0:4999]))

print( "N=50: ", var(gaussian[0:49]))
print( "N=500: ", var(gaussian[0:499]))
print( "N=5000: ", var(gaussian[0:4999]))



# surface and contours
# first create a single histogram

# the histogram of the data with histtype='step'
n, bins, patches = hist(gaussian[0:49], 20, normed=1, histtype='stepfilled')
setp(patches, 'facecolor', 'white')

# add a line showing the expected distribution
y = norm.pdf(linspace(-4.0, 4.0, 100), 0.0, 1.0)
l = ax.plot(linspace(-4.0, 4.0, 100), y, 'k--', linewidth=1.0)


# axis adjustment
ax.set_xlabel(r'X')
ax.set_xlim(-4.,4.)
ax.xaxis.set_major_locator(LinearLocator(3))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.xaxis.grid(False)

ax.set_ylabel(r'P(X)')
ax.set_ylim(0.,0.5)
ax.yaxis.set_major_locator(LinearLocator(3))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.grid(False)

title(r'(A)', fontsize=fontsize, color=label_color)
ttl = ax.title
ttl.set_position([1.0, 1.0])





#############
ax2 = fig.add_subplot(1,3,2)

n2, bins2, patches2 = hist(gaussian[0:499], 20, normed=1, histtype='stepfilled')
setp(patches2, 'facecolor', 'white')

# add a line showing the expected distribution
y = norm.pdf(linspace(-4.0, 4.0, 100), 0.0, 1.0)
l = ax2.plot(linspace(-4.0, 4.0, 100), y, 'k--', linewidth=1.0)


# axis adjustment
ax2.set_xlabel(r'X')
ax2.set_xlim(-4.,4.)
ax2.xaxis.set_major_locator(LinearLocator(3))
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax2.xaxis.grid(False)

ax2.set_ylabel(r' ')
ax2.set_ylim(0.,0.5)
ax2.yaxis.set_major_locator(LinearLocator(0))
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax2.yaxis.grid(False)

title(r'(B)', fontsize=fontsize, color=label_color)
ttl = ax2.title
ttl.set_position([1.0, 1.0])

ax3 = fig.add_subplot(1,3,3)

n3, bins3, patches3 = hist(gaussian[0:4999], 20, normed=1, histtype='stepfilled')
setp(patches3, 'facecolor', 'white')

# add a line showing the expected distribution
y = norm.pdf(linspace(-4.0, 4.0, 100), 0.0, 1.0)
l = ax3.plot(linspace(-4.0, 4.0, 100), y, 'k--', linewidth=1.0)


# axis adjustment
ax3.set_xlabel(r'X')
ax3.set_xlim(-4.,4.)
ax3.xaxis.set_major_locator(LinearLocator(3))
ax3.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax3.xaxis.grid(False)

ax3.set_ylabel(r' ')
ax3.set_ylim(0.,0.5)
ax3.yaxis.set_major_locator(LinearLocator(0))
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax3.yaxis.grid(False)

title(r'(C)', fontsize=fontsize, color=label_color)
ttl = ax3.title
ttl.set_position([1.0, 1.0])

savefig('./random_histogram.tiff', facecolor=bg_color, edgecolor='none', dpi=default_dpi)