#!/usr/bin/python
from defaults import *
from pylab import *
import matplotlib.animation as animation
from matplotlib.collections import LineCollection

# load data
ini = open('../fish_ini.dat', 'r')
nparticles = int(str.strip(ini.readline()))
time = np.loadtxt('../fish_position.dat', usecols=[0], unpack=True)
start = 0
end = len(time) - 1
datawidth = 4
span = end-start+1
fpsx = zeros( (nparticles, span) ) # particle position
fpsy = zeros( (nparticles, span) ) # particle position
threshold = 7.0
res_time = zeros((nparticles,span))
mres = zeros( (nparticles, 2))
min_alpha = 0.1
max_alpha = 0.4
fade_window = 40
alpha_fade = linspace(0.0, 1.0, fade_window)



for ii in range(0, nparticles):
    xcol = ii*datawidth + 2
    ycol = ii*datawidth + 3
    zcol = ii*datawidth + 4
    fpsx[ii,:], fpsy[ii,:] = np.loadtxt('../fish_position.dat', usecols=(xcol, ycol), unpack=True)

# Calculate Domain Data
vert_x, vert_y = np.loadtxt('../mesh_node.dat', skiprows=1, usecols=(1,2), unpack=True)
xmax = max(vert_x[:])
ymax = max(vert_y[:])
xmin = min(vert_x[:])
ymin = min(vert_y[:])


# calculate residence time
for ii in range(0, nparticles):
    for jj in range(0, end-1):
        for kk in range(jj+1, end):
            dx = fpsx[ii, kk] - fpsx[ii, jj]
            dy = fpsy[ii, kk] - fpsy[ii, jj]
            distance = sqrt(dx**2.0+dy**2.0)
            if distance > threshold:
                res_time[ii,jj] = time[kk-1] - time[jj]
                break
                
    mres[ii,0] = min(res_time[ii,:])
    mres[ii,1] = max(res_time[ii,:])
avg_res = np.mean(np.mean(res_time[:,:]))
min_res = min(mres[:,0])
res_range = max(mres[:,1]) - min_res
        
# plot path and end markers
for ii in range (0, nparticles):
    
    points = np.array([fpsx[ii,0:fade_window+1],fpsy[ii,0:fade_window+1]]).transpose().reshape(-1,1,2)
    segs = np.concatenate([points[:-1],points[1:]],axis=1) # set up a list of segments
    trajectory = LineCollection(segs, zorder=9)
    if ii == 0:
        ptrj = [trajectory]
    else:
        ptrj.append(trajectory)
    ax.add_collection(ptrj[ii])

xyz = zeros((nparticles,2))
scat = plt.scatter(fpsx[:,:], fpsy[:,:], s=25, color='gray', zorder=10, edgecolors='face') # end markers
clock = plt.text(0.0, 0.0,"{:5.2f}".format(time[start+ii])+' hours',color=label_color, fontsize=fontsize*(2./3.), zorder=10) # start markers


clock.set_position( (0.01*(xmax-xmin)+xmin, 0.01*(ymax-ymin)+ymin ) )


ax.set_xlabel(r'easting (m)')
ax.set_xlim(xmin,xmax)
ax.xaxis.set_major_locator(mticker.MultipleLocator(xmax))
ax.xaxis.set_minor_locator(mticker.MultipleLocator(xmax))
ax.xaxis.grid(False)

ax.set_ylabel(r'northing (m)')
ax.set_ylim(ymin,ymax)
ax.yaxis.set_major_locator(mticker.MultipleLocator(ymax))
ax.yaxis.set_minor_locator(mticker.MultipleLocator(ymax))
ax.yaxis.grid(False)

set_tick_sizes(ax, 0, 0)

title(r'fish particle trajectories', fontsize=fontsize, color=label_color)
ttl = ax.title
ttl.set_position([0.5, 1.05])

def animate(ii):
    if (ii < fade_window):
        t_window = ii
    else: 
        t_window = fade_window
        
    # set marker position to current particle position
    xyz[:,1] = fpsy[:, ii+t_window]
    xyz[:,0] = fpsx[:, ii+t_window]
    scat.set_offsets(xyz)
    clock.set_text("{:8.2f}".format(time[start+ii])+' hours, '+"{:5.1f}".format(float(ii*100)/float(end-start))+'%')
    
    for jj in range(0, nparticles):
        color_seq=zeros((fade_window, 4), float)
        color_seq[:,0] = sqrt((res_time[jj, ii:ii+fade_window] - min_res) / res_range) # normalizes
        color_seq[:,2] = ones(fade_window) - color_seq[:,0]
        color_seq[:,3] = 0.75*exp(sqrt(alpha_fade) - 1.0)
        ptrj[jj].set_color(color_seq)
        points = np.array([fpsx[jj, ii:ii+t_window+1],fpsy[jj, ii:ii+t_window+1]]).transpose().reshape(-1,1,2)
        ptrj[jj].set_segments(np.concatenate([points[:-1],points[1:]],axis=1))
        
    
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, frames=(end-fade_window), interval=25, blit=False, repeat=True)

plt.show()
