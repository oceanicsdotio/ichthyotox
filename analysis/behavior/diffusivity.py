#!/usr/bin/python
from matplotlib.pyplot import MultipleLocator, FormatStrFormatter, figure, savefig
from numpy import zeros, loadtxt, std, mean, arange
from sys import argv
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi


hpadding = 0.1
vpadding = 0.05
marginWidth = 6.5
fheight = 4.5

if __name__ == "__main__":
    data_dir = argv[1]

    # load data
    with open(f'{data_dir}/100/fish_ini.dat', 'r', encoding="utf8") as ini:
        nfish = int(str.strip(ini.readline()))
    time = loadtxt(f'{data_dir}/100/fish_position.dat', usecols=[0], unpack=True)
    time = time/24.0
    start = 0
    end = len(time) - 1
    dwidth = 4
    data = zeros((nfish*dwidth+1,end-start+1))
    fpsx = zeros((nfish,end-start+1))
    fpsy = zeros((nfish,end-start+1))
    dx = zeros((nfish,end-start+1))
    dy = zeros((nfish,end-start+1))
    disp_x = zeros(nfish)
    disp_y = zeros(nfish)
    displacement = zeros((nfish,end-start+1))
    residence = zeros((nfish,end-start+1))

    # Figure and subplots
    fig = figure(facecolor=bg_color, figsize=(marginWidth, fheight)) #Change this
    fig.subplots_adjust(top=1.0-vpadding, bottom=vpadding+0.1, left=hpadding, right=1.0-hpadding/3., hspace=0.20)
    ax = []
    ax.append(fig.add_subplot(1,1,1))

    ax[0].patch.set_facecolor(bg_color)
    ax[0].spines['top'].set_color(overlay_color)
    ax[0].spines['bottom'].set_color(overlay_color)
    ax[0].spines['left'].set_color(overlay_color)
    ax[0].spines['right'].set_color(overlay_color)
    ax[0].xaxis.label.set_color(label_color)
    ax[0].yaxis.label.set_color(label_color)
    ax[0].tick_params(axis='x', colors=label_color)
    ax[0].tick_params(axis='y', colors=label_color)


    # experiment phases (days)
    phaselength = 10
    nphases = 2*30//phaselength - 1
    phase_start = range(0, nphases*phaselength//2, phaselength//2)
    phase_end = range(phaselength, nphases*phaselength//2+phaselength, phaselength//2)
    print (phase_start)
    print (phase_end)
    window = 40 # steps = 4 hours
    width = 0.25

    msd_std_a = zeros((nphases,2))
    msd_value_a = zeros((nphases,2))
    msd_std_b = zeros((nphases,2))
    msd_value_b = zeros((nphases,2))
    msd_std_c = zeros((nphases,2))
    msd_value_c = zeros((nphases,2))
    msd_std_d = zeros((nphases,2))
    msd_value_d = zeros((nphases,2))

    # load experiment data
    data = loadtxt(f'{data_dir}/100/fish_position.dat', unpack=True) # Experiment A
    for ii in range(0, nfish):
        xcol = ii*dwidth + 2
        ycol = ii*dwidth + 3
        fpsx[ii,:] = data[xcol,:]
        fpsy[ii,:] = data[ycol,:]
    # calculate individual displacement over one hour time window
    dx[:,1:end] = fpsx[:,1:end] - fpsx[:,0:end-1]
    dy[:,1:end] = fpsy[:,1:end] - fpsy[:,0:end-1]
    dx = dx - 500.0*(dx > 200)
    dx = dx + 500.0*(dx < -200)
    dy = dy - 500.0*(dy > 200)
    dy = dy + 500.0*(dy < -200)
    displacement[:,:] = 0.0
    for tt in range(start+window, end): # loop of values to calculate
        disp_x[:] = 0.0
        disp_y[:] = 0.0
        for uu in range(tt-window, tt): # loop of steps in calc (40)
            disp_x[:] = disp_x[:] + dx[:,uu]
            disp_y[:] = disp_y[:] + dy[:,uu]
        displacement[:,tt] = displacement[:,tt] + disp_x**2 + disp_y**2 # displacement at end of one hour
        displacement[:,tt] = 0.25*displacement[:,tt] # average square displacement per hour
    # calculate summary stats
    for ii in range(0, nphases):
        msd_value_a[ii,0] = mean(mean(displacement[:,phase_start[ii]*240:phase_end[ii]*240], axis=0))
        msd_std_a[ii,0] = std(mean(displacement[:,phase_start[ii]*240:phase_end[ii]*240], axis=0))
    #rects_a = ax[0].bar(arange(nphases), msd_value_a[:,0], width, color='black', edgecolor='none', yerr=msd_std_a[:,0], ecolor='black')
    rects_a = ax[0].bar(arange(nphases), msd_value_a[:,0], width, color='black', edgecolor='none', label='Control (A)')
    print ("A:", msd_value_a[:,0])
    print ("A:", msd_std_a[:,0])



    ############################# load experiment data #######################################
    data = loadtxt(f'{data_dir}/101/fish_position.dat', unpack=True) # experiment B
    for ii in range(0, nfish):
        xcol = ii*dwidth + 2
        ycol = ii*dwidth + 3
        fpsx[ii,:] = data[xcol,:]
        fpsy[ii,:] = data[ycol,:]
    # calculate individual displacement over one hour time window
    dx[:,1:end] = fpsx[:,1:end] - fpsx[:,0:end-1]
    dy[:,1:end] = fpsy[:,1:end] - fpsy[:,0:end-1]
    dx = dx - 500.0*(dx > 200)
    dx = dx + 500.0*(dx < -200)
    dy = dy - 500.0*(dy > 200)
    dy = dy + 500.0*(dy < -200)
    displacement[:,:] = 0.0
    for tt in range(start+window, end): # loop of values to calculate
        disp_x[:] = 0.0
        disp_y[:] = 0.0
        for uu in range(tt-window, tt): # loop of steps in calc (40)
            disp_x[:] = disp_x[:] + dx[:,uu]
            disp_y[:] = disp_y[:] + dy[:,uu]
        displacement[:,tt] = displacement[:,tt] + disp_x**2 + disp_y**2 # displacement at end of one hour
        displacement[:,tt] = 0.25*displacement[:,tt] # average square displacement per hour
    # calculate summary stats
    for ii in range(0, nphases):
        msd_std_b[ii,0] = std(mean(displacement[:,phase_start[ii]*240:phase_end[ii]*240], axis=0))
        msd_value_b[ii,0] = mean(mean(displacement[:,phase_start[ii]*240:phase_end[ii]*240], axis=0))
    #rects_b = ax[0].bar(arange(nphases)+width, msd_value_b[:,0], width, color='red', edgecolor='none', yerr=msd_std_b[:,0], ecolor='red')
    rects_b = ax[0].bar(arange(nphases)+width, msd_value_b[:,0], width, color='red', edgecolor='none', label='Formation (B)')
    print ("B:", msd_value_b[:,0])
    print ("B:", msd_std_b[:,0])

    # load experiment data
    data = loadtxt(f'{data_dir}/102/fish_position.dat', unpack=True) # experiment C
    for ii in range(0, nfish):
        xcol = ii*dwidth + 2
        ycol = ii*dwidth + 3
        fpsx[ii,:] = data[xcol,:]
        fpsy[ii,:] = data[ycol,:]
    # calculate individual displacement over one hour time window
    dx[:,1:end] = fpsx[:,1:end] - fpsx[:,0:end-1]
    dy[:,1:end] = fpsy[:,1:end] - fpsy[:,0:end-1]
    dx = dx - 500.0*(dx > 200)
    dx = dx + 500.0*(dx < -200)
    dy = dy - 500.0*(dy > 200)
    dy = dy + 500.0*(dy < -200)
    displacement[:,:] = 0.0
    for tt in range(start+window, end): # loop of values to calculate
        disp_x[:] = 0.0
        disp_y[:] = 0.0
        for uu in range(tt-window, tt): # loop of steps in calc (10)
            disp_x[:] = disp_x[:] + dx[:,uu]
            disp_y[:] = disp_y[:] + dy[:,uu]
        displacement[:,tt] = displacement[:,tt] + disp_x**2 + disp_y**2 # displacement at end of one hour
        displacement[:,tt] = 0.25*displacement[:,tt]# average square displacement per hour
    # calculate summary stats
    for ii in range(0, nphases):
        msd_std_c[ii,0] = std(mean(displacement[0:99,phase_start[ii]*240:phase_end[ii]*240], axis=0))
        msd_value_c[ii,0] = mean(mean(displacement[0:99,phase_start[ii]*240:phase_end[ii]*240], axis=0))
        msd_std_c[ii,1] = std(mean(displacement[100:199,phase_start[ii]*240:phase_end[ii]*240], axis=0))
        msd_value_c[ii,1] = mean(mean(displacement[100:199,phase_start[ii]*240:phase_end[ii]*240], axis=0))
    #rects_c = ax[0].bar(arange(nphases)+2*width, msd_value_c[:,0], width/2, color='green', edgecolor='none', yerr=msd_std_c[:,0], ecolor='green')
    #rects_c = ax[0].bar(arange(nphases)+2.5*width, msd_value_c[:,1], width/2, color='green', edgecolor='none', yerr=msd_std_c[:,1], ecolor='green')

    rects_c = ax[0].bar(arange(nphases)+2*width, msd_value_c[:,0], width/2, color='green', edgecolor='none', label='Intensification, surface (C)')
    rects_c = ax[0].bar(arange(nphases)+2.5*width, msd_value_c[:,1], width/2, color=[0.0,1.0,0.0,0.5], edgecolor='none', label='Intensification, demersal')
    print ("C1:", msd_value_c[:,0])
    print ("C1:", msd_std_c[:,0])
    print ("C2:", msd_value_c[:,1])
    print ("C2:", msd_std_c[:,1])


    # load experiment data
    data = loadtxt(f'{data_dir}/103/fish_position.dat', unpack=True) # experiment D
    for ii in range(0, nfish):
        xcol = ii*dwidth + 2
        ycol = ii*dwidth + 3
        fpsx[ii,:] = data[xcol,:]
        fpsy[ii,:] = data[ycol,:]
    # calculate individual displacement over one hour time window
    dx[:,1:end] = fpsx[:,1:end] - fpsx[:,0:end-1]
    dy[:,1:end] = fpsy[:,1:end] - fpsy[:,0:end-1]
    dx = dx - 500.0*(dx > 200)
    dx = dx + 500.0*(dx < -200)
    dy = dy - 500.0*(dy > 200)
    dy = dy + 500.0*(dy < -200)
    displacement[:,:] = 0.0
    for tt in range(start+window, end): # loop of values to calculate
        disp_x[:] = 0.0
        disp_y[:] = 0.0
        for uu in range(tt-window, tt): # loop of steps in calc (10)
            disp_x[:] = disp_x[:] + dx[:,uu]
            disp_y[:] = disp_y[:] + dy[:,uu]
        displacement[:,tt] = displacement[:,tt] + disp_x**2 + disp_y**2 # displacement at end of one hour
        displacement[:,tt] = 0.25*displacement[:,tt] # average square displacement per hour
    # calculate summary stats
    for ii in range(0, nphases):
        msd_std_d[ii,0] = std(mean(displacement[0:100,phase_start[ii]*240:phase_end[ii]*240], axis=0))
        msd_value_d[ii,0] = mean(mean(displacement[0:100,phase_start[ii]*240:phase_end[ii]*240], axis=0))
        msd_std_d[ii,1] = std(mean(displacement[100:199,phase_start[ii]*240:phase_end[ii]*240],axis=0))
        msd_value_d[ii,1] = mean(mean(displacement[100:199,phase_start[ii]*240:phase_end[ii]*240],axis=0))
    #rects_d = ax[0].bar(arange(nphases)+3*width, msd_value_d[:,0], width/2, color='blue', edgecolor='none', yerr=msd_std_d[:,0], ecolor='blue')
    #rects_d = ax[0].bar(arange(nphases)+3.5*width, msd_value_d[:,1], width/2, color='blue', edgecolor='none', yerr=msd_std_d[:,1], ecolor='blue')

    rects_d = ax[0].bar(arange(nphases)+3*width, msd_value_d[:,0], width/2, color='blue', edgecolor='none', label='Decline, surface (D)')
    rects_d = ax[0].bar(arange(nphases)+3.5*width, msd_value_d[:,1], width/2, color=[0.0,0.0,1.0,0.5], edgecolor='none', label='Decline, demersal')
    print ("D1:", msd_value_d[:,0])
    print ("D1:", msd_std_d[:,0])
    print ("D2:", msd_value_d[:,1])
    print( "D2:", msd_std_d[:,1])

    # plot path and end markers
    ax[0].set_xlabel(r'Phase')
    ax[0].set_xlim()
    #ax[0].xaxis.set_major_locator(MultipleLocator())
    ax[0].xaxis.grid(False)

    ax[0].set_ylabel(r'Diffusivity (m^2/hr)')
    ax[0].set_ylim(0,800)
    ax[0].yaxis.set_major_locator(MultipleLocator(100))
    ax[0].yaxis.grid(False)

    legend = ax[0].legend(loc='upper left')
    legend.get_frame().set_facecolor('none')
    legend.get_frame().set_edgecolor('none')
    text = legend.get_texts()
    text[0].set_color('black')
    text[1].set_color('red')
    text[2].set_color('green')
    text[3].set_color([0.0,1.0,0.0,0.5])
    text[4].set_color('blue')
    text[5].set_color([0.0,0.0,1.0,0.5])


    savefig('./figures/behavior/diffusivity.png', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
