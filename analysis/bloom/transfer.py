#!/usr/bin/python
from matplotlib.pyplot import figure, savefig
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from pylab import LinearLocator, FormatStrFormatter, linspace
from analysis.defaults import bg_color, overlay_color, label_color, default_dpi

if __name__ == "__main__":
    fontsize  = 10
    linewidth = 1
    uniformPadding = 0.05
    hpad=0.05
    vpad=0.05
    marginWidth = 6.5
    default_alpha = 1.0
    lineRGBA = [1.0, 0.0, 0.0, default_alpha]
    style = ['-',':']
    div=150

    fig = figure(facecolor=bg_color, figsize=(marginWidth, 2.*marginWidth/3.)) #Change this
    fig.subplots_adjust(top=1.0-uniformPadding, bottom=uniformPadding, left=uniformPadding, right=1.0-uniformPadding/2.)
    
    # surface subplots
    ax: list[Axes3D] = []
    ax.append(fig.add_subplot(2,2,1,projection='3d'))
    ax.append(fig.add_subplot(2,2,2,projection='3d'))
    ax.append(fig.add_subplot(2,2,3,projection='3d'))
    ax.append(fig.add_subplot(2,2,4,projection='3d'))

    ax[0].patch.set_facecolor(bg_color)
    ax[0].spines['top'].set_color(overlay_color)
    ax[0].spines['bottom'].set_color(overlay_color)
    ax[0].spines['left'].set_color(overlay_color)
    ax[0].spines['right'].set_color(overlay_color)
    ax[0].xaxis.label.set_color(label_color)
    ax[0].yaxis.label.set_color(label_color)
    ax[0].tick_params(axis='x', colors=label_color)
    ax[0].tick_params(axis='y', colors=label_color)

    protein = 1.0
    irradiance = linspace(0.,650.,div)
    temperature = linspace(20.0,30.0,div)
    tp1 = linspace(0.0,1.0,div)
    tp2 = linspace(0.0,1.0,div)
    ratio = linspace(0.0,4.0,div)
    protein = linspace(100.0,100.0,div)
    carbs = protein*ratio
    maxxes = np.zeros(4)

    xx = linspace(1.0,1.0,div)

    #FIXATION
    irrad_r = irradiance/250.
    light_limit = 1.98 * irrad_r / (irrad_r**2. - 0.02*irrad_r + 1.)
    llgr, carbsgr = np.meshgrid(light_limit, carbs, sparse=True)
    irradiance, rgr = np.meshgrid(irradiance, ratio, sparse=True)
    fixation = 0.02622 * (4.*protein - carbsgr) * llgr
    maxxes[0] = np.max(abs(fixation[:,:]))
    ax[0].plot_surface(irradiance, np.log(rgr), abs(fixation), color=bg_color, edgecolor=overlay_color, linewidth=1, antialiased=True, shade=False)

    #SYNTHESIS
    tp1 = (temperature/25. * (-0.1*(temperature - 35.))**(3./25.))**4. # error here
    for ii in range(0,div):
        if (temperature[ii] > 35.0):
            tp1[ii]= 0.0
    tp1gr, carbsgr = np.meshgrid(tp1, carbs, sparse=True)
    temperature, rgr = np.meshgrid(temperature, ratio, sparse=True)
    synthesis = 0.05 * carbsgr * tp1gr
    maxxes[1] = np.max(abs(synthesis[:,:]))
    ax[1].plot_surface(temperature, np.log(rgr), abs(synthesis), color=bg_color, edgecolor=overlay_color, linewidth=1, antialiased=True, shade=False)


    tmp_planes = ax[0].zaxis._PLANES 
    ax[0].zaxis._PLANES = ( tmp_planes[2], tmp_planes[3], tmp_planes[0], tmp_planes[1], tmp_planes[4], tmp_planes[5] )


    # RESPIRATION
    tp2 = 0.286*np.exp(.05*temperature - 0.15)
    tp2gr, progr = np.meshgrid(tp2, protein, sparse=True)
    respiration = 0.01*carbsgr*tp1gr - 0.004*protein*tp2gr
    maxxes[2] = np.max(fixation[:,:])
    ax[2].plot_surface(temperature, np.log(rgr), abs(respiration), color=bg_color, edgecolor=overlay_color, linewidth=1, antialiased=True, shade=False)

    # EXCRETION
    excretion = .1*(0.004*carbsgr + 0.05*progr)*tp2gr
    ax[3].plot_surface(temperature, np.log(rgr), abs(excretion), color=bg_color, edgecolor=overlay_color, linewidth=1, antialiased=True, shade=False)

    sublabels = [['I', r'ln(C/M)', r'| fixation |'], [r'T', r'ln(C/M)', r'| synthesis |'], [r'T', r'ln(C/M)', r'| respiration| '], [r'T', r'ln(C/M)', r'| excretion |']]
    subtitles = [r'(A)', r'(B)', r'(C)', r'(D)']
    for ii in range(0,4):
        # axis adjustment
        ax[ii].zaxis._PLANES = ( tmp_planes[2], tmp_planes[3], tmp_planes[0], tmp_planes[1], tmp_planes[4], tmp_planes[5] )

        ax[ii].view_init(elev = 30., azim=135.0)
        if (ii > 0):
            ax[ii].view_init(elev = 30., azim=135.0+180.0)
        
        
        ax[ii].set_xlabel(sublabels[ii][0]); ax[ii].set_xlim()
        ax[ii].xaxis.set_major_locator(LinearLocator(2))
        ax[ii].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        ax[ii].xaxis.grid(False)
        ax[ii].xaxis.set_rotate_label(False)
        ax[ii].xaxis.set_pane_color(bg_color)

        ax[ii].set_ylabel(sublabels[ii][1]); ax[ii].set_ylim()
        ax[ii].yaxis.set_major_locator(LinearLocator(2))
        ax[ii].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        ax[ii].yaxis.grid(False)
        ax[ii].yaxis.set_rotate_label(False)
        ax[ii].yaxis.set_pane_color(bg_color)

        ax[ii].set_zlabel(sublabels[ii][2])
        ax[ii].autoscale(enable=True, axis='z', tight=True)
        #ax[ii].set_zlim(0.0,)
        ax[ii].zaxis.set_major_locator(LinearLocator(2))
        ax[ii].zaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax[ii].zaxis.grid(False)
        ax[ii].zaxis.label.set_rotation(360.0)
        ax[ii].zaxis.set_pane_color(bg_color)
        
        
        ax[ii].set_title(subtitles[ii], fontsize=fontsize, color=label_color)
        ax[ii].title.set_position([1.0, 1.0])
        
        # ax[ii].set_frame_on(False)

    ax[0].set_xlim(0.0,650.0)
    ax[1].set_xlim(20.,30.); ax[1].set_ylim()
    ax[2].set_xlim(20.,30.); ax[2].set_ylim()
    ax[3].set_xlim(20.,30.); ax[3].set_ylim()


    fig.tight_layout()
    savefig('./figures/bloom/carbon/transfer.png', facecolor=bg_color, edgecolor='none', dpi=default_dpi)
    