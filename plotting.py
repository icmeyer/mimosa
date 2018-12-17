import re
import numpy as np
from matplotlib import collections as mc
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib

from tools import get_trailing_numbers

# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def plotlines(lines='', circles = '', length = 1.26, ax=[]):
    if ax == []: # Need a new axis for plot
        fig, ax = plt.subplots()
    else:
        fig = []

    if lines:
        lc = mc.LineCollection(lines[0], colors=lines[1],  linewidths=1)
        ax.add_collection(lc)
        # lc = mc.LineCollection(lines[0], lines[1])
    if circles:
        for circle in circles:
            pltcircle = plt.Circle(circle[0], circle[1], color=circle[2],
                                   fill = False)
            ax.add_artist(pltcircle)
    # ax.autoscale()
    ax.set_xlim([-0.1,length+0.1])
    ax.set_ylim([-0.1,length+0.1])
    ax.margins(0.1)
    ax.set_aspect(1.0)
    return fig, ax

def plot_from_rays(rays, regions, MATERIALS, length=1.26):
    fuelcolors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                  'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    modcolors = ['xkcd:aqua','xkcd:azure','xkcd:blue','xkcd:darkblue',
                 'xkcd:lightblue','xkcd:navy']
    linesegs = []
    linecols = []
    for ray in rays:
        for segment in ray.segments:
            linesegs.append([segment.r0,segment.r1])
            mat = regions[segment.region].mat
            if mat == 'mod':
                linecols.append(modcolors[segment.region%6])
            elif mat == 'fuel':
                linecols.append(fuelcolors[segment.region%10])
    lines = [linesegs, linecols]
    fig, ax = plotlines(lines=lines,length=length)
    fig.show()

def plot_k(iterations, ks, title):
    plt.scatter(iterations, ks)
    plt.xlabel('Iteration')
    plt.ylabel('k')
    plt.title(title)
    plt.savefig('ks.png')
    plt.show()

def plot_flux(e_groups, regions, adjoint=False):
    flux_fig, ax = plt.subplots()
    legend_names = []
    e_groups_width = []
    for i in range(len(e_groups)-1):
        e_groups_width.append(e_groups[i+1] - e_groups[i])

    for region in regions:
        # group_vals = region.phi/e_groups_width
        # flux = np.insert(group_vals[::-1],0,0)
        flux = np.insert(region.phi[::-1],0,0)
        plt.step(e_groups, flux)
        legend_names.append((region.mat + ' Region '+str(region.uid)))
    flux_fig.legend(legend_names)
    plt.xlabel('Energy (eV)')
    plt.ylabel('Forward Flux cm$^{-2}$')
    # ax.set_yscale('log')
    ax.set_xscale('log')

    if adjoint:
        a_fig, ax = plt.subplots()
        legend_names = []
        for region in regions:
            flux = np.insert(region.a_phi[::-1],0,0)
            # flux = np.insert(region.a_phi[::-1]/e_groups_width[::-1],0,0)
            plt.step(e_groups, flux)
            legend_names.append((region.mat + ' Region '+str(region.uid)))
        a_fig.legend(legend_names)
        plt.xlabel('Energy (eV)')
        plt.ylabel('Adjoint Flux cm$^{-2}$')
        # ax.set_yscale('log')
        ax.set_xscale('log')

def plot_flux_on_geometry(ngroup, regions, rays, length, e_group, adjoint=False):
    fluxes = []
    for region in regions:
        if adjoint:
            fluxes.append(region.a_phi[e_group])
        else:
            fluxes.append(region.phi[e_group])
        # fluxes.append(region.phi[ngroup-1])
    max_flux = np.amax(fluxes)
    cmap = matplotlib.cm.YlGnBu_r
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=max_flux)

    linesegs = []
    linecols = []
    for ray in rays:
        for segment in ray.segments:
            linesegs.append([segment.r0,segment.r1])
            if adjoint:
                phi = regions[segment.region].a_phi[e_group]
            else:
                phi = regions[segment.region].phi[e_group]
            linecols.append(cmap(norm(phi)))
    lines = [linesegs, linecols]
    fog_fig, ax = plotlines(lines=lines,length=length)

    fog_fig.subplots_adjust(right = 0.8)
    cbar_ax = fog_fig.add_axes([0.85, 0.15, 0.06, 0.7])
    cbl = matplotlib.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm)
    if adjoint:
        cbl.set_label('Adjoint Flux')
    else:
        cbl.set_label('Flux')
    # fig.colorbar(ims[0], cax=cbar_ax)
    fog_fig.show()

def plot_all_flux_on_geometry(ngroup, regions, rays, length, adjoint=False):
    if ngroup == 10:
        x_plots = 3
        y_plots = 4
        tot_plots = int(x_plots * y_plots)
        main_fig, main_axis = plt.subplots(x_plots,y_plots) 
        main_fig.tight_layout()
    elif ngroup == 2:
        x_plots = 1
        y_plots = 2
        tot_plots = int(x_plots * y_plots)
        main_fig, main_axis = plt.subplots(x_plots,y_plots) 
        main_fig.tight_layout()
    elif ngroup == 1:
        x_plots = 1
        y_plots = 1
        tot_plots = int(x_plots * y_plots)
        main_fig, main_axis = plt.subplots(x_plots,y_plots) 
        main_fig.tight_layout()

    for g in range(tot_plots):
        main_axis = np.reshape(main_axis, [tot_plots, 1])
        axis = main_axis[g][0]
        if g> ngroup-1:
            axis.axis('off')
        else:
            axis.set_title('Group '+str(g+1))
            fluxes = []
            for region in regions:
                if adjoint:
                    fluxes.append(region.a_phi[g])
                else:
                    fluxes.append(region.phi[g])
                # fluxes.append(region.phi[ngroup-1])
            max_flux = np.amax(fluxes)
            cmap = matplotlib.cm.YlGnBu_r
            norm = matplotlib.colors.Normalize(vmin=0.0, vmax=max_flux)

            linesegs = []
            linecols = []
            for ray in rays:
                for segment in ray.segments:
                    linesegs.append([segment.r0,segment.r1])
                    if adjoint:
                        phi = regions[segment.region].a_phi[g]
                    else:
                        phi = regions[segment.region].phi[g]
                    linecols.append(cmap(norm(phi)))
            lines = [linesegs, linecols]
            plotlines(lines=lines, length=length, ax=axis)

    main_fig.subplots_adjust(right = 0.8)
    cbar_ax = main_fig.add_axes([0.85, 0.15, 0.06, 0.7])
    cbl = matplotlib.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm)
    if adjoint:
        cbl.set_label('Adjoint Flux')
    else:
        cbl.set_label('Forward Flux')
    # fig.colorbar(ims[0], cax=cbar_ax)
    main_fig.show()

def pert_plot(pert_dict, scatter_pert_ks, nuf_pert_ks, k0):
    scat_predicted_ks = 1/(1/k0 + pert_dict['scatter'])
    scat_fig, ax = plt.subplots()
    plt.scatter(pert_dict['percent']*100, [k0]*3)
    plt.scatter(pert_dict['percent']*100, scat_predicted_ks)
    plt.scatter(pert_dict['percent']*100, scatter_pert_ks)
    scat_fig.legend(['Original', 'Theory', 'True Perturbation'])
    plt.xlabel('Percent Increase')
    plt.ylabel('$k_{inf}$')
    plt.title('$\Sigma_s$')

    nuf_predicted_ks = 1/(1/k0 + pert_dict['nuf'])
    nuf_fig, ax = plt.subplots()
    plt.scatter(pert_dict['percent']*100, [k0]*3)
    plt.scatter(pert_dict['percent']*100, nuf_predicted_ks)
    plt.scatter(pert_dict['percent']*100, nuf_pert_ks)
    nuf_fig.legend(['Original','Theory', 'True Perturbation'])
    plt.xlabel('Percent Increase')
    plt.ylabel('$k_{inf}$')
    plt.title('nu-$\Sigma_f$')

    plt.show()

if __name__ == '__main__':
    segments = [[(0, 1), (1, 1)], [(2, 3), (-3, 2)], [(0, 2), (2, 3)]]
    colors = np.array([(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)])
    lines = [segments, colors]
    circles = [[(0,0),1,'b'],[(1,1),2,'r']]
    plotlines(lines, circles = circles)
    # plotlines(lines, c)
    
