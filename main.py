# This file will drive the simulation
import re
import numpy as np
import time
import sys
from numpy.random import random_sample as rand
from math import pi

from surface import XPlane, YPlane, Circle, intersection
from tools import *
from plotting import *
from region import Region
from ray import Ray, make_segments
from physics import calc_q, ray_contributions, normalize_phi
from materials import import_xs
from perturbation import *
from header import *
np.random.seed(42)


def main(n_rays, surfaces, regions, limits, ngroup, plot=False,
        cutoff_length=300, deadzone=10, super_regions=[],
        super_surfaces=[], pert=[], k_guess = 0, rays=[]):
    """ Run MOC and write outputs to file 

    Parameters
    ----------
    limits : list
        [xmin, xmax, ymin, ymax]

    """
    # Import materials library
    MATERIALS = import_xs(ngroup, pert=pert)

    # Assign regions and surfaces for ray tracing if no heirarchal
    # data is given 
    if super_regions == []:
        super_regions = regions
    if super_surfaces == []:
        super_surfaces = surfaces


    physics = True
    perturbation = False
    quiet = False

    start = time.perf_counter()
    print(header)
    if rays == []:
        print('Laying down tracks')
        track_time_0 = time.perf_counter()
        all_track_length = 0
        all_active_length = 0

        one_percent = np.max([np.floor(n_rays/100),1])
        # Initiate rays and fill each with segments
        for i in range(n_rays):
            if i%one_percent == 0:
                sys.stdout.write("\r{0:.0%}".format(i/n_rays))
                sys.stdout.flush()
            rstart = np.zeros(2)
            rstart[0] = rand()*(limits[1]-limits[0])+limits[0]
            rstart[1] = rand()*(limits[3]-limits[2])+limits[2]
            polar = (2*rand()-1)*pi/2
            theta = rand()*2*pi
            ray_init = Ray(r=rstart, theta=theta, varphi=polar)
            ray = make_segments(ray_init, surfaces, regions, 
                                cutoff_length=cutoff_length, deadzone=deadzone,
                                super_surfaces = super_surfaces, 
                                super_regions = super_regions)
            all_track_length += ray.length
            all_active_length += ray.active_length
            rays.append(ray)

        for region in regions:
            # Assign volumes
            region.vol = region.tot_track_length/all_active_length
            print('Region vol:', region.mat, region.vol)

        print('Tracks laid and volume calculated')
        track_time_1 = time.perf_counter()
        print('Track laying time: ', track_time_1-track_time_0)

    else: 
        all_active_length = 0
        for ray in rays:
            all_active_length += ray.active_length

    if physics:
        print('Begin physics')
        counter = 0
        #Initial k, a_k, and q guess
        if k_guess == 0:
            k = 1
            a_k = 1
        else:
            k = k_guess
            a_k = k_guess
        print('Calculating initial q')
        fission_source_old, a_fission_source_old, k, a_k = calc_q(regions, ngroup, k, a_k, pert=pert)

        ks = [k]
        a_ks = [a_k]
        converged = False
        phi_conv_flag = False
        print('Begin iterations')
        # while counter < 3:
        # while not converged and counter < 2:
        # while counter < 300:
        while not converged and counter < 300:
            # normalize_phi(regions, ngroup)
            #Print out flux in each region
            # for region in regions:
            #     print(counter, 'Flux in region', region.uid, region.mat, region.phi)
            counter += 1
            if not quiet:
                print('------Iterations: ', counter, ' k = ', k, ' a_k = ', a_k)
            rays = ray_contributions(rays, ngroup, regions, pert=pert)

            #Update phi and set counter to 0 for next iteration
            for region in regions:
                sigma_t = MATERIALS[region.mat]['total']
                vol = region.vol
                term = (1/vol/sigma_t)
                # term = (1/sigma_t)
                region.phi = (term*region.tracks_phi/all_active_length
                              + 4*pi*region.q)
                region.a_phi = (term*region.tracks_a_phi/all_active_length
                              + 4*pi*region.a_q)

                # Zero out phi counters
                region.tracks_phi = np.zeros(region.phi.shape)
                region.tracks_a_phi = np.zeros(region.phi.shape)
                region.q_phi = np.zeros(region.phi.shape)

            fission_source_new, a_fission_source_new, k, a_k = calc_q(regions, ngroup, k, a_k, 
                                                update_k=True, old_fission_source=fission_source_old, old_a_fission_source=a_fission_source_old, pert=pert)

            # print('fission_sources', fission_source_old, a_fission_source_old)

            fission_source_old = fission_source_new
            a_fission_source_old = a_fission_source_new

            ks.append(k)
            a_ks.append(a_k)
            converged = checktol(ks[counter-1], k, tol=1e-6) & \
                        checktol(a_ks[counter-1], a_k, tol=1e-5) & \
                        phi_conv_flag

            if counter==1:
                old_a_phi = collect_phi(regions, adjoint=True)
            else:
                phi_conv_flag, new_a_phi = check_phi_convergence(regions,
                                                                 old_a_phi,
                                                                 adjoint=True)
                old_a_phi = new_a_phi


            # for region in regions:
            #     print(region.mat)
            #     print(region.phi)
            #     print(region.a_phi)
        
        print('k = ', k, ' and a_k = ', a_k, ' after ', counter, 'iterations')
        end = time.perf_counter()
        elapsed_time = end - start
        segments = 0
        for ray in rays:
            for segment in ray.segments:
                segments += 1

        print('Elapsed time:               ', elapsed_time)
        try:
            print('Time per segment per group: ', elapsed_time/(segments*ngroup))
        except:
            print('Time per segment per group: n/a')

    normalize_phi(regions, ngroup, adjoint=True)
    if plot:
        ktitle ='k = '+str(k)+' Rays ='+str(n_rays)
        length = np.amax([limits[1]-limits[0],limits[3]-limits[2]])
        print('Plotting tracks')
        plot_from_rays(rays, regions, MATERIALS, length = length)
        # plot_k(np.arange(counter+1),ks, ktitle)
        print('Plotting forward flux')
        plot_all_flux_on_geometry(ngroup, regions, rays, length)
        print('Plotting adjoint flux')
        plot_all_flux_on_geometry(ngroup, regions, rays, length, adjoint=True)
        e_group = 0
        # plot_flux_on_geometry(ngroup, regions, rays, length, e_group, adjoint=True)
        # e_group = ngroup
        # plot_flux_on_geometry(ngroup, regions, rays, length, e_group, adjoint=True)
        # e_group = 1
        # plot_flux_on_geometry(ngroup, regions, rays, length, e_group, adjoint=False)
        if ngroup == 10:
            energy_groups = [0.0, 0.058, 0.14, 0.28, 0.625, 4.0, 10.0, 40.0, 5530.0, 821e3, 20e6]
            plot_flux(energy_groups, regions, adjoint = True)
        elif ngroup == 20 or ngroup == 40:
            energy_groups = np.geomspace(1e-4,20e6,ngroup+1)
            energy_groups[0] = 0
            plot_flux(energy_groups, regions, adjoint = True)
        plt.show()

    if perturbation:
        pert_dict = calc_perturbation(regions, MATERIALS, k, a_k)

            
    # Zero out ray lengths and volumes so region objects can be reused
    # for region in regions:
    #     region.tot_track_length = 0
    #     region.vol = 0

    return k, a_k, regions, rays

# Helpful snippet below for checking for negative values

# if any(region.phi < 0):
#     print(region.mat)
#     print('q', region.q)
#     print('phi', region.phi)
#     print('tracks_phi', region.tracks_phi)
#     print('tracks_phi normalized', region.tracks_phi/all_track_length)
#     print('q_phi', region.q_phi)
#     raise ValueError('whattttt')
