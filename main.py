NGROUP = 10
header = """
  _  _  __  _  _   __   ____   __   
 ( \/ )(  )( \/ ) /  \ / ___) / _\  
 / \/ \ )( / \/ \(  O )\___ \/    \ 
 \_)(_/(__)\_)(_/ \__/ (____/\_/\_/ 
 
 """
#MOC Is Magic, Onerous, Silly, and Awesome
                                    
# This file will drive the simulation
import re
import numpy as np
import time
from numpy.random import random_sample as rand
from math import pi

from surface import XPlane, YPlane, Circle
from tools import normalize, intersection, checktol
from plotting import *
from region import Region
from ray import Ray, make_segments
from physics import calc_q, ray_contributions, normalize_phi
from materials import MATERIALS
np.random.seed(42)

def main(n_rays, surfaces, regions, limits, ngroup, plot=False, physics=True,
        cutoff_length=300, deadzone=50):
    """ Run MOC and write outputs to file 

    Parameters
    ----------
    limits : list
        [xmin, xmax, ymin, ymax]

    """
    start = time.perf_counter()
    print(header)
    rays = []
    print('Laying down tracks')
    all_track_length = 0
    all_active_length = 0
    # Initiate rays and fill each with segments
    for i in range(n_rays):
        rstart = np.zeros(2)
        rstart[0] = rand()*(limits[1]-limits[0])+limits[0]
        rstart[1] = rand()*(limits[3]-limits[2])+limits[2]
        polar = (2*rand()-1)*pi/2
        theta = rand()*2*pi
        ray_init = Ray(r=rstart, theta=theta, varphi=polar)
        ray = make_segments(ray_init, surfaces, regions, cutoff_length=cutoff_length, deadzone=deadzone)
        all_track_length += ray.length
        all_active_length += ray.active_length
        rays.append(ray)

    for region in regions:
        # Assign volumes
        region.vol = region.tot_track_length/all_active_length
        print('Region vol:', region.mat, region.vol)

    print('Tracks laid and volume calculated')

    if physics:
        print('Begin physics')
        counter = 0
        #Initial k, a_k, and q guess
        k = 1
        a_k = 1
        print('Calculating initial q')
        fission_source_old, k, a_k = calc_q(regions, ngroup, k, a_k)

        ks = [k]
        a_ks = [a_k]
        converged = False
        print('Begin iterations')
        # while counter < 3:
        while not converged and counter < 500:
            normalize_phi(regions, ngroup)
            #Print out flux in each region
            # for region in regions:
            #     print(counter, 'Flux in region', region.uid, region.mat, region.phi)
            counter += 1
            print('Iterations: ', counter, ' k = ', k)
            rays = ray_contributions(rays, ngroup, regions)

            #Update phi and set counter to 0 for next iteration
            for region in regions:
                sigma_t = MATERIALS[region.mat]['total']
                vol = region.vol
                term = (1/vol/sigma_t)
                region.phi = (term*region.tracks_phi/all_active_length
                              + 4*pi*region.q)
                region.a_phi = (term*region.tracks_a_phi/all_active_length
                              + 4*pi*region.a_q)

                # Zero out phi counters
                region.tracks_phi = np.zeros(region.phi.shape)
                region.tracks_a_phi = np.zeros(region.phi.shape)
                region.q_phi = np.zeros(region.phi.shape)

            fission_source_new, k, a_k = calc_q(regions, ngroup, k, a_k, 
                                                update_k=True, old_fission_source=fission_source_old)
            fission_source_old = fission_source_new

            ks.append(k)
            a_ks.append(a_k)
            converged = checktol(ks[counter-1], k, tol=1e-5) &  \
                        checktol(a_ks[counter-1], a_k, tol=1e-5)

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

    if plot:
        ktitle ='k = '+str(k)+' Rays ='+str(n_rays)
        print('Plotting tracks')
        length = np.amax([limits[1]-limits[0],limits[3]-limits[2]])
        plot_from_rays(rays, regions, MATERIALS, length = length)
        plot_k(np.arange(counter+1),ks, ktitle)
        if ngroup == 10:
            energy_groups = [0.0, 0.058, 0.14, 0.28, 0.625, 4.0, 10.0, 40.0, 5530.0, 821e3, 20e6]
            plot_flux(energy_groups, regions, adjoint = True)

    return k, a_k, regions

# Helpful snippet below for checking for negative values

# if any(region.phi < 0):
#     print(region.mat)
#     print('q', region.q)
#     print('phi', region.phi)
#     print('tracks_phi', region.tracks_phi)
#     print('tracks_phi normalized', region.tracks_phi/all_track_length)
#     print('q_phi', region.q_phi)
#     raise ValueError('whattttt')
