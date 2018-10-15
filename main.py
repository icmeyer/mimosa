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
from numpy.random import random_sample as rand
from math import pi

from surface import XPlane, YPlane, Circle
from tools import normalize, intersection, get_trailing_numbers, checktol
from plotting import *
from region import Region
from ray import Ray, make_segments
from physics import calc_q, ray_contributions, MATERIALS, normalize_phi
np.random.seed(42)

def main(n_rays, surfaces, regions, length, ngroup, plot=False, physics=False):
    """ Run MOC and write outputs to file 

    Parameters
    ----------
    n_rays : int
        number of rays to simulate
    deadzone : int
        Amount of deadzone to use in simulation
    
    """
    print(header)
    rays = []
    print('Laying down tracks')
    all_track_length = 0
    all_active_length = 0
    # Initiate rays and fill each with segments
    for i in range(n_rays):
        rstart = np.array([rand(),rand()])*length
        polar = (2*rand()-1)*pi/2
        theta = rand()*2*pi
        ray_init = Ray(r=rstart, theta=theta, varphi=polar)
        ray = make_segments(ray_init, surfaces, regions, cutoff_length=300)
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
        #Initial k and q guess
        k = 1
        print('Calculating initial q')
        normalize_phi(regions, ngroup)
        fission_source_old, k = calc_q(regions, ngroup, k)

        ks = [k]
        converged = False
        # while counter < 2:
        print('Begin iterations')
        while not converged:
            counter += 1
            print('Iterations: ', counter, ' k = ', k)
            rays = ray_contributions(rays, ngroup, regions)

            #Update phi and set counter to 0 for next iteration
            for region in regions:
                sigma_t = MATERIALS[region.mat]['total']
                vol = region.vol
                term = (1/vol/sigma_t)
                region.phi = (term*region.tracks_phi/all_active_length
                              + region.q_phi)

                # Zero out phi counters
                region.tracks_phi = np.zeros(region.phi.shape)
                region.q_phi = np.zeros(region.phi.shape)

            fission_source_new, k = calc_q(regions, ngroup, k, update_k=True, 
                                           old_fission_source = fission_source_old)
            fission_source_old = fission_source_new

            ks.append(k)
            converged = checktol(ks[counter-1], k, tol=1e-5)
            
            #Print out flux in each region
            for region in regions:
                print('Flux in region', region.uid, ' ', region.phi)
        
        print('k = ', k, ' after ', counter, 'iterations')

    if plot:
        ktitle ='k='+str(k)+' n_rays='+str(n_rays)
        print('Plotting tracks')
        plot_from_rays(rays, regions, MATERIALS, length = 3*1.26)
        plot_k(np.arange(counter+1),ks, ktitle)

# Helpful snippet below for checking for negative values

# if any(region.phi < 0):
#     print(region.mat)
#     print('q', region.q)
#     print('phi', region.phi)
#     print('tracks_phi', region.tracks_phi)
#     print('tracks_phi normalized', region.tracks_phi/all_track_length)
#     print('q_phi', region.q_phi)
#     raise ValueError('whattttt')
