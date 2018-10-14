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
from physics_new import calc_q, ray_contributions, MATERIALS, normalize_phi
np.random.seed(42)

def main(n_rays, surfaces, regions, length, plot=False, physics=False):
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
    for i in range(n_rays):
        # print('Ray '+str(i)+'/'+str(n_rays))
        rstart = np.array([rand(),rand()])*length
        polar = (2*rand()-1)*pi/2
        theta = rand()*2*pi
        ray_init = Ray(r=rstart, theta=theta, varphi=polar, ngroup = 2)
        ray = make_segments(ray_init, surfaces, regions, deadzone=50)
        all_track_length += ray.length
        rays.append(ray)

    for region in regions:
        # Assign volumes
        region.vol = region.tot_track_length/all_track_length

    print('Tracks laid and volume calculated')
    print(' ')

    if physics:
        print('Begin physics')
        counter = 0
        #Initial k and q guess
        k = 1
        print('Calculating initial q')
        # regions, fission_source_old = calc_q(regions, rays[0].ngroup, k)
        normalize_phi(regions, rays[0].ngroup)
        fission_source_old, k = calc_q(regions, rays[0].ngroup, k)
        print('Begin iterations')
        ks = [k]
        while counter < 50:
            counter += 1
            print('Iteration: ', counter, ' k =', k)

            rays = ray_contributions(rays, regions)

            #Update phi and set counter to 0 for next iteration
            for region in regions:
                # Lecture 5 method
                region.phi = region.tracks_phi/all_track_length + region.q_phi
                # Lecture 7 method
                # region.phi = region.iter_phi + region.q/MATERIALS[region.mat]['total]
                # if any(region.phi < 0):
                #     print(region.mat)
                #     print('q', region.q)
                #     print('phi', region.phi)
                #     print('tracks_phi', region.tracks_phi)
                #     print('tracks_phi normalized', region.tracks_phi/all_track_length)
                #     print('q_phi', region.q_phi)
                #     raise ValueError('whattttt')

                region.tracks_phi = np.zeros(region.phi.shape)
                region.q_phi = np.zeros(region.phi.shape)

            fission_source_new, k = calc_q(regions, rays[0].ngroup, k, 
                                           update_k=True, 
                                           old_fission_source = fission_source_old)
            # print(fission_source_new/fission_source_old)
            fission_source_old = fission_source_new

            ks.append(k)
        
        print('k = ', k, ' after ', counter, 'iterations')

    if plot:
        print('Plotting tracks')
        plot_from_rays(rays, regions, MATERIALS, length = 3*1.26)
        plt.plot(np.arange(counter+1), ks)
        plt.xlabel('Iteration')
        plt.ylabel('k')
        plt.show()
