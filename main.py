# This file will drive the simulation
import re
import numpy as np
from numpy.random import random_sample as rand
from math import pi

from surface import XPlane, YPlane, Circle
from tools import normalize, intersection, get_trailing_numbers
from plotting import *
from region import Region
from ray import Ray, make_segments
np.random.seed(42)

def main(n_rays, surfaces, regions, pitch, plot=False):
    """ Run MOC and write outputs file 

    Parameters
    ----------
    n_rays : int
        number of rays to simulate
    deadzone : int
        Amount of deadzone to use in simulation

    Returns
    -------
    
    """
    print('Start of simulation!')

    rays = []
    for i in range(n_rays):
        print('Ray '+str(i)+'/'+str(n_rays))
        rstart = np.array([rand(),rand()])*pitch
        polar = (2*rand()-1)*pi/2
        theta = rand()*2*pi
        ray_init = Ray(r=rstart, theta=theta, varphi=polar)
        ray = make_segments(ray_init, surfaces, regions)
        rays.append(ray)

    if plot:
        plot_from_rays(rays)
