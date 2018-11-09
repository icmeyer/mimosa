import os
import sys
mimosa_dir = os.path.dirname(os.getcwd())
sys.path.insert(0, mimosa_dir)

import numpy as np
import copy
from main import main 

from surface import XPlane, YPlane, Circle
from tools import normalize, intersection
from plotting import *
from region import Region
from ray import Ray, make_segments

pitch = 1.26
circle = Circle(surface_id=1, boundary_type='transmission', x0=pitch/2,
             y0=pitch/2, R=0.39218) #0.39218
left = XPlane(surface_id=2, boundary_type='reflection', x0=0)
right = XPlane(surface_id=3, boundary_type='reflection', x0=pitch)
top = YPlane(surface_id=4, boundary_type='reflection', y0=pitch)
bottom = YPlane(surface_id=5, boundary_type='reflection', y0=0)

surfaces = [circle, left, right, top, bottom]

ngroup = 10
fuelphiguess = np.ones([ngroup,])
modphiguess = fuelphiguess*0.1
moderator = Region([left, right, top, bottom, circle],[1, -1, -1, 1, 1],
                    uid=0, mat='mod', phi=fuelphiguess)
fuel = Region([circle], [-1], uid=1, mat='fuel', phi=fuelphiguess)
regions = [moderator, fuel]

n_rays = 100
# k, regions = main(n_rays, surfaces, regions, pitch, ngroup, plot=False)

# Ray Length sensitivity analysis
ks_ray_length =  []
lengths = []
for i in range(100):
    try:
        cutoff = 300*i/100+1
        k, regions_trash = main(n_rays, surfaces, copy.deepcopy(regions), pitch, ngroup, cutoff_length=cutoff, deadzone=0)
        lengths.append(cutoff)
        ks_ray_length.append(k)
    except:
        print('oh well')
data = np.vstack((lengths, ks_ray_length)).T
np.savetxt('./sensitivity/cutoff_length_sensitivity_data', data)

ks_dead_zone =  []
lengths = []
# Dead zone sensitivity analysis
for i in range(50):
    try:
        dz = 50*i/50
        k, regions_trash = main(n_rays, surfaces, copy.deepcopy(regions), pitch, ngroup, cutoff_length=50, deadzone=dz)
        lengths.append(dz)
        ks_dead_zone.append(k)
    except:
        print('oh well')
data = np.vstack((lengths, ks_dead_zone)).T
np.savetxt('./sensitivity/dead_zone_sensitivity_data', data)

ks_nrays =  []
nrays_list = []
# Dead zone sensitivity analysis
for i in range(1000):
    try:
        n_rays = i
        k, regions_trash = main(n_rays, surfaces, copy.deepcopy(regions), pitch, ngroup, cutoff_length=300, deadzone=50)
        nrays_list.append(n_rays)
        ks_nrays.append(k)
    except:
        print('oh well')
data = np.vstack((nrays_list, ks_nrays)).T
np.savetxt('./sensitivity/nrays_sensitivity_data', data)
