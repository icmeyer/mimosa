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

length = 1.26
height = 0.1
limits = [0, length, 0, height]
left = XPlane(surface_id=2, boundary_type='reflection', x0=0)
right = XPlane(surface_id=3, boundary_type='reflection', x0=length)
top = YPlane(surface_id=4, boundary_type='reflection', y0=height)
bottom = YPlane(surface_id=5, boundary_type='reflection', y0=0)

surfaces = [left, right, top, bottom]

ngroup = 10
fuelphiguess = np.ones([ngroup,])
slab = Region([left, right, top, bottom],[1, -1, -1, 1],
               uid=0, mat='fuel', phi=fuelphiguess)
regions = [slab]

n_rays = 10
# k, regions = main(n_rays, surfaces, regions, pitch, ngroup, plot=False)

lengths = []

cutoff = 100
k, a_k, regions_trash = main(n_rays, surfaces, copy.deepcopy(regions), limits, ngroup,plot=True, cutoff_length=cutoff, deadzone=10)
# k, a_k, regions_trash = main(n_rays, surfaces, copy.deepcopy(regions), pitch, ngroup, cutoff_length=cutoff, deadzone=10)
lengths.append(cutoff)
 
