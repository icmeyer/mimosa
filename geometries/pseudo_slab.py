import os
import sys
mimosa_dir = os.path.dirname(os.getcwd())
sys.path.insert(0, mimosa_dir)

import numpy as np
import copy
from main import main 

from surface import XPlane, YPlane, Circle
from tools import normalize
from plotting import *
from region import Region, what_region
from ray import Ray, make_segments

length = 1.26
height = 0.2
limits = [0, length, 0, height]
left = XPlane(surface_id=0, boundary_type='vacuum', x0=0)
right = XPlane(surface_id=1, boundary_type='vacuum', x0=length)
# left = XPlane(surface_id=0, boundary_type= 'reflection', x0=0)
# right = XPlane(surface_id=1, boundary_type='reflection', x0=length)
top = YPlane(surface_id=2, boundary_type='reflection', y0=height)
bottom = YPlane(surface_id=3, boundary_type='reflection', y0=0)

surfaces = [left, right, top, bottom]
regions = []

ngroup = 1
refine_slab = 30
delta = length/refine_slab
fuelguess = np.ones([ngroup,])
for i in range(refine_slab):
    surf = XPlane(surface_id=4+i, boundary_type = 'transmission', x0 = (i+1)*delta)
    surfaces.append(surf)
    if i == 0 :
        print(left.x0)
        print(surf.x0)
        slab_bit = Region([left, surf, top, bottom], [1, -1, -1, 1], uid=i,
                          mat='fuel', phi=fuelguess)
    elif i == refine_slab-1:
        surf1 = surfaces[4+i - 1]
        print(surf1.x0)
        print(right.x0)
        slab_bit = Region([surf1, right, top, bottom], [1, -1, -1, 1], uid=i,
                          mat='fuel', phi=fuelguess)
    else:
        surf1 = surfaces[4+i - 1]
        surf2 = surf
        print(surf1.x0)
        print(surf2.x0)
        slab_bit = Region([surf1, surf2, top, bottom], [1, -1, -1, 1], uid=i,
                          mat='fuel', phi=fuelguess)
    regions.append(slab_bit)

for region in regions:
    print(region)
    for surface in region.surfaces:
        print(surface)
        print(surface.boundary_type)


# grid = 10
# xs = np.linspace(0,length,grid)
# ys = np.linspace(0,height,grid)
# for i in range(grid):
#     for j in range(grid):
#         print([xs[i],ys[j]])
#         print(what_region([xs[i],ys[j]], regions))

n_rays = 30
# k, regions = main(n_rays, surfaces, regions, pitch, ngroup, plot=False)

lengths = []

cutoff = 200
k, a_k, regions_trash = main(n_rays, surfaces, copy.deepcopy(regions), limits, ngroup,plot=True, cutoff_length=cutoff, deadzone=10)
# k, a_k, regions_trash = main(n_rays, surfaces, copy.deepcopy(regions), pitch, ngroup, cutoff_length=cutoff, deadzone=10)
lengths.append(cutoff)
 
