import os
import sys
mimosa_dir = os.path.dirname(os.getcwd())
sys.path.insert(0, mimosa_dir)

import numpy as np
from main import main 

from surface import SuperSurface, XPlane, YPlane, Circle
from tools import normalize
from plotting import *
from region import SuperRegion, Region, what_region
from ray import Ray, make_segments

lattice = 2
fuel_refine_level = 4
ngroup = 10
pitch = 1.26
length = 1.26*lattice
limits = [0, length, 0, length]
# Surfaces starting with 100 will be circles
# Surfaces starting with 200 will be planes
surfaces = []
super_surfaces = []
circ_counter = 100
xcoords = []
for i in range(lattice):
    xcoords.append(i*1.26+1.26/2)
ycoords = xcoords[::-1]

# Biggest radius first
radii = np.linspace(0,0.39218,fuel_refine_level+1)[:0:-1]

for i in range(lattice):
    for j in range(lattice):
        sub_surfaces = []
        for level in range(fuel_refine_level):
            circ_counter += 1
            circle = Circle(surface_id=circ_counter, boundary_type='transmission',
                            x0=xcoords[j], y0=ycoords[i], R=radii[level]) #0.39218
            surfaces.append(circle)
            if level == 0:
                super_surface_ext = circle
            sub_surfaces.append(circle)
        super_surface = SuperSurface(circ_counter, surfaces, super_surface_ext,
                                     boundary_type='transmission')
        super_surfaces.append(super_surface)

# Regions: current limitation, uid of region must be its order in the list
regions = []
super_regions = []
region_counter = 0
fuel_phi_guess = np.ones([ngroup,])
for i in range(lattice*lattice):
    sub_regions = []
    outer = fuel_refine_level*i
    inner = outer+fuel_refine_level-1
    # Make donuts
    for level in range(fuel_refine_level - 1):
        fuel = Region([surfaces[outer+level], surfaces[outer+level+1]], [-1, 1], uid=region_counter, mat='fuel',phi=fuel_phi_guess)
        region_counter += 1
        regions.append(fuel)
        sub_regions.append(fuel)
    # Make donut hole
    fuel = Region([surfaces[inner]], [-1], uid=region_counter, mat='fuel',phi=fuel_phi_guess)
    region_counter += 1
    regions.append(fuel)
    sub_regions.append(fuel)
    super_region = SuperRegion(sub_regions, [surfaces[outer]], [-1], uid=i)
    super_regions.append(super_region)

planes_only = []
# top is planes_only[0] bottomw is planes_only[9]
# left is planes_only[10] bottomw is planes_only[19]
plane_count = 0
for i in range(lattice+1):
    if i==0 or i==lattice:
        bt = 'reflection'
    else:
        bt = 'transmission'
    y = 1.26*lattice-1.26*i
    planes_only.append(YPlane(surface_id=plane_count, boundary_type=bt, y0=y))
    plane_count += 1
for i in range(lattice+1):
    if i==0 or i==lattice:
        bt = 'reflection'
    else:
        bt = 'transmission'
    x = 1.26*i
    planes_only.append(XPlane(surface_id=plane_count, boundary_type=bt, x0=x))
    plane_count += 1
surfaces += planes_only

mod_regions = []
mod_phi_guess = np.ones([ngroup,])*0.1
for i in range(lattice):
    for j in range(lattice):
        top = planes_only[i]
        bottom = planes_only[i+1]
        left = planes_only[lattice+1+j]
        right = planes_only[lattice+1+j+1]
        circle_num = fuel_refine_level*(i*lattice+j)
        mod = Region([left, right, top, bottom, surfaces[circle_num]],
                     [1, -1, -1, 1, 1], uid=region_counter, mat='mod',
                     phi = mod_phi_guess)
        # print('Position', circle_num)
        # print('x',left.x0, right.x0)
        # print('y',bottom.y0, top.y0)
        print(surfaces[circle_num].x0, surfaces[circle_num].y0, surfaces[circle_num].r)
        mod_regions.append(mod)
        region_counter += 1

regions += mod_regions
super_regions += mod_regions


# n = 40
# ys = np.linspace(0,2.52,n)
# for i in range(n):
#     r = (1.89,ys[i])
#     print(what_region(r,regions))

n_rays = 20
main(n_rays, surfaces, regions, limits, ngroup, super_regions = super_regions,
      super_surfaces = [], plot=True)
# main(n_rays, surfaces, regions, limits, ngroup, super_regions = [],
#      super_surfaces = [], plot=True)

