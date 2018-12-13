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

ngroup = 10
pitch = 1.26
length = 1.26*3
limits = [0, length, 0, length]
# Surfaces starting with 100 will be circles
# Surfaces starting with 200 will be planes
surfaces = []
super_surfaces = []
circ_counter = 100
xcoords = [2*1.26+1.26/2, 1*1.26+1.26/2, 0*1.26+1.26/2]
ycoords = xcoords

fuel_refine_level = 5
# Biggest radius first
radii = np.linspace(0,0.39218,fuel_refine_level+1)[:0:-1]

for i in range(3):
    for j in range(3):
        sub_surfaces = []
        for level in range(fuel_refine_level):
            circ_counter += 1
            circle = Circle(surface_id=circ_counter, boundary_type='transmission',
                            x0=xcoords[i], y0=ycoords[j], R=radii[level]) #0.39218
            surfaces.append(circle)
            if level == 0:
                super_surface_ext = circle
            sub_surfaces.append(circle)
        super_surface = SuperSurface(circ_counter, surfaces, super_surface_ext,
                                     boundary_type='transmission')
        super_surfaces.append(super_surface)


left = XPlane(surface_id=201, boundary_type='reflection', x0=0)
xplane1 = XPlane(surface_id=202, boundary_type='transmission', x0=pitch)
xplane2 = XPlane(surface_id=203, boundary_type='transmission', x0=2*pitch)
right = XPlane(surface_id=204, boundary_type='reflection', x0=length)
top = YPlane(surface_id=205, boundary_type='reflection', y0=length)
yplane1 = YPlane(surface_id=206, boundary_type='transmission', y0=2*pitch)
yplane2 = YPlane(surface_id=207, boundary_type='transmission', y0=pitch)
bottom = YPlane(surface_id=208, boundary_type='reflection', y0=0)


surfaces += [left, xplane1, xplane2, right, top, yplane1, yplane2, bottom]
super_surfaces += [left, xplane1, xplane2, right, top, yplane1, yplane2, bottom]
# Regions: current limitation, uid of region must be its order in the list
regions = []
super_regions = []
region_counter = 0
fuel_phi_guess = np.ones([ngroup,])
for i in range(9):
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

mod_phi_guess = np.ones([ngroup,])*0.1
mod1 = Region([left, xplane1, top, yplane1, surfaces[fuel_refine_level*0]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
mod2 = Region([xplane1, xplane2, top, yplane1, surfaces[fuel_refine_level*1]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
mod3 = Region([xplane2, right, top, yplane1, surfaces[fuel_refine_level*2]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
mod4 = Region([left, xplane1, yplane1, yplane2, surfaces[fuel_refine_level*3]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
mod5 = Region([xplane1, xplane2, yplane1, yplane2, surfaces[fuel_refine_level*4]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
mod6 = Region([xplane2, right, yplane1, yplane2, surfaces[fuel_refine_level*5]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
mod7 = Region([left, xplane1, yplane2, bottom, surfaces[fuel_refine_level*6]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
mod8 = Region([xplane1, xplane2, yplane2, bottom, surfaces[fuel_refine_level*7]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
mod9 = Region([xplane2, right, yplane2, bottom, surfaces[fuel_refine_level*8]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
regions += [mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9]
super_regions += [mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9]

n_rays = 300
main(n_rays, surfaces, regions, limits, ngroup, super_regions = super_regions,
     super_surfaces = [], plot=True)
# main(n_rays, surfaces, regions, limits, ngroup, super_regions = [],
#      super_surfaces = [], plot=True)

