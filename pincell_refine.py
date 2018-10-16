import numpy as np
from main import main 

from surface import XPlane, YPlane, Circle
from tools import normalize, intersection
from plotting import *
from region import Region, what_region
from ray import Ray, make_segments

ngroup = 10
pitch = 1.26
# Surfaces starting with 100 will be circles
# Surfaces starting with 200 will be planes
surfaces = []
circ_counter = 100
xcoord = pitch/2
ycoord = xcoord

fuel_refine_level = 6
# Biggest radius first
radii = np.linspace(0,0.39218,fuel_refine_level+1)[:0:-1]

for level in range(fuel_refine_level):
    circ_counter += 1
    circle = Circle(surface_id=circ_counter, boundary_type='transmission',
                    x0=xcoord, y0=ycoord, R=radii[level]) #0.39218
    surfaces.append(circle)

left = XPlane(surface_id=201, boundary_type='reflection', x0=0)
right = XPlane(surface_id=204, boundary_type='reflection', x0=pitch)
top = YPlane(surface_id=205, boundary_type='reflection', y0=pitch)
bottom = YPlane(surface_id=208, boundary_type='reflection', y0=0)

surfaces += [left, right, top, bottom]
# Regions: current limitation, uid of region must be its order in the list
regions = []
region_counter = 0
fuel_phi_guess = np.ones([ngroup,])

outer = fuel_refine_level*0
inner = outer+fuel_refine_level-1
# Make donuts
for level in range(fuel_refine_level - 1):
    fuel = Region([surfaces[outer+level], surfaces[outer+level+1]], [-1, 1], uid=region_counter, mat='fuel',phi=fuel_phi_guess)
    region_counter += 1
    regions.append(fuel)
# Make donut hole
fuel = Region([surfaces[inner]], [-1], uid=region_counter, mat='fuel',phi=fuel_phi_guess)
region_counter += 1
regions.append(fuel)

mod_phi_guess = np.ones([ngroup,])*0.1
mod1 = Region([left, right, top, bottom, surfaces[fuel_refine_level*0]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
regions += [mod1]

n_rays = 100
main(n_rays, surfaces, regions, pitch, ngroup, plot=True, physics=True)

