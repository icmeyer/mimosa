import numpy as np
from main import main 

from surface import XPlane, YPlane, Circle
from tools import normalize, intersection
from plotting import *
from region import Region
from ray import Ray, make_segments

pitch = 1.26
length = 1.26*3
# Surfaces starting with 100 will be circles
# Surfaces starting with 200 will be planes
surfaces = []
circ_counter = 100
xcoords = [2*1.26+1.26/2, 1*1.26+1.26/2, 0*1.26+1.26/2]
ycoords = xcoords

for i in range(3):
    for j in range(3):
        circ_counter += 1
        circle = Circle(surface_id=circ_counter, boundary_type='transmission',
                        x0=xcoords[i], y0=ycoords[j], R=0.39218) #0.39218
        surfaces.append(circle)

left = XPlane(surface_id=201, boundary_type='reflection', x0=0)
xplane1 = XPlane(surface_id=202, boundary_type='transmission', x0=pitch)
xplane2 = XPlane(surface_id=203, boundary_type='transmission', x0=2*pitch)
right = XPlane(surface_id=204, boundary_type='reflection', x0=length)
top = YPlane(surface_id=205, boundary_type='reflection', y0=length)
yplane1 = YPlane(surface_id=206, boundary_type='transmission', y0=2*pitch)
yplane2 = YPlane(surface_id=207, boundary_type='transmission', y0=pitch)
bottom = YPlane(surface_id=208, boundary_type='reflection', y0=0)

surfaces += [left, xplane1, xplane2, right, top, yplane1, yplane2, bottom]
# Regions
regions = []
region_counter = 0
for i in range(9):
    fuel = Region([surfaces[i]], [-1], uid=region_counter, mat='fuel',phi=[1,0.1])
    region_counter += 1
    regions.append(fuel)
mod1 = Region([left, xplane1, top, yplane1, surfaces[1]],[1, -1, -1, 1, 1],
                    uid=9, mat='mod', phi = [0.1, 0.1])
mod2 = Region([xplane1, xplane2, top, yplane1, surfaces[1]],[1, -1, -1, 1, 1],
                    uid=10, mat='mod', phi = [0.1, 0.1])
mod3 = Region([xplane2, right, top, yplane1, surfaces[1]],[1, -1, -1, 1, 1],
                    uid=11, mat='mod', phi = [0.1, 0.1])
mod4 = Region([left, xplane1, yplane1, yplane2, surfaces[1]],[1, -1, -1, 1, 1],
                    uid=12, mat='mod', phi = [0.1, 0.1])
mod5 = Region([xplane1, xplane2, yplane1, yplane2, surfaces[1]],[1, -1, -1, 1, 1],
                    uid=13, mat='mod', phi = [0.1, 0.1])
mod6 = Region([xplane2, right, yplane1, yplane2, surfaces[1]],[1, -1, -1, 1, 1],
                    uid=14, mat='mod', phi = [0.1, 0.1])
mod7 = Region([left, xplane1, yplane2, bottom, surfaces[1]],[1, -1, -1, 1, 1],
                    uid=15, mat='mod', phi = [0.1, 0.1])
mod8 = Region([xplane1, xplane2, yplane2, bottom, surfaces[1]],[1, -1, -1, 1, 1],
                    uid=16, mat='mod', phi = [0.1, 0.1])
mod9 = Region([xplane2, right, yplane2, bottom, surfaces[1]],[1, -1, -1, 1, 1],
                    uid=17, mat='mod', phi = [0.1, 0.1])
regions += [mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9]

n_rays = 20
main(n_rays, surfaces, regions, length, plot=True, physics=True)

