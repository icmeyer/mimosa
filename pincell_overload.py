import numpy as np
from main import main 

from surface import XPlane, YPlane, Circle
from tools import normalize, intersection
from plotting import *
from region import Region
from ray import Ray, make_segments

pitch = 1.26
circle1 = Circle(surface_id=1, boundary_type='transmission', x0=pitch/2,
             y0=pitch/2, R=0.39218) #0.39218
circle2 = Circle(surface_id=1, boundary_type='transmission', x0=pitch/2,
             y0=pitch/2, R=0.29218) #0.39218
circle3 = Circle(surface_id=1, boundary_type='transmission', x0=pitch/2,
             y0=pitch/2, R=0.19218) #0.39218
circle4 = Circle(surface_id=1, boundary_type='transmission', x0=0.1,
             y0=0.1, R=0.09) #0.39218
left = XPlane(surface_id=2, boundary_type='reflection', x0=0)
right = XPlane(surface_id=3, boundary_type='reflection', x0=pitch)
top = YPlane(surface_id=4, boundary_type='reflection', y0=pitch)
bottom = YPlane(surface_id=5, boundary_type='reflection', y0=0)


surfaces = [left, right, top, bottom]
circles = []
for i in range(50):
    circle = Circle(surface_id=1, boundary_type='transmission', x0=pitch/2,
             y0=pitch/2, R=0.5*(i/50))
    circles.append(circle)
    surfaces.append(circle)
regions = []
for i in range(49):
    i = i+1
    regions.append(Region([circles[i-1],circles[i]], [1, -1], name='fuel'+str(i)))
regions.append(Region([circles[0]], [-1], name='fuel0'))
moderator = Region([left, right, top, bottom, circles[-1]], [1, -1, -1, 1, 1], name = 'mod') 
regions.append(moderator)

n_rays = 50
main(n_rays, surfaces, regions, pitch, plot=True)

# linesegs = []
# linecols = []
# if False:
#     for segment in rays[0].segments:
#         linesegs.append([segment.r0,segment.r1])
#         if segment.region == 'mod':
#             linecols.append('b')
#         elif segment.region == 'fuel':
#             linecols.append('g')
#         elif segment.region == None:
#             linecols.append('r')
#     lines = [linesegs, linecols]
#     plotlines(lines=lines, circles= '')
