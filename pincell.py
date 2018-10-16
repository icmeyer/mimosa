import numpy as np
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
                    uid=0, mat='mod', phi=modphiguess)
fuel = Region([circle], [-1], uid=1, mat='fuel', phi=fuelphiguess)
regions = [moderator, fuel]

n_rays = 1000
main(n_rays, surfaces, regions, pitch, ngroup, plot=True, physics=True)

