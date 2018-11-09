import numpy as np
import sys
from math import pi
sys.path.insert(0, '/home/icmeyer/22.212/moc') #add directory above to path

from my_extension import intersection_c
from ray import Ray
from surface import XPlane


ray_init = Ray(r=np.array([1,1]), theta= pi, varphi=0)
ray_init = Ray(r=np.array([1,2]), theta= pi, varphi=0)
xplane = XPlane(x0=1, boundary_type = 'reflective', surface_id=1)

out = intersection_c(ray_init.r[0], ray_init.r[1], ray_init.u[0], ray_init.u[1], xplane)


