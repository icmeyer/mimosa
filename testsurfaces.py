import numpy as np 

from surface import *
from tools import normalize, intersection
from plotting import * 

circle_test = Circle(surface_id=1, boundary_type='transmission', x0=0.5,
                     y0=0.5, R=0.2)
plane_test =  XPlane(surface_id=2, boundary_type='transmission', x0=1)
test_r = np.asarray([0,0])
test_u = normalize(np.array([1,1]))
circle_out = intersection(test_r, test_u, circle_test)
print(circle_out)
plan_out = intersection(test_r, test_u, plane_test)

