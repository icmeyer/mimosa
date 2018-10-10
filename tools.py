import numpy as np
import cmath
import re

def normalize(x):
    return x/np.linalg.norm(x)

def intersection(r, u, surface):
    """ Calculates the intersection of a ray with a given surface

    Parameters
    ----------
    r : 2-tuple
        Position coordinates
    u : 2-tuple
        unit vector
    surface : Surface
        Surface object to check for intersection with

    Returns
    -------
    intersect : n-tuple
        All points of intersection on the surface
    """
    if surface.type == 'circle':
        radius = surface.r
        d_to_center = r-[surface.x0, surface.y0]
        quad1 = -(np.dot(u, d_to_center))
        test = -np.linalg.norm(d_to_center)**2
        test2 = radius**2
        test3 = (quad1)**2
        radicand = ((quad1)**2-np.linalg.norm(d_to_center)**2+radius**2)
        quad2 = cmath.sqrt(radicand)
        d1 = quad1+quad2
        d2 = quad1-quad2
        r1 = r + np.dot(d1,u)
        r2 = r + np.dot(d2,u)

        r1 = real_reduce(r1)
        r2 = real_reduce(r2)
        d1 = real_reduce(d1)
        d2 = real_reduce(d2)

        return [[r1, d1],[r2, d2]]

    elif surface.type == 'x-plane':
        p0 = np.array([surface.x0,0])
        normal = np.array([1,0])
        d = real_reduce(np.dot((p0-r),normal)/np.dot(u,normal))
        r1 = real_reduce(r + np.dot(d,u))
        return [[r1, d]]

    elif surface.type == 'y-plane':
        p0 = np.array([0,surface.y0])
        normal = np.array([0,1])
        d = real_reduce(np.dot((p0-r),normal)/np.dot(u,normal))
        r1 = real_reduce(r + np.dot(d,u))
        return [[r1, d]]

def get_trailing_numbers(s, zero=False):
    m = re.search(r'\d+$', s)
    if m:
        return int(m.group())
    elif zero:
        return 0
    else:
        return None

def real_reduce(num):
    """ Reduce complex numbers to reals if imaginary part is 0 """
    return num.real if (num.imag == 0).all() else num
