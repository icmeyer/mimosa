import numpy as np
import cmath

from region import what_region
from tools import intersection

def make_segments(ray, surfaces, regions, cutoff_length = 300):
    """Make segments using a ray with an initial location and direction

    Parameters
    ----------
    ray : Ray object
        An initialized ray object
    surfaces : list of Surface objects
        List of the surfaces used in the problem
    regions : list of Region objects
        List of regions used in the problem
    cutoff_length : float, opt
        How far to track the ray, default is 200 cm

    Returns
    -------
    ray : Ray object
        A ray object with populated segments attribute
    """
    while ray.length < cutoff_length:
        # Find current region by iterating through surfaces
        region = what_region(ray.r, regions)
        r1 = []
        d_to_beat = cutoff_length

        # Find nearest surface to intersect with
        for surface in surfaces:
            intersect_out = intersection(ray.r, ray.u, surface)
            for loc_dist_pair in intersect_out:
                r1 = loc_dist_pair[0]
                d = loc_dist_pair[1]
                if d > 0 and (np.imag(r1) == 0).all():
                    if d < d_to_beat:
                        d_to_beat = d
                        r = r1
                        boundary = surface.boundary_type
                        col_surface = surface
        # print(ray.r, r, region)
        ray.segments.append(Segment(ray.r, r, region))
        ray.length += np.linalg.norm(r - r1)
        ray.r = r
        d_to_beat = cutoff_length

        # Change r and u in ray and find new surface
        if boundary == 'transmission':
            ray.u = ray.u
        elif boundary == 'reflection' and col_surface.type == 'x-plane':
            ray.u = np.array([-ray.u[0], ray.u[1]])
        elif boundary == 'reflection' and col_surface.type == 'y-plane':
            ray.u = np.array([ray.u[0], -ray.u[1]])
        else:
            raise TypeError('Unkown boundary condition for surface ' + 
                            col_surface.name + 'with ' + boundary)
        # Move ray forward a small bit to ensure location in new region
        smudge = 1e-8
        ray.r += ray.u*smudge

    return ray


class Ray():
    """An object with a position and direction that stores all segements of its
    transport through the geometry

    Parameters
    ----------
    r : 2-tuple
        Position coordinates
    theta : float
        Azimuthal angle
    varphi : float
        Polar angle

    Attributes
    ----------
    r : 2-tuple
        Position coordinates
    u : 2-tuple
        unit vector
    varphi : float
        Polar angle
    region : in
        id of the current region
    segments : list of Segment
        all of the segments of the given ray
    length : float
        length of the total ray
    """
    def __init__(self, r, theta, varphi):
        self.r = r
        self.u = np.array([np.cos(theta), np.sin(theta)])
        self.varphi = varphi
        self.region = None
        self.segments = []
        self.length = 0

class Segment():
    """One segment of a ray

    Parameters
    ----------
    r0 : 2-tuple
        Initial position of segment
    r1 : 2-tuple
        Final position of segment
    region_name : str
        Name of region that the segment is in 
    #TODO: add psi0 and psi1

    Attributes
    ----------
    r0 : 2-tuple
        Initial position of segment
    r1 : 2-tuple
        Final position of segment
    region_name : str
        Name of region that the segment is in 
    #TODO: add psi0 and psi1
    """
    def __init__(self, r0, r1, region_name):
        self.r0 = r0
        self.r1 = r1
        self.region = region_name


    


        
