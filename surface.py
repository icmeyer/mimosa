import numpy as np
import cmath

from tools import real_reduce

class SuperSurface():
    """
    A set of surfaces, will return closest intersection in case
    of a ray passing through outermost surface, otherwise only
    returns the imaginary

    Parameters
    ----------
    surfaces : list
        List of surfaces contained within the set
    ext_surface : list
        Exterior surface, first to be evaluated for intersection
    surface_id : int
        Unique id for surface

    Attributes
    ----------
    surfaces : list
        List of surfaces contained within the set
    type : str
        Type of the surface
    """
    def __init__(self, surface_id, surfaces, ext_surface, 
                 boundary_type):
        self.id = surface_id
        self.surfaces = surfaces
        self.ext_surface = ext_surface
        self.boundary_type = boundary_type
        self.type = 'super_surface'

    def intersect(cls, r, u):
        """
        Parameters
        ----------
        r : tuple

        Returns
        -------
        n-tuple
            The closest real intersection in the set of surfaces
        """
        # First check if there are any intersections with the outer surface
        # If not we won't evaluate any more surfaces
        d_to_beat = np.Inf
        intersect_out = intersection(r, u, cls.ext_surface)
        for loc_dist_pair in intersect_out:
            r1 = loc_dist_pair[0]
            d = loc_dist_pair[1]
            if d > 0 and (np.imag(r1) == 0).all(): #Has a real collision
                # Find nearest surface to intersect with
                for surface in cls.surfaces:
                    intersect_out = intersection(r, u, surface)
                    for loc_dist_pair in intersect_out:
                        r1test = loc_dist_pair[0]
                        dtest = loc_dist_pair[1]
                        if dtest > 0 and (np.imag(r1test) == 0).all():
                            # If new intersection is closer, make this
                            # collision surface
                            if dtest < d_to_beat:
                                d_to_beat = dtest
                                dout = dtest
                                rout = r1test
                return [[rout, d]]
            else: # Not a real collision
                return [[r1, d]]


class Surface(object):
    """
    The surface class is used to define the surfaces of the problem. 
    It will contain methods for evaluating intersections at surfaces.
    Only three surface types are currently nece

    Parameters
    ----------
    surface_id : int
        Unique integer used to define surface
    boundary_type : str
        Type of boundary [transmission, reflection, vacuum]
    name : str, optional
        Name of surface

    Attributes
    ----------
    type : str
        Type of the surface
    """
    def __init__(self, surface_id, boundary_type, name =''):
        self.id = surface_id
        self.boundary_type = boundary_type
        self.name = name

class XPlane(Surface):
    """
    The xplane represents a plane perpindicular to the x surface
    
    Parameters
    ----------
    surface_id : int
        Unique integer used to define surface
    boundary_type : str
        Type of boundary [transmission, reflection, vacuum]
    x0 : float
        x-coordinate of plane
    name : str, optional
        Name of surface

    Attributes
    ----------
    id : int
        Unique integer used to define surface
    boundary_type : str
        Type of boundary [transmission, reflection, vacuum]
    name : str
        Name of surface
    type : str
        Type of the surface
    """
    type = 'x-plane'

    def __init__(self, surface_id, boundary_type, x0, name=''):
        super().__init__(surface_id, boundary_type, name=name)
        self.x0 = x0

    def evaluate(self, point):
        """Evaluate the surface equation at a given point.

        Parameters
        ----------
        point : 3-tuple of float
            (x,y,x)

        Returns
        -------
        float
            x - x_0
        """
        return point[0] - self.x0

class YPlane(Surface):
    """
    The yplane represents a plane perpindicular to the y surface
    
    Parameters
    ----------
    surface_id : int
        Unique integer used to define surface
    boundary_type : str
        Type of boundary [transmission, reflection, vacuum]
    y0 : float
        y-coordinate of plane
    name : str, optional
        Name of surface

    Attributes
    ----------
    id : int
        Unique integer used to define surface
    boundary_type : str
        Type of boundary [transmission, reflection, vacuum]
    name : str
        Name of surface
    type : str
        Type of the surface
    """
    type = 'y-plane'

    def __init__(self, surface_id, boundary_type, y0, name=''):
        super().__init__(surface_id, boundary_type, name=name)
        self.y0 = y0

    def evaluate(self, point):
        """Evaluate the surface equation at a given point.

        Parameters
        ----------
        point : 2-tuple of float
            (x,y)

        Returns
        -------
        float
            y - y_0
        """
        return point[1] - self.y0

class Circle(Surface):
    """ A circle defined by (x-x0)^2 + (y-y0)^2 = R^2

    Parameters
    ----------
    surface_id : int
        Unique integer used to define surface
    boundary_type : str
        Type of boundary [transmission, reflection, vacuum]
    x0 : float
    y0 : float
    R : float
    name : str, optional
        Name of surface

    Attributes
    ----------
    id : int
        Unique integer used to define surface
    boundary_type : str
        Type of boundary [transmission, reflection, vacuum]
    x0 : float
    y0 : float
    r : float
        radius of the circle
    name : str, optional
        Name of surface
    type : str
        type of the surface
    """
    type = 'circle'

    def __init__(self, surface_id, boundary_type, x0, y0, R, name=''):
        super().__init__(surface_id, boundary_type, name=name)
        self.x0 = x0
        self.y0 = y0
        self.r = R

    def evaluate(self, point):
        """Evaluate the surface equation at a given point.

        Parameters
        ----------
        point : 2-tuple of float
            (x,y)

        Returns
        -------
        float
            x - x_0
        """
        x = point[0] - self.x0
        y = point[1] - self.y0
        return x**2 + y**2 - self.r**2


def intersection(r, u, surface):
    """ Calculates the intersection of a ray with a given surface, or
    in the case of a SuperSurface, returns the closest intersecting
    point within that set

    Parameters
    ----------
    r : 2-tuple
        Position coordinates
    u : 2-tuple
        unit vector
    surface : Surface or SuperSurface
        Surface object to check for intersection with

    Returns
    -------
    intersect : n-tuple
        All points of intersection on the surface
    """
    if surface.type == 'circle':
        radius = surface.r
        d_to_center = np.array(r)-np.array([surface.x0, surface.y0])
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

    elif surface.type == 'super_surface':
        return surface.intersect(r, u)

