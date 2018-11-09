import numpy as np

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






    
