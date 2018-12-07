import numpy as np

class SuperRegion():
    """A class to contain a list of regions which can be queried in the same
    way. Introduces minor heirarchal concepts to speed up raytracing.

    Parameters
    ----------
    regions : list of Region objects

    Attributes
    ----------
    regions : list of Region objects
        Regions contained within the super region
    surfaces : list of Surface objects
        surfaces that define the exterior of the super region
    orientations : list of -1,1
        On which side of each surface the region is located
    uid : int
        Unique id of the super region
    """
    def __init__(self, regions, surfaces, orientations, uid):
        self.regions = regions
        self.surfaces = surfaces
        self.orientations = orientations
        self.uid = uid

    def evaluate(cls, r):
        """Check if r is located within the region
        
        Parameters
        ----------
        r : 2-tuple
            coordinates of point

        Returns
        -------
        bool
            Boolean on whether or not r is located in this super region
        """
        # Loop through all surfaces defining this super region
        # If none of the evaluations fail, point is in region
        for idx, surface in enumerate(cls.surfaces):
            sgn = cls.orientations[idx]
            if sgn*surface.evaluate(r) < 0:
                return False
        return True

    def get_uid(cls, r):
        """ Return uid of region containing r

        Parameters
        ----------
        r : 2-tuple
            Coordinates of point

        Returns
        -------
        int
            uid of region containing r
        """
        for region in cls.regions:
            if region.evaluate(r):
                return region.uid

class Region():
    """The Region class describes the region using CSG

    Parameters
    ----------
    surfaces : list of Surface objects
        The surfaces that define the region
    orientations : list of -1,1
        On which side of each surface the region is located
    mat : int, opt
        id of material contained within surface
    name : str, opt
        name of the region

    Attributes
    ----------
    surfaces : list of Surface objects
        The surfaces that define the region
    orientations : list of -1,1
        On which side of each surface the region is located
    uid : int
        Unique id of the region
    mat : str, opt
        name of material contained within surface
    """
    def __init__(self, surfaces, orientations, uid, mat = [], phi=0, a_phi=0, 
                 vol=0):
        self.surfaces = surfaces
        self.orientations = orientations
        self.uid = uid
        self.mat = mat
        self.phi = phi
        self.a_phi = phi
        self.vol = vol

        self.q = np.zeros([len(phi),])
        self.a_q = np.zeros([len(phi),])
        self.tracks_phi = np.zeros([len(phi),])
        self.tracks_a_phi = np.zeros([len(phi),])
        self.q_phi = np.zeros([len(phi),])
        self.a_q_phi = np.zeros([len(phi),])
        self.tot_track_length = 0

    def evaluate(cls, r):
        """Check if r is located within the region
        
        Parameters
        ----------
        r : 2-tuple
            coordinates of point

        Returns
        -------
        int 
            uid of the sub region in which the point is contained
        """
        # Loop through all surfaces defining this region
        # If none of the evaluations fail, point is in region
        for idx, surface in enumerate(cls.surfaces):
            sgn = cls.orientations[idx]
            # print('region: ', cls.uid)
            # print('surface: ', surface.id)
            # print('surface: ', surface.name)
            # print('r from center: ', [r[0]-1.26/2,r[1]-1.26/2])
            # print('orientation: ', sgn)
            # print('evaluation: ', surface.evaluate(r), 'r: ',r)
            if sgn*surface.evaluate(r) < 0:
                return False
        # print("in region", cls.name)
        return True

    def get_uid(cls, r):
        """ Return uid by form of a function. Allows variable to be called 
        in same method as super region

        Parameters
        ----------
        r : 2-tuple
            Coordinates of point (un-used)

        Returns
        -------
        int
            uid of region
        """
        return cls.uid

def what_region(r, regions):
    """ Find out what region the point is in 
    Parameters
    ----------
    r : 2-tuple
        coordinates
    regions : list of Region objects
        The regions being used in this problem 

    Returns
    -------
    str 
        name of the region 
    """
    for region in regions:
        if region.evaluate(r):
            return region.get_uid(r)
        
    
