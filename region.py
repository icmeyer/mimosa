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
    name : str, opt
        name of the region
    mat : int, opt
        id of material contained within surface
    """
    def __init__(self, surfaces, orientations, name = '', mat = []):
        self.surfaces = surfaces
        self.orientations = orientations
        self.name = name
        self.mat = mat

    def evaluate(cls, r):
        """Check if r is located within the region
        
        Parameters
        ----------
        r : 2-tuple
            coordinates of point

        Returns
        -------
        bool
            Boolean on whether or not r is located in this region
        """
        # Loop through all surfaces defining this region
        # If none of the evaluations fail, point is in region
        for idx, surface in enumerate(cls.surfaces):
            sgn = cls.orientations[idx]
            # print('region: ', cls.name)
            # print('surface: ', surface.name)
            # print('r from center: ', [r[0]-1.26/2,r[1]-1.26/2])
            # print('orientation: ', sgn)
            # print('evaluation: ', surface.evaluate(r))
            if sgn*surface.evaluate(r) < 0:
                return False
        # print("in region", cls.name)
        return True

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
            return region.name
        
    
