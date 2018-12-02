import numpy as np
from math import pi
import copy

from materials import MATERIALS

def normalize_phi(regions, ngroup):
    """Normalize the phi of the problem."""
    phi_sum = 0
    for region in regions:
        for group in range(ngroup):
            phi_sum += region.phi[group]
    for region in regions:
        for group in range(ngroup):
            region.phi[group] = region.phi[group]/phi_sum

def normalize_q(regions, ngroup):
    """Normalize the q of the problem. Currently unused."""
    q_sum = 0
    for region in regions:
        for group in range(ngroup):
            q_sum += region.q[group]
    for region in regions:
        for group in range(ngroup):
            region.q[group] = region.q[group]/q_sum

def calc_q(regions, ngroup, k, a_k, update_k=False, old_fission_source=0):
    """Updates the q's of the regions and returns a total q

    Parameters
    ----------
    regions : list of Region objects
    ngroup : int
        Number of energy groups
    k : float
        Current eigenvalue
    update_k : bool, opt
        Flag whether or not to calculate a new k
    old_fission_source : float, opt
        Sum of fission source from last q calculation

    Returns
    -------
    normalize_fission_source : float
        Fission source of this calculation normalized by k
    k : float
        Eigenvalue calculated from new and old source
    """

    fission_source = 0
    a_fission_source = 0 
    for idx, region in enumerate(regions):
        # Get materials properties for this region
        scatter = MATERIALS[region.mat]['scatter']
        a_scatter = scatter.transpose()
        sigma_t = MATERIALS[region.mat]['total']
        nuf = MATERIALS[region.mat]['nufission']
        chi = MATERIALS[region.mat]['chi']

        phi = region.phi
        a_phi = region.a_phi
        region_fission_source = 0
        region_a_fission_source = 0
        region.q = np.zeros([len(phi),])
        region.a_q = np.zeros([len(phi),])

        # Energy 'loop'
        # Calculate region fission source
        # Total tracker
        fission_source += np.dot(nuf, phi)
        a_fission_source += np.dot(nuf, a_phi)
        # Region tracker
        region_fission_source += np.dot(nuf, phi)
        region_a_fission_source += np.dot(nuf, a_phi)

        # Distribute fission source using xi
        q_fission = chi*region_fission_source/k
        a_q_fission = chi*region_a_fission_source/a_k
        # scatter is organized by [group out, group in]
        q_scatter = np.matmul(scatter,phi)
        a_q_scatter = np.matmul(a_scatter,phi)

        reduction = (1/4/pi/sigma_t)
        region.q += reduction*(q_scatter + q_fission)
        region.a_q += reduction*(a_q_scatter + a_q_fission)

    if update_k:
        k = fission_source/old_fission_source
        a_k = a_fission_source/old_fission_source
    # normalize_q(regions, ngroup)
    return fission_source/k, k, a_k


def ray_contributions(rays, ngroup, regions):
    """ Update the phi using method of characteristics along
    the stored segments

    Parameters
    ---------
    rays : list of Ray objects
        Ray objects populated with segments
    ngroup : int
        Number of energy groups
    regions : list of Region objects

    Returns
    -------
    rays : list of Ray objects
    """
    for ray in rays:
        # Calculate initial psi
        psi = copy.deepcopy(regions[ray.segments[0].region].q)
        a_psi = copy.deepcopy(regions[ray.segments[0].region].a_q)

        for segment in ray.segments:
            d = segment.d
            region = regions[segment.region]
            sigma_t = MATERIALS[region.mat]['total']
            
            #Energy "loop" performed by vectors here
            q = region.q
            a_q = region.a_q
            tau = sigma_t*d
            exp_term = (1-np.exp(-tau))
            delta_psi = (psi - q)*exp_term
            delta_a_psi = (a_psi - a_q)*exp_term

            if segment.active:
                region.tracks_phi += 4*pi*delta_psi
                region.tracks_a_phi += 4*pi*delta_a_psi
                
            psi -= delta_psi
            
    return rays
