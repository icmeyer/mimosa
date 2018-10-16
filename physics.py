import numpy as np
from math import pi

from materials import MATERIALS

def normalize_phi(regions, ngroup):
    """Normalize the phi of the problem. Currently unused."""
    phi_sum = 0
    for region in regions:
        for group in range(ngroup):
            phi_sum += region.q[group]
    for region in regions:
        for group in range(ngroup):
            region.q[group] = region.q[group]/phi_sum

def normalize_q(regions, ngroup):
    """Normalize the q of the problem. Currently unused."""
    q_sum = 0
    for region in regions:
        for group in range(ngroup):
            q_sum += region.q[group]
    for region in regions:
        for group in range(ngroup):
            region.q[group] = region.q[group]/q_sum

def calc_q(regions, ngroup, k, update_k=False, old_fission_source=0):
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
    for idx, region in enumerate(regions):
        scatter = MATERIALS[region.mat]['scatter']
        phi = region.phi
        region_fission_source = 0
        region.q = np.zeros([len(phi),])

        for group in range(ngroup):
            # Calculate region fission source
            nuf = MATERIALS[region.mat]['nufission'][group]
            fission_source += nuf*phi[group]
            region_fission_source += nuf*phi[group]

            # Caclculate source from scattering
            for group_prime in range(ngroup):
                region.q[group] += scatter[group, group_prime]*phi[group_prime]

        #Distribute fission source using xi
        for group in range(ngroup):
            chi = MATERIALS[region.mat]['chi'][group]
            region.q[group] += chi*region_fission_source/k
        #     ####MIRIAM BIT
        #     if idx == 0:
        #         region.q[group] = 0
        #     elif idx == 1 and group == 1:
        #         region.q[group] = 0
        # print('region', idx, region.mat, region.q)


    if update_k:
        k = fission_source/old_fission_source
    # normalize_q(regions, ngroup)
    return fission_source/k, k


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

    for region in regions:
        for group in range(ngroup):
            sigma_t = MATERIALS[region.mat]['total'][group]
            # region.q_phi[group] += (4*pi/sigma_t)*region.q[group]
            region.q_phi[group] += (1/sigma_t)*region.q[group]

    for ray in rays:
        # Calculate initial psi
        region = regions[ray.segments[0].region]
        sigma_t = MATERIALS[region.mat]['total']
        psi = regions[ray.segments[0].region].q/(4*pi*sigma_t)

        for segment in ray.segments:
            d = segment.d
            region = regions[segment.region]
            sigma_t = MATERIALS[region.mat]['total']
            
            #Energy "loop"
            q = region.q
            # Calculate delta_psi
            tau = sigma_t*d
            delta_psi = (psi - (q/4/pi/sigma_t))*(1-np.exp(-tau))

            if segment.active:
                region.tracks_phi += 4*delta_psi
                
            psi -= delta_psi
            
            # WORKING
            # for group in range(ngroup):
            #     sigma_t = MATERIALS[region.mat]['total'][group]
            #     q = region.q[group]
            #     # Calculate delta_psi
            #     tau = sigma_t*d
            #     delta_psi = (psi[group] - (q/4/pi/sigma_t))*(1-np.exp(-tau))

            #     if segment.active:
            #         region.tracks_phi[group] += 4*delta_psi
            #         
            #     psi[group] -= delta_psi
            # END WORKING
        
    return rays
