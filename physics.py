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
        sigma_t = MATERIALS[region.mat]['total']
        reduction = (1/4/pi/sigma_t)
        phi = region.phi
        region_fission_source = 0
        region.q = np.zeros([len(phi),])

        # Energy 'loop'
        # Calculate region fission source
        nuf = MATERIALS[region.mat]['nufission']
        fission_source += np.dot(nuf,phi)
        region_fission_source += np.dot(nuf,phi)

        # scatter is organized by [group out, group in]
        region.q += reduction*np.matmul(scatter,phi)
        # region.q += np.matmul(scatter,phi)

        # Distribute fission source using xi
        chi = MATERIALS[region.mat]['chi']
        region.q += reduction*chi*region_fission_source/k
        # region.q += chi*region_fission_source/k
        # print('region q reduced_q')
        # print([region.mat, region.q/reduction, region.q])


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
    for ray in rays:
        # Calculate initial psi
        region = regions[ray.segments[0].region]
        sigma_t = MATERIALS[region.mat]['total']
        reduction = (1/4/pi/sigma_t)
        # psi = regions[ray.segments[0].region].q*reduction
        psi = copy.deepcopy(regions[ray.segments[0].region].q)
        # print('psi q region r reduction')
        # print([psi,regions[ray.segments[0].region].q, ray.segments[0].region, ray.segments[0].r0, reduction])

        for segment in ray.segments:
            d = segment.d
            region = regions[segment.region]
            sigma_t = MATERIALS[region.mat]['total']
            reduction = (1/4/pi/sigma_t)
            
            #Energy "loop"
            q = region.q
            # Calculate delta_psi
            tau = sigma_t*d
            # delta_psi = (psi - q*reduction)*(1-np.exp(-tau))
            delta_psi = (psi - q)*(1-np.exp(-tau))
            # print('delta_psi')
            # print(delta_psi)
            # print('reduction')
            # print(reduction)
            # print('q times reduction')
            # print(q*reduction)


            if segment.active:
                region.tracks_phi += 4*pi*delta_psi
                
            psi -= delta_psi
            # print('delta_psi')
            # print(delta_psi)
            # print('psi')
            # print(psi)
            
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
