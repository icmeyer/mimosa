import numpy as np
from math import pi
import copy

from materials import import_xs

def normalize_phi(regions, ngroup, adjoint=False):
    """Normalize the phi of the problem."""
    phi_sum = 0
    a_phi_sum = 0
    for region in regions:
        for group in range(ngroup):
            phi_sum += region.phi[group]
            if adjoint:
                a_phi_sum += region.a_phi[group]

    for region in regions:
        for group in range(ngroup):
            region.phi[group] = region.phi[group]/phi_sum
            if adjoint:
                region.a_phi[group] = region.a_phi[group]/a_phi_sum

def normalize_q(regions, ngroup):
    """Normalize the q of the problem. Currently unused."""
    q_sum = 0
    for region in regions:
        for group in range(ngroup):
            q_sum += region.q[group]
    for region in regions:
        for group in range(ngroup):
            region.q[group] = region.q[group]/q_sum

def calc_q(regions, ngroup, k, a_k, update_k=False, old_fission_source=0,
           old_a_fission_source=0, pert=[]):
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
    MATERIALS = import_xs(ngroup, pert=pert)

    fission_source = 0
    a_fission_source = 0 
    tot_phi = 0
    tot_absorption = 0
    a_tot_absorption = 0
    tot_fission_source = 0
    tot_fission_source_test = 0
    a_tot_fission_source = 0
    for idx, region in enumerate(regions):
        # Get materials properties for this region
        scatter = MATERIALS[region.mat]['scatter']
        a_scatter = scatter.transpose()
        sigma_t = MATERIALS[region.mat]['total']
        nuf = MATERIALS[region.mat]['nufission']
        chi = MATERIALS[region.mat]['chi']
        absor = MATERIALS[region.mat]['absorption']

        phi = region.phi
        tot_phi += np.sum(region.phi)*region.vol
        a_phi = region.a_phi
        region.q = np.zeros([len(phi),])
        region.a_q = np.zeros([len(phi),])

        # Energy 'loop'
        # Calculate region fission source
        
        # Region tracker
        region_fission_source = np.dot(nuf, phi)
        region_a_fission_source = np.dot(chi, a_phi)
        region_absorption = np.dot(absor, phi)
        a_region_absorption = np.dot(absor, a_phi)

        # Distribute fission source using chi
        q_fission_uncorrected = chi*region_fission_source
        q_fission = q_fission_uncorrected/k
        # print('forward: chi, nuf, fission_source', chi, nuf, phi, q_fission)
        a_q_fission_uncorrected = nuf*region_a_fission_source
        a_q_fission = a_q_fission_uncorrected/a_k
        # print('adjoint: chi, nuf, fission_source', chi, nuf, a_phi, a_q_fission)

        # Total tracker
        fission_source += np.sum(q_fission_uncorrected)
        tot_fission_source += region.vol * np.sum(q_fission_uncorrected)
        tot_fission_source_test += np.sum(q_fission_uncorrected)
        a_tot_fission_source += region.vol * np.sum(a_q_fission_uncorrected)
        a_fission_source += np.sum(a_q_fission_uncorrected)

        # print('fission source in q calc:', fission_source, a_fission_source)
        tot_absorption += np.sum(region_absorption)*region.vol
        a_tot_absorption += np.sum(a_region_absorption)*region.vol
        # a_tot_absorption += np.sum(a_region_absorption)
        # print('Values the same?', np.dot(nuf,phi), np.sum(q_fission)*k)
        # scatter is organized by [group out, group in]
        if ngroup > 1:
            q_scatter = np.matmul(scatter, phi)
            a_q_scatter = np.matmul(a_scatter, a_phi)
        else:
            q_scatter = scatter*phi
            a_q_scatter = a_scatter*phi

        reduction = (1/4/pi/sigma_t)
        region.q += reduction*(q_scatter + q_fission)
        region.a_q += reduction*(a_q_scatter + a_q_fission)
        # print('qs', region.q, region.a_q)
        # print('scatter sources', q_scatter, a_q_scatter)
        # print('fission sources', q_fission, a_q_fission)
        # print('q_fission sources before k', q_fission, a_q_fission)
        # print('ks', k, a_k)

    k_from_balance = tot_fission_source/tot_absorption
    a_k_from_balance = a_tot_fission_source/a_tot_absorption
    if update_k:
        # print('fission source at k update - old', old_fission_source, old_a_fission_source)
        # print('fission source at k update - new', fission_source, a_fission_source)
        k = fission_source/old_fission_source
        a_k = a_fission_source/old_a_fission_source
        # print('forward k/balance', k, k_from_balance)
        # print('forward balance factors:', tot_fission_source, tot_absorption)
        # print('adjoint a_k/balance', a_k, a_k_from_balance)
        # print('adjoint balance factors:', a_tot_fission_source, a_tot_absorption)

    # normalize_q(regions, ngroup)
    return fission_source/k, a_fission_source/a_k, k_from_balance, a_k


def ray_contributions(rays, ngroup, regions, pert=[]):
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
    MATERIALS = import_xs(ngroup, pert=pert)

    for ray in rays:
        # Calculate initial psi
        psi = copy.deepcopy(regions[ray.segments[0].region].q)
        # a_psi = copy.deepcopy(regions[ray.segments[0].region].a_q)
        a_psi = copy.deepcopy(regions[ray.segments[0].region].q)

        for segment in ray.segments:
            if segment.from_vacuum:
                psi = 0
                a_psi = 0

            d = segment.d
            region = regions[segment.region]
            sigma_t = MATERIALS[region.mat]['total']
            
            #Energy "loop" performed by vectors here
            q = region.q
            a_q = region.a_q
            tau = sigma_t*d
            exp_term = (1-np.exp(-tau))
            # a_exp_term = (1-np.exp(-tau[::-1]))
            delta_psi = (psi - q)*exp_term
            delta_a_psi = (a_psi - a_q)*exp_term
            # print('tau', tau)
            # print('delta terms', delta_psi, delta_a_psi)

            if segment.active:
                region.tracks_phi += 4*pi*delta_psi
                region.tracks_a_phi += 4*pi*delta_a_psi

                
            psi -= delta_psi
            a_psi -= delta_a_psi
            # print('psis', psi, a_psi)
            
    return rays
