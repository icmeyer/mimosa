import numpy as np
from math import pi

# MATERIALS = {'fuel': {'total': [0.37, 0.66], 
#                       'nufission': [0.01, 0.38],
#                       'scatter': np.array([[0.24, 0.0],[0.072, 1.35]]), 
#                       'xi': [1.0, 0] },
#              'mod':  {'total': [0.68, 2.2], 
#                       'nufission': [0.0, 0.0],
#                       'scatter': np.array([[0.30, 0.0],[0.001, 0.42]]), 

def import_xs(folder):
   """Create a MATERIALS dictionary using files output by OpenMC"""
   fuel_file = '_cell_1'
   mod_file = '_cell_0'
   MATERIALS = {'fuel': {'total': np.loadtxt(folder+'total'+fuel_file),
                         'nufission': np.loadtxt(folder+'nufission'+fuel_file),
                         'scatter': np.loadtxt(folder+'scattering'+fuel_file),
                         'chi': np.loadtxt(folder+'chi'+fuel_file)
                         },
                'mod': {'total': np.loadtxt(folder+'total'+mod_file),
                        'nufission': np.loadtxt(folder+'nufission'+mod_file),
                        'scatter': np.loadtxt(folder+'scattering'+mod_file),
                        'chi': np.loadtxt(folder+'chi'+mod_file)
                       }
               }
   return MATERIALS

MATERIALS = import_xs('./make_xs/xs/')

def normalize_phi(regions, ngroup):
    phi_sum = 0
    for region in regions:
        for group in range(ngroup):
            phi_sum += region.q[group]
    for region in regions:
        for group in range(ngroup):
            region.q[group] = region.q[group]/phi_sum

def normalize_q(regions, ngroup):
    q_sum = 0
    for region in regions:
        for group in range(ngroup):
            q_sum += region.q[group]
    for region in regions:
        for group in range(ngroup):
            region.q[group] = region.q[group]/q_sum

def calc_q(regions, ngroup, k, update_k=False, old_fission_source=0):
    """Updates the q's of the regions and returns a total q"""
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
        # if any(region.q < 0):
        #     print(region.mat)
        #     print('q', region.q)
        #     print('phi', region.phi)
        #     raise ValueError('whattttt')

    if update_k:
        k = fission_source/old_fission_source

    normalize_q(regions, ngroup)

    return fission_source/k, k


def ray_contributions(rays, regions):
    # Phi from source term
    for region in regions:
        # print('before and after for region', region.mat)
        for group in range(rays[0].ngroup):
            sigma_t = MATERIALS[region.mat]['total'][group]
            # region.q_phi[group] += (4*pi/sigma_t)*region.q[group]
            region.q_phi[group] += (1/sigma_t)*region.q[group]

    for ray in rays:
        # Calculate initial psi
        psi = regions[ray.segments[0].region].q/(4*pi)
        for segment in ray.segments:
            d = segment.d
            for group in range(ray.ngroup):
                region = regions[segment.region]
                vol = region.vol
                sigma_t = MATERIALS[region.mat]['total'][group]
                q = region.q[group]
                # Calculate delta_psi
                tau = sigma_t*d
                ###OLD
                delta_psi = (psi[group] - (q/4/pi/sigma_t)*(1-np.exp(-tau)))
                ###NEW
                # psi_1 = psi[group]*np.exp(-tau)+q/sigma_t*(1-np.exp(-tau))
                # delta_psi = psi_1 - psi[group]
                # psi[group] = psi_1

                if segment.active:
                    # Do I need a sin term here?
                    region.tracks_phi[group] += (1/vol/sigma_t)*delta_psi
                    # sin_term = np.sin(ray.varphi)
                    # region.iter_phi[group] += (4*pi/vol/sigma_t)*delta_psi*d*sin_term
                    # region.iter_phi[group] -= delta_psi/(sigma_t*vol)
                    # if any(region.iter_phi < 0):
                    #     print(region.mat)
                    #     print('group', group)
                    #     print('psi', psi)
                    #     print('d', segment.d)
                    #     print('mu', ray.mu)
                    #     print('delta_psi', delta_psi)
                    #     print('q', region.q)
                    #     print('phi', region.iter_phi)
                    #     raise ValueError('whattttt')
                    
                psi -= delta_psi
                # if any(psi < 0):
                #     print('psi', psi)
                #     print('q', q)
                #     print('sigma_t', sigma_t)
                #     print('region', region.mat)
                #     raise ValueError('whattttt')
        
    return rays


