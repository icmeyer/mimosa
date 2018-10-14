import numpy as np
from math import pi

# MATERIALS = {'fuel': {'total': [0.37, 0.66], 
#                       'nufission': [0.01, 0.38],
#                       'scatter': np.array([[0.24, 0.0],[0.072, 1.35]]), 
#                       'xi': [1.0, 0] },
#              'mod':  {'total': [0.68, 2.2], 
#                       'nufission': [0.0, 0.0],
#                       'scatter': np.array([[0.30, 0.0],[0.001, 0.42]]), 
#                       'xi': [0.0, 0.0]} }
def import_xs(folder):
    """Create a MATERIALS dictionary using files output by OpenMC"""
    fuel_file = '_cell_1'
    mod_file = '_cell_0'
    MATERIALS = {'fuel': {'total': np.loadtxt(folder+'total'+fuel_file),
                          'nufission': np.loadtxt(folder+'nufission'+fuel_file),
                          'scatter': np.loadtxt(folder+'scatter'+fuel_file),
                          'xi': np.loadtxt(folder+'xi'+fuel_file)
                          },
                 'mod': {'total': np.loadtxt(folder+'total'+mod_file),
                         'nufission': np.loadtxt(folder+'nufission'+mod_file),
                         'scatter': np.loadtxt(folder+'scatter'+mod_file),
                         'xi': np.loadtxt(folder+'xi'+mod_file)
                        }
                }
    return MATERIALS

MATERIALS = import_xs('./make_xs/xs/')

def calc_q(regions, ngroup, k):
    """Updates the q's of the regions and returns a total q"""
    total_source = 0
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
                region.q[group] += (1/4/pi)*scatter[group, group_prime]*phi[group_prime]
                print('source after scatter', region.q)
        #Distribute fission source using xi
        for group in range(ngroup):
            xi = MATERIALS[region.mat]['xi'][group]
            region.q[group] += (1/4/pi)*region_fission_source
            print('source after fission dist', region.q)
            total_source += region.q[group]
            print('total_source', total_source)

    #Normalize phi and q by q
    checksum = 0
    for region in regions:
        for group in range(ngroup):
            region.q[group] = region.q[group]/total_source
            checksum += region.q[group]
            region.phi[group] = region.phi[group]/total_source
        print("This should be 1", checksum)

    return fission_source


def ray_contributions(rays, regions):
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
                delta_psi = (psi[group] - q/sigma_t)*(1-np.exp(-tau))
                ###NEW
                # psi_1 = psi[group]*np.exp(-tau)+q/sigma_t*(1-np.exp(-tau))
                # delta_psi = psi_1 - psi[group]
                # psi[group] = psi_1

                if segment.active:
                    # Do I need a sin term here?
                    region.iter_phi[group] += delta_psi/(vol)
                psi -= delta_psi
        
    return rays


