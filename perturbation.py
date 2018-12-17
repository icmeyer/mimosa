import numpy as np

def calc_perturbation(regions, MATERIALS, k, a_k):
    """
    Calculate the change in reactivity using manufactured uncertainties
    """
    pert_dict = {'percent': np.array([0.01, 0.03, 0.05]), 'nuf' : [], 'scatter': []}
    print('------------Perturbation------------')
    top = 0
    bottom = 0
    for n in [0.01, 0.03, 0.05]:
        for region in regions:
            d_sigma_s = MATERIALS[region.mat]['scatter']*n
            nuf = MATERIALS[region.mat]['nufission']
            chi = MATERIALS[region.mat]['chi']
            phi = region.phi
            a_phi = region.a_phi
            top += np.inner(region.a_phi, np.matmul(d_sigma_s, phi))
            bottom += np.inner(a_phi, chi*np.dot(nuf, phi))
        d_lamb = top/bottom
        print('Sigma_s', n, d_lamb)
        pert_dict['scatter'].append(d_lamb)
            
    top = 0
    bottom = 0
    for n in [0.01, 0.03, 0.05]:
        for region in regions:
            d_nuf = MATERIALS[region.mat]['nufission']*n
            nuf = MATERIALS[region.mat]['nufission']
            chi = MATERIALS[region.mat]['chi']
            phi = region.phi
            a_phi = region.a_phi
            top += np.inner(region.a_phi, -(1/k)*chi*np.dot(d_nuf,phi))
            bottom += np.inner(a_phi, chi*np.dot(nuf, phi))
        d_lamb = top/bottom
        print('nuf', n, d_lamb)
        pert_dict['nuf'].append(d_lamb)
    return pert_dict
