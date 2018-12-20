import os
import sys
import copy
mimosa_dir = os.path.dirname(os.getcwd())
sys.path.insert(0, mimosa_dir)

import numpy as np
from main import main 

from surface import SuperSurface, XPlane, YPlane, Circle
from tools import normalize
from plotting import *
from region import SuperRegion, Region, what_region
from ray import Ray, make_segments
from physics import calc_q
from materials import import_xs

ngroup = 10
pitch = 1.26
length = 1.26*3
limits = [0, length, 0, length]
# Surfaces starting with 100 will be circles
# Surfaces starting with 200 will be planes
surfaces = []
super_surfaces = []
circ_counter = 100
xcoords = [2*1.26+1.26/2, 1*1.26+1.26/2, 0*1.26+1.26/2]
ycoords = xcoords

fuel_refine_level = 3
# Biggest radius first
radii = np.linspace(0,0.39218,fuel_refine_level+1)[:0:-1]
print(radii)

for i in range(3):
    for j in range(3):
        sub_surfaces = []
        for level in range(fuel_refine_level):
            circ_counter += 1
            circle = Circle(surface_id=circ_counter, boundary_type='transmission',
                            x0=xcoords[i], y0=ycoords[j], R=radii[level]) #0.39218
            surfaces.append(circle)
            if level == 0:
                super_surface_ext = circle
            sub_surfaces.append(circle)
        super_surface = SuperSurface(circ_counter, surfaces, super_surface_ext,
                                     boundary_type='transmission')
        super_surfaces.append(super_surface)


left = XPlane(surface_id=201, boundary_type='reflection', x0=0)
xplane1 = XPlane(surface_id=202, boundary_type='transmission', x0=pitch)
xplane2 = XPlane(surface_id=203, boundary_type='transmission', x0=2*pitch)
right = XPlane(surface_id=204, boundary_type='reflection', x0=length)
top = YPlane(surface_id=205, boundary_type='reflection', y0=length)
yplane1 = YPlane(surface_id=206, boundary_type='transmission', y0=2*pitch)
yplane2 = YPlane(surface_id=207, boundary_type='transmission', y0=pitch)
bottom = YPlane(surface_id=208, boundary_type='reflection', y0=0)


surfaces += [left, xplane1, xplane2, right, top, yplane1, yplane2, bottom]
super_surfaces += [left, xplane1, xplane2, right, top, yplane1, yplane2, bottom]
# Regions: current limitation, uid of region must be its order in the list
regions = []
super_regions = []
region_counter = 0
fuel_phi_guess = np.ones([ngroup,])
for i in range(9):
    if i == 4:
        mat = 'mod'
    else:
        mat = 'fuel'
    sub_regions = []
    outer = fuel_refine_level*i
    inner = outer+fuel_refine_level-1
    # Make donuts
    for level in range(fuel_refine_level - 1):
        fuel = Region([surfaces[outer+level], surfaces[outer+level+1]], [-1, 1], uid=region_counter, mat=mat,phi=fuel_phi_guess)
        region_counter += 1
        regions.append(fuel)
        sub_regions.append(fuel)
    # Make donut hole
    fuel = Region([surfaces[inner]], [-1], uid=region_counter, mat=mat,phi=fuel_phi_guess)
    region_counter += 1
    regions.append(fuel)
    sub_regions.append(fuel)
    super_region = SuperRegion(sub_regions, [surfaces[outer]], [-1], uid=i)
    super_regions.append(super_region)

mod_phi_guess = np.ones([ngroup,])*0.1
mod1 = Region([left, xplane1, top, yplane1, surfaces[fuel_refine_level*0]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
mod2 = Region([xplane1, xplane2, top, yplane1, surfaces[fuel_refine_level*1]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
mod3 = Region([xplane2, right, top, yplane1, surfaces[fuel_refine_level*2]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
mod4 = Region([left, xplane1, yplane1, yplane2, surfaces[fuel_refine_level*3]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
mod5 = Region([xplane1, xplane2, yplane1, yplane2, surfaces[fuel_refine_level*4]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
mod6 = Region([xplane2, right, yplane1, yplane2, surfaces[fuel_refine_level*5]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
mod7 = Region([left, xplane1, yplane2, bottom, surfaces[fuel_refine_level*6]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
mod8 = Region([xplane1, xplane2, yplane2, bottom, surfaces[fuel_refine_level*7]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
mod9 = Region([xplane2, right, yplane2, bottom, surfaces[fuel_refine_level*8]],[1, -1, -1, 1, 1],
                    uid=region_counter, mat='mod', phi = mod_phi_guess)
region_counter += 1
regions += [mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9]
super_regions += [mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9]

for region in regions:
    print(region, region.mat)

n_rays = 50
# main(n_rays, surfaces, regions, limits, ngroup, super_regions = super_regions,
#      super_surfaces = [], plot=False)
k0, a_k0, regions, rays = main(n_rays, surfaces, regions, 
                                limits, ngroup, 
                                super_regions = copy.deepcopy(super_regions), 
                                super_surfaces = [], plot=False)

# Check that the using the same rays and regions produces the same results
# k0, a_k0, regions2, rays = main(n_rays, surfaces, copy.deepcopy(regions), 
#                                 limits, ngroup, 
#                                 super_regions = copy.deepcopy(super_regions), 
#                                 super_surfaces = [], plot=False, k_guess = k0,
#                                 rays = rays)
# for i in range(len(regions)):
#     region1 = regions[i]
#     region2 = regions2[i]
#     print(region1.phi)
#     print(region2.phi)
#     print(np.allclose(region1.phi, region2.phi))


# Save the fluxes from first calculation
for region in regions:
    region.phi_0 = region.phi
    region.a_phi_0 = region.a_phi

# Build Sensitivity Vector
index = 0


# Calculate denominator
bottom = 0
MATERIALS = import_xs(ngroup, pert = [])
for region in regions:
    nuf = MATERIALS[region.mat]['nufission']
    chi = MATERIALS[region.mat]['chi']
    phi = region.phi
    a_phi = region.a_phi
    bottom += np.inner(a_phi, chi*np.dot(nuf, phi))

# Bottom should be in a loop over all regions and energy groups
total_sens = len(regions)*ngroup*2
sens_vec = np.zeros(total_sens)
sens_vec_counter = 0
region_list = [0,1,2,3,4,5,6,12,13,14,27,28,29,30,31,32,33,34,35]
total_runs = len(region_list)*10*2

for region_num in region_list:
    region_pert = regions[region_num]
    for i in range(ngroup):
        print('Sensitivity Progress:',sens_vec_counter,'/',total_runs)
        print(sens_vec_counter<5)
        material = region_pert.mat
        # Use perturbed material
        if material == 'fuel':
            region_pert.mat = 'fuel1' 
        elif material == 'mod':
            region_pert.mat = 'mod1'
        pert_vec = np.zeros(ngroup)
        pert_vec[i] = 0.05

        # Finding absorption sensitivites
        pert = ['absorption', pert_vec] # Perturbation instructions
        # Run perturbed simulation
        pert_k, a_k, regions, rays = main(n_rays, surfaces, regions, 
                                        limits, ngroup, 
                                        super_regions = copy.deepcopy(super_regions), 
                                        super_surfaces = [], plot=False, pert=pert, 
                                        k_guess = k0, rays = rays)
        top = 0
        for region2 in regions:
            a_phi = region2.a_phi_0
            phi = region2.phi_0
            phi_pert = region2.phi
            trans_pert = phi_pert - phi
            if region_pert.mat[0] == 'f':
                nuf = MATERIALS['fuel']['nufission']
                chi = MATERIALS['fuel']['chi']
            elif region_pert.mat[0] == 'm':
                nuf = MATERIALS['mod']['nufission']
                chi = MATERIALS['mod']['chi']
            nuf_pert = MATERIALS[region_pert.mat]['nufission']
            chi_pert = MATERIALS[region_pert.mat]['chi']
            fiss_pert = (1/k0)*(chi*np.dot(nuf_pert, phi_pert) - chi*np.dot(nuf, phi))
            pert_val = MATERIALS['pert_val']
            # top += np.inner(a_phi, (trans_pert - fiss_pert)/pert_val )
            top += np.inner(a_phi, (trans_pert - fiss_pert))

        d_lambda = top/bottom
        print(d_lambda)
        sens_vec[region_num*10*2+i*2] = d_lambda
        sens_vec_counter += 1

        # Finding nufission sensitivites
        pert = ['nuf', pert_vec] # Perturbation instructions
        # Run perturbed simulation
        pert_k, a_k, regions, rays = main(n_rays, surfaces, regions, 
                                        limits, ngroup, 
                                        super_regions = copy.deepcopy(super_regions), 
                                        super_surfaces = [], plot=False, pert=pert, 
                                        k_guess = k0, rays = rays)
        top = 0
        for region2 in regions:
            a_phi = region2.a_phi_0
            phi = region2.phi_0
            phi_pert = region2.phi
            trans_pert = phi_pert - phi
            if region_pert.mat[0] == 'f':
                nuf = MATERIALS['fuel']['nufission']
                chi = MATERIALS['fuel']['chi']
            elif region_pert.mat[0] == 'm':
                nuf = MATERIALS['mod']['nufission']
                chi = MATERIALS['mod']['chi']
            nuf_pert = MATERIALS[region_pert.mat]['nufission']
            chi_pert = MATERIALS[region_pert.mat]['chi']
            fiss_pert = (1/k0)*(chi*np.dot(nuf_pert, phi_pert) - chi*np.dot(nuf, phi))
            pert_val = MATERIALS['pert_val']
            # top += np.inner(a_phi, (trans_pert - fiss_pert)/pert_val )
            top += np.inner(a_phi, (trans_pert - fiss_pert))

        d_lambda = top/bottom
        print(d_lambda)
        sens_vec[region_num*10*2+i*2+1] = d_lambda
        sens_vec_counter += 1
        # Reset to unperturbed material
        if material == 'fuel1':
            region_pert.mat = 'fuel' 
        elif material == 'mod1':
            region_pert.mat = 'mod'

np.savetxt('./sensitivity_vector_before_symmetry', sens_vec)
# Apply symmetry
# Corner pins
pinstart = 6
sens_vec[pinstart*10*2:(pinstart+3)*10*2] = sens_vec[0*10*2:3*10*2]
pinstart = 18
sens_vec[pinstart*10*2:(pinstart+3)*10*2] = sens_vec[0*10*2:3*10*2]
pinstart = 24
sens_vec[pinstart*10*2:(pinstart+3)*10*2] = sens_vec[0*10*2:3*10*2]
#Side pins
pinstart = 9
sens_vec[pinstart*10*2:(pinstart+3)*10*2] = sens_vec[3*10*2:6*10*2]
pinstart = 15
sens_vec[pinstart*10*2:(pinstart+3)*10*2] = sens_vec[3*10*2:6*10*2]
pinstart = 21
sens_vec[pinstart*10*2:(pinstart+3)*10*2] = sens_vec[3*10*2:6*10*2]
np.savetxt('./sensitivity_vector', sens_vec)

