from math import pi
import numpy as np

import openmc
import openmc.mgxs
import openmc.model


# # Materials and Geometry

# ## Geometry parameters

# In[2]:


radius_fuel = 0.3922
radius_gap = 0.4001
radius_clad = 0.4572
pitch = 1.26


# ## Basic materials
# 
# The basic materials for this problem are 3.2% enriched UO2, Zircaloy-4, and water with 650 ppm boron

# In[3]:


uo2 = openmc.Material(name='fuel')
uo2.add_element('U', 1, enrichment=3.2)
uo2.add_element('O', 2)
uo2.set_density('g/cc', 10.341)


# In[4]:


zirc = openmc.Material()
zirc.add_element('O', 0.125, 'wo')
zirc.add_element('Cr', 0.10, 'wo')
zirc.add_element('Fe', 0.21, 'wo')
zirc.add_element('Zr', 98.115, 'wo')
zirc.add_element('Sn', 1.45, 'wo')
zirc.set_density('g/cm3', 6.55)


# In[5]:


water = openmc.model.borated_water(650)


# ## Homogenized water/cladding
# 
# We want to make some comparisons to textbook slowing-down theory which only works for two-region problems so we're going to smear the gap, cladding, and water together into one material.  We want to make sure we conserve the overall mass of each materail in the problem when we do this.

# In[6]:


area_gap = pi * (radius_gap**2 - radius_fuel**2)
area_clad = pi * (radius_clad**2 - radius_gap**2)
area_channel = pitch**2 - pi * radius_clad**2
area_homog = area_gap + area_clad + area_channel


# In[7]:


homogenized = openmc.Material(name='moderator')
homogenized.add_s_alpha_beta('c_H_in_H2O')

for name, frac in (zirc.get_nuclide_atom_densities().values()):
    homogenized.add_nuclide(name, frac * area_clad / area_homog)

for name, frac in (water.get_nuclide_atom_densities().values()):
    homogenized.add_nuclide(name, frac * area_channel / area_homog)


# In[8]:


materials = openmc.Materials([uo2, homogenized])


# ## Surfaces

# In[9]:


rfo = openmc.ZCylinder(R=radius_fuel)
xy_box = openmc.model.get_rectangular_prism(pitch, pitch, boundary_type='reflective')
z0 = openmc.ZPlane(z0=-10, boundary_type='reflective')
z1 = openmc.ZPlane(z0=10, boundary_type='reflective')


# ## Cells, etc.

# In[10]:


fuel = openmc.Cell(name='fuel', fill=uo2)
fuel.region = -rfo & +z0 & -z1
mod = openmc.Cell(name='moderator', fill=homogenized)
mod.region = +rfo & xy_box & +z0 & -z1
root = openmc.Universe(cells=(fuel, mod))
geometry = openmc.Geometry(root)


# # Settings

# In[11]:


settings = openmc.Settings()


# In[12]:


#TODO: how many batches/particles are needed for good cross sections?
settings.batches = 60
settings.inactive = 10
settings.particles = 1000


# In[13]:


# Set the initial source to a flat distribution born only in the fuel.
space = openmc.stats.Box((-radius_fuel, -radius_fuel, 0),
     (radius_fuel, radius_fuel, 0), only_fissionable=True)
settings.source = openmc.Source(space=space)


# # MGXS Tallies

# In[14]:


# Define the 8-group energy bounds.
ngroup = 1
if ngroup==10:
    e_groups = [0.0, 0.058, 0.14, 0.28, 0.625, 4.0, 10.0, 40.0, 5530.0, 821e3, 20e6]
elif ngroup==2:
    e_groups = [0.0, 1000, 20e6]
elif ngroup==1:
    e_groups = [0.0, 20e6]
groups = openmc.mgxs.EnergyGroups(e_groups)


# In[15]:


# Instantiate an MGXS library.
mgxs_lib = openmc.mgxs.Library(geometry)
mgxs_lib.energy_groups = groups


# In[16]:


# Don't apply any anisotropic scattering corrections.
mgxs_lib.correction = None


# In[17]:


# Set the desired MGXS data.
mgxs_lib.mgxs_types = ('total', 'scatter matrix', 'chi', 'nu-fission')


# In[18]:


# Define the domain and build the library.
mgxs_lib.domain_type = 'cell'
mgxs_lib.domains = geometry.get_all_material_cells().values()
mgxs_lib.build_library()


# In[19]:


# Add the tallies.
tallies = openmc.Tallies()
mgxs_lib.add_to_tallies_file(tallies)


# In[20]:


# Export to XML and run.
materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()
tallies.export_to_xml()
openmc.run()


# In[21]:


# Load the statepoint and the MGXS results.
sp_file = 'statepoint.{}.h5'.format(settings.batches)
sp = openmc.StatePoint(sp_file)
mgxs_lib.load_from_statepoint(sp)


# In[22]:


# Pick out the fuel and moderator domains.
fuel = mgxs_lib.domains[0]
moderator = mgxs_lib.domains[1]
assert fuel.name == 'fuel'
assert moderator.name == 'moderator'

#Save the cross-sections
df = mgxs_lib.get_mgxs(fuel, 'total').get_pandas_dataframe()
total_mat = df['mean'].values
np.savetxt('./xs/total_cell_1', total_mat)

df = mgxs_lib.get_mgxs(fuel, 'scatter matrix').get_pandas_dataframe()
scatter_col = df['mean'].values
scatter_mat = np.reshape(scatter_col, [ngroup,ngroup]).T
np.savetxt('./xs/scatter_cell_1', scatter_mat)

df = mgxs_lib.get_mgxs(fuel, 'chi').get_pandas_dataframe()
chi_mat = df['mean'].values
np.savetxt('./xs/chi_cell_1', chi_mat)

df = mgxs_lib.get_mgxs(fuel, 'nu-fission').get_pandas_dataframe()
nufission_mat = df['mean'].values
np.savetxt('./xs/nufission_cell_1', nufission_mat)

df = mgxs_lib.get_mgxs(moderator, 'total').get_pandas_dataframe()
total_mat = df['mean'].values
np.savetxt('./xs/total_cell_0', total_mat)

df = mgxs_lib.get_mgxs(moderator, 'scatter matrix').get_pandas_dataframe()
scatter_col = df['mean'].values
scatter_mat = np.reshape(scatter_col, [ngroup,ngroup]).T
np.savetxt('./xs/scatter_cell_0', scatter_mat)

df = mgxs_lib.get_mgxs(moderator, 'chi').get_pandas_dataframe()
chi_mat = df['mean'].values
np.savetxt('./xs/chi_cell_0', chi_mat)

df = mgxs_lib.get_mgxs(moderator, 'nu-fission').get_pandas_dataframe()
nufission_mat = df['mean'].values
np.savetxt('./xs/nufission_cell_0', nufission_mat)


