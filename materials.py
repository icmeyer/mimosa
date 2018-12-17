import os 
import numpy as np

def import_xs(ngroup, pert = []):
   """Create a MATERIALS dictionary using files output by an OpenMC script
   Parameters
   ----------
   ngroup : int
       Number of energy groups

   pert : 2-tuple
      ['string of reaction', percent of perturbation]
   
   Returns
   -------
   MATERIALS : dict
       dictionary containing the cross section information for materials in the
       problem
   """

   mimosa_dir = os.path.dirname(os.path.realpath(__file__))
   folder = mimosa_dir+'/make_xs/xs_'+str(ngroup)+'group/'
   fuel_file = '_cell_1'
   mod_file = '_cell_0'
   MATERIALS = {'fuel': {'total': np.loadtxt(folder+'total'+fuel_file),
                         'nufission': np.loadtxt(folder+'nufission'+fuel_file),
                         'scatter': np.loadtxt(folder+'scatter'+fuel_file),
                         'chi': np.loadtxt(folder+'chi'+fuel_file),
                         'absorption': np.loadtxt(folder+'absorption'+fuel_file),
                         'cov' : {'total': [], 'nufission': [], 'chi': []}
                         },
                'mod': {'total': np.loadtxt(folder+'total'+mod_file),
                        'nufission': np.loadtxt(folder+'nufission'+mod_file),
                        'scatter': np.loadtxt(folder+'scatter'+mod_file),
                        'chi': np.loadtxt(folder+'chi'+mod_file),
                        'absorption': np.loadtxt(folder+'absorption'+mod_file),
                        'cov' : {'total': [], 'nufission': [], 'chi': []}
                        }
               }

   # Perturn cross-sections based on `pert` variable
   if pert != []:
       factor = 1 + pert[1]
       for mat in MATERIALS: 
           for xs in MATERIALS[mat]:
               if xs == 'scatter' and pert[0] == 'scatter':
                   MATERIALS[mat]['total'] -= np.diag(MATERIALS[mat][xs])
                   new_scat_diag = factor * np.diag(MATERIALS[mat][xs])
                   for i in range(ngroup):
                       MATERIALS[mat][xs][i,i] = new_scat_diag[i]
                   MATERIALS[mat]['total'] += np.diag(MATERIALS[mat][xs])
               if xs == 'nufission' and pert[0] == 'nuf':
                   MATERIALS[mat]['total'] -= MATERIALS[mat]['nufission']
                   MATERIALS[mat]['nufission'] = factor*MATERIALS[mat]['nufission']
                   MATERIALS[mat]['total'] += MATERIALS[mat]['nufission']

   return MATERIALS
