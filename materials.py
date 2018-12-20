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
                        },
               'fuel1': {'total': np.loadtxt(folder+'total'+fuel_file),
                         'nufission': np.loadtxt(folder+'nufission'+fuel_file),
                         'scatter': np.loadtxt(folder+'scatter'+fuel_file),
                         'chi': np.loadtxt(folder+'chi'+fuel_file),
                         'absorption': np.loadtxt(folder+'absorption'+fuel_file),
                         'cov' : {'total': [], 'nufission': [], 'chi': []}
                         },
                'mod1': {'total': np.loadtxt(folder+'total'+mod_file),
                        'nufission': np.loadtxt(folder+'nufission'+mod_file),
                        'scatter': np.loadtxt(folder+'scatter'+mod_file),
                        'chi': np.loadtxt(folder+'chi'+mod_file),
                        'absorption': np.loadtxt(folder+'absorption'+mod_file),
                        'cov' : {'total': [], 'nufission': [], 'chi': []}
                        },
                'pert_val': 0 
               }

   # Perturn cross-sections based on `pert` variable
   if pert != []:
       factor = np.ones(ngroup) + pert[1]
       for mat in ['fuel','mod']: 
           pert_mat = mat+str(1)
           for xs in ['absorption', 'nufission']:
               if xs == 'absorption' and pert[0] == 'absorption':
                   MATERIALS[pert_mat]['total'] -= MATERIALS[pert_mat]['absorption']
                   MATERIALS['pert_val'] = np.sum(pert[1]*MATERIALS[pert_mat]['absorption'])
                   MATERIALS[pert_mat]['absorption'] = factor*MATERIALS[pert_mat]['absorption']
                   MATERIALS[pert_mat]['total'] += MATERIALS[pert_mat]['absorption']
               if xs == 'nufission' and pert[0] == 'nuf':
                   MATERIALS[pert_mat]['total'] -= MATERIALS[pert_mat]['nufission']
                   MATERIALS['pert_val'] = np.sum(pert[1]*MATERIALS[pert_mat]['nufission'])
                   MATERIALS[pert_mat]['nufission'] = factor*MATERIALS[pert_mat]['nufission']
                   MATERIALS[pert_mat]['total'] += MATERIALS[pert_mat]['nufission']

   return MATERIALS
