import os 
import numpy as np

def import_xs(ngroup):
   """Create a MATERIALS dictionary using files output by an OpenMC script
   Parameters
   ----------
   ngroup : int
       Number of energy groups
   
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
                         'cov' : {'total': [], 'nufission': [], 'chi': []}
                         },
                'mod': {'total': np.loadtxt(folder+'total'+mod_file),
                        'nufission': np.loadtxt(folder+'nufission'+mod_file),
                        'scatter': np.loadtxt(folder+'scatter'+mod_file),
                        'chi': np.loadtxt(folder+'chi'+mod_file),
                        'cov' : {'total': [], 'nufission': [], 'chi': []}
                        }
               }
   # Add covariance data, assume 5%
   # for mat in MATERIALS: 
   #     for xs in MATERIALS[mat]:
   #         if xs != 'scatter':
   #             MATERIALS[mat][cov][xs] = np.diag(0.05*MATERIALS[mat][xs])

   return MATERIALS
