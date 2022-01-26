
import numpy as np

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from flares_utility import analyse



filename = analyse.flares_master_file+'/flares_highz_v3_nosed.hdf5' # if you don't have an environment variable set you need to give this the full path to the master file instead

a = analyse.analyse(filename, default_tags = False)





# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []
quantities.append({'path': 'Galaxy', 'dataset': 'Mstar', 'name': None, 'log10': True})
quantities.append({'path': 'Galaxy/SFR', 'dataset': '100Myr', 'name': 'SFR', 'log10': True})



# --- get quantities (and weights and deltas)

for tag in a.tags:
    D = a.get_datasets(tag, quantities)
    print(tag, len(D['Mstar']), len(D['SFR']))
