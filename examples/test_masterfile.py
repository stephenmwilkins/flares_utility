
import numpy as np

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from flares_utility import analyse


filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_highz.hdf5'

a = analyse.analyse(filename, default_tags = False)





# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []
quantities.append({'path': 'Galaxy', 'dataset': 'Mstar', 'name': None, 'log10': True})
quantities.append({'path': 'Galaxy/SFR_total', 'dataset': 'SFR_100_Myr', 'name': 'SFR', 'log10': True})



# --- get quantities (and weights and deltas)

for tag in a.tags:
    D = a.get_datasets(tag, quantities)
    print(tag, len(D['Mstar']), len(D['SFR']))
