
import numpy as np
import matplotlib.pyplot as plt

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from flares_utility import analyse

z = 11

filename = analyse.flares_master_file+'/flares_highz_v3_nosed.hdf5' # if you don't have an environment variable set you need to give this the full path to the master file instead

a = analyse.analyse(filename, default_tags = False)

# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []
quantities.append({'path': 'Galaxy/SFR', 'dataset': '100Myr', 'name': 'SFR', 'log10': True}) # SFR averaged over 100 Myr
quantities.append({'path': 'Galaxy', 'dataset': 'Mstar', 'name': None, 'log10': True}) # Stellar mass




# --- get quantities (and weights and deltas)
D = a.get_datasets(a.tag_from_zed[z], quantities)

# ----------------------------------------------
# define new quantities
D['log10sSFR'] = D['log10SFR']-D['log10Mstar']+9 # +9 converts to Gyr

# ----------------------------------------------
# define selection
s = D['log10Mstar']>8.0 # only select galaxies with Mstar>10^9

# ----------------------------------------------
# Print number of galaxies meeting the selection
print(f"Total number of galaxies: {len(D['log10Mstar'][s])}")


# ----------------------------------------------
# Make a basic scatter plot

fig = plt.figure(figsize = (4, 4))

left  = 0.2
height = 0.75
bottom = 0.2
width = 0.75

ax = fig.add_axes((left, bottom, width, height))

ax.scatter(D['log10Mstar'][s], D['log10sSFR'][s], s=1, alpha=1)

ax.set_xlabel(r'$\rm M_{\star}/M_{\odot}$')
ax.set_ylabel(r'$\rm sSFR/Gyr^{-1}$')

plt.show()
