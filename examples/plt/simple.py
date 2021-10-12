
import numpy as np
import matplotlib.cm as cm

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from flares_utility import analyse
from flares_utility import plt as plots
from flares_utility import limits as limits_
from flares_utility import labels as labels_


# ----------------------------------------------------------------------
# --- load analysis object
filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_highz.hdf5'
a = analyse.analyse(filename, default_tags = False)

# --- define the redshift (used below to select the tag)
z = 10
log10Mstar_limit = 9.0


# --- define quantities to read in [not those for the corner plot, that's done later]
quantities = []
quantities.append({'path': 'Galaxy/SFR_total', 'dataset': 'SFR_100_Myr', 'name': 'SFR', 'log10': True}) # SFR averaged over 100 Myr
quantities.append({'path': 'Galaxy', 'dataset': 'Mstar', 'name': None, 'log10': True}) # Stellar mass

# --- get quantities (and weights and deltas)
D = a.get_datasets(a.tag_from_zed[z], quantities)

# ----------------------------------------------
# define new quantities
D['log10sSFR'] = np.log10(D['SFR'])-np.log10(D['Mstar'])+9

# ----------------------------------------------
# define selection
s = D['log10Mstar']>log10Mstar_limit # only select galaxies meeting the criteria



# ----------------------------------------------
# ----------------------------------------------
# plot with colour bar


# --- get default limits and modify them to match the selection range
limits = limits_.limits
limits['log10Mstar_30'] = [log10Mstar_limit, 10.9]

# --- get default labels and modify them as required
labels = labels_.labels
labels['log10sSFR'] = 'specific\ star\ formation\ rate/Gyr^{-1}' # this is just to demonstrate how to use this


x = 'log10Mstar'
y = 'log10sSFR'

# --- make plot with colour bar plot
fig = plots.simple(D, x, y, s, limits = limits, labels = labels)


fig.savefig(f'figs/simple.pdf')
