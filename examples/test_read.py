
import numpy as np

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from flares_utility import analyse


filename = analyse.flares_master_file+'/flares_highz_v3_nosed.hdf5' # if you don't have an environment variable set you need to give this the full path to the master file instead

a = analyse.analyse(filename, default_tags = False)

print(a.ids) # print simulation ids
print(a.tags) # print tag/snapshot ids
print(a.zeds) # print corresponding redshifts
print(a.deltas) # print delta


# ----------------------------------------------------------------------
# --- define parameters and tag

tag = a.tags[-1] # get last (lowest redshift) tag
z = a.zed_from_tag[tag] # get redshift of that tag
print(tag, z, a.tag_from_zed[z]) # print, but also shows how to the tag from zed


# ----------------------------------------------------------------------
# --- list datasets (specifically for the 1st sim/tag)
a.list_datasets()


# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []
quantities.append({'path': 'Galaxy', 'dataset': 'Mstar', 'name': None, 'log10': True})

# --- get quantities (and weights and deltas)
D = a.get_datasets(tag, quantities)

print(D.keys())
print(len(D['Mstar']))

P = a.get_particle_datasets(tag)

print(P.keys()) # print the quantities we've extracted
print(len(P['S_Age'])) # print the number of galaxes
print(P['S_Age'][0]) # print the ages of the star particles in the first galaxy
