


# --- this example shows how to get and use star particle information

import numpy as np

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from flares_utility import analyse

filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_highz.hdf5'

a = analyse.analyse(filename, default_tags = False)



# ----------------------------------------------------------------------
# --- define parameters and tag

tag = a.tags[0] # get first (highest redshift) tag


# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []
quantities.append({'path': 'Galaxy', 'dataset': 'Mstar', 'name': None, 'log10': True})

# --- get quantities (and weights and deltas)
D = a.get_datasets(tag, quantities)

print('the integrated galaxy properties that have been extracted:', D.keys())
print('the TOTAL number of galaxies:', len(D['Mstar']))
s = D['log10Mstar']>8.
print('the number of galaxies with M*>1E8:', len(D['Mstar'][s]))

P = a.get_particle_datasets(tag)

print('the particle datasets:', P.keys()) # print the quantities we've extracted

# P is dictionary, with an element for each quantity.
# S_Age: the ages of the star particles
# S_Mass: the masses of the star particles
# S_MassInitial: the initial masses of the star Particles
# S_Z: the metallicity of the star particles (i.e. the mass fraction in metals)
# For each quantity there is a list of arrays, with one entry for each galaxy.

# This each array is a difference length we can't **currently** slice P like this: P['S_Age'][s]


for i in np.arange(len(s))[s]: # this loops over the indices of galaxies for which s is True (i.e. log10Mstar>8)
    print('-'*30)
    print('-'*20, i)
    print(f"integrated stellar mass: {D['log10Mstar'][i]:.2f}")
    print(f"summed particle stellar mass: {np.log10(np.sum(P['S_Mass'][i])):.2f}")
    print(f"number of star particles: {len(P['S_Mass'][i])}")
    print(f"oldest star particle/Myr: {np.max(P['S_Age'][i]):.0f}")
    print(f"mean star particle age/Myr: {np.mean(P['S_Age'][i]):.0f}")
    print(f"median star particle age/Myr: {np.median(P['S_Age'][i]):.0f}")
