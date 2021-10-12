

# Calculate the ionising luminosity of a FLARES galaxy using synthobs



import sys
import os

import numpy as np

import synthobs
from synthobs.sed import models

import flare
import flare.filters

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from flares_utility import analyse





filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_highz.hdf5'

a = analyse.analyse(filename, default_tags = False)

z = 10
tag = a.tag_from_zed[z]

i = 0 # galaxy index


P = a.get_particle_datasets(tag)

print('Total number of galaxies:', len(P['S_Age']))

print('-'*20)
print(f'Particles for galaxy {i}:')
print('Ages:', P['S_Age'][i])
print('Initial Masses:', P['S_MassInitial'][i])
print('Metallicities:', P['S_Z'][i])




# --- initialise SED grid ---
#  this can take a long time do don't do it for every object
model = models.define_model('BPASSv2.2.1.binary/ModSalpeter_300') # DEFINE SED GRID -

# --- calculate the ionising photon luminosity
log10Q = models.generate_log10Q(model, P['S_MassInitial'][i],  P['S_Age'][i],  P['S_Z'][i])

print('-'*20)
print(f'log10(M*): {np.log10(np.sum(P["S_MassInitial"][i]))}')
print(f'log10Q (direct): {log10Q}')

# --- calculate the full SEDs
o = models.generate_SED(model, P['S_MassInitial'][i], P['S_Age'][i], P['S_Z'][i], fesc = 1.0) # looking for the pure stellar luminosities so fesc=1

# --- calculate the ionising photon luminosity from integrating the SED
log10Q_SED = o.stellar.return_log10Q()
print(f'log10Q (SED): {log10Q_SED}')


# --- calculate the pure stellar UV luminosity
f = 'FAKE.TH.FUV'
F = flare.filters.add_filters([f], new_lam = model.lam) # --- define the filters. FAKE.FAKE are just top-hat filters using for extracting rest-frame quantities.
model.create_Lnu_grid(F) # --- create new L grid for each filter. In units of erg/s/Hz
Lnu = models.generate_Lnu(model, F, P['S_MassInitial'][i], P['S_Age'][i], P['S_Z'][i], fesc = 1.0)
log10LFUV = np.log10(Lnu[f]) # erg/s/Hz

print(f'log10LFUV: {log10LFUV}')
print(f'ionising photon production efficiecy: {log10Q-log10LFUV}')
