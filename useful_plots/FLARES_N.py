
import numpy as np

import h5py


import flare.photom as phot

import flares_utility.analyse as analyse



# --- FLARES

flares = analyse.analyse(f'{analyse.flares_master_file}')

print(flares.tags)


quantities = [{'path': 'Galaxy/', 'dataset': 'Mstar', 'name': None, 'log10': True}]


limits = [7., 8., 9., 10., 11.]

for z, tag in zip(flares.zeds, flares.tags):


    X_ = flares.get_datasets(tag, quantities)
    X = X_['log10'+quantities[0]['dataset']]

    print(f'N[{z}]={{ '+', '.join([f'{limit} : {np.sum(X>limit)}' for limit in limits])+'}')




    # d = ', '.join([f'{p:.2e}' for p in phi])
    #
    # print(f'self.phi[{z}] = [{d}]')
