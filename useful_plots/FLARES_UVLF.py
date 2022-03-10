
import numpy as np

import h5py


import flare.photom as phot

import flares_utility.analyse as analyse



X_limits = [27., 31.]
binw = 0.2
bin_edges = np.arange(*X_limits, binw)
bin_centres = bin_edges[:-1]+binw/2


# --- FLARES

flares = analyse.analyse(f'{analyse.flares_master_file}')

print(flares.tags)


# flares.list_datasets()

V = (4./3) * np.pi * (flares.radius)**3 # Mpc^3


for z, tag in zip(flares.zeds, flares.tags):


    ## ---- get data
    phi = np.zeros(len(bin_centres))
    N = np.zeros(len(bin_centres))

    # X = flares.load_dataset(tag, 'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI/', 'FUV')
    X = flares.load_dataset(tag, 'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/Intrinsic/', 'FUV')

    for i, (sim, w) in enumerate(zip(flares.sims, flares.weights)):

        x = np.log10(np.array(X[sim]))
        x = x[x>0.0]

        N_temp, _ = np.histogram(x, bins = bin_edges)

        N += N_temp

        phi_temp = (N_temp / V) / binw

        phi += phi_temp * w


    d = ', '.join([f'{p:.2e}' for p in phi])

    print(f'self.phi[{z}] = [{d}]')
