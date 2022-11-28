
import numpy as np
from astropy.table import Table, Column, vstack
import astropy.units as u
import h5py

import flares_utility.analyse as analyse





redshifts = [5,6,7,8,9,10,11,12,13,14,15]



x, x_, x_unit, range, binw, conversion = 'log10Mstar', ('Galaxy/Mstar_aperture', '30'), 'dex(Msun)', [8.,12.], 0.2, 10
# x, x_, x_unit, range, binw, conversion = 'log10SFR', ('Galaxy/SFR_aperture/30', '50Myr'), 'dex(Msun yr^-1)', [-0.5, 3.5], 0.2, 0.0

bin_edges =  np.arange(*range, binw)
bin_centres = 0.5*(bin_edges[:-1]+bin_edges[1:])


# --- FLARES
flares = analyse.analyse('/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_noparticlesed_v2.hdf5', default_tags = False)

print(flares.tags)
# flares.list_datasets()

V = (4./3) * np.pi * (flares.radius)**3 # Mpc^3

tables = []

for z in redshifts:

    tag = flares.tag_from_zed[z]

    ## ---- get data
    phi = np.zeros(len(bin_centres))
    N = np.zeros(len(bin_centres))

    Xs = flares.load_dataset(tag, *x_)
    print(z, '-'*40)

    for i, (sim, w) in enumerate(zip(flares.sims, flares.weights)):

        X = np.log10(np.array(Xs[sim])) + conversion
        X = X[X>-99.]

        if len(X)>0:

            print(sim, np.min(X), np.max(X))

            N_temp, _ = np.histogram(X, bins = bin_edges)

            N += N_temp

            phi_temp = (N_temp / V) / binw

            phi += phi_temp * w


    s = N>0
    if np.sum(s)>0:
        log10phi = np.log10(phi)

        t = Table()
        t.add_column(Column(data = z*np.ones(np.sum(s)), name = 'z'))
        t.add_column(Column(data = np.round(bin_centres[s],2), name = x, unit = x_unit))
        t.add_column(Column(data = N[s].astype(int), name = 'N'))
        t.add_column(Column(data = np.round(log10phi[s],2), name = 'log10phi', unit = 'dex(Mpc^-3 dex^-1)'))
        tables.append(t)


table = vstack(tables)

table.meta['x'] = x
table.meta['y'] = 'log10phi'
table.meta['name'] = 'FLARES'
table.meta['redshifts'] = list(set(table['z']))
table.meta['type'] = 'binned'
table.meta['references'] = ['2021MNRAS.500.2127L', '2022arXiv220409431W']

table.write(f'tables/flares_{x}DF.ecsv', format = 'ascii.ecsv', overwrite=True)
table.write(f'/Users/stephenwilkins/Dropbox/Research/modules/flare_data/flags_data/data/DistributionFunctions/{x[5:]}/models/binned/flares.ecsv', format = 'ascii.ecsv', overwrite=True)
# table.write(f'models/binned/flares.ecsv', format = 'ascii.ecsv', overwrite=True)
