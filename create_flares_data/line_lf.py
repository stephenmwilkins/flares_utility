
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.lines import Line2D

from astropy.table import Table, Column

import cmasher as cmr

import h5py

import flare.plt as fplt
import flare.photom as phot

import flares_utility.colors
import flares_utility.analyse as analyse

import flags_data.distribution_functions as df

tags = ['005_z010p000', '006_z009p000', '007_z008p000',
        '008_z007p000', '009_z006p000', '010_z005p000']
redshifts = [10, 9, 8, 7, 6, 5]


filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_noparticlesed_v3.hdf5'

flares = analyse.analyse(filename, default_tags=False)


lines = [
    ('OIII5007', ['OIII5007']),
    ('OIII', ('OIII5007', 'OIII4959')),
]

binw = 0.1
V = (4./3) * np.pi * (flares.radius)**3  # total volume of FLARES Mpc^3

for line in lines:

    line_name, lines_ = line

    print('-'*50, line_name)

    for lum_type, lt in zip(['DustModelI', 'Intrinsic'], ['', '_intrinsic']):

        # ---
        quantities = []
        quantities.append({'path': 'Galaxy/Mstar_aperture',
                          'dataset': f'30', 'name': 'Mstar', 'log10': True})

        for line_ in lines_:
            quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/{lum_type}/{line_}',
                              'dataset': f'Luminosity', 'name': line_, 'log10': False})

        t = Table()

        z_ = np.array([])
        bin_centres_ = np.array([])
        N_ = np.array([])
        phi_ = np.array([])

        for z in redshifts:

            tag = flares.tag_from_zed[z]

            D = flares.get_datasets(tag, quantities)
            L = np.zeros(len(D['log10Mstar']))

            for line_ in lines_:
                L += D[line_]

            log10Mstar = D['log10Mstar']
            w = D['weight']
            log10L = np.log10(L)

            # identify galaxies around our imposed stellar mass limit
            s = np.fabs(log10Mstar - 8.) < 0.1

            mn = np.round(np.percentile(log10L[s], 95), 1)
            mx = np.round(np.max(log10L), 1)

            bin_edges = np.arange(mn, mx, binw)
            bin_centres = bin_edges[:-1] + binw/2

            N, _ = np.histogram(log10L, bins=bin_edges)
            Neff, _ = np.histogram(log10L, bins=bin_edges, weights=w)

            print(z, mn, mx, N)

            phi = Neff/V/binw

            # --- stack

            z_ = np.hstack((z_, np.ones(len(N))*z))
            bin_centres_ = np.hstack((bin_centres_, bin_centres))
            N_ = np.hstack((N_, N))
            phi_ = np.hstack((phi_, phi))

        t.add_column(Column(data=z_, name='z'))
        t.add_column(Column(data=np.round(bin_centres_, 2), name='log10L', unit='dex(erg/s)',
                     description=rf'log_{{10}}(L_{{ {line_name} }}/erg\ s^{{-1}}'))
        t.add_column(Column(data=N_.astype(int), name='N'))
        t.add_column(Column(data=np.round(np.log10(phi_), 2),
                     name='log10phi', unit='dex(Mpc^-3 dex^-1)', description=rf'log_{{10}}(\phi/Mpc^{{-3}}\ dex^{{-1}})'))

        t.meta['name'] = 'FLARES'
        t.meta['references'] = []
        t.meta['lines'] = lines_

        t.write(
            f'/Users/stephenwilkins/Dropbox/Research/projects/flares/flares_data/flares_data/data/distribution_functions/{line_name}{lt}.ecsv', format='ascii.ecsv', overwrite=True)
