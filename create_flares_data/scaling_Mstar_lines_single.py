
from synthesizer.sed import Lnu_to_M
import numpy as np
from astropy.table import Table, Column

import flares_utility.analyse as analyse
import flares_utility.labels as labels
import flares_utility.stats as stats

from unyt import unyt_quantity

tags = ['005_z010p000', '006_z009p000', '007_z008p000',
        '008_z007p000', '009_z006p000', '010_z005p000']
redshifts = [10, 9, 8, 7, 6, 5]


filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_noparticlesed_v3.hdf5'

flares = analyse.analyse(filename, default_tags=False)


lines = [
    'HI4861',
    'OIII4959',
    'OIII5007',
    ['OIII4959', 'OIII5007']
]


binw = 0.2
quantiles = [0.022, 0.158, 0.50, 0.842, 0.978]
percentiles = [np.round(q*100, 1) for q in quantiles]

for line in lines:

    # ---
    quantities = []
    quantities.append({'path': 'Galaxy/Mstar_aperture',
                       'dataset': f'30', 'name': 'Mstar', 'log10': True})
    quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/{line}',
                       'dataset': 'Luminosity', 'name': 'L', 'log10': True})
    quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/{line}',
                       'dataset': 'EW', 'name': 'EW', 'log10': True})

    t = Table()

    z_ = np.array([])
    bin_centres_ = np.array([])
    N_ = np.array([])
    lum_ = {}
    ew_ = {}

    for q in quantiles:
        lum_[q] = np.array([])
        ew_[q] = np.array([])

    for z in redshifts:

        tag = flares.tag_from_zed[z]

        D = flares.get_datasets(tag, quantities)
        X = D['log10Mstar']
        w = D['weight']
        log10L = D['log10L']
        log10EW = D['log10EW']

        bins = np.arange(8., np.max(X), binw)
        bin_centres = 0.5*(bins[:-1]+bins[1:])

        lum = stats.binned_weighted_quantile(X, log10L, w, bins, quantiles)
        ew = stats.binned_weighted_quantile(X, log10EW, w, bins, quantiles)
        N, _ = np.histogram(X, bins=bins)

        # --- make stacks
        z_ = np.hstack((z_, np.ones(len(N))*z))
        bin_centres_ = np.hstack((bin_centres_, bin_centres))
        N_ = np.hstack((N_, N))
        for i, q in enumerate(quantiles):
            lum_[q] = np.hstack((lum_[q], np.array(lum[:, i])))
            ew_[q] = np.hstack((ew_[q], np.array(ew[:, i])))

    t.add_column(Column(data=z_, name='z'))
    t.add_column(Column(data=np.round(bin_centres_, 2), name='log10Mstar', unit='dex(solMass)',
                 description=r'log_{{10}}(M_{\star}/M_{\odot})'))
    t.add_column(Column(data=N_.astype(int), name='N'))

    for q in quantiles:
        t.add_column(Column(data=np.round(
            lum_[q], 2), name=f'{line}_log10L_P{q*100:.1f}', unit='dex(erg/s)', description=rf'log_{{10}}(L_{{ {line} }}/erg\ s^{{-1}})'))
        t.add_column(Column(data=np.round(
            ew_[q], 2), name=f'{line}_EW_P{q*100:.1f}', unit='dex(AA)', description=rf'log_{{10}}(EW_{{ {line},0 }}/\AA)'))

    t.meta['name'] = 'FLARES'
    t.meta['references'] = []
    t.meta['percentiles'] = percentiles
    t.meta['description'] = f'predictions for the relationship betwee stellar mass and {line} luminosity and equivalent width'

    t.write(
        f'/Users/stephenwilkins/Dropbox/Research/projects/flares/flares_data/flares_data/data/scaling_relations/log10Mstar/{line}.ecsv', format='ascii.ecsv', overwrite=True)
