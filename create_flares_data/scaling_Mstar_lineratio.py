

# NOT WORKING YET


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


line_ratio, lines = 'R3', ['OIII5007', 'HI4861']


binw = 0.2
quantiles = [0.022, 0.158, 0.50, 0.842, 0.978]
percentiles = [np.round(q*100, 1) for q in quantiles]


# ---
quantities = []
quantities.append({'path': 'Galaxy/Mstar_aperture',
                   'dataset': f'30', 'name': 'Mstar', 'log10': True})

for i, line in enumerate(lines):
    quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/{line}',
                       'dataset': 'Luminosity', 'name': f'L{i}', 'log10': False})

t = Table()

z_ = np.array([])
bin_centres_ = np.array([])
N_ = np.array([])
lr_ = {}  #  line_ratio

for q in quantiles:
    lr_[q] = np.array([])

for z in redshifts:

    tag = flares.tag_from_zed[z]

    D = flares.get_datasets(tag, quantities)
    X = D['log10Mstar']
    w = D['weight']

    R = D['L0']/D['L1']

    s = D['log10Mstar'] > 8
    print(z, np.median(R[s]))

    bins = np.arange(8., np.max(X), binw)
    bin_centres = 0.5*(bins[:-1]+bins[1:])

    lr = stats.binned_weighted_quantile(X, R, w, bins, quantiles)
    N, _ = np.histogram(X, bins=bins)

    # --- make stacks
    z_ = np.hstack((z_, np.ones(len(N))*z))
    bin_centres_ = np.hstack((bin_centres_, bin_centres))
    N_ = np.hstack((N_, N))
    for i, q in enumerate(quantiles):
        lr_[q] = np.hstack((lr_[q], np.array(lr[:, i])))

t.add_column(Column(data=z_, name='z'))
t.add_column(Column(data=np.round(bin_centres_, 2), name='log10Mstar', unit='dex(solMass)',
             description=r'log_{{10}}(M_{\star}/M_{\odot})'))
t.add_column(Column(data=N_.astype(int), name='N'))

for q in quantiles:
    t.add_column(Column(data=np.round(
        lr_[q], 2), name=f'{line_ratio}_P{q*100:.1f}', description=rf' {line_ratio} ({lines[0]}/{lines[1]}) line ratio'))

t.meta['name'] = 'FLARES'
t.meta['references'] = []
t.meta['percentiles'] = percentiles
t.meta[
    'description'] = f'predictions for the relationship betwee stellar mass and {line_ratio} ({lines[0]}/{lines[1]}) line ratio'

t.write(
    f'/Users/stephenwilkins/Dropbox/Research/projects/flares/flares_data/flares_data/data/scaling_relations/log10Mstar/{line_ratio}.ecsv', format='ascii.ecsv', overwrite=True)
