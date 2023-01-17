
from synthesizer.sed import log10L_to_M
import numpy as np
from astropy.table import Table, Column
import flares_utility.analyse as analyse

tags = ['005_z010p000', '006_z009p000', '007_z008p000',
        '008_z007p000', '009_z006p000', '010_z005p000']
redshifts = [10, 9, 8, 7, 6, 5]


filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_noparticlesed_v3.hdf5'

flares = analyse.analyse(filename, default_tags=False)


binw = 0.1
V = (4./3) * np.pi * (flares.radius)**3  # total volume of FLARES Mpc^3


# ---
quantities = []
quantities.append({'path': 'Galaxy/Mstar_aperture',
                  'dataset': f'30', 'name': 'Mstar', 'log10': True})


t = Table()

z_ = np.array([])
bin_centres_ = np.array([])
N_ = np.array([])
phi_ = np.array([])

for z in redshifts:

    tag = flares.tag_from_zed[z]

    D = flares.get_datasets(tag, quantities)

    log10Mstar = D['log10Mstar']
    w = D['weight']

    bin_edges = np.arange(8., 11.5, binw)
    bin_centres = bin_edges[:-1] + binw/2

    N, _ = np.histogram(log10Mstar, bins=bin_edges)
    Neff, _ = np.histogram(log10Mstar, bins=bin_edges, weights=w)

    phi = Neff/V/binw

    z_ = np.hstack((z_, np.ones(len(N))*z))
    bin_centres_ = np.hstack((bin_centres_, bin_centres))
    N_ = np.hstack((N_, N))
    phi_ = np.hstack((phi_, phi))

    # --- find the bin where the number of galaxies falls below 5
    # s = N > 5

    # --- stack
    # z_ = np.hstack((z_, np.ones(len(N[s]))*z))
    # bin_centres_ = np.hstack((bin_centres_, bin_centres[s]))
    # N_ = np.hstack((N_, N[s]))
    # phi_ = np.hstack((phi_, phi[s]))


t.add_column(Column(data=z_, name='z'))
t.add_column(Column(data=np.round(bin_centres_, 2), name='log10Mstar', unit='dex(Msol)',
             description=r'log_{10}(M_{\star}/M_{\odot})'))
t.add_column(Column(data=N_.astype(int), name='N'))
t.add_column(Column(data=np.round(np.log10(phi_), 2),
             name='log10phi', unit='dex(Mpc^-3 dex^-1)', description=rf'log_{{10}}(\phi/Mpc^{{-3}}\ dex^{{-1}})'))

t.meta['name'] = 'FLARES'
t.meta['references'] = []


t.write(
    f'/Users/stephenwilkins/Dropbox/Research/projects/flares/flares_data/flares_data/data/distribution_functions/log10Mstar.ecsv', format='ascii.ecsv', overwrite=True)
