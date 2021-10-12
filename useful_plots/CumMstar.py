
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

import cmasher as cmr

import h5py

from flares_utility import analyse

import flare.plt




Mstar_limits = [8., 11]

fig = plt.figure(figsize = (3.5,3.5))

left  = 0.15
bottom = 0.15
height = 0.8
width = 0.8

ax = fig.add_axes((left, bottom, width, height))


redshifts = [5,6,7,8,9,10]

snap = {}
snap['EAGLE'] = {5:'008_z005p037',6: '006_z005p971',7:'005_z007p050',8:'004_z008p075',9:'003_z008p988',10:'002_z009p993'}
snap['FLARES'] = {5:'010_z005p000',6: '009_z006p000',7:'008_z007p000',8:'007_z008p000',9:'006_z009p000',10:'005_z010p000'}


EAGLE = h5py.File('/Users/stephenwilkins/Dropbox/Research/data/simulations/FLARES/EAGLE_REF_sp_info.hdf5', 'r')

fl = analyse.analyse('/Users/stephenwilkins/Dropbox/Research/data/simulations/FLARES/flares_no_particlesed.hdf5')

quantities = [{'path': 'Galaxy', 'dataset': 'Mstar_30', 'name': None, 'log10': True}]



for sim, ls in zip(['EAGLE','FLARES'],[':','-']):

    for z, c in zip(redshifts, cmr.take_cmap_colors('cmr.gem_r', len(redshifts))):

        if sim == 'EAGLE':
            log10Mstar = np.log10(EAGLE[snap['EAGLE'][z]]['Galaxy']['Mstar_30'][:]) + 10.
            label = None

        if sim == 'FLARES':
            D = fl.get_datasets(snap['FLARES'][z], quantities, return_weights = False)
            log10Mstar = D['log10Mstar_30']
            label = rf'$\rm z={z}$'

        N = []

        M = np.arange(*Mstar_limits, 0.1)

        for log10Mlimit in M:
            N.append(np.log10(len(log10Mstar[log10Mstar>log10Mlimit])))

        ax.plot(M, N, ls = ls, c=c, label = label)



ax.legend(labelspacing=0.05,fontsize=8)
ax.set_xlim(Mstar_limits)
# ax.set_ylim([-12., -1.51])


ax.set_xlabel(r'$\rm \log_{10}(M_{\star}/M_{\odot})$')
ax.set_ylabel(r'$\rm\log_{10}[N(>M_{\star})]$')


fig.savefig(f'figures/CumMstar.pdf')


fig.clf()
