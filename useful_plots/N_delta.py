
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

import cmasher as cmr

import h5py

import flare.plt as fplt
import flare.photom as phot

import flares_utility.analyse as analyse













filename = analyse.flares_master_file+'/flares_highz_v3_nosed.hdf5'
flares = analyse.analyse(filename, default_tags = False)

# redshifts = [11,12,13,14,15]
# redshifts = [11,13,15]
redshifts = flares.zeds

print(redshifts)


V = (4./3) * np.pi * (flares.radius)**3 # Mpc^3

xlim = [-0.3, 0.35]
ylim = [-0.75, 1.7]


# --- all redshifts on one figure


fig, ax = fplt.simple()

for z, c in zip(redshifts, cmr.take_cmap_colors('cmr.gem_r', len(redshifts))):

    tag = flares.tag_from_zed[z]

    X = flares.load_dataset(tag, f'Galaxy', 'Mstar')


    N = []
    for sim in flares.sims:
        x = np.log10(np.array(X[sim])) +10.
        N.append(np.sum(x>8.))

    N = np.array(N)
    N_avg = np.sum(N*flares.weights)/np.sum(flares.weights)

    ax.axhline(np.log10(N_avg), c=c, alpha = 0.3)
    ax.scatter(np.log10(1+flares.deltas), np.log10(N), s=5, c=[c], label = rf'$\rm z={z}$')


ax.set_xlim(xlim)
ax.set_ylim(ylim)

ax.legend(fontsize=8)

ax.set_xlabel(r'$\rm\log_{10}(1+\delta)$')
ax.set_ylabel(r'$\rm \log_{10}[N(M_{\star}>10^{8}\ M_{\odot})]$')

fig.savefig(f'figs/N_delta.pdf')

fig.clf()



# --- all redshifts on one figure but normalised by average density


fig, ax = fplt.simple()

for z, c in zip(redshifts, cmr.take_cmap_colors('cmr.gem_r', len(redshifts))):

    tag = flares.tag_from_zed[z]

    X = flares.load_dataset(tag, f'Galaxy', 'Mstar')


    N = []
    for sim in flares.sims:
        x = np.log10(np.array(X[sim])) +10.
        N.append(np.sum(x>8.))

    N = np.array(N)
    N_avg = np.sum(N*flares.weights)/np.sum(flares.weights)

    ax.scatter(np.log10(1+flares.deltas), np.log10(N/N_avg), s=5, c=[c], label = rf'$\rm z={z}$')


ax.set_xlim(xlim)
ax.set_ylim(ylim)

ax.legend(fontsize=8)

ax.set_xlabel(r'$\rm\log_{10}(1+\delta)$')
ax.set_ylabel(r'$\rm \log_{10}[N(M_{\star}>10^{8}\ M_{\odot})/\bar{N}]$')

fig.savefig(f'figs/N_delta_mean.pdf')

fig.clf()



# --- individual redshifts

for z, c in zip(redshifts, cmr.take_cmap_colors('cmr.gem_r', len(redshifts))):

    fig, ax = fplt.simple()

    tag = flares.tag_from_zed[z]

    X = flares.load_dataset(tag, f'Galaxy', 'Mstar')


    N = []
    for sim in flares.sims:
        x = np.log10(np.array(X[sim])) +10.
        N.append(np.sum(x>8.))

    N = np.array(N)
    N_avg = np.sum(N*flares.weights)/np.sum(flares.weights)

    ax.axhline(np.log10(N_avg), c=c, alpha = 0.3)
    ax.scatter(np.log10(1+flares.deltas), np.log10(N), s=5, c=[c], label = rf'$\rm z={z}$')


    ax.set_ylim([-0.75, 1.7])

    ax.legend(fontsize=8)

    ax.set_xlabel(r'$\rm\log_{10}(1+\delta)$')
    ax.set_ylabel(r'$\rm \log_{10}[N(M_{\star}>10^{8}\ M_{\odot})]$')

    fig.savefig(f'figs/N_delta_{z}.pdf')

    fig.clf()
