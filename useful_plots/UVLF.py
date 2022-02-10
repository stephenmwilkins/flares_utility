
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

import cmasher as cmr

import h5py

import flare.plt as fplt
import flare.photom as phot

import flares_utility.analyse as analyse



X_limits = [28, 30.5]

fig = plt.figure(figsize = (3.5,3.5))

left  = 0.15
bottom = 0.15
height = 0.8
width = 0.8

ax = fig.add_axes((left, bottom, width, height))


redshifts = [5,6,7,8,9,10]








binw = 0.1
bin_edges = np.arange(*X_limits, binw)
bin_centres = bin_edges[:-1]+binw/2



# --- EAGLE REF

EAGLE = h5py.File('/Users/stephenwilkins/Dropbox/Research/data/simulations/FLARES/EAGLE_REF_sp_info.hdf5', 'r')
snap = {5:'008_z005p037',6: '006_z005p971',7:'005_z007p050',8:'004_z008p075',9:'003_z008p988',10:'002_z009p993'}
V = 100**3

for z, c in zip(redshifts, cmr.take_cmap_colors('cmr.gem_r', len(redshifts))):

    X = np.log10(EAGLE[snap[z]]['Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI/FUV'][:])

    hist, _ = np.histogram(X, bins = bin_edges)

    phi = (hist / V) / binw

    ax.plot(bin_centres, np.log10(phi), ls = '--', c=c)




# --- FLARES

flares = analyse.analyse('/Users/stephenwilkins/Dropbox/Research/data/simulations/FLARES/flares_no_particlesed.hdf5', default_tags = False)

# flares.list_datasets()

V = (4./3) * np.pi * (flares.radius)**3 # Mpc^3


for z, c in zip(redshifts, cmr.take_cmap_colors('cmr.gem_r', len(flares.zeds))):


    tag = flares.tag_from_zed[z]

    ## ---- get data
    phi = np.zeros(len(bin_centres))
    N = np.zeros(len(bin_centres))

    X = flares.load_dataset(tag, 'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI/', 'FUV')


    for i, (sim, w) in enumerate(zip(flares.sims, flares.weights)):

        x = np.log10(np.array(X[sim]))
        x = x[x>0.0]

        N_temp, _ = np.histogram(x, bins = bin_edges)

        N += N_temp

        phi_temp = (N_temp / V) / binw

        phi += phi_temp * w

    ax.plot(bin_centres, np.log10(phi), ls = '-', c=c, label = rf'$\rm z={z}$')



ax.legend(labelspacing=0.05,fontsize=8)
ax.set_xlim(X_limits)
# ax.set_ylim([-12., -1.51])


ax.set_xlabel(r'$\rm \log_{10}(L_{FUV}/erg\ s^{-1}\ Hz^{-1})$')
ax.set_ylabel(r'$\rm\log_{10}[\phi/Mpc^{-3}\ dex^{-1}]$')


fig.savefig(f'figures/UVLF.pdf')


fig.clf()
