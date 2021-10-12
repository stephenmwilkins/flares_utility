
import numpy as np

from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

import cmasher as cmr

from astropy.cosmology import Planck15 as cosmo

import flare.plt as fplt
import flare.photom as phot

import flares_utility.analyse as analyse



redshifts = [5,7,9]


Y_limits = [5.5, 11.75]

X_limits = [27.5, 31]
binw = 0.1
bin_edges = np.arange(X_limits[0]-0.5, X_limits[1], binw)
bin_centres = bin_edges[:-1]+binw/0.5


Surveys = {}

# --- Euclid
Surveys['Euclid'] = {}
Surveys['Euclid']['Deep'] = {'area': 40*3600, 'Hlimit': 26.0, 'public': True}
Surveys['Euclid']['Wide'] = {'area': 18000*3600, 'Hlimit': 24.0, 'public': True}

# --- Roman
Surveys['Roman'] = {}
Surveys['Roman']['HLS'] = {'area': 2000*3600, 'Hlimit': 26.9, 'public': True}

# --- Webb
Surveys['Webb'] = {}
Surveys['Webb']['COSMOS'] = {'area': 0.6*3600, 'Hlimit': 27.6, 'public': True}







fig = plt.figure(figsize = (3.5,3.5))

left  = 0.15
bottom = 0.15
height = 0.8
width = 0.8

ax = fig.add_axes((left, bottom, width, height))


for z, c in zip(redshifts, cmr.take_cmap_colors('cmr.gem_r', len(redshifts))):

    v1, v2 = cosmo.comoving_volume([z-0.5, z+0.5]).to('Mpc3').value

    for Survey in Surveys.keys():
        for subSurvey, SS in Surveys[Survey].items():
            log10L_limit = np.log10(phot.m_to_lum(SS['Hlimit'],z, cosmo=cosmo))# absolute magnitude at z'=z
            log10V = np.log10((v2 - v1) * (SS['area']/(41253.*3600))) # co-moving volume z'=z-0.5, z+0.5

            # ax.fill_between([M_limit, -30],[0,0],[log10V, log10V],color='k',alpha=0.1)
            ax.plot([log10L_limit, 50], [log10V]*2, color = c, lw = 3, alpha = 0.25)
            if z==5.: ax.text(X_limits[1]-0.05, log10V+0.1, rf'$\rm {Survey}-{subSurvey}$', c='k', fontsize=7, ha = 'right')



for l, ls in zip([100, 250], ['--','-.']):

    log10V = 3*np.log10(l)
    ax.axhline(log10V, color = 'k', ls = '-.', lw = 1, alpha = 0.5)
    ax.text(X_limits[0]+0.05, log10V+0.1, rf'$\rm ({l}\ Mpc)^{{3}}$', c='k', fontsize=7)







# --- FLARES



flares = analyse.analyse('/Users/stephenwilkins/Dropbox/Research/data/simulations/FLARES/flares_no_particlesed.hdf5', default_tags = False)

# flares.list_datasets()

V = (4./3) * np.pi * (flares.radius)**3 # Mpc^3


for z, tag, c in zip(flares.zeds, flares.tags, cmr.take_cmap_colors('cmr.gem_r', len(flares.zeds))):

    ## ---- get data
    phi = np.zeros(len(bin_centres))
    N = np.zeros(len(bin_centres))
    V_total = 0.

    X = flares.load_dataset(tag, 'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI/', 'FUV')

    for i, (sim, w) in enumerate(zip(flares.sims, flares.weights)):

        x = np.log10(np.array(X[sim]))
        x = x[x>0.0]

        N_temp, _ = np.histogram(x, bins = bin_edges)

        N += N_temp

        phi_temp = (N_temp / V) / binw

        phi += phi_temp * w

    V_eff = N/(phi*binw)

    ax.plot(bin_centres, np.log10(V_eff), ls = '-', c=c, label = rf'$\rm z={z}$')


ax.legend(loc='upper left', fontsize=7, labelspacing=0.05, title=r'$\rm\bf FLARES-1$')

ax.set_xlim(X_limits)
ax.set_ylim(Y_limits)

ax.set_xlabel(r'$\rm \log_{10}(L_{FUV}/erg\ s^{-1}\ Hz^{-1})$')
ax.set_ylabel(r'$\rm\log_{10}(V_{eff}/cMpc^{3})$')

fig.savefig(f'figures/EffectiveVolume.pdf')
fig.clf()
