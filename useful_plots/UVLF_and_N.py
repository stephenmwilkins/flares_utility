
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

import cmasher as cmr

import h5py

import flare

cosmo = flare.default_cosmo()

import flare.plt as fplt
import flare.photom as phot

import flares_utility.analyse as analyse

import flare.surveys

X_limits = [28, 30.5]

fig = plt.figure(figsize = (3.5,3.5))

left  = 0.15
bottom = 0.15
height = 0.6
width = 0.8

ax = fig.add_axes((left, bottom, width, height))

axN = fig.add_axes((left, bottom+height, width, 0.2))

z = 5
c = 'k'
include_EAGLE = True
include_FLARES = True
include_surveys = False
include_surveys = ['Euclid']




binw = 0.1
bin_edges = np.arange(*X_limits, binw)
bin_centres = bin_edges[:-1]+binw/0.5




if include_surveys:

    # --- surveys

    v1, v2 = cosmo.comoving_volume([z-0.5, z+0.5]).to('Mpc3')

    i = -1
    imax = 5

    for survey_name in include_surveys:

        survey = flare.surveys.surveys[survey_name]

        for field_name, field in survey.fields.items():

            i += 1
            cs = cm.plasma(i/imax)

            print(field)

            flux = phot.m_to_flux(field.depths_mag[field.depth_reference_filter])
            lum = phot.flux_to_L(flux, cosmo, z)
            # M = phot.lum_to_M(lum)

            low_lim = np.log10(1. / ((v2 - v1) * (field.area/(41253.*3600))).value)

            # ax.axhline(-low_lim, c='k', alpha=0.1, lw=10)
            # ax.plot([10, mass_limit-1.0], [low_lim]*2, c='k', alpha=0.1, lw=10, zorder = 3)

            ax.fill_between([np.log10(lum), 50], [low_lim, low_lim], [10,10], alpha=0.1, label = fr'$\rm {{\bf {survey_name} }}/{field_name}$', color=cs)







# --- EAGLE REF

if include_EAGLE:

    EAGLE = h5py.File('/Users/stephenwilkins/Dropbox/Research/data/simulations/FLARES/EAGLE_REF_sp_info.hdf5', 'r')
    snap = {5:'008_z005p037',6: '006_z005p971',7:'005_z007p050',8:'004_z008p075',9:'003_z008p988',10:'002_z009p993'}
    V = 100**3

    X = np.log10(EAGLE[snap[z]]['Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI/FUV'][:])

    N, _ = np.histogram(X, bins = bin_edges)

    phi = (N / V) / binw

    ax.plot(bin_centres, np.log10(phi), ls = '--', c=c, alpha = 0.8)

    # --- plot the number of galaxies in each bin
    # axN.plot(bin_centres, np.log10(N), ls = '--', c=c, alpha = 0.8, label = rf'$\rm EAGLE $')

    # --- plot the cumulative number of galaxies in each bin
    C = []
    for bc in bin_centres:
        C.append(np.log10(len(X[X>bc])))
    axN.plot(bin_centres, C, ls = '--', c=c, label = rf'$\rm EAGLE $')



# --- FLARES

if include_FLARES:

    flares = analyse.analyse('/Users/stephenwilkins/Dropbox/Research/data/simulations/FLARES/flares_no_particlesed.hdf5', default_tags = False)

    # flares.list_datasets()

    V = (4./3) * np.pi * (flares.radius)**3 # Mpc^3

    tag = flares.tag_from_zed[z]

    ## ---- get data
    phi = np.zeros(len(bin_centres))
    N = np.zeros(len(bin_centres))

    X = flares.load_dataset(tag, 'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI/', 'FUV')

    x_ = np.array([])

    for i, (sim, w) in enumerate(zip(flares.sims, flares.weights)):

        x = np.log10(np.array(X[sim]))
        x = x[x>0.0]

        N_temp, _ = np.histogram(x, bins = bin_edges)

        N += N_temp

        phi_temp = (N_temp / V) / binw

        phi += phi_temp * w

        x_ = np.append(x_, x)


    ax.plot(bin_centres, np.log10(phi), ls = '-', c=c, alpha = 0.8)

    # --- plot the number of galaxies in each bin
    # axN.plot(bin_centres, np.log10(N), ls = '-', c=c, alpha = 0.8, label = rf'$\rm FLARES $')

    # --- plot the cumulative number of galaxies in each bin
    C = []
    for bc in bin_centres:
        C.append(np.log10(len(x_[x_>bc])))
    axN.plot(bin_centres, C, ls = '-', c=c, label = rf'$\rm FLARES $')




ax.legend(loc='lower left', labelspacing=0.05,fontsize=8)
axN.legend(labelspacing=0.05,fontsize=8)

axN.set_xticks([])

ax.set_xlim(X_limits)
ax.set_ylim(top = np.max(np.log10(phi))+0.7)
axN.set_xlim(X_limits)

ax.set_xlabel(r'$\rm \log_{10}(L_{FUV}/erg\ s^{-1}\ Hz^{-1})$')
ax.set_ylabel(r'$\rm\log_{10}[\phi/Mpc^{-3}\ dex^{-1}]$')


ax.text(0.85, 0.9, rf'$\rm z={z}$', fontsize = 14, transform=ax.transAxes)

axN.set_ylabel(r'$\rm\log_{10}[N(>L)]$')


fig.savefig(f'figures/UVLF_and_N.pdf')


fig.clf()
