
import numpy as np

from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

import cmasher as cmr

from astropy.cosmology import Planck15 as cosmo

import flare.plt as fplt
import flare.photom as phot

import flares
import flares_analysis as fa

import h5py



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

fl = flares.flares('/Users/stephenwilkins/Dropbox/Research/data/simulations/FLARES/flares_no_particlesed.hdf5', sim_type='FLARES')

snap = {5:'010_z005p000',6: '009_z006p000',7:'008_z007p000',8:'007_z008p000',9:'006_z009p000',10:'005_z010p000'}

dat = np.loadtxt('/Users/stephenwilkins/Dropbox/Research/modules/flares/weight_files/weights_grid.txt', skiprows=1, delimiter=',')
weights = dat[:,8]
index = dat[:,0]

X = fl.load_dataset('FUV',arr_type='Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI')

R = 14./0.6777 # Mpc
V = (4./3) * np.pi * (R)**3 # Mpc^3

for z, c in zip(redshifts, cmr.take_cmap_colors('cmr.gem_r', len(redshifts))):

    ## ---- get data
    phi = np.zeros(len(bin_centres))
    N = np.zeros(len(bin_centres))
    V_total = 0.

    for i,halo in enumerate(fl.halos):

        w = weights[np.where(["%02d"%i == halo for i in index])[0]]

        x = X[halo][snap[z]]
        x = x[x>0.0]

        V_total += V

        N_temp, _ = np.histogram(np.log10(x), bins = bin_edges)

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
