
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

import cmasher as cmr

import flares
import flares_analysis as fa

import h5py


log10L_limits = [27, 32]

fig = plt.figure(figsize = (3.5,3.5))

left  = 0.15
bottom = 0.15
height = 0.8
width = 0.8

ax = fig.add_axes((left, bottom, width, height))


redshifts = [5,6,7,8,9,10]








binw = 0.1
bin_edges = np.arange(*log10L_limits, binw)
bin_centres = bin_edges[:-1]+binw/0.5



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


    ax.plot(bin_centres, np.log10(phi), ls = '-', c=c)





ax.legend(labelspacing=0.05,fontsize=8)
ax.set_xlim(log10L_limits)
# ax.set_ylim([-12., -1.51])


ax.set_xlabel(r'$\rm \log_{10}(L_{FUV}/erg\ s^{-1}\ Hz^{-1})$')
ax.set_ylabel(r'$\rm\log_{10}[\phi/Mpc^{-3}\ dex^{-1}]$')


fig.savefig(f'figures/LF.pdf')


fig.clf()
