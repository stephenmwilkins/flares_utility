
import numpy as np

from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

from astropy.cosmology import Planck15 as cosmo

import FLARE.plt as fplt
import FLARE.photom as phot








fig = plt.figure(figsize = (3.5,3.5))

left  = 0.15
bottom = 0.15
height = 0.8
width = 0.8

ax = fig.add_axes((left, bottom, width, height))


for i,z in enumerate([5.0, 6.0, 7.0, 8.0, 9.0, 10.0]):

    c = cm.viridis(i/6)

    for simID, ls, alpha, label in zip(['ref','flares'],['--','-'], [0.1, 0.2],['EAGLE\ REF', 'FLARES']):

        sim = ascii.read(f"data/lf_{simID}.csv", format='csv', fast_reader=False)

        M = sim[f'M{z}']
        phi = sim[f'phi{z}']
        err = sim[f'err{z}']

        if simID == 'ref':
            binw = M[1] - M[0]
            N = phi*binw*100**3
        else:
            N = sim[f'num{z}']



        V_eff = (N/phi)

        # ax.fill_between(M[s], np.log10(phi[s]-err[s]), np.log10(phi[s]+err[s]), color=c, alpha=alpha)

        if simID == 'flares':
            ax.plot(M, np.log10(V_eff), c=c, ls=ls, lw=1, label = rf'$\rm z={z}$')
        # else:
        #     ax.plot(M, np.log10(V_eff), c=c, ls=ls, lw=1)



# ax.axhline(3*np.log10(6.4E3), color = 'k', ls = '--', lw = 1, alpha = 0.5, label = r'$\rm (6.4\ Gpc)^3$')
#
# ax.axhline(3*np.log10(3.2E3), color = 'k', ls = ':', lw = 1, alpha = 0.5, label = r'$\rm (3.2\ Gpc)^3$')

# ax.axhline(3*np.log10(400/0.7), color = 'k', ls = '-.', lw = 1, alpha = 0.5, label = r'$\rm (571\ Mpc)^3$')

# ax.axhline(3*np.log10(300), color = 'k', ls = '-.', lw = 1, alpha = 0.5, label = r'$\rm (300\ Mpc)^3$')

ax.axhline(np.log10(40*(4/3)*np.pi*(14/0.7)**3), color = 'k', ls = '-', lw = 3, alpha = 0.25, label = r'$\rm FLARES\ total\ simulated$')

ax.axhline(3*np.log10(100), color = 'k', ls = '-', lw = 1, alpha = 0.5, label = r'$\rm EAGLE\ REF\ (100\ Mpc)^3$')




ax.legend(loc='upper left', fontsize=7)


ax.set_xlim([-17.5, -25.])
ax.set_ylim([5., 10.])


ax.set_xlabel(r'$\rm M_{FUV}$')
ax.set_ylabel(r'$\rm\log_{10}(V_{eff}/cMpc^{3})$')


fig.savefig(f'figures/effective_volume.pdf')
fig.clf()
