
import numpy as np

from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

from astropy.cosmology import Planck15 as cosmo

import FLARE.plt as fplt
import FLARE.photom as phot

# plt.style.use('simple')






fig = plt.figure(figsize = (3.5,3.5))

left  = 0.15
bottom = 0.15
height = 0.6
width = 0.8

ax = fig.add_axes((left, bottom, width, height))

axH = fig.add_axes((left, bottom+height, width, 0.2))

for z in [10., 6.]:

    for simID, ls, alpha, c, label in zip(['ref','flares'],['--','-'], [0.1, 0.2], ['0.5','k'],['EAGLE\ REF', 'FLARES']):

        sim = ascii.read(f"data/lf_{simID}.csv", format='csv', fast_reader=False)

        M = sim[f'M{z}']
        phi = sim[f'phi{z}']
        err = sim[f'err{z}']

        if simID == 'ref':
            binw = M[1] - M[0]
            N = phi*binw*100**3
        else:
            N = sim[f'num{z}']

        s = err<phi
        s = phi>0

        print((N/phi)**(1/3))

        # ax.fill_between(M[s], np.log10(phi[s]-err[s]), np.log10(phi[s]+err[s]), color=c, alpha=alpha)
        ax.plot(M[s], np.log10(phi[s]), c=c, ls=ls, lw=1)

        axH.plot(M, np.log10(N), ls=ls, lw=1, c=c, label = rf'$\rm {label}$')

    # --- surveys

    Surveys = {}

    # --- Euclid
    Surveys['Euclid'] = {}
    Surveys['Euclid']['Deep'] = {'area': 40*3600, 'Hlimit': 26.0, 'public': True}
    Surveys['Euclid']['Wide'] = {'area': 18000*3600, 'Hlimit': 24.0, 'public': True}

    # --- Roman
    Surveys['Roman'] = {}
    Surveys['Roman']['HLS'] = {'area': 2000*3600, 'Hlimit': 26.9, 'public': True}
    # Surveys['Roman']['SN-Wide'] = {'area': 18.5*3600, 'Hlimit': 28.7, 'public': True}
    # Surveys['Roman']['SN-Deep'] = {'area': 8.5*3600, 'Hlimit': 29.6, 'public': True}


    # --- Webb
    # Surveys['Webb'] = {}
    # Surveys['Webb']['JADES-Deep'] = {'area': 46., 'Hlimit': 30.7, 'public': False}
    # Surveys['Webb']['JADES-Medium'] = {'area': 144., 'Hlimit': 29.7, 'public': False}

    v1, v2 = cosmo.comoving_volume([z-0.5, z+0.5]).to('Mpc3')

    for Survey in Surveys.keys():

        for subSurvey, SS in Surveys[Survey].items():


            flux = phot.m_to_flux(SS['Hlimit'])
            lum = phot.flux_to_L(flux, cosmo, z)
            M = phot.lum_to_M(lum)



            low_lim = np.log10(1. / ((v2 - v1) * (SS['area']/(41253.*3600))).value)

            print(low_lim)

            # ax.axhline(-low_lim, c='k', alpha=0.1, lw=10)
            # ax.plot([10, mass_limit-1.0], [low_lim]*2, c='k', alpha=0.1, lw=10, zorder = 3)

            ax.fill_between([M, -30], [low_lim, low_lim], alpha=0.1, label = fr'$\rm {{\bf {Survey} }}/{subSurvey}$')




ax.axhline(-3*np.log10(3.2E3), color = 'k', ls = ':', lw = 1, alpha = 0.5, label = r'$\rm (3.2\ Gpc)^3$')
ax.axhline(-3*np.log10(6.4E3), color = 'k', ls = '-.', lw = 1, alpha = 0.5, label = r'$\rm (6.4\ Gpc)^3$')

ax.legend(loc='center left', fontsize=8)
axH.legend(loc='upper right', fontsize=8)

ax.set_xlim([-17.5, -25.])
ax.set_ylim([-12., -2.01])

axH.set_xlim([-17.5, -25.])
axH.set_ylim([-0.1, 4.])
axH.set_xticks([])


ax.set_xlabel(r'$\rm M_{FUV}$')
ax.set_ylabel(r'$\rm\log_{10}(\phi/cMpc^{-3}\, mag^{-1})$')

axH.set_ylabel(r'$\rm\log_{10}(N_{i})$')


ax.text(-23.5, -3, fr'$\rm z={z}$', fontsize = 12)


fig.savefig(f'figures/LF_2.pdf')
fig.clf()
