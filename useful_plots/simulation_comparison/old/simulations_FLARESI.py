import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

from astropy.cosmology import Planck15 as cosmo

# plt.style.use('simple')

import FLARE.plt as fplt


fig = plt.figure(figsize = (3.5,3.5))

left  = 0.1
bottom = 0.1
height = 0.75
width = 0.75

ax = fig.add_axes((left, bottom, width, height))

# axN = ax.twinx()


include_FLARES = True
add_current_surveys = False
add_future_surveys = False



for i in range(15):
    ax.plot([10, 2], [4+i, -3+i], lw = 1, c='k', alpha = 0.025, zorder = -1)



simulations = {}

simulations['Technicolor Dawn'] = {'size': 12/0.7, 'DM_mass': 1.3725E6, 'resimulation': False,  'RT': True}
simulations['GIMIC'] = {'size': (4*(4/3)*np.pi*(18/0.7)**3 + (4/3)*np.pi*(25/0.7)**3)**(1/3), 'DM_mass': 9.3E6, 'resimulation': 500/0.7,  'RT': False}
simulations['EAGLE-Ref'] = {'size': 100, 'DM_mass': 9.7E6, 'resimulation': False,  'RT': False}
simulations['EAGLE-Recal'] = {'size': 25, 'DM_mass': 9.7E5, 'resimulation': False,  'RT': False}
simulations['CROC'] = {'size': 30/0.7, 'DM_mass': 7E6, 'resimulation': False,  'RT': True}
simulations['CoDA'] = {'size': 91, 'DM_mass': 7E6, 'resimulation': False,  'RT': True}
simulations['Illustris-TNG100'] = {'size': 110.7, 'DM_mass': 5E6, 'resimulation': False,  'RT': False}
simulations['Illustris-TNG300'] = {'size': 302.6, 'DM_mass': 3.98E7/0.7, 'resimulation': False,  'RT': False}
simulations['Renaissance'] = {'size': 8.3, 'DM_mass': 3E4, 'resimulation': 40,  'RT': True}
simulations['Katz+17'] = {'size': 10/0.7, 'DM_mass': 6.5E6, 'resimulation': False,  'RT': True}
simulations['SPHINX-5'] = {'size': 5/0.7, 'DM_mass': 3.1E4, 'resimulation': False,  'RT': True}
simulations['SPHINX-10'] = {'size': 10/0.7, 'DM_mass': 2.5E5, 'resimulation': False,  'RT': True}
simulations['SIMBA-100'] = {'size': 100/0.7, 'DM_mass': 9.6E7, 'resimulation': False,  'RT': False}
simulations['SIMBA-50'] = {'size': 50/0.7, 'DM_mass': 1.2E7, 'resimulation': False,  'RT': False}
simulations['SIMBA-25'] = {'size': 25/0.7, 'DM_mass': 1.5E6, 'resimulation': False,  'RT': False}
simulations['Bluetides'] = {'size': 400/0.7, 'DM_mass': 1.7E7, 'resimulation': False, 'RT': False}

simulations['FLARES-1'] = {'size': (40*(4/3)*np.pi*(14/0.7)**3)**(1/3), 'DM_mass': 9.7E6, 'resimulation': 3200,  'RT': False}
# simulations['FLARES-HR'] = {'size': (10*(4/3)*np.pi*14**3)**(1/3), 'DM_mass': 9.7E5, 'resimulation': 3200,  'RT': False}



all_simulations = list(simulations.keys())

my_simulations = ['FLARES-1']

other_simulations = [x for x in all_simulations if x not in my_simulations]

c = lambda x: cm.RdYlBu((x-8.5)/3.5)

i = 1
for simulation_name in other_simulations:

    simulation = simulations[simulation_name]
    res = np.log10(simulation['DM_mass'])
    volume = 3.*np.log10(simulation['size'])
    elements = np.log10(2*(7000)**3) - res + volume + 7.23 - 8.27 # approximate based on Bluetides scaling
    print(simulation_name, res, volume, elements)

    if simulation['RT']:
        marker = 's'
    else:
        marker = 'o'

    color = c(elements)
    # color = '0.9'

    if include_FLARES:
        alpha = 0.3
    else:
        alpha = 1.0


    if simulation['resimulation']:
        parent_volume = 3.*np.log10(simulation['resimulation'])
        ax.plot([res]*2, [volume, parent_volume], c='0.8', lw=3, zorder = 1, alpha = alpha)
        ax.scatter(res, parent_volume, s=30, marker=marker, c='0.8', zorder = 1, alpha = alpha)


    ax.scatter(res, volume, c=[color], s=100, marker=marker, edgecolors='k',lw=1, zorder = 2, alpha = alpha)

    ax.text(res, volume, i, fontsize = 5, ha = 'center', va = 'center', alpha = alpha)
    ax.text(9.8, 10-i*0.4, rf'$\rm{{\bf{i}}}:{simulation_name}$', fontsize = 8, ha = 'left', va = 'center', alpha = 0.8)

    i += 1


# --- add surveys

if add_current_surveys:

    Surveys = {}

    Surveys['VISTA'] = {}
    Surveys['VISTA']['UltraVISTA'] = {'area': 1.5*3600, 'Hlimit': 25.7, 'public': True}
    Surveys['VISTA']['VIDEO'] = {'area': 12*3600, 'Hlimit': 24.0, 'public': True}

    Surveys['Hubble'] = {}
    Surveys['Hubble']['XDF'] = {'area': 4.7, 'Hlimit': 29.4, 'public': True}
    Surveys['Hubble']['CANDELS-DEEP'] = {'area': 64.5, 'Hlimit': 27.5, 'public': True}
    Surveys['Hubble']['CANDELS-WIDE'] = {'area': 34.2+151.2+151.9+150.7, 'Hlimit': 26.8, 'public': True}

    z = 8.

    v1, v2 = cosmo.comoving_volume([z-0.5, z+0.5]).to('Mpc3')

    for Survey in Surveys.keys():

        for subSurvey, SS in Surveys[Survey].items():

            mass_limit = 9.0-0.4*(SS['Hlimit']-26.)

            low_lim = -np.log10(1. / ((v2 - v1) * (SS['area']/(41253.*3600))).value)

            print(low_lim)

            a1 = 0.1
            a2 = 1.0

            if add_future_surveys:
                a1 = 0.02
                a2 = 0.2

            # ax.axhline(-low_lim, c='k', alpha=0.1, lw=10)
            ax.plot([10, mass_limit-1.0], [low_lim]*2, c='k', alpha=a1, lw=10, zorder = 3)

            ax.text(9.8, low_lim-0.05, fr'$\rm {{\bf {Survey}}}/{subSurvey}$', size=8, va = 'center', ha='left', color='k',zorder = 4, alpha=a2)



if add_future_surveys:

    Surveys = {}

    # --- Euclid
    Surveys['Euclid'] = {}
    Surveys['Euclid']['Deep'] = {'area': 40*3600, 'Hlimit': 26.0, 'public': True}
    Surveys['Euclid']['Wide'] = {'area': 18000*3600, 'Hlimit': 24.0, 'public': True}

    # --- Webb
    Surveys['Webb'] = {}
    Surveys['Webb']['JADES-Deep'] = {'area': 46., 'Hlimit': 30.7, 'public': False}
    Surveys['Webb']['JADES-Medium'] = {'area': 144., 'Hlimit': 29.7, 'public': False}


    z = 8.

    v1, v2 = cosmo.comoving_volume([z-0.5, z+0.5]).to('Mpc3')

    for Survey in Surveys.keys():

        for subSurvey, SS in Surveys[Survey].items():

            mass_limit = 9.0-0.4*(SS['Hlimit']-26.)

            low_lim = -np.log10(1. / ((v2 - v1) * (SS['area']/(41253.*3600))).value)

            print(low_lim)

            # ax.axhline(-low_lim, c='k', alpha=0.1, lw=10)
            ax.plot([10, mass_limit-1.0], [low_lim]*2, c='k', alpha=0.1, lw=10, zorder = 3)

            ax.text(9.8, low_lim-0.05, fr'$\rm {{\bf {Survey}}}/{subSurvey}$', size=8, va = 'center', ha='left', color='k',zorder = 4)



# --- my simulations

if include_FLARES:
    for simulation_name in my_simulations:

        simulation = simulations[simulation_name]
        res = np.log10(simulation['DM_mass'])
        volume = 3.*np.log10(simulation['size'])
        elements = np.log10(2*(7000)**3) - res + volume + 7.23 - 8.27 # approximate based on Bluetides scaling
        print(simulation_name, res, volume, elements)

        marker = 'o'
        color = c(10.38) # FLARES-I number
        # color = '0.9'
        alpha = 1.0

        if simulation['resimulation']:
            parent_volume = 3.*np.log10(simulation['resimulation'])
            ax.plot([res]*2, [volume, parent_volume], c='0.8', lw=3, alpha = alpha,zorder = 5)
            ax.scatter(res, parent_volume, s=30, marker=marker, c='0.8', alpha = alpha,zorder = 6)

        ax.scatter(res, volume, c=[color], s=100, marker=marker, edgecolors='k',lw=1, alpha = alpha,zorder = 6)

        ax.text(res-0.2, volume-0.05, rf'$\rm\bf {simulation_name}$', fontsize = 9, ha = 'left', va = 'center', alpha = alpha,zorder = 6)













ax.set_xlim([10., 4.])
ax.set_ylim([2., 11.5])

ax.set_xlabel(r'$\rm\log_{10}(m_{DM}/M_{\odot})$')
ax.set_ylabel(r'$\rm\log_{10}(volume/cMpc^{3})$')

def forward(x):
    return x + 1.0

def inverse(x):
    return x - 1.0

secax = ax.secondary_xaxis('top', functions=(forward, inverse))
secax.set_xlabel(r'$\rm\log_{10}(M_{\star}/M_{\odot})$')


ax1 = fig.add_axes([0.85, 0.1, 0.03, 0.75])

cmap = mpl.cm.RdYlBu
norm = mpl.colors.Normalize(vmin=8.5, vmax=12.0)
cb = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='vertical')

ax1.set_ylabel(r'$\rm\log_{10}(resolution\ elements)$')

ax1.tick_params(axis='both', which='major', labelsize=4)


# ax.legend(loc='upper right', fontsize=6, labelspacing=0.4)

fig.savefig(f'figures/simulations_FLARESI.pdf')

fig.clf()
