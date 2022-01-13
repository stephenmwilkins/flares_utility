import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

import cmasher as cmr

from astropy.cosmology import Planck15 as cosmo

from scipy import interpolate


import flare.plt as fplt

# plt.style.use('simple')






fig = plt.figure(figsize = (3.5,4.5))

left  = 0.1
bottom = 0.1
height = 0.8
width = 0.85

ax = fig.add_axes((left, bottom, width, height))

cax = fig.add_axes((left, bottom+height, width, 0.02))


norm = mpl.colors.Normalize(vmin=0, vmax=10)
cmap = cmr.cosmic


for i in range(15):
    ax.plot([10, 2], [4+i, -3+i], lw = 1, c='k', alpha = 0.025, zorder = -1)



simulations = {}

# simulations['Technicolor Dawn'] = {'size': 12/0.7, 'DM_mass': 1.3725E6, 'resimulation': False,  'RT': True}
# simulations['CROC'] = {'size': 30/0.7, 'DM_mass': 7E6, 'resimulation': False,  'RT': True}
# simulations['CoDA'] = {'size': 91, 'DM_mass': 7E6, 'resimulation': False,  'RT': True}
# simulations['Renaissance'] = {'size': 8.3, 'DM_mass': 3E4, 'resimulation': 40,  'RT': True}
# simulations['Katz+17'] = {'size': 10/0.7, 'DM_mass': 6.5E6, 'resimulation': False,  'RT': True}
# simulations['SPHINX-5'] = {'size': 5/0.7, 'DM_mass': 3.1E4, 'resimulation': False,  'RT': True}
# simulations['SPHINX-10'] = {'size': 10/0.7, 'DM_mass': 2.5E5, 'resimulation': False,  'RT': True}



# simulations['Bahamas'] = {'zoom': False, 'size': 400/0.7, 'm_g': 8E8/0.7, 'resimulation': False,  'RT': False, 'redshift_end': 0.0}


simulations['EAGLE-Ref'] = {'zoom': False, 'size': 100, 'm_g': 1.81E6,  'RT': False, 'complete': True, 'redshift_end': 0.0, 'label': True}
simulations['EAGLE-Recal'] = {'zoom': False, 'size': 25, 'm_g': 2.26E5, 'RT': False, 'complete': True, 'redshift_end': 0.0, 'label': True}

simulations['Illustris-TNG50'] = {'zoom': False, 'size': 51.7, 'm_g': 8.5E4,  'RT': False, 'complete': True, 'redshift_end': 0.0, 'label': True}
simulations['Illustris-TNG100'] = {'zoom': False, 'size': 110.7, 'm_g': 1.4E6,  'RT': False, 'complete': True, 'redshift_end': 0.0, 'label': False}
simulations['Illustris-TNG300'] = {'zoom': False, 'size': 302.6, 'm_g': 1.1E7,  'RT': False, 'complete': True, 'redshift_end': 0.0, 'label': True}

simulations['Simba-100'] = {'zoom': False, 'size': 100/0.7, 'm_g': 1.82E7,   'RT': False, 'complete': True, 'redshift_end': 0.0, 'label': True}
simulations['Simba-50'] = {'zoom': False, 'size': 50/0.7, 'm_g': 2.28E6,  'RT': False, 'complete': True, 'redshift_end': 1.0, 'label': True}
simulations['Simba-25'] = {'zoom': False, 'size': 25/0.7, 'm_g': 2.85E5, 'RT': False, 'complete': True, 'redshift_end': 2.0, 'label': True}

simulations['Bluetides'] = {'zoom': False, 'size': 400/0.7, 'm_g': 2.36E6/0.7,  'RT': False, 'complete': True, 'redshift_end': 7.0, 'label': False}

simulations['ASTRID'] = {'zoom': False, 'size': 250/0.7, 'm_g': 1.27E6/0.7,  'RT': False, 'complete': True, 'redshift_end': 3.0, 'label': False}

simulations['FLARES-1'] = {'zoom': True, 'parent': 3200, 'size': [np.array([8.0, 8.5, 9.0, 9.5, 10., 10.5, 11.0]),np.array([6.5, 6.7, 6.8, 6.9, 7.1, 7.5, 8.0])], 'm_g': 1.81E6,  'RT': False, 'complete': True, 'redshift_end': 5.0, 'label': True}
simulations['FLARES-2-SD'] = {'zoom': True, 'parent': 6400, 'size': [np.array([8.0, 8.5, 9.0, 9.5, 10., 10.5, 11.0]),np.array([6.5, 6.7, 6.8, 6.9, 7.1, 7.5, 8.0])+1], 'm_g': 1E6,  'RT': False, 'complete': False, 'redshift_end': 5.0, 'label': True}
simulations['FLARES-2-HD'] = {'zoom': True, 'parent': 6400, 'size': [np.array([7.0, 8.0, 8.5, 9.0, 9.5, 10., 10.5, 11.0]),np.array([6.2, 6.5, 6.7, 6.8, 6.9, 7.1, 7.5, 8.0])], 'm_g': 1E5,  'RT': False, 'complete': False, 'redshift_end': 5.0, 'label': True}
simulations['COLIBRE-100'] = {'zoom': False, 'size': 100, 'm_g': 1E5,  'RT': False, 'complete': False, 'redshift_end': 0.0, 'label': True}
simulations['COLIBRE-250'] = {'zoom': False, 'size': 250, 'm_g': 1E6,  'RT': False, 'complete': False, 'redshift_end': 0.0, 'label': False}




j = 1
for i, (simulation_name, simulation) in enumerate(simulations.items()):

    s = simulation

    marker = 'o'

    if s['complete']:
        alpha = 1.0
    else:
        alpha = 0.3

    c = cmap(norm(s['redshift_end']))
    c = 'k'

    if s['zoom']:

        x = np.linspace(s['size'][0][0], s['size'][0][-1], 100)
        # y = np.interp(x, s['size'][0], s['size'][1])
        f = interpolate.interp1d(s['size'][0], s['size'][1], kind='cubic')
        y = f(x)

        norm_ = mpl.colors.Normalize(vmin=7, vmax=11)
        cmap_ = cmr.bubblegum

        c_ = cmap_(norm_(x))

        ax.plot([np.log10(s['m_g'])]*2, [y[0], 3*np.log10(s['parent'])-0.1],c='k',alpha=0.05, zorder = 0, lw = 5, solid_capstyle='round')
        ax.scatter([np.log10(s['m_g'])]*100, y, color=c_, s=10, zorder = 1)


        # ax.plot([np.log10(s['m_g'])]*2, np.log10(np.array(s['size'])), c=c, lw=1, zorder = 1, alpha = alpha)
        #

        label_loc = np.mean([y[0],y[-1]])
        label_loc = y[0]

        ax.text(np.log10(s['m_g'])+0.1, label_loc, rf'$\rm\bf {simulation_name}$', c='k', rotation = 90, fontsize = 8, ha = 'center', va = 'bottom', alpha = alpha)

    else:
        ax.scatter(np.log10(s['m_g']), 3*np.log10(s['size']), c=[c], s=20, lw=0, marker=marker, zorder = 2, alpha = alpha)

        # ax.text(np.log10(s['m_g'])-0.1, 3*np.log10(s['size'])-0.05, rf'$\rm {simulation_name}$', fontsize = 8, ha = 'left', va = 'center', alpha = alpha)
        # ax.text(np.log10(s['m_g'])-0.1, 3*np.log10(s['size']), simulation_name, fontsize = 8, ha = 'left', va = 'center', alpha = alpha)

        if s['label']:
            ax.text(np.log10(s['m_g']), 3*np.log10(s['size'])-0.2, rf'$\rm {simulation_name}$', fontsize = 7, ha = 'center', va = 'center', alpha = alpha, c=c)
        else:
            ax.text(np.log10(s['m_g'])-0.05, 3*np.log10(s['size'])-0.15, rf'$\rm {j}$', fontsize = 7, ha = 'left', va = 'center', alpha = alpha, c=c)
            ax.text(7.9, 5.-j*0.2, rf'$\rm {j}: {simulation_name}$', fontsize = 7, ha = 'left', va = 'center', alpha = alpha, c=c)
            j += 1

    #
    # ax.text(res, volume, i, fontsize = 5, ha = 'center', va = 'center', alpha = alpha)

    # ax.text(3.8, 9-i*0.4, rf'$\rm{{\bf{i}}}:{simulation_name}$', fontsize = 8, ha = 'left', va = 'center', alpha = alpha)
    # ax.text(3.8, 9-i*0.4, rf'$\rm{{\bf{i}}}:{simulation_name}$', fontsize = 8, ha = 'left', va = 'center', alpha = alpha)





add_future_surveys = True

if add_future_surveys:

    Surveys = {}

    # --- Euclid
    Surveys['Euclid'] = {}
    Surveys['Euclid']['Deep'] = {'area': 40*3600, 'Hlimit': 26.0, 'public': True}
    Surveys['Euclid']['Wide'] = {'area': 18000*3600, 'Hlimit': 24.0, 'public': True}

    # --- Webb
    Surveys['Webb'] = {}
    Surveys['Webb']['COSMOS-Web'] = {'area': 0.6*3600, 'Hlimit': 30.7, 'public': False}
    # Surveys['Webb']['NGDEEP'] = {'area': 8., 'Hlimit': 29.7, 'public': False}


    z = 7.

    v1, v2 = cosmo.comoving_volume([z-0.5, z+0.5]).to('Mpc3')

    for Survey in Surveys.keys():

        for subSurvey, SS in Surveys[Survey].items():

            volume = -np.log10(1. / ((v2 - v1) * (SS['area']/(41253.*3600))).value)


            print(volume)
            ax.axhline(volume, c='k', alpha=0.05, lw=10)

            ax.text(7.9, volume-0.025, fr'$\rm {{\bf {Survey}}}/{subSurvey}$', size=7, va = 'center', ha='left', color='k',zorder = 4, alpha=0.5)







ax.set_xlim([8., 4.])
ax.set_ylim([3.8, 11.5])

ax.set_xlabel(r'$\rm\log_{10}(m_{g}/M_{\odot})$')
ax.set_ylabel(r'$\rm\log_{10}(volume/cMpc^{3})$')


cmapper = cm.ScalarMappable(norm=norm_, cmap=cmap_)
cmapper.set_array([])
cbar = fig.colorbar(cmapper, cax=cax, orientation='horizontal')
cbar.set_label(r'$\rm log_{10}(M_{\star}/M_{\odot})$', fontsize = 7)
# cax.xaxis.tick_top()
# cax.xaxis.set_label_position(r'$\rm log_{10}(M_{\star}/M_{\odot})$')
cbar.ax.xaxis.set_ticks_position('top')
cbar.ax.xaxis.set_label_position('top')
cbar.ax.tick_params(labelsize=6)

fig.savefig(f'figures/simulation_comparison.pdf')


fig.clf()
