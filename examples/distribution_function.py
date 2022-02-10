
import sys
import os

import numpy as np

import matplotlib.pyplot as plt

# --- this is necessary if the module isn't "installed" or in the PythonPath
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import flares_utility.analyse as analyse



redshift = 11



left  = 0.25
bottom = 0.25
height = 0.7
width = 0.7


# --- initialise a basic plot
fig_N = plt.figure(figsize = (3.5,3.5))
ax_N = fig_N.add_axes((left, bottom, width, height))

# --- initialise a SECOND basic plot
fig_phi = plt.figure(figsize = (3.5,3.5))
ax_phi = fig_phi.add_axes((left, bottom, width, height))




# --- open FLARES data
filename = analyse.flares_master_file+'/flares_highz_v3_nosed.hdf5' # if you don't have an environment variable set you need to give this the full path to the master file instead
flares = analyse.analyse(filename, default_tags = False)
tag = flares.tag_from_zed[redshift] # determine the tag of the redshift we're interested in

# --- define the volume of each FLARES simulation
V = (4./3) * np.pi * (flares.radius)**3 # Mpc^3



# --- define the bins we're going to use
X_limits = [8., 11] # x-range to consider
binw = 0.25
bin_edges = np.arange(*X_limits, binw)
bin_centres = bin_edges[:-1]+binw/2

# --- make empty arrays
phi = np.zeros(len(bin_centres)) # \phi - the number of galaxies per unit volume
N = np.zeros(len(bin_centres)) # N - the number of galaxies

# --- select relevant dataset
# X = flares.load_dataset(tag, 'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI/', 'FUV') # --- far-UV luminosity
X = flares.load_dataset(tag, 'Galaxy', 'Mstar') # -- stellar mass

# --- loop over simulations
for i, (sim, w) in enumerate(zip(flares.sims, flares.weights)):

    x = np.log10(np.array(X[sim])) + 10
    x = x[x>0.0] # only select objects with non-zero L, mass etc.

    N_temp, _ = np.histogram(x, bins = bin_edges) # make temporary histogram

    N += N_temp # add tempoary histogram to number count

    phi += (N_temp / V) / binw  * w # add tempoary histogram to running \phi accounting for volume, bin_width, and simulation weight


ax_N.plot(bin_centres, np.log10(N)) # plot
ax_phi.plot(bin_centres, np.log10(phi)) # plot

ax_N.set_xlim(X_limits)
ax_phi.set_xlim(X_limits)

ax_N.set_ylabel(r'$\rm\log_{10}N$')
ax_phi.set_ylabel(r'$\rm\log_{10}[\phi/Mpc^{-3}\ dex^{-1}]$')

for ax in [ax_N, ax_phi]:
    ax.set_xlabel(r'$\rm \log_{10}(M_{\star}/M_{\odot})$')

fig_N.savefig(f'figs/distribution_function_N.pdf')
fig_phi.savefig(f'figs/distribution_function_phi.pdf')
