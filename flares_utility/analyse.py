
import numpy as np

import h5py

import astropy.constants as constants
import astropy.units as units

import os
this_dir, this_filename = os.path.split(__file__)


print('initialised cosmology')

from astropy.cosmology import WMAP9

cosmo = WMAP9




flares_master_file = os.environ['FLARES_MASTER']
# FLARES_MASTER=/cosma7/data/dp004/dc-payy1/my_files/flares_pipeline/data/flares.hdf5


# --- apertures and averaging timescales
apertures = [1, 3, 5, 10, 30, 50, 100]
averaging_timescales = [1,5,10,20,40,50,100,200]


# --- default scalings to the units specified in labels.py
scalings = {}
scalings['Mstar'] = 1E10
scalings['Mstar_total'] = 1E10
for aperture in apertures:
    scalings[f'Mstar_{aperture}'] = 1E10
scalings['BH_Mass'] = 1E10
scalings['age'] = 1000 # Gyr -> Myr
scalings['S_Mass'] = 1E10
scalings['S_MassInitial'] = 1E10
scalings['S_Age'] = 1000 # Gyr -> Myr

# converting MBHacc units to M_sol/yr
h = 0.6777  # Hubble parameter
BH_Mdot_scaling = h * 6.445909132449984E23  # g/s
BH_Mdot_scaling /= constants.M_sun.to('g').value  # convert to M_sol/s
BH_Mdot_scaling *= units.yr.to('s')  # convert to M_sol/yr
scalings['BH_Mdot'] = BH_Mdot_scaling




class analyse_flares:

    def __init__(self, fname, default_tags = True):

        self.fname =fname

        self.radius = 14./h # cMpc
        self.x, self.y, self.z, self.deltas, self.sigmas, self.weights = np.loadtxt(this_dir+'/data/weights_grid.txt', skiprows=1, unpack = True, usecols = (2,3,4,6,7,8), delimiter=',')

        self.ids = np.array(['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39'])
        self.sims = self.ids

        if default_tags:
            self.tags = np.array(['000_z015p000','001_z014p000','002_z013p000','003_z012p000','004_z011p000','005_z010p000','006_z009p000','007_z008p000','008_z007p000','009_z006p000','010_z005p000','011_z004p770'])
        else:
            with h5py.File(self.fname,'r') as f:
                self.tags = list(f['00'].keys())
                self.tags.pop() # get rid of header entry



        self.zeds = np.array([float(tag[5:].replace('p','.')) for tag in self.tags])

        self.zed_from_tag = {key:item for key, item in zip(self.tags, self.zeds)} # get redshift from tag
        self.tag_from_zed = {key:item for item, key in zip(self.tags, self.zeds)} # get tag from redshift

        self.apertures = apertures
        self.averaging_timescales = averaging_timescales


    def list_datasets(self):
        with h5py.File(self.fname,'r') as f:
            f[f'{self.sims[0]}/{self.tags[0]}'].visit(print)


    def load_single_dataset(self, sim, tag, path, dataset):

        with h5py.File(self.fname,'r') as f:
            try:
                k = f'{sim}/{tag}/{path}/{dataset}'
                out = f[k][:]
            except KeyError:
                out = []

        return out


    def load_dataset(self, tag, path, dataset):

        out = {}
        for sim in self.ids:
            out[sim] = self.load_single_dataset(sim, tag, path, dataset)
        return out


    def get_datasets(self, tag, quantities, apply_scalings = True, return_weights = True):


        # --- create dictionary of quantities

        Q = {}
        d = {}
        D = {}
        for q in quantities:

            if not q['name']:
                qid = q['dataset']
            else:
                qid = q['name']

            d[qid] = self.load_dataset(tag, q['path'], q['dataset'])
            D[qid] = np.array([])
            Q[qid] = q

        # --- read in weights
        if return_weights:
            D['weight'] = np.array([])

        D['delta'] = np.array([])

        for ii, sim in enumerate(self.ids):
            for qid in Q.keys():
                D[qid] = np.append(D[qid], d[qid][sim])
            if return_weights: D['weight'] = np.append(D['weight'], np.ones(np.shape(d[qid][sim]))*self.weights[ii])
            D['delta'] = np.append(D['delta'], np.ones(np.shape(d[qid][sim]))*self.deltas[ii])

        # --- apply standard scaling
        if apply_scalings:
            for qid, q in Q.items():
                if qid in scalings:
                    D[qid] *= scalings[qid]

        # --- create logged versions of the quantities
        for qid, q in Q.items():
            if q['log10']:
                D['log10'+qid] = np.log10(D[qid])

        return D


    def load_particles(self, sim, tag, q, return_dict = True):

        length_array = self.load_single_dataset(sim, tag, 'Galaxy', q[0]+'_Length')

        # find beginning:end indexes for each galaxy

        begin = np.zeros(len(length_array), dtype = np.int64)
        end = np.zeros(len(length_array), dtype = np.int64)
        begin[1:] = np.cumsum(length_array)[:-1]
        end = np.cumsum(length_array)


        d = self.load_single_dataset(sim, tag, 'Particle', q)


        if return_dict: # return a dictionary of arrays where the key is the galaxy index
            out = {}
            for i in np.arange(len(length_array)): # loop through gals
                out[i] =  d[begin[i]:end[i]]

        else: # otherwise return a list of array
            out = []
            for i in np.arange(len(length_array)): # loop through gals
                out.append(d[begin[i]:end[i]])

        # print(out)
        # print(len(out))
        return out


    def load_aperture_mask(self, sim, tag, particle_type = 'Star', aperture = '30', return_dict = True):

        length_array = self.load_single_dataset(sim, tag, 'Galaxy', particle_type[0]+'_Length')

        # find beginning:end indexes for each galaxy

        begin = np.zeros(len(length_array), dtype = np.int64)
        end = np.zeros(len(length_array), dtype = np.int64)
        begin[1:] = np.cumsum(length_array)[:-1]
        end = np.cumsum(length_array)


        d = self.load_single_dataset(sim, tag, 'Particle', f'Apertures/{particle_type}/{aperture}')

        if return_dict: # return a dictionary of arrays where the key is the galaxy index
            out = {}
            for i in np.arange(len(length_array)): # loop through gals
                out[i] =  d[begin[i]:end[i]]

        else: # otherwise return a list of array
            out = []
            for i in np.arange(len(length_array)): # loop through gals
                out.append(d[begin[i]:end[i]])

        # print(out)
        # print(len(out))
        return out




    def get_particle_datasets(self, tag, quantities = ['S_Age', 'S_Mass', 'S_MassInitial', 'S_Z'], aperture = '30', apply_scalings = True):

        P = {}

        for qid in quantities:
            P[qid] = []

        for ii, sim in enumerate(self.ids):

            if aperture:
                s = self.load_aperture_mask(sim, tag, aperture = aperture, return_dict = False)

            for qid in quantities:

                p = self.load_particles(sim, tag, qid, return_dict = False)

                if aperture:
                    for iii, p_ in enumerate(p):
                        p[iii] = p_[s[iii]]

                if apply_scalings and (qid in list(scalings.keys())):
                    for iii, p_ in enumerate(p):
                        p[iii] = p_*scalings[qid]

                P[qid] += p






        # D['pweight'] = []
        # for v,w in zip(D[ks[0]],D['weight']):
        #     D['pweight'].append(w*np.ones(v.shape))
        # D['pweight'] = np.array(D['pweight'])

        return P


analyse = analyse_flares
