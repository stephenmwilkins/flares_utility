


# --- extract a single simulation and snapshot. Useful for testing or giving to students.


filedir = '/cosma7/data/dp004/dc-payy1/my_files/flares_pipeline/data/'

import h5py

sim = '00'
snap = '008_z007p000'

fs = h5py.File(filedir + 'flares.hdf5', 'r')
fd = h5py.File(f'flares_{sim}_{snap}.hdf5', 'w')
fs.copy(f'{sim}/{snap}', fd)
fd.flush()
