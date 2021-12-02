import h5py
import flares


# this strips the particle data and SED data. Probably overkill.



filedir = '/cosma7/data/dp004/dc-payy1/my_files/flares_pipeline/data/'

fl = flares.flares(filedir + 'flares.hdf5', sim_type='FLARES')

fs = h5py.File(filedir + 'flares.hdf5', 'r')
fd = h5py.File('flares_nosed.hdf5', 'w')

for sim in fl.halos:
    fd.create_group(f'{sim}')
    for tag in fl.tags:
        print(sim, tag)
        fs.copy(f'{sim}/{tag}', fd[f'{sim}/'])
        del fd[f'{sim}/{tag}/Galaxy/BPASS_2.2.1/Chabrier300/SED']


fd.flush()
