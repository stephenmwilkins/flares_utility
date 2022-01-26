import h5py
import flares


# this strips the particle data and SED data. Probably overkill.

# tags = fl.tags
tags = ['000_z015p000', '001_z014p000', '002_z013p000', '003_z012p000', '004_z011p000', '005_z010p000']

# in_dir = '/cosma7/data/dp004/dc-payy1/my_files/flares_pipeline/data/'
in_dir = '/cosma/home/dp004/dc-wilk2/data/flare/simulations/flares/'

# filename = 'flares'
filename = 'flares_highz_v3'

# out_dir
out_dir = in_dir

fl = flares.flares(in_dir + filename +'.hdf5', sim_type='FLARES')

fs = h5py.File(in_dir + filename +'.hdf5', 'r')
fd = h5py.File(out_dir + filename + '_nosed.hdf5', 'w')



for sim in fl.halos:
    fd.create_group(f'{sim}')
    for tag in tags:
        print(sim, tag)
        fs.copy(f'{sim}/{tag}', fd[f'{sim}/'])
        del fd[f'{sim}/{tag}/Galaxy/BPASS_2.2.1/Chabrier300/SED']


fd.flush()
