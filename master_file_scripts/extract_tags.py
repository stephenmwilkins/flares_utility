

import h5py
# import flares_utility


import flares





# in_file = '/cosma7/data/dp004/dc-payy1/my_files/flares_pipeline/data/flares.hdf5' # Aswin's version
in_file = '/cosma7/data/dp004/dc-love2/codes/flares/data/flares.hdf5' # Chris's version
in_file = '/cosma7/data/dp004/dc-payy1/my_files/flares_pipeline/data/flares.hdf5' #Aswin's newer version

out_file = '/cosma/home/dp004/dc-wilk2/data/flare/simulations/flares/flares_highz_v3.hdf5'

fl = flares.flares(in_file, sim_type='FLARES')


fin = h5py.File(in_file, 'r')

# tags = fin['00'].keys()
tags = ['000_z015p000', '001_z014p000', '002_z013p000', '003_z012p000', '004_z011p000', '005_z010p000']

fout = h5py.File(out_file, 'w')

for sim in fl.halos:
    fout.create_group(f'{sim}')
    for tag in tags:
        print(sim, tag)
        fout.create_group(f'{sim}/{tag}')
        fin.copy(f'{sim}/{tag}/Galaxy', fout[f'{sim}/{tag}'])
        fin.copy(f'{sim}/{tag}/Particle', fout[f'{sim}/{tag}'])
        # del fd[f'{sim}/{tag}/Galaxy/BPASS_2.2.1/Chabrier300/SED']

fout.flush()
