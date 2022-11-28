import h5py
import flares


# this strips the particle data and SED data. Probably overkill.


# indir = '/cosma7/data/dp004/dc-payy1/my_files/flares_pipeline/data/'
indir = '/cosma7/data/dp004/dc-seey1/data/flares/steve/'

outdir = '/cosma/home/dp004/dc-wilk2/data/flare/simulations/flares/'


fl = flares.flares(indir + 'flares.hdf5', sim_type='FLARES')

fs = h5py.File(indir + 'flares.hdf5', 'r')
fd = h5py.File(outdir + 'flares_noparticlesed_v3.hdf5', 'w')

for sim in fl.halos:
    fd.create_group(f'{sim}')
    for tag in fl.tags:
        print(sim, tag)
        fd.create_group(f'{sim}/{tag}')
        fs.copy(f'{sim}/{tag}/Galaxy', fd[f'{sim}/{tag}'])
        del fd[f'{sim}/{tag}/Galaxy/BPASS_2.2.1/Chabrier300/SED']


fd.flush()
