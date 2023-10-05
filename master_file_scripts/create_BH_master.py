import h5py


# this strips the full SED and the non-BH particle data

indir = '/cosma7old/data/dp004/dc-payy1/my_files/flares_pipeline/data/'
# indir = '/cosma7/data/dp004/dc-seey1/data/flares/steve/'

outdir = '/cosma/home/dp004/dc-wilk2/scratch/'



fs = h5py.File(indir + 'flares.hdf5', 'r')
fd = h5py.File(outdir + 'flares_bh.hdf5', 'w')

# testing
sims = ['00']
tags = ['008_z007p000']

# sims = []
# tags = ['000_z015p000', '001_z014p000', '002_z013p000', '003_z012p000', '004_z011p000', '005_z010p000']

for sim in sims:
    fd.create_group(f'{sim}')
    for tag in tags:
        print(sim, tag)
        fd.create_group(f'{sim}/{tag}')
        fs.copy(f'{sim}/{tag}/Galaxy', fd[f'{sim}/{tag}'])
        del fd[f'{sim}/{tag}/Galaxy/BPASS_2.2.1/Chabrier300/SED']

        fd.create_group(f'{sim}/{tag}/Particle')
        print(len(fd['{sim}/{tag}/Particle/BH_ID']))
        for k in ['BH_Age', 'BH_Coordinates', 'BH_ID', 'BH_Index', 'BH_Mass', 'BH_Mdot', 'BH_los', 'BH_sml']:
            fs.copy(f'{sim}/{tag}/Particle/{k}', fd[f'{sim}/{tag}/Particle'])


fd.flush()
