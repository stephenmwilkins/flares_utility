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

sims = list(fs.keys())
print(sims)
tags = list(fs[sims[0]].keys())
print(tags)

for sim in sims:
    fd.create_group(f'{sim}')
    for tag in tags:
        print(sim, tag)

        fd.create_group(f'{sim}/{tag}')

        fd.create_group(f'{sim}/{tag}/Galaxy')
        for k in ['Mbh', 'Mstar', 'Mdm', 'BH_Length']:
            fs.copy(f'{sim}/{tag}/Galaxy/{k}', fd[f'{sim}/{tag}/Galaxy'])

        fd.create_group(f'{sim}/{tag}/Particle')
        for k in ['BH_Age', 'BH_Coordinates', 'BH_ID', 'BH_Index', 'BH_Mass', 'BH_Mdot', 'BH_los', 'BH_sml']:
            fs.copy(f'{sim}/{tag}/Particle/{k}', fd[f'{sim}/{tag}/Particle'])


fd.flush()
