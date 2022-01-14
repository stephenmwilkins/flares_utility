
import h5py

f = '/cosma7/data/dp004/dc-love2/codes/flares/data/flares.hdf5' # current cosma main version
# f = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_highz.hdf5'


x1 = 'Galaxy/Mstar'
x2 = 'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI/FUV'
x2 = 'Galaxy/SFR_aperture/SFR_30/SFR_50_Myr'

with h5py.File(f, mode='r') as hf:

    sims = hf.keys()
    tags = hf[list(sims)[0]].keys()

    for sim in sims:
        for tag in tags:

            try:
                Ngal_1 = len(hf[sim][tag][x1][:])
            except:
                Ngal_1 = 'FAILED'

            try:
                Ngal_2 = len(hf[sim][tag][x2][:])
            except:
                Ngal_2 = 'FAILED'

            if Ngal_2 == 'FAILED' or Ngal_1 == 'FAILED':
                print(sim, tag, Ngal_1, Ngal_2)
            else:
                if Ngal_1-Ngal_2>0:
                    print(sim, tag, Ngal_1, Ngal_2, Ngal_1-Ngal_2, '*'*10)
                else:
                    print(sim, tag, Ngal_1, Ngal_2, Ngal_1-Ngal_2)
