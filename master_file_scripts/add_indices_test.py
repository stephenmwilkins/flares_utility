

import h5py
import flares
import numpy as np

printname = lambda name: print(name)


filedir = '/cosma7/data/dp004/dc-payy1/my_files/flares_pipeline/data/'
filename = filedir + 'flares.hdf5'

#  --- testing on single sim/snapshot
sim = '00'
snap = '008_z007p000'
filedir = ''
filename = filedir + f'flares_{sim}_{snap}.hdf5'


f = h5py.File(filename, 'a')

g = f'{snap}/Galaxy/BPASS_2.2.1/Chabrier300/Indices'
if not g in f:
    f.create_group(g)




f.visit(printname)


# 008_z007p000/Galaxy/BPASS_2.2.1/Chabrier300/SED
# 008_z007p000/Galaxy/BPASS_2.2.1/Chabrier300/SED/DustModelI
# 008_z007p000/Galaxy/BPASS_2.2.1/Chabrier300/SED/Intrinsic
# 008_z007p000/Galaxy/BPASS_2.2.1/Chabrier300/SED/No_ISM
# 008_z007p000/Galaxy/BPASS_2.2.1/Chabrier300/SED/Pure_Stellar
# 008_z007p000/Galaxy/BPASS_2.2.1/Chabrier300/SED/Wavelength


# -----------------------------------------------
# D4000


add_D4000 = False
if add_D4000:

    h = f'{snap}/Galaxy/BPASS_2.2.1/Chabrier300'
    g = f'{h}/Indices/D4000'

    f.create_group(g)

    lam = np.array(f[f'{h}/SED/Wavelength'][:])

    s1 = ((lam>=3750)&(lam<3950)).nonzero()[0]
    s2 = ((lam>=4050)&(lam<4250)).nonzero()[0]

    for spec_type in ['DustModelI','Intrinsic','No_ISM','Pure_Stellar']:
        print(spec_type)
        fnu = np.array(f[f'{h}/SED/{spec_type}'][:])
        f[f'{g}/{spec_type}'] = np.sum(fnu[:, s2], axis=1)/np.sum(fnu[:, s1], axis=1)



# Bolometric Luminosity

add_bolometric = True
if add_bolometric:

    h = f'{snap}/Galaxy/BPASS_2.2.1/Chabrier300'
    g = f'{h}/Indices/Lbol'

    del f[g]

    f.create_group(g)

    lam = np.array(f[f'{h}/SED/Wavelength'][:])

    for spec_type in ['DustModelI','Intrinsic','No_ISM','Pure_Stellar']:
        fnu = np.array(f[f'{h}/SED/{spec_type}'][:])
        flam = fnu*3E8/(lam*1E-10)**2
        f[f'{g}/{spec_type}'] = np.sum(flam, axis=1)*1E-10
        print(spec_type, np.median(np.log10(f[f'{g}/{spec_type}'])))


# UV continuum slope (Beta)

add_UVC = False
if add_UVC:

    h = f'{snap}/Galaxy/BPASS_2.2.1/Chabrier300'
    g = f'{h}/Indices/beta'
    f.create_group(g)

    for spec_type in ['DustModelI','Intrinsic','No_ISM','Pure_Stellar']:
        FUV = f[f'{h}/Luminosity/{spec_type}/FUV'][:]
        NUV = f[f'{h}/Luminosity/{spec_type}/NUV'][:]
        beta = (np.log10(FUV/NUV)/np.log10(1500/2500))-2.0
        print(spec_type, np.median(beta))
        f[f'{g}/{spec_type}'] = beta
