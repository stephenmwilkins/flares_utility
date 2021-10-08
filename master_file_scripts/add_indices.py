

import h5py
import flares
import numpy as np


add_D4000 = True
add_BB = True
add_bolometric = True
add_UVC = True

filedir = '/cosma7/data/dp004/dc-payy1/my_files/flares_pipeline/data/'
filename = filedir + 'flares.hdf5'

fl = flares.flares(filedir + 'flares.hdf5', sim_type='FLARES')
f = h5py.File(filedir + 'flares.hdf5', 'a')


# loop over simulations
for sim in fl.halos:

    # loop over tags
    for tag in fl.tags:

        # create indices groups if it doesn't already exist

        model_tag = f'{sim}/{tag}/Galaxy/BPASS_2.2.1/Chabrier300'
        if not f'{model_tag}/Indices' in f:
            f.create_group(f'{model_tag}/Indices')



        # -----------------------------------------------
        # D4000
        if add_D4000:

            indice_tag = f'{model_tag}/Indices/D4000'

            f.create_group(indice_tag)

            lam = np.array(f[f'{model_tag}/SED/Wavelength'][:])

            s1 = ((lam>=3750)&(lam<3950)).nonzero()[0]
            s2 = ((lam>=4050)&(lam<4250)).nonzero()[0]

            for spec_type in ['DustModelI','Intrinsic','No_ISM','Pure_Stellar']:
                fnu = np.array(f[f'{model_tag}/SED/{spec_type}'][:])
                f[f'{indice_tag}/{spec_type}'] = np.sum(fnu[:, s2], axis=1)/np.sum(fnu[:, s1], axis=1)


        # -----------------------------------------------
        # 3400 - 3600 # 4150 - 4250
        if add_BB_Wilkins:

            indice_tag = f'{model_tag}/Indices/BB_Wilkins'

            f.create_group(indice_tag)

            lam = np.array(f[f'{model_tag}/SED/Wavelength'][:])

            s1 = ((lam>=3400)&(lam<3600)).nonzero()[0]
            s2 = ((lam>=4150)&(lam<4250)).nonzero()[0]

            for spec_type in ['DustModelI','Intrinsic','No_ISM','Pure_Stellar']:
                fnu = np.array(f[f'{model_tag}/SED/{spec_type}'][:])
                f[f'{indice_tag}/{spec_type}'] = (np.sum(fnu[:, s2], axis=1)/np.sum(fnu[:, s1], axis=1))/(len(lam[s2])/len(lam[s1]))


        # -----------------------------------------------
        # Balmer Break using method by Binggeli (https://arxiv.org/pdf/1908.11393.pdf). Just the flux at 4200 / 3500.
        if add_BB_Binggeli:

            indice_tag = f'{model_tag}/Indices/BB_Binggeli'

            f.create_group(indice_tag)

            lam = np.array(f[f'{model_tag}/SED/Wavelength'][:])

            for spec_type in ['DustModelI','Intrinsic','No_ISM','Pure_Stellar']:
                fnu = np.array(f[f'{model_tag}/SED/{spec_type}'][:])
                f[f'{indice_tag}/{spec_type}'] = fnu[:, 4199]/fnu[:, 3499] # I believe this should be correct since d\lam = 1\AA




        # -----------------------------------------------
        # Bolometric Luminosity
        if add_bolometric:

            indice_tag = f'{model_tag}/Indices/Lbol'

            f.create_group(indice_tag)

            lam = np.array(f[f'{model_tag}/SED/Wavelength'][:])

            for spec_type in ['DustModelI','Intrinsic','No_ISM','Pure_Stellar']:
                fnu = np.array(f[f'{model_tag}/SED/{spec_type}'][:])
                flam = fnu*3E8/(lam*1E-10)**2
                f[f'{indice_tag}/{spec_type}'] = np.sum(flam, axis=1)*1E-10


        # -----------------------------------------------
        # UV continuum slope (Beta)

        if add_UVC:

            g = f'{model_tag}/Indices/beta'

            f.create_group(indice_tag)

            for spec_type in ['DustModelI','Intrinsic','No_ISM','Pure_Stellar']:
                FUV = f[f'{model_tag}/Luminosity/{spec_type}/FUV'][:]
                NUV = f[f'{model_tag}/Luminosity/{spec_type}/NUV'][:]
                beta = (np.log10(FUV/NUV)/np.log10(1500/2500))-2.0
                print(spec_type, np.median(beta))
                f[f'{indice_tag}/{spec_type}'] = beta
