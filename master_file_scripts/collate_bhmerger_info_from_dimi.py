

import numpy as np
import h5py



data_path = '/cosma7/data/dp004/dc-irod1/FLARES/LISA_data'
output_path = '/cosma/home/dp004/dc-wilk2/data/flares/flares-1'


redshifts = [4.77, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00]
tags = ['011_z004p770', '010_z005p000', '009_z006p000', '008_z007p000', '007_z008p000', '006_z009p000', '005_z010p000',
        '004_z011p000', '003_z012p000', '002_z013p000', '001_z014p000', '000_z015p000']
regions = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17',
           '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35',
           '36', '37', '38', '39']



# open hdf5 output file
with h5py.File(f'{output_path}/bhmergers.h5', 'w') as hf:

    for region in regions:

        # Load the remaining data saved in the read-in block of code. #
        a = np.load(data_path + '/times_' + str(region) + '.npy')
        ids_primary = np.load(data_path + '/ids_primary_' + str(region) + '.npy')
        ids_secondary = np.load(data_path + '/ids_secondary_' + str(region) + '.npy')
        masses_primary = np.load(data_path + '/masses_primary_' + str(region) + '.npy')
        masses_secondary = np.load(data_path + '/masses_secondary_' + str(region) + '.npy')

        # Convert the data into arrays of floats. #
        a = np.array(a, dtype=float)
        ids_primary = np.array(ids_primary, dtype=int)
        ids_secondary = np.array(ids_secondary, dtype=int)
        masses_primary = np.array(masses_primary, dtype=float)
        masses_secondary = np.array(masses_secondary, dtype=float)

        hf[f'{region}/a'] = a
        hf[f'{region}/ids_primary'] = ids_primary
        hf[f'{region}/ids_secondary'] = ids_secondary
        hf[f'{region}/masses_primary'] = masses_primary
        hf[f'{region}/masses_secondary'] = masses_secondary


