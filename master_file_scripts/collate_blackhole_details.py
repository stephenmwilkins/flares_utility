

import h5py
import numpy as np
import pandas as pd
from astropy.cosmology import Planck13 as cosmo


halos = np.array([
    '00','01','02','03','04',
    '05','06','07','08','09',
    '10','11','12','13','14',
    '15','16','17','18','19',
    '20','21','22','23','24',
    '25','26','27','28','29',
    '30','31','32','33','34',
    '35','36','37','38','39'])


output_file = f'/cosma/home/dp004/dc-wilk2/data/flares/flares-1/blackhole_details.h5'

with h5py.File(output_file, 'w') as hf:

    for sim in halos:

        filename = f'/cosma7/data/dp004/dc-love2/codes/flares_passive/analysis/data/blackhole_details_h{sim}.csv'

        details = pd.read_csv(filename, float_precision='round_trip')

        PIDs = set(details['PID'])
        print(sim, len(PIDs))

        for PID in PIDs:  
            mask = details['PID'] == PID
            bh_history = details.loc[mask].sort_values('Time')
            hf[f'{sim}/{PID}/a'] = bh_history['Time']
            hf[f'{sim}/{PID}/z'] = (1. / bh_history['Time']) - 1
            hf[f'{sim}/{PID}/aou'] = cosmo.age(bh_history['z']).value
            for col in ['BH_Particle_Mass', 'BH_Subgrid_Mass', 'Mdot']:
                hf[f'{PID}/{col}'] = bh_history[col]

