

import h5py
import numpy as np
import pandas as pd

sim = '00'
filename = f'/Users/sw376/Dropbox/Research/data/simulations/flares/blackhole_details_h{sim}.csv'

# details = ascii.read(filename)
details = pd.read_csv(filename, float_precision='round_trip')

PIDs = set(details['PID'])
print(len(PIDs))

output_file = f'/Users/sw376/Dropbox/Research/data/simulations/flares/blackhole_details_h{sim}.h5'

with h5py.File(output_file, 'w') as hf:
    for PID in PIDs:
        print(PID)
        mask = details['PID'] == PID
        bh_history = details.loc[mask].sort_values('Time')
        hf[f'{PID}/a'] = bh_history['Time']
        for col in ['BH_Particle_Mass', 'BH_Subgrid_Mass', 'Mdot']:
            hf[f'{PID}/{col}'] = bh_history[col]

