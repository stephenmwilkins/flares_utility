
import numpy as np


# --- default plot limits

limits = {}
limits['log10Mstar_30'] = [8.5,11.49]

sfr_limits = [0.,3]
limits['log10SFR_10'] = sfr_limits
limits['log10SFR_inst_30'] = sfr_limits

limits['log10BH_Mass'] = [5.1,9.9]
limits['log10BH_Mdot'] = [-2.9,1.9]

limits['log10Z'] = [-3.19,-1.51]

limits['age'] = [0,400]
limits['log10age'] = [1.4, 2.9]
limits['timescale'] = [0,500]


limits['log10sSFR'] = [-0.9,1.9]
limits['log10sSFR_inst'] = limits['log10sSFR']
limits['log10sSFR_10'] = limits['log10sSFR']

limits['log10HbetaEW'] = [0.01,2.4]
limits['AFUV'] = [0,3.9]
limits['log10FUV'] = [28.1,30.1]

limits['log10M*/LV'] = limits['log10VMTOL'] = [-20.4, -19.1]

limits['BD'] = [2, 5]
limits['log10BD'] = [0.,1.5]

limits['beta'] = [-2.7,-1.1]
limits['log10BB'] = [-0.19,0.49]
limits['log10HbetaEW'] = [1.26,2.74]
limits['log10HI4861_EW'] = limits['log10HbetaEW']
limits['log10OIIIHbEW'] = [1.26, 3.74]


limits['beta_Intrinsic'] = limits['beta_Pure_Stellar'] = [-2.9, -1.7]
limits['log10BB_Intrinsic'] = limits['log10BB_Pure_Stellar'] = limits['log10BB']
limits['log10HbetaEW_Intrinsic'] = [0.01,2.4]


limits['log10SFR50/200'] = [-0.9,0.9]
limits['log10SFR10/50'] = [-0.9,0.9]
limits['log10SFRinst/10'] = [-0.9,0.9]


# limits['log10lambda'] = [0.21, 1.49]
limits['log10lambda'] = np.array([0.21, 1.49])-3 # Myr
limits['nvariance'] = [0., 2.9]
limits['skew'] = [-1.9, 3.9]
limits['nkurtosis'] = [-1.9, 19]
