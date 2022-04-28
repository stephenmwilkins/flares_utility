


log10x = lambda x: rf'\log_{{10}}({x})'


# quantities

quantities = {}
quantities['Mstar'] = quantities['Mstar_30'] = r'M_{\star}'

for f in ['FUV','NUV','V']:
    quantities[f] = rf'L_{{ {f} }}'
    quantities['L'+f] = rf'L_{{ {f} }}'

quantities['beta'] = r'\beta'
# quantities['CIIIEW'] = r'EW_0([CIII],CIII]\lambda\lambda 1907,1909\AA)'
quantities['CIIIEW'] = r'EW_0([CIII],CIII])'
# units

units = {}
units['Mstar'] = units['Mstar_30'] = r'M_{\odot}'

for f in ['FUV','NUV','V']:
    units[f] = r'erg\ s^{{-1}}\ Hz^{{-1}}'
    units['L'+f] = r'erg\ s^{{-1}}\ Hz^{{-1}}'

units['beta'] = None
units['CIIIEW'] =r'\AA'


# composite labels

labels = {}
for k in list(quantities.keys()):
    if units[k]:
        labels[k] = f'{quantities[k]}/{units[k]}'
    else:
        labels[k] = f'{quantities[k]}'

    quantities['log10'+k] = rf'\log_{{10}}({quantities[k]})'
    labels['log10'+k] = rf'\log_{{10}}({labels[k]})'
    units['log10'+k] = units[k]







# labels['log10Mstar_30'] = rf'\log_{{10}}(M_{{\star}}/ {units['Mstar']} )'
# labels['log10Mstar'] = rf'\log_{{10}}(M_{{\star}}/ {units['Mstar']} )'

# for f in ['FUV','NUV','V']:
#     labels[f'log10{f}'] = rf'\log_{{10}}(L_{{ {f} }}/{units["FUV"]})'
#     labels[f'log10{f}_DustModelI'] = labels['log10FUV']
#     labels[f'log10{f}_Intrinsic'] = rf'\log_{{10}}(L_{{ {f} }}^{{intrinsic}}/erg\ s^{{-1}}\ Hz^{{-1}})'



labels['OIIIHbEW'] = r'EW([OIII+H\beta])/\AA'

labels['log10BB'] = r'\log_{10}(L_{4200}/L_{3500})'
labels['log10BB_Intrinsic'] = rf"{labels['log10BB']}^{{\rm int}}"
labels['log10BB_Pure_Stellar'] = rf"{labels['log10BB']}^{{\star}}"

labels['beta'] = r'\beta'
labels['beta_DustModelI'] = labels['beta']
labels['beta_Intrinsic'] = rf"{labels['beta']}^{{\rm int}}"
labels['beta_Pure_Stellar'] = rf"{labels['beta']}^{{\star}}"

labels['log10HbetaEW'] = r'\log_{10}(H\beta\ EW/\AA)'
labels['log10HbetaEW_Intrinsic'] = rf"{labels['log10HbetaEW']}^{{\rm int}}"

labels['log10SFR_inst_30'] = r'\log_{10}({\rm SFR_{inst}}/{\rm M_{\odot}\ yr^{-1})}'
labels['log10SFR_10'] = r'\log_{10}({\rm SFR_{10}}/{\rm M_{\odot}\ yr^{-1})}'
labels['log10SFR'] = r'\log_{10}({\rm SFR}/{\rm M_{\odot}\ yr^{-1})}'





labels['log10sSFR'] = r'\log_{10}({\rm sSFR}/{\rm Gyr^{-1})}'
labels['log10sSFR_inst'] = r'\log_{10}({\rm sSFR_{inst}}/{\rm Gyr^{-1})}'
labels['log10sSFR_10'] = r'\log_{10}({\rm sSFR_{10}}/{\rm Gyr^{-1})}'

labels['age'] = r'age/{\rm Myr}'
labels['log10age'] = log10x(labels['age'])

labels['MassWeightedStellarAge'] = r'age/{\rm Myr}'
labels['log10MassWeightedStellarAge'] = log10x(labels['MassWeightedStellarAge'])


labels['timescale'] = r'SF\ timescale\ \tau_{SF}/{\rm Myr}'


labels['MassWeightedGasZ'] = r'Z_{g}'
labels['log10MassWeightedGasZ'] = log10x(labels['MassWeightedGasZ'])

labels['Z'] = r'Z'
labels['log10Z'] = log10x(labels['Z'])

labels['log10SFRinst/10'] = r'log_{10}(SFR_{inst}/SFR_{10})'
labels['log10SFR10/50'] = r'log_{10}(SFR_{10}/SFR_{50})'
labels['log10SFR50/200'] = r'log_{10}(SFR_{50}/SFR_{200})'

labels['R'] = labels['log10SFR50/200']
labels['R2'] = labels['log10SFR10/50']

labels['lambda'] = r'\lambda/Myr^{-1}'
labels['log10lambda'] = rf'log_{{10}}{labels["lambda"]}'

labels['skew'] = r'skew'
labels['nvariance'] = r'(\lambda\sigma)^2'
labels['nkurtosis'] = r'excess\ kurtosis'




labels['AFUV'] = r'A_{FUV}'
labels['log10BH_Mass'] = r'\log_{10}({\rm M_{\bullet}}/{\rm M_{\odot})}'
labels['log10BH_Mdot'] = r'\log_{10}({\rm M_{\bullet}}/{\rm M_{\odot}\ yr^{-1})}'
# labels['log10VMTOL'] = labels['log10M*/LV'] = r'\log_{10}[({\rm M_{\star}}/L_{V})/({\rm M_{\odot}/erg\ s^{-1}\ Hz^{-1})]'
labels['log10VMTOL'] = labels['log10M*/LV'] = r'\log_{10}({\rm M_{\star}}/L_{V})'


labels['BD'] = r'L_{H\alpha}/L_{H\beta}'
labels['log10BD'] = log10x(labels['BD'])
