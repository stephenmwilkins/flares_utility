

import numpy as np

from astropy.io import ascii
from astropy.table import Table, Column

from . import stats
from . import limits as limits_
from . import labels as labels_

limits = limits_.limits
labels = labels_.labels


def simple(D, x, y, labels = labels, limits = limits, bins = 50, formats = {'x': '%.1f', 'y': '%.2f'}, print = True):

    if type(bins) is int:
        bins = np.linspace(*limits[x], bins)

    bincen = (bins[:-1]+bins[1:])/2.

    out = stats.binned_weighted_quantile(D[x],D[y], D['weight'], bins, [0.16,0.50,0.84])

    N, _ = np.histogram(D[x], bins = bins)

    s = N>5

    formats_ = {'x': formats['x'], 'P16': formats['y'], 'P50': formats['y'], 'P84': formats['y']}

    table = {'x': bincen[s], 'P16': out.T[0][s], 'P50': out.T[1][s], 'P84': out.T[2][s]}
    if print:
        ascii.write(table, Writer=ascii.Latex, formats = formats_)
    else:
        return table



def redshift(D_, zeds, x, y, filename = None, latex = True, labels = labels, limits = limits, bins = 50, formats = {'x': '%.1f', 'y': '%.2f'}, percentiles = np.array([2.2, 15.8,50.,84.2,97.8 ])):

    if type(bins) is int:
        bins = np.linspace(*limits[x], bins)

    bincen = (bins[:-1]+bins[1:])/2.

    table = {x: bincen}
    formats_ = {x: formats['x']}

    header1 = f' ${labels_.quantities[x]}$'
    header2 = f' ${labels_.units[x]}$'

    for z in zeds:

        header1 += rf' & \multicolumn{{ {len(percentiles)} }}{{c}}{{ $z={z:.0f}$ }}'
        header2 += rf' & '+' & '.join([rf' P$_{{ {percentile*100:.0f} }} $' for percentile in percentiles])

        D = D_[z]

        P = stats.binned_weighted_quantile(D[x],D[y], D['weight'], bins, percentiles/100.)

        N, _ = np.histogram(D[x], bins = bins)

        s = N>5

        for i, percentile in enumerate(percentiles):

            formats_[f'{y}_z{z}_P{percentile}'] = formats['y']
            table[f'{y}_z{z}_P{percentile}'] = P.T[i]
            table[f'{y}_z{z}_P{percentile}'][~s] = None

    if latex:
        print(header1 + r'\\')
        print(header2 + r'\\')
        ascii.write(table, Writer=ascii.Latex, formats = formats_)

    if filename:
        ascii.write(table, filename, formats = formats_, overwrite = True)



def binned_single_z(D, x, y, limits = limits, bins = 50, percentiles = np.array([2.2, 15.8,50.,84.2,97.8 ])):

    if type(bins) is int:
        bins = np.linspace(*limits[x], bins)

    bincen = (bins[:-1]+bins[1:])/2.
    bincen = np.round(bincen, 2)

    t = Table()




    N, _ = np.histogram(D[x], bins = bins)
    P = stats.binned_weighted_quantile(D[x],D[y], D['weight'], bins, percentiles/100.)

    s = N>1


    if x[:5] == 'log10':
        x_ = x[5:]
        x_log10 = True
    else:
        x_ = x
        x_log10 = False

    if y[:5] == 'log10':
        y_ = y[5:]
        y_log10 = True
    else:
        y_ = y
        y_log10 = False

    if x_log10:
        unit = f'dex({labels_.unit[x_]})'
    else:
        unit = labels_.unit[x_]

    col = Column(data=bincen[s], name = x, unit = unit )
    t.add_column(col)

    t['N'] = N[s]

    for i, percentile in enumerate(percentiles):

        if y_log10:
            unit = f'dex({labels_.unit[y_]})'
        else:
            unit = labels_.unit[y_]

        col = Column(data = np.round(P.T[i][s], 3), name = f'{y}_P{percentile}', unit = unit)
        t.add_column(col)


    return t



# def standard(D, x, y, redshifts, limits = limits, bins = 50, percentiles = np.array([2.2, 15.8,50.,84.2,97.8 ])):
#
#     if type(bins) is int:
#         bins = np.linspace(*limits[x], bins)
#
#     bincen = (bins[:-1]+bins[1:])/2.
#     bincen = np.round(bincen, 2)
#
#     t = {}
#
#     for z in redshifts:
#
#         t[z] = Table()
#
#         P = stats.binned_weighted_quantile(D[z][x],D[z][y], D['weight'], bins, percentiles/100.)
#
#         N, _ = np.histogram(D[z][x], bins = bins)
#
#         s = N>1
#
#         t[z][x] = bincen[s]
#         t[z]['N'] = N[s]
#
#         for i, percentile in enumerate(percentiles):
#             t[z][f'{y}_P{percentile}'] = np.round(P.T[i], 3)
#
#     return t
