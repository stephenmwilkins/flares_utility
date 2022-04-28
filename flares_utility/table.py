

import numpy as np

from astropy.io import ascii

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
