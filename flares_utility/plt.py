

from . import colors as colors_
from . import labels as labels_
from . import limits as limits_
from . import stats
import flare.plt as fplt
from flare.photom import log10lum_to_M
import cmasher as cmr
from scipy.stats import pearsonr
import numpy as np
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.patheffects as pe
mpl.use('Agg')


limits = limits_.limits
labels = labels_.labels


def fancy(x): return r'$\rm '+x.replace(' ', '\ ')+'$'
def ml(x): return r'$\rm '+x+'$'


default_quantiles = [0.022, 0.158, 0.842, 0.978]

nbins = {}


def get_bins_N(X, limits=limits, bins=20, cuts_=[10, 1]):

    # --- if integer number of bins provided decide the bins based on the limits, otherwise assume bin_edges
    if type(bins) == int:
        bins = np.linspace(*limits, bins)

    bincen = (bins[:-1]+bins[1:])/2.

    N, _ = np.histogram(X, bins=bins)

    i = np.array(range(len(N)))

    cuts = [np.ones(len(N))]*len(cuts_)

    for i, cut in enumerate(cuts_):
        cuts[i] = N > cut

    # ss = i[N<1]
    # if len(ss)>0:
    #     Ncut1 = i[i<ss[0]]
    # else:
    #     Ncut1 = i
    #
    # ss = i[N<10]
    # if len(ss)>0:
    #     Ncut2 = i[i<ss[0]]
    # else:
    #     Ncut2 = i

    return N, bins, bincen, cuts


def measure_weighted_median(X, Y, w, limits=limits, bins=20, weighted=True, cuts_=[10, 1]):

    N, bins, bincen, cuts = get_bins_N(X, limits=limits, bins=bins, cuts_=cuts_)

    if not weighted:
        w = np.ones(len(X))

    med = stats.binned_weighted_quantile(X, Y, w, bins, [0.5])

    return N, bins, bincen, cuts, med


def weighted_median(ax, X, Y, w, limits=limits, bins=20, weighted=True, c='k', label=None, ls=None, lw=1, cuts_=[10, 1], outline=True):

    N, bins, bincen, cuts, med = measure_weighted_median(
        X, Y, w, limits=limits, bins=bins, weighted=weighted, cuts_=cuts_)

    if not ls:
        if outline:
            ax.plot(bincen[cuts[1]], med[cuts[1]], c=c, ls='--', lw=lw,
                    path_effects=[pe.withStroke(linewidth=3, foreground='white')], zorder=3)  # >1, <10
            ax.plot(bincen[cuts[0]], med[cuts[0]], c=c, ls='-', lw=lw, label=label,
                    path_effects=[pe.withStroke(linewidth=3, foreground='white')], zorder=3)  # >= 10
        else:
            ax.plot(bincen[cuts[1]], med[cuts[1]], c=c, ls='--', lw=lw, zorder=3)  # >1, <10
            ax.plot(bincen[cuts[0]], med[cuts[0]], c=c, ls='-',
                    lw=lw, label=label, zorder=3)  # >= 10
    else:
        ax.plot(bincen, med, c=c, ls=ls, lw=lw, label=label, zorder=3)  # >1, <10

    return bincen[cuts[1]], med[cuts[1]], bincen[cuts[0]], med[cuts[0]]


def weighted_range(ax, X, Y, w, limits=limits, bins=20, weighted=True, quantiles=default_quantiles, c='k'):

    N, bins, bincen, cuts = get_bins_N(X, limits=limits, bins=bins)

    if not weighted:
        w = np.ones(len(X))

    P = stats.binned_weighted_quantile(X, Y, w, bins, quantiles)

    if len(quantiles) == 2:
        ax.fill_between(bincen[cuts[0]], P[:, 0][cuts[0]], P[:, -1][cuts[0]], color=c, alpha=0.1)
    if len(quantiles) == 4:
        ax.fill_between(bincen[cuts[0]], P[:, 1][cuts[0]], P[:, -2]
                        [cuts[0]], color=c, alpha=0.1, lw=0)
        ax.fill_between(bincen[cuts[0]], P[:, 0][cuts[0]], P[:, -1]
                        [cuts[0]], color=c, alpha=0.05, lw=0)


def simple_wcbar(D, x, y, z, s=None, labels=labels, limits=limits,  cmap=cm.viridis, add_weighted_median=True):

    #
    # # --- if no selection provided select all galaxies
    # if not s:
    #     s = D[x] == D[x]

    # --- if no limits provided base limits on selected data ranges
    for v in [x, y, z]:
        if v not in limits.keys():
            limits[v] = [np.min(D[v][s]), np.max(D[v][s])]

    # --- if no labels provided just use the name
    for v in [x, y, z]:
        if v not in labels.keys():
            labels[v] = v

    # --- get template figure from flare.plt
    fig, ax, cax = fplt.simple_wcbar()

    # --- define colour scale
    norm = mpl.colors.Normalize(vmin=limits[z][0], vmax=limits[z][1])

    # --- plot
    ax.scatter(D[x][s], D[y][s], s=1, alpha=0.5, c=cmap(norm(D[z][s])))

    # --- weighted median Lines

    if add_weighted_median:
        bins = np.linspace(*limits[x], 20)
        bincen = (bins[:-1]+bins[1:])/2.
        out = stats.binned_weighted_quantile(
            D[x][s], D[y][s], D['weight'][s], bins, [0.84, 0.50, 0.16])

        ax.plot(bincen, out[:, 1], c='k', ls='-')
        # ax.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color='k', alpha = 0.2)

    ax.set_xlim(limits[x])
    ax.set_ylim(limits[y])

    ax.set_ylabel(rf'$\rm {labels[y]}$', fontsize=9)
    ax.set_xlabel(rf'$\rm {labels[x]}$', fontsize=9)

    # --- add colourbar

    cmapper = cm.ScalarMappable(norm=norm, cmap=cmap)
    cmapper.set_array([])
    cbar = fig.colorbar(cmapper, cax=cax, orientation='vertical')
    cbar.set_label(rf'$\rm {labels[z]} $')

    return fig, ax, cax


def simple_wcbar_connected(D, x, y1, y2, z, s=None, labels=labels, limits=limits,  cmap=cm.viridis):
    """ Same as above but you provide two quantities (e.g. intrinsic and dust attenuated beta) to compare the effect """

    # --- if no limits provided base limits on selected data ranges
    for v in [x, y1, y2, z]:
        if v not in limits.keys():
            limits[v] = [np.min(D[v][s]), np.max(D[v][s])]

    # --- if no labels provided just use the name
    for v in [x, y1, y2, z]:
        if v not in labels.keys():
            labels[v] = v

    # --- get template figure from flare.plt
    fig, ax, cax = fplt.simple_wcbar()

    # --- define colour scale
    norm = mpl.colors.Normalize(vmin=limits[z][0], vmax=limits[z][1])

    # --- plot

    for x_, y1_, y2_, z_ in zip(D[x][s], D[y1][s], D[y2][s], D[z][s]):
        ax.plot([x_]*2, [y1_, y2_], c=cmap(norm(z_)), lw=1)

    # print(limits[x])
    #
    # ax.set_xlim(limits[x])
    # ax.set_ylim(limits[y1])

    ax.set_ylabel(rf'$\rm {labels[y1]}$', fontsize=9)
    ax.set_xlabel(rf'$\rm {labels[x]}$', fontsize=9)

    # --- add colourbar

    cmapper = cm.ScalarMappable(norm=norm, cmap=cmap)
    cmapper.set_array([])
    cbar = fig.colorbar(cmapper, cax=cax, orientation='vertical')
    cbar.set_label(rf'$\rm {labels[z]} $')

    return fig, ax, cax


def simple_whist(D, x, y, s=None, labels=labels, limits=limits, base_size=3.5, hist_bins=20, bins=20, full_width=True, add_weighted_median=True, weighted_median_outline=False, add_correlation_coefficient=False, add_weighted_range=True, weighted=True, quantiles=default_quantiles):

    left = 0.15
    height = 0.80
    bottom = 0.15
    width = 0.6
    hwidth = 0.15

    fig = plt.figure(figsize=(base_size, base_size*width/height))

    ax = fig.add_axes((left, bottom, width, height))
    hax = fig.add_axes([left+width, bottom, hwidth, height])

    hax.axis('off')

    # --- if no limits provided base limits on selected data ranges
    for v in [x, y]:
        if v not in limits.keys():
            limits[v] = [np.min(D[v][s]), np.max(D[v][s])]

    # --- if no labels provided just use the name
    for v in [x, y]:
        if v not in labels.keys():
            labels[v] = v

    X = D[x][s]
    Y = D[y][s]
    w = D['weight'][s]

    hax.hist(Y, bins=hist_bins, orientation='horizontal',
             color='k', histtype=u'step', fill=False, density=True)

    hax.set_xticks([])
    hax.set_yticks([])

    # --- weighted median Lines

    # --- weighted median Lines
    if add_weighted_range:
        weighted_range(ax, X, Y, w, limits=limits[x],
                       bins=bins, weighted=weighted, quantiles=quantiles)

    if add_weighted_median:
        if weighted_median_outline:
            b1, m1, b2, m2 = weighted_median(
                ax, X, Y, w, limits=limits[x], bins=bins, weighted=weighted, lw=2, c='w')
        b1, m1, b2, m2 = weighted_median(
            ax, X, Y, w, limits=limits[x], bins=bins, weighted=weighted)

    ax.set_xlim(limits[x])
    ax.set_ylim(limits[y])

    print(labels[y])

    ax.set_ylabel(rf'$\rm {labels[y]}$', fontsize=9)
    ax.set_xlabel(rf'$\rm {labels[x]}$', fontsize=9)

    return fig, ax, hax


def simple_wcbar_whist(D, x, y, z, s=None, labels=labels, limits=limits,  cmap=cm.viridis, base_size=3.5, hist_bins=20, alpha=0.5, bins=20, full_width=True, add_weighted_median=True, weighted_median_outline=False, add_correlation_coefficient=False, weighted=True):

    left = 0.15
    height = 0.70
    bottom = 0.15
    width = 0.6
    hwidth = 0.15

    fig = plt.figure(figsize=(base_size, base_size*width/height))

    ax = fig.add_axes((left, bottom, width, height))
    hax = fig.add_axes([left+width, bottom, hwidth, height])
    cax = fig.add_axes([left, bottom+height, width, 0.03])

    # --- if no limits provided base limits on selected data ranges
    for v in [x, y, z]:
        if v not in limits.keys():
            limits[v] = [np.min(D[v][s]), np.max(D[v][s])]

    # --- if no labels provided just use the name
    for v in [x, y, z]:
        if v not in labels.keys():
            labels[v] = v

    # --- define colour scale
    norm = mpl.colors.Normalize(vmin=limits[z][0], vmax=limits[z][1])

    X = D[x][s]
    Y = D[y][s]
    Z = D[z][s]
    w = D['weight'][s]

    # --- plot
    ax.scatter(X, Y, s=1, alpha=alpha, c=cmap(norm(Z)))

    hax.hist(Y, bins=hist_bins, orientation='horizontal',
             color='k', histtype=u'step', fill=False, density=True)

    hax.set_xticks([])
    hax.set_yticks([])

    # --- weighted median Lines

    if add_weighted_median:
        if weighted_median_outline:
            b1, m1, b2, m2 = weighted_median(
                ax, X, Y, w, limits=limits[x], bins=bins, weighted=weighted, lw=2, c='w')
        b1, m1, b2, m2 = weighted_median(
            ax, X, Y, w, limits=limits[x], bins=bins, weighted=weighted)

    ax.set_xlim(limits[x])
    ax.set_ylim(limits[y])

    ax.set_ylabel(rf'$\rm {labels[y]}$', fontsize=9)
    ax.set_xlabel(rf'$\rm {labels[x]}$', fontsize=9)

    # --- add colourbar

    cmapper = cm.ScalarMappable(norm=norm, cmap=cmap)
    cmapper.set_array([])
    cbar = fig.colorbar(cmapper, cax=cax, orientation='horizontal')
    cbar.set_label(rf'$\rm {labels[z]} $', fontsize=8)
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')

    return fig, ax, cax, hax


def simple(D, x, y, s=None, labels=labels, limits=limits, nbins=nbins, add_weighted_median=True):

    # --- if no limits provided base limits on selected data ranges
    for v in [x, y]:
        if v not in limits.keys():
            limits[v] = [np.min(D[v][s]), np.max(D[v][s])]

    # --- if no labels provided just use the name
    for v in [x, y]:
        if v not in labels.keys():
            labels[v] = v

    # --- if no bins provided just use the name
    for v in [x, y]:
        if v not in nbins.keys():
            nbins[v] = 25

    # --- get template figure from flare.plt
    fig, ax = fplt.simple()

    # --- plot
    ax.scatter(D[x][s], D[y][s], s=1, alpha=0.5, c='k')

    # --- weighted median Lines

    if add_weighted_median:
        bins = np.linspace(*limits[x], nbins[x])
        bincen = (bins[:-1]+bins[1:])/2.
        out = stats.binned_weighted_quantile(
            D[x][s], D[y][s], D['weight'][s], bins, [0.84, 0.50, 0.16])

        ax.plot(bincen, out[:, 1], c='k', ls='-')
        # ax.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color='k', alpha = 0.2)

    ax.set_xlim(limits[x])
    ax.set_ylim(limits[y])

    ax.set_ylabel(rf'$\rm {labels[y]}$', fontsize=9)
    ax.set_xlabel(rf'$\rm {labels[x]}$', fontsize=9)

    return fig


def linear(D, properties, s, labels=labels, limits=limits, scatter_colour_quantity=False, scatter_cmap=None, bins=50, full_width=True):

    # --- if no limits provided base limits on selected data ranges
    for v in properties:
        if v not in limits.keys():
            limits[v] = [np.min(D[v][s]), np.max(D[v][s])]

    # --- if no labels provided just use the name
    for v in properties:
        if v not in labels.keys():
            labels[v] = v

    if scatter_colour_quantity:
        norm = mpl.colors.Normalize(
            vmin=limits[scatter_colour_quantity][0], vmax=limits[scatter_colour_quantity][1])
        cmap = scatter_cmap

    N = len(properties)-1

    if full_width:
        left = 0.1
        top = 0.9
        bottom = 0.2
        right = 0.9
    else:
        left = 0.15
        top = 0.9
        bottom = 0.2
        right = 0.85

    panel_width = (right-left)/N
    panel_height = top-bottom

    if full_width:
        fig, axes = plt.subplots(1, N, figsize=(7, 7/(panel_height/panel_width)), sharey=True)
    else:
        fig, axes = plt.subplots(1, N, figsize=(3.5, 3.5/(panel_height/panel_width)), sharey=True)

    plt.subplots_adjust(left=left, top=top, bottom=bottom, right=right, wspace=0.0, hspace=0.0)

    y = properties[0]

    for ax, x in zip(axes, properties[1:]):

        # --- scatter plot here

        if scatter_colour_quantity:
            ax.scatter(D[x][s], D[y][s], s=1, alpha=0.5,
                       c=cmap(norm(D[scatter_colour_quantity][s])))
        else:
            ax.scatter(D[x][s], D[y][s], s=1, alpha=0.5, c='k')

        # --- weighted median Lines

        bins = np.linspace(*limits[x], 20)
        bincen = (bins[:-1]+bins[1:])/2.
        out = stats.binned_weighted_quantile(
            D[x][s], D[y][s], D['weight'][s], bins, [0.84, 0.50, 0.16])

        ax.plot(bincen, out[:, 1], c='k', ls='-')
        # ax.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color='k', alpha = 0.2)

        ax.set_xlim(limits[x])
        ax.set_ylim(limits[y])

        ax.set_xlabel(rf'$\rm {labels[x]}$', fontsize=8)

    axes[0].set_ylabel(rf'$\rm {labels[y]}$', fontsize=8)

    # --- add colourbar

    if scatter_colour_quantity:

        cmapper = cm.ScalarMappable(norm=norm, cmap=cmap)
        cmapper.set_array([])

        cax = fig.add_axes([right, bottom, 0.015, top-bottom])
        fig.colorbar(cmapper, cax=cax, orientation='vertical')
        cax.set_ylabel(rf'$\rm {labels[scatter_colour_quantity]} $', fontsize=7)
        cax.tick_params(axis='y', labelsize=6)

    return fig


def linear2(D, y, properties, s, labels=labels, limits=limits, scatter_colour_quantity=False, scatter_cmap=None, bins=20, full_width=True, add_weighted_median=True, weighted_median_outline=False, add_correlation_coefficient=False, weighted=True):
    """ an improved version of the above, using elements from linear_mcol """

    z = scatter_colour_quantity

    if not scatter_cmap:
        if scatter_colour_quantity in colors_.cmap.keys():
            scatter_cmap = colors_.cmap[z]
        else:
            scatter_cmap = 'plasma'

    # --- if no limits provided base limits on selected data ranges
    for v in properties:
        if v not in limits.keys():
            limits[v] = [np.min(D[v][s]), np.max(D[v][s])]

    # --- if no labels provided just use the name
    for v in properties:
        if v not in labels.keys():
            labels[v] = v

    if scatter_colour_quantity:
        norm = mpl.colors.Normalize(
            vmin=limits[scatter_colour_quantity][0], vmax=limits[scatter_colour_quantity][1])
        cmap = scatter_cmap

    Nx = len(properties)

    left = 0.1
    top = 0.9
    bottom = 0.2
    right = 0.9

    panel_width = (right-left)/Nx
    panel_height = (top-bottom)

    fig, axes = plt.subplots(1, Nx, figsize=(7, 2), sharey='row')

    plt.subplots_adjust(left=left, top=top, bottom=bottom, right=right, wspace=0.0, hspace=0.0)

    for i, (x, ax) in enumerate(zip(properties, axes)):

        # --- scatter plot here

        X = D[x][s]
        Y = D[y][s]
        w = D['weight'][s]

        if scatter_colour_quantity:
            ax.scatter(X, Y, s=1, alpha=0.5, c=cmap(norm(D[scatter_colour_quantity][s])))

        # --- weighted median Lines

        if add_weighted_median:
            if weighted_median_outline:
                b1, m1, b2, m2 = weighted_median(
                    ax, X, Y, w, limits=limits[x], bins=bins, weighted=weighted, lw=2, c='w')
            b1, m1, b2, m2 = weighted_median(
                ax, X, Y, w, limits=limits[x], bins=bins, weighted=weighted)

        if add_correlation_coefficient:

            # s2 = (~np.isnan(Y))&(~np.isnan(X))&(~np.isnan(Y))&(~np.isnan(X))
            s2 = (np.isfinite(X)) & (np.isfinite(Y))
            r, p = pearsonr(X[s2], Y[s2])

            ax.text(0.1, 0.9, rf'$\rm r={r:.2f}$', horizontalalignment='left', verticalalignment='center',
                    transform=ax.transAxes, fontsize=8, path_effects=[pe.withStroke(linewidth=4, foreground='white')])

        ax.set_xlim(limits[x])
        ax.set_ylim(limits[y])

        ax.set_xlabel(rf'$\rm {labels[x]}$', fontsize=8)

        if i == 0:
            ax.set_ylabel(rf'$\rm {labels[y]}$', fontsize=8)

    # --- add colourbar

    if scatter_colour_quantity:

        cmapper = cm.ScalarMappable(norm=norm, cmap=cmap)
        cmapper.set_array([])

        cax = fig.add_axes([right, bottom, 0.015, top-bottom])
        fig.colorbar(cmapper, cax=cax, orientation='vertical')
        cax.set_ylabel(rf'$\rm {labels[scatter_colour_quantity]} $', fontsize=7)
        cax.tick_params(axis='y', labelsize=6)

    return fig


def vlinear(D, x, properties, s, labels=labels, limits=limits, scatter_colour_quantity=False, scatter_cmap=None, bins=20, full_width=True, add_weighted_median=True, weighted_median_outline=False, add_correlation_coefficient=False, weighted=True):
    """ as above but stacked vertically """

    z = scatter_colour_quantity

    if not scatter_cmap:
        if scatter_colour_quantity in colors_.cmap.keys():
            scatter_cmap = colors_.cmap[z]
        else:
            scatter_cmap = 'plasma'

    # --- if no limits provided base limits on selected data ranges
    for v in properties:
        if v not in limits.keys():
            limits[v] = [np.min(D[v][s]), np.max(D[v][s])]

    # --- if no labels provided just use the name
    for v in properties:
        if v not in labels.keys():
            labels[v] = v

    if scatter_colour_quantity:
        norm = mpl.colors.Normalize(
            vmin=limits[scatter_colour_quantity][0], vmax=limits[scatter_colour_quantity][1])
        cmap = scatter_cmap

    Nx = len(properties)

    left = 0.15
    top = 0.90
    bottom = 0.075
    right = 0.95

    panel_width = (right-left)/Nx
    panel_height = (top-bottom)

    fig, axes = plt.subplots(Nx, 1, figsize=(3.5, np.max([Nx*2, 7.])), sharex=True)

    plt.subplots_adjust(left=left, top=top, bottom=bottom, right=right, wspace=0.0, hspace=0.0)

    for i, (y, ax) in enumerate(zip(properties, axes)):

        # --- scatter plot here

        X = D[x][s]
        Y = D[y][s]
        w = D['weight'][s]

        if scatter_colour_quantity:
            ax.scatter(X, Y, s=1, alpha=0.5, c=cmap(norm(D[scatter_colour_quantity][s])))

        # --- weighted median Lines

        if add_weighted_median:
            if weighted_median_outline:
                b1, m1, b2, m2 = weighted_median(
                    ax, X, Y, w, limits=limits[x], bins=bins, weighted=weighted, lw=2, c='w')
            b1, m1, b2, m2 = weighted_median(
                ax, X, Y, w, limits=limits[x], bins=bins, weighted=weighted)

        if add_correlation_coefficient:

            # s2 = (~np.isnan(Y))&(~np.isnan(X))&(~np.isnan(Y))&(~np.isnan(X))
            s2 = (np.isfinite(X)) & (np.isfinite(Y))
            r, p = pearsonr(X[s2], Y[s2])

            ax.text(0.1, 0.9, rf'$\rm r={r:.2f}$', horizontalalignment='left', verticalalignment='center',
                    transform=ax.transAxes, fontsize=8, path_effects=[pe.withStroke(linewidth=4, foreground='white')])

        ax.set_xlim(limits[x])
        ax.set_ylim(limits[y])

        ax.set_ylabel(rf'$\rm {labels[y]}$', fontsize=8)

    axes[-1].set_xlabel(rf'$\rm {labels[x]}$', fontsize=8)

    # --- add colourbar

    if scatter_colour_quantity:

        cmapper = cm.ScalarMappable(norm=norm, cmap=cmap)
        cmapper.set_array([])

        cax = fig.add_axes([left, top, right-left, 0.02])
        fig.colorbar(cmapper, cax=cax, orientation='horizontal')
        cax.set_xlabel(rf'$\rm {labels[scatter_colour_quantity]} $', fontsize=7)
        cax.tick_params(axis='x', labelsize=6)
        cax.xaxis.set_label_position('top')
        cax.xaxis.tick_top()

    return fig, axes


def linear_mcol(D, diagnostics, properties, s, labels=labels, limits=limits, scatter_colour_quantity=False, scatter_cmap=None, bins=20, full_width=True, add_weighted_median=True, weighted_median_outline=False, add_correlation_coefficient=False, weighted=True):

    z = scatter_colour_quantity

    if not scatter_cmap:
        if scatter_colour_quantity in colors_.cmap.keys():
            scatter_cmap = colors_.cmap[z]
        else:
            scatter_cmap = 'plasma'

    # --- if no limits provided base limits on selected data ranges
    for v in properties:
        if v not in limits.keys():
            limits[v] = [np.min(D[v][s]), np.max(D[v][s])]

    # --- if no labels provided just use the name
    for v in properties:
        if v not in labels.keys():
            labels[v] = v

    if scatter_colour_quantity:
        norm = mpl.colors.Normalize(
            vmin=limits[scatter_colour_quantity][0], vmax=limits[scatter_colour_quantity][1])
        cmap = scatter_cmap

    Ny = len(diagnostics)
    Nx = len(properties)

    if full_width:
        left = 0.1
        top = 0.9
        bottom = 0.1
        right = 0.9
    else:
        left = 0.15
        top = 0.95
        bottom = 0.125
        right = 0.85

    panel_width = (right-left)/Nx
    panel_height = (top-bottom)/Ny

    if full_width:
        fig, axes = plt.subplots(Ny, Nx, figsize=(7, 7/(panel_height/panel_width)), sharey='row')
    else:
        fig, axes = plt.subplots(Ny, Nx, figsize=(
            3.5, 3.5/(panel_height/panel_width)), sharey='row')

    plt.subplots_adjust(left=left, top=top, bottom=bottom, right=right, wspace=0.0, hspace=0.0)

    for j, y in enumerate(diagnostics):

        for i, x in enumerate(properties):

            ax = axes[j, i]

            # --- scatter plot here

            X = D[x][s]
            Y = D[y][s]
            w = D['weight'][s]

            if scatter_colour_quantity:
                ax.scatter(X, Y, s=1, alpha=0.5, c=cmap(norm(D[scatter_colour_quantity][s])))

            # --- weighted median Lines

            if add_weighted_median:
                if weighted_median_outline:
                    b1, m1, b2, m2 = weighted_median(
                        ax, X, Y, w, limits=limits[x], bins=bins, weighted=weighted, lw=2, c='w')
                b1, m1, b2, m2 = weighted_median(
                    ax, X, Y, w, limits=limits[x], bins=bins, weighted=weighted)

            if add_correlation_coefficient:

                # s2 = (~np.isnan(Y))&(~np.isnan(X))&(~np.isnan(Y))&(~np.isnan(X))
                s2 = (np.isfinite(X)) & (np.isfinite(Y))
                r, p = pearsonr(X[s2], Y[s2])

                ax.text(0.1, 0.9, rf'$\rm r={r:.2f}$', horizontalalignment='left', verticalalignment='center',
                        transform=ax.transAxes, fontsize=8, path_effects=[pe.withStroke(linewidth=4, foreground='white')])

            ax.set_xlim(limits[x])
            ax.set_ylim(limits[y])

            ax.set_xlabel(rf'$\rm {labels[x]}$', fontsize=8)

            if i == 0:
                ax.set_ylabel(rf'$\rm {labels[y]}$', fontsize=8)

    # axes[0].set_ylabel(rf'$\rm {labels[y]}$', fontsize = 8)

    # --- add colourbar

    if scatter_colour_quantity:

        cmapper = cm.ScalarMappable(norm=norm, cmap=cmap)
        cmapper.set_array([])

        cax = fig.add_axes([right, bottom, 0.015, top-bottom])
        fig.colorbar(cmapper, cax=cax, orientation='vertical')
        cax.set_ylabel(rf'$\rm {labels[scatter_colour_quantity]} $', fontsize=7)
        cax.tick_params(axis='y', labelsize=6)

    return fig


def zevo(D, zeds, x, y, s, labels=labels, limits=limits, bins=20, colors=colors_.redshift, add_weighted_range=False, fig_size=(3.5, 3.5)):

    # this is actually repeated above but with different definitions.

    fig, ax = fplt.simple(fig_size)

    for z, c in zip(zeds, colors):

        X = D[z][x][s[z]]
        Y = D[z][y][s[z]]
        w = D[z]['weight'][s[z]]

        if add_weighted_range:
            weighted_range(ax, X, Y, w, limits=limits[x],
                           bins=bins, weighted=weighted, quantiles=quantiles)

        weighted_median(ax, X, Y, w, limits=limits[x], bins=bins, c=c, label=rf'$\rm z={z:.0f}$')

    ax.set_xlim(limits[x])
    ax.set_ylim(limits[y])

    ax.set_ylabel(rf'$\rm {labels[y]}$', fontsize=9)
    ax.set_xlabel(rf'$\rm {labels[x]}$', fontsize=9)

    ax.legend(fontsize=7)

    return fig, ax


def zevo_whist(D, zeds, x, y, s, labels=labels, limits=limits, nbins=nbins, bins=50, zevo_cmap='cmr.gem_r', hist_bins=30):

    # this is actually repeated above but with different definitions.

    left = 0.15
    width = 0.6
    hwidth = 0.2
    height = 0.75
    bottom = 0.2

    fig = plt.figure(figsize=(3.5, 3.0))

    ax = fig.add_axes((left, bottom, width, height))
    hax = fig.add_axes([left+width, bottom, hwidth, height])

    colors = cmr.take_cmap_colors(zevo_cmap, len(zeds))

    for z, c in zip(zeds, colors):

        # --- weighted median Lines

        if type(bins) is not np.ndarray:
            bins = np.linspace(*limits[x], bins)

        bincen = (bins[:-1]+bins[1:])/2.
        out = stats.binned_weighted_quantile(
            D[z][x][s[z]], D[z][y][s[z]], D[z]['weight'][s[z]], bins, [0.84, 0.50, 0.16])

        N, be = np.histogram(D[z][x][s[z]], bins=bins)

        ax.plot(bincen, out[:, 1], c=c, ls=':', lw=1)
        ax.plot(bincen[N > 10], out[:, 1][N > 10], c=c, ls='-', lw=1)

        hax.hist(D[z][y][s[z]], bins=hist_bins, orientation='horizontal',
                 color=c, histtype=u'step', fill=False, density=True)

    handles = []
    for z, c in zip(zeds, colors):
        handles.append(Line2D([0], [0], label=rf'$\rm z={z:.0f}$', color=c, lw=1))
    hax.legend(handles=handles, fontsize=7, labelspacing=0.1)

    ax.set_xlim(limits[x])
    ax.set_ylim(limits[y])
    hax.set_ylim(limits[y])

    ax.set_ylabel(rf'$\rm {labels[y]}$', fontsize=9)
    ax.set_xlabel(rf'$\rm {labels[x]}$', fontsize=9)

    # hax.legend(fontsize=7, labelspacing=0.1)

    hax.set_yticks([])
    hax.set_xticks([])

    return fig, ax, hax


def linear_redshift(D, zeds, x, y, s, labels=labels, limits=limits, nbins=nbins, scatter=True, scatter_colour_quantity=False, scatter_cmap=None, bins=20, rows=1, add_weighted_median=True, add_weighted_range=False, lowz=False, add_zevo=False, weighted=True, quantiles=default_quantiles, single_column=False, ylabel_fontsize=9):

    z = scatter_colour_quantity
    if not scatter_cmap:
        if scatter_colour_quantity in colors_.cmap.keys():
            scatter_cmap = colors_.cmap[z]
        else:
            scatter_cmap = 'plasma'

    # --- if no limits provided base limits on selected data ranges
    for v in [x, y]:
        if v not in limits.keys():
            limits[v] = [np.min(D[zeds[-1]][v][s[zeds[-1]]]), np.max(D[zeds[-1]][v][s[zeds[-1]]])]

    # --- if no labels provided just use the name
    for v in [x, y]:
        if v not in labels.keys():
            labels[v] = v

    if scatter_colour_quantity:
        norm = mpl.colors.Normalize(
            vmin=limits[scatter_colour_quantity][0], vmax=limits[scatter_colour_quantity][1])
        cmap = scatter_cmap

    Npanels = len(zeds)

    if add_zevo:
        Npanels += 1

    if single_column:
        fig, axes, (left, right, top, bottom) = fplt.multiple_small(
            Npanels, rows=rows, flatten=False)
    else:
        fig, axes, (left, right, top, bottom) = fplt.multiple(
            Npanels, rows=rows, flatten=False)

    if add_zevo:
        ax_zevo = axes[j, -1]

    for ax, z in zip(axes.flatten(), zeds):

        X = D[z][x][s[z]]
        Y = D[z][y][s[z]]
        w = D[z]['weight'][s[z]]

        # --- scatter plot here
        if scatter:
            if scatter_colour_quantity:
                ax.scatter(X, Y, s=1, alpha=0.5, c=cmap(norm(D[z][scatter_colour_quantity][s[z]])))
            else:
                ax.scatter(X, Y, s=1, alpha=0.1, c='k')

        # --- weighted median Lines
        if add_weighted_range:
            weighted_range(ax, X, Y, w, limits=limits[x],
                           bins=bins, weighted=weighted, quantiles=quantiles)

        if add_weighted_median:
            b1, m1, b2, m2 = weighted_median(
                ax, X, Y, w, limits=limits[x], bins=bins, weighted=weighted)

        ax.set_xlim(limits[x])
        ax.set_ylim(limits[y])

        if rows == 1:
            ax.text(0.5, 1.02, rf'$\rm z={z:.0f}$', horizontalalignment='center',
                    verticalalignment='bottom', transform=ax.transAxes, fontsize=7)
        if rows == 2:
            ax.text(0.5, 1.01, rf'$\rm z={z:.0f}$', horizontalalignment='center',
                    verticalalignment='bottom', transform=ax.transAxes, fontsize=8)
        if rows == 3:
            ax.text(0.5, 1.01, rf'$\rm z={z:.0f}$', horizontalalignment='center',
                    verticalalignment='bottom', transform=ax.transAxes, fontsize=8)

    # --- add final redshift trend to all plots
    if lowz and add_weighted_median:
        for ax in axes.flatten()[:-1]:
            ax.plot(b1, m1, c='0.6', ls='--', lw=2, zorder=2)
            ax.plot(b2, m2, c='0.6', ls='-', lw=2, zorder=2)

    if rows > 1:
        for i in range(rows):
            axes[i, 0].set_ylabel(rf'$\rm {labels[y]}$', fontsize=ylabel_fontsize)
    else:
        axes[0].set_ylabel(rf'$\rm {labels[y]}$', fontsize=ylabel_fontsize)

    # --- add colourbar

    if scatter_colour_quantity:

        cmapper = cm.ScalarMappable(norm=norm, cmap=cmap)
        cmapper.set_array([])

        cax = fig.add_axes([right, bottom, 0.015, top-bottom])
        fig.colorbar(cmapper, cax=cax, orientation='vertical')
        cax.set_ylabel(rf'$\rm {labels[scatter_colour_quantity]} $', fontsize=7)
        if rows == 2:
            cax.set_ylabel(rf'$\rm {labels[scatter_colour_quantity]} $', fontsize=8)
        cax.tick_params(axis='y', labelsize=6)

    fig.text(left+(right-left)/2, 0.04, rf'$\rm {labels[x]}$', ha='center', fontsize=9)

    return fig, axes


def linear_redshift_dual(D, zeds, x, y, z, s, labels=labels, limits=limits, scatter_cmap=None, bins=20, weighted=True, quantiles=default_quantiles):

    if not scatter_cmap:
        if z in colors_.cmap.keys():
            scatter_cmap = colors_.cmap[z]
        else:
            scatter_cmap = 'plasma'

    # --- if no limits provided base limits on selected data ranges
    for v in [x, y, z]:
        if v not in limits.keys():
            limits[v] = [np.min(D[zeds[-1]][v][s[zeds[-1]]]), np.max(D[zeds[-1]][v][s[zeds[-1]]])]

    # --- if no labels provided just use the name
    for v in [x, y]:
        if v not in labels.keys():
            labels[v] = v

    norm = mpl.colors.Normalize(vmin=limits[z][0], vmax=limits[z][1])
    cmap = scatter_cmap

    Npanels = len(zeds)

    left = 0.1
    top = 0.90
    bottom = 0.15
    right = 0.9
    hspace = 0.0
    panel_width = (right-left)/int(Npanels)
    panel_height = (top-bottom-hspace*2)/2
    fig, axes = plt.subplots(2, Npanels, figsize=(
        7, 7/(panel_height/panel_width)), sharey=True, sharex=True)
    plt.subplots_adjust(left=left, top=top, bottom=bottom, right=right, wspace=0.0, hspace=hspace)

    print(axes.shape)

    for i, zed in enumerate(zeds):

        X = D[zed][x][s[zed]]
        Y = D[zed][y][s[zed]]
        Z = D[zed][z][s[zed]]
        w = D[zed]['weight'][s[zed]]

        axes[0, i].scatter(X, Y, s=1, alpha=0.5, c=cmap(norm(Z)))
        b1, m1, b2, m2 = weighted_median(
            axes[0, i], X, Y, w, limits=limits[x], bins=bins, weighted=weighted, lw=2, c='w')
        b1, m1, b2, m2 = weighted_median(
            axes[0, i], X, Y, w, limits=limits[x], bins=bins, weighted=weighted)
        axes[0, i].set_xlim(limits[x])
        axes[0, i].set_ylim(limits[y])

        weighted_range(axes[1, i], X, Y, w, limits=limits[x], bins=bins,
                       weighted=weighted, quantiles=quantiles)
        b1, m1, b2, m2 = weighted_median(
            axes[1, i], X, Y, w, limits=limits[x], bins=bins, weighted=weighted)
        axes[1, i].set_xlim(limits[x])
        axes[1, i].set_ylim(limits[y])

        axes[0, i].text(0.5, 1.01, rf'$\rm z={zed:.0f}$', horizontalalignment='center',
                        verticalalignment='bottom', transform=axes[0, i].transAxes, fontsize=8)

    #
    # # --- add final redshift trend to all plots
    # if lowz and add_weighted_median:
    #     for ax in axes[:-1]:
    #         ax.plot(b1, m1, c='k', ls = '--', lw=2, alpha=0.2)
    #         ax.plot(b2, m2, c='k', ls = '-', lw=2, alpha=0.2)

    # axes[0].set_ylabel(rf'$\rm {labels[y]}$', fontsize = 9)

    # --- add colourbar

    cmapper = cm.ScalarMappable(norm=norm, cmap=cmap)
    cmapper.set_array([])

    cax = fig.add_axes([right, bottom+panel_height+hspace/2, 0.015, panel_height])
    fig.colorbar(cmapper, cax=cax, orientation='vertical')
    cax.set_ylabel(rf'$\rm {labels[z]} $', fontsize=7)
    # if rows==2: cax.set_ylabel(rf'$\rm {labels[z]} $', fontsize = 8)
    cax.tick_params(axis='y', labelsize=6)

    fig.text(left+(right-left)/2, 0.04, rf'$\rm {labels[x]}$', ha='center', fontsize=9)
    fig.text(0.025, bottom+(top-bottom)/2,
             rf'$\rm {labels[y]}$', ha='center', va='center', fontsize=9, rotation=90)

    return fig, axes


def linear_redshift_dual_wmag(D, zeds, x, y, z, s, labels=labels, limits=limits, scatter_cmap=None, bins=20, weighted=True, quantiles=default_quantiles):

    if not scatter_cmap:
        if z in colors_.cmap.keys():
            scatter_cmap = colors_.cmap[z]
        else:
            scatter_cmap = 'plasma'

    # --- if no limits provided base limits on selected data ranges
    for v in [x, y, z]:
        if v not in limits.keys():
            limits[v] = [np.min(D[zeds[-1]][v][s[zeds[-1]]]), np.max(D[zeds[-1]][v][s[zeds[-1]]])]

    # --- if no labels provided just use the name
    for v in [x, y]:
        if v not in labels.keys():
            labels[v] = v

    norm = mpl.colors.Normalize(vmin=limits[z][0], vmax=limits[z][1])
    cmap = scatter_cmap

    Npanels = len(zeds)

    left = 0.1
    top = 0.8
    bottom = 0.15
    right = 0.9
    hspace = 0.0
    panel_width = (right-left)/int(Npanels)
    panel_height = (top-bottom-hspace*2)/2
    fig, axes = plt.subplots(2, Npanels, figsize=(
        7, 7/(panel_height/panel_width)), sharey=True, sharex=True)
    plt.subplots_adjust(left=left, top=top, bottom=bottom, right=right, wspace=0.0, hspace=hspace)

    print(axes.shape)

    for i, zed in enumerate(zeds):

        X = D[zed][x][s[zed]]
        Y = D[zed][y][s[zed]]
        Z = D[zed][z][s[zed]]
        w = D[zed]['weight'][s[zed]]

        axes[0, i].scatter(X, Y, s=1, alpha=0.5, c=cmap(norm(Z)))
        b1, m1, b2, m2 = weighted_median(
            axes[0, i], X, Y, w, limits=limits[x], bins=bins, weighted=weighted, lw=2, c='w')
        b1, m1, b2, m2 = weighted_median(
            axes[0, i], X, Y, w, limits=limits[x], bins=bins, weighted=weighted)
        axes[0, i].set_xlim(limits[x])
        axes[0, i].set_ylim(limits[y])

        weighted_range(axes[1, i], X, Y, w, limits=limits[x], bins=bins,
                       weighted=weighted, quantiles=quantiles)
        b1, m1, b2, m2 = weighted_median(
            axes[1, i], X, Y, w, limits=limits[x], bins=bins, weighted=weighted, outline=False)
        axes[1, i].set_xlim(limits[x])
        axes[1, i].set_ylim(limits[y])

        axes[0, i].text(0.5, 1.41, rf'$\rm z={zed:.0f}$', horizontalalignment='center',
                        verticalalignment='bottom', transform=axes[0, i].transAxes, fontsize=8)

        axm = axes[0, i].twiny()
        axm.set_xlim(log10lum_to_M(np.array(limits[x])))
        print(log10lum_to_M(np.array(limits[x])))
        # axm.set_xticks([-20,-22,-24])

    #
    # # --- add final redshift trend to all plots
    # if lowz and add_weighted_median:
    #     for ax in axes[:-1]:
    #         ax.plot(b1, m1, c='k', ls = '--', lw=2, alpha=0.2)
    #         ax.plot(b2, m2, c='k', ls = '-', lw=2, alpha=0.2)

    # axes[0].set_ylabel(rf'$\rm {labels[y]}$', fontsize = 9)

    # --- add colourbar

    cmapper = cm.ScalarMappable(norm=norm, cmap=cmap)
    cmapper.set_array([])

    cax = fig.add_axes([right, bottom+panel_height+hspace/2, 0.015, panel_height])
    fig.colorbar(cmapper, cax=cax, orientation='vertical')
    cax.set_ylabel(rf'$\rm {labels[z]} $', fontsize=7)
    # if rows==2: cax.set_ylabel(rf'$\rm {labels[z]} $', fontsize = 8)
    cax.tick_params(axis='y', labelsize=6)

    fig.text(left+(right-left)/2, 0.04, rf'$\rm {labels[x]}$', ha='center', fontsize=9)
    fig.text(0.025, bottom+(top-bottom)/2,
             rf'$\rm {labels[y]}$', ha='center', va='center', fontsize=9, rotation=90)

    fig.text(left+(right-left)/2, 0.88, r'$\rm M_{FUV}$', ha='center', fontsize=9)

    return fig, axes


def linear_redshift_density(D, zeds, x, y, s, labels=labels, limits=limits, nbins=20):

    # --- if no limits provided base limits on selected data ranges
    for v in [x, y]:
        if v not in limits.keys():
            limits[v] = [np.min(D[zeds[-1]][v][s[zeds[-1]]]), np.max(D[zeds[-1]][v][s[zeds[-1]]])]

    # --- if no labels provided just use the name
    for v in [x, y]:
        if v not in labels.keys():
            labels[v] = v

    N = len(zeds)

    left = 0.075
    top = 0.9
    bottom = 0.25
    right = 0.85
    panel_width = (right-left)/N
    panel_height = top-bottom
    fig, axes = plt.subplots(1, N, figsize=(7, 7/(panel_height/panel_width)), sharey=True)
    plt.subplots_adjust(left=left, top=top, bottom=bottom, right=right, wspace=0.0, hspace=0.0)

    ldelta_bins = np.array([-0.3, -0.2, -0.1, 0.1, 0.2, 0.3])

    norm = mpl.colors.Normalize(vmin=-0.3, vmax=0.3)

    for ax, z in zip(axes, zeds):

        for ldelta in zip(ldelta_bins[:-1][::-1], ldelta_bins[1:][::-1]):

            c = colors_.density_cmap(norm(np.mean(ldelta)))

            sd = s[z] & (D[z]['ldelta'] > ldelta[0]) & (D[z]['ldelta'] < ldelta[1])

            # --- weighted median Lines

            X = D[z][x][sd]
            Y = D[z][y][sd]
            w = D[z]['weight'][sd]

            _ = weighted_median(
                ax, X, Y, w, limits=limits[x], bins=20, weighted=True, c=c, label=rf'$\rm [{ldelta[0]:.2f},{ldelta[1]:.2f})$')

        ax.set_xlim(limits[x])
        ax.set_ylim(limits[y])

        ax.text(0.5, 1.02, rf'$\rm z={z:.0f}$', horizontalalignment='center',
                verticalalignment='bottom', transform=ax.transAxes, fontsize=7)

    axes[-1].legend(title=r'$\rm \log_{{10}}(1+\delta_{14})$', title_fontsize=8,
                    fontsize=7, labelspacing=0.0, bbox_to_anchor=(1, 0, 2, 1), loc='center left')

    axes[0].set_ylabel(rf'$\rm {labels[y]}$', fontsize=9)

    fig.text(left+(right-left)/2, 0.04, rf'$\rm {labels[x]}$', ha='center', fontsize=9)

    return fig, axes


def linear_redshift_comparison(D, zeds, x, y_, s, labels=labels, limits=limits, bins=20, colors='k', line_styles=['-', '--', '-.', ':'], ylabel=None, ylims=None):

    # --- if no limits provided base limits on selected data ranges
    for v in [x]+y_:
        if v not in limits.keys():
            limits[v] = [np.min(D[zeds[-1]][v][s[zeds[-1]]]), np.max(D[zeds[-1]][v][s[zeds[-1]]])]

    # --- if no labels provided just use the name
    for v in [x]+y_:
        if v not in labels.keys():
            labels[v] = v

    N = len(zeds)

    left = 0.075
    top = 0.9
    bottom = 0.25
    right = 0.95
    panel_width = (right-left)/N
    panel_height = top-bottom
    fig, axes = plt.subplots(1, N, figsize=(7, 7/(panel_height/panel_width)), sharey=True)
    plt.subplots_adjust(left=left, top=top, bottom=bottom, right=right, wspace=0.0, hspace=0.0)

    if type(colors) == str:

        colors = [colors]*len(y_)

    for ax, z in zip(axes, zeds):

        sd = s[z]

        for y, c, ls in zip(y_, colors, line_styles):

            # --- weighted median Lines

            sy = ~np.isnan(D[z][y])

            X = D[z][x][sd & sy]
            Y = D[z][y][sd & sy]
            w = D[z]['weight'][sd & sy]

            _ = weighted_median(
                ax, X, Y, w, limits=limits[x], bins=bins, weighted=True, c=c, lw=1, ls=ls, label=rf'$\rm {labels[y]}$')

        ax.set_xlim(limits[x])

        if ylims:
            ax.set_ylim(ylims)

        ax.text(0.5, 1.02, rf'$\rm z={z:.0f}$', horizontalalignment='center',
                verticalalignment='bottom', transform=ax.transAxes, fontsize=7)

    axes[0].legend(fontsize=6, labelspacing=0.1, loc='upper left')

    axes[0].set_ylabel(rf'$\rm {ylabel}$', fontsize=9)

    fig.text(left+(right-left)/2, 0.04, rf'$\rm {labels[x]}$', ha='center', fontsize=9)

    return fig, axes


def linear_redshift_mcol(D, zeds, x, properties, s, labels=labels, limits=limits, scatter=False, scatter_colour_quantity=False, scatter_cmap=None, add_weighted_median=True, add_weighted_range=False, bins=20, weighted=True, quantiles=default_quantiles, add_linear_fit=False, height=1, add_zevo=False, zevo_cmap=None, fontsize=9, master_label=None):

    # --- if no limits provided base limits on selected data ranges
    for v in [x]+properties:
        if v not in limits.keys():
            limits[v] = [np.min(D[zeds[-1]][v][s[zeds[-1]]]), np.max(D[zeds[-1]][v][s[zeds[-1]]])]

    # --- if no labels provided just use the name
    for v in [x]+properties:
        if v not in labels.keys():
            labels[v] = v

    if scatter_colour_quantity:
        norm = mpl.colors.Normalize(
            vmin=limits[scatter_colour_quantity][0], vmax=limits[scatter_colour_quantity][1])
        cmap = scatter_cmap

    Np = len(properties)
    N = len(zeds)

    if add_zevo:
        N += 1

    fig, axes, (left, right, top, bottom) = fplt.multiplerows(N, Np, height=height)

    if master_label:
        fig.text(0.02, (bottom+top)*0.5, rf'$\rm {master_label}$',
                 fontsize=9, va='center', ha='center', rotation=90)

    if zevo_cmap:
        colors = cmr.take_cmap_colors(zevo_cmap, len(zeds))

    for j, y in enumerate(properties):

        if add_zevo:
            ax_zevo = axes[j, -1]

        for i, z in enumerate(zeds):

            ax = axes[j, i]

            if zevo_cmap:
                c = colors[i]
            else:
                c = 'k'

            # --- weighted median Lines

            X = D[z][x][s[z]]
            Y = D[z][y][s[z]]
            w = D[z]['weight'][s[z]]

            # --- scatter plot here
            if scatter:
                if scatter_colour_quantity:
                    ax.scatter(X, Y, s=1, alpha=0.5, c=cmap(
                        norm(D[z][scatter_colour_quantity][s[z]])))
                else:
                    ax.scatter(X, Y, s=1, alpha=0.1, c='k')

            # --- weighted median Lines

            if add_weighted_range:
                weighted_range(
                    ax, X, Y, w, limits=limits[x], bins=bins, weighted=weighted, quantiles=quantiles, c=c)

            if add_weighted_median:
                b1, m1, b2, m2 = weighted_median(
                    ax, X, Y, w, limits=limits[x], bins=bins, weighted=weighted, c=c)

            if add_zevo:
                ax_zevo.plot(b1, m1, c=c, ls='--', lw=1)
                ax_zevo.plot(b2, m2, c=c, ls='-', lw=1)

            # --- linear first

            if add_linear_fit:

                # fit_p = np.polyfit(D[z][x][s[z]],D[z][y][s[z]], 1, w = D[z]['weight'][s[z]])

                x_ = D[z][x][s[z]]
                y_ = D[z][y][s[z]]
                w_ = D[z]['weight'][s[z]]

                s_ = (~np.isnan(x_)) & (~np.isnan(y_)) & (~np.isinf(y_))  # capture NaNs

                fit_p = np.polyfit(x_[s_], y_[s_], 1, w=w_[s_])
                def fit(x): return fit_p[1] + fit_p[0]*x
                ax.plot(limits[x], fit(np.array(limits[x])), lw=1, c='k', ls='--')

            ax.set_xlim(limits[x])
            ax.set_ylim(limits[y])

        axes[j, 0].set_ylabel(rf'$\rm {labels[y]}$', fontsize=fontsize)

    for i, z in enumerate(zeds):
        axes[0, i].text(0.5, 1.02, rf'$\rm z={z:.0f}$', horizontalalignment='center',
                        verticalalignment='bottom', transform=axes[0, i].transAxes, fontsize=7, color='k')

    # --- add colourbar

    if scatter_colour_quantity:

        cmapper = cm.ScalarMappable(norm=norm, cmap=cmap)
        cmapper.set_array([])

        cax = fig.add_axes([right, bottom, 0.015, top-bottom])
        fig.colorbar(cmapper, cax=cax, orientation='vertical')
        cax.set_ylabel(rf'$\rm {labels[scatter_colour_quantity]} $', fontsize=7)
        cax.tick_params(axis='y', labelsize=6)

    fig.text(left+(right-left)/2, 0.04, rf'$\rm {labels[x]}$', ha='center', fontsize=9)

    return fig, axes


def corner_plot(D, properties, s, labels=labels, limits=limits, scatter_colour_quantity=False, scatter_cmap=None, bins=50):
    return corner_whist(D, properties, s, labels=labels, limits=limits, scatter_colour_quantity=scatter_colour_quantity, scatter_cmap=scatter_cmap, bins=bins)


def corner_whist(D, properties, s, labels=labels, limits=limits, scatter_colour_quantity=False, scatter_cmap=None, bins=50):

    # --- if no limits provided base limits on selected data ranges
    for v in properties:
        if v not in limits.keys():
            limits[v] = [np.min(D[v][s]), np.max(D[v][s])]

    # --- if no labels provided just use the name
    for v in properties:
        if v not in labels.keys():
            labels[v] = v

    if scatter_colour_quantity:
        norm = mpl.colors.Normalize(
            vmin=limits[scatter_colour_quantity][0], vmax=limits[scatter_colour_quantity][1])
        cmap = scatter_cmap

    N = len(properties)

    fig, axes = plt.subplots(N, N, figsize=(7, 7))
    plt.subplots_adjust(left=0.1, top=0.9, bottom=0.1, right=0.9, wspace=0.02, hspace=0.02)

    for i in np.arange(N):
        for j in np.arange(N):
            axes[i, j].set_axis_off()

    for i, x in enumerate(properties):
        for j, y in enumerate(properties[1:][::-1]):

            jj = N-1-j
            ii = i

            ax = axes[jj, ii]

            if j+i < (N-1):
                ax.set_axis_on()

                # --- scatter plot here

                if scatter_colour_quantity:
                    ax.scatter(D[x][s], D[y][s], s=1, alpha=0.5,
                               c=cmap(norm(D[scatter_colour_quantity][s])))
                else:
                    ax.scatter(D[x][s], D[y][s], s=1, alpha=0.5, c='k')

                # --- weighted median Lines

                bins = np.linspace(*limits[x], 20)
                bincen = (bins[:-1]+bins[1:])/2.
                out = stats.binned_weighted_quantile(
                    D[x][s], D[y][s], D['weight'][s], bins, [0.84, 0.50, 0.16])

                ax.plot(bincen, out[:, 1], c='k', ls='-')
                # ax.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color='k', alpha = 0.2)

                ax.set_xlim(limits[x])
                ax.set_ylim(limits[y])

            if i == 0:  # first column
                ax.set_ylabel(rf'$\rm {labels[y]}$', fontsize=7)
            else:
                ax.yaxis.set_ticklabels([])

            if j == 0:  # first row
                ax.set_xlabel(rf'$\rm {labels[x]}$', fontsize=7)
            else:
                ax.xaxis.set_ticklabels([])

            # ax.text(0.5, 0.5, f'x{i}-y{j}', transform = ax.transAxes)

        # --- histograms

        ax = axes[ii, ii]
        ax.set_axis_on()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])

        X = D[x][s]

        H, bin_edges = np.histogram(X, bins=bins, range=limits[x])
        Hw, bin_edges = np.histogram(X, bins=bins, range=limits[x], weights=D['weight'][s])

        Hw *= np.max(H)/np.max(Hw)

        bin_centres = bin_edges[:-1] + (bin_edges[1]-bin_edges[0])*0.5

        ax.fill_between(bin_centres, H*0.0, H, color='0.9')
        ax.plot(bin_centres, Hw, c='0.7', lw=1)

        ax.set_ylim([0.0, np.max(H)*1.2])

    # --- add colourbar

    if scatter_colour_quantity:

        cmapper = cm.ScalarMappable(norm=norm, cmap=cmap)
        cmapper.set_array([])

        cax = fig.add_axes([0.25, 0.87, 0.5, 0.015])
        fig.colorbar(cmapper, cax=cax, orientation='horizontal')
        cax.set_xlabel(rf'$\rm {labels[scatter_colour_quantity]} $')

    return fig


def corner(D, properties, s, labels=labels, limits=limits, scatter_colour_quantity=False, scatter_cmap=None, bins=20, full_width=True, add_weighted_range=False, add_weighted_median=True, weighted=True):

    # --- if no limits provided base limits on selected data ranges
    for v in properties:
        if v not in limits.keys():
            limits[v] = [np.min(D[v][s]), np.max(D[v][s])]

    # --- if no labels provided just use the name
    for v in properties:
        if v not in labels.keys():
            labels[v] = v

    if scatter_colour_quantity:
        norm = mpl.colors.Normalize(
            vmin=limits[scatter_colour_quantity][0], vmax=limits[scatter_colour_quantity][1])
        cmap = scatter_cmap

    N = len(properties)-1

    if full_width:
        fig, axes = plt.subplots(N, N, figsize=(7, 7))
    else:
        fig, axes = plt.subplots(N, N, figsize=(3.5, 3.5))

    if full_width:
        left = 0.1
        right = 0.9
        bottom = 0.1
        top = 0.9
    else:
        left = 0.15
        right = 0.95
        bottom = 0.1
        top = 0.9

    plt.subplots_adjust(left=left, top=top, bottom=bottom, right=right, wspace=0.0, hspace=0.0)

    for i in np.arange(N):
        for j in np.arange(N):
            axes[i, j].set_axis_off()

    for i, x in enumerate(properties[:-1]):
        for j, y in enumerate(properties[1:][::-1]):

            jj = N-1-j
            ii = i

            ax = axes[jj, ii]

            if j+i < N:
                ax.set_axis_on()

                X = D[x][s]
                Y = D[y][s]
                w = D['weight'][s]

                # --- scatter plot here

                if scatter_colour_quantity:
                    ax.scatter(X, Y, s=3, alpha=0.5, c=cmap(
                        norm(D[scatter_colour_quantity][s])), lw=0)

                # --- weighted median Lines
                if add_weighted_range:
                    weighted_range(
                        ax, X, Y, w, limits=limits[x], bins=bins, weighted=weighted, quantiles=quantiles)

                if add_weighted_median:
                    b1, m1, b2, m2 = weighted_median(
                        ax, X, Y, w, limits=limits[x], bins=bins, weighted=weighted)
                # ax.text(0.5, 0.5, f'{ii}_{jj}',transform=ax.transAxes)

                ax.set_xlim(limits[x])
                ax.set_ylim(limits[y])

            if i == 0:  # first column
                ax.set_ylabel(rf'$\rm {labels[y]}$', fontsize=7)
            else:
                ax.yaxis.set_ticklabels([])

            if j == 0:  # first row
                ax.set_xlabel(rf'$\rm {labels[x]}$', fontsize=7)
            else:
                ax.xaxis.set_ticklabels([])

            # ax.text(0.5, 0.5, f'x{i}-y{j}', transform = ax.transAxes)

    # --- add colourbar

    if scatter_colour_quantity:

        cmapper = cm.ScalarMappable(norm=norm, cmap=cmap)
        cmapper.set_array([])

        width = right-left
        height = top - bottom

        cax = fig.add_axes([left+width/2 + 0.025, bottom+height/2 + 0.2, width/2 - 0.025, 0.025])
        fig.colorbar(cmapper, cax=cax, orientation='horizontal')
        cax.set_xlabel(rf'$\rm {labels[scatter_colour_quantity]} $', fontsize=6)
        cax.xaxis.tick_top()
        cax.xaxis.set_label_position('top')
        cax.tick_params(axis='x', labelsize=5)

    return fig


def corner3(D, properties, s, cmaps, labels=labels, limits=limits, bins=50, full_width=True):

    # --- if no limits provided base limits on selected data ranges
    for v in properties:
        if v not in limits.keys():
            limits[v] = [np.min(D[v][s]), np.max(D[v][s])]

    # --- if no labels provided just use the name
    for v in properties:
        if v not in labels.keys():
            labels[v] = v

    N = len(properties)-1

    if full_width:
        left = 0.1
        right = 0.9
        bottom = 0.1
        top = 0.9
        fig, axes = plt.subplots(N, N, figsize=(7, 7))
    else:
        left = 0.15
        right = 0.95
        bottom = 0.1
        top = 0.9
        fig, axes = plt.subplots(N, N, figsize=(3.5, 3.5))

    plt.subplots_adjust(left=left, top=top, bottom=bottom, right=right, wspace=0.0, hspace=0.0)

    for i in np.arange(N):
        for j in np.arange(N):
            axes[i, j].set_axis_off()

    for i, x in enumerate(properties[:-1]):
        for j, y in enumerate(properties[1:][::-1]):

            jj = N-1-j
            ii = i

            ax = axes[jj, ii]

            if j+i < N:
                ax.set_axis_on()

                z = list(set(properties) ^ set([x, y]))[0]

                norm = mpl.colors.Normalize(vmin=limits[z][0], vmax=limits[z][1])
                cmap = cmaps[z]

                # --- scatter plot here

                ax.scatter(D[x][s], D[y][s], s=1, alpha=0.5, c=cmap(norm(D[z][s])))

                # --- weighted median Lines

                bins = np.linspace(*limits[x], 20)
                bincen = (bins[:-1]+bins[1:])/2.
                out = stats.binned_weighted_quantile(
                    D[x][s], D[y][s], D['weight'][s], bins, [0.84, 0.50, 0.16])

                ax.plot(bincen, out[:, 1], c='k', ls='-')
                # ax.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color='k', alpha = 0.2)

                # ax.text(0.5, 0.5, f'{ii}_{jj}',transform=ax.transAxes)

                ax.set_xlim(limits[x])
                ax.set_ylim(limits[y])

            if i == 0:  # first column
                ax.set_ylabel(rf'$\rm {labels[y]}$', fontsize=7)
            else:
                ax.yaxis.set_ticklabels([])

            if j == 0:  # first row
                ax.set_xlabel(rf'$\rm {labels[x]}$', fontsize=7)
            else:
                ax.xaxis.set_ticklabels([])

            # ax.text(0.5, 0.5, f'x{i}-y{j}', transform = ax.transAxes)

    # --- add colourbar

    for i, z in enumerate(properties):

        norm = mpl.colors.Normalize(vmin=limits[z][0], vmax=limits[z][1])
        cmap = cmaps[z]

        cmapper = cm.ScalarMappable(norm=norm, cmap=cmap)
        cmapper.set_array([])

        width = right - left
        height = top - bottom

        cax = fig.add_axes([left+width/2 + 0.025, bottom+height/2 + i*width /
                            (2*len(properties)) + 0.025, width/2 - 0.025, 0.015])
        fig.colorbar(cmapper, cax=cax, orientation='horizontal')
        cax.set_xlabel(rf'$\rm {labels[z]} $', fontsize=6)
        cax.xaxis.tick_top()
        cax.xaxis.set_label_position('top')
        cax.tick_params(axis='x', labelsize=5)

    return fig, axes
