
import cmasher as cmr

cmap = {}
cmap['log10Mstar'] = cmap['log10Mstar_30'] = cmr.get_sub_cmap('cmr.rainforest', 0.15, 0.85)
cmap['log10FUV'] = cmap['log10LFUV'] = cmr.get_sub_cmap('cmr.sapphire', 0.3, 1.0)

redshift_cmap = 'cmr.infinity'
redshift = cmr.take_cmap_colors(redshift_cmap, 6, cmap_range=(0.1, 0.9))
frontier_redshift = cmr.take_cmap_colors('cmr.gem_r', 6)


density_cmap = cmr.get_sub_cmap('cmr.infinity_s_r', 0.2, 0.8)
