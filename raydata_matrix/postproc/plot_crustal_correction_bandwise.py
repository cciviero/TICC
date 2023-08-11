#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  plot_crustal_correction_bandwise.py
#   Purpose:   simple plot for overview of common corrections for
#              all stations/events, one band
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.basemap import Basemap
import numpy as np
import sys
import os
import glob
import numpy as np

rdrm_path = '/mnt/seismodata/MT/SEA-SEIS_TOMO/rdrm/rdrm_IA_21.10.2021'
bands = ['band01', 'band02', 'band03', 'band04', 'band05', 'band06', 'band07', 'band08']

# cmap = plt.get_cmap('jet_r')
# colors = cmap(np.linspace(0, 1.0, len(bands)))

for band in bands:
    print ('Working on %s' % band)
    counter = 0
    # bandwise will be plotted, therefore for new band we can delete the last contend
    band_array = np.array([])
    longitude = np.array([])
    latitude = np.array([])

    path_2_ell_ccorr = glob.glob(os.path.join(rdrm_path, '*' + band + '*', '*ell_ccor*'))
    for file in path_2_ell_ccorr:
        ccorr_arr = np.loadtxt(file, delimiter=',', dtype=object)
        # 4th column in the file is the tau/common correction column

        longitude = np.append(longitude, ccorr_arr[:, 1].astype(float))
        latitude = np.append(latitude, ccorr_arr[:, 0].astype(float))
        band_array = np.append(band_array, ccorr_arr[:, 4].astype(float))

    plt.figure(figsize=(15, 12))
    plt.ioff()

    mymap = Basemap(projection='robin', lon_0=0, lat_0=0)
    # import ipdb; ipdb.set_trace()

    x, y = mymap(longitude, latitude)
    cmap = plt.cm.get_cmap('jet')
    norm = matplotlib.colors.BoundaryNorm(np.arange(-2, 2, 0.3), cmap.N)
    #mymap.scatter(x, y, c=band_array, cmap=plt.cm.get_cmap('seismic', 256), s=30, zorder=10,
    #               edgecolors='none', vmin=-2, vmax=2, alpha=0.3)
    mymap.scatter(x, y, c=band_array, cmap=cmap, norm=norm, s=30, zorder=10,
                                 edgecolors='none', vmin=-2, vmax=2, alpha=0.5)

    mymap.drawcoastlines(color='grey', linewidth=2)

    cbar = plt.colorbar(orientation='horizontal', shrink=0.5)
    cbar.ax.tick_params(labelsize=10)
    name_part1 = os.path.basename(rdrm_path)
    save_name = os.path.join(rdrm_path, name_part1 + '_' + band + '.png')
    plt.savefig(save_name)
    counter += 1
    plt.clf()
    # import ipdb; ipdb.set_trace()

print('Finished saving... ')