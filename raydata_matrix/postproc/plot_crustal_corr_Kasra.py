#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  plot_crustal_corr.py
#   Purpose:   simple plot for overview of common corrections for one
#              event, all stations, one band
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import sys
import os

plt.figure(figsize=(15, 12))
ell_ccor_fio = sys.argv[1]
index_2_plot = int(sys.argv[2])

ccorr_arr = np.loadtxt(ell_ccor_fio, delimiter=',', dtype=object)
print 'Passed loadtxt'

mymap = Basemap(projection='robin', lon_0=0, lat_0=0)

x, y = mymap(ccorr_arr[:, 1], ccorr_arr[:, 0])

mymap.scatter(x, y, c=ccorr_arr[:, index_2_plot], cmap=plt.cm.get_cmap('seismic', 256), s=30, zorder=10, edgecolors='none',
              vmin=-2, vmax=2)

mymap.drawcoastlines(color='black', linewidth=5)

cbar = plt.colorbar(orientation='horizontal', shrink=0.5)
cbar.ax.tick_params(labelsize=10)
plt.ion()
plt.show()
