#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  compare_crustal_plot_Kasra.py
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

def plot_map(ccorr_arr_1, ccorr_arr_2=[], arr_indx=4, id=111, plot_bar=False):
    plt.subplot(id)
    mymap = Basemap(projection='robin', lon_0=0, lat_0=0)
    x, y = mymap(ccorr_arr_1[:, 1], ccorr_arr_1[:, 0])

    if len(ccorr_arr_2) == 0:
        mymap.scatter(x, y, c=ccorr_arr_1[:, arr_indx],
                      cmap=plt.cm.get_cmap('seismic', 256), s=5, zorder=10, edgecolors='none',
                      vmin=-2, vmax=2)
    else:
        mymap.scatter(x, y, c=ccorr_arr_1[:, arr_indx] - ccorr_arr_2[:, arr_indx],
                     cmap=plt.cm.get_cmap('seismic', 256), s=5, zorder=10, edgecolors='none',
                      vmin=-2, vmax=2)
    mymap.drawcoastlines(color='black')

    if plot_bar:
        cbar = plt.colorbar(orientation='horizontal', shrink=0.5 )
        cbar.ax.tick_params(labelsize=10)

plt.figure()
ell_ccor_fio_1 = sys.argv[1]
ell_ccor_fio_2 = sys.argv[2]
ind_el_crust = int(sys.argv[3])

ccorr_arr_1 = np.loadtxt(ell_ccor_fio_1, delimiter=',')
ccorr_arr_2 = np.loadtxt(ell_ccor_fio_2, delimiter=',')

plot_map(ccorr_arr_1, [], ind_el_crust, 221, plot_bar=True)
plot_map(ccorr_arr_2, [], ind_el_crust, 222, plot_bar=True)
plot_map(ccorr_arr_1, ccorr_arr_2, ind_el_crust, 223, plot_bar=True)

plt.subplot(224)
diff_arr = ccorr_arr_1[:, ind_el_crust] - ccorr_arr_2[:, ind_el_crust]
plt.scatter(range(len(diff_arr)), abs(diff_arr), s=2)
plt.axhline(np.mean(abs(diff_arr)), ls='--')
print "Mean: %s" % np.mean(abs(diff_arr))
print "STD: %s" % np.std(abs(diff_arr))

plt.ion()
plt.show()
