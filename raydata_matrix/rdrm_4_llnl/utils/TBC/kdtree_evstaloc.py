#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  kdtree_evstaloc.py
#   Purpose:   find the neighbors of stations and events (for discretization)
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python modules will be imported in this part.
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import os
from pylab import cm as pltcm
from scipy import spatial

# ------------------- INPUT -----------------------------
loc_file = '/home/hosseini/FFINVERSION/INVERSION/Programs/pyfiles/' \
           'pyray_data_matrix/OUTPUTS/P_Pdiff/evsta_locs/evsta_loc.txt'
gc_dist = 1.0
# ------------------- END INPUT -------------------------

loc_arr = np.loadtxt(loc_file, comments='#', skiprows=False)
kd_loc = spatial.KDTree(loc_arr[:, 0:2])
count_neigh = []
counter = 1
loc_size = len(loc_arr)
for loc in loc_arr:
    print '%s/%s ' % (counter, loc_size),
    neighbors = kd_loc.query_ball_point(loc[0:2], gc_dist)
    count_neigh.append([loc[0], loc[1], loc[2], len(neighbors), neighbors])
    counter += 1

plt.ion()
plt.figure()
m = Basemap(projection='robin', lon_0=0.0, lat_0=0.0, resolution='c')
#m.fillcontinents()
m.drawcoastlines(zorder=20)
m.drawparallels(np.arange(-90., 120., 30.))
m.drawmeridians(np.arange(0., 420., 60.))
m.drawmapboundary()

cmap = pltcm.get_cmap('bone_r', 8)
for neigh in count_neigh:
    x, y = m(neigh[1], neigh[0])
    m.scatter(x, y, c=neigh[3], s=100, marker="o", edgecolor="none",
            vmin=0, vmax=8, zorder=10, cmap=cmap,
              alpha=1.0)
cbar = plt.colorbar(orientation='vertical',
                    ticks=range(0, 10, 2), shrink=0.8)
cbar.ax.tick_params(labelsize=24)

plt.figure()
m = Basemap(projection='robin', lon_0=0.0, lat_0=0.0, resolution='c')
#m.fillcontinents()
m.drawcoastlines(zorder=20)
m.drawparallels(np.arange(-90., 120., 30.))
m.drawmeridians(np.arange(0., 420., 60.))
m.drawmapboundary()

high_neigh_fio = open(os.path.join(os.path.dirname(loc_file),
                                   'evsta_high_neighbor.txt'), 'w')
cmap = pltcm.get_cmap('bone_r', 20)
for neigh in count_neigh:
    if neigh[3] < 8:
        continue
    high_neigh_fio.writelines('%s  %s  %s\n' % (neigh[0], neigh[1], neigh[2]))
    x, y = m(neigh[1], neigh[0])
    m.scatter(x, y, c=neigh[3], s=100, marker="o", edgecolor="none",
              vmin=0, vmax=20, zorder=10, cmap=cmap,
              alpha=1.0)
high_neigh_fio.close()
cbar = plt.colorbar(orientation='vertical',
                    ticks=range(0, 22, 2), shrink=0.8)
cbar.ax.tick_params(labelsize=24)
plt.show()

plt.figure()
for i in range(len(count_neigh)):
    plt.plot(i, count_neigh[i][3], 'o', c='k')
plt.xticks(size='large', weight='bold')
plt.yticks(size='large', weight='bold')
plt.xlabel('Station Number', size='x-large', weight='bold')
plt.ylabel('Number of Neighbours', size='x-large', weight='bold')
plt.show()

raw_input('Press enter...')
