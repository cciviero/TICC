#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  TomoModelViz.py
#   Purpose:   Visualization tool for tomographic models
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GNU Lesser General Public License, Version 3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
import sys

from utils import read_tomo_file
from utils import plot_plates, plot_hotspots, plot_great_circle, plot_stations
from utils import depth_profile, vertical_profile, plot_cross_section

import warnings
warnings.filterwarnings("ignore")

# ---------------- INPUTS
# '../models/kasra_june_2016/solx.02.dlnVp.iter4'
# '../models/pri05/solx.PRI-P05'
# '../models/pri05/solx.PRI-S05'
# HMSL-P06_percent.nc
# HMSL-S06_percent.nc
# GYPSUM_P_percent.nc
# GYPSUM_S_percent.nc
# spani_P.npy
# spani_S.npy
#tomo_file = '../models/nc/spani_S.npy'
tomo_file = '../models/maria_aug_2016/homo_400_rhum_rum_P_v0103/solx.05.dlnVp.homo_400'
# tomo_file = '../models/maria_aug_2016/homo_400_rhum_rum_P_v0103/columndensity.1.homo_400'

# netcdf format
netcdf_format = False
dvp_dvs = 'dvp'

# min and max for the colorbar
vmin = -0.02
vmax = 0.02

# vmin = 0
# vmax = 700.
# plot_mode:
# default
# poles
plot_mode = 'default'
# min and max radius for depth profiles
min_r = 3482
max_r = 6371
# depth steps in km for depth profiles
step_km = 20
# number of neighbours to be considered for interpolation
num_neigh = 4
# Projection
map_proj = "merc"
lon_0 = 55
llcrnrlat = -40
urcrnrlat = 0
llcrnrlon = 0
urcrnrlon = 90

# max depth to plot plates:
dplate = 800
num_lat_lon_cuts = 720
# remove mean (depth slices)
remove_mean = False

# cmap:
# 1. choose from the available colormaps
# 2. 'custom': change the function manually
cmap=plt.cm.hot_r
cmap = 'seismic_r'
eradius = 6371.
# sigma for Gauss filter applied to interpolated values (False: if no filter)
sigma_gauss_filter = 1.0
# ----------------

if 'custom' in cmap:
    from utils import custom_div_cmap
    # cmap = custom_div_cmap(numcolors=11, mincol='white', midcol='yellow', mid2col='red', maxcol='black')
                           # mincol='red', midcol='white', maxcol='blue'
    cmap = custom_div_cmap(numcolors=11, mincol='red', midcol='white', maxcol='blue')

    # brown-white-teal
    # cmap = custom_div_cmap(numcolors=11,
    #                        mincol='#8c510a', maxcol='#01665e')

print("===> reading %s" % tomo_file)
tomo_rd = read_tomo_file(tomo_file,
                         dvp_dvs=dvp_dvs,
                         netcdf_format=netcdf_format,
                         eradius=eradius)
print("===> creating kd-tree for the tomo model")
kd_tree = cKDTree(data=tomo_rd[:, 0:3], leafsize=10)
print 'Done!'
# XXX
# cut_type = raw_input("\nSelect one of the following options:\n"
#                      "1. vertical cross-sections\n")
cut_type = "1"

# Vertical cross-sections
if "1" in cut_type:
    fig = plt.figure(figsize=(15, 10))
    plt.ion()
    plt.subplot(2, 2, 1)

    m = Basemap(projection=map_proj, lon_0=lon_0, lat_0=0, llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon)
    m.drawcoastlines(zorder=100)

    hdepth = float(raw_input("Depth for the depth-profile:  "))
    hlonlats_arr, hxyz = depth_profile(hdepth, eradius=eradius,
                                       cuts=num_lat_lon_cuts)

    print("===> querying kd-tree for depth profile")
    nps = kd_tree.query(hxyz, k=num_neigh)
    print 'again done...'

    # interpolation
    inter_values = np.zeros([np.shape(hxyz)[0], 4])
    inter_values[:, 0:3] = hxyz[:, 0:3]
    for in1 in range(np.shape(nps)[1]):
        inter_values[in1][3] = \
            np.sum(1.0 / nps[0][in1] * tomo_rd[nps[1][in1]][:, 3]) / \
            np.sum(1.0 / nps[0][in1])
    reshape_interpolate = np.reshape(inter_values[:, 3],
                                     [num_lat_lon_cuts, -1])
    if sigma_gauss_filter:
        reshape_interpolate = gaussian_filter(reshape_interpolate,
                                              sigma=sigma_gauss_filter,
                                              order=0)

    hlons_grd, hlats_grd = np.meshgrid(hlonlats_arr[:, 0], hlonlats_arr[:, 1])
    if remove_mean:
        reshape_mean = np.mean(reshape_interpolate.T)
    else:
        reshape_mean = 0
    m.pcolormesh(hlons_grd, hlats_grd, reshape_interpolate.T - reshape_mean,
                 latlon=True, zorder=10, cmap=cmap,
                 vmin=vmin, vmax=vmax)
    plt.title("Radius: %s -- Depth: %s -- Removed Mean: %.4f%%" %
              (eradius - hdepth, hdepth, reshape_mean*100.))

    plot_hotspots(m)
    plot_stations(m)
    if hdepth <= dplate:
        plot_plates(m)
    plt.show()

    plt.subplot(2, 2, 2)
    if 'poles' in plot_mode:
        m2 = Basemap(projection='nplaea', boundinglat=50, lon_0=0)
        m2.drawcoastlines(zorder=100)
        m2.pcolormesh(hlons_grd, hlats_grd,
                      reshape_interpolate.T - reshape_mean,
                      latlon=True, zorder=10, cmap=cmap,
                      vmin=vmin, vmax=vmax)
        plot_hotspots(m2)
        if hdepth <= dplate:
            plot_plates(m2)

        plt.subplot(224)
        m3 = Basemap(projection='splaea', boundinglat=-50, lon_0=180)
        m3.drawcoastlines(zorder=100)
        m3.pcolormesh(hlons_grd, hlats_grd,
                      reshape_interpolate.T - reshape_mean,
                      latlon=True, zorder=10, cmap=cmap,
                      vmin=vmin, vmax=vmax)
        plot_hotspots(m3)
        if hdepth <= dplate:
            plot_plates(m3)

    else:
        m2 = Basemap(projection=map_proj, lon_0=lon_0, lat_0=0, llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat,
                    llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon)
        m2.drawcoastlines(zorder=100)
        m2.pcolormesh(hlons_grd, hlats_grd,
                      reshape_interpolate.T - np.mean(reshape_interpolate.T),
                      latlon=True, zorder=10, cmap=cmap,
                      vmin=vmin, vmax=vmax)
        plot_hotspots(m2)
        plot_stations(m2)
        if hdepth <= dplate:
            plot_plates(m2)
        plt.title("Mean: %.4f %%" % (np.mean(reshape_interpolate.T)*100))

    path_num = 1
    while True:
        print(20*"-")
        onscreen = raw_input("select on the screen? (y: press enter)\n")
        if 'n' in onscreen.lower():
            latlon_inps = raw_input("Enter min_lat/max_lat/min_lon/max_lon:\n")
            lls = latlon_inps.split('/')
            min_lat_v, max_lat_v, min_lon_v, max_lon_v = \
                float(lls[0]), float(lls[1]), float(lls[2]), float(lls[3])
        elif 'q' in onscreen.lower():
            'EXIT!'
            plt.close('all')
            sys.exit()
        else:
            ncuts = plt.ginput(n=2, timeout=-10)
            min_lon_v, min_lat_v = m(ncuts[0][0], ncuts[0][1], inverse=True)
            max_lon_v, max_lat_v = m(ncuts[1][0], ncuts[1][1], inverse=True)
            if 'poles' in plot_mode:
                map_sel = raw_input("which map? (1, 2, 3)")
                if map_sel == '2':
                    min_lon_v, min_lat_v = m2(ncuts[0][0], ncuts[0][1],
                                              inverse=True)
                    max_lon_v, max_lat_v = m2(ncuts[1][0], ncuts[1][1],
                                              inverse=True)
                if map_sel == '3':
                    min_lon_v, min_lat_v = m3(ncuts[0][0], ncuts[0][1],
                                              inverse=True)
                    max_lon_v, max_lat_v = m3(ncuts[1][0], ncuts[1][1],
                                              inverse=True)

        vxyz, vlonlats_arr, vdepths, dist, num_samples = \
            vertical_profile(min_lat_v, max_lat_v, min_lon_v,
                             max_lon_v, min_r, max_r, step_km, eradius)

        # print("===> Querying kd-tree for the vertical profile")
        nps = kd_tree.query(vxyz, k=num_neigh)

        # interpolation
        inter_values = np.zeros([np.shape(vxyz)[0], 4])
        inter_values[:, 0:3] = vxyz[:, 0:3]
        for in1 in range(np.shape(nps)[1]):
            inter_values[in1][3] = \
                np.sum(1.0 / nps[0][in1] * tomo_rd[nps[1][in1]][:, 3]) / \
                np.sum(1.0 / nps[0][in1])
        reshape_interpolate = np.reshape(
            inter_values[:, 3],
            [len(np.arange(min_r, max_r + step_km, step_km)), -1])

        if sigma_gauss_filter:
            reshape_interpolate = gaussian_filter(reshape_interpolate,
                                                  sigma=sigma_gauss_filter,
                                                  order=0)

        plot_cross_section(fig, reshape_interpolate,
                           dist, num_samples,
                           vdepths, vmin, vmax,
                           cmap=cmap,
                           subp=223)

        if 'poles' in plot_mode:
            plt.subplot(2, 2, 4)
            x3, y3 = plot_great_circle(m3, vlonlats_arr)
            plt.annotate('  %s' % path_num, (x3, y3), zorder=100, size=10)
            plt.show()

        plt.subplot(2, 2, 1)
        x3, y3 = plot_great_circle(m, vlonlats_arr)
        plt.annotate('  %s' % path_num, (x3, y3), zorder=100, size=10)

        plt.subplot(2, 2, 2)
        x3, y3 = plot_great_circle(m2, vlonlats_arr)
        plt.annotate('  %s' % path_num, (x3, y3), zorder=100, size=10)
        plt.show()

        print "%s: %s/%s/%s/%s" % (path_num,
                                   min_lat_v, max_lat_v,
                                   min_lon_v, max_lon_v)
        path_num += 1
