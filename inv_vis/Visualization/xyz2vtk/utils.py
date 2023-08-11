#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  utils.py
#   Purpose:   utility codes
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GNU Lesser General Public License, Version 3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as floating_axes
import mpl_toolkits.axisartist.angle_helper as angle_helper
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from mpl_toolkits.basemap.pyproj import Geod
from netCDF4 import Dataset
import numpy as np
from obspy.geodetics import locations2degrees
import os
import pickle

# ---------------- read_tomo_file


def read_tomo_file(tomo_file, dvp_dvs='dvp', netcdf_format=False,
                   eradius=6371.):
    if not netcdf_format:
        tomo_rd = np.loadtxt(tomo_file, dtype=np.float, comments='#')
    else:
        tfio_base = os.path.basename(tomo_file)
        tfio_dir = os.path.dirname(tomo_file)
        tfio_npy = tfio_base.split('.')[0] + '.npy'
        add_npy = os.path.join(tfio_dir, tfio_npy)
        try:
            tomo_rd = np.load(add_npy, "r")
        except Exception, e:
            tomo_nc = Dataset(tomo_file, "r", format='netcdf4')
            tomo_rd = []
            # XXX This is sooo inefficient!
            for dp_i in range(len(tomo_nc['depth'])):
                print "Depth: %s" % tomo_nc['depth'][dp_i]
                for lat_i in range(len(tomo_nc['latitude'])):
                    for lon_i in range(len(tomo_nc['longitude'])):
                        dp = tomo_nc['depth'][dp_i]
                        lat = tomo_nc['latitude'][lat_i]
                        lon = tomo_nc['longitude'][lon_i]
                        x = (eradius - dp) * \
                            np.cos(lat*np.pi/180.) * np.cos(lon*np.pi/180.)
                        y = (eradius - dp) * \
                            np.cos(lat*np.pi/180.) * np.sin(lon*np.pi/180.)
                        z = (eradius - dp) * np.sin(lat*np.pi/180.)
                        tomo_rd.append(
                            [x, y, z,
                             tomo_nc[dvp_dvs][dp_i, lat_i, lon_i]/100.])
            tomo_rd = np.array(tomo_rd)
            np.save(add_npy, tomo_rd)
    return tomo_rd


# ---------------- setup_axes


def setup_axes(fig, rect, theta, radius):
    # PolarAxes.PolarTransform takes radian. However, we want our coordinate
    # system in degree
    tr = Affine2D() + PolarAxes.PolarTransform()

    # Find grid values appropriate for the coordinate (degree).
    # The argument is an approximate number of grids.
    grid_locator1 = angle_helper.LocatorD(2)

    # And also use an appropriate formatter:
    tick_formatter1 = angle_helper.FormatterDMS()

    # set up number of ticks for the r-axis
    grid_locator2 = MaxNLocator(4)

    # the extremes are passed to the function
    grid_helper = floating_axes.GridHelperCurveLinear(
        tr,
        extremes=(theta[0], theta[1], radius[0], radius[1]),
        grid_locator1=grid_locator1,
        grid_locator2=grid_locator2,
        tick_formatter1=tick_formatter1,
        tick_formatter2=None,
    )

    ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
    fig.delaxes(plt.subplot(rect))
    fig.add_subplot(ax1)

    # adjust axis
    # the axis artist lets you call axis with
    # "bottom", "top", "left", "right"
    # # ax1.axis["left"].set_axis_direction("bottom")
    # # ax1.axis["right"].set_axis_direction("top")
    # # ax1.axis["bottom"].set_visible(False)
    # # ax1.axis["top"].set_axis_direction("bottom")
    # # ax1.axis["top"].toggle(ticklabels=True, label=True)
    # # ax1.axis["top"].major_ticklabels.set_axis_direction("top")
    # # ax1.axis["top"].label.set_axis_direction("top")
    # # ax1.axis["left"].label.set_text("R")
    # # ax1.axis["top"].label.set_text(ur"$\theta$")

    ax1.axis["bottom"].set_visible(False)
    ax1.axis["top"].set_visible(False)
    ax1.axis["left"].set_visible(False)
    ax1.axis["right"].set_visible(False)

    # create a parasite axes
    aux_ax = ax1.get_aux_axes(tr)

    aux_ax.patch = ax1.patch  # for aux_ax to have a clip path as in ax
    ax1.patch.zorder = 0.9  # but this has a side effect that the patch is
    # drawn twice, and possibly over some other
    # artists. So, we decrease the zorder a bit to
    # prevent this.

    return ax1, aux_ax

# ---------------- plot_plates


def plot_plates(m, base_path='.', lon360=False, **kwargs):
    for bound in ['ridge', 'transform', 'trench']:
        name, segs = pickle.load(
            open(os.path.join(base_path, 'plates', '%s.pkl' % bound), 'rb'))
        ind_nan, = np.nonzero(np.isnan(segs[:, 0]))
        segs[ind_nan, 0] = 0
        segs[ind_nan, 1] = 0
        if lon360:
            segs[:, 0] = (segs[:, 0] + 360) % 360.0
        x, y = m(segs[:, 0], segs[:, 1])
        x[ind_nan] = np.nan
        y[ind_nan] = np.nan
        dx = np.abs(x[1:] - x[:-1])
        ind_jump, = np.nonzero(dx > 1000000)
        x[ind_jump] = np.nan
        if kwargs:
            m.plot(x, y, '-', zorder=100, **kwargs)
        else:
            m.plot(x, y, '-', zorder=100)

# ---------------- plot_hotspots


def plot_hotspots(m, base_path='.', lon360=False, **kwargs):
    hotspots = pickle.load(
        open(os.path.join(base_path, 'plates', 'hotspots.pkl'), 'rb'))
    if lon360:
        hotspots[:, 0] = (hotspots[:, 0] + 360) % 360.0
    x, y = m(hotspots[:, 0], hotspots[:, 1])
    if kwargs:
        m.scatter(x, y, zorder=400, **kwargs)
    else:
        m.scatter(x, y, zorder=400, c='g', s=30)

# ---------------- plot_stations


def plot_stations(m, base_path='.', lon360=False, **kwargs):
    list_stations = ['120s.txt', 'DEPAS.txt', 'GEOMAR.txt', 'important_loctions.txt', 'INSU.txt', 'LA_LAND.txt', 'LAND.txt', 'SWIR.txt']

    for l_stas in list_stations:
        latlon = np.loadtxt(os.path.join(base_path, 'coords', l_stas), usecols=(0, 1), comments='#')
        x, y = m(latlon[:, 1], latlon[:, 0])
        if kwargs:
            m.scatter(x, y, zorder=200, **kwargs)
        else:
            m.scatter(x, y, marker='v', zorder=200, c='gray', s=25)

# ---------------- depth_profile


def depth_profile(depth, min_lat=-89.9, max_lat=89.9,
                  min_lon=-179.9, max_lon=179.9,
                  eradius=6371., cuts=720):

    lats_degree = np.linspace(min_lat, max_lat, cuts)
    lons_degree = np.linspace(min_lon, max_lon, cuts)
    lonlats_arr = np.vstack([lons_degree, lats_degree]).T

    lats = lats_degree*np.pi/180.
    lons = lons_degree*np.pi/180.

    x = np.array([])
    y = np.array([])
    z = np.array([])
    for i in range(len(lons)):
        x = np.append(x, (eradius - depth) * np.cos(lats) * np.cos(lons[i]))
        y = np.append(y, (eradius - depth) * np.cos(lats) * np.sin(lons[i]))
        z = np.append(z, (eradius - depth) * np.sin(lats))
    xyz = np.vstack([x, y, z]).T

    return lonlats_arr, xyz

# ---------------- vertical_profile


def vertical_profile(min_lat=-90, max_lat=90, min_lon=-180, max_lon=180,
                     min_r=3482, max_r=6371, step_km=20, eradius=6371):

    g = Geod(a=eradius, b=eradius)
    dist = locations2degrees(min_lat, min_lon, max_lat, max_lon)
    print("Distance: %s" % dist)

    num_samples = int(dist)*10
    lonlats = g.npts(min_lon, min_lat, max_lon, max_lat, num_samples)
    lonlats_arr = np.array(lonlats)
    lonlats_arr = np.insert(lonlats_arr, 0, [min_lon, min_lat], axis=0)
    lonlats_arr = np.insert(lonlats_arr, len(lonlats_arr), [max_lon, max_lat],
                            axis=0)

    lats = lonlats_arr[:, 1]*np.pi/180.
    lons = lonlats_arr[:, 0]*np.pi/180.

    x = np.array([])
    y = np.array([])
    z = np.array([])
    depths = np.array([])
    for rv in np.arange(min_r, max_r + step_km, step_km):
        x = np.append(x, rv * np.cos(lats) * np.cos(lons))
        y = np.append(y, rv * np.cos(lats) * np.sin(lons))
        z = np.append(z, rv * np.sin(lats))
        depths = np.append(depths, rv)
    xyz = np.vstack([x, y, z]).T
    return xyz, lonlats_arr, depths, dist, num_samples

# ---------------- plot_cross_section


def plot_cross_section(fig, data, dist, num_samples, depths,
                       vmin, vmax, cmap, subp=223):
    ax1, aux_ax1 = setup_axes(fig, subp,
                              theta=[(90.-dist/2.)*np.pi/180.,
                                     (90.+dist/2.)*np.pi/180.],
                              radius=[3482./6371, 1.])

    theta, rad = np.meshgrid(np.linspace((90+dist/2.)*np.pi/180.,
                                         (90-dist/2.)*np.pi/180.,
                                         num_samples+2),
                             depths/max(depths))

    vplt = aux_ax1.pcolormesh(theta, rad, data,
                              cmap=cmap,
                              vmin=vmin, vmax=vmax)

    cbar = fig.colorbar(vplt,
                        orientation='horizontal', fraction=0.046,
                        ticks=np.linspace(vmin, vmax, 5))
    cbar.ax.set_xticklabels(np.linspace(vmin, vmax, 5)*100)
    # cbar.ax.set_xticklabels(np.linspace(vmin, vmax, 5))

    aux_ax1.scatter((90.+dist/2.)*np.pi/180., 1.,
                    c='red', s=40, edgecolor=None, clip_on=False,
                    zorder=10)
    aux_ax1.scatter((90.-dist/2.)*np.pi/180., 1.,
                    c='blue', s=40, edgecolor=None, clip_on=False,
                    zorder=10)
    ls = ['-', '-', '--', '--', '-', '-']
    i = 0
    for dl in [410, 660, 1000, 2000, 0, 2889]:
        aux_ax1.plot(
            np.linspace((90.-dist/2.), (90.+dist/2.), num_samples)*np.pi/180.,
            [(6371.-dl)/6371]*num_samples, c='k', ls=ls[i])
        i += 1

# ---------------- plot_great_cricle


def plot_great_circle(m, lonlats_arr):
    x1, y1 = m(lonlats_arr[:, 0], lonlats_arr[:, 1])
    x2, y2 = m(lonlats_arr[0, 0], lonlats_arr[0, 1])
    x3, y3 = m(lonlats_arr[-1, 0], lonlats_arr[-1, 1])

    m.scatter(x1, y1, color="k", s=1, zorder=10)
    m.scatter(x2, y2, color="red", s=40, zorder=10)
    m.scatter(x3, y3, color="blue", s=40, zorder=10)
    return x3, y3

# ---------------- custom_div_cmap


def custom_div_cmap(numcolors=11, name='custom_div_cmap',
                    mincol='blue', midcol='white', mid2col='red', maxcol='red'):
    """ Create a custom diverging colormap with three colors

    Default is blue to white to red with 11 colors.  Colors can be specified
    in any way understandable by matplotlib.colors.ColorConverter.to_rgb()
    """

    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list(name=name,
                                             colors=[mincol, midcol, mid2col, maxcol],
                                             N=numcolors)
    return cmap
