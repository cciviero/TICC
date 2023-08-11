#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  utility_codes.py
#   Purpose:   Collection of utility tools
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
from datetime import datetime
import glob

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from scipy.spatial import cKDTree

try:
    from obspy.imaging.beachball import beach as Beach
except Exception as e:
    from obspy.imaging.beachball import Beach
from obspy import read_inventory

import os
import sys
import socket

EARTH_RADIUS = 6371009

# ------------------ no_inp_file_exit ---------------------------


def no_inp_file_exit():
    """
    if no argument is given by the user
    :return:
    """
    print("Serial:   python ffprocessing.py /path/to/input")
    print("Parallel: mpirun -np 4 python ffprocessing.py /path/to/input")
    print("version:  python ffprocessing.py version")
    sys.exit()

# ------------------ print_ffproc_version ---------------------------


def print_ffproc_version():
    """
    ffprocessing version!
    :return:
    """
    print(17*"=")
    print("Version: 0.1.3.0")
    print(17*"=")
    sys.exit()

# ------------------ SphericalNearestNeighbour ---------------------------


class SphericalNearestNeighbour:
    """
    Spherical nearest neighbour queries using scipy's fast kd-tree
    implementation.
    """
    def __init__(self, lat, lon, el_dp):
        cart_data = self.spherical2cartesian(lat, lon, el_dp)
        self.kd_tree = cKDTree(data=cart_data, leafsize=10)

    def query(self, lat, lon, el_dp, k=1):
        points = self.spherical2cartesian(lat, lon, el_dp)
        d, i = self.kd_tree.query(points, k=k)
        return d, i

    def query_pairs(self, maximum_distance):
        return self.kd_tree.query_pairs(maximum_distance)

    def spherical2cartesian(self, lat, lon, el_dp):
        """
        Converts a list of :class:`~obspy.fdsn.download_status.Station`
        objects to an array of shape(len(list), 3) containing x/y/z in meters.
        """
        shape = len(lat)
        r = EARTH_RADIUS + el_dp
        # Convert data from lat/lng to x/y/z.
        colat = 90.0 - lat
        cart_data = np.empty((shape, 3), dtype=np.float64)
        cart_data[:, 0] = r * np.sin(np.deg2rad(colat)) * \
        np.cos(np.deg2rad(lon))
        cart_data[:, 1] = r * np.sin(np.deg2rad(colat)) * \
        np.sin(np.deg2rad(lon))
        cart_data[:, 2] = r * np.cos(np.deg2rad(colat))
        return cart_data

# ------------------ geocen_array ---------------------------


def geocen_array(arr):
    """
    Calculate geocentric latitudes on an array of latitudes
    :param arr:
    :return:
    """
    fac = 0.993305621334896

    # ----------------- First for station:
    colat = 90.0 - arr
    colat[abs(colat) < 1.0e-5] = np.sign(colat[abs(colat) < 1.0e-5])*1.0e-5
    # arg = colat*rpd
    colat *= np.pi/180.
    colat_sin = np.sin(colat)
    colat_sin[colat_sin < 1.0e-30] = 1.0e-30
    # geocen=pi2-atan(fac*cos(arg)/(max(1.0e-30,sin(arg))))
    geocen_colat = np.pi/2. - np.arctan(fac*np.cos(colat)/colat_sin)
    geocen_colat = geocen_colat*180./np.pi
    geocen_lat = 90.0 - geocen_colat

    return geocen_lat

# ------------------ plot_events ---------------------------


def plot_events(all_events, add_save, lat_0=0, lon_0=0, proj='robin'):
    """
    Plot all the events on the selected projection
    :param all_events:
    :param add_save:
    :param lat_0:
    :param lon_0:
    :param proj:
    :return:
    """
    plt.ioff()
    plt.figure(figsize=(30, 20))
    mymap = Basemap(projection=proj, lat_0=round(lat_0, 1), lon_0=round(lon_0, 1))
    mymap.fillcontinents()
    # mymap.drawcoastlines()

    for ev in range(len(all_events.lat)):
        try:
            x, y = mymap(all_events.lon[ev], all_events.lat[ev])
            focmecs = [all_events.mrr[ev],
                       all_events.mtt[ev],
                       all_events.mpp[ev],
                       all_events.mrt[ev],
                       all_events.mrp[ev],
                       all_events.mtp[ev]]

            ax = plt.gca()
            b = Beach(focmecs, xy=(x, y), width=6e5, linewidth=1)
            b.set_zorder(100)
            ax.add_collection(b)
        except Exception as e:

            cprint('utility_codes.py',
                   '[WARNING][PLOT]', bc.red,
                   "Event: %s -- Warning: %s" % (all_events.name[ev], e))
    plt.savefig(os.path.join(add_save, 'all_events.png'), format='png',
                bbox_inches='tight')
    plt.clf()
    plt.ioff()
    plt.close()

# ------------------ plot_ttime ---------------------------


def plot_ttime(all_stations, add_save):
    """
    Plot travel-time of selected phases (all stations)
    :param all_stations:
    :param add_save:
    :return:
    """
    plt.ioff()
    plt.figure(figsize=(20, 10))
    plt.scatter(all_stations.dist, all_stations.tt_ph,
                marker='o', color='k', edgecolors='k',
                alpha=0.5)
    plt.xlabel('Distance (deg)', size=24, weight='bold')
    plt.ylabel('Travel Time (sec)', size=24, weight='bold')
    plt.xticks(size=18, weight='bold')
    plt.yticks(size=18, weight='bold')

    if not os.path.isdir(add_save):
        os.makedirs(add_save)
    plt.savefig(os.path.join(add_save, 'tt_phase.png'),
                format='png', bbox_inches='tight')
    plt.clf()
    plt.ioff()
    plt.close()

# ------------------ plot_event_stations ---------------------------


def plot_event_stations(all_events, all_stations, ev, inp_ffproc, proj='aeqd'):
    """
    Plot one event and all found stations
    :param all_events:
    :param all_stations:
    :param ev:
    :param inp_ffproc:
    :param proj:
    :return:
    """
    ev_lat = all_events.lat[ev]
    ev_lon = all_events.lon[ev]

    plt.ioff()
    plt.figure(figsize=(10, 10))
    
    try:
        mymap = Basemap(projection=proj, lat_0=round(ev_lat, 1), lon_0=round(ev_lon, 1))
    except Exception as e:
        mymap = Basemap(projection=proj, lat_0=round(ev_lat, 2), lon_0=round(ev_lon, 2))

    # mymap.fillcontinents()
    mymap.drawcoastlines()

    x, y = mymap(all_stations.lon, all_stations.lat)
    mymap.scatter(x, y, c='b', edgecolor='none', marker='v', zorder=20, s=100)
    # creating epicentral-distance-circles
    radii = np.arange(30., 180., 30.) * 111.
    for r in radii:
        equi(mymap, ev_lon, ev_lat, r,
             lw=1., color='k', linestyle='--')
    equi(mymap, ev_lon, ev_lat, 90.*111,
         lw=1., color='k', linestyle='-')
    equi(mymap, ev_lon, ev_lat, inp_ffproc.ph_min_epi*111,
         lw=1., color='r', linestyle='-')
    equi(mymap, ev_lon, ev_lat, inp_ffproc.ph_max_epi*111,
         lw=1., color='r', linestyle='-')

    x, y = mymap(ev_lon, ev_lat)
    focmecs = [all_events.mrr[ev],
               all_events.mtt[ev],
               all_events.mpp[ev],
               all_events.mrt[ev],
               all_events.mrp[ev],
               all_events.mtp[ev]]

    try:
        ax = plt.gca()
        b = Beach(focmecs, xy=(x, y), width=int(2e6), linewidth=1, facecolor='r')
        b.set_zorder(50)
        ax.add_collection(b)
    except Exception as e:
        mymap.scatter(x, y, c='r', edgecolor='none', marker='o', zorder=20, s=100)
    plt.savefig(os.path.join(inp_ffproc.output_dir,
                             all_events.name[ev],
                             'all_stations_%s.png' % all_events.name[ev]),
                format='png', bbox_inches='tight')
    plt.clf()
    plt.ioff()
    plt.close()

# ------------------ plot_cut_wave ---------------------------


def plot_cut_wave(all_stations, add_save, fig_name, waveforms, ev_datetime):
    """
    Plot waveforms (synthetic/real/...)
    :param all_stations:
    :param add_save:
    :param fig_name:
    :param waveforms:
    :param ev_datetime:
    :return:
    """
    plt.ioff()
    fig, ax1 = plt.subplots()
    fig.set_size_inches(18, 15)

    sta_collector = []
    y_ticks = []
    y_label = []
    preset_arr = []
    for sta in range(len(waveforms)):
        tr = waveforms[sta]
        y_ticks.append(all_stations.dist[sta])
        y_label.append(all_stations.name[sta].split('.')[1])

        preset = all_stations.tt_ph[sta] - (tr.stats.starttime - ev_datetime)
        preset_arr.append(preset)
        try:
            ax1.plot(
                np.linspace(-preset,
                            (tr.stats.npts-1)/tr.stats.sampling_rate - preset,
                            tr.stats.npts),
                tr.data/max(abs(tr.data)) + all_stations.dist[sta],
                color='k', lw=1)
        except Exception as exp:
            import idpb; ipdb.set_trace()
        sta_collector.append([all_stations.dist[sta],
                              tr.data / max(abs(tr.data)),
                              preset])
    ax1.axvline(x=0, lw=2, linestyle='--', color='r')
    ax1.set_xlabel('Time (sec)', size=24, weight='bold')
    ax1.set_ylabel('Epicentral Distance (deg)', size=24, weight='bold')
    ax1.set_xticks(ax1.get_xticks())
    ax1.set_xticklabels(ax1.get_xticks(), size=18, weight='bold')
    ax1.set_yticks(ax1.get_yticks())
    ax1.set_yticklabels(ax1.get_yticks(), size=18, weight='bold')

    if len(sta_collector) < 1:
        cprint('utility_codes.py',
               '[OUTPUT][PLOT]', bc.red,
               '[OUTPUT] No phases/stations for this event!')

    preset_arr = np.array(preset_arr)
    preset_max = np.max(preset_arr)
    plt.xlim(xmin=-1*preset_max)
    plt.savefig(os.path.join(add_save, fig_name + '.png'),
                format='png', bbox_inches='tight')
    plt.clf()
    plt.ioff()
    plt.close()

    # ---------------- pcolormesh figure
    plt.ioff()
    fig, ax1 = plt.subplots()
    fig.set_size_inches(18, 15)

    sta_collector = sorted(sta_collector, key=lambda stas: stas[0])
    # making the x axis
    x = np.linspace(-1*preset_max,
                    (tr.stats.npts-1)/tr.stats.sampling_rate - preset_max,
                    tr.stats.npts)
    # creating the y axis
    y = []
    sta_collector_data = []
    for sta in range(len(sta_collector)):
        y.append(sta_collector[sta][0])
        sta_collector_data.append(sta_collector[sta][1])
    y = np.array(y)
    z = np.array(sta_collector_data)
    x, y = np.meshgrid(x, y)

    len_stripe_half = (np.max(y) - np.min(y))/len(z)/2.
    if len_stripe_half == 0:
        len_stripe_half = 0.5
    for i in range(len(x)):
        # to give the band a fixed width
        y1 = y[i] - len_stripe_half
        y2 = y[i] + len_stripe_half
        y_axis = np.array((y1, y2), dtype=float)

        # to fill the 'gaps' with white areas
        l = [0] * len(z[i])
        z_axis = np.array((z[i], l), dtype=float)
        offset = -1*sta_collector[i][2] - x[i][0]
        pcp = ax1.pcolormesh(x[i] + offset, y_axis, z_axis,
                             cmap='seismic_r',
                             vmin=-1., vmax=1.)

    ax1.axvline(x=0, lw=2, linestyle='--', color='k')
    ax1.set_xlabel('Time (sec)', size=24, weight='bold')
    ax1.set_ylabel('Epicentral Distance (deg)', size=24, weight='bold')
    ax1.set_xticks(ax1.get_xticks().tolist())
    ax1.set_yticks(ax1.get_yticks().tolist())
    ax1.set_xticklabels(ax1.get_xticks(), size=18, weight='bold')
    ax1.set_yticklabels(ax1.get_yticks(), size=18, weight='bold')

    if len(sta_collector) < 1:
        cprint('utilit_codes.py', '[OUTPUT][PLOT]',
               bc.lred, 'No phases/stations for this event!')

    plt.colorbar(pcp)
    plt.xlim(xmin=-1*preset_max)
    plt.savefig(os.path.join(add_save, "pcolor_" + fig_name + ".png"),
                format='png', bbox_inches='tight')
    plt.clf()
    plt.ioff()
    plt.close()

# ------------------- plot_filt_compare_groups --------------------


def plot_filt_compare_groups(all_stations, realBP, synBP,
                             cut_realBP, cut_synBP,
                             inp_ffproc, add_save, fig_name):
    """
    Plot both real and synthetic waveforms over each other after correcting
    for the time shift
    :param all_stations:
    :param realBP:
    :param synBP:
    :param cut_realBP:
    :param cut_synBP:
    :param inp_ffproc:
    :param add_save:
    :param fig_name:
    :return:
    """
    plt.ioff()
    for band in range(inp_ffproc.check_min_bands, inp_ffproc.check_max_bands):
        fig = plt.figure(frameon=False)
        fig.set_size_inches(20, 15)

        p_counter = 0
        fig_plt_counter = 0
        loc_tick = []
        txt_tick = []
        for sta in range(len(all_stations.name)):
            plt_opacity = 1.
            if not all_stations.cc_step2[band, sta] >= inp_ffproc.check_min_cc:
                plt_opacity = 0.3
            if all_stations.clip_step2[band, sta] > 0:
                plt_opacity = 0.3
            real_start_time = - inp_ffproc.ph_preset \
                              - all_stations.final_time_shift[band, sta]
            len_real = (len(realBP[:, sta, band])-1)/all_stations.srate[sta]
            pp1, = plt.plot(
                np.linspace(real_start_time,
                            real_start_time + len_real,
                            len(realBP[:, sta, band])),
                realBP[:, sta, band]/max(abs(realBP[:, sta, band])) +
                p_counter*2,
                color='k', lw=3, alpha=plt_opacity)

            # ------------------ REAL (NOT SHIFTED)
            real_start_time = - inp_ffproc.ph_preset
            pp2, = plt.plot(
                np.linspace(real_start_time,
                            real_start_time + len_real,
                            len(realBP[:, sta, band])),
                realBP[:, sta, band]/max(abs(realBP[:, sta, band])) +
                p_counter*2,
                color='k', ls='--', lw=2, alpha=plt_opacity)

            # ------------------ SYN
            syn_start_time = - inp_ffproc.ph_preset
            len_syn = (len(synBP[:, sta, band])-1)/all_stations.srate[sta]
            pp3, = plt.plot(
                np.linspace(syn_start_time,
                            syn_start_time + len_syn,
                            len(synBP[:, sta,  band])),
                synBP[:, sta, band]/max(abs(synBP[:, sta, band])) +
                p_counter*2,
                color='b', lw=3, alpha=plt_opacity)

            plt.text(syn_start_time, p_counter*2+1.0,
                     '%s' % all_stations.name[sta],
                     fontsize=8.5, weight='bold', color='red')
            plt.text(syn_start_time, p_counter*2+0.80,
                     'CC: %4.2f -- Clip: %i'
                     % (all_stations.cc_step2[band, sta],
                        all_stations.clip_step2[band, sta]),
                     fontsize=8, weight='bold')
            plt.text(syn_start_time, p_counter*2+0.60,
                     'dT: %4.2f -- A: %4.2f'
                     % (all_stations.final_time_shift[band, sta],
                        all_stations.final_amp[band, sta]),
                     fontsize=8, weight='bold')
            plt.text(syn_start_time, p_counter*2+0.40,
                     'SNR: %4.2f'
                     % (all_stations.snr[band, sta]),
                     fontsize=8, weight='bold')

            # ------------------ CUT SYN
            len_cut_syn = (len(cut_synBP[band][sta])-1)/all_stations.srate[sta]
            syn_start_time = - inp_ffproc.ph_preset
            syn_start_time += \
                (all_stations.cut_sample_syn[band][sta]-1) / \
                inp_ffproc.ph_sampling_rate
            pp4, = plt.plot(
                np.linspace(syn_start_time,
                            syn_start_time + len_cut_syn,
                            len(cut_synBP[band][sta])),
                cut_synBP[band][sta]/max(abs(synBP[:, sta, band])) +
                p_counter*2,
                color='r', lw=3, alpha=plt_opacity)

            loc_tick.append(p_counter*2)
            txt_tick.append(round(all_stations.dist[sta], 2))
            p_counter += 1
            if np.mod(p_counter, 10) < 1 or \
                    (sta == (len(all_stations.name) - 1)):
                plt.title('Band: %s' % (band+1), size=24, weight='bold')
                plt.xlabel('Time (sec)', size=24, weight='bold')
                plt.ylabel('Epicentral Distance (deg)', size=24, weight='bold')
                plt.xticks(size=18, weight='bold')
                plt.yticks(loc_tick, txt_tick, size=18, weight='bold')
                plt.xlim(-inp_ffproc.ph_preset)
                plt.ylim(-1, p_counter*2)
                plt.axvline(0, 0, 1, c='black', ls='--', lw=2.5)
                plt.legend([pp1, pp2, pp4],
                           ['real (shifted)',
                            'real (not-shifted)',
                            'syn (measured window)'])
                plt.savefig(
                    os.path.join(add_save,
                                 '%04i_%s_' % (fig_plt_counter, band+1) +
                                 fig_name + ".png"),
                    format='png', bbox_inches='tight')

                plt.clf()
                plt.ioff()
                plt.close()
                fig_plt_counter += p_counter
                p_counter = 0
                # set everything back after 10 smgr plotted in one figure
                loc_tick = []
                txt_tick = []
                fig = plt.figure(frameon=False)
                fig.set_size_inches(20, 15)
                plt.ioff()
            if sta == (len(all_stations.name) - 1):
                plt.clf()
                plt.close()

# ------------------- plot_measured_glob --------------------


def plot_measured_glob(all_stations, all_events, ev, inp_ffproc,
                       band=0, proj='aeqd', center_lat=False,
                       center_lon=False, vmin=-3, vmax=3):
    """
    projection of measured travel-time and cross-correlation factor
    over the glob
    :param all_stations:
    :param all_events:
    :param ev:
    :param inp_ffproc:
    :param band:
    :param proj:
    :param center_lat:
    :param center_lon:
    :param vmin:
    :param vmax:
    :return:
    """
    # Some useful latitude and longitudes
    # US:
    # center_lat = 38.5; center_lon = -115.

    # ======================= Time Shift
    plt.ioff()
    plt.figure(figsize=(10, 10))
    if not center_lat or center_lon:
        center_lat = all_events.lat[ev]
        center_lon = all_events.lon[ev]
    mymap = Basemap(projection=proj, lat_0=round(center_lat, 1), lon_0=round(center_lon, 1))
    mymap.drawcoastlines(color='black', linewidth=1)

    # In case that only good stations should be plotted:
    # sel_stas = \
    #     (all_stations.cc_step2[band, :] >= inp_ffproc.check_min_cc) * \
    #     (all_stations.clip_step2[band, :] < 1)

    sel_stas = (all_stations.cc_step2[band, :] >= -100)
    lats = all_stations.lat[sel_stas]
    lons = all_stations.lon[sel_stas]

    time_shift = all_stations.final_time_shift[band, sel_stas]
    x, y = mymap(lons, lats)
    # removing the median from the measured travel times
    cprint('utility_codes.py', '[PLOTTING]', bc.yellow,
           'Removing the median of the measurements!')
    c_plot = time_shift - np.median(time_shift)
    mymap.scatter(x, y, c=c_plot.astype(np.float),
                  edgecolor='none',
                  zorder=20, s=20, vmin=vmin, vmax=vmax)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=24)

    # creating circles with the same epicentral distances
    radii = np.arange(30., 180., 30.) * 111.
    for r in radii:
        equi(mymap, all_events.lon[ev], all_events.lat[ev], r,
             lw=1., color='k',  linestyle='--')
    equi(mymap, all_events.lon[ev], all_events.lat[ev], 90.*111,
         lw=1., color='k', linestyle='-')

    # plotting the event
    x_b, y_b = mymap(all_events.lon[ev], all_events.lat[ev])
    focmecs = [all_events.mrr[ev],
               all_events.mtt[ev],
               all_events.mpp[ev],
               all_events.mrt[ev],
               all_events.mrp[ev],
               all_events.mtp[ev]]

    try:
        ax = plt.gca()
        b = Beach(focmecs, xy=(x_b, y_b), width=int(2e6), linewidth=1)
        b.set_zorder(50)
        ax.add_collection(b)
    except Exception as e:
        mymap.scatter(x_b, y_b, c='r',
                      edgecolor='none',
                      zorder=20, s=20)
    plt.title('dT', weight='bold', size=18)

    plt.savefig(os.path.join(inp_ffproc.output_dir, all_events.name[ev],
                             'dt_band%02i.png' % (band+1)),
                format='png', bbox_inches='tight')
    plt.clf()
    plt.ioff()
    plt.close()

    # ======================= CC
    plt.ioff()
    plt.figure(figsize=(10, 10))
    mymap = Basemap(projection=proj, lat_0=round(center_lat, 1), lon_0=round(center_lon, 1))
    mymap.drawcoastlines(color='black', linewidth=1)

    lats = all_stations.lat
    lons = all_stations.lon
    cc_calc = all_stations.cc_step2[band, :]

    x, y = mymap(lons, lats)
    mymap.scatter(x, y, c=cc_calc,
                  edgecolor='none',
                  zorder=20, s=20, vmin=0, vmax=1.)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=24)

    # creating circles with the same epicentral distances
    radii = np.arange(30., 180., 30.) * 111.
    for r in radii:
        equi(mymap, all_events.lon[ev], all_events.lat[ev], r,
             lw=1.,  color='k',  linestyle='--')
    equi(mymap, all_events.lon[ev], all_events.lat[ev], 90.*111,
         lw=1., color='k', linestyle='-')

    try:
        ax = plt.gca()
        b = Beach(focmecs, xy=(x_b, y_b), width=int(2e6), linewidth=1)
        b.set_zorder(50)
        ax.add_collection(b)
    except Exception as e:
        mymap.scatter(x_b, y_b, c='r',
                      edgecolor='none',
                      zorder=20, s=20)
    plt.title('CC factor', weight='bold', size=18)

    plt.savefig(os.path.join(inp_ffproc.output_dir, all_events.name[ev],
                             'cc_band%02i.png' % (band+1)),
                format='png', bbox_inches='tight')
    plt.clf()
    plt.ioff()
    plt.close()

    # ======================= AMPLITUDE
    plt.ioff()
    plt.figure(figsize=(10, 10))
    mymap = Basemap(projection=proj, lat_0=round(center_lat, 1), lon_0=round(center_lon, 1))
    mymap.drawcoastlines(color='black', linewidth=1)

    lats = all_stations.lat
    lons = all_stations.lon
    amp_calc = all_stations.final_amp[band, :]

    x, y = mymap(lons, lats)
    mymap.scatter(x, y, c=amp_calc,
                  edgecolor='none',
                  zorder=20, s=20, vmin=0, vmax=1)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=24)

    # creating circles with the same epicentral distances
    radii = np.arange(30., 180., 30.) * 111.
    for r in radii:
        equi(mymap, all_events.lon[ev], all_events.lat[ev], r,
             lw=1.,  color='k',  linestyle='--')
    equi(mymap, all_events.lon[ev], all_events.lat[ev], 90.*111,
         lw=1., color='k', linestyle='-')

    try:
        ax = plt.gca()
        b = Beach(focmecs, xy=(x_b, y_b), width=int(2e6), linewidth=1)
        b.set_zorder(50)
        ax.add_collection(b)
    except Exception as e:
        mymap.scatter(x_b, y_b, c='r',
                      edgecolor='none',
                      zorder=20, s=20)
    plt.title('Amplitude', weight='bold', size=18)

    plt.savefig(os.path.join(inp_ffproc.output_dir, all_events.name[ev],
                             'amp_band%02i.png' % (band+1)),
                format='png', bbox_inches='tight')
    plt.clf()
    plt.ioff()
    plt.close()
 
# ------------------- plot_spec_gabor_filter_compare --------------------


def plot_spec_gabor_filter_compare(all_stations, real_waveforms, realBP, synBP,
                           event_time, cut_realBP, cut_synBP,
                           inp_ffproc, add_save, fig_name):
    """
    plot the result of the gabor filter in frequency domain beside the filtered
    smgr/syn
    :param all_stations:
    :param realBP:
    :param synBP:
    :param cut_realBP:
    :param cut_synBP:
    :param inp_ffproc:
    :param add_save:
    :param fig_name:
    :return:
    """

    plt.ioff()
    for band in range(inp_ffproc.check_min_bands, inp_ffproc.check_max_bands):
        progress_bar(band, inp_ffproc.check_max_bands)
        f, (ax1, ax2) = plt.subplots(1, 2, figsize=(30, 15), sharey=True)

        p_counter = 0
        fig_plt_counter = 0
        sta_collector = []
        ph_shifts = []
        # for epicentral distance on yaxis
        loc_tick = []
        txt_tick = []
        list_counter = 0
        for sta in range(len(all_stations.name)):
            ph_shifts.append([all_stations.dist[sta], list_counter])
            list_counter += 1

        plt_counter = 0
        for sta in range(len(all_stations.name)):
            plt_opacity = 1.
            if not all_stations.cc_step2[band, sta] >= inp_ffproc.check_min_cc:
                plt_opacity = 0.3
            if all_stations.clip_step2[band, sta] > 0:
                plt_opacity = 0.3
            # ------------------ REAL (SHIFTED)
            real_start_time = - inp_ffproc.ph_preset \
                              - all_stations.final_time_shift[band, sta]
            len_real = (len(realBP[:, sta, band]) - 1) / all_stations.srate[sta]
            pp1, = ax1.plot(np.linspace(real_start_time,
                            real_start_time + len_real,
                            len(realBP[:, sta, band])),
                            realBP[:, sta, band] / max(abs(realBP[:, sta, band])) + p_counter * 2,
                            color='k', lw=3, label='real (shifted)',
                            alpha=plt_opacity)

            # ------------------ Create spectrogram out of shifted real signal
            # original real singal to spectrogram
            freqs_real = np.fft.rfftfreq(real_waveforms[sta].stats.npts) / real_waveforms[sta].stats.delta
            fft_real = np.fft.rfft(real_waveforms[sta].data)
            sp1, = ax2.semilogx(freqs_real, abs(fft_real)/max(abs(fft_real)) + p_counter * 2,
                            color=(57/250., 106/250., 177/250.), lw=4, label='original trace', alpha=0.8)

            # cut out noise window from real signal
            noise_window = real_waveforms[sta].copy()

            window = event_time + all_stations.tt_ph[sta] - real_waveforms[sta].stats.starttime
            # XXX fixed values for cutting the noise window depending on the first arrival
            noise_window.trim(real_waveforms[sta].stats.starttime + 10,
                              real_waveforms[sta].stats.starttime + window - 10)
            freqs_noise = np.fft.rfftfreq(noise_window.stats.npts) / noise_window.stats.delta
            fft_noise= np.fft.rfft(noise_window.data)
            sp2, = ax2.semilogx(freqs_noise, abs(fft_noise) / max(abs(fft_noise)) + p_counter * 2,
                                color=(204/250., 37/250., 41/250.), lw=3, label='noise window')

            # bandpass filtered signal
            freqs_band = np.fft.rfftfreq(len(realBP[:, sta, band])) / real_waveforms[sta].stats.delta
            fft_band = np.fft.rfft(realBP[:, sta, band])
            sp3, = ax2.semilogx(freqs_band, abs(fft_band)/max(abs(fft_band)) + p_counter * 2,
                            color=(62/250., 150/250., 81/250.), lw=3, label='filtred trace')

            # gabor filter
            freqs_gabor = np.fft.rfftfreq(len(inp_ffproc.lgb_filt_IR[:, band])) / real_waveforms[sta].stats.delta
            fft_gabor = np.fft.rfft(inp_ffproc.lgb_filt_IR[:, band])
            sp4, = ax2.semilogx(freqs_gabor, abs(fft_gabor)/max(abs(fft_gabor)) + p_counter * 2,
                            color=(107/250., 76/250., 154/250.), lw=4, label='gabor filter')

            # ------------------ Create CEPSTRUM out of shifted real signal
            # #import ipdb;
            # #ipdb.set_trace()
            #
            # # original real signal to cepstrum
            # spectrum_tr = np.fft.rfft(real_waveforms[sta].data)
            # log_spectrum = np.log(spectrum_tr)
            # cepstrum = np.fft.irfft(log_spectrum)
            # cp1, = ax3.plot(np.linspace(-inp_ffproc.ph_preset,
            #         (real_waveforms[sta].stats.npts - 1) /
            #                             real_waveforms[sta].stats.sampling_rate - inp_ffproc.ph_preset,
            #                             real_waveforms[sta].stats.npts)[:-1],
            #                 abs(cepstrum) + p_counter * 2,
            #                 color=(57/250., 106/250., 177/250.), lw=3)
            #
            # # bandpass filtered signal to cepstrum
            # spectrum_filter = np.fft.rfft(realBP[:, sta, band])
            # log_spectrum_filter = np.log(spectrum_filter)
            # cepstrum_filter = np. .irfft(log_spectrum_filter)
            # cp2, = ax3.plot(np.linspace(real_start_time,
            #                 real_start_time + len_real,
            #                 len(realBP[:, sta, band]))[:-1],
            #                 abs(cepstrum_filter) + p_counter * 2,
            #                 color=(62 / 250., 150 / 250., 81 / 250.), lw=2, alpha=0.8)

            # ------------------ REAL (NOT SHIFTED)
            real_start_time = - inp_ffproc.ph_preset
            pp2, = ax1.plot(
                np.linspace(real_start_time,
                            real_start_time + len_real,
                            len(realBP[:, sta, band])),
                realBP[:, sta, band] / max(abs(realBP[:, sta, band])) +
                p_counter * 2,
                color='k', ls='--', lw=2, label='real (not-shifted)',
                alpha=plt_opacity)

            # ------------------ SYN
            syn_start_time = - inp_ffproc.ph_preset
            len_syn = (len(synBP[:, sta, band]) - 1) / all_stations.srate[sta]
            pp3, = ax1.plot(
                np.linspace(syn_start_time,
                            syn_start_time + len_syn,
                            len(synBP[:, sta, band])),
                synBP[:, sta, band] / max(abs(synBP[:, sta, band])) +
                p_counter * 2,
                color='b', lw=3, alpha=plt_opacity)

            ax1.text(syn_start_time, p_counter*2+1.0,
                     '%s' % all_stations.name[sta],
                     fontsize=8.5, weight='bold')
            ax1.text(syn_start_time, p_counter*2+0.80,
                     'CC: %4.2f -- Clip: %i'
                     % (all_stations.cc_step2[band, sta],
                        all_stations.clip_step2[band, sta]),
                     fontsize=8)
            ax1.text(syn_start_time, p_counter*2+0.60,
                     'dT: %4.2f -- A: %4.2f'
                     % (all_stations.final_time_shift[band, sta],
                        all_stations.final_amp[band, sta]),
                     fontsize=8, weight='bold', color='red')
            ax1.text(syn_start_time, p_counter*2+0.40,
                     'SNR: %4.2f'
                     % (all_stations.snr[band, sta]),
                     fontsize=8)

            # ------------------ CUT SYN
            len_cut_syn = (len(cut_synBP[band][sta]) - 1) / all_stations.srate[sta]
            syn_start_time = - inp_ffproc.ph_preset
            syn_start_time += \
                (all_stations.cut_sample_syn[band][sta] - 1) / \
                inp_ffproc.ph_sampling_rate
            pp4, = ax1.plot(
                np.linspace(syn_start_time,
                            syn_start_time + len_cut_syn,
                            len(cut_synBP[band][sta])),
                cut_synBP[band][sta] / max(abs(synBP[:, sta, band])) +
                p_counter * 2,
                color='r', lw=3, label='measured window', alpha=plt_opacity)

            sta_collector.append([all_stations.dist[sta],
                                  all_stations.tt_ph[sta]])
            loc_tick.append(plt_counter * 2)

            txt_tick.append(round(all_stations.dist[sta], 2))
            plt_counter += 1
            p_counter += 1

            if np.mod(plt_counter, 10) < 1 or \
                    (sta == (len(all_stations.name) - 1)):
                # ax1 settings
                ax1.set_title('Band: %s' % (band + 1), size=24, weight='bold')
                ax1.set_xlabel('Time (sec)', size=24, weight='bold')
                ax1.set_ylabel('Epicentral Distance (deg)', size=24, weight='bold')
                x1_tick_ticks = ax1.get_xticks()
                ax1.set_xticks(x1_tick_ticks)
                ax1.set_xticklabels(x1_tick_ticks, size=18, weight='bold')
                # creating y ticks based on previously created ph_shifts
                ax1.set_yticks(loc_tick)
                ax1.set_yticklabels(txt_tick, size=18, weight='bold')
                # XXX this is hard coded
                ax1.set_xlim(- inp_ffproc.ph_preset, 100)
                ax1.set_ylim(-1, p_counter * 2)

                ax1.axvline(0, 0, 1, c='black', ls='--', lw=2.5)
                ax1.legend([pp1, pp2, pp4], ['real (shifted)', 'real (not-shifted)',
                                             'syn (measured window)'])

                # ax2 settings
                ax2.set_title('Spectrum (normalised gain)', size=24, weight='bold')
                x2_tick_ticks = ax2.get_xticks()
                ax2.set_xticks(x2_tick_ticks)
                ax2.set_xticklabels(x2_tick_ticks, size=18, weight='bold')
                ax2.set_xlim(1/150., 5)
                ax2.set_xlabel('Frequency (Hz)', size=24, weight='bold')

                ax2.legend([sp1, sp2, sp3, sp4],
                           ['original trace', 'noise window', 'filtered trace',
                            'gabor filter'])


                plt.savefig(
                    os.path.join(add_save,
                                 '%04i_%s_' % (fig_plt_counter, band + 1) +
                                 fig_name + ".png"),
                    format='png', bbox_inches='tight')

                plt.clf()
                plt.ioff()
                plt.close()
                p_counter = 0
                fig_plt_counter += plt_counter
                # set everything back after 10 smgr plotted in one figure
                plt_counter = 0
                loc_tick = []
                txt_tick = []
                f, (ax1, ax2) = plt.subplots(1, 2, figsize=(30, 15), sharey=True)
                #f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(30, 15), sharey=True)

                plt.ioff()
            if sta == (len(all_stations.name) - 1):
                plt.clf()
                plt.ioff()

# ------------------- shoot --------------------


def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    """
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.

    EPS = 0.00000000005
    if np.abs(np.cos(glat1)) < EPS and not np.abs(np.sin(faz)) < EPS:
        print("Only N-S courses are meaningful, starting at a pole!")

    a = 6378.13/1.852
    f = 1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if cf == 0:
        b = 0.
    else:
        b = 2. * np.arctan2(tu, cf)

    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while np.abs(y - c) > EPS:
        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
             d / 4. - cz) * sy * d + tu

    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi

    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)

    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi

    return glon2, glat2, baz

# ------------------- equi --------------------


def equi(m, centerlon, centerlat, radius, *args, **kwargs):
    """
    Plotting the lines using shoot function
    :param m:
    :param centerlon:
    :param centerlat:
    :param radius:
    :param args:
    :param kwargs:
    :return:
    """
    glon1 = centerlon
    glat1 = centerlat
    X = []
    Y = []
    for azimuth in range(0, 360):
        glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])
    X, Y = m(X, Y)
    plt.plot(X, Y, **kwargs)

# ------------------- progress_bar --------------------


def progress_bar(curr_indx, total_number):
    """
    Showing the progress in a loop
    :param curr_indx:
    :param total_number:
    :return:
    """
    sys.stdout.write('\r')
    sys.stdout.write("[%-100s] %d%%"
                     % ('='*int(100.*(curr_indx+1)/total_number),
                        100.*(curr_indx+1)/total_number))
    sys.stdout.flush()

# ------------------- bc --------------------


class bc:
    lgrey = '\033[1;90m'
    grey = '\033[90m'     # broing information
    yellow = '\033[93m'   # FYI
    orange = '\033[0;33m'  # Warning

    lred = '\033[1;31m'  # there is smoke
    red = '\033[91m'    # fire!
    dred = '\033[2;31m'  # Everything is on fire

    lblue = '\033[1;34m'
    blue = '\033[94m'
    dblue = '\033[2;34m'

    lgreen = '\033[1;32m'  # all is normal
    green = '\033[92m'    # something else
    dgreen = '\033[2;32m'  # even more interesting

    lmagenta = '\033[1;35m'
    magenta = '\033[95m'  # for title
    dmagenta = '\033[2;35m'

    cyan = '\033[96m'    # system time
    white = '\033[97m'   # final time

    black = '\033[0;30m'

    end = '\033[0m'
    bold = '\033[1m'
    under = '\033[4m'

# ------------------- get_time --------------------


def get_time():
    time = datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S')
    return time

# ------------------- cprint --------------------


def cprint(code, type_info, bc_color, text):
    """
    Instead the logging function which does not work with parallel mpi4py
    codes
    :param code:
    :param type_info:
    :param bc_color:
    :param text:
    :return:
    """
    ho_nam = socket.gethostname().split('.')[0]

    print(bc.green + get_time() + bc.end, \
          bc.magenta + ho_nam + bc.end, \
          bc.cyan + code + bc.end, \
          bc.bold + bc.grey + type_info + bc.end, \
          bc_color + text + bc.end)

# ------------------- check_exitus_file --------------------


def check_exitus_file(ev_name, inp_ffproc):
    """
    Check if the event was already processed to the end
    :param ev_name:
    :param inp_ffproc:
    :return:
    """
    file_exitus = glob.glob(os.path.join(inp_ffproc.output_dir,
                                         ev_name,
                                         "exitus_%s.txt" % ev_name))
    if file_exitus:
        cprint('utility_codes.py', '[EVENT]', bc.dmagenta,
               '%s exists. CONTINUE!' % ev_name)
        return True
    else:
        return False

# ------------------- exitus_file --------------------


def exitus_file(all_events, ev, all_stfs, inp_ffproc, time):
    """
    writes a file with processing time information at the end of every event
    :param all_events:
    :param ev:
    :param inp_ffproc:
    :param time:
    :return:
    """
    ev_name = all_events.name[ev]
    with open(os.path.join(inp_ffproc.output_dir, ev_name,
                           "exitus_%s.txt" % ev_name), 'w') as f:
        f.write('%s min PROCESSING TIME FOR EVENT: %s \n' % (str(time), ev_name))

        if inp_ffproc.int_syn:
            f.write('\n ATTENTION!!! \n int_syn = %s \n' % inp_ffproc.int_syn)
        if inp_ffproc.int_real:
            f.write('\n ATTENTION!!! \n int_real = %s \n' % inp_ffproc.int_real)

        f.write('\n--------------------------------------------------------------\n')
        f.write('Input summary from your input_ffprocessing.ini file:\n')
        f.write('--------------------------------------------------------------\n\n')
        f.write('ph_phase = %s\n' % inp_ffproc.ph_phase)
        f.write('ph_channel = %s\n' % inp_ffproc.ph_channel)
        f.write('ph_min_epi = %s \n' % inp_ffproc.ph_min_epi)
        f.write('ph_max_epi = %s \n' % inp_ffproc.ph_max_epi)
        f.write('ph_bg = %s\n' % inp_ffproc.ph_bg)
        f.write('ph_sampling_rate = %s \n' % inp_ffproc.ph_sampling_rate)
        f.write('ph_preset = %s\n' % inp_ffproc.ph_preset)
        f.write('ph_offset = %s\n' % inp_ffproc.ph_offset)
        f.write('ph_static_preset = %s\n' % inp_ffproc.ph_static_preset)
        f.write('ph_static_offset = %s\n' % inp_ffproc.ph_static_offset)
        f.write('mmeant_clip_time1 = %s\n' % inp_ffproc.mmeant_clip_time1)
        f.write('mmeant_clip_time2 = %s\n' % inp_ffproc.mmeant_clip_time2)
        f.write('real mode = %s\n' % inp_ffproc.real_mode)
        f.write('syn mode = %s\n' % inp_ffproc.syn_mode)
        f.write('evproc mode = %s\n' % inp_ffproc.evproc_mode)
        f.write('stf mode = %s\n' % inp_ffproc.stf_mode)

        f.write('\n--------------------------------------------------------------\n')
        f.write('Gabor filter:\n')
        f.write('--------------------------------------------------------------\n\n')
        f.write('lgb_filt_pmax = %s \n' % inp_ffproc.lgb_filt_pmax)
        f.write('lgb_filt_nscale = %s \n' % inp_ffproc.lgb_filt_nscale)
        f.write('lgb_filt_fmult = %s \n' % inp_ffproc.lgb_filt_fmult)
        f.write('lgb_filt_sigmaIfc = %s \n' % inp_ffproc.lgb_filt_sigmaIfc)
        f.write('lgb_filt_npad = %s \n' % inp_ffproc.lgb_filt_npad)
        f.write('lgb_filt_energy_frac = %s \n' % inp_ffproc.lgb_filt_energy_frac)
        f.write('lgb_filt_nlambda = %s \n' % inp_ffproc.lgb_filt_nlambda)
        f.write('selected window lengths: %s \n' % inp_ffproc.lgb_filt_duration)
        f.write('dominant periods: %s \n' % inp_ffproc.lgb_filt_center_period)

        f.write('\n--------------------------------------------------------------\n')
        f.write('STF information:\n')
        f.write('--------------------------------------------------------------\n\n')
        f.write('address_stf: %s\n' % all_stfs.address_stf)
        f.write('qual_stf: %s' % inp_ffproc.stf_min_qual)

        f.write('\n--------------------------------------------------------------\n')
        f.write('Paths defined by you:\n')
        f.write('--------------------------------------------------------------\n\n')
        f.write('real_path = %s \n' % inp_ffproc.real_path)
        f.write('syn_path = %s \n' % inp_ffproc.syn_path)
        f.write('evproc_path = %s \n' % inp_ffproc.evproc_path)
        f.write('output_dir = %s \n' % inp_ffproc.output_dir)

        f.close()


# ------------------- create_stxml_files --------------------

def create_stxml_files(tr, corr_sei, corr_hyd, real_mode):

    '''

    instr_corr_sps = 250
    instr_corr_instr = ['Trillium', 'HTI', 'HELGA']  # for 0 - seismometer, 1 - hydrophone, 2 - hydrophone
    save_files = True    # save instrument corrected files and resp files
    path_to_icf = '/Volumes/myDODO/Test_OBSTB_data/TestSave'   # where to save instrument corrected files and resp files
    '''

    # find corresponding base station correction file
    # import ipdb; ipdb.set_trace()
    # structure of '/output/STATION/YYYY-MM-DD/processed' or '/output/STATION/YYYY-MM-DD/resp'

    if tr.stats.channel in ['Z', '1', '2', 'E', 'N', 'X', 'Y']:
        inv = read_inventory(corr_sei)
        inv[0]._code = tr.stats.network
        inv[0][0]._code = tr.stats.station
        inv[0][0][0]._code = tr.stats.channel
        inv[0][0][0].location_code = tr.stats.location
        inv[0][0][0].start_date = tr.stats.starttime
        inv[0][0][0].end_date = tr.stats.endtime

        return inv

    elif tr.stats.channel == 'H':
        inv = read_inventory(corr_hyd)
        inv[0]._code = tr.stats.network
        inv[0][0]._code = tr.stats.station
        inv[0][0][0]._code = tr.stats.channel
        inv[0][0][0].location_code = tr.stats.location
        inv[0][0][0].start_date = tr.stats.starttime
        inv[0][0][0].end_date = tr.stats.endtime

        return inv

    elif len(tr.stats.channel) == 3:
        
        # XXX need to improve this list? xXX
        if tr.stats.channel in ['HHZ', 'HHX', 'HHY', 
                                'BHZ', 'BH1', 'BH2', 
                                'BHE', 'BHN', 'BHX', 'BHY']:
            # import ipdb; ipdb.set_trace()
            if real_mode == 'repo':
                # RESP.IA.IA014..HHN
                inv = read_inventory(os.path.join(corr_sei, 'RESP.%s.%s..%s' % (tr.stats.network, tr.stats.station, tr.stats.channel)))
            else:
                inv = read_inventory(corr_sei)
        elif tr.stats.channel == 'BDH':
            inv = read_inventory(corr_hyd)
        else:
            sys.exit('ERROR(2)!')
        inv[0]._code = tr.stats.network
        inv[0][0]._code = tr.stats.station
        inv[0][0][0]._code = tr.stats.channel
        inv[0][0][0].location_code = tr.stats.location
        inv[0][0][0].start_date = tr.stats.starttime
        inv[0][0][0].end_date = tr.stats.endtime

        return inv

    else:
        # import ipdb; ipdb.set_trace()
        sys.exit('[SUBFUNCTION.PY] [create_stxml_files] This channel is not implemented yet. EXIT!')

