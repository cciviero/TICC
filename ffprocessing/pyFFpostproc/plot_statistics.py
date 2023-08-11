# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  pyffpp_src.py
#   Purpose:   pyffpp: Python Finite Frequency Post Proc
#   Author:    Maria Tsekhmistrenko
#   Email:     mariat@earth.ox.ac.uk
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")

import ConfigParser
import glob
import copy
import sys
import os

import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from obspy.taup import TauPyModel
from obspy.geodetics import gps2dist_azimuth

import time
import pickle

from collections import Counter
from operator import itemgetter

# -----------------------------------------------------------------------
# ------------------------------ INPUT ----------------------------------
# -----------------------------------------------------------------------
input_path = '/mnt/seismodata2/MT/SEA-SEIS/ffproc/SEA-SEIS_18.01.2021'
# path to event_list_pickle in EVENTS-INFO folder
original_path = '/mnt/home_geo/mariat/Codes/_INPUT_/neic_2018_2020'
# original_path = '/mnt/seismodata/MT/ICELAND'

cc_factor = 0.70
ts_factor = 20
clip_factor = 0
snr_factor = 0

plot_cc_hist = True
plot_ts_hist = True
plot_piechart = False
plot_snr_hist = True
plot_net_hist = False
plot_map_hc = True

# -----------------------------------------------------------------------
# --------------------------- FUNCTIONS ---------------------------------
# -----------------------------------------------------------------------


# ------------------ read output files genereted from pyffproc  ------------------
def read_ampstt(path, event_name, band):
    """
    In order to extract quickly the information from ffproc.ampstt.band0?
    generated by pyffporc
    :param path: path to the folder where the events are stored going
                 to the ampstt files
    :param event_name: the name of the event, e.g. '20130924_35'
    :param band: which band should be read
    :return: the whole file as a np.array[line -> stations / column -> information]
    How to read the station list:
    0 idx
    1 grp
    2 stla
    3 stlo
    4 stazie
    5 stnam
    6 xc_coeff
    7 Tobs
    8 dT
    9 sigma_dT
    10 not_used
    11 not_used
    12 A
    13 sigma_A
    14 not_used
    15 not_used
    16 tB_smgr
    17 tB_mfi
    18 winlen
    19 clip_taumax
    not_used
    SNR
    [ts_step1    ts_step2    cc_step1    cc_step2    clip_step1    clip_step2]
    """

    dt = np.dtype([('idx', int), ('grp', int), ('stlat', float), ('stlon', float), ('epd', float),
                   ('name', 'S10'), ('cc', float), ('tobs', float), ('dt', float), ('sigma_t', float),
                   ('N/A_1', 'S10'), ('N/A_2', 'S10'), ('f_amp', float), ('sigma_a', float), ('N/A_3', 'S10'),
                   ('N/A_4', 'S10'), ('tB_smgr', float), ('tB_mfi', float),
                   ('winlen', float), ('clip', int), ('N/A_5', 'S10'), ('snr', float)])

    ampstt_dir = os.path.join(path, event_name, 'outfiles', 'ffproc.ampstt.band%02d' % band)
    try:
        station_list = np.loadtxt(ampstt_dir, dtype=dt, ndmin=1)
    except Exception, e:
        print '%s:This is not an ffproc.ampstt file' % ampstt_dir
        return None

    if np.shape(station_list) == ():
        station_list = np.array([station_list])

    return station_list


#  ------------------ read kernel_filter.txt ------------------
def read_kernel_filter(path, event_name):
    '''
    :param path:
    :param event_name:
    :return:
    '''
    dt = np.dtype([('band', 'S10'), ('type', 'S10'), ('dominant_period', float), ('sigma', float)])
    # import ipdb; ipdb.set_trace()
    kernel_filter_dir = os.path.join(path, event_name, 'outfiles', 'kernel_filter.txt')
    band_info = np.loadtxt(kernel_filter_dir, dtype=dt, skiprows=1, usecols=(0, 1, 2, 3))

    return band_info

#  ------------------ read event pkl ------------------
def read_event_pkl(path):
    # import ipdb; ipdb.set_trace()

    ev_load = open(os.path.join(path, 'EVENTS-INFO', 'event_list_pickle'), 'r')
    ev_pkl = pickle.load(ev_load)
    ev_load.close()

    return ev_pkl


# -----------------------------------------------------------------------
# -------------------------------- MAIN ---------------------------------
# -----------------------------------------------------------------------

path_events = glob.glob(os.path.join(input_path, '*.?'))
all_events = []
for i in path_events:
    join = i.split('/')[-1]
    all_events.append(join)
#import ipdb; ipdb.set_trace()
inp_band = read_kernel_filter(input_path, all_events[0])
bands_sec = inp_band['dominant_period']

# create output directory according to the input parameters
output_dir = os.path.join(input_path, 'stat_plots_%s_%s_%s/'
                          % (cc_factor, ts_factor, clip_factor))
if not os.path.isdir(output_dir):
    os.makedirs(os.path.join(output_dir))

# this list will contain the information from all bands (various lengths)
# therefore different in length
final_cc = []
final_ts = []
final_amp = []
final_mag = []
final_azi = []
final_name = []
final_stlo = []
final_stla = []
final_snr = []

check_all_events = copy.copy(all_events)

# this counter counts the overall amount for all events, bands and
# stations that are fitting the criteria
count_criteria_true = 0
count_criteria_false = 0

for i, value in enumerate(bands_sec):
    band = i + 1

    cc_array = np.array([])
    t_shift_array = np.array([])
    sta_name_array = np.array([])
    ampl_array = np.array([])
    mag_array = np.array([])
    azi_array = np.array([])

    snr_array = np.array([])
    stla_array = np.array([])
    stlo_array = np.array([])

    remove_bad_events = copy.copy(all_events)

    # this counter counts only the band but all events and stations that
    # fit the criteria
    count_stations = 0

    string = 'Working on Band %s' % band
    print '%s' % (len(string) * '=')
    print '%s' % string
    print '%s' % (len(string) * '=')

    for event in check_all_events:
        # import ipdb; ipdb.set_trace()
        info_pkl = read_event_pkl(original_path)
        indx = [i for i, ev in enumerate(info_pkl) if ev['event_id'] == event][0]
        mag = info_pkl[indx]['magnitude']

        print '%s' % event
        station_list = read_ampstt(input_path, event, band)
        # import ipdb; ipdb.set_trace()

        if station_list is None:
            remove_bad_events.remove(event)
            print '%s will be removed from list!' % event
            continue

        original_length = len(station_list)

        filter_crit = station_list['cc'] >= cc_factor
        filter_crit *= station_list['clip'] == clip_factor
        filter_crit *= abs(station_list['dt']) <= ts_factor

        station_list = station_list[filter_crit]

        updated_length = len(station_list)
        sorted_out = original_length - updated_length

        count_criteria_true += updated_length
        count_criteria_false += sorted_out
        count_stations += original_length

        cc_array = np.append(cc_array, station_list['cc'])
        t_shift_array = np.append(t_shift_array, station_list['dt'])
        sta_name_array = np.append(sta_name_array, station_list['name'])
        ampl_array = np.append(ampl_array, station_list['f_amp'])
        azi_array = np.append(azi_array, station_list['epd'])
        mag_array = np.append(mag_array, [mag]*updated_length)
        snr_array = np.append(snr_array, station_list['snr'])
        stla_array = np.append(stla_array, station_list['stlat'])
        stlo_array = np.append(stlo_array, station_list['stlon'])

    final_cc.append(cc_array)
    final_ts.append(t_shift_array)
    final_amp.append(ampl_array)
    final_mag.append(mag_array)
    final_azi.append(azi_array)
    final_name.append(sta_name_array)
    final_snr.append(snr_array)
    final_stla.append(stla_array)
    final_stlo.append(stlo_array)

print '\n'
# I don't know why i defined this?
bands = bands_sec
# -----------------------------------------------------------------------
# plot CC histogram
if plot_cc_hist:
    print 'plotting cc histogram'

    label = 'cross correlation [%]'
    save_label = 'cc'
    plt.ioff()
    fig = plt.figure(frameon=False)
    fig.set_size_inches(30, 20)

    # for the specific band different colors and linewidths
    cmap = plt.get_cmap('jet_r')
    # import ipdb; ipdb.set_trace()
    colors = cmap(np.linspace(0, 1.0, len(bands)))
    line_width = range(len(bands) + 1)[::-1][0:-1]
    band_freq = bands

    bins_bin = len(np.arange(cc_factor, 1.0, 0.03)) - 1
    bins_range = [cc_factor, 1.0]
    plt.xlim(cc_factor - 0.05, 1.0)
    plt.axvline(x=cc_factor, color='r', lw=3, linestyle='--')

    # import ipdb; ipdb.set_trace()
    for idx, band in enumerate(bands):
        # import ipdb; ipdb.set_trace()
        y, binEdges = np.histogram(final_cc[idx].astype(np.float), bins=bins_bin, range=bins_range)
        bin_centers = 0.5 * (binEdges[1:] + binEdges[:-1])

        # import ipdb; ipdb.set_trace()
        plt.plot(bin_centers, y, color=colors[idx], lw=line_width[idx],
                 label='B0%s - %2.2fs: %s data points' % (band, band,
                                                          len(final_cc[idx])))

    plt.xlabel('%s' % label, size=30, weight='bold')
    plt.ylabel('number of measurements', size=30, weight='bold')
    plt.title('Histogram of Measurements\nCriteria: cc>=%s -- ts <=%s -- clip=%s'
              % (cc_factor, ts_factor, clip_factor), size=40, weight='bold')
    plt.legend(fontsize=25, loc=3)

    plt.grid(True)

    # to put the text relative to the x/y axis
    plt.annotate('Number of all measurements pro band: %s\n'
                 'Number of good measurements fitting the criteria (all bands): %s\n'
                 'Number of rejected measurements (not fitting criteria): %s'
                 % (count_stations, count_criteria_true, count_criteria_false),
                 xy=(0.02, 0.90), size=23, xycoords='axes fraction',
                 bbox=dict(boxstyle="round", facecolor='red', alpha=0.8))

    plt.xticks(weight='bold', size=18)
    plt.yticks(weight='bold', size=18)

    plt.savefig(output_dir + '%s_%s_%s_%s.png' % (save_label, cc_factor, ts_factor, clip_factor), dpi=500)

# -----------------------------------------------------------------------
# plot time shift histogram
if plot_ts_hist:
    print 'plotting time shift histogram'

    label = 'time shift [s]'
    save_label = 'time_shift'

    plt.ioff()
    fig = plt.figure(frameon=False)
    fig.set_size_inches(30, 20)

    # for the specific band different colors and linewidths
    cmap = plt.get_cmap('jet_r')
    # import ipdb; ipdb.set_trace()
    colors = cmap(np.linspace(0, 1.0, len(bands)))
    line_width = range(len(bands) + 1)[::-1][0:-1]

    bins_bin = len(np.arange(-ts_factor, ts_factor, 1))
    bins_range = [-ts_factor, ts_factor]
    plt.axvline(x=0, color='r', lw=3, linestyle='--')
    plt.axvline(x=0.5, color='r', lw=2, linestyle='--')
    plt.axvline(x=-0.5, color='r', lw=2, linestyle='--')
    plt.axvline(x=1, color='r', lw=1, linestyle='--')
    plt.axvline(x=-1, color='r', lw=1, linestyle='--')
    plt.axvline(x=-2, color='orange', lw=1, linestyle='--')
    plt.axvline(x=2, color='orange', lw=1, linestyle='--')

    # import ipdb; ipdb.set_trace()
    for idx, band in enumerate(bands):
        # import ipdb; ipdb.set_trace()
        y, binEdges = np.histogram(final_ts[idx].astype(np.float), bins=bins_bin, range=bins_range)
        bin_centers = 0.5 * (binEdges[1:] + binEdges[:-1])

        # import ipdb; ipdb.set_trace()
        plt.plot(bin_centers, y, color=colors[idx], lw=line_width[idx],
                 label='B0%s - %2.2fs: %s data points' % (band, band,
                                                          len(final_ts[idx])))

    plt.xlabel('%s' % label, size=30, weight='bold')
    plt.ylabel('number of measurements', size=30, weight='bold')
    plt.title('Histogram of Measurements\nCriteria: cc>=%s -- ts <=%s -- clip=%s'
              % (cc_factor, ts_factor, clip_factor), size=40, weight='bold')
    plt.legend(fontsize=25, loc=3)

    plt.grid(True)

    # to put the text relative to the x/y axis
    plt.annotate('Number of all measurements pro band: %s\n'
                 'Number of good measurements fitting the criteria (all bands): %s\n'
                 'Number of rejected measurements (not fitting criteria): %s'
                 % (count_stations, count_criteria_true, count_criteria_false),
                 xy=(0.02, 0.90), size=23, xycoords='axes fraction',
                 bbox=dict(boxstyle="round", facecolor='red', alpha=0.8))

    plt.xticks(weight='bold', size=18)
    plt.yticks(weight='bold', size=18)

    plt.savefig(output_dir + '%s_%s_%s_%s.png' % (save_label, cc_factor, ts_factor, clip_factor), dpi=500)

# -----------------------------------------------------------------------
# plot pie chart
if plot_piechart:
    print 'plotting pie chart'

    net_array = ([])

    plt.ioff()
    fig = plt.figure(frameon=False)
    fig.set_size_inches(30, 20)

    for band in final_name:
        unique_stations = set(band)
        len_uniq_stat = len(unique_stations)
        for station in unique_stations:
            net = station.split('.')[0]
            net_array = np.append(net_array, net)

    count = Counter(net_array)
    # list_of_networks = sorted(count)
    total_number = sum(count.values())

    sort_count = np.array(sorted(count.items(), key=itemgetter(1)))

    labels = sort_count[:, 0]
    sizes = sort_count[:, 1].astype(np.float) * 100 / total_number

    cmap = plt.get_cmap('jet')
    colors = cmap(np.linspace(0, 1.0, len(labels)))
    plt.rcParams['font.size'] = 9.0

    explode = []
    for i in range(len(labels)):
        if labels[i] == 'YV' or labels[i] == 'XV' or labels[i] == 'ZE' or labels[i] == 'X1' \
                or labels[i] == '3E':
            explode.append(0.2)
        else:
            explode.append(0)

    #plt.pie(sizes, explode=explode, labels=labels, colors=colors,
    #        autopct='%1.1f%%', shadow=False, startangle=90)
    plt.pie(sizes, explode=explode, colors=colors,
            shadow=False, startangle=90)      
    # Set aspect ratio to be equal so that pie is drawn as a circle.
    plt.axis('equal')
    plt.savefig(os.path.join(output_dir, 'pie_chart_network_contribution.png'), dpi=500)

# -----------------------------------------------------------------------
# plot network histogram
if plot_net_hist:

    label = 'Signal to Noise ratio histogram'
    save_label = 'snr_hist'

    net_array = np.array([])

    plt.ioff()
    fig = plt.figure(frameon=False)
    fig.set_size_inches(30, 20)

    # get the unique network names
    for band in final_name:
        unique_stations = set(band)
        len_uniq_stat = len(unique_stations)
        for station in unique_stations:
            net = station.split('.')[0]
            net_array = np.append(net_array, net)

    count = Counter(net_array)
    # list_of_networks = sorted(count)
    sort_count = np.array(sorted(count.items(), key=itemgetter(1)))
    labels = sort_count[:, 0]

    x = []
    y = []

    count = 0

    cmap = plt.get_cmap('jet_r')
    colors = cmap(np.linspace(0, 1.0, len(bands)))
    plt.rcParams['font.size'] = 9.0

    shift = np.arange(-0.25, 0.25, 0.07)

    # get for each band/net the snr

    for net in labels:
        count += 1
        y = []
        for i, name_band in enumerate(final_name):
            temp = []
            for j, name in enumerate(name_band):
                if net in name:
                    temp.append(final_snr[i][j])
            y.append(np.mean(temp))
        plt.bar(count + shift, y[::-1], width=0.25, color=colors[::-1])
        x.append(count)

    plt.xticks(x, labels)
    plt.xticks(weight='bold', size=18)
    plt.yticks(weight='bold', size=18)

    plt.savefig(output_dir + '%s_%s_%s_%s.png' % (save_label, cc_factor, ts_factor, clip_factor), dpi=500)

# -----------------------------------------------------------------------
# plot cc histogram
if plot_cc_hist:

    label = 'CC network wise'
    save_label = 'cc_hist'

    net_array = np.array([])

    plt.ioff()
    fig = plt.figure(frameon=False)
    fig.set_size_inches(30, 20)

    # get the unique network names
    for band in final_name:
        unique_stations = set(band)
        len_uniq_stat = len(unique_stations)
        for station in unique_stations:
            net = station.split('.')[0]
            net_array = np.append(net_array, net)

    count = Counter(net_array)
    # list_of_networks = sorted(count)
    sort_count = np.array(sorted(count.items(), key=itemgetter(1)))
    labels = sort_count[:, 0]

    x = []
    y = []

    count = 0

    cmap = plt.get_cmap('jet_r')
    colors = cmap(np.linspace(0, 1.0, len(bands)))
    plt.rcParams['font.size'] = 9.0

    shift = np.arange(-0.25, 0.25, 0.07)

    # get for each band/net the snr

    for net in labels:
        count += 1
        y = []
        for i, name_band in enumerate(final_name):
            temp = []
            for j, name in enumerate(name_band):
                if net in name:
                    temp.append(final_cc[i][j])
            y.append(np.mean(temp))
        plt.bar(count + shift, y[::-1], width=0.25, color=colors[::-1])
        x.append(count)

    plt.ylim(0.7, 1.01)
    plt.xticks(x, labels)
    plt.xticks(weight='bold', size=18)
    plt.yticks(weight='bold', size=18)

    plt.savefig(output_dir + '%s_%s_%s_%s.png' % (save_label, cc_factor, ts_factor, clip_factor), dpi=500)


# -----------------------------------------------------------------------
# plot SNR histogram
if plot_snr_hist:
    print 'plot snr histogram'

    label = 'Signal to Noise ratio'
    save_label = 'snr'

    plt.ioff()
    fig = plt.figure(frameon=False)
    fig.set_size_inches(30, 20)

    # for the specific band different colors and linewidths
    cmap = plt.get_cmap('jet_r')
    # import ipdb; ipdb.set_trace()
    colors = cmap(np.linspace(0, 1.0, len(bands)))
    line_width = range(len(bands) + 1)[::-1][0:-1]

    max_snr = max([max(i) for i in final_snr])
    bins_bin = len(np.arange(-1, max_snr, 10)) - 1
    if bins_bin <= 0:
        bins_bin = 1
    bins_range = [-1, max_snr]
    plt.xlim(-1, max_snr)
    plt.axvline(x=1, color='r', lw=3, linestyle='--')

    # import ipdb; ipdb.set_trace()
    for idx, band in enumerate(bands):

        try:
        # import ipdb; ipdb.set_trace()
            y, binEdges = np.histogram(final_snr[idx].astype(np.float), bins=bins_bin, range=bins_range)
            bin_centers = 0.5 * (binEdges[1:] + binEdges[:-1])

            # import ipdb; ipdb.set_trace()
            plt.plot(bin_centers, y, color=colors[idx], lw=line_width[idx])
        except Exception, exp:
            import ipdb; ipdb.set_trace()
            print exp
            pass


    plt.xlabel('%s' % label, size=30, weight='bold')
    plt.ylabel('number of measurements', size=30, weight='bold')
    plt.title('Histogram of Measurements\nCriteria: cc>=%s -- ts <=%s -- clip=%s'
              % (cc_factor, ts_factor, clip_factor), size=40, weight='bold')
    plt.grid(True)

    plt.xticks(weight='bold', size=18)
    plt.yticks(weight='bold', size=18)

    plt.savefig(output_dir + '%s_%s_%s_%s.png' % (save_label, cc_factor, ts_factor, clip_factor), dpi=500)

# -----------------------------------------------------------------------
# plot map hitcount
if plot_map_hc:
    print 'plot map hitcount for each station'

    label = 'Hit Count'
    save_label = 'map_hc'
    sta_array = np.array([])
    coord_lat = np.array([])
    coord_lon = np.array([])

    # count how many times one stations comes up
    for i, band in enumerate(final_name):
        for j, station in enumerate(band):
            sta_array = np.append(sta_array, station)


    count = Counter(sta_array)

    sort_count = np.array(sorted(count.items(), key=itemgetter(1)))

    labels = sort_count[:, 0]
    numbers = sort_count[:, 1].astype(float)
    max_count = np.mean(numbers)
    min_count = numbers[0]

    for sta in labels:
        count = 0
        for i, name_band in enumerate(final_name):
            for j, name in enumerate(name_band):
                if sta in name and count < 1:
                    coord_lat = np.append(coord_lat, final_stla[i][j])
                    coord_lon = np.append(coord_lon, final_stlo[i][j])
                    count += 1
                else:
                    pass

    #import ipdb; ipdb.set_trace()

    plt.ioff()
    fig = plt.figure(frameon=False)
    fig.set_size_inches(30, 20)


    mymap = Basemap(projection='robin', lon_0=0, lat_0=0)

    x, y = mymap(coord_lon, coord_lat)

    mymap.scatter(x, y, c=numbers, cmap=plt.cm.get_cmap('jet', 256), s=150, zorder=10,
                  edgecolors='none',
                  vmin=min_count, vmax=max_count)

    mymap.drawcoastlines(color='black', linewidth=3)

    cbar = plt.colorbar(orientation='horizontal', shrink=0.5)
    cbar.ax.tick_params(labelsize=10)
    plt.savefig(output_dir + '%s_%s_%s_%s.png' % (save_label, cc_factor, ts_factor, clip_factor), dpi=500)