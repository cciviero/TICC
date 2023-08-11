#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  statioin_plotter.py
#   Purpose:   simple plotting tool to load and plot waveforms of desired
#              event, station, channels for comparison. Saves plot.
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

import matplotlib.pylab as plt
import glob
import os
import pickle
import numpy as np

from obspy import read
from obspy.taup import TauPyModel
from obspy.geodetics import gps2DistAzimuth
from obspy.geodetics import locations2degrees

# -----------------------------------------------------------------------
# ------------------------------- INPUT ---------------------------------
# -----------------------------------------------------------------------

path_2_dataset = '/Users/maria/PhD/Codes/__TEST_DATA/The17'
event = '20130416_104420.a'
# event = '20130722_070142.a'
station = 'RR49'

filter_min = 30
filter_max = 20

ph_bg = 'iasp91'
phase = 'P'


# -----------------------------------------------------------------------
# ---------------------------- Functions --------------------------------
# -----------------------------------------------------------------------

def read_event(event_path):
    '''

    :param event_path:
    :return:
    '''

    ev_load = open(os.path.join(event_path, 'info', 'event.pkl'), 'r')
    ev_pkl = pickle.load(ev_load)
    ev_load.close()
    # import ipdb; ipdb.set_trace()

    ev_lat = ev_pkl['latitude']
    ev_lon = ev_pkl['longitude']
    ev_dep = ev_pkl['depth']
    ev_origt = ev_pkl['datetime']

    return ev_lat, ev_lon, ev_dep, ev_origt


def calculate_arrival_time(ev_lat, ev_lon, ev_dep, sta_lat, sta_lon, stel, bg_model, req_phase):
    # import ipdb; ipdb.set_trace()

    (dist, azi, bazi) = gps2DistAzimuth(ev_lat, ev_lon, sta_lat, sta_lon)

    dist_deg = locations2degrees(ev_lat,
                                 ev_lon,
                                 sta_lat,
                                 sta_lon)
    print dist_deg, dist, azi, bazi

    model = TauPyModel(model=bg_model)
    try:
        print req_phase, ev_dep, stel
        # import ipdb; ipdb.set_trace()
        # arrivals = model.get_travel_times(source_depth_in_km=ev_dep,
        #                        distance_in_degree=azi,
        #                        phase_list=[req_phase],
        #                        receiver_depth_in_km=stel)
        # print arrivals

        arrivals = model.get_travel_times(source_depth_in_km=ev_dep,
                                    distance_in_degree = dist_deg,
                                    phase_list = [req_phase])
        ttime = arrivals[0].time
    except Exception, exp:
        print exp
        import ipdb;
        ipdb.set_trace()
        tt = False

    return ttime

def read_inv_dmt(path_2_dataset, event, tr):
    """
    read station_event file from DMT directories
    :param inp_ffproc:
    :param all_events:
    :param ev:
    :param tr:
    :return:
    """
    inv_info = np.loadtxt(os.path.join(path_2_dataset, event,
                                       'info',
                                       'station_event'),
                          dtype='object', delimiter=',')
    inv_info_net = inv_info[inv_info[:, 0] == tr.stats.network]
    inv_info_sta = inv_info_net[inv_info_net[:, 1] == tr.stats.station]
    inv_info_loc = inv_info_sta[inv_info_sta[:, 2] == tr.stats.location]
    inv_info_cha = inv_info_sta[np.where(inv_info_loc[:, 3] == tr.stats.channel)]
    if len(inv_info_sta) == 0:
        return False, False, False, False
    # import ipdb; ipdb.set_trace()
    stla = inv_info_cha[0][4]
    stlo = inv_info_cha[0][5]
    stel = inv_info_cha[0][6]
    stdp = inv_info_cha[0][7]

    return stla, stlo, stdp, stel
# -----------------------------------------------------------------------
# -------------------------------- MAIN ---------------------------------
# -----------------------------------------------------------------------

event_path = os.path.join(path_2_dataset, event)

station_path = glob.glob(os.path.join(event_path, 'processed', '*' + station + '*'))

ev_lat, ev_lon, ev_dep, ev_origt = read_event(event_path)
tr = read(station_path[0])[0]
stla, stlo, stdp, stel = read_inv_dmt(path_2_dataset, event, tr)
print stla, stlo, stdp, stel
print eval(stla), eval(stlo), eval(stdp), eval(stel)
tt = calculate_arrival_time(ev_lat, ev_lon, ev_dep, eval(stla), eval(stlo), eval(stel), ph_bg, phase)
print tt
if not tt == False:
    # import ipdb;
    # ipdb.set_trace()
    vline = ev_origt + tt
    diff = vline - tr.stats.starttime

plt.figure(figsize=(15,10))
plt.ion()
for i in range(len(station_path)):

    tr = read(station_path[i])[0]
    tr.filter('bandpass', freqmin=1./filter_min, freqmax=1./filter_max)

    cha = tr.stats.channel
    x_axis = tr.times()
    y_axis = tr.data
    plt.subplot(len(station_path), 1, i+1)
    plt.plot(x_axis, y_axis)
    plt.axvline(diff, 0, 1, c='r', ls='--', lw=2)

    plt.title(cha, weight='bold', size=12)

plt.suptitle('%s - %s - Filter: %s - %s Hz' % (event, station, filter_min, filter_max))
plt.show()

raw_input('press enter')