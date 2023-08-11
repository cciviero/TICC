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

# import ConfigParser
import glob
import copy
import sys
import os

import numpy as np
import pandas as pd

import matplotlib as mpl
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
input_path = '/mnt/seismodata2/MT/SEA-SEIS/ffproc/SEA-SEIS_rotation_testing_v3'
#input_path = '/mnt/seismodata2/MT/_OUTPUT_/TICC_ffproc_testing_rotation'
# path to event_list_pickle in EVENTS-INFO folder
# original_path = '/Users/mary/Codes/_INPUT_/neic_2018_2020'

cc_factor = 0.70

output_path = '/mnt/seismodata2/MT/SEA-SEIS/ffproc/SEA-SEIS_rotation_testing_v3/rotation'
#output_path = '/mnt/seismodata2/MT/_OUTPUT_/TICC_ffproc_testing_rotation/rotation'
if not os.path.isdir(output_path):
    os.makedirs(os.path.join(output_path))


# -----------------------------------------------------------------------
# --------------------------- FUNCTIONS ---------------------------------
# -----------------------------------------------------------------------


# ------------------ read output files genereted from pyffproc  ------------------
def read_pm(path, station, event_name):
    """
    In order to extract quickly the information from ffproc.ampstt.band0?
    generated by pyffporc
    :param path: path to the folder where the events are stored going
                 to the ampstt files
    :param event_name: the name of the event, e.g. '20130924_35'

    :return: the whole file as a np.array[line -> stations / column -> information]
    How to read the station list:
    0 - station name
    1 - band
    2 - stalat
    3 - stalon
    4 - evlat
    5 - evlon
    6 - xc_coeff
    7 - AZI catalouge
    8 - BAZI catalouge
    9 - GCD
    10 - AZI AB
    11 - AZI BA
    12 - azi_pm
    13 - inci_pm
    14 - err_azi_pm
    15 - err_inci_pm
    16 - azi_flinn
    17 - inci_flinn
    18 - rectilinearity_flinn
    19 - planarity_flinn
    20 - orientation
    """

    ampstt_dir = os.path.join(path, event_name, 'outfiles', 'ffproc.pm.%s' % station)
    try:
        station_list = pd.read_csv(ampstt_dir, delimiter="\t", comment="#", header=None).to_numpy()
    except Exception, e:
        # import ipdb; ipdb.set_trace()
        print '%s:This is not an ffproc.ampstt file' % ampstt_dir
        return None
    if np.shape(station_list) == ():
        station_list = np.array([station_list])
    return station_list

# -----------------------------------------------------------------------
# -------------------------------- MAIN ---------------------------------
# -----------------------------------------------------------------------

# find all unique stations

all_events = glob.glob(os.path.join(input_path, '*.a'))
#import ipdb; ipdb.set_trace()
all_stations = np.array([])
for event in all_events:
    station_file = os.path.join(event, 'outfiles/ffproc.receivers')
    station_list = np.loadtxt(station_file, dtype=object, usecols=2, comments='#')
    all_stations = np.append(all_stations, [i.split('.')[1] for i in station_list])

uniq_stations = np.unique(all_stations)

for station in uniq_stations:
    collector = np.array([])
    events_name = np.array([])
    events_counter = np.array([])
    counter = 0
    for event in all_events:

        station_file = os.path.join(event, 'outfiles/ffproc.pm.%s' % station)
        try:
            pm_info = read_pm(station_file, station, event)
            if counter == 0:
                collector = pm_info
            else:
                collector = np.r_[collector, pm_info]
            events_name = np.append(events_name, [os.path.basename(event).split('.')[0]] * 8)
            events_counter = np.append(events_counter, [counter] * 8)
            counter += 1
        except Exception, exp:
            print exp, 'skipping this one..'
            # import ipdb; ipdb.set_trace()
            continue

    # create one figure with all event measurements
    # import ipdb; ipdb.set_trace()
    try:
        if cc_factor > 0.0:
            filter_list = collector[:, 6] >= cc_factor
            # not good here XXX change as soon as possible XXX
            collector = collector[filter_list]
            events_name = events_name[filter_list]
            events_counter = events_counter[filter_list]
    except Exception, exp:
        continue
    # ==== SAVE collector to file ====
    # import ipdb; ipdb.set_trace()
    array_to_export = np.c_[events_name, collector]
    np.savetxt(os.path.join(output_path, 'pm_information_%s.txt' % station), array_to_export, delimiter="\t\t", fmt="%s")


    # ==== FIGURE ====
    plt.figure(figsize=(20, 10))

    # Color handling
    cmap = plt.cm.jet
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # create the new map
    cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
    # define the bins and normalize
    bounds = np.linspace(0, 1, 11)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    plt.subplot(3, 1, 1)
    plt.title('BAZ measured | mean: %s | std: %s' % (np.mean(collector[:, 12]), np.std(collector[:, 12])))
    plt.scatter(events_counter, collector[:, 12], c=collector[:, 6], s=(collector[:, 6]*1000).astype(int),
                cmap=cmap, norm=norm, label='calculated BAZ')
    plt.ylim(0, 180)
    plt.xlabel('events')
    plt.xlim(0, int(events_counter[-1]))

    plt.subplot(3, 1,  2)
    plt.title('Orientation | mean: %s | std: %s' % (np.mean(collector[:, -1]), np.std(collector[:, -1])))
    plt.scatter(events_counter, collector[:, -1], c=collector[:, 6],
                cmap=cmap, norm=norm, s=(collector[:, 6]*1000).astype(int))
    plt.xlim(0, int(events_counter[-1]))
    plt.xlabel('events')

    plt.subplot(3, 1, 3)
    plt.title('BAZ expected')
    plt.scatter(events_counter, collector[:, 10], s=150, marker='x', c='red', alpha=0.3, label='BAZ')
    plt.xlim(0, int(events_counter[-1]))

    locs, labels = plt.xticks()
    plt.xticks(np.unique(events_counter), np.unique(events_name), rotation=90)  # Set locations and labels
   
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, '%s_orientation.png' % station), dpi=300)

print 'Done.'