#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  create_station_list.py
#   Purpose:   from obspyDMT created database extract unique station list
#   Author:    Maria Tsekhmistrenko
#   Email:     mariat@earth.ox.ac.uk
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------


import glob
import os
import time
import numpy as np

# -------------------------------------------------------------------
# --------------------------- INPUT ---------------------------------
# -------------------------------------------------------------------

data_path = '/mnt/seismodata/MT/ICELAND-TOMO/ice_data'
folder_type = '*_*.*'

# -------------------------------------------------------------------
# --------------------------- INPUT ---------------------------------
# -------------------------------------------------------------------

start_time = time.time()

print '\n====>>>> Reading events info\n'
all_events_path = glob.glob(os.path.join(data_path, folder_type))
all_events_path.sort()


# e.g.: YV, MAID, 00, HHE, -21.0797, 55.3831, 2169.0, 0.0, RESIF, 20121111_011238.a,23.005,95.885,13.7,6.8,110,
dt = np.dtype([('network',  'S10'), ('station', 'S10'), ('location', 'S10'), ('channel', 'S10'), ('stla', float),
               ('stlo', float), ('stel', float), ('stdp', float), ('source', 'S10'), ('ev_id', 'S20'),
               ('evla', float), ('evlo', float), ('evdp', float), ('mag', float)])

station_list_comp = np.array([], dtype=dt)
calc_length = 0
for event in all_events_path:
    # print event
    try:
        sta_ev = np.loadtxt(os.path.join(event, 'info', 'station_event'), dtype=dt, delimiter=',',
                            usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13))
        # only get the vertical component for now...
        filter1 = sta_ev['channel'] == 'BHZ'
        filter2 = sta_ev['channel'] == 'HHZ'
        update = sta_ev[filter1 + filter2]

        # print 'Length is: %s ' % len(update)

        station_list_comp = np.append(station_list_comp, update)

        calc_length += len(update)

    except Exception, e:

        print event, ':', e

uniq_sta_list, uniq_sta_index = np.unique(station_list_comp['station'], return_index=True)

print '\nUnique stations:', len(uniq_sta_list)

outfile_fio = open(os.path.join(data_path, 'station_file.txt'), 'w')
outfile_fio.writelines('# stla\tstlo\tnet.stat.loc.chan\tstel\tstdp\n')

for index in uniq_sta_index:
    outfile_fio.writelines('%20s %20s %5s.%s.%s.%s %15s %10s\n'
                           % (station_list_comp['stla'][index],
                              station_list_comp['stlo'][index],
                              station_list_comp['network'][index],
                              station_list_comp['station'][index],
                              station_list_comp['location'][index],
                              station_list_comp['channel'][index],
                              station_list_comp['stel'][index],
                              station_list_comp['stdp'][index]))


outfile_fio.close()
print('\n====>>>> DONE. Time lapsed: %s seconds' % (round(time.time() - start_time, 2)))




