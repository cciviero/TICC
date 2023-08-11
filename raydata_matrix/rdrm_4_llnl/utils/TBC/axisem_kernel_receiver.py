#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  axisem_kernel_receiver.py
#   Purpose:   Generate input files (receiver.dat) for AXISEM
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import datetime
import numpy as np

import output_reader as outread

# ------------------- INPUT -----------------------------
#master_dir = '/home/hosseini/Work/Scripts/gitHUB/MEASUREMENTS/P_measure_1_sec_LAMBDA_1-5_32_100'
#phase = 'P'
master_dir = '/home/hosseini/Work/Scripts/gitHUB/MEASUREMENTS/Pdiff_measure_1_sec_LAMBDA_1-5_90_180'
phase = 'Pdiff'
req_bands = ['01', '02', '03', '04']
all_events = ['0274.2009.273.a']
#all_events = True

############## CRITERIA ###############
# WARNING: all the intervals are min <= ... < max
min_depth = -10
max_depth = 2000

min_xcorr = 0.8
max_xcorr = 1.01

#min_epi = 32
#max_epi = 85.01
min_epi = 97.0
max_epi = 180.01

check_clip = True
#######################################

selected_events_add = './info/selected_events_indexed.txt'
# -------------------------------------------------------

# ======================= PRINTING INPUTS
print '==============INPUT==================='
print 'Event DIR: %s' % master_dir
print 'Requested Phase: %s' % phase
print 'Requested Band: %s' % req_bands
if all_events == True:
    print 'Requested events: ALL'
else:
    print 'Requested events: %s' % all_events
print '============CRITERIA=============='
print '%s <= Depth < %s' % (min_depth, max_depth)
print '%s <= xcorr < %s' % (min_xcorr, max_xcorr)
print '%s <= epicentral < %s' % (min_epi, max_epi)
print 'Check clip: %s' % check_clip
print '==============END INPUT==================='
# ======================= END PRINTING INPUTS

print '\n\n================================================================='
print '!!!!!!!!!!!!!!!!!!!! PROGRAM START !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print '================================================================='

t_start = datetime.datetime.now()
print 'Start time: %s' % t_start

print '\n======>> reading the event information and filter them %s <= depth < %s' % (min_depth, max_depth)
passed_event_adds = outread.event_filter(master_dir, selected_events_add=selected_events_add, all_events=all_events,
                                         min_dp=min_depth, max_dp=max_depth)

print '\n======>> reading ALL FFM outputs'
all_output_files = []
for i in range(len(passed_event_adds)):
    output_files = np.array([])
    for rb in range(len(req_bands)):
        req_band = 'band%s' % req_bands[rb]
        tmp_output_file = outread.output_file_reader(evinfo=passed_event_adds[i], req_band=req_band)
        if tmp_output_file.size == 0:
            continue
        if output_files.size == 0:
            output_files = np.ones((tmp_output_file.shape[0], tmp_output_file.shape[1], len(req_bands)),
                                   dtype='object')
        output_files[:, :, rb] = tmp_output_file
    if not output_files.size == 0:
        all_output_files.append(output_files)

print '\n========================='
print 'All events   : %s' % len(passed_event_adds)
print 'Used events  : %s' % len(all_output_files)
print 'Missed events: %s' % (len(passed_event_adds) - len(all_output_files))
print '========================='

print '\n======>> filtering FFM outputs:'
print '%s <= xcorr < %s' % (min_xcorr, max_xcorr)
print '%s <= epicentral < %s' % (min_epi, max_epi)
print 'check clip: %s' % check_clip

all_filt_evsta = []
for i in range(len(all_output_files)):
    tmp_filt_array = outread.array_station_filter_mark(all_output_files[i], min_xcorr=min_xcorr, max_xcorr=max_xcorr,
                                                       min_epi=min_epi, max_epi=max_epi, check_clip=check_clip)
    all_filt_evsta.append(tmp_filt_array)
for i in range(len(all_filt_evsta)):
    outread.axi_kernel_receiver_writer(all_filt_evsta[i])
