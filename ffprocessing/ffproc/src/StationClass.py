#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  StationClass.py
#   Purpose:   Station class for ffprocessing
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import fnmatch
import glob
import math
import numpy as np
from obspy import read, Stream, UTCDateTime
try:
    from obspy.geodetics import locations2degrees
except Exception, e:
    from obspy.core.util import locations2degrees
try:
    from obspy.geodetics.base import gps2dist_azimuth as gps2DistAzimuth
except Exception, e:
    try:
        from obspy.geodetics import gps2DistAzimuth
    except Exception, e:
        from obspy.core.util import gps2DistAzimuth
from obspy.signal.interpolation import lanczos_interpolation
import os
import sys

from src.pres2dis import pres2dis
from src.utility_codes import bc, cprint
from src.utility_codes import geocen_array
from src.utility_codes import create_stxml_files
# XXX to not show the sampling warning...why?
# in our processing scheme, all the real waveforms are filtered (below 4Hz)
# and the synthetic waveforms have usually lower sampling rates as
# 10Hz (default); therefore, everything should be fine. XXX
suppress_sample_warn = True

# ------------------ read_data_header ---------------------------


def read_data_header(inp_ffproc, all_events, ev):
    """
    read the data header to create station class,
    this function only decides which method of reading should be used
    :param inp_ffproc:
    :param all_events:
    :param ev:
    :return:
    """
    cprint('StationClass.py', '[DATA]', bc.lgreen, 'Read the headers')

    if inp_ffproc.real_mode == 'read':
        station_cl = read_data_header_cat(inp_ffproc=inp_ffproc,
                                          all_events=all_events,
                                          ev=ev)

    elif inp_ffproc.real_mode == 'ami':
        station_cl = read_data_header_ami(inp_ffproc=inp_ffproc,
                                          all_events=all_events,
                                          ev=ev)
    elif inp_ffproc.real_mode == 'obs':
        station_cl = read_data_header_obs(inp_ffproc=inp_ffproc,
                                          all_events=all_events,
                                          ev=ev)
    elif inp_ffproc.real_mode == 'repo':
        station_cl = read_data_header_repo(inp_ffproc=inp_ffproc,
                                          all_events=all_events,
                                          ev=ev)
    else:
        cprint('StationClass.py', '[ERROR] [DATA]', bc.dred,
               '>>%s<< is not implemented! EXIT'
               % inp_ffproc.real_mode)
        import ipdb; ipdb.set_trace()
        sys.exit()
    return station_cl

# ------------------ read_data ---------------------------


def read_data(all_events, all_stations, inp_ffproc, ev):
    """
    reading mode for real data,
    this function only decides which method of reading should be used
    :param all_events:
    :param all_stations:
    :param inp_ffproc:
    :param ev:
    :return:
    """

    if inp_ffproc.real_mode == 'read':
        cprint('StationClass.py', '[DATA]', bc.lgreen,
               'Reading real waveforms.')
        real_waveforms = \
            read_real_cat(all_events, all_stations, inp_ffproc, ev)
    elif inp_ffproc.real_mode == 'obs':
        cprint('StationClass.py', '[DATA]', bc.lgreen,
               'Reading **raw** OBS waveforms.')
        real_waveforms = \
            read_raw_obs(all_events, all_stations, inp_ffproc, ev)
    elif inp_ffproc.real_mode == 'repo':
        cprint('StationClass.py', '[DATA]', bc.lgreen,
               'Reading **raw** OBS waveforms.')
        real_waveforms = \
            read_raw_repo(all_events, all_stations, inp_ffproc, ev)
    else:
        import ipdb; ipdb.set_trace()
        cprint('StationClass.py', '[ERROR] [DATA]', bc.dred,
               '>>%s<< is not implemented! EXIT'
               % inp_ffproc.real_mode)
        sys.exit()
    return real_waveforms

# ------------------ read_syn ---------------------------


def read_syn(all_events, all_stations, all_stfs, inp_ffproc, ev):
    """
    reading mode for synthetics,
    this function decides which method of reading should be used
    :param all_events:
    :param all_stations:
    :param inp_ffproc:
    :param ev:
    :return:
    """
    if inp_ffproc.syn_mode == 'read':
        cprint('StationClass.py', '[SYNTH] [YSPEC]', bc.lgreen,
               'Reading YSPEC synthetic waveforms!')
        raw_input("YSPEC (read) has not tested/cleaned-up!")
        syn_waveforms = read_syn_yspec(all_events=all_events,
                                       all_stations=all_stations,
                                       inp_ffproc=inp_ffproc,
                                       ev=ev)
    elif inp_ffproc.syn_mode == 'specfem':
        cprint('StationClass.py', '[SYNTH] [SPECFEM]', bc.lgreen,
               'Reading SPECFEM synthetic waveforms!')
        raw_input("SPECFEM has not tested/cleaned-up!")
        syn_waveforms = read_syn_specfem(all_events=all_events,
                                         all_stations=all_stations,
                                         inp_ffproc=inp_ffproc,
                                         ev=ev)
    elif inp_ffproc.syn_mode == 'instaseis':
        cprint('StationClass.py', '[SYNTH] [INSTASEIS]', bc.lgreen,
               'Instaseis synthetic waveforms!')
        syn_waveforms = read_syn_instaseis(all_events=all_events,
                                           all_stations=all_stations,
                                           all_stfs=all_stfs,
                                           inp_ffproc=inp_ffproc,
                                           ev=ev)
    else:
        cprint('StationClass.py', '[ERROR] [SYNTH]', bc.dred,
               '%s has not implemented!' % inp_ffproc.syn_mode)
        sys.exit()
    return syn_waveforms

# ------------------ read_data_header_cat ---------------------------


def read_data_header_cat(inp_ffproc, all_events, ev):
    """
    read data header from catalog and create a StationClass
    :param inp_ffproc:
    :param all_events:
    :param ev:
    :return:
    """
    station_cl = StationClass()
    # Loop over requested channels
    for cha in inp_ffproc.ph_channel:
        # import ipdb; ipdb.set_trace()
        real_data_path = glob.glob(os.path.join(inp_ffproc.real_path,
                                                all_events.name[ev],
                                                inp_ffproc.real_name_format,
                                                '*%s' % cha))
        real_data_path.sort()
        cprint('StationClass.py', '[DATA]', bc.lgreen,
               'Creating station class for %s.' % cha)
        list_stas = []
        # import ipdb;
        # ipdb.set_trace()
        for sta in real_data_path:
            # print sta
            try:
                try:
                    tr = read(sta, format='SAC', headonly=True)[0]
                    stla = tr.stats.sac.stla
                    stlo = tr.stats.sac.stlo
                    stdp = tr.stats.sac.stdp
                    stel = tr.stats.sac.stel
                except Exception, e:
                    tr = read(sta, format='MSEED', headonly=True)[0]
                    stla, stlo, stdp, stel = \
                        read_inv_dmt(inp_ffproc, all_events, ev, tr)
                    if not stla:
                        raise Exception
                name = '%s.%s.%s.%s' % (tr.stats.network, tr.stats.station,
                                        tr.stats.location, tr.stats.channel)
                geocen_stla = geocen_array(np.array([float(stla)]))[0]
                dist_deg = locations2degrees(all_events.geocen_lat[ev],
                                             all_events.lon[ev],
                                             geocen_stla,
                                             float(stlo))
                list_stas.append([dist_deg, sta, name, stla, stlo, stdp, stel,
                                  all_events.geocen_lat[ev], all_events.lon[ev]])
                list_stas = sorted(list_stas, key=lambda stas: stas[0])
            except Exception, e:
                cprint('read_data_header_cat', '[WARNING] [REAL]', bc.orange,
                       'Can not read %s, %s' % (sta, e))
                # import ipdb; ipdb.set_trace()
                continue

        for l_sta in list_stas:
            station_cl.addsta(l_sta[1], l_sta[2], l_sta[3], l_sta[4],
                              l_sta[5], l_sta[6], l_sta[7], l_sta[8])

    return station_cl


# ------------------ read_data_header_ami ---------------------------


def read_data_header_ami(inp_ffproc, all_events, ev):
    """
    read data header from catalog and create a StationClass
    :param inp_ffproc:
    :param all_events:
    :param ev:
    :return:
    """
    station_cl = StationClass()
    # Loop over requested channels
    for cha in inp_ffproc.ph_channel:
        # import ipdb; ipdb.set_trace()
        real_data_path = glob.glob(os.path.join(inp_ffproc.real_path,
                                                all_events.name[ev],
                                                inp_ffproc.real_name_format,
                                                '*%s' % cha))
        real_data_path.sort()
        cprint('StationClass.py', '[DATA]', bc.lgreen,
               'Creating station class.')
        list_stas = []
        # import ipdb;
        # ipdb.set_trace()
        for sta in real_data_path:
            try:
                try:
                    tr = read(sta, format='SAC', headonly=True)[0]
                    stla = tr.stats.sac.stla
                    stlo = tr.stats.sac.stlo
                    stdp = tr.stats.sac.stdp
                    stel = tr.stats.sac.stel
                except Exception, e:
                    tr = read(sta, format='MSEED', headonly=True)[0]
                    stla, stlo, stdp, stel = \
                        read_inv_dmt(inp_ffproc, all_events, ev, tr)
                    if not stla:
                        raise Exception
                name = '%s.%s.%s.%s' % (tr.stats.network, tr.stats.station,
                                        tr.stats.location, tr.stats.channel)
                geocen_stla = geocen_array(np.array([float(stla)]))[0]
                dist_deg = locations2degrees(all_events.geocen_lat[ev],
                                             all_events.lon[ev],
                                             geocen_stla,
                                             float(stlo))
                list_stas.append([dist_deg, sta, name, stla, stlo, stdp, stel,
                                  all_events.geocen_lat[ev], all_events.lon[ev]])
                list_stas = sorted(list_stas, key=lambda stas: stas[0])
            except Exception, e:
                cprint('read_data_header_ami', '[WARNING] [REAL]', bc.orange,
                       'Can not read %s, %s' % (sta, e))
                # import ipdb; ipdb.set_trace()
                continue

        for l_sta in list_stas:
            station_cl.addsta(l_sta[1], l_sta[2], l_sta[3], l_sta[4],
                              l_sta[5], l_sta[6], l_sta[7], l_sta[8])

    return station_cl


# ------------------ read_data_header_obs ---------------------------

def read_data_header_obs(inp_ffproc, all_events, ev):
    """
    read obs station information from file
    following information is necessary
    real_path, name, lat, lon, dp, el, evla, evlo
    :param inp_ffproc:
    :param all_events:
    :param ev:
    :return:
    """
    station_cl = StationClass()
    net = 'AI'
    cprint('StationClass.py', '[read_data_header_obs]', bc.lred,
            'ATTENTION - for OBS data the network name %s is in this part'
            'hard coded at this moment.' % net)

    # Loop over requested channels
    for cha in inp_ffproc.ph_channel:
        # this is only valid if we have OBS data
        real_data_path = glob.glob(inp_ffproc.real_name_format)
        cprint('StationClass.py', '[DATA]', bc.lgreen,
               'Creating station class.')
        # lat | lon | station_name | depth | elevation | NR | station code | station name
        sta_info = np.loadtxt(real_data_path[0], dtype=object)
        # import ipdb; ipdb.set_trace()
        list_stas = []

        for sta in sta_info:
            try:
                # import ipdb; ipdb.set_trace()
                name = '%s.%s.%s.%s' % (net, sta[6], '00', cha)
                geocen_stla = geocen_array(np.array([float(sta[0])]))[0]
                dist_deg = locations2degrees(all_events.geocen_lat[ev],
                                             all_events.lon[ev],
                                             geocen_stla,
                                             float(sta[1]))

                # XXX need better name for mseed folder and cha name XXX
                mseed_file_name = '%04d-%02d-%02d-%s.%s' % (all_events.year[ev], all_events.datetime[ev].month,
                                                            all_events.datetime[ev].day, 'HH'+cha[-1], 'mseed')

                sta_path = os.path.join(inp_ffproc.real_path, sta[6], 'mseed', mseed_file_name)
                cprint('read_data_header_obs', '[ATTENTION] [REAL]', bc.dgreen, 'This path is used for OBS data: %s' % sta_path)
                if glob.glob(sta_path) == []:
                    cprint('read_data_header_obs', '[WARNING] [REAL]', bc.orange, 'Path does not exisit: %s' % sta_path)
                    # if the station did not record this long - no reason to safe this information in the first place
                    continue
                list_stas.append([dist_deg, sta_path, name, sta[0], sta[1], sta[3], sta[4],
                                  all_events.geocen_lat[ev], all_events.lon[ev]])

                list_stas = sorted(list_stas, key=lambda stas: stas[0])

            except Exception, e:
                cprint('read_data_header_obs', '[WARNING] [REAL]', bc.orange,
                       'Can not read %s, %s' % (sta, e))
                continue

        for l_sta in list_stas:
            station_cl.addsta(l_sta[1], l_sta[2], l_sta[3], l_sta[4],
                              l_sta[5], l_sta[6], l_sta[7], l_sta[8])

    return station_cl


# ------------------ read_data_header_repo ---------------------------

def read_data_header_repo(inp_ffproc, all_events, ev):
    """
    read obs station information from file
    following information is necessary
    real_path, name, lat, lon, dp, el, evla, evlo
    :param inp_ffproc:
    :param all_events:
    :param ev:
    :return:
    """
    # import ipdb; ipdb.set_trace()
    station_cl = StationClass()

    # Loop over requested channels
    for cha in inp_ffproc.ph_channel:
        # this is only valid if we have OBS data
        real_data_path = glob.glob(inp_ffproc.real_name_format)
        cprint('StationClass.py', '[DATA]', bc.lgreen,
               'Creating station class.')
        # 0 - NET| 1 - STA| 2 - LAT| 3 - LON| 4 - ELEV 
        sta_info = np.loadtxt(real_data_path[0], delimiter=',', dtype=object)
        # import ipdb; ipdb.set_trace()
        list_stas = []

        for sta in sta_info:
            try:
                # import ipdb; ipdb.set_trace()
                name = '%s.%s.%s.%s' % (sta[0], sta[1], '', cha)
                geocen_stla = geocen_array(np.array([float(sta[2])]))[0]
                dist_deg = locations2degrees(all_events.geocen_lat[ev],
                                             all_events.lon[ev],
                                             geocen_stla,
                                             float(sta[3]))

                # XXX need better name for mseed folder and cha name XXX
                # IA.IA001..HHZ.D.2018.049
                mseed_file_name = '%s.%s..%s.D.%04d.%03d' % (sta[0], sta[1], cha, all_events.year[ev], 
                                                             all_events.datetime[ev].julday)
                # '/mnt/REPO/IA/MINISEED/YEAR/NET/STATION/CHANNEL.D/mseed_file_name'
                sta_path = os.path.join(inp_ffproc.real_path, str(all_events.year[ev]), sta[0], sta[1], cha + '.D', mseed_file_name)
                cprint('read_data_header_repo', '[ATTENTION] [REAL]', bc.dgreen, 'This path is used for REPO data: %s' % sta_path)
                if glob.glob(sta_path) == []:
                    cprint('read_data_header_repo', '[WARNING] [REAL]', bc.orange, 'Path does not exisit: %s' % sta_path)
                    # if the station did not record this long - no reason to safe this information in the first place
                    continue

                # dist_deg, sta_path, name, lat | lon |  depth | elevation,all_events.geocen_lat[ev], all_events.lon[ev]
                list_stas.append([dist_deg, sta_path, name, eval(sta[2]), eval(sta[3]), -eval(sta[4]), eval(sta[4]),
                                  all_events.geocen_lat[ev], all_events.lon[ev]])

                list_stas = sorted(list_stas, key=lambda stas: stas[0])

            except Exception, e:
                cprint('read_data_header_repo', '[WARNING] [REAL]', bc.orange,
                       'Can not read %s, %s' % (sta, e))
                continue

        for l_sta in list_stas:
            # addsta(self, real_path, name, lat, lon, dp, el, evla, evlo)
            
            station_cl.addsta(l_sta[1], l_sta[2], l_sta[3], l_sta[4],
                              l_sta[5], l_sta[6], l_sta[7], l_sta[8])

    return station_cl


# ------------------ read_inv_dmt ---------------------------


def read_inv_dmt(inp_ffproc, all_events, ev, tr):
    """
    read station_event file from DMT directories
    :param inp_ffproc:
    :param all_events:
    :param ev:
    :param tr:
    :return:
    """
    try:
        inv_info = np.loadtxt(os.path.join(inp_ffproc.real_path,
                                           all_events.name[ev],
                                           'info',
                                           'station_event'),
                              dtype='object', delimiter=',')
        inv_info_net = inv_info[inv_info[:, 0] == tr.stats.network]
        inv_info_sta = inv_info_net[inv_info_net[:, 1] == tr.stats.station]
        inv_info_loc = inv_info_sta[inv_info_sta[:, 2] == tr.stats.location]
        inv_info_cha = inv_info_loc[inv_info_loc[:, 3] == tr.stats.channel]
        if len(inv_info_sta) == 0:
            return False, False, False, False
        stla = inv_info_cha[0][4]
        stlo = inv_info_cha[0][5]
        stel = inv_info_cha[0][6]
        stdp = inv_info_cha[0][7]
    except Exception, exp:
        print exp
        import ipdb;ipdb.set_trace()
    return stla, stlo, stdp, stel

# ------------------ read_real_cat ---------------------------


def read_real_cat(all_events, all_stations, inp_ffproc, ev, lancz_a=12):
    """
    reading real data from catalog
    :param all_events:
    :param all_stations:
    :param inp_ffproc:
    :param ev:
    :param lancz_a:
    :return:
    """
    if inp_ffproc.filt_mode == 'log-gabor':
        max_filt_duration = inp_ffproc.lgb_filt_duration[0]
    else:
        cprint('StationClass.py', '[ERROR] [FILTER]', bc.dred,
               '%s filter mode has not implemented' % inp_ffproc.filt_mode)
        sys.exit()

    real_waveforms = Stream()
    # accurate duration of the seismogram
    acc_time = inp_ffproc.ph_preset + max_filt_duration + inp_ffproc.ph_offset
    # acc_dur: number of samples
    acc_dur = (inp_ffproc.ph_sampling_rate * acc_time) + 1
    acc_dur = int(acc_dur)

    #import ipdb; ipdb.set_trace()
    for cha in range(len(all_stations.name)):
        cha_id = all_stations.name[cha].split('.')
        real_name = '%s.%s.%s.%s' % (cha_id[0],
                                     cha_id[1],
                                     cha_id[2],
                                     cha_id[3])
        real_data_path = all_stations.real_path[cha]
        try:
            start_time = all_events.datetime[ev] + all_stations.tt_ph[cha]
            len_tr = max_filt_duration + inp_ffproc.ph_offset + 10.
            try:
                tr = read(real_data_path,
                          starttime=start_time-inp_ffproc.ph_preset,
                          endtime=start_time+len_tr,
                          format='SAC')[0]
            except Exception, e:
                tr = read(real_data_path,
                          starttime=start_time-inp_ffproc.ph_preset,
                          endtime=start_time+len_tr)[0]

            starttime_diff = tr.stats.starttime - \
                             (start_time - inp_ffproc.ph_preset)
            if not starttime_diff <= tr.stats.delta:
                cprint('StationClass.py', '[WARNING] [START-TIME]', bc.orange,
                       'starttime of the trace differs from the requested '
                       'starttime.')
                raise Exception

            if not suppress_sample_warn:
                if tr.stats.sampling_rate > inp_ffproc.ph_sampling_rate:
                    cprint('StationClass.py', '[WARNING] [SAMPLING]',
                           bc.orange,
                           'sampling rate of %s (%s) is higher than the '
                           'requested sampling rate. '
                           'This means down-sampling which is NOT properly '
                           'implemented yet.'
                           % (real_name, tr.stats.sampling_rate))

            starttime = tr.stats.starttime.timestamp
            endtime = tr.stats.endtime.timestamp
            dt = 1./inp_ffproc.ph_sampling_rate
            tr.data = lanczos_interpolation(
                np.require(tr.data, requirements=["C"]),
                starttime, tr.stats.delta, starttime, dt,
                int(math.floor((endtime - starttime) / dt)) + 1,
                a=lancz_a, window="blackman")

            tr.stats.starttime = UTCDateTime(starttime)
            tr.stats.sampling_rate = inp_ffproc.ph_sampling_rate

            tr.data = tr.data[:acc_dur]
            if not len(tr.data) == acc_dur:
                cprint('StationClass.py', '[ERROR] [DATA]', bc.dred,
                       "Length of the data: %s\nRequested length: %s"
                       % (len(tr.data), acc_dur))
                raise Exception
            if len(np.where(np.isnan(tr.data))[0]) > 0:
                raise Exception
            all_stations.real_rdable = \
                np.append(all_stations.real_rdable, True)
            # filling in the starttime of the real data relative to event time
            all_stations.real_start = np.append(all_stations.real_start,
                                                tr.stats.starttime -
                                                all_events.datetime[ev])
            real_waveforms.append(tr)
        except Exception, e:
            import ipdb;ipdb.set_trace()
            cprint('StationClass.py', '[WARNING] [REAL]', bc.orange,
                   'Can not read %s, %s' % (real_data_path, e))
            all_stations.real_rdable = \
                np.append(all_stations.real_rdable, False)
            all_stations.real_start = \
                np.append(all_stations.real_start, False)
    return real_waveforms


# ------------------ read_real_cat ---------------------------


def read_raw_obs(all_events, all_stations, inp_ffproc, ev, lancz_a=12):
    """
    reading raw obs data which still need to be instrument corrected on the fly
    :param all_events:
    :param all_stations:
    :param inp_ffproc:
    :param ev:
    :param lancz_a:
    :return:
    """
    # import ipdb; ipdb.set_trace()
    if inp_ffproc.filt_mode == 'log-gabor':
        max_filt_duration = inp_ffproc.lgb_filt_duration[0]
    else:
        cprint('StationClass.py', '[ERROR] [FILTER]', bc.dred,
               '%s filter mode has not implemented' % inp_ffproc.filt_mode)
        sys.exit()

    real_waveforms = Stream()
    # accurate duration of the seismogram
    acc_time = inp_ffproc.ph_preset + max_filt_duration + inp_ffproc.ph_offset
    # acc_dur: number of samples
    acc_dur = (inp_ffproc.ph_sampling_rate * acc_time) + 1
    acc_dur = int(acc_dur)

    # import ipdb; ipdb.set_trace()
    for cha in range(len(all_stations.name)):
        # cha_id = all_stations.name[cha].split('.')
        # real_name = '%s.%s.%s.%s' % (cha_id[0],
        #                              cha_id[1],
        #                              cha_id[2],
        #                              cha_id[3])

        real_data_path = all_stations.real_path[cha]
        # import ipdb; ipdb.set_trace()
        try:
            start_time = all_events.datetime[ev] + all_stations.tt_ph[cha]
            len_tr = max_filt_duration + inp_ffproc.ph_offset + 10.
            try:
                tr = read(real_data_path,
                          starttime=start_time-inp_ffproc.ph_preset,
                          endtime=start_time+len_tr,
                          format='SAC')[0]
            except Exception, e:
                tr = read(real_data_path,
                          starttime=start_time-inp_ffproc.ph_preset,
                          endtime=start_time+len_tr)[0]

            # resample - like the STF resampling
            starttime = start_time-inp_ffproc.ph_preset
            endtime = all_events.start_time+len_tr
            dt = 1. / inp_ffproc.ph_sampling_rate
            tr.data = lanczos_interpolation(np.require(tr.data, requirements=["C"]),
                                            starttime,
                                            tr.stats.delta,
                                            starttime,
                                            dt,
                                            int(math.floor((endtime - starttime) / dt)) + 1,
                                            a=lancz_a,
                                            window="blackman")

            starttime_diff = tr.stats.starttime - \
                             (start_time - inp_ffproc.ph_preset)

            if not starttime_diff <= tr.stats.delta:
                cprint('StationClass.py', '[WARNING] [START-TIME]', bc.orange,
                       'starttime of the trace differs from the requested '
                       'starttime.')
                raise Exception

            if not suppress_sample_warn:
                if tr.stats.sampling_rate > inp_ffproc.ph_sampling_rate:
                    cprint('StationClass.py', '[WARNING] [SAMPLING]',
                           bc.orange,
                           'sampling rate of %s (%s) is higher than the '
                           'requested sampling rate. '
                           'This means down-sampling which is NOT properly '
                           'implemented yet.'
                           % (all_stations.real_path[cha], tr.stats.sampling_rate))

            starttime = tr.stats.starttime.timestamp
            endtime = tr.stats.endtime.timestamp
            dt = 1./inp_ffproc.ph_sampling_rate
            tr.data = lanczos_interpolation(
                np.require(tr.data, requirements=["C"]),
                starttime, tr.stats.delta, starttime, dt,
                int(math.floor((endtime - starttime) / dt)) + 1,
                a=lancz_a, window="blackman")

            tr.stats.starttime = UTCDateTime(starttime)
            tr.stats.sampling_rate = inp_ffproc.ph_sampling_rate

            tr.data = tr.data[:acc_dur]
            if not len(tr.data) == acc_dur:
                cprint('StationClass.py', '[ERROR] [DATA]', bc.dred,
                       "Length of the data: %s\nRequested length: %s"
                       % (len(tr.data), acc_dur))
                raise Exception
            if len(np.where(np.isnan(tr.data))[0]) > 0:
                raise Exception
            all_stations.real_rdable = \
                np.append(all_stations.real_rdable, True)
            # filling in the starttime of the real data relative to event time
            all_stations.real_start = np.append(all_stations.real_start,
                                                tr.stats.starttime -
                                                all_events.datetime[ev])

            # XXX here the instrument correction for the real smgr has to happen
            # import ipdb; ipdb.set_trace()
            # import ipdb; ipdb.set_trace()
            cprint('StationClass.py', '[DATA] [read_raw_obs]', bc.lgreen,
                  'INSTRUMENT correction of OBS waveforms.')
            inv = create_stxml_files(tr, inp_ffproc.correction_seismometer, 
                                     inp_ffproc.correction_hydrophone, inp_ffproc.real_mode)
            tr.attach_response(inv)
            tr = tr.remove_response(output=inp_ffproc.trace_unit, pre_filt=inp.trace_pre_filt,
                                              water_level=inp_ffproc.trace_waterlevel,
                                              zero_mean=True,
                                              taper=True,
                                              taper_fraction=inp_ffproc.trace_taper)

            real_waveforms.append(tr)

        except Exception, e:
            cprint('StationClass.py', '[WARNING] [REAL]', bc.orange,
                   'Can not read %s, %s' % (real_data_path, e))
            all_stations.real_rdable = \
                np.append(all_stations.real_rdable, False)
            all_stations.real_start = \
                np.append(all_stations.real_start, False)
    return real_waveforms

# ------------------ read_real_cat ---------------------------


def read_raw_repo(all_events, all_stations, inp_ffproc, ev, lancz_a=12):
    """
    reading raw obs data which still need to be instrument corrected on the fly
    :param all_events:
    :param all_stations:
    :param inp_ffproc:
    :param ev:
    :param lancz_a:
    :return:
    """
    # import ipdb; ipdb.set_trace()
    if inp_ffproc.filt_mode == 'log-gabor':
        max_filt_duration = inp_ffproc.lgb_filt_duration[0]
    else:
        cprint('StationClass.py', '[ERROR] [FILTER]', bc.dred,
               '%s filter mode has not implemented' % inp_ffproc.filt_mode)
        sys.exit()

    real_waveforms = Stream()
    # accurate duration of the seismogram
    acc_time = inp_ffproc.ph_preset + max_filt_duration + inp_ffproc.ph_offset
    # acc_dur: number of samples
    acc_dur = (inp_ffproc.ph_sampling_rate * acc_time) + 1
    acc_dur = int(acc_dur)

    # import ipdb; ipdb.set_trace()
    for cha in range(len(all_stations.name)):
        real_data_path = all_stations.real_path[cha]

        try:
            start_time = all_events.datetime[ev] + all_stations.tt_ph[cha]
            len_tr = max_filt_duration + inp_ffproc.ph_offset + 10.
            try:
                tr = read(real_data_path,
                          starttime=start_time-inp_ffproc.ph_preset,
                          endtime=start_time+len_tr,
                          format='SAC')[0]
            except Exception, e:
                tr = read(real_data_path,
                          starttime=start_time-inp_ffproc.ph_preset,
                          endtime=start_time+len_tr)[0]

            starttime_diff = tr.stats.starttime - \
                             (start_time - inp_ffproc.ph_preset)

            if not starttime_diff <= tr.stats.delta:
                cprint('StationClass.py', '[WARNING] [START-TIME]', bc.orange,
                       'starttime of the trace differs from the requested '
                       'starttime.')
                raise Exception

            if not suppress_sample_warn:
                if tr.stats.sampling_rate > inp_ffproc.ph_sampling_rate:
                    cprint('StationClass.py', '[WARNING] [SAMPLING]',
                           bc.orange,
                           'sampling rate of %s (%s) is higher than the '
                           'requested sampling rate. '
                           'This means down-sampling which is NOT properly '
                           'implemented yet.'
                           % (all_stations.real_path[cha], tr.stats.sampling_rate))

            starttime = tr.stats.starttime.timestamp
            endtime = tr.stats.endtime.timestamp
            dt = 1./inp_ffproc.ph_sampling_rate
            tr.data = lanczos_interpolation(
                np.require(tr.data, requirements=["C"]),
                starttime, tr.stats.delta, starttime, dt,
                int(math.floor((endtime - starttime) / dt)) + 1,
                a=lancz_a, window="blackman")

            tr.stats.starttime = UTCDateTime(starttime)
            tr.stats.sampling_rate = inp_ffproc.ph_sampling_rate

            tr.data = tr.data[:acc_dur]
            if not len(tr.data) == acc_dur:
                cprint('StationClass.py', '[ERROR] [DATA]', bc.dred,
                       "Length of the data: %s\nRequested length: %s"
                       % (len(tr.data), acc_dur))
                raise Exception
            if len(np.where(np.isnan(tr.data))[0]) > 0:
                raise Exception
            all_stations.real_rdable = \
                np.append(all_stations.real_rdable, True)
            # filling in the starttime of the real data relative to event time
            all_stations.real_start = np.append(all_stations.real_start,
                                                tr.stats.starttime -
                                                all_events.datetime[ev])

            # XXX here the instrument correction for the real smgr has to happen
            # import ipdb; ipdb.set_trace()
            # import ipdb; ipdb.set_trace()
            cprint('StationClass.py', '[DATA] [read_raw_repo]', bc.lgreen,
                  'INSTRUMENT correction of REPO waveforms.')
            inv = create_stxml_files(tr, inp_ffproc.correction_seismometer, 
                                     inp_ffproc.correction_hydrophone, inp_ffproc.real_mode)
            tr.attach_response(inv)
            tr = tr.remove_response(output=inp_ffproc.trace_unit,
                                              water_level=inp_ffproc.trace_waterlevel,
                                              zero_mean=True,
                                              taper=True,
                                              taper_fraction=inp_ffproc.trace_taper)

            real_waveforms.append(tr)

        except Exception, e:
            cprint('StationClass.py', '[WARNING] [REAL]', bc.orange,
                   'Can not read %s, %s' % (real_data_path, e))
            all_stations.real_rdable = \
                np.append(all_stations.real_rdable, False)
            all_stations.real_start = \
                np.append(all_stations.real_start, False)
    return real_waveforms



# ------------------ read_syn_yspec ---------------------------


def read_syn_yspec(all_events, all_stations, inp_ffproc, ev, lancz_a=20):
    """
    reading already saved YSPEC synthetic waveforms
    :param all_events:
    :param all_stations:
    :param inp_ffproc:
    :param ev:
    :param lancz_a:
    :return:
    """
    if inp_ffproc.filt_mode == 'log-gabor':
        max_filt_duration = inp_ffproc.lgb_filt_duration[0]
    else:
        cprint('StationClass.py', '[ERROR] [FILTER]', bc.dred,
               '%s filter mode has not implemented! EXIT!'
               % inp_ffproc.filt_mode)
        sys.exit()

    cprint('StationClass.py', '[YSPEC]', bc.lgreen, 'Read the synthetic data')

    syn_waveforms = Stream()
    # accurate duration of the seismogram
    acc_time = inp_ffproc.ph_preset + max_filt_duration + inp_ffproc.ph_offset
    # acc_dur: number of samples
    acc_dur = (inp_ffproc.ph_sampling_rate * acc_time) + 1
    acc_dur = int(acc_dur)

    for cha in range(len(all_stations.name)):
        cha_id = all_stations.name[cha].split('.')
        # Special name for YSPEC synthetic waveforms
        syn_name = 'grf.%s.%s.%s.x00.%s' % (cha_id[0],
                                            cha_id[1],
                                            cha_id[2],
                                            cha_id[3])
        syn_data_path = os.path.join(inp_ffproc.syn_path,
                                     all_events.name[ev],
                                     inp_ffproc.syn_name_format,
                                     syn_name)
        all_stations.syn_path = np.append(all_stations.syn_path, syn_data_path)
        try:
            tr_orig = read(syn_data_path,
                           headonly=True,
                           format='SAC')[0]
            if not abs(all_events.datetime[ev]-tr_orig.stats.starttime) < 0.1:
                cprint('StationClass.py', '[ERROR] [YSPEC]', bc.dred,
                       '[YSPEC] the starttime of the trace differs '
                       'from the event datetime! EXIT!')
                sys.exit()
            start_time = all_events.datetime[ev] + all_stations.tt_ph[cha]
            len_tr = max_filt_duration + inp_ffproc.ph_offset + 10.
            tr = read(syn_data_path,
                      starttime=start_time-inp_ffproc.ph_preset,
                      endtime=start_time+len_tr,
                      format='SAC')[0]
            # re-sample YSPEC data
            if not suppress_sample_warn:
                if tr.stats.sampling_rate > inp_ffproc.ph_sampling_rate:
                    cprint('StationClass.py', '[WARNING] [SAMPLING]',
                           bc.orange,
                           'sampling rate of %s (%s) is higher than the '
                           'requested sampling rate. This means down-sampling '
                           'which is NOT properly implemented yet.'
                           % (syn_name, tr.stats.sampling_rate))

            starttime = tr.stats.starttime.timestamp
            endtime = tr.stats.endtime.timestamp
            dt = 1./inp_ffproc.ph_sampling_rate
            tr.data = lanczos_interpolation(
                tr.data, starttime, tr.stats.delta, starttime, dt,
                int(math.floor((endtime - starttime) / dt)) + 1,
                a=lancz_a, window="blackman")

            tr.stats.sampling_rate = inp_ffproc.ph_sampling_rate
            tr.data = tr.data[:acc_dur]
            if not len(tr.data) == acc_dur:
                cprint('StationClass.py', '[ERROR] [DATA]', bc.dred,
                       "Length of the data: %s\nRequested length: %s"
                       % (len(tr.data), acc_dur))
                raise Exception
            all_stations.srate = np.append(all_stations.srate,
                                           tr.stats.sampling_rate)
            all_stations.syn_rdable = np.append(all_stations.syn_rdable, True)
            # filling in the starttime of the synthetic relative to event time
            all_stations.syn_start = np.append(all_stations.syn_start,
                                               tr.stats.starttime -
                                               all_events.datetime[ev])
            syn_waveforms.append(tr)
        except Exception, e:
            cprint('StationClass.py', '[WARNING] [YSPEC]', bc.orange,
                   'Can not read %s, %s' % (syn_data_path, e))
            all_stations.srate = np.append(all_stations.srate, False)
            all_stations.syn_rdable = np.append(all_stations.syn_rdable, False)
            all_stations.syn_start = np.append(all_stations.syn_start, False)
    return syn_waveforms

# ------------------ read_syn_specfem ---------------------------


def read_syn_specfem(all_events, all_stations, inp_ffproc, ev, lancz_a=20):
    """
    reading already saved SPECFEM synthetic waveforms
    :param all_events:
    :param all_stations:
    :param inp_ffproc:
    :param ev:
    :param lancz_a:
    :return:
    """
    if inp_ffproc.filt_mode == 'log-gabor':
        max_filt_duration = inp_ffproc.lgb_filt_duration[0]
    else:
        cprint('StationClass.py', '[ERROR] [FILTER]', bc.dred,
               '%s filter mode has not implemented! EXIT'
               % inp_ffproc.filt_mode)
        sys.exit()

    cprint('StationClass.py', '[SPECFEM]', bc.lgreen,
           'Read the synthetic data')
    syn_waveforms = Stream()
    # accurate duration of the seismogram
    acc_time = inp_ffproc.ph_preset + max_filt_duration + inp_ffproc.ph_offset
    # acc_dur: number of samples
    acc_dur = (inp_ffproc.ph_sampling_rate * acc_time) + 1
    acc_dur = int(acc_dur)

    for cha in range(len(all_stations.name)):
        cha_id = all_stations.name[cha].split('.')
        # Special name for SPECEM synthetic waveforms
        # Note that this might be replaced with MX for low res synthetics
        syn_name = 'dis.%s.%s..BH%s' % (cha_id[0],
                                        cha_id[1],
                                        cha_id[3].split('BX')[1])
        syn_data_path = os.path.join(inp_ffproc.syn_path,
                                     all_events.name[ev],
                                     inp_ffproc.syn_name_format,
                                     syn_name)
        all_stations.syn_path = np.append(all_stations.syn_path, syn_data_path)
        try:
            try:
                tr_orig = read(syn_data_path,
                               headonly=True,
                               format='SAC')[0]
            except Exception, e:
                tr_orig = read(syn_data_path,
                               headonly=True)[0]

            start_time = all_events.datetime[ev] + all_stations.tt_ph[cha]
            len_tr = max_filt_duration + inp_ffproc.ph_offset + 10.
            try:
                tr = read(syn_data_path,
                          starttime=start_time-inp_ffproc.ph_preset,
                          endtime=start_time+len_tr,
                          format='SAC')[0]
            except Exception, e:
                tr = read(syn_data_path,
                          starttime=start_time-inp_ffproc.ph_preset,
                          endtime=start_time+len_tr)[0]

            # resample SPECFEM data
            if not suppress_sample_warn:
                if tr.stats.sampling_rate > inp_ffproc.ph_sampling_rate:
                    cprint('StationClass.py', '[WARNING] [SAMPLING]',
                           bc.orange, 'sampling rate of %s (%s) is higher '
                                      'than the requested sampling rate. '
                                      'This means down-sampling which is NOT '
                                      'properly implemented yet.'
                           % (syn_name, tr.stats.sampling_rate))

            starttime = tr.stats.starttime.timestamp
            endtime = tr.stats.endtime.timestamp
            dt = 1./inp_ffproc.ph_sampling_rate
            tr.data = lanczos_interpolation(
                tr.data, starttime, tr.stats.delta, starttime, dt,
                int(math.floor((endtime - starttime) / dt)) + 1,
                a=lancz_a, window="blackman")

            tr.stats.sampling_rate = inp_ffproc.ph_sampling_rate
            tr.data = tr.data[:acc_dur]
            if not len(tr.data) == acc_dur:
                cprint('StationClass.py', '[ERROR] [DATA]', bc.dred,
                       "Length of the data: %s\n"
                       "Requested length: %s" % (len(tr.data), acc_dur))
                raise Exception
            all_stations.srate = np.append(all_stations.srate,
                                           tr.stats.sampling_rate)
            all_stations.syn_rdable = np.append(all_stations.syn_rdable, True)
            # filling in the starttime of the synthetic relative to event time
            all_stations.syn_start = np.append(all_stations.syn_start,
                                               tr.stats.starttime -
                                               all_events.datetime[ev])
            syn_waveforms.append(tr)
        except Exception, e:
            cprint('StationClass.py', '[WARNING] [SPECFEM]', bc.orange,
                   'Can not read %s, %s' % (syn_data_path, e))
            all_stations.srate = np.append(all_stations.srate, False)
            all_stations.syn_rdable = np.append(all_stations.syn_rdable, False)
            all_stations.syn_start = np.append(all_stations.syn_start, False)
    return syn_waveforms

# ------------------ read_syn_instaseis ---------------------------


def read_syn_instaseis(all_events, all_stations, all_stfs, inp_ffproc, ev,
                       lancz_a=12):
    """
    Creating instaseis synthetic waveforms
    :param all_events:
    :param all_stations:
    :param inp_ffproc:
    :param ev:
    :param lancz_a:
    :return:
    """
    import instaseis
    if inp_ffproc.filt_mode == 'log-gabor':
        max_filt_duration = inp_ffproc.lgb_filt_duration[0]
    else:
        cprint('StationClass.py', '[ERROR] [FILTER]', bc.dred,
               '%s filter mode has not implemented.'
               % inp_ffproc.filt_mode)
        sys.exit()

    cprint('StationClass.py', '[INSTASEIS]', bc.lgreen,
           'Read the synthetic data.')

    db_insta = instaseis.open_db(inp_ffproc.syn_path)
    # XXX these can be removed in case that I find a good solution for STF
    # convolution and deconvolution!
    all_stations.db_slip = db_insta.info.slip
    all_stations.db_sampling_rate = db_insta.info.sampling_rate
    all_stations.db_src_shift = db_insta.info.src_shift

    syn_waveforms = Stream()
    # accurate duration of the seismogram
    acc_time = inp_ffproc.ph_preset + max_filt_duration + inp_ffproc.ph_offset
    # acc_dur: number of samples
    acc_dur = (inp_ffproc.ph_sampling_rate * acc_time) + 1
    acc_dur = int(acc_dur)

    all_stfs.stfs.interpolate(sampling_rate=db_insta.info.sampling_rate,
                              interpolation='lanczos')

    add_save = os.path.join(inp_ffproc.output_dir, all_events.name[ev],
                            'outfiles')

    if not os.path.isdir(add_save):
        os.makedirs(add_save)
    stfk_fio = open(os.path.join(add_save, 'stf.dat'), 'w')
    stfk_fio.writelines('%s\n' % all_stfs.address_stf[0])
    stfk_fio.writelines('%i\n' % len(all_stfs.stfs))
    for istf in range(len(all_stfs.stfs)):
        stfk_fio.writelines('stf%02i  %s  %s  %s\n'
                            % (istf+1, len(all_stfs.stfs[istf]),
                            1./db_insta.info.sampling_rate,
                               all_stfs.taxs[istf][0]))
        for ic in range(len(all_stfs.stfs[istf].data)):
            stfk_fio.writelines('%s\n' % all_stfs.stfs[istf].data[ic])

    cprint('StationClass.py', '[INSTASEIS]', bc.green,
           "# Stations to be extracted from Instaseis database: %s"
           % len(all_stations.name))
    for cha in range(len(all_stations.name)):
        if cha > 0 and np.mod(cha, 100) == 0:
            cprint('StationClass.py', '[INSTASEIS]', bc.green,
                   "#extracted stations: %s/%s" % (cha, len(all_stations.name)))
        try:

            cha_id = all_stations.name[cha].split('.')
                
            receiver = instaseis.Receiver(
                latitude=all_stations.geocen_lat[cha],
                longitude=all_stations.lon[cha],
                network=cha_id[0],
                station=cha_id[1])

            if cha_id[3][-1] == 'Z' or cha_id[3][-1] == 'H':
                tr_comp = 'Z'
            elif cha_id[3][-1] == 'T':
                tr_comp = 'T'
            elif cha_id[3][-1] == 'R':
                tr_comp = 'R'
            # the next two are for the case of particle motion really XXX just in case
            # check if there any other possibilities where this could go wrong XXX
            #cprint('StationClass.py', '[WARNING] [CHANNELS]', bc.red, 
            #       'Check how channels X, Y, 1, 2, N, E are related!!!')
            
            elif cha_id[3] == 'BHX' or cha_id[3] == 'BH2':
                tr_comp = 'E'
            elif cha_id[3] == 'BHY' or cha_id[3] == 'BH1':
                tr_comp = 'N'
            elif cha_id[3] == 'BHE':
                tr_comp = 'E'
            elif cha_id[3] == 'BHN':
                tr_comp = 'N'
            else:
                import ipdb; ipdb.set_trace()
                sys.exit("This channel is not supported: %s" % cha_id[3])
            source = instaseis.Source(
                          latitude=all_events.geocen_lat[ev],
                          longitude=all_events.lon[ev],
                          depth_in_m=all_events.inv_dp[ev]*1000.,
                          m_rr=all_events.mrr[ev],
                          m_tt=all_events.mtt[ev],
                          m_pp=all_events.mpp[ev],
                          m_rt=all_events.mrt[ev],
                          m_rp=all_events.mrp[ev],
                          m_tp=all_events.mtp[ev],
                          origin_time=all_events.datetime[ev],
                sliprate=all_stfs.stfs[all_stfs.grps[cha]-1].data,
                dt=all_stfs.stfs[all_stfs.grps[cha]-1].stats.delta,
                time_shift=all_stfs.taxs[all_stfs.grps[cha]-1][0])

            # remove_source_shift means to cut the waveform according to
            # the source shift AND change the datetime. This does not
            # have any effect on the waveform measurement.
            tr = db_insta.get_seismograms(
                source, receiver,
                components=tr_comp,
                kind='displacement',
                remove_source_shift=False,
                reconvolve_stf=True,
                dt=1./inp_ffproc.ph_sampling_rate)[0]
            starttime = tr.stats.starttime.timestamp
            endtime = tr.stats.endtime.timestamp
            dt = 1./inp_ffproc.ph_sampling_rate
            tr.data = lanczos_interpolation(
                np.require(tr.data, requirements=["C"]),
                starttime, tr.stats.delta, starttime, dt,
                int(math.floor((endtime - starttime) / dt)) + 1,
                a=lancz_a, window="blackman")
            tr.stats.starttime = UTCDateTime(starttime)
            tr.stats.sampling_rate = inp_ffproc.ph_sampling_rate

            start_time = all_events.datetime[ev] + all_stations.tt_ph[cha]
            len_tr = max_filt_duration + inp_ffproc.ph_offset + 10.
            tr = tr.slice(starttime=start_time-inp_ffproc.ph_preset,
                          endtime=start_time+len_tr)
            tr.data = tr.data[:acc_dur]
            if not len(tr.data) == acc_dur:
                cprint('StationClass.py', '[ERROR] [DATA]', bc.dred,
                       "Length of the data: %s\n"
                       "Requested length: %s" % (len(tr.data), acc_dur))
                raise Exception
            if len(np.where(np.isnan(tr.data))[0]) > 0:
                raise Exception
            all_stations.srate = np.append(all_stations.srate,
                                           tr.stats.sampling_rate)
            all_stations.syn_rdable = np.append(all_stations.syn_rdable, True)
            tr.stats.network = cha_id[0]
            tr.stats.station = cha_id[1]
            tr.stats.location = cha_id[2]
            if inp_ffproc.proc_hydrophone and cha_id[3] == 'BDH':
                # XXX Check the name 
                tr.stats.channel = 'BDH'
            else:
                tr.stats.channel = cha_id[3]
            # filling in the starttime of the synthetic relative to event time
            all_stations.syn_start = np.append(all_stations.syn_start,
                                               tr.stats.starttime -
                                               all_events.datetime[ev])
            syn_waveforms.append(tr)
        except Exception, e:
            cprint('StationClass.py', '[WARNING] [INSTASEIS]', bc.orange,
                   'Can not read %s, %s' % (all_stations.name[cha], e))
            all_stations.srate = np.append(all_stations.srate, False)
            all_stations.syn_rdable = np.append(all_stations.syn_rdable, False)
            all_stations.syn_start = np.append(all_stations.syn_start, False)
            # XXX XXX 
            # import ipdb; ipdb.set_trace()
    return syn_waveforms

# ------------------ filter_stations ---------------------------


def filter_stations(inp_ffproc, all_stations):
    """
    filter stations based on the specified criteria
    :param inp_ffproc:
    :param all_stations:
    :return:
    """
    cprint('StationClass.py', '[DATA]', bc.lgreen,
           'Filtering %s stations' % all_stations.number_stations)
    min_lon, max_lon, min_lat, max_lat = inp_ffproc.real_rect.split('/')

    # import ipdb; ipdb.set_trace()
    indxs = all_stations.dist >= inp_ffproc.ph_min_epi
    indxs = np.multiply(indxs, all_stations.dist <= inp_ffproc.ph_max_epi)
    indxs = np.multiply(indxs, all_stations.lon >= float(min_lon))
    indxs = np.multiply(indxs, all_stations.lon <= float(max_lon))
    indxs = np.multiply(indxs, all_stations.lat >= float(min_lat))
    indxs = np.multiply(indxs, all_stations.lat <= float(max_lat))
    indxs = np.multiply(indxs, all_stations.azi >= inp_ffproc.real_min_azi)
    indxs = np.multiply(indxs, all_stations.azi <= inp_ffproc.real_max_azi)
    indxs_ones = np.zeros(len(all_stations.name))

    if inp_ffproc.real_id_state == 0:
        # import ipdb; ipdb.set_trace()
        indxs_name = np.where(np.in1d(
            all_stations.name, fnmatch.filter(all_stations.name,
                                              inp_ffproc.real_id)))[0]
        indxs_ones[indxs_name] = 1
        indxs = np.multiply(indxs, indxs_ones.astype(np.bool))
        all_stations.rmsta(indxs)

    if inp_ffproc.real_id_state == 1:
        for id in inp_ffproc.real_id_list:
            indxs_name = np.where(np.in1d(
                all_stations.name, fnmatch.filter(all_stations.name, id)))[0]
            indxs_ones[indxs_name] = 1
        indxs = np.multiply(indxs, indxs_ones.astype(np.bool))
        all_stations.rmsta(indxs)

    if inp_ffproc.real_id_state == 2:
        indxs_ones = np.ones(len(all_stations.name))
        for id in inp_ffproc.real_id_list:
            indxs_name = np.where(np.in1d(
                all_stations.name, fnmatch.filter(all_stations.name, id)))[0]
            indxs_ones[indxs_name] = 0
        indxs = np.multiply(indxs, indxs_ones.astype(np.bool))
        all_stations.rmsta(indxs)

    cprint('StationClass.py', '[DATA]', bc.dgreen,
           '# stations: %s (after filtering)' % all_stations.number_stations)

# ------------------ StationClass ---------------------------


class StationClass:
    def __init__(self):
        self.number_stations = 0

        self.nam = np.array([])
        self.loc = np.array([])
        self.net = np.array([])
        self.cha = np.array([])

        self.name = np.array([])
        self.lat = np.array([])
        self.geocen_lat = np.array([])
        self.lon = np.array([])
        self.dp = np.array([])
        self.el = np.array([])
        self.dist = np.array([])
        self.azi = np.array([])
        self.bazi = np.array([])
        self.tt_ph = np.array([])
        self.srate = np.array([])
        self.syn_rdable = np.array([])
        self.syn_path = np.array([])
        self.syn_start = np.array([])
        self.real_rdable = np.array([])
        self.real_path = np.array([])
        self.real_start = np.array([])
        self.ss_shift_step1 = np.array([])
        self.ss_shift_step2 = np.array([])
        self.cc_step1 = np.array([])
        self.cc_step2 = np.array([])
        self.clip_step1 = np.array([])
        self.clip_step2 = np.array([])
        self.final_time_shift = np.array([])
        self.final_amp = np.array([])
        self.final_two_sigma = np.array([])
        self.cut_sample_real = np.array([])
        self.cut_sample_syn = np.array([])
        self.length_bp = np.array([])
        self.snr = np.array([])

    def addsta(self, real_path, name, lat, lon, dp, el, evla, evlo):
        """
        Add station to the StationClass
        :param real_path:
        :param name:
        :param lat:
        :param lon:
        :param dp:
        :param el:
        :param evla:
        :param evlo:
        :return:
        """
        self.name = np.append(self.name, name)
        net, nam, loc, cha = name.split('.')
        self.nam = np.append(self.nam, nam)
        self.loc = np.append(self.loc, loc)
        self.net = np.append(self.net, net)
        self.cha = np.append(self.cha, cha)

        self.lat = np.append(self.lat, float(lat))
        geocen_stla = geocen_array(np.array([float(lat)]))[0]
        self.geocen_lat = np.append(self.geocen_lat, geocen_stla)
        self.lon = np.append(self.lon, float(lon))
        self.dp = np.append(self.dp, float(dp))
        self.el = np.append(self.el, float(el))
        (dist_m, azi, bazi) = gps2DistAzimuth(evla, evlo,
                                              float(geocen_stla), float(lon))
        dist_deg = locations2degrees(evla, evlo, float(geocen_stla), float(lon))
        self.dist = np.append(self.dist, float(dist_deg))
        self.azi = np.append(self.azi, float(azi))
        self.bazi = np.append(self.bazi, float(bazi))
        self.real_path = np.append(self.real_path, real_path)
        self.number_stations += 1

    def add_ttime(self, req_phase, evdp, bg_model):
        """
        Add theoretical arrival times to the StationClass
        :param req_phase:
        :param evdp:
        :param bg_model:
        :return:
        """
        if ',' in req_phase:
            list_req_phase = []
            for r_ph in req_phase.split(','):
                list_req_phase.append(r_ph.strip())
            req_phase = tuple(list_req_phase)
        else:
            req_phase = tuple([req_phase])

        cprint('StationClass.py', '[PHASE]', bc.green,
               'finding %s arrival times of %s stations'
               % (req_phase, self.number_stations))

        try:
            from obspy.taup import tau
            tau_bg = tau.TauPyModel(model=bg_model)
            for sta in range(len(self.lat)):
                self.tt_calc(self.dist[sta], req_phase, evdp,
                             bg_model, tau_bg)
        except Exception, e:
            cprint('StationClass.py', '[WARINING]', bc.orange,
                   'from obspy.taup import tau ---> Failed'
                   '\n WARNING: %s\n'
                   'getTravelTime will be used instead - '
                   'only the first phase will be used: %s' % (e, req_phase[0]))
            req_phase = req_phase[0]
            for sta in range(len(self.lat)):
                self.tt_calc(self.dist[sta], req_phase, evdp,
                             bg_model, tau_bg=False)

    def rmsta(self, indxs):
        """
        Remove stations from the StationClass
        :param indxs:
        :return:
        """
        self.name = self.name[indxs]

        self.nam = self.nam[indxs]
        self.loc = self.loc[indxs]
        self.net = self.net[indxs]
        self.cha = self.cha[indxs]

        self.lat = self.lat[indxs]
        self.geocen_lat = self.geocen_lat[indxs]
        self.lon = self.lon[indxs]
        self.dp = self.dp[indxs]
        self.el = self.el[indxs]
        self.dist = self.dist[indxs]
        self.azi = self.azi[indxs]
        self.bazi = self.bazi[indxs]
        if len(self.tt_ph) > 0:
            self.tt_ph = self.tt_ph[indxs]
        if len(self.srate) > 0:
            self.srate = self.srate[indxs]
        if len(self.syn_rdable) > 0:
            self.syn_rdable = self.syn_rdable[indxs]
        if len(self.syn_path) > 0:
            self.syn_path = self.syn_path[indxs]
        if len(self.real_path) > 0:
            self.real_path = self.real_path[indxs]
        if len(self.syn_start) > 0:
            self.syn_start = self.syn_start[indxs]
        if len(self.real_rdable) > 0:
            self.real_rdable = self.real_rdable[indxs]
        if len(self.real_start) > 0:
            self.real_start = self.real_start[indxs]
        if len(self.ss_shift_step1) > 0:
            self.ss_shift_step1 = self.ss_shift_step1[indxs]
        if len(self.ss_shift_step2) > 0:
            self.ss_shift_step2 = self.ss_shift_step2[indxs]
        if len(self.cc_step1) > 0:
            self.cc_step1 = self.cc_step1[indxs]
        if len(self.cc_step2) > 0:
            self.cc_step2 = self.cc_step2[indxs]
        if len(self.clip_step1) > 0:
            self.clip_step1 = self.clip_step1[indxs]
        if len(self.clip_step2) > 0:
            self.clip_step2 = self.clip_step2[indxs]
        if len(self.final_time_shift) > 0:
            self.final_time_shift = self.final_time_shift[indxs]
        if len(self.final_amp) > 0:
            self.final_amp = self.final_amp[indxs]
        if len(self.final_two_sigma) > 0:
            self.final_two_sigma = self.final_two_sigma[indxs]
        if len(self.cut_sample_real) > 0:
            self.cut_sample_real = self.cut_sample_real[indxs]
        if len(self.cut_sample_syn) > 0:
            self.cut_sample_syn = self.cut_sample_syn[indxs]
        if len(self.length_bp) > 0:
            self.length_bp = self.length_bp[indxs]
        if len(self.snr) > 0:
            self.snr = self.snr[indxs]
        self.number_stations -= len(indxs[indxs == False])

    def tt_calc(self, dist, req_phase, evdp, bg_model, tau_bg):
        """
        Calculate the theoretical arrival times
        :param dist:
        :param req_phase:
        :param evdp:
        :param bg_model:
        :param tau_bg:
        :return:
        """
        if not tau_bg:
            from obspy.taup.taup import getTravelTimes
            try:
                tt = [_i for _i in
                      getTravelTimes(dist, evdp, bg_model) if
                      req_phase == _i['phase_name']][0]['time']
            except Exception, e:
                cprint('StationClass.py', '[ATTENTION] [PHASE]', bc.orange,
                       '%s can not be found in %s degree distance! %s'
                       % (req_phase, round(dist, 2), e))
                tt = False
            self.tt_ph = np.append(self.tt_ph, tt)
        else:
            try:
                tt = tau_bg.get_travel_times(evdp, dist,
                                             phase_list=req_phase)[0].time
                # tt = round(tt, 1)
            except Exception, e:
                cprint('StationClass.py', '[ATTENTION] [PHASE]', bc.orange,
                       '%s can not be found in %s degree distance! %s'
                       % (req_phase, round(dist, 2), e))
                tt = False
            self.tt_ph = np.append(self.tt_ph, tt)

# Using TauP from the terminal...maybe not the best idea, time-wise!
# taup_process = subprocess.Popen(['taup_time', '-mod', bg_model,
#                                  '-time',
#                                  '-h', str(evdp),
#                                  '-ph', req_phase,
#                                  '-deg', str(dist)],
#                                 stdout=subprocess.PIPE)

# tt_raw = taup_process.communicate()[0]
# try:
#     tt = tt_raw.split('\n')[0].split()[-1]
#     tt = float(tt)
# except Exception, e:
#     logging.warn('%s can not be found in %s degree distance! %s'
#                  % (req_phase, dist, e))
#     tt = False


