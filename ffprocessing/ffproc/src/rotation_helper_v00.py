#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  rotation_helper.py
#   Purpose:   ?H? and ?H? components to R and T
#   Author:    Maria Tsekhmistrenko
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import matplotlib.image as image
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import numpy as np
from obspy.signal.rotate import rotate_ne_rt, rotate2zne
from obspy import read, read_inventory, Stream, UTCDateTime
from obspy.io.xseed import Parser
import glob
import os
import sys
import math

from obspy import Stream
from obspy.signal.interpolation import lanczos_interpolation
from obspy.signal.polarization import flinn, particle_motion_odr
from src.utility_codes import bc, cprint, progress_bar


suppress_sample_warn = True

# ------------- particle motion -----------------

def calc_part_mot(all_stations, real_waveforms, realBP, synBP,
                  inp_ffproc, all_events, ev):

    cprint('rotation_helper.py', '[PARTICLE MOTION]', bc.blue,
           'Analysing particle motion for this event for all stations.')

    # making a special folder for these plots in the event folder
    pm_output = os.path.join(inp_ffproc.output_dir, all_events.name[ev], 'particle_motion')
    if not os.path.isdir(pm_output):
        os.mkdir(pm_output)
    # ===== 1. Find all the station channels
    for indx, station in enumerate(all_stations.nam):
        print indx, station
        # XXX hard coded for now XXX
        nr_of_streams = all_stations.number_stations/3
        indx_cha = np.where(all_stations.nam == all_stations.nam[indx])
        #import ipdb; ipdb.set_trace()
        # XXX order of traces depends on input file ph_channel flag...         
        for ic in indx_cha[0]:
            if all_stations.cha[ic][-1:] == 'Z':
                trace_z = real_waveforms[ic] 
            if all_stations.cha[ic][-1:] == 'X':
                trace_x = real_waveforms[ic]
            if all_stations.cha[ic][-1:] == 'Y':
                trace_y = real_waveforms[ic]

        st_pm = Stream(traces=[trace_z, trace_x, trace_y])

        # ===== 2. read in the horizontal components and do instrument correction
        # ===== 3. plot and calculate partical motion 


        azi, inci, err_azi, err_inci = particle_motion_odr(st_pm)
        azi_fl, inci_fl, rectilinearity_fl, planarity_fl = flinn(st_pm)
        
        #plt.figure(figsize=(16, 20), constrained_layout=True)
        #gs = fig.add_gridspec(2, 2)

        #plt.subplot()
        #fig, (ax0, ax1, a2, a3) = plt.subplots(2, 2)
        # fig, axs = plt.subplots(2, 2)
        import ipdb; ipdb.set_trace()
        fig = plt.figure()
        # set up the axes for the first plot
        ax0 = fig.add_subplot(2, 2, 1, projection='3d')
 
        ax0.plot(trace_x.data, trace_y.data, trace_z.data, 'gray')
        
        plt.suptitle('BAZ cat: %s | BAZ calc %s' % (all_stations.bazi[indx], (360+azi) % 360))

        plt.savefig(os.path.join(pm_output, 'pm_%s.png' % station))
        
        if indx == nr_of_streams -1:
            return

# ------------------ NE_2_RT --------------------


def NE_2_RT(all_stations, real_waveforms, inp_ffproc,
            all_events, ev, plot_RT=True):
    """
    Function requests T and R component and sorts out real_waveforms and
    all_station data and information to make it consistent
    with the rest of the dataset
    :param all_stations:
    :param real_waveforms:
    :param inp_ffproc:
    :param station_pairs:
    :param ev:
    :param plot_RT:
    :return:
    """
    cprint('rotation_helper.py', '[S PHASE MODE]', bc.blue,
           'Will change from NE to RT system '
           'and continue with %s component for measurements!' % inp_ffproc.ph_horizontal)

    # filter_crit: flags those stations that do not have 3 components
    # Note that, we need 3 components for the rotation
    TR_comp, filter_crit = getTR(all_stations, real_waveforms, inp_ffproc,
                         all_events, ev, plot_RT=True)

    for indx, t_station in enumerate(TR_comp):
        # here the new channel name takes care of the difference between eg BHE and HHE
        new_channel = '%sH%s' % (all_stations.name[indx].split('.')[-1][0], inp_ffproc.ph_horizontal)
        real_waveforms[indx].stats['channel'] = new_channel
        if not filter_crit[indx]:
            real_waveforms[indx].stats['channel'] = 'XXX'
        real_waveforms[indx].data = t_station
        all_stations.name[indx] = '%s.%s.%s.%s' \
                               % (real_waveforms[indx].stats.network,
                                  real_waveforms[indx].stats.station,
                                  real_waveforms[indx].stats.location,
                                  new_channel)
    # remove the unwanted traces
    all_stations.rmsta(filter_crit)
    if not len(real_waveforms) < 1:
        for tr in real_waveforms.select(component="X"):
            real_waveforms.remove(tr)


# ------------------ get_inv --------------------

def get_inv(inp_ffproc, events_name, station_name):
    '''
    :param inp_ffproc:
    :param station_name:
    :return: azimuth, dip
    '''
    # XXX instead of hard-coding "resp" here, maybe we can have an option in the input file
    # XXX STXML is fine because we used that also in the previous version of DMT
    path_xml = glob.glob(os.path.join(inp_ffproc.real_path, events_name, 'resp', 'STXML.' + station_name))[0]

    if os.path.exists(path_xml):
        inventory = read_inventory(path_xml)
        return inventory[0][0][0].azimuth, inventory[0][0][0].dip

    else:
        try:
            path_less = glob.glob(os.path.join(inp_ffproc.real_path, events_name, 'resp', 'DATALESS.' + station_name))
            inventory = Parser(path_less)
            azimuth = inventory.stations[0][1].azimuth
            dip = inventory.stations[0][1].dip
            return azimuth, dip
        except:
            cprint('rotation_helper.py', '[get_inv] [ERROR]', bc.red, 'NO response file exists or cannot read file. '
                                                                      'NEEDS checking! EXIT')
            # theoretically this should never happen, because we would not have a processed trace in the first trace
            # to check for the response file
            sys.exit(2)

# ------------------ plot_rot --------------------


def plot_rot(R, T, Zcomp, tr_N, tr_E, inp_ffproc, all_events, ev, tr_1, tr_2):
    '''

    :param R: radial
    :param T: transversal
    :param tr_N:
    :param tr_E:
    :param inp_ffproc:
    :param ev: event_name
    :param tr_1: optional
    :param tr_2: optional
    :return:
    '''
    # import ipdb; ipdb.set_trace()
    station_name = '%s.%s.%s' \
                   % (Zcomp.stats['network'],
                      Zcomp.stats['station'],
                      Zcomp.stats['location'])

    time_ax = Zcomp.times()

    output_dir = os.path.join(inp_ffproc.output_dir, all_events.name[ev],
                              'rotate_NE_2_RT')

    file_name = '%s_rotate_NE_RT.png' % station_name

    if not os.path.isdir(output_dir):
        os.makedirs(os.path.join(output_dir))

    plt.figure()
    plt.ioff()

    plt.subplot(4, 1, 1)
    if len(tr_1) > 1:
        plt.plot(time_ax,
                 tr_1 / max(abs(T)),
                 color=(204 / 250., 37 / 250., 41 / 250.),
                 lw=3, alpha=0.5, label='1')

    plt.plot(time_ax,
             tr_N / max(abs(T)),
             color=(204/250., 37/250., 41/250.),
             lw=2, label='N')

    plt.legend(fontsize=8)
    plt.ylim(-1, 1)

    plt.subplot(4, 1, 2)
    if len(tr_2) > 1:
        plt.plot(time_ax,
                 tr_2 / max(abs(T)),
                 color=(57/250., 106/250., 177/250.),
                 lw=3, alpha=0.5, label='2')

    plt.plot(time_ax,
             tr_E / max(abs(T)),
             color=(57/250., 106/250., 177/250.),
             lw=2, alpha=1, label='E')


    plt.legend(fontsize=8)
    plt.ylim(-1, 1)

    plt.subplot(4, 1, 3)
    plt.plot(time_ax, R / max(abs(T)),
             color=(218/250., 124/250., 48/250.),
             lw=2, label='R')

    plt.legend(fontsize=8)
    plt.ylim(-1, 1)

    plt.subplot(4, 1, 4)
    plt.plot(time_ax, T / max(abs(T)),
             color=(63/250., 150/250., 81/250.),
             lw=2, label='T')

    plt.legend(fontsize=8)
    plt.xlabel('[sec]', weight='bold')
    plt.ylim(-1, 1)

    plt.suptitle('%s: Rotation NE -> RT (normalized for T)' % station_name,
                 weight='bold', fontsize=15)
    plt.savefig(os.path.join(output_dir, file_name), format='png')
    plt.clf()
    plt.close()

# ------------------ getTR --------------------


def getTR(all_stations, real_waveforms, inp_ffproc, all_events, ev, plot_RT):
    """
    get T-component or R-component depending on user input
    :param all_stations:
    :param real_waveforms:
    :param inp_ffproc:
    :param ev:
    :param plot_RT:
    :return:
    """
    T_comp = []
    R_comp = []
    filter_crit = np.ones(all_stations.number_stations, dtype=bool)
    path_2_S = os.path.join(inp_ffproc.real_path, all_events.name[ev], inp_ffproc.real_name_format)

    for indx, station in enumerate(all_stations.name):
        staName_split = station.split('.')
        print 'Working on: %s.%s.%s' %(staName_split[0], staName_split[1], staName_split[2])

        # =====>>>>> see if we have 1 and 2 components
        sta_nam_1 = station[:-1] + '1'
        sta_nam_2 = station[:-1] + '2'
        path_tr_1 = os.path.join(path_2_S, sta_nam_1)
        path_tr_2 = os.path.join(path_2_S, sta_nam_2)

        if not (os.path.exists(path_tr_1) or os.path.exists(path_tr_2)):
            pass
        else:
            try:
                tr_1 = read_horz_comp(path_tr_1, all_events, ev, all_stations, indx, inp_ffproc)
                tr_2 = read_horz_comp(path_tr_2, all_events, ev, all_stations, indx, inp_ffproc)
                # extract azimuth and dip information from xml/dataless files
                az_0, dip_0 = get_inv(inp_ffproc, all_events.name[ev], station)
                az_1, dip_1 = get_inv(inp_ffproc, all_events.name[ev], sta_nam_1)
                az_2, dip_2 = get_inv(inp_ffproc, all_events.name[ev], sta_nam_2)
                if not len(tr_1.data) == len(tr_2.data):
                    cprint('rotation_helper.py',
                           '[get_t_comp] [ERROR]', bc.red,
                           '1 and 2 components should have the same length! \nEXIT!')
                    sys.exit(2)
                # rotate the 1 and 2 component to N and E
                # zcomp need azimuth and dip values from the metadata
                Zcomp, Ncomp, Ecomp = rotate2zne(real_waveforms[indx], az_0, dip_0,
                                                 tr_1, az_1, dip_1,
                                                 tr_2, az_2, dip_2)
                R, T = rotate_ne_rt(Ncomp,
                                    Ecomp,
                                    all_stations.bazi[indx])
                T_comp.append(T)
                R_comp.append(R)
                filter_crit[indx] = True

                if plot_RT:
                    plot_rot(R, T, real_waveforms[indx], Ncomp, Ecomp, inp_ffproc, all_events, ev, tr_1.data, tr_2.data)
                continue

            except Exception, exp:
                filter_crit[indx] = False
                R_comp.append(np.zeros(len(real_waveforms[indx])))
                T_comp.append(np.zeros(len(real_waveforms[indx])))
                cprint('rotation_helper.py', '[INFORMATION]', bc.yellow, 'Skipping this station (%s)! Continue.' % exp)
                continue

        # =====>>>>> see if we have N and E components
        sta_nam_n = station[:-1] + 'N'
        sta_nam_e = station[:-1] + 'E'
        path_tr_n = os.path.join(path_2_S, sta_nam_n)
        path_tr_e = os.path.join(path_2_S, sta_nam_e)

        if not (os.path.exists(path_tr_n) or os.path.exists(path_tr_e)):
            filter_crit[indx] = False
            R_comp.append(np.zeros(len(real_waveforms[indx])))
            T_comp.append(np.zeros(len(real_waveforms[indx])))
            cprint('rotation_helper.py', '[INFORMATION]', bc.orange,
                   '%s does not have all horizontal components! Continue.' % station)
            continue
        else:
            try:
                tr_n = read_horz_comp(path_tr_n, all_events, ev, all_stations, indx, inp_ffproc)
                tr_e = read_horz_comp(path_tr_e, all_events, ev, all_stations, indx, inp_ffproc)

                if not len(tr_e.data) == len(tr_n.data):
                    cprint('rotation_helper.py',
                           '[get_t_comp] [ERROR]', bc.red,
                           'N and E components should have the same length! \nEXIT!')
                    sys.exit(2)

                R, T = rotate_ne_rt(tr_n.data,
                                    tr_e.data,
                                    all_stations.bazi[indx])
                T_comp.append(T)
                R_comp.append(R)
                filter_crit[indx] = True
                if plot_RT:
                    plot_rot(R, T, real_waveforms[indx], tr_n.data, tr_e.data, inp_ffproc, all_events, ev, [0], [0])

            except Exception, exp:
                filter_crit[indx] = False
                R_comp.append(np.zeros(len(real_waveforms[indx])))
                T_comp.append(np.zeros(len(real_waveforms[indx])))
                cprint('rotation_helper.py', '[INFORMATION]', bc.yellow, 'Skipping this station (%s)! Continue.' % exp)
                continue

    if inp_ffproc.ph_horizontal == 'T':
        return T_comp, filter_crit
    elif inp_ffproc.ph_horizontal == 'R':
        return R_comp, filter_crit
    else:
        cprint('rotation_helper.py', '[ERROR]', bc.dred, 'This component is not implemted or does not exist for '
                                                         'the horizontal measruements. EXIT!')
        sys.exit(33)
# ------------------ read_horz_comp --------------------


def read_horz_comp(path_tr, all_events, ev, all_stations, indx, inp_ffproc, lancz_a=12):
    """
    modified accordingly from StationClass/ read_real_cat
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
        cprint('rotation_helper.py', '[ERROR] [FILTER]', bc.dred,
               '%s filter mode has not implemented' % inp_ffproc.filt_mode)
        sys.exit()
    # accurate duration of the seismogram
    acc_time = inp_ffproc.ph_preset + max_filt_duration + inp_ffproc.ph_offset
    # acc_dur: number of samples
    acc_dur = (inp_ffproc.ph_sampling_rate * acc_time) + 1
    acc_dur = int(acc_dur)

    cha_id = path_tr.split('.')
    real_name = '%s.%s.%s.%s' % (cha_id[0],
                                 cha_id[1],
                                 cha_id[2],
                                 cha_id[3])
    try:
        start_time = all_events.datetime[ev] + all_stations.tt_ph[indx]
        len_tr = max_filt_duration + inp_ffproc.ph_offset + 10.
        try:
            tr = read(path_tr,
                      starttime=start_time-inp_ffproc.ph_preset,
                      endtime=start_time+len_tr,
                      format='SAC')[0]
        except Exception, e:
            tr = read(path_tr,
                      starttime=start_time-inp_ffproc.ph_preset,
                      endtime=start_time+len_tr)[0]

        starttime_diff = tr.stats.starttime - \
                         (start_time - inp_ffproc.ph_preset)
        if not starttime_diff <= tr.stats.delta:
            cprint('rotation_helper.py', '[WARNING] [START-TIME]', bc.orange,
                   'starttime of the trace differs from the requested '
                   'starttime.')
            raise Exception

        if not suppress_sample_warn:
            if tr.stats.sampling_rate > inp_ffproc.ph_sampling_rate:
                cprint('rotation_helper.py', '[WARNING] [SAMPLING]',
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
            cprint('rotation_helper.py', '[ERROR] [DATA]', bc.dred,
                   "Length of the data: %s\nRequested length: %s"
                   % (len(tr.data), acc_dur))
            raise Exception
        if len(np.where(np.isnan(tr.data))[0]) > 0:
            raise Exception
        return tr
    except Exception, e:
        cprint('rotation_helper.py', '[WARNING] [REAL]', bc.orange,
               '%s' % (e))