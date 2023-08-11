#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  STFutils.py
#   License:   GPLv3
#   Info:      most of the functions are recycled from existing codes,
#              like pyffproc or pyFFPostProc(refact)
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

import ConfigParser
import sys

from datetime import datetime
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

try:
    from obspy.imaging.beachball import beach as Beach
except Exception, e:
    from obspy.imaging.beachball import Beach
import os
from scipy.spatial import cKDTree
import sys
import socket
import math

import pickle
from obspy import UTCDateTime, read, Stream, Trace
from netCDF4 import Dataset
from obspy.clients.iris import Client
from obspy.signal.interpolation import lanczos_interpolation


from STFhelpers import *


# ------------------ read Event information from database


def read_events_dmt(input_file):
    """
    reads from the database downloaded with DMT the event information
    :param input_file:
    :return:
    """
    event_cl = EventClass()
    all_events_rd = glob.glob(os.path.join(input_file.path_database, input_file.database_format))
    all_events_rd.sort()

    for event_rd in all_events_rd:
        try:
            # Pickle file is generated for each event separately and should
            # contain all the information that we need at this step.
            ev_load = open(os.path.join(event_rd, 'info', 'event.pkl'), 'r')
            ev_pkl = pickle.load(ev_load)
            ev_load.close()
        except Exception, exc:
            cprint('STFutils.py', '[read_events_dmt]', bc.lred,
                   'event.pkl can not be open in %s - %s. CONTINUE!'
                   % (event_rd, exc))
            continue

        ev_date = ev_pkl['datetime']
        mrr = ev_pkl['focal_mechanism'][0]
        mtt = ev_pkl['focal_mechanism'][1]
        mpp = ev_pkl['focal_mechanism'][2]
        mrt = ev_pkl['focal_mechanism'][3]
        mrp = ev_pkl['focal_mechanism'][4]
        mtp = ev_pkl['focal_mechanism'][5]
        ev_name = ev_pkl['event_id']

        scalar_moment_magnitude = \
            np.sqrt(
                np.power(float(mrr), 2) +
                np.power(float(mtt), 2) +
                np.power(float(mpp), 2) +
                np.power(float(mrt), 2) +
                np.power(float(mrp), 2) +
                np.power(float(mtp), 2)) / np.sqrt(2)
        moment_magnitude = 2. / 3. * (np.log10(scalar_moment_magnitude) - 9.1)

        # calculate Flinn region
        client = Client()
        flinn_unformat = client.flinnengdahl(lat=ev_pkl['latitude'], lon=ev_pkl['longitude'],
                                    rtype="code")
        flinn = "{0:0=4d}".format(flinn_unformat)

        # get also month and day:
        date = UTCDateTime(year=ev_date.year, julday=ev_date.julday)
        month = date.month
        day = date.day

        qual = 2

        try:
            ev_half_duration = float(ev_pkl['source_duration'][1])/2.
        except Exception, error:
            cprint('STFutils.py', '[read_events_dmt] [WARNING SOURCE DURATION]', bc.lred,
                   'is not set! Half-duration is set to 2.5sec!')
            ev_half_duration = mag_duration(moment_magnitude)

        # create stf from dmt information
        # needed: number_grps, address_stf, grpaff_arr_all, stfs, stf_taxs
        number_grps, address_stf, grpaff_arr_all, stfs, stf_taxs = stf_read_dmt(input_file, ev_half_duration)

        event_cl.addev(ev_name, ev_date.year, ev_date.julday, month, day, ev_date.hour,
                       ev_date.minute, ev_date.second, ev_date.microsecond,
                       ev_pkl['latitude'], ev_pkl['longitude'],
                       ev_pkl['depth'], ev_pkl['depth'],
                       scalar_moment_magnitude, ev_half_duration,
                       mrr, mtt, mpp, mrt, mrp, mtp,
                       ev_pkl['magnitude'], qual,
                       flinn,
                       number_grps, address_stf, grpaff_arr_all, stfs, stf_taxs, 'DMT',
                       'DMT')

    stf_database = netCDF_writer(event_cl, input_file)
    # import ipdb;
    # ipdb.set_trace()
    return stf_database

# ------------------ netCDF_writer


def netCDF_writer(event_cl, input_file):
    '''

    :param event_cl:
    :param input_file:
    :return:
    '''

    # check/create output file folder and
    output_nc = os.path.join(input_file.path_out)
    if not os.path.isdir(output_nc):
        os.makedirs(output_nc)
    file_name = os.path.join(output_nc, '%s.nc' % input_file.nc_name)

    # create a stf's dataset (root group)

    stfs_ds = Dataset(file_name, "w", format="NETCDF4")

    # create groups for each band
    for indx, event in enumerate(event_cl.name):
        group_name = '%s.%02i.%02i-%02i:%02i:%02i-%04i' % (event_cl.year[indx], event_cl.month[indx],
                                               event_cl.day[indx], event_cl.hr[indx],
                                               event_cl.minu[indx], event_cl.sec[indx],
                                               event_cl.flinn[indx])

        # import ipdb; ipdb.set_trace()

        cprint('STFutils.py', '[netCDF_writer]', bc.grey, '...adding: %s' % group_name)

        # group name should be uniform throughout the netCDF
        event_group = stfs_ds.createGroup('%s' % group_name)

        # assign attributes for comparing later to other subgroups (SCARDEC, LOCAL)
        event_group.year = event_cl.year[indx]
        event_group.month = event_cl.month[indx]
        event_group.day = event_cl.day[indx]
        event_group.hour = event_cl.hr[indx]
        event_group.minute = event_cl.minu[indx]
        event_group.sec = event_cl.sec[indx]
        event_group.lon = event_cl.lon[indx]
        event_group.lat = event_cl.lat[indx]

        # event_subgroup = event_group.createGroup('%s' % event_cl.catalog[indx])
        create_nc_subgroup(event_group, event_cl.catalog[indx], event_cl, indx)

    # import ipdb; ipdb.set_trace()
    return output_nc


# ------------------ update netCDF


def update(stf_input, stf_database):
    '''

    :param stf_input:
    :param stf_database:
    :return:
    '''

    priority = ['database', 'scardec', 'localstf']
    priority.remove(stf_input.in_ref)
    for update in priority:
        if update == 'scardec':
            cprint('STFutils.py', '[update]', bc.green, 'Updating with SCARDEC database between the years: %s and %s'
                   % (stf_input.up_mindate, stf_input.up_maxdate))
            volcgrp = update_scardec(stf_input, stf_database)
        elif update == 'localstf':
            cprint('STFutils.py', '[update]', bc.green, 'Updating with local STFs database between the years: %s and %s'
                   % (stf_input.up_mindate, stf_input.up_maxdate))
            volcgrp = update_localstf(stf_input, stf_database)
        else:
            print update
            cprint('STFutils.py', '[update] [PROCESS ERROR]', bc.red,
                   'Not implemented yet. Work harder!')


# ------------------ update_localstf


def update_localstf(stf_input, stf_database):
    '''

    :param stf_input:
    :param stf_database:
    :return:
    '''

    cprint('STFutils.py', '[update_localstf]', bc.lgreen, ' ')

    # =====>>>>> read and write permission
    volcgrp = Dataset(stf_database, 'r+', format="NETCDF4")

    localstf_path = glob.glob(os.path.join(stf_input.path_localstf, stf_input.localstf_format))

    # =====>>>>> read source information from scardec
    localstf_cl = read_localstf(stf_input, localstf_path)

    # =====>>>>> compare to what we have in existing database
    cprint('STFutils.py', '[update_localstf]', bc.lgreen, ' ')

    volgroup = compare(localstf_cl, volcgrp, stf_input)
    # import ipdb;
    # ipdb.set_trace()

    return volcgrp


# ------------------ udpate_scardec


def update_scardec(stf_input, stf_database):
    '''

    :param stf_input:
    :param stf_database:
    :return:
    '''

    # =====>>>>> read and write permission
    volcgrp = Dataset(stf_database, 'r+', format="NETCDF4")

    scardec_path = glob.glob(os.path.join(stf_input.path_scardec, stf_input.scardec_format))

    # =====>>>>> read source information from scardec
    scardec_cl = read_scardec(stf_input, scardec_path)

    # =====>>>>> compare to what we have in existing database
    cprint('STFutils.py', '[update_scardec]', bc.lgreen, ' ')

    volcgrp = compare(scardec_cl, volcgrp, stf_input)

    return volcgrp


# ------------------ read_localstf


def read_localstf(stf_input, localstf_path):

    cprint('STFutils.py', '[read_localstf]', bc.lgreen, ' ')
    localstf_cl = EventClass()

    # import ipdb; ipdb.set_trace()

    for indx, path in enumerate(localstf_path):
        # import ipdb; ipdb.set_trace()
        stfs_glob = glob.glob(os.path.join(path, 'data', 'stf.0?'))
        stfs_glob.sort()
        
        # =====>>>>> read quality information
        try:
            postproc_qual = glob.glob(os.path.join(path, 'outfiles', 'postproc.qual'))
            qual = np.loadtxt(postproc_qual[0])
        except Exception, exp:
            cprint('STFutils.py', '[read_localstf] [NO postproc.qual]', bc.orange, '%s: Quality set to 1.'
                   % os.path.basename(path))
            qual = 1

        # =====>>>>> read source information
        try:
            ampinv_source_fio = open(os.path.join(path, 'outfiles',
                                     'ampinv.source'), 'r')
            f_source = ampinv_source_fio.readlines()
            ev_year, ev_julianday, ev_hr, ev_min, ev_sec, ev_msec = f_source[1].split()

            if stf_input.up_mindate <= eval(ev_year) <= stf_input.up_maxdate is False:
                print stf_input.up_mindate, '<=', ev_year, '<=', stf_input.up_maxdate
                continue

            # get also month and day:
            date = UTCDateTime(year=eval(ev_year), julday=eval(ev_julianday))
            ev_month = date.month
            ev_day = date.day

            ev_lat, ev_lon, catalog_depth, inverted_depth = f_source[3].split()
            scalar_moment, tau = f_source[5].split()

            # calculate Flinn region
            client = Client()
            flinn = client.flinnengdahl(lat=eval(ev_lat), lon=eval(ev_lon), rtype="code")
            flinn = "{0:0=4d}".format(flinn)
        except Exception, exp:
            cprint('STFutils.py', '[read_localstf] [NO ampinv.source]', bc.dred, '%s: No source information. '
                                                                                 'Skipping this event. Continue!'
                   % os.path.basename(path))
            continue

        try:
            mrr, mtt, mpp, mrt, mrp, mtp = f_source[13].split()
            event_catalog = 'inverted'
        except Exception, exc:
            cprint('STFutils.py', '[read_localstf]', bc.orange,
                   'Inverted information does not exist: %s\n'
                   'Read the information from catalog!' % exc)
            mrr, mtt, mpp, mrt, mrp, mtp = f_source[7].split()
            event_catalog = 'catalog'
        try:
            fio_mag = open(os.path.join(path, 'README'), 'r')
            f_mag = fio_mag.readlines()
            ev_mag = f_mag[1].split()[1]
            if not isinstance(eval(ev_mag), float):
                # print ev_mag
                raise Exception
        except Exception, exc:
            cprint('STFutils.py', '[read_localstf]', bc.orange,
                   'Magnitude can not be read: %s \n'
                   'Magnitude will be calculated!' % exc)
            scalar_moment_magnitude = \
                np.sqrt(
                    np.power(float(mrr), 2) +
                    np.power(float(mtt), 2) +
                    np.power(float(mpp), 2) +
                    np.power(float(mrt), 2) +
                    np.power(float(mrp), 2) +
                    np.power(float(mtp), 2)) / np.sqrt(2)
            moment_magnitude = \
                2. / 3. * (np.log10(scalar_moment_magnitude) - 9.1)
            scalar_moment = scalar_moment_magnitude
            ev_mag = moment_magnitude

        # =====>>>>> read group affiliation
        # header: idx_ai  grpaff_ai idx_ppc grpaff_ppc stla   stlo   dist  azi   stnam
        dt = np.dtype([('grpaff_idx_ai', int), ('grpaff_ai', int), ('grpaff_idx_ppc', int), ('grpaff_ppc', int),
                       ('grpaff_stla', float), ('grpaff_stlo', float), ('grpaff_dist', float), ('grpaff_azi', float),
                       ('grpaff_stnam', 'S20')])
        grpaff_arr_all = np.loadtxt(os.path.join(path, 'outfiles', 'ampinv.grpaff'),
                                    dtype=dt, skiprows=1)

        traces = []
        taxs = []
        number_grps = 0
        for stftmp in stfs_glob:
            # =====>>>>> default just one group
            number_grps += 1

            # =====>>>>> get stf
            stf = read(stftmp, format='SAC')[0]

            # =====>>>>> resampling STF
            starttime = stf.stats.starttime.timestamp
            endtime = stf.stats.endtime.timestamp
            dt = 1. / stf_input.ph_sampling_rate
            stf.data = lanczos_interpolation(
                np.require(stf.data, requirements=["C"]),
                starttime, stf.stats.delta, starttime, dt,
                int(math.floor((endtime - starttime) / dt)) + 1,
                a=stf_input.lancz_a, window="blackman")

            stf.stats.starttime = UTCDateTime(starttime)
            stf.stats.sampling_rate = stf_input.ph_sampling_rate

            # =====>>>>> get time axis
            stf_taxs = Trace(np.linspace(stf.stats.sac.b, stf.stats.sac.b + (stf.stats.npts - 1) / stf.stats.sampling_rate,
                                   stf.stats.npts))

            traces.append(stf)
            taxs.append(stf_taxs)

        ev_name = '%04i.%02i.%02i-%02i:%02i:%02i' \
                  % (eval(ev_year), ev_month, ev_day, eval(ev_hr), eval(ev_min), eval(ev_sec))

        localstf_cl.addev(ev_name, ev_year, ev_julianday, ev_month,
                      ev_day, ev_hr, ev_min, ev_sec,
                      ev_msec,
                      ev_lat, ev_lon, catalog_depth, inverted_depth,
                      scalar_moment, tau,
                      mrr, mtt, mpp, mrt, mrp, mtp,
                      ev_mag, float(qual),
                      flinn,
                      number_grps, path, grpaff_arr_all, traces, taxs, event_catalog,
                      'LOCAL')


    return localstf_cl


# ------------------ read_scardec


def read_scardec(stf_input, scardec_path):
    '''

    :param stf_input:
    :param scardec_path:
    :return:
    '''

    scardec_cl = EventClass()

    # how to read the files (from SCARDEC webpage)
    # 1st line : YYYY MM DD HH MM SS'.0' Latitude Longitude [origin time and epicentral location from NEIC]
    # 2nd line : Depth(km) M0(N.m) Mw strike1(°) dip1(°) rake1(°) strike2(°) dip2(°) rake2(°) [all from SCARDEC]
    # All the other lines are the temporal STF, with format : time(s), moment rate(N.m/s)

    for indx, path in enumerate(scardec_path):
        # np.loadtxt(os.path.join(scardec_path, '*moy*', )
        path_average = glob.glob(os.path.join(path, '*moy*'))[0]
        path_optimum = glob.glob(os.path.join(path, '*opt*'))[0]

        # headers are *always* identical, therefore I'm just reading one header for this stf
        with open(path_average) as myfile:
            # first two lines are the header I need the information
            head = [next(myfile).splitlines() for x in xrange(2)]

        line = head[0][0]
        year = int(line.split(' ')[0])

        if stf_input.up_mindate <= year <= stf_input.up_maxdate:
            line_1 = head[0][0]
            line_1 = line_1.split(' ')
            line_2 = head[1][0]
            line_2 = line_2.split(' ')

            # here I choose line 18 because all stf start time=0 at line 19
            # XXX change that to read from the negative values! therefore start from 3rd line!
            average_stf = np.loadtxt(path_average, skiprows=2)
            optimum_stf = np.loadtxt(path_optimum, skiprows=2)

            # stf to traces to be consistent with reading local stfs
            # 14 is sampling rate of the SCARDEC STF's
            tr_data_av = average_stf[:, 1]
            tr_eng_norm_av = tr_data_av/np.sum(tr_data_av)*14

            tr_data_opt = optimum_stf[:, 1]
            tr_eng_norm_opt = tr_data_opt/np.sum(tr_data_opt)*14

            tr_average = Trace(tr_eng_norm_av)
            tr_optimum = Trace(tr_eng_norm_opt)
            traces = [tr_average, tr_optimum]

            # x axis as list
            taxs_average = Trace(average_stf[:, 0])
            taxs_optimum = Trace(optimum_stf[:, 0])
            taxs = [taxs_average, taxs_optimum]

            ev_year = int(line_1[0])
            ev_month = int(line_1[1])
            ev_day = int(line_1[2])
            ev_hour = int(line_1[3])
            ev_minute = int(line_1[4])
            ev_sec = int(float(line_1[5]))

            ev_microsec = 0

            datetime = UTCDateTime(year=ev_year, month=ev_month, day=ev_day,
                                   hour=ev_hour, minute=ev_minute, second=ev_sec)

            ev_julday = "{0:0=4d}".format(datetime.julday)

            ev_name = '%04i.%02i.%02i-%02i:%02i:%02i' % (ev_year, ev_month, ev_day, ev_hour, ev_minute, ev_sec)

            ev_lat = float(line_1[6])
            ev_lon = float(line_1[7])

            ev_dp = float(line_2[0])
            ev_mom = float(line_2[1])
            ev_mag = float(line_2[2])

            ev_strike1 = int(line_2[3])
            ev_dip1 = int(line_2[4])
            ev_rake1 = int(line_2[5])
            moment_1 = sdr2mt(ev_strike1, ev_dip1, ev_rake1, ev_mom)
            # print '2 - %s ' % year

            ev_qual = 3.5

            # calculate Flinn region
            client = Client()
            try:
                flinn_unformat = client.flinnengdahl(lat=ev_lat, lon=ev_lon,
                                                     rtype="code")
                flinn = "{0:0=4d}".format(flinn_unformat)
            except Exception, exp:
                import ipdb; ipdb.set_trace()

            ev_halfduration = hd_calculator(average_stf[:, 0], average_stf[:, 1])
            # print '3 - %s ' % year

            # scardec stfs have ALWAYS two stfs (average and optimum, but no grouping)
            scardec_cl.addev(ev_name, ev_year, ev_julday, ev_month, ev_day, ev_hour,
                             ev_minute, ev_sec, ev_microsec, ev_lat, ev_lon,
                             ev_dp, ev_dp, ev_mom, ev_halfduration,
                             moment_1[0], moment_1[1], moment_1[2],
                             moment_1[3], moment_1[4], moment_1[5],
                             ev_mag, ev_qual, flinn,
                             2, path, False, traces, taxs,
                             'SCARDEC', 'SCARDEC')

        progress_bar(curr_indx=indx, total_number=len(scardec_path))

    print '\n'
    cprint('STFutils.py', '[read_scardec]', bc.green, 'Finished reading SCARDEC library!')
    return scardec_cl


# ------------------ plot

def plot_it(volcgrp):

    plt.figure()
    plt.ioff()

    # XXXXX Stopped here to plt mag ag half duration/tau
    # Q: how to access all the content at once without for loop? possible?exit
    mag_l = []
    tau_l = []

    mag_s = []
    tau_s = []

    mag_d = []
    tau_d = []

    cmap = plt.get_cmap('jet_r')
    colors = cmap(np.linspace(0, 1.0, 3))

    col_l = []
    col_s = []
    col_d = []

    for group in volcgrp.groups:
        print bc.green + group + bc.end
        for subgroup in volcgrp.groups[group].groups:
            print bc.blue + subgroup + bc.end
            if subgroup == 'LOCAL':
                col_l.append(colors[0])
                tau_l.append(volcgrp.groups[group].groups[subgroup].variables['tau'][0])
                mag_l.append(volcgrp.groups[group].groups[subgroup].variables['mag'][0])
            elif subgroup == 'SCARDEC':
                col_s.append(colors[1])
                tau_s.append(volcgrp.groups[group].groups[subgroup].variables['tau'][0])
                mag_s.append(volcgrp.groups[group].groups[subgroup].variables['mag'][0])
            else:
                col_d.append(colors[2])
                tau_d.append(volcgrp.groups[group].groups[subgroup].variables['tau'][0])
                mag_d.append(volcgrp.groups[group].groups[subgroup].variables['mag'][0])

    l = plt.scatter(tau_l, mag_l, c=col_l, lw=0, alpha=0.7)
    s = plt.scatter(tau_s, mag_s, c=col_s, lw=0, alpha=0.7)
    d = plt.scatter(tau_d, mag_d, c=col_d, lw=0, alpha=0.7)
    plt.legend((l, s, d), ('LOCAL', 'SCARDEC', 'DATABASE'), scatterpoints=1, loc='lower right', fontsize=10)
    plt.xlabel('Halfduration [sec]', weight='bold', size=15)
    plt.ylabel('Magnitude', weight='bold', size=15)
    plt.xticks(weight='bold')
    plt.yticks(weight='bold')
    plt.xlim(0, 100)
    plt.savefig('mag_vs_tau.png')