#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  EventClass.py
#   Purpose:   Event class for ffprocessing
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
import glob
import numpy as np
from obspy import UTCDateTime
import os
import pickle
import sys
from src.utility_codes import bc, cprint
from src.utility_codes import geocen_array

# necessary for ami cmt formatted files
import fortranformat as ff
# ------------------ read_events ---------------------------


def read_events(inp_ffproc):
    """
    reading all the events, this function decides which method of reading
    should be used.
    :param inp_ffproc:
    :return:
    """
    cprint('EventClass.py', '[EVENT]', bc.lgreen,
           'Reading event mode: %s' % inp_ffproc.evproc_mode)
    if inp_ffproc.evproc_mode == 'read':
        # read_events_cat is based on our original data archive
        eventcl = read_events_cat(inp_ffproc=inp_ffproc)
    elif inp_ffproc.evproc_mode == 'read_dmt':
        # read_events_dmt is based on DMT archive
        eventcl = read_events_dmt(inp_ffproc=inp_ffproc)
    elif inp_ffproc.evproc_mode == 'read_ami':
        # read_events_ami is based on AMI folder
        eventcl = read_events_ami(inp_ffproc=inp_ffproc)
    elif inp_ffproc.evproc_mode == 'read_obs':
        # read_events_obs is reading DMT EVENTS-INFO folder
        # because no event information is stored in the mseed folders of the OBS data
        eventcl = read_events_obs(inp_ffproc=inp_ffproc)
    else:
        cprint('EventClass.py', '[ERROR] [READ EVENTS]', bc.dred,
               '%s has not implemented!' % inp_ffproc.evproc_mode)
        sys.exit()

    return eventcl

# ------------------ read_events_cat ---------------------------


def read_events_cat(inp_ffproc):
    """
    read ampinv.source files and create EventClass
    :param inp_ffproc:
    :return:
    """
    event_cl = EventClass()
    all_events_rd = glob.glob(os.path.join(inp_ffproc.evproc_path,
                                           inp_ffproc.evproc_name_format))
    all_events_rd.sort()
    # import ipdb; ipdb.set_trace()
    for event_rd in all_events_rd:
        try:
            ampinv_source_fio = open(os.path.join(event_rd, 'outfiles',
                                                  'ampinv.source'), 'r')
        except Exception, exc:
            cprint('EventClass.py', '[WARNING] [AMPINV.SOURCE]', bc.orange,
                   'ampinv.source can not be open in %s\n%s.'
                   '\n====>>>> Going to use obspyDMT event.pkl information '
                   'for this event.\n' % (event_rd, exc))
            # MT: snippet of read_events_dmt function
            try:
                # Pickle file is generated for each event separately and should
                # contain all the information that we need at this step.
                ev = os.path.basename(event_rd)
                ev_load = open(os.path.join(inp_ffproc.real_path, ev, 'info', 'event.pkl'), 'r')
                ev_pkl = pickle.load(ev_load)
                ev_load.close()
            except Exception, exc:
                cprint('EventClass.py', '[PICKLE]', bc.lred,
                       'event.pkl can not be open in %s - %s. CONTINUE!\n'
                       % (os.path.join(inp_ffproc.real_path, ev, 'info', 'event.pkl'), exc))
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
            try:
                ev_half_duration = float(ev_pkl['source_duration'][1]) / 2.
            except Exception, error:
                cprint('EventClass.py', '[WARNING][SOURCE DURATION]', bc.lred,
                       'is not set! Half-duration is set to 2.5sec!')
                # ev_half_duration = 2.5
                ev_half_duration = mag_duration(moment_magnitude)

            event_cl.addev(ev_name, ev_date.year, ev_date.julday, ev_date.hour,
                           ev_date.minute, ev_date.second, ev_date.microsecond,
                           ev_pkl['latitude'], ev_pkl['longitude'],
                           ev_pkl['depth'], ev_pkl['depth'],
                           scalar_moment_magnitude, ev_half_duration,
                           mrr, mtt, mpp, mrt, mrp, mtp,
                           ev_pkl['magnitude'], 'DMT')
            continue
        # MT: back to original function
        f_source = ampinv_source_fio.readlines()
        ev_year, ev_julianday, ev_hr, ev_min, ev_sec, ev_msec = \
            f_source[1].split()
        evlat, evlon, catalog_depth, inverted_depth = f_source[3].split()
        # tausrc: adjusted rupture duration
        scalar_moment, tausrc = f_source[5].split()
        # if the source file does not contain the inverted source info,
        # the data will be read based on NEIC, HARVARD catalogs
        # However, this is different from reading the information directly
        # from NEIC and GCMT
        try:
            mrr, mtt, mpp, mrt, mrp, mtp = f_source[13].split()
            event_catalog = 'inverted'
        except Exception, exc:
            cprint('EventClass.py', '[INVERTED]', bc.orange,
                   'Inverted information does not exist: %s -- %s ---> '
                   'Read the information from catalog!' % (exc, event_rd))
            mrr, mtt, mpp, mrt, mrp, mtp = f_source[7].split()
            event_catalog = 'catalog'
        try:
            fio_mag = open(os.path.join(event_rd, 'README'), 'r')
            f_mag = fio_mag.readlines()
            ev_mag = float(f_mag[1].split()[1])
        except Exception, exc:
            cprint('EventClass.py', '[MAGNITUDE]', bc.orange,
                   'Magnitude can not be read: %s -- %s\n'
                   'Magnitude will be calculated!' % (exc, event_rd))
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

        ev_name = os.path.basename(event_rd)
        event_cl.addev(ev_name, ev_year, ev_julianday, ev_hr, ev_min, ev_sec,
                       int(ev_msec)*1000, evlat, evlon, catalog_depth,
                       inverted_depth, scalar_moment, float(tausrc)/2.,
                       mrr, mtt, mpp, mrt, mrp, mtp, ev_mag,
                       event_catalog)
    return event_cl

# ------------------ read_events_dmt ---------------------------


def read_events_dmt(inp_ffproc):
    """
    create EventClass from archives generated by DMT
    :param inp_ffproc:
    :return:
    """
    event_cl = EventClass()
    all_events_rd = glob.glob(os.path.join(inp_ffproc.evproc_path,
                                           inp_ffproc.evproc_name_format))
    all_events_rd.sort()

    # import ipdb; ipdb.set_trace()

    for event_rd in all_events_rd:
        try:
            # Pickle file is generated for each event separately and should
            # contain all the information that we need at this step.
            ev_load = open(os.path.join(event_rd, 'info', 'event.pkl'), 'r')
            ev_pkl = pickle.load(ev_load)
            ev_load.close()
        except Exception, exc:
            cprint('EventClass.py', '[PICKLE]', bc.lred,
                   'event.pkl can not be open in %s - %s. CONTINUE!'
                   % (event_rd, exc))
            continue
        #import ipdb; ipdb.set_trace()
        ev_date = ev_pkl['datetime']

        try:
            mrr = ev_pkl['focal_mechanism'][0]
            mtt = ev_pkl['focal_mechanism'][1]
            mpp = ev_pkl['focal_mechanism'][2]
            mrt = ev_pkl['focal_mechanism'][3]
            mrp = ev_pkl['focal_mechanism'][4]
            mtp = ev_pkl['focal_mechanism'][5]
        except Exception, exc:
            # import ipdb; ipdb.set_trace()
            cprint('EventClass.py', '[PICKLE]', bc.lred,
                   'event.pkl in %s does not contain moment tensor information - %s. CONTINUE!'
                   % (event_rd, exc))
            event_cl = read_events_obs(inp_ffproc)
            return event_cl

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

        try:
            ev_half_duration = float(ev_pkl['source_duration'][1])/2.
        except Exception, error:
            cprint('EventClass.py', '[WARNING][SOURCE DURATION]', bc.lred,
                   'is not set! Half-duration is set to 2.5sec!')
            #ev_half_duration = 2.5
            ev_half_duration = mag_duration(moment_magnitude)

        event_cl.addev(ev_name, ev_date.year, ev_date.julday, ev_date.hour,
                       ev_date.minute, ev_date.second, ev_date.microsecond,
                       ev_pkl['latitude'], ev_pkl['longitude'],
                       ev_pkl['depth'], ev_pkl['depth'],
                       scalar_moment_magnitude, ev_half_duration,
                       mrr, mtt, mpp, mrt, mrp, mtp,
                       ev_pkl['magnitude'], 'DMT')
    return event_cl

# ------------------ read_events_obs ---------------------------


def read_events_obs(inp_ffproc):
    """
    create EventClass from archives generated by DMT
    :param inp_ffproc:
    :return:
    """
    event_cl = EventClass()

    ev_load = open(os.path.join(inp_ffproc.evproc_path, 'EVENTS-INFO', 'event_list_pickle'), 'r')
    ev_all_pkl = pickle.load(ev_load)
    ev_load.close()

    if len(ev_all_pkl) < 1:
        cprint('EventClass.py', '[read_events_obs]', bc.lred,
               'event_list_pickle does not contain any events! [EXIT]')
        sys.exit()


    for ev_pkl in ev_all_pkl:

        #import ipdb; ipdb.set_trace()
        ev_date = ev_pkl['datetime']
        try:
            mrr = ev_pkl['focal_mechanism'][0]
            mtt = ev_pkl['focal_mechanism'][1]
            mpp = ev_pkl['focal_mechanism'][2]
            mrt = ev_pkl['focal_mechanism'][3]
            mrp = ev_pkl['focal_mechanism'][4]
            mtp = ev_pkl['focal_mechanism'][5]
        except Exception as exc:
            cprint('EventClass.py', '[read_events_obs]', bc.lred,
                   'event_list_pickle for %s des not contain moment tensor information - %s. CONTINUE!'
                   % (ev_pkl['datetime'], exc))
            continue

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

        try:
            ev_half_duration = float(ev_pkl['source_duration'][1])/2.

        except Exception as exc:
            # cprint('EventClass.py', '[WARNING][SOURCE DURATION]', bc.lred,
            #        'is not set! Half-duration is set to 2.5sec!')
            # ev_half_duration = 2.5
            # XXX maybe this needs a proper documentation in the event_cl? XXX
            ev_half_duration = mag_duration(moment_magnitude)

        event_cl.addev(ev_name, ev_date.year, ev_date.julday, ev_date.hour,
                       ev_date.minute, ev_date.second, ev_date.microsecond,
                       ev_pkl['latitude'], ev_pkl['longitude'],
                       ev_pkl['depth'], ev_pkl['depth'],
                       scalar_moment_magnitude, ev_half_duration,
                       mrr, mtt, mpp, mrt, mrp, mtp,
                       ev_pkl['magnitude'], 'DMT')

    return event_cl


# ------------------ read_events_ami ---------------------------


def read_events_ami(inp_ffproc):
    """
    reads the event information from the AMI structure in DIAS
    :param inp_ffproc:
    :return:
    """
    event_cl = EventClass()

    # since we are reading a global dataset starting 1979 the amount of paths is too high and takes a long time with glob.glob
    all_events_rd = []
    for year in inp_ffproc.evproc_year_list:
        all_events_rd.extend(glob.glob(os.path.join(inp_ffproc.evproc_path,
                                                    year, '*',
                                                    inp_ffproc.evproc_name_format)))

    all_events_rd.sort()

    for event_rd in all_events_rd:
        try:
            # cmt file is generated for each event separately and should contain the information that we need at this step.
            ev_info_file = glob.glob(os.path.join(event_rd, 'cmt*'))[0]

            # mt.harvard.bin or individual files in DIAS have this specific form
            # 0 - year, 1 - jday, 2- hour, 3- min, 4 - sec, 5 - lat, 6 - lon, 7 - dep, 8 - M0, 9 - tau, 10 - scale,
            # 11 - Mrr, 12 - Mtt, 13 - Mpp, 14 - Mrt, 15 - Mrp, 16 - Mtp,
            # 17 - strike1, 18 - dip1, 19 - rake1, 20 - strike2, 21 - dip2, 22 - rake2
            #       M0 => scalar moment tensor
            #       tau => source duration
            #       scale => scaling factor to estimate tau (see note 2 ind explanation of the ndk format)
            # format(i4,i4,2i3,f5.1,f6.2,f7.2,f6.1,e11.2,f5.1,e11.2,6f6.3,1x,i3, 1x, i2, i4,1x, i3, 1x, i2, i4)
            # XXX more efficient way to read this file? np.loadtxt/np.genfromtxt gives error message with dtype

            reader = ff.FortranRecordReader('(i4,i4,2i3,f5.1,f6.2,f7.2,f6.1,e11.2,f5.1,e11.2,6f6.3,1x,i3, 1x, i2, i4,1x, i3, 1x, i2, i4)')
            file = open(ev_info_file, 'r')
            cmt_line = file.readline().replace('\n', '')
            cmt_list = reader.read(cmt_line)

            if len(cmt_list) < 23:
                cprint('EventClass.py', '[ERROR AMI CMT FORMAT]', bc.lred,
                   'Problem reading the cmt file. IPDB initiated.')
                import ipdb; ipdb.set_trace()

        except Exception, exc:
            cprint('EventClass.py', '[AMI CMT FORMAT]', bc.lred,
                   'cmt%s can not be open in %s - %s. CONTINUE!'
                   % (os.path.basename(event_rd), event_rd, exc))
            continue

        sec_int = int(cmt_list[4])
        sec_frac = int((cmt_list[4] - sec_int)*1000000)

        datetime = UTCDateTime(year=cmt_list[0], julday=cmt_list[1], hour=cmt_list[2], minute=cmt_list[3], second=sec_int, microsecond=sec_frac)

        ev_date = '%04i/%02i/%02i-%02i:%02i:%02i' \
                          % (datetime.year,
                             datetime.month,
                             datetime.day,
                             datetime.hour,
                             datetime.minute,
                             datetime.second)
        #import ipdb; ipdb.set_trace()
        try:
            mrr = cmt_list[11] * cmt_list[10]
            mtt = cmt_list[12] * cmt_list[10]
            mpp = cmt_list[13] * cmt_list[10]
            mrt = cmt_list[14] * cmt_list[10]
            mrp = cmt_list[15] * cmt_list[10]
            mtp = cmt_list[16] * cmt_list[10]
            
        except Exception, exc:
            cprint('EventClass.py', '[AMI CMT FORMAT]', bc.lred,
                   'cmt%s can not read moment tensor information in %s - %s. CONTINUE!'
                   % (os.path.basename(event_rd), event_rd, exc))
            continue

        ev_name = os.path.basename(event_rd)

        scalar_moment_magnitude = \
            np.sqrt(
                np.power(float(mrr), 2) +
                np.power(float(mtt), 2) +
                np.power(float(mpp), 2) +
                np.power(float(mrt), 2) +
                np.power(float(mrp), 2) +
                np.power(float(mtp), 2)) / np.sqrt(2)
        moment_magnitude = 2. / 3. * (np.log10(scalar_moment_magnitude) - 9.1)
        #import ipdb; ipdb.set_trace()
        try:
            ev_half_duration = float(cmt_list[9])/2.
        except Exception, error:
            cprint('EventClass.py', '[WARNING][SOURCE DURATION]', bc.lred,
                   'is not set! Half-duration is set to 2.5sec!')
            #ev_half_duration = 2.5
            ev_half_duration = mag_duration(moment_magnitude)

        event_cl.addev(ev_name, datetime.year, datetime.julday, datetime.hour,
                       datetime.minute, datetime.second, datetime.microsecond,
                       cmt_list[5], cmt_list[6],
                       cmt_list[7], cmt_list[7],
                       scalar_moment_magnitude, ev_half_duration,
                       mrr, mtt, mpp, mrt, mrp, mtp,
                       moment_magnitude, 'AMI')
    #print 'Done with', ev_name
    #import ipdb; ipdb.set_trace()
    return event_cl


# ------------------ filter_events ---------------------------


def filter_events(inp_ffproc, all_events):
    """
    filter the events based on specified criteria in input file
    :param inp_ffproc:
    :param all_events:
    :return:
    """
    #import ipdb; ipdb.set_trace()
    cprint('EventClass.py', '[EVENT]', bc.lmagenta,
           'Filtering %s events' % all_events.number_events)
    min_lon, max_lon, min_lat, max_lat = inp_ffproc.evproc_rect.split('/')

    indxs = all_events.mag >= inp_ffproc.evproc_min_mag
    indxs = np.multiply(indxs, all_events.mag <= inp_ffproc.evproc_max_mag)
    indxs = np.multiply(indxs,
                        all_events.inv_dp >= inp_ffproc.evproc_min_depth)
    indxs = np.multiply(indxs,
                        all_events.inv_dp <= inp_ffproc.evproc_max_depth)
    indxs = np.multiply(indxs, all_events.year >= inp_ffproc.evproc_min_year)
    indxs = np.multiply(indxs, all_events.year <= inp_ffproc.evproc_max_year)
    indxs = np.multiply(indxs, all_events.lon >= float(min_lon))
    indxs = np.multiply(indxs, all_events.lon <= float(max_lon))
    indxs = np.multiply(indxs, all_events.lat >= float(min_lat))
    indxs = np.multiply(indxs, all_events.lat <= float(max_lat))
    all_events.rmev(indxs)
    cprint('EventClass.py', '[EVENT]', bc.dblue,
           'Number of events: %s (after filtering)' % all_events.number_events)

# ------------------ mag_duration ---------------------------


def mag_duration(mag, type_curve=1):
    """
    calculate the source duration out of magnitude
    type_curve can be 1, 2, 3:
    1: 2005-2014
    2: 1976-1990
    3: 1976-2014
    :param mag:
    :param type_curve:
    :return:
    """
    if type_curve == 1:
        half_duration = 0.00272*np.exp(1.134*mag)
    elif type_curve == 2:
        half_duration = 0.00804*np.exp(1.025*mag)
    elif type_curve == 3:
        half_duration = 0.00392*np.exp(1.101*mag)
    else:
        sys.exit('%s Type for magnitude to source duration conversion is not '
                 'implemented' % type_curve)
    return round(half_duration, 3)

# ------------------ EventClass ---------------------------


class EventClass:
    def __init__(self):
        self.number_events = 0
        self.name = np.array([])
        self.year = np.array([])
        self.julday = np.array([])
        self.hr = np.array([])
        self.minu = np.array([])
        self.sec = np.array([])
        self.msec = np.array([])
        self.datetime = np.array([])
        self.lat = np.array([])
        self.geocen_lat = np.array([])
        self.lon = np.array([])
        self.cat_dp = np.array([])
        self.inv_dp = np.array([])
        self.scmom = np.array([])
        self.tau = np.array([])
        self.mrr = np.array([])
        self.mtt = np.array([])
        self.mpp = np.array([])
        self.mrt = np.array([])
        self.mrp = np.array([])
        self.mtp = np.array([])
        self.mag = np.array([])
        self.catalog = np.array([])

    def addev(self, name, year, julday, hr, minu, sec, msec, lat, lon, cat_dp,
              inv_dp, scmom, tau, mrr, mtt, mpp, mrt, mrp, mtp, mag, catalog):
        """
        Add event to the EventClass
        :param name:
        :param year:
        :param julday:
        :param hr:
        :param minu:
        :param sec:
        :param msec:
        :param lat:
        :param lon:
        :param cat_dp:
        :param inv_dp:
        :param scmom:
        :param tau:
        :param mrr:
        :param mtt:
        :param mpp:
        :param mrt:
        :param mrp:
        :param mtp:
        :param mag:
        :param catalog:
        :return:
        """
        self.name = np.append(self.name, name)
        self.year = np.append(self.year, int(year)).astype(int)
        self.julday = np.append(self.julday, int(julday)).astype(int)
        self.hr = np.append(self.hr, int(hr)).astype(int)
        self.minu = np.append(self.minu, int(minu)).astype(int)
        self.sec = np.append(self.sec, int(sec)).astype(int)
        self.msec = np.append(self.msec, int(msec)).astype(int)
        # Avoid wrong msec reported in the processed data
        if int(msec) == 1e6:    
            sec = int(sec) + 1
            msec = 0
        ev_datetime = UTCDateTime(year=int(year),
                                  julday=int(julday),
                                  hour=int(hr),
                                  minute=int(minu),
                                  second=int(sec),
                                  microsecond=int(msec))
        self.datetime = np.append(self.datetime, ev_datetime)
        self.lat = np.append(self.lat, float(lat))
        self.geocen_lat = np.append(self.geocen_lat, geocen_array(np.array([float(lat)])))
        self.lon = np.append(self.lon, float(lon))
        self.cat_dp = np.append(self.cat_dp, abs(float(cat_dp)))
        self.inv_dp = np.append(self.inv_dp, abs(float(inv_dp)))
        self.scmom = np.append(self.scmom, float(scmom))
        self.tau = np.append(self.tau, float(tau))
        self.mrr = np.append(self.mrr, float(mrr))
        self.mtt = np.append(self.mtt, float(mtt))
        self.mpp = np.append(self.mpp, float(mpp))
        self.mrt = np.append(self.mrt, float(mrt))
        self.mrp = np.append(self.mrp, float(mrp))
        self.mtp = np.append(self.mtp, float(mtp))
        self.mag = np.append(self.mag, float(mag))
        self.catalog = np.append(self.catalog, catalog)
        self.number_events += 1

    def rmev(self, indxs):
        """
        Remove events from the EventClass ---> in case that the event(s)
        do(es) not meet the criteria
        :param indxs:
        :return:
        """
        self.name = self.name[indxs]
        self.year = self.year[indxs]
        self.julday = self.julday[indxs]
        self.hr = self.hr[indxs]
        self.minu = self.minu[indxs]
        self.sec = self.sec[indxs]
        self.msec = self.msec[indxs]
        self.datetime = self.datetime[indxs]
        self.lat = self.lat[indxs]
        self.geocen_lat = self.geocen_lat[indxs]
        self.lon = self.lon[indxs]
        self.cat_dp = self.cat_dp[indxs]
        self.inv_dp = self.inv_dp[indxs]
        self.scmom = self.scmom[indxs]
        self.tau = self.tau[indxs]
        self.mrr = self.mrr[indxs]
        self.mtt = self.mtt[indxs]
        self.mpp = self.mpp[indxs]
        self.mrt = self.mrt[indxs]
        self.mrp = self.mrp[indxs]
        self.mtp = self.mtp[indxs]
        self.mag = self.mag[indxs]
        self.catalog = self.catalog[indxs]
        self.number_events -= len(indxs[indxs == False])
