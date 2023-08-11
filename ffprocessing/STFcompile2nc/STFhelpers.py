#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  STFhelper.py
#   License:   GPLv3
#   Info:      helper functions for STFutils.py
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

import ConfigParser
from datetime import datetime
import numpy as np
import sys
import socket
from netCDF4 import Dataset
from obspy import UTCDateTime
from scipy.integrate import trapz, cumtrapz
from obspy.core import Stream, Trace

# ------------------ input_reader

class STFreadit:
    def __init__(self, inp_file):

        cprint('STFhelpers.py', '[STFreadit]', bc.bold + bc.yellow, 'Creating input class!')
        config_ini = ConfigParser.ConfigParser()
        config_ini.read(inp_file)

        # --------------- paths
        self.path_database = eval(config_ini.get('paths', 'path_database'))
        self.database_format = eval(config_ini.get('paths', 'database_format'))
        self.path_localstf = eval(config_ini.get('paths', 'path_localstf'))
        self.localstf_format = eval(config_ini.get('paths', 'localstf_format'))
        self.path_scardec = eval(config_ini.get('paths', 'path_scardec'))
        self.scardec_format = eval(config_ini.get('paths', 'scardec_format'))

        # --------------- threshold
        self.distance = float(config_ini.get('threshold', 'distance'))
        self.time = float(config_ini.get('threshold', 'time'))

        # --------------- initialise
        self.in_ref = eval(config_ini.get('initialise', 'in_ref'))
        self.in_mindate = eval(config_ini.get('initialise', 'in_mindate'))
        self.in_maxdate = eval(config_ini.get('initialise', 'in_maxdate'))

        # --------------- update
        self.up_ref = eval(config_ini.get('update', 'up_ref'))
        self.up_mindate = eval(config_ini.get('update', 'up_mindate'))
        self.up_maxdate = eval(config_ini.get('update', 'up_maxdate'))

        # --------------- best_choice
        self.qual_min = eval(config_ini.get('best_choice', 'qual_min'))

        # --------------- localstf
        self.ph_sampling_rate = eval(config_ini.get('localstf', 'ph_sampling_rate'))
        self.lancz_a = eval(config_ini.get('localstf', 'lancz_a'))

        # --------------- output
        self.path_out = eval(config_ini.get('output', 'path_out'))
        self.nc_name = eval(config_ini.get('output', 'nc_name'))


# ------------------ EventClass


class EventClass:
    def __init__(self):
        self.number_events = 0
        self.name = np.array([])
        self.year = np.array([])
        self.julday = np.array([])
        self.month = np.array([])
        self.day = np.array([])
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
        self.qual = np.array([])
        self.flinn = np.array([])
        self.number_grps = np.array([])
        self.address_stf = np.array([])
        self.grpaff_arr_all = []
        self.stfs = []
        self.stf_taxs = []
        self.event_catalog = np.array([])
        self.catalog = np.array([])

    def addev(self, name, year, julday, month, day, hr, minu,
              sec, msec, lat, lon, cat_dp,
              inv_dp, scmom, tau, mrr, mtt, mpp, mrt, mrp,
              mtp, mag, qual, flinn,
              number_grps, address_stf, grpaff_arr_all, stfs, stf_taxs, event_catalog,
              catalog):
        """
        Add event to the EventClass
        :param name:
        :param year:
        :param julday:
        :param month:
        :param day:
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
        :param qual:
        :param flinn:
        :param number_grps:
        :param address_stf:
        :param grpaff_arr_all:
        :param stfs:
        :param stf_taxs:
        :param catalog:
        :return:
        """
        self.name = np.append(self.name, str(name))
        self.year = np.append(self.year, int(year)).astype(int)
        self.julday = np.append(self.julday, int(julday)).astype(int)
        self.month = np.append(self.month, int(month)).astype(int)
        self.day = np.append(self.day, int(day)).astype(int)
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
        self.geocen_lat = np.append(self.geocen_lat,
                                    geocen_array(np.array([float(lat)])))

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
        self.qual = np.append(self.qual, float(qual))
        self.flinn = np.append(self.flinn, int(flinn))

        self.number_grps = np.append(self.number_grps, int(number_grps))
        self.address_stf = np.append(self.address_stf, str(address_stf))

        self.grpaff_arr_all.append(grpaff_arr_all)
        self.stfs.append(Stream(traces=stfs))
        self.stf_taxs.append(Stream(traces=stf_taxs))

        self.event_catalog = np.append(self.event_catalog, str(event_catalog))
        self.catalog = np.append(self.catalog, catalog)
        self.number_events += 1

    def rmev(self, indxs):
        """
        Remove events from the EventClass when
        do(es) not meet the criteria
        :param indxs:
        :return:
        """
        self.name = self.name[indxs]

        self.year = self.year[indxs]
        self.julday = self.julday[indxs]
        self.month = self.month[indxs]
        self.day = self.day[indxs]
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
        self.qual = self.qual[indxs]
        self.flinn = self.flinn[indxs]

        #
        self.number_grps = self.number_grps[indxs]
        self.address_stf = self.address_stf[indxs]
        # import ipdb; ipdb.set_trace()
        self.grpaff_arr_all = np.array(self.grpaff_arr_all)[indxs]
        self.stfs = np.array(self.stfs)[indxs]
        self.stf_taxs = np.array(self.stf_taxs)[indxs]

        self.event_catalog = self.event_catalog[indxs]
        self.catalog = self.catalog[indxs]
        self.number_events -= len(indxs[indxs == False])


# ------------------ create_nc_subgroup


def create_nc_subgroup(nc_group, name_subgroup, event_cl, indx):
    '''

    :param nc_group:
    :param name_subgroup:
    :param class2fill:
    :return:
    '''
    print name_subgroup
    event_subgroup = nc_group.createGroup('%s' % name_subgroup)
    # create dimension for input
    try:
        event_subgroup.createDimension('single_value', 1)
        event_subgroup.createDimension('unlimited', None)
    except Exception, exp:
        import ipdb; ipdb.set_trace()
    # event_out.createDimension(name, dimension.isunlimited())

    # create the variables

    name = event_subgroup.createVariable('name', 'S20', ('single_value',))

    year = event_subgroup.createVariable('year', 'i4', ('single_value',))
    julday = event_subgroup.createVariable('julday', 'i4', ('single_value',))
    month = event_subgroup.createVariable('month', 'i4', ('single_value',))
    day = event_subgroup.createVariable('day', 'i4', ('single_value',))
    hr = event_subgroup.createVariable('hr', 'f4', ('single_value',))
    minu = event_subgroup.createVariable('minu', 'f4', ('single_value',))
    sec = event_subgroup.createVariable('sec', 'f4', ('single_value',))
    msec = event_subgroup.createVariable('msec', 'f4', ('single_value',))

    lat = event_subgroup.createVariable('lat', 'f4', ('single_value',))
    geocen_lat = event_subgroup.createVariable('geocen_lat', 'f4', ('single_value',))
    lon = event_subgroup.createVariable('lon', 'f4', ('single_value',))
    cat_dp = event_subgroup.createVariable('cat_dp', 'f4', ('single_value',))
    inv_dp = event_subgroup.createVariable('inv_dp', 'f4', ('single_value',))

    scmom = event_subgroup.createVariable('scmom', 'f4', ('single_value',))
    tau = event_subgroup.createVariable('tau', 'f4', ('single_value',))
    mrr = event_subgroup.createVariable('mrr', 'f4', ('single_value',))
    mtt = event_subgroup.createVariable('mtt', 'f4', ('single_value',))
    mpp = event_subgroup.createVariable('mpp', 'f4', ('single_value',))
    mrt = event_subgroup.createVariable('mrt', 'f4', ('single_value',))
    mrp = event_subgroup.createVariable('mrp', 'f4', ('single_value',))
    mtp = event_subgroup.createVariable('mtp', 'f4', ('single_value',))
    mag = event_subgroup.createVariable('mag', 'f4', ('single_value',))

    qual = event_subgroup.createVariable('qual', 'f4', ('single_value',))
    flinn = event_subgroup.createVariable('flinn', 'f4', ('single_value',))

    number_grps = event_subgroup.createVariable('number_grps', 'f4', ('single_value',))
    address_stf = event_subgroup.createVariable('address_stf', 'S200', ('single_value',))

    if event_cl.grpaff_arr_all[indx] == 0:
        grpaff_arr_all = event_subgroup.createVariable('grpaff_arr_all', 'S20', ('single_value',))
        grpaff_arr_all[0] = str(event_cl.grpaff_arr_all[indx])

    else:
        grpaff_arr_all = event_subgroup.createVariable('grpaff_arr_all', 'S20', ('single_value',))
        grpaff_idx_ai = event_subgroup.createVariable('grpaff_idx_ai', 'i4', ('unlimited',))
        grpaff_ai = event_subgroup.createVariable('grpaff_ai', 'i4', ('unlimited',))
        grpaff_idx_ppc = event_subgroup.createVariable('grpaff_idx_ppc', 'i4', ('unlimited',))
        grpaff_ppc = event_subgroup.createVariable('grpaff_ppc', 'i4', ('unlimited',))
        grpaff_stla = event_subgroup.createVariable('grpaff_stla', 'f4', ('unlimited',))
        grpaff_stlo = event_subgroup.createVariable('grpaff_stlo', 'f4', ('unlimited',))
        grpaff_dist = event_subgroup.createVariable('grpaff_dist', 'f4', ('unlimited',))
        grpaff_azi = event_subgroup.createVariable('grpaff_azi', 'f4', ('unlimited',))
        grpaff_stnam = event_subgroup.createVariable('grpaff_stnam', 'S20', ('unlimited',))

        grpaff_arr_all[0] = str(event_cl.grpaff_arr_all[indx])
        grpaff_idx_ai[:] = event_cl.grpaff_arr_all[indx]['grpaff_idx_ai']
        grpaff_ai[:] = event_cl.grpaff_arr_all[indx]['grpaff_ai']
        grpaff_idx_ppc[:] = event_cl.grpaff_arr_all[indx]['grpaff_idx_ppc']
        grpaff_ppc[:] = event_cl.grpaff_arr_all[indx]['grpaff_ppc']
        grpaff_stla[:] = event_cl.grpaff_arr_all[indx]['grpaff_stla']
        grpaff_stlo[:] = event_cl.grpaff_arr_all[indx]['grpaff_stlo']
        grpaff_dist[:] = event_cl.grpaff_arr_all[indx]['grpaff_dist']
        grpaff_azi[:] = event_cl.grpaff_arr_all[indx]['grpaff_azi']
        grpaff_stnam[:] = event_cl.grpaff_arr_all[indx]['grpaff_stnam']

    if event_cl.number_grps[indx] == 1:
        stfs_1 = event_subgroup.createVariable('stfs_1', 'f4', ('unlimited',))
        stf_taxs_1 = event_subgroup.createVariable('stf_taxs_1', 'f4', ('unlimited',))

        stfs_1[:] = event_cl.stfs[indx][0].data
        stf_taxs_1[:] = event_cl.stf_taxs[indx][0].data

    elif event_cl.number_grps[indx] > 1:
        for stf_grp in range(int(event_cl.number_grps[indx])):
            grp_stf_name = 'stfs_%s' % (stf_grp + 1)
            grp_taxs_name = 'stf_taxs_%s' % (stf_grp + 1)

            grp_stf_name = event_subgroup.createVariable(grp_stf_name, 'f4', ('unlimited',))
            grp_taxs_name = event_subgroup.createVariable(grp_taxs_name, 'f4', ('unlimited',))

            grp_stf_name[:] = event_cl.stfs[indx][stf_grp].data
            grp_taxs_name[:] = event_cl.stf_taxs[indx][stf_grp].data

    else:
        cprint('STFhelpers.py', '[create_nc_subgroup]', bc.red + bc.bold,
               'Major problem. This should not happen: NO STF.')

    event_catalog = event_subgroup.createVariable('event_catalog', 'S20', ('single_value',))
    catalog = event_subgroup.createVariable('catalog', 'S20', ('single_value',))


    # assigning the values to the variables

    name[0] = str(event_cl.name[indx][:])
    try:
        year[:] = event_cl.year[indx]
    except Exception, exp:
        import ipdb; ipdb.set_trace()

    julday[:] = event_cl.julday[indx]
    month[:] = event_cl.month[indx]
    day[:] = event_cl.day[indx]
    hr[:] = event_cl.hr[indx]
    minu[:] = event_cl.minu[indx]
    sec[:] = event_cl.sec[indx]
    msec[:] = event_cl.msec[indx]

    lat[:] = event_cl.lat[indx]
    geocen_lat[:] = event_cl.geocen_lat[indx]
    lon[:] = event_cl.lon[indx]
    cat_dp[:] = event_cl.cat_dp[indx]
    inv_dp[:] = event_cl.inv_dp[indx]

    scmom[:] = event_cl.scmom[indx]
    tau[:] = event_cl.tau[indx]
    mrr[:] = event_cl.mrr[indx]
    mtt[:] = event_cl.mtt[indx]
    mpp[:] = event_cl.mpp[indx]
    mrt[:] = event_cl.mrt[indx]
    mrp[:] = event_cl.mrp[indx]
    mtp[:] = event_cl.mtp[indx]
    mag[:] = event_cl.mag[indx]

    qual[:] = event_cl.qual[indx]
    flinn[:] = event_cl.flinn[indx]

    number_grps[:] = event_cl.number_grps[indx]
    address_stf[0] = str(event_cl.address_stf[indx])

    event_catalog[0] = str(event_cl.event_catalog[indx][:])
    catalog[0] = str(event_cl.catalog[indx][:])


# ------------------ compare

def compare(event_cl, volcgrp, stf_input):
    # we might have two lists now of different lengths, all of them need to be checked...
    print event_cl.number_events
    # print event_cl.name
    for group in volcgrp.groups.values():
        year_volcgrp = group.year
        month_volcgrp = group.month
        day_volcgrp = group.day
        hour_volcgrp = group.hour
        minu_volcgrp = group.minute
        sec_volcgrp = group.sec

        lon_volcgrp = group.lon
        lat_volcgrp = group.lat

        # creating a filter list determining
        # XXX why not just use the date time for the filter, therefore the search would reduce to 3 lines
        # see the assimilate code for stf (works good)
        filter_list = event_cl.year == year_volcgrp
        filter_list *= event_cl.month == month_volcgrp
        filter_list *= event_cl.day == day_volcgrp

        true_list = np.where(filter_list == True)[0]

        # import ipdb; ipdb.set_trace()


        # finer search
        # i index v value
        try:
            i_lat, v_lat = min(enumerate(event_cl.lat[filter_list]),
                               key=lambda x: abs(x[1] - lat_volcgrp))
            i_lon, v_lon = min(enumerate(event_cl.lon[filter_list]),
                               key=lambda x: abs(x[1] - lon_volcgrp))
            i_hr, v_hr = min(enumerate(event_cl.hr[filter_list]),
                             key=lambda x: abs(x[1] - hour_volcgrp))
            i_minu, v_minu = min(enumerate(event_cl.minu[filter_list]),
                                 key=lambda x: abs(x[1] - minu_volcgrp))
            i_sec, v_sec = min(enumerate(event_cl.sec[filter_list]),
                               key=lambda x: abs(x[1] - sec_volcgrp))
        except Exception, exp:
            cprint('STFhelpers.py', '[EXCEPTION] [compare.py]', bc.orange, '%s - %s. Continue.' % (group.name, exp))
            continue

        # =====>>>>> if in the range of threshold add as group
        if abs(lat_volcgrp - v_lat) > stf_input.distance or abs(lon_volcgrp - v_lon) > stf_input.distance:
            continue
        if abs(sec_volcgrp - v_sec) > stf_input.time:
            continue
        if i_lat == i_lon == i_hr == i_minu == i_sec:

            # print group.name
            try:
                ev_indx = true_list[i_lat]
                cprint('STFhelpers.py', '[compare]', bc.grey, '...updating existing netDCF event %s with %s'
                       % (bc.blue + group.name + bc.end, bc.green + event_cl.name[ev_indx] + bc.end))
            except Exception, exp:
                import ipdb;
                ipdb.set_trace()

            remover = np.ones([len(filter_list)], dtype=bool)
            remover[ev_indx] = False
            # import ipdb;
            # ipdb.set_trace()
            name_subgroup = event_cl.catalog[ev_indx]
            create_nc_subgroup(group, name_subgroup, event_cl, ev_indx)
            # import ipdb;
            # ipdb.set_trace()
            event_cl.rmev(remover)
        else:
            continue

    print event_cl.number_events
    # print event_cl.name

    # =====>>>>> otherwise create new event group
    for indx, event in enumerate(event_cl.name):
        group_name = '%s.%02i.%02i-%02i:%02i:%02i-%04i' \
                     % (event_cl.year[indx], event_cl.month[indx],
                        event_cl.day[indx], event_cl.hr[indx],
                        event_cl.minu[indx], event_cl.sec[indx],
                        event_cl.flinn[indx])

        cprint('STFhelpers.py', '[compare]', bc.grey, '...adding new event to netCDF: %s' % group_name)

        # group name should be uniform throughout the netCDF
        event_group = volcgrp.createGroup('%s' % group_name)

        # assign attributes for comaring later to other subgroups (SCARDEC, LOCAL)
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

    return volcgrp


# ------------------ stf_read_dmt


def stf_read_dmt(stf_input, tau):
    """
    Read STF(s) from the DMT directory
    :param all_events:
    :param inp_ffproc:
    :param ev:
    :return:
    """
    stf = Trace()
    stf_tax = np.arange(0, 100, 1. / stf_input.ph_sampling_rate)
    # Attention: tau is half duration
    stf_data = np.exp(-(stf_tax / tau) ** 2) * \
               (stf_tax / tau) ** 2
    stf_data = stf_data / np.sum(stf_data) * stf_input.ph_sampling_rate

    stf.data = stf_data
    stf.stats.sampling_rate = stf_input.ph_sampling_rate
    
    address_stf = 'DMT'
    number_grps = 1
    grpaff_arr_all = False
    stfs = [stf]
    stf_taxs = [stf_tax]
    return number_grps, address_stf, grpaff_arr_all, stfs, stf_taxs


# ------------------ sdr2mt


def sdr2mt(strike, dip, rake, M0):
    '''
    compare to https://github.com/g2e/seizmo/blob/master/cmt/sdr2mt.m

    convert strike-dip-rake to HRV mt
    f = strike, d = dip, l = rake

    Mrr =  Mzz =  Mo sin 2d sin l
    Mtt =  Mxx = -Mo (sin d cos l sin 2f +     sin 2d sin l (sin f)^2 )
    Mpp =  Myy =  Mo (sin d cos l sin 2f -     sin 2d sin l (cos f)^2 )
    Mrt =  Mxz = -Mo (cos d cos l cos f  +     cos 2d sin l sin f )
    Mrp = -Myz =  Mo (cos d cos l sin f  -     cos 2d sin l cos f )
    Mtp = -Mxy = -Mo (sin d cos l cos 2f + 0.5 sin 2d sin l sin 2f )

    :param strike:
    :param dip:
    :param rake:
    :return: mt
    '''

    # get it in radian
    strike_rad = np.deg2rad(strike)
    dip_rad = np.deg2rad(dip)
    rake_rad = np.deg2rad(rake)

    Mrr = M0*(np.sin(2 * dip_rad) * np.sin(rake_rad))

    Mtt = -M0*(np.sin(dip_rad) * np.cos(rake_rad) * np.sin(2 * strike_rad) + \
          np.sin(2 * dip_rad) * np.sin(rake_rad) * np.sin(strike_rad) * np.sin(strike_rad))

    Mpp = M0*(np.sin(dip_rad) * np.cos(rake_rad) * np.sin(2 * strike_rad) - \
          np.sin(2 * dip_rad) * np.sin(rake_rad) * np.sin(strike_rad) * np.sin(strike_rad))

    Mrt = -M0*(np.cos(dip_rad) * np.cos(rake_rad) * np.cos(strike_rad) + \
          np.cos(2 * dip_rad) * np.sin(rake_rad) * np.sin(strike_rad))

    Mrp = M0*(np.cos(dip_rad) * np.cos(rake_rad) * np.sin(strike_rad) - \
          np.cos(2 * dip_rad) * np.sin(rake_rad) * np.cos(strike_rad))

    Mtp = -M0*(np.sin(dip_rad) * np.cos(rake_rad) * np.cos(2 * strike_rad) + \
          0.5 * np.sin(2 * dip_rad) * np.sin(rake_rad) * np.cos(2 * strike_rad))

    return [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]


# ------------------ half_duration from graph


def hd_calculator(x, y):
    '''
    calculates the cumulative integral and decides when it's half -> half duration
    :param x:
    :param y:
    :return:
    '''

    # initial : scalar, optional
    # If given, uses this value as the first value in the returned result.
    # Typically this value should be 0. Default is None, which means no
    # value at x[0] is returned and res has one element less than
    # y along the axis of integration.

    y_int = cumtrapz(y, x, initial=0)
    y_sum = trapz(y, x)
    y_halfsum = y_sum / 2

    # find closest value to half of th integral
    indx, value = min(enumerate(y_int), key=lambda x: abs(x[1] - y_halfsum))
    half_duration = x[indx]

    return half_duration


# ------------------ mag_duration


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


# ------------------ geocen_array


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


# ------------------ bc


class bc:
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


# ------------------ get_time


def get_time():
    time = datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S')
    return time


# ------------------ cprint


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

    print bc.dgreen + get_time() + bc.end, \
        bc.dmagenta + ho_nam + bc.end, \
        bc.cyan + code + bc.end, \
        bc.bold + bc.white + type_info + bc.end, \
        bc_color + text + bc.end


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
