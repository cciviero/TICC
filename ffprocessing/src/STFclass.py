#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  STFclass.py
#   Purpose:   Source Time Function class
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
from glob import glob
import matplotlib.pyplot as plt
import math
from mpl_toolkits.basemap import Basemap
import numpy as np
from obspy import read, read_events, UTCDateTime
from obspy.core import Stream, Trace
from obspy.imaging.beachball import beach
from obspy.signal.interpolation import lanczos_interpolation
from obspy.clients.iris import Client
import os
import sys

from src.utility_codes import SphericalNearestNeighbour, bc, cprint

# ------------------ conv_stf ---------------------------


def conv_stf(address_stf, all_stations, inp_ffproc, all_events, ev):
    """
    Decide which method of STF reading should be used and convolve the STF
    with the synthetic waveforms
    :param address_stf:
    :param all_stations:
    :param inp_ffproc:
    :param all_events:
    :param ev:
    :return:
    """
    if inp_ffproc.stf_mode == 'read':
        cprint('STFclass.py', '[STF]', bc.lgreen, 'Read archived STFs!')
        all_stfs = conv_stf_cat(address_stf=address_stf,
                                all_stations=all_stations,
                                inp_ffproc=inp_ffproc,
                                all_events=all_events,
                                ev=ev)
    elif inp_ffproc.stf_mode == 'read_dmt':
        cprint('STFclass.py', '[STF]', bc.lgreen, 'Read STFs in DMT archives!')
        all_stfs = conv_stf_dmt(all_stations=all_stations,
                                inp_ffproc=inp_ffproc,
                                all_events=all_events,
                                ev=ev)
    elif inp_ffproc.stf_mode == 'read_scardec':
        cprint('STFclass.py', '[STF]', bc.lgreen, 'Read STFs in SCARDEC archives!')
        all_stfs = conv_stf_scardec(address_stf=address_stf, all_stations=all_stations,
                                inp_ffproc=inp_ffproc,
                                all_events=all_events,
                                ev=ev)

    elif inp_ffproc.stf_mode == 'none':
        cprint('STFclass.py', '[STF]', bc.dgreen, 'No STF convolution!')
        all_stfs = False
    else:
        cprint('STFclass.py', '[ERROR] [STF]', bc.dred,
               '%s has not implemented!' % inp_ffproc.stf_mode)
        sys.exit()
    return all_stfs

# ------------------ conv_stf_cat ---------------------------


def conv_stf_cat(address_stf, all_stations, inp_ffproc, all_events, ev):
    """
    Reading STF(s) and assigning the groups
    :param address_stf:
    :param all_stations:
    :param inp_ffproc:
    :param all_events:
    :param ev:
    :return:
    """
    all_stfs = STFclass()
    try:
        postproc_qual = glob(os.path.join(address_stf, 'outfiles',
                                          'postproc.qual'))
        check_qual = np.loadtxt(postproc_qual[0])
    except Exception as e:
        cprint('STFclass.py', '[NO quality assigned]', bc.lred,
               'Continue with the inverted STF!')
        check_qual = 10.0

    path_stf = glob(os.path.join(address_stf, 'data', 'stf.0?'))
    if check_qual >= inp_ffproc.stf_min_qual and path_stf:
        cprint('STFclass.py', '[STF] [QUAL]', bc.orange, '%s' % check_qual)
        all_stfs.stf_read(address_stf=address_stf,
                          inp_ffproc=inp_ffproc,
                          add_save=os.path.join(inp_ffproc.output_dir,
                                                all_events.name[ev]))

    elif check_qual < inp_ffproc.stf_min_qual or check_qual == 10.0 or path_stf == []:
        cprint('STFclass.py', '[STF] [BAD INPUT]', bc.lmagenta,
               'STF has a low quality, no qualitiy, or STF(s) do(es) not exist at all.\n Please check your '
               'STF(s). A STF from halfduration (check: obspyDMT) will be created instead.')
        all_stfs.stf_read_dmt(all_events=all_events,
                              inp_ffproc=inp_ffproc,
                              ev=ev)

    all_stfs.stf_grouping(all_stations=all_stations,
                          all_events=all_events,
                          ev=ev,
                          add_save=os.path.join(inp_ffproc.output_dir,
                                                all_events.name[ev]))
    # ========================== OLD
    # syn_waveforms = all_stfs.stf_convolve(
    #     syn_waveforms=syn_waveforms,
    #     all_stations=all_stations,
    #     add_save=os.path.join(inp_ffproc.output_dir, all_events.name[ev]),
    #     evname=all_events.name[ev])

    # cutting the Matched filters based on the longest IR duration, preset
    # and offset
    # ==========================

    return all_stfs

# ------------------ conv_stf_dmt ---------------------------


def conv_stf_dmt(all_stations, inp_ffproc, all_events, ev):
    """
    Reading the STF(s), creating dirac STF, convolve the STF and
    cut the waveforms
    :param all_stations:
    :param syn_waveforms:
    :param inp_ffproc:
    :param all_events:
    :param ev:
    :return:
    """
    all_stfs = STFclass()
    all_stfs.stf_read_dmt(all_events=all_events,
                          inp_ffproc=inp_ffproc,
                          ev=ev)
    #all_stfs.stf_dirac(srate=inp_ffproc.ph_sampling_rate,
    #                   syn_nsample=syn_waveforms[0].stats.npts)
    all_stfs.stf_grouping(all_stations=all_stations,
                          all_events=all_events,
                          ev=ev,
                          add_save=os.path.join(inp_ffproc.output_dir,
                                                all_events.name[ev]))

    # syn_waveforms = all_stfs.stf_convolve(
    #     syn_waveforms=syn_waveforms,
    #     add_save=os.path.join(inp_ffproc.output_dir, all_events.name[ev]),
    #     evname=all_events.name[ev])
    # # cutting the Matched filters based on the longest IR duration, preset
    # # and offset
    # for sta in range(len(syn_waveforms)):
    #     start_time = syn_waveforms[sta].stats.starttime
    #     len_tr = inp_ffproc.ph_preset + \
    #              inp_ffproc.lgb_filt_duration[0] + \
    #              inp_ffproc.ph_offset
    #     syn_waveforms[sta] = syn_waveforms[sta].slice(
    #         starttime=start_time,
    #         endtime=start_time+len_tr)
    #     all_stations.syn_start[sta] = syn_waveforms[sta].stats.starttime - \
    #                                   all_events.datetime[ev]
    # return syn_waveforms, all_stfs
    return all_stfs


# ------------------ conv_stf_scardec -----------------------


def conv_stf_scardec(address_stf, all_stations, inp_ffproc, all_events, ev):
    '''

    :param all_stations:
    :param inp_ffproc:
    :param all_events:
    :param ev:
    :return:
    '''
    all_stfs = STFclass()

    # find event in scardec database
    
    event_name = '%s%02d%02d_%02d%02d*' % (all_events.datetime[ev].year,
                                             all_events.datetime[ev].month,
                                             all_events.datetime[ev].day,
                                             all_events.datetime[ev].hour,
                                             all_events.datetime[ev].minute)
    # event_name = '%s%02d%02d_%02d%02d%02d' % (all_events.datetime[ev].year,
    #                                          all_events.datetime[ev].month,
    #                                          all_events.datetime[ev].day,
    #                                          all_events.datetime[ev].hour,
    #                                          all_events.datetime[ev].minute,
    #                                          all_events.datetime[ev].second)
    path_2_scardec_ev = glob(os.path.join(inp_ffproc.stf_scardec_db, '*'+event_name+'*'))
    # import ipdb; ipdb.set_trace()
    if len(path_2_scardec_ev) == 0:
        # if not found scardec event then go back to archive STF
        # cprint('STFclass.py', '[NO SCARDEC STF]', bc.lred, f'No corresponding SCARDEC STF exists. Event {event_name} will be excluded. Continue.')
        # return None

        cprint('STFclass.py', '[NO SCARDEC STF]', bc.lred, 'No corresponding SCARDEC STF exists. '
                                                           'Will go back to archived STFs! Continue.')
        all_stfs = conv_stf_cat(address_stf=address_stf,
                     all_stations=all_stations,
                     inp_ffproc=inp_ffproc,
                     all_events=all_events,
                     ev=ev)
    else:
        all_stfs.stf_read_scardec(path_2_scardec_ev[0], all_events, inp_ffproc, ev)

        # change moment tensor etc... accordingly to SCARDEC
        path_2_scardec_stf = glob(os.path.join(path_2_scardec_ev[0], 'fctoptsource_*'))[0]
        scardec_cat_info = read_scardec_header(path_2_scardec_stf)

        # ev_year, ev_month, ev_day, ev_hour, ev_minute,
        # ev_sec, ev_microsec, ev_julday, ev_lat, ev_lon,
        # ev_dp, ev_mom, ev_mag, moment_1, ev_qual,
        # flinn

        all_events.inv_dp[ev] = scardec_cat_info[10]
        # moment_1 =  [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]
        all_events.mrr[ev] = scardec_cat_info[13][0]
        all_events.mtt[ev] = scardec_cat_info[13][1]
        all_events.mpp[ev] = scardec_cat_info[13][2]
        all_events.mrt[ev] = scardec_cat_info[13][3]
        all_events.mrp[ev] = scardec_cat_info[13][4]
        all_events.mtp[ev] = scardec_cat_info[13][5]

    all_stfs.stf_grouping(all_stations=all_stations,
                          all_events=all_events,
                          ev=ev,
                          add_save=os.path.join(inp_ffproc.output_dir,
                                                all_events.name[ev]))

    return all_stfs

# ------------------ cut_syn_only ---------------------------


def cut_syn_only(all_stations, syn_waveforms, inp_ffproc, all_events, ev):
    """
    Just cuts the synthetic waveform without any STF convolution
    :param all_stations:
    :param syn_waveforms:
    :param inp_ffproc:
    :param all_events:
    :param ev:
    :return:
    """
    # dummy stf's such that all works anyways
    all_stfs = STFclass()
    all_stfs.stf_read_dmt(all_events=all_events,
                          inp_ffproc=inp_ffproc,
                          ev=ev)
    all_stfs.stf_dirac(srate=inp_ffproc.ph_sampling_rate,
                       syn_nsample=syn_waveforms[0].stats.npts)
    all_stfs.stf_grouping(all_stations=all_stations,
                          all_events=all_events,
                          ev=ev,
                          add_save=os.path.join(inp_ffproc.output_dir,
                                                all_events.name[ev]))
    # cutting the Matched filters based on the longest IR duration, preset
    # and offset
    for sta in range(len(syn_waveforms)):
        start_time = syn_waveforms[sta].stats.starttime
        len_tr = inp_ffproc.ph_preset + \
                 inp_ffproc.lgb_filt_duration[0] + \
                 inp_ffproc.ph_offset
        syn_waveforms[sta] = syn_waveforms[sta].slice(
            starttime=start_time,
            endtime=start_time+len_tr)
        all_stations.syn_start[sta] = syn_waveforms[sta].stats.starttime - \
                                      all_events.datetime[ev]
    return syn_waveforms, all_stfs


def read_scardec_header(path_2_scardec_stf):

    with open(path_2_scardec_stf) as myfile:
        # first two lines are the header I need the information
        head = [next(myfile).splitlines() for x in range(2)]

    line = head[0][0]
    year = int(line.split(' ')[0])
    line_1 = head[0][0]
    line_1 = line_1.split(' ')
    line_2 = head[1][0]
    line_2 = line_2.split(' ')
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

    ev_lat = float(line_1[6])
    ev_lon = float(line_1[7])

    ev_dp = float(line_2[0])
    ev_mom = float(line_2[1])
    ev_mag = float(line_2[2])

    ev_strike1 = int(line_2[3])
    ev_dip1 = int(line_2[4])
    ev_rake1 = int(line_2[5])
    moment_1 = sdr2mt(ev_strike1, ev_dip1, ev_rake1, ev_mom)

    ev_qual = 3.5

    # calculate Flinn region
    client = Client()
    try:
        flinn_unformat = client.flinnengdahl(lat=ev_lat, lon=ev_lon,
                                             rtype="code")
        flinn = "{0:0=4d}".format(flinn_unformat)
    except Exception as exp:
        sys.exit(12)
    output_list = [ev_year, ev_month, ev_day, ev_hour, ev_minute, ev_sec,
                   ev_microsec, ev_julday, ev_lat, ev_lon, ev_dp, ev_mom, ev_mag, moment_1, ev_qual, flinn]
    return output_list

# ------------------ STFclass ---------------------------


class STFclass:
    def __init__(self):
        self.number_grps = 0
        self.address_stf = []
        self.grpaff_arr_all = False
        self.stfs = Stream()
        self.taxs = []
        self.grps = []
        self.dirac = False
        self.dirac_time = False

    def stf_read(self, address_stf, inp_ffproc, add_save=False, lancz_a=12):
        """
        Read STF(s) from the directory (similar to the original method)
        :param address_stf:
        :param inp_ffproc:
        :param add_save:
        :param lancz_a:
        :return:
        """
        stf_counter = 1
        stfs_glob = glob(os.path.join(address_stf, 'data', 'stf.0?'))
        stfs_glob.sort()
        for stftmp in stfs_glob:
            stf = read(stftmp, format='SAC')[0]
            # resampling STF
            starttime = stf.stats.starttime.timestamp
            endtime = stf.stats.endtime.timestamp
            dt = 1./inp_ffproc.ph_sampling_rate
            stf.data = lanczos_interpolation(
                np.require(stf.data, requirements=["C"]),
                starttime, stf.stats.delta, starttime, dt,
                int(math.floor((endtime - starttime) / dt)) + 1,
                a=lancz_a, window="blackman")

            stf.stats.starttime = UTCDateTime(starttime)
            stf.stats.sampling_rate = inp_ffproc.ph_sampling_rate

            if inp_ffproc.kernel_stf_gauss and stf_counter == 1:
                cum_stf = np.cumsum(abs(stf.data))
                cum_stf_norm = cum_stf/cum_stf[-1]
                stfk_half_duration = \
                    np.where(cum_stf_norm > 0.9)[0][0]/2. / \
                    stf.stats.sampling_rate

                stfk_tax = np.linspace(0,
                                       0 +
                                       (stf.stats.npts - 1) /
                                       stf.stats.sampling_rate,
                                       stf.stats.npts)

                # XXX kernel code can not incorporate this?
                # stfk_tax = np.linspace(stf.stats.sac.b,
                #                        stf.stats.sac.b +
                #                        (stf.stats.npts - 1) /
                #                        stf.stats.sampling_rate,
                #                        stf.stats.npts)

                stfk_data = np.exp(-(stfk_tax/stfk_half_duration)**2) * \
                                  (stfk_tax/stfk_half_duration)**2
                # XXX kernel code: energy normalized?
                stfk_data = stfk_data/np.sum(stfk_data)*stf.stats.sampling_rate
                stfk_fio = open(os.path.join(add_save, 'stf.dat'), 'w')
                stfk_fio.writelines('%s  %s\n' % (len(stfk_data),
                                    1./inp_ffproc.ph_sampling_rate))
                for ic in range(len(stfk_data)):
                    stfk_fio.writelines('%s\n' % stfk_data[ic])
                stfk_fio.close()

                stf.data = stfk_data
                self.address_stf.append(stftmp)
                self.taxs.append(stfk_tax)
                self.stfs.append(stf)
                self.number_grps = 1
            elif inp_ffproc.kernel_stf_gauss and stf_counter > 1:
                continue
            else:
                self.address_stf.append(stftmp)
                self.taxs.append(np.linspace(stf.stats.sac.b,
                                             stf.stats.sac.b +
                                             (stf.stats.npts - 1) /
                                             stf.stats.sampling_rate,
                                             stf.stats.npts))
                self.stfs.append(stf)
                self.number_grps += 1
            stf_counter += 1

        if self.number_grps > 1:
            try:
                # loading the file that contains the grouping information
                self.grpaff_arr_all = np.loadtxt(os.path.join(address_stf,
                                                              'outfiles',
                                                              'ampinv.grpaff'),
                                                 dtype='S', skiprows=1)
            except Exception as e:
                self.number_grps = 1

    def stf_read_dmt(self, all_events, inp_ffproc, ev):
        """
        Read STF(s) from the DMT directory
        :param all_events:
        :param inp_ffproc:
        :param ev:
        :return:
        """
        stf = Trace()
        stf_tax = np.arange(0, 100, 1./inp_ffproc.ph_sampling_rate)
        # Attention: tau is half duration
        stf_data = np.exp(-(stf_tax/all_events.tau[ev])**2) * \
                         (stf_tax/all_events.tau[ev])**2
        stf_data = stf_data/np.sum(stf_data)*inp_ffproc.ph_sampling_rate

        stf.data = stf_data
        stf.stats.sampling_rate = inp_ffproc.ph_sampling_rate
        self.address_stf.append('DMT')
        self.taxs.append(stf_tax)
        self.stfs.append(stf)
        self.number_grps = 1

    def stf_read_scardec(self, path_2_scardec_ev, all_events, inp_ffproc, ev, lancz_a=12):
        '''

        :param path_2_scardec_ev:
        :param all_events:
        :param inp_ffproc:
        :param ev:
        :return:
        '''

        # read stf metadata information for fctOPTsource
        path_2_scardec_stf = glob(os.path.join(path_2_scardec_ev, 'fctoptsource_*'))[0]

        # read stf graph and convert to trace
        # skip the first two lines which are more conveniently read with obspy
        opt_stf = np.loadtxt(path_2_scardec_stf, skiprows=2)
        tr_data_opt = opt_stf[:, 1]

        # normalise
        tr_eng_norm_opt = tr_data_opt / np.sum(tr_data_opt) * 14
        opt_tr = Trace(tr_eng_norm_opt)

        # resampling STF (SCARDEC sampling rate is 14)
        opt_tr.stats.sampling_rate = 14

        # starttime is eventtime
        starttime = all_events.datetime[ev].timestamp + opt_stf[:,0][0]
        # endttime is time of stf + eventtime
        endtime = all_events.datetime[ev].timestamp + opt_stf[:,0][-1]

        dt = 1. / inp_ffproc.ph_sampling_rate
        # lanczos_interpolation(data, old_start, old_dt, new_start, new_dt, new_npts,
        # a, window='lanczos', *args, **kwargs)[source]
        opt_tr.data = lanczos_interpolation(np.require(opt_tr.data, requirements=["C"]),
                                            starttime,
                                            opt_tr.stats.delta,
                                            starttime,
                                            dt,
                                            int(math.floor((endtime - starttime) / dt)) + 1,
                                            a=lancz_a,
                                            window="blackman")

        opt_tr.stats.starttime = UTCDateTime(starttime)
        opt_tr.stats.sampling_rate = inp_ffproc.ph_sampling_rate

        self.address_stf.append('SCARDEC')
        self.taxs.append(np.linspace(opt_stf[:,0][0],
                                     opt_stf[:, 0][0] + (opt_tr.stats.npts - 1) / opt_tr.stats.sampling_rate,
                                     opt_tr.stats.npts))
        self.stfs.append(opt_tr)
        self.number_grps = 1

    def stf_dirac(self, srate, syn_nsample):
        """
        creating Dirac STF
        :param srate:
        :param syn_nsample:
        :return:
        """
        self.dirac_time = np.linspace(0, syn_nsample - 1, syn_nsample)
        self.dirac_time /= srate
        self.dirac = np.zeros(syn_nsample)
        self.dirac[0] = 1

    def stf_grouping(self, all_stations, all_events, ev, add_save):
        """
        Finding the STF groups that each seismogram belongs to
        :param all_stations:
        :param all_events:
        :param ev:
        :param add_save:
        :return:
        """
        if self.number_grps == 1:
            self.grps = np.ones(all_stations.number_stations, dtype='int')
        elif self.number_grps > 1:
            lats = self.grpaff_arr_all[:, 4].astype(np.float)
            lons = self.grpaff_arr_all[:, 5].astype(np.float)
            els = np.zeros(len(lats))
            grps = self.grpaff_arr_all[:, 1].astype(np.int)
            aff_orig_sta = SphericalNearestNeighbour(lats, lons, els)
            self.grps = \
                grps[aff_orig_sta.query(all_stations.lat,
                                        all_stations.lon,
                                        np.zeros(len(all_stations.lat)))[1]]
        else:
            cprint('STFclass.py', '[ERROR] [STF]', bc.dred,
                   'Something is wrong with grouping. EXIT!')
            sys.exit()
        # Plotting the group affiliations
        self.plot_grpaff(all_events=all_events,
                         ev=ev,
                         all_stations=all_stations,
                         add_save=add_save)

    def stf_convolve(self, syn_waveforms, all_stations, add_save, evname):
        """
        convolve SYN*STF to generate Matched Filters
        :param syn_waveforms:
        :param all_stations:
        :param add_save:
        :param evname:
        :return:
        """
        for sta in range(len(syn_waveforms)):
            # First time sample of the STF
            tb_stf = self.taxs[self.grps[sta] - 1][0]
            syn_npts = syn_waveforms[sta].stats.npts
            srate = syn_waveforms[sta].stats.sampling_rate
            syn_start = syn_waveforms[sta].stats.starttime
            syn_waveforms[sta].data = \
                np.convolve(syn_waveforms[sta].data,
                            self.stfs[self.grps[sta] - 1])
            syn_waveforms[sta].stats.starttime = syn_start + tb_stf
            start_time = syn_waveforms[sta].stats.starttime - tb_stf
            syn_waveforms[sta] = \
                syn_waveforms[sta].slice(
                    starttime=start_time,
                    endtime=start_time + (syn_npts-1)/srate)
        return syn_waveforms

        # Plotting the spectrum
        # plt.subplot(211)
        # plt.plot(syn_waveforms[sta].times(), syn_waveforms[sta].data, 'k')
        # plt.subplot(212)
        # from src.potential_scripts.plot_spectrum import spectrum_calc
        # freqs, spec_tr = spectrum_calc(syn_waveforms[sta])
        # plt.loglog(freqs, spec_tr, 'k')
        # freqs, spec_tr = spectrum_calc(self.stfs[self.grps[sta] - 1])
        # plt.loglog(freqs, spec_tr, 'b--')

        # plt.subplot(211)
        # plt.plot(syn_waveforms[sta].times(), syn_waveforms[sta].data, 'r')
        # plt.subplot(212)
        # freqs, spec_tr = spectrum_calc(syn_waveforms[sta])
        # plt.loglog(freqs, spec_tr, c='r')

    def plot_grpaff(self, all_events, ev, all_stations, add_save):
        """
        Plot group affiliations
        :param all_events:
        :param ev:
        :param all_stations:
        :param add_save:
        :return:
        """
        from matplotlib import gridspec
        plt.ion()
        fig = plt.figure(frameon=False)
        fig.set_size_inches(15, 10)
        gs = gridspec.GridSpec(4, 4)
        fig.add_subplot(gs[0:3, 0:4])
        try:
            mymap = Basemap(projection='aeqd',
                            lat_0=round(all_events.lat[ev], 1),
                            lon_0=round(all_events.lon[ev], 1))
        except Exception as e:
            mymap = Basemap(projection='aeqd',
                            lat_0=round(all_events.lat[ev], 2),
                            lon_0=round(all_events.lon[ev], 2))
        mymap.fillcontinents()
        mymap.drawcoastlines()

        x, y = mymap(all_events.lon[ev], all_events.lat[ev])
        focmecs = [all_events.mrr[ev],
                   all_events.mtt[ev],
                   all_events.mpp[ev],
                   all_events.mrt[ev],
                   all_events.mrp[ev],
                   all_events.mtp[ev]]

        try:
            ax = plt.gca()
            b = beach(focmecs, xy=(x, y),
                      width=int(2e6), linewidth=1,
                      facecolor='r')
            b.set_zorder(50)
            ax.add_collection(b)
        except Exception as e:
            mymap.scatter(x, y, c='r', edgecolor='none', marker='o', zorder=20, s=100) 

        # better color?
        colors = ['b', 'r', 'g', 'c', 'm', 'y', 'k', 'b', 'r', 'g']
        cprint('STFclass.py', '[STF]', bc.dmagenta,
               '#STF(s): %s' % self.number_grps)
        for grp in range(self.number_grps):
            x, y = mymap(all_stations.lon[self.grps == (grp + 1)],
                         all_stations.lat[self.grps == (grp + 1)])
            cols = [colors[grp]]*len(x)
            mymap.scatter(x, y, c=cols, marker='v', s=100, zorder=10)
        plt.title('STF Groups', size=18, weight='bold')

        ax_stf = fig.add_subplot(gs[3, 1:3])
        for grps in range(self.number_grps):
            ax_stf.plot(self.taxs[grps], self.stfs[grps],
                        c=colors[grps], lw=3,
                        label='Group %s' % (grps+1))
        ax_stf.set_xlabel('Time [sec]', size=15, weight='bold')
        ax_stf.set_ylabel('Energy', size=15, weight='bold')
        ax_stf.legend()
        plt.savefig(os.path.join(
            add_save, 'stfs_grp_%s.png' % all_events.name[ev]),
            format='png', bbox_inches='tight')
        plt.clf()
        plt.ioff()
        plt.close()


#  ---------------- SCARDEC helpers --------------------------

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