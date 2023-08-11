#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  utility_codes.py
#   Purpose:   Collection of utility codes for raydata_raymatrix.py
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

# Required Python modules will be imported in this part.

import configparser
import glob
import math as m_math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
# from mpl_toolkits.basemap import Basemap
import numpy as np
from obspy import UTCDateTime
try:
    from obspy.geodetics import locations2degrees
except Exception as e:
    from obspy.core.util import locations2degrees
#from obspy.taup import getTravelTimes
from obspy.taup import tau
from datetime import datetime
import os
import socket
import subprocess
import shutil
import sys
import time
# import pyvtk as pvtk

# ###################### GENERAL INPUT ##################################
# Transparency of the scatterer plots
global alpha_plt
alpha_plt = 0.1

# ###################### InpClass #############################


class InpClass():
    def __init__(self, inp):
        Config_ini = configparser.ConfigParser()
        Config_ini.read(inp)
        # general
        self.plot_statistics = \
            eval(Config_ini.get('section_general', 'plot_statistics'))
        self.run_raydata = \
            eval(Config_ini.get('section_general', 'run_raydata'))
        self.run_raymatrix = \
            eval(Config_ini.get('section_general', 'run_raymatrix'))
        self.parallel_exec = \
            eval(Config_ini.get('section_general', 'parallel_exec'))
        self.np_req = \
            int(Config_ini.get('section_general', 'np_req'))
        self.selected_events_add = \
            eval(Config_ini.get('section_general', 'selected_events_add'))
        # phase
        self.phase = eval(Config_ini.get('section_phase', 'phase'))
        self.events_dir = eval(Config_ini.get('section_phase', 'events_dir'))
        self.min_epi = float(Config_ini.get('section_phase', 'min_epi'))
        self.max_epi = float(Config_ini.get('section_phase', 'max_epi'))
        self.req_bands = eval(Config_ini.get('section_phase', 'req_bands'))
        self.all_events = eval(Config_ini.get('section_phase', 'all_events'))
        # raydata
        self.input_file_name_part = \
            eval(Config_ini.get('section_raydata', 'input_file_name_part'))
        self.twinned = eval(Config_ini.get('section_raydata', 'twinned'))
        self.max_num_arrival = \
            int(Config_ini.get('section_raydata', 'max_num_arrival'))
        self.delay_wrt_first_arrival = \
            float(Config_ini.get('section_raydata', 'delay_wrt_first_arrival'))
        self.corr_io_list = \
            eval(Config_ini.get('section_raydata', 'corr_io_list'))
        self.crust_corr = \
            eval(Config_ini.get('section_raydata', 'crust_corr'))
        self.sigma_fix = \
            eval(Config_ini.get('section_raydata', 'sigma_fix'))
        self.sigma_cc = \
            eval(Config_ini.get('section_raydata', 'sigma_cc'))
        self.sigma_factor = \
            eval(Config_ini.get('section_raydata', 'sigma_factor'))
        # raymatrix
        self.vp_vs_Qs = eval(Config_ini.get('section_raymatrix', 'vp_vs_Qs'))
        self.kernel_quad_km = \
            float(Config_ini.get('section_raymatrix', 'kernel_quad_km'))
        self.vertex_file = \
            eval(Config_ini.get('section_raymatrix', 'vertex_file'))
        self.facet_file = \
            eval(Config_ini.get('section_raymatrix', 'facet_file'))
        # criteria
        self.min_depth = float(Config_ini.get('section_criteria', 'min_depth'))
        self.max_depth = float(Config_ini.get('section_criteria', 'max_depth'))
        self.min_xcorr = float(Config_ini.get('section_criteria', 'min_xcorr'))
        self.max_xcorr = float(Config_ini.get('section_criteria', 'max_xcorr'))
        self.bg_model = eval(Config_ini.get('section_criteria', 'bg_model'))
        self.check_clip = \
            eval(Config_ini.get('section_criteria', 'check_clip'))
        # check_plot
        self.check_selections = \
            eval(Config_ini.get('section_check_plot', 'check_selections'))
        self.pickle_filt_array_quit = \
            eval(Config_ini.get('section_check_plot',
                                'pickle_filt_array_quit'))
        self.pickle_filt_array = \
            eval(Config_ini.get('section_check_plot', 'pickle_filt_array'))
        self.pickle_filt_array_quit_corr = \
            eval(Config_ini.get('section_check_plot',
                                'pickle_filt_array_quit_corr'))
        self.write_staev_loc = \
            eval(Config_ini.get('section_check_plot', 'write_staev_loc'))
        self.abs_vtk = eval(Config_ini.get('section_check_plot', 'abs_vtk'))
        self.plot_stas_proj = \
            eval(Config_ini.get('section_check_plot', 'plot_stas_proj'))
        self.plot_stas_proj_type = \
            eval(Config_ini.get('section_check_plot', 'plot_stas_proj_type'))
        self.plot_lat_0 = \
            eval(Config_ini.get('section_check_plot', 'plot_lat_0'))
        self.plot_lon_0 = \
            eval(Config_ini.get('section_check_plot', 'plot_lon_0'))
        self.plot_fixed_median = \
            eval(Config_ini.get('section_check_plot', 'plot_fixed_median'))
        self.num_sta_thresh_median = \
            int(Config_ini.get('section_check_plot', 'num_sta_thresh_median'))
        self.plot_selected_stas = \
            eval(Config_ini.get('section_check_plot', 'plot_selected_stas'))
        # section_output
        self.out_path = eval(Config_ini.get('section_output', 'out_path'))
        # import ipdb; ipdb.set_trace()
        # create if not exist
        if not os.path.isdir(self.out_path):
            os.makedirs(self.out_path)

# ###################### print_inp #############################


def print_inp(inp, req_band):
    """
    printing the input class
    :param inp:
    :return:
    """
    cprint('utility_codes.py', '[INPUT INFORMATION]', bc.white,
           'Event DIR: %s' % inp.events_dir)
    cprint('utility_codes.py', '[INPUT INFORMATION]', bc.white,
           'Output DIR: %s' % inp.out_path)
    cprint('utility_codes.py', '[INPUT INFORMATION]', bc.white,
           'Requested Phase: %s' % inp.phase)
    cprint('utility_codes.py', '[INPUT INFORMATION]', bc.white,
           'Requested Band: %s' % req_band)

    if inp.all_events == True:
        cprint('utility_codes.py', '[INPUT INFORMATION]', bc.white,
               'Requested events: ALL')
    else:
        cprint('utility_codes.py', '[INPUT INFORMATION]', bc.white,
               'Requested events: %s' % inp.all_events)

    if inp.run_raydata:
        cprint('utility_codes.py', '[INPUT RAYDATA]', bc.white,
               'INPUT name: %s' % inp.input_file_name_part)
        cprint('utility_codes.py', '[INPUT RAYDATA]', bc.white,
               'Twinned: %s' % inp.twinned)
        cprint('utility_codes.py', '[INPUT RAYDATA]', bc.white,
               'Maximum number of arrivals: %s' % inp.max_num_arrival)
        cprint('utility_codes.py', '[INPUT RAYDATA]', bc.white,
               'Delay wrt the first arrival: %s' % inp.delay_wrt_first_arrival)
        if inp.corr_io_list[0] == 1:
            cprint('utility_codes.py', '[INPUT RAYDATA]', bc.white,
                   'Correction for Ellipticity: %s' % (bc.green + 'Yes' + bc.end))
        else:
            cprint('utility_codes.py', '[INPUT RAYDATA]', bc.white,
                   'Correction for Ellipticity: %s' % (bc.red + 'No' + bc.end))
        if inp.corr_io_list[1] == 1:
            cprint('utility_codes.py', '[INPUT RAYDATA]', bc.white,
                   'Correction for Crustal time: %s' % (bc.green + 'Yes' + bc.end))
        else:
            cprint('utility_codes.py', '[INPUT RAYDATA]', bc.white,
                   'Correction for Crustal time: %s' % (bc.red + 'No' + bc.end))
        if inp.corr_io_list[2] == 1:
            cprint('utility_codes.py', '[INPUT RAYDATA]', bc.white,
                   'Correction for Elevation: %s' % (bc.green + 'Yes' + bc.end))
        else:
            cprint('utility_codes.py', '[INPUT RAYDATA]', bc.white,
                   'Correction for Elevation: %s' % (bc.red + 'No' + bc.end))

    if inp.run_raymatrix:
        cprint('utility_codes.py', '[INPUT RAYMATRIX]', bc.white,'Vp, Vs, Qs : %s' % inp.vp_vs_Qs)
        cprint('utility_codes.py', '[INPUT RAYMATRIX]', bc.white,'Kernel Quadrature (KM): %s' % inp.kernel_quad_km)
        cprint('utility_codes.py', '[INPUT RAYMATRIX]', bc.white,'Vertex file: %s' % inp.vertex_file)
        cprint('utility_codes.py', '[INPUT RAYMATRIX]', bc.white,'Facet file: %s' % inp.facet_file)

    if inp.parallel_exec:
        cprint('utility_codes.py', '[INPUT PARALLEL]', bc.green, 'YES')
        print ('#Processes: %s' % inp.np_req)
    else:
        cprint('utility_codes.py', '[INPUT PARALLEL]', bc.green, 'NO')

    cprint('utility_codes.py', '[INPUT CRITERIA]', bc.white, '%s <= Depth < %s' % (inp.min_depth, inp.max_depth))
    cprint('utility_codes.py', '[INPUT CRITERIA]', bc.white, '%s <= xcorr < %s' % (inp.min_xcorr, inp.max_xcorr))
    cprint('utility_codes.py', '[INPUT CRITERIA]', bc.white, '%s <= epicentral < %s' % (inp.min_epi, inp.max_epi))
    cprint('utility_codes.py', '[INPUT CRITERIA]', bc.white, 'Check clip: %s' % inp.check_clip)
    cprint('utility_codes.py', '[INPUT CRITERIA]', bc.white, 'Background model: %s' % inp.bg_model)
    cprint('utility_codes.py', '[END INPUT]', bc.green, '\n\n')

# ###################### event_filter #############################


def event_filter(inp):
    """
    Reading and filtering events based on the criteria defined by the user
    :param inp:
    :return:
    """
    # import ipdb; ipdb.set_trace()
    cprint('utility_codes.py', '[event_filter]', bc.green, 'Reading the event information and '
                                                   'filter them %s <= depth < %s' % (inp.min_depth, inp.max_depth))
    if not os.path.isdir(inp.events_dir):
        cprint('utility_codes.py', '[event_filter]', bc.bred, '%s is not a valid directory' % inp.events_dir)
        return False, False
    selected_events = np.loadtxt(fname=inp.selected_events_add,
                                 dtype='S', comments='#', delimiter=',')
    event_adds = []
    for i in range(len(selected_events[:, 1])):
        if inp.all_events != True:
            if not selected_events[:, 1][i] in inp.all_events:
                continue
        event_adds.append(
            [os.path.join(inp.events_dir, selected_events[:, 1][i]),
             selected_events[:, 0][i]])
    passed_event_adds = []
    for i in range(len(event_adds)):
        add_flag = False
        try:
            fio_source = open(os.path.join(event_adds[i][0], 'outfiles',
                                           'ffproc.source'), 'r')
            f_source = fio_source.readlines()
            ev_year, ev_julianday, ev_hr, ev_min, ev_sec, ev_msec = \
                f_source[1].split()
            evlat, evlon, catalog_depth, inverted_depth = f_source[3].split()
            try:
                mrr, mtt, mpp, mrt, mrp, mtp = f_source[13].split()
            except Exception as e:
                print ('WARNING: inverted moment tensor was not found: %s' \
                      % event_adds[i][0])
                mrr, mtt, mpp, mrt, mrp, mtp = f_source[7].split()
            ev_date = UTCDateTime(year=int(ev_year), julday=int(ev_julianday))
            ev_date_str = '%4s%3s' % (ev_date.year, ev_date.julday)
            ev_date_str = ev_date_str.replace(' ', '0')
            ev_time = '%2s%2s%2s' % (ev_hr, ev_min, ev_sec)
            ev_time = ev_time.replace(' ', '0')
            ev_id = event_adds[i][1]
            # Check for depth
            if inp.min_depth <= float(inverted_depth) < inp.max_depth:
                add_flag = True
            if add_flag:
                passed_event_adds.append([event_adds[i][0], ev_date_str,
                                          ev_time, ev_id, evlat, evlon,
                                          inverted_depth,
                                          mrr, mtt, mpp, mrt, mrp, mtp])
        except Exception as e:
            cprint('utility_codes.py', '[event_filter]', bc.red, 'ERROR: %s' % e)
    # Sort the events by depth and then pass it
    return sorted(passed_event_adds,
                  key=lambda passed_event_adds: float(passed_event_adds[6]))

# ###################### output_file_reader #############################


def output_file_reader(evinfo, inp, req_band='band01'):
    """
    This function reads one output file and pass all the values
    ATTENTION: there is not filtering at this stage, it reads all...
    :param evinfo:
    :param req_band:
    :return:
    """
    add_output = evinfo[0]
    empty_array = np.array([])
    import ipdb; ipdb.set_trace()
    if not os.path.isfile(os.path.join(add_output, 'outfiles',
                                       'ffproc.ampstt.%s' % req_band)):
        print ("%s is not found!" % os.path.join(add_output, 'outfiles',
                                                'ffproc.ampstt.%s' % req_band))
        return empty_array
    if not os.path.isfile(os.path.join(add_output, 'outfiles',
                                       'ffproc.receivers')):
        print ("%s is not found!" % os.path.join(add_output, 'outfiles',
                                                'ffproc.receivers'))
        return empty_array
    # reads always only the default 22 columns, what comes after that does not matter
    rd_output = np.loadtxt(os.path.join(add_output, 'outfiles',
                                        'ffproc.ampstt.%s' % req_band),
                           dtype='S', comments='#', usecols=range(0, 22))
    # import ipdb; ipdb.set_trace()
    sta_output = np.loadtxt(os.path.join(add_output, 'outfiles',
                                         'ffproc.receivers'),
                            dtype='S', comments='#')

    # To avoid any problem in case that there is only one src-rcv available:
    if np.size(rd_output) == 22:
        rd_output = np.reshape(rd_output, [1, 22])
    if np.size(sta_output) == 8:
        sta_output = np.reshape(sta_output, [1, 8])

    # import ipdb;
    # ipdb.set_trace()
    # these lines change the uncertainty assigned previously in pyffproc
    # if it seems that the inversion does not gives out a 'proper' L-curve
    if inp.sigma_fix:
        for indx, cc in enumerate(rd_output[:, 6]):
            if eval(cc) < inp.sigma_cc:
                rd_output[indx, 9] = rd_output[indx, 9].astype(np.float) * inp.sigma_factor

    # Generate an empty array with size of: rd_output_rows + 9
    new_col_cr = np.empty([np.shape(rd_output)[0], 9], dtype=object)
    new_col_cr[:, 0] = (sta_output[:, 5].astype(np.float) -
                        sta_output[:, 6].astype(np.float))/1000.
    new_col_cr[:, 1] = os.path.basename(add_output)
    new_col_cr[:, 2] = evinfo[1]
    new_col_cr[:, 3] = evinfo[2]
    new_col_cr[:, 4] = evinfo[3]
    new_col_cr[:, 5] = evinfo[4]
    new_col_cr[:, 6] = evinfo[5]
    new_col_cr[:, 7] = evinfo[6]
    new_col_cr[:, 8] = req_band
    # import ipdb; ipdb.set_trace()
    # Now append the rd_output_rows X 9 array to the original rd_output
    # import ipdb; ipdb.set_trace()
    output_sta_evname = np.append(rd_output, new_col_cr, 1)

    return output_sta_evname

# ###################### array_station_filter #############################


def array_station_filter(passed_array, inp):
    """
    Filters the stations in an array based on the required inputs
    :param passed_array:
    :param inp:
    :return:
    """
    cprint('utility_codes.py', '[filtering FFM outputs]', bc.white, '')
    cprint('utility_codes.py', '%s <= xcorr < %s' % (inp.min_xcorr, inp.max_xcorr), bc.white, '')
    cprint('utility_codes.py', '%s <= epicentral < %s' % (inp.min_epi, inp.max_epi), bc.white, '')
    cprint('utility_codes.py', 'check clip: %s' % inp.check_clip, bc.white, '')
    print ('\n')

    cprint('utility_codes.py', '=================================', bc.white, '')
    cprint('utility_codes.py', '#event-station pairs (before): %s' % len(passed_array), bc.white, '')

    empty_array = np.array([])
    # import ipdb; ipdb.set_trace()
    # --------------- XCORRELATION ---------------
    pass_stas_corr_1 = passed_array[
        passed_array[:, 6].astype(np.float) >= inp.min_xcorr]
    if not pass_stas_corr_1.size == 0:
        pass_stas_final = pass_stas_corr_1[
            pass_stas_corr_1[:, 6].astype(np.float) < inp.max_xcorr]
    else:
        return empty_array
    # --------------- EPICENTRAL ---------------
    passed_stas_epi_1 = pass_stas_final[
        pass_stas_final[:, 4].astype(np.float) >= inp.min_epi]
    if not passed_stas_epi_1.size == 0:
        passed_stas_final = passed_stas_epi_1[
            passed_stas_epi_1[:, 4].astype(np.float) < inp.max_epi]
    else:
        return empty_array
    # --------------- CHECK CLIPS ---------------
    if inp.check_clip:
        passed_stas_final = passed_stas_final[
            passed_stas_final[:, 19].astype(np.float) < 0.1]

    if not passed_stas_final.size == 0:
        cprint('utility_codes.py', '#event-station pairs (after):  %s' % len(passed_stas_final), bc.white, '')
        cprint('utility_codes.py', '=================================', bc.white, '')
        return passed_stas_final
    else:
        cprint('utility_codes.py', '#event-station pairs (after):  0', bc.white, '')
        cprint('utility_codes.py', '=================================', bc.white, '')
        return empty_array

# ###################### array_fill_in #############################


def array_fill_in(filt_array, req_band_prev, output_dir):
    srcrcv_previous = \
        np.load(os.path.join(output_dir, 'srcrcv_%s.npy' % req_band_prev))
    srcrcv_id = filt_array[:, 0] + filt_array[:, 26]
    print ('Length source-receiver pairs (before): %s' % len(filt_array))
    decision = np.in1d(srcrcv_id, srcrcv_previous)
    filt_array = filt_array[decision != True]
    filt_array = np.reshape(filt_array, (-1, 31))
    print ('Length source-receiver pairs (after): %s' % len(filt_array))
    return filt_array

# ###################### array_fill_in_writer #############################


def array_fill_in_writer(filt_array, req_band, output_dir):
    # import ipdb; ipdb.set_trace()
    filt_array = np.reshape(filt_array, (-1, 31))
    srcrcv_id = filt_array[:, 0] + filt_array[:, 26]
    np.save(os.path.join(output_dir, 'srcrcv_%s.npy' % req_band),
            srcrcv_id)

# ###################### array_station_filter_mark ############################


def array_station_filter_mark(all_output_files, min_xcorr=-100, max_xcorr=100,
                              min_epi=0., max_epi=360., check_clip=True):
    """
    Filters the stations in an array based on the required inputs
    """
    # --------------- XCORRELATION ---------------
    all_output_files[:, 10, :][all_output_files[:, 6, :].astype(np.float) <
                               min_xcorr] = -1
    all_output_files[:, 10, :][all_output_files[:, 6, :].astype(np.float) >=
                               max_xcorr] = -1
    # --------------- EPICENTRAL ---------------
    all_output_files[:, 10, :][all_output_files[:, 4, :].astype(np.float) <
                               min_epi] = -1
    all_output_files[:, 10, :][all_output_files[:, 4, :].astype(np.float) >=
                               max_epi] = -1
    # --------------- CHECK CLIPS ---------------
    if check_clip:
        all_output_files[:, 10, :] \
            [all_output_files[:, 19, :].astype(np.float) > 0.1] = -1

    return all_output_files

# ###################### check_selection #############################


def check_selection(filt_array, inp, output_dir):
    """
    Check whether the selection procedure worked well
    :param filt_array:
    :param inp:
    :return:
    """
    print ('\n======>> check the selected event-station pairs')
    global alpha_plt

    if len(filt_array) < 100:
        alpha_plt = 1.0
    elif 100 <= len(filt_array) < 10000:
        alpha_plt = 0.2
    elif 10000 <= len(filt_array) < 50000:
        alpha_plt = 0.1
    elif 50000 <= len(filt_array) < 200000:
        alpha_plt = 0.01
    else:
        alpha_plt = 0.002

    plt.ion()
    plt.figure()
    plt.subplot(2, 2, 1)
    plt.scatter(filt_array[:, 2].astype(np.float),
                filt_array[:, 6].astype(np.float), c='r', edgecolors='none',
                alpha=alpha_plt)
    plt.xlabel('Station latitude', size='large', weight='bold')
    plt.ylabel('xcorr factor', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.hlines(inp.min_xcorr, np.min(filt_array[:, 2].astype(np.float)),
               np.max(filt_array[:, 2].astype(np.float)), 'k',
               linestyles='--')
    plt.hlines(inp.max_xcorr, np.min(filt_array[:, 2].astype(np.float)),
               np.max(filt_array[:, 2].astype(np.float)), 'k',
               linestyles='--')
    plt.vlines(50.0, inp.min_xcorr, inp.max_xcorr, 'k', linestyles='-')
    plt.vlines(26.0, inp.min_xcorr, inp.max_xcorr, 'k', linestyles='-')
    plt.vlines(60.0, inp.min_xcorr, inp.max_xcorr, 'k', linestyles='dotted')
    plt.vlines(34.0, inp.min_xcorr, inp.max_xcorr, 'k', linestyles='dotted')

    plt.subplot(2, 2, 2)
    plt.scatter(filt_array[:, 3].astype(np.float),
                filt_array[:, 6].astype(np.float), c='r', edgecolors='none',
                alpha=alpha_plt)
    plt.xlabel('Station longitude', size='large', weight='bold')
    plt.ylabel('xcorr factor', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.hlines(inp.min_xcorr, np.min(filt_array[:, 3].astype(np.float)),
               np.max(filt_array[:, 3].astype(np.float)), 'k',
               linestyles='--')
    plt.hlines(inp.max_xcorr, np.min(filt_array[:, 3].astype(np.float)),
               np.max(filt_array[:, 3].astype(np.float)), 'k',
               linestyles='--')
    plt.vlines(-127.0, inp.min_xcorr, inp.max_xcorr, 'k', linestyles='-')
    plt.vlines(-60.0, inp.min_xcorr, inp.max_xcorr, 'k', linestyles='-')
    plt.vlines(-12.0, inp.min_xcorr, inp.max_xcorr, 'k', linestyles='dotted')
    plt.vlines(41.0, inp.min_xcorr, inp.max_xcorr, 'k', linestyles='dotted')

    plt.subplot(2, 2, 3)
    plt.scatter(filt_array[:, 4].astype(np.float),
                filt_array[:, 6].astype(np.float), c='b', edgecolors='none',
                alpha=alpha_plt)
    plt.xlabel('Epicentral distance', size='large', weight='bold')
    plt.ylabel('xcorr factor', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.vlines(inp.min_epi, inp.min_xcorr, inp.max_xcorr, 'k', linestyles='--')
    plt.vlines(inp.max_epi, inp.min_xcorr, inp.max_xcorr, 'k', linestyles='--')

    plt.subplot(2, 2, 4)
    plt.scatter(filt_array[:, 19].astype(np.float),
                filt_array[:, 6].astype(np.float), c='r', edgecolors='none',
                alpha=alpha_plt)
    plt.xlabel('Clip value', size='large', weight='bold')
    plt.ylabel('xcorr factor', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')

    # Importing as a MATLAB object for further processing
    #py2mat(filt_array[:, 2].astype(np.float), 'stlat', 'RESULTS/stlat')
    #py2mat(filt_array[:, 3].astype(np.float), 'stlon', 'RESULTS/stlon')
    #py2mat(filt_array[:, 4].astype(np.float), 'epi_dist',
    #       'RESULTS/epi_dist')
    #py2mat(filt_array[:, 6].astype(np.float), 'xcorr', 'RESULTS/xcorr')

    """
    ----------------------
    Useful matlab commands
    ----------------------
    load('stlat.mat')
    load('stlon.mat')
    load('epi_dist.mat')
    load('xcorr.mat')

    dat = [epi_dist; xcorr];
    dat = dat';
    %n = log(hist3(dat, [123, 20])); % default is to 10x10 bins
    n = hist3(dat, [123, 20]); % default is to 10x10 bins
    n1 = n';
    n1(size(n,2) + 1, size(n,1) + 1) = 0;
    xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
    yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,2)+1);
    h = pcolor(xb,yb,n1);
    caxis([0, 1500])
    colormap(jet(15));

    shading flat
    """

    plt.tight_layout()
    plt.savefig('RESULTS/check_selection.png')

# ###################### unique_rows #############################


def unique_rows(data):
    """
    Make a unique array based on the array's row
    """
    data_mid = np.ascontiguousarray(data).view(
        np.dtype((np.void, data.dtype.itemsize * data.shape[1])))
    _, idx = np.unique(data_mid, return_index=True)
    unique_data = data[idx]
    return unique_data

# ###################### plot_sta_ev_unique #############################


def plot_sta_ev_unique(sta_info_uniq, sta_info_str_arr, ev_info_uniq,
                       ev_info_str_arr, phase, input_file_name_part,
                       req_band, num_color_grp=11):
    """
    Plot all stations and events which have been passed the criteria
    The criteria for checking the numbers is only based on lat, lon and
    el/depth. This means that we can have 4 repetitions for only two events
    because of the stations channel
    :param sta_info_uniq:
    :param sta_info_str_arr:
    :param ev_info_uniq:
    :param ev_info_str_arr:
    :param phase:
    :param input_file_name_part:
    :param req_band:
    :param num_color_grp:
    :return:
    """
    print ('Plotting events and stations...')
    plt.ion()
    plt.figure()
    plt.subplot2grid((4, 4), (0, 0), colspan=3, rowspan=3)

    m = Basemap(projection='robin', lon_0=0.0, lat_0=0.0, resolution='c')
    m.fillcontinents()
    m.drawparallels(np.arange(-90., 120., 30.))
    m.drawmeridians(np.arange(0., 420., 60.))
    m.drawmapboundary()

    # STATIONS
    col_values = range(num_color_grp+4)
    cm_cmap = plt.get_cmap('Greens')
    cNorm_cmap = colors.Normalize(vmin=0, vmax=col_values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm_cmap, cmap=cm_cmap)

    x_ev, y_ev = m(-360, 0)
    for i in range(num_color_grp-1):
        m.scatter(x_ev, y_ev, 100, color=scalarMap.to_rgba(col_values[i+4]),
                  marker="v", edgecolor="none", zorder=0,
                  label='%s-%s' % (100*i, 100*(i+1)))
    m.scatter(x_ev, y_ev, 100, color=scalarMap.to_rgba(col_values[14]),
              marker="v", edgecolor="none", zorder=0,
              label='>=%s' % ((num_color_grp-1)*100))

    sta_fio = open(os.path.join(
        os.path.curdir, 'RESULTS/%s_staloc_%s_%s.txt'
                        % (phase, input_file_name_part, req_band)), 'w')
    num_unique_sta = []
    for i in range(len(sta_info_uniq)):
        try:
            x, y = m(float(sta_info_uniq[i, 1]), float(sta_info_uniq[i, 0]))
            sta_lat_lon_str = '%s_%s_%s' % (round(sta_info_uniq[i, 0], 3),
                                            round(sta_info_uniq[i, 1], 3),
                                            round(sta_info_uniq[i, 2], 3))
            sta_size = len(sta_info_str_arr[sta_info_str_arr[:] ==
                                            sta_lat_lon_str])
            num_unique_sta.append(sta_size)
            sta_fio.writelines('%s   %s   %s \n' % (sta_info_uniq[i, 0],
                                                    sta_info_uniq[i, 1],
                                                    sta_info_uniq[i, 2]))
            if 0 <= sta_size < 100:
                sta_color = scalarMap.to_rgba(col_values[4])
            elif 100 <= sta_size < 200:
                sta_color = scalarMap.to_rgba(col_values[5])
            elif 200 <= sta_size < 300:
                sta_color = scalarMap.to_rgba(col_values[6])
            elif 300 <= sta_size < 400:
                sta_color = scalarMap.to_rgba(col_values[7])
            elif 400 <= sta_size < 500:
                sta_color = scalarMap.to_rgba(col_values[8])
            elif 500 <= sta_size < 600:
                sta_color = scalarMap.to_rgba(col_values[9])
            elif 600 <= sta_size < 700:
                sta_color = scalarMap.to_rgba(col_values[10])
            elif 700 <= sta_size < 800:
                sta_color = scalarMap.to_rgba(col_values[11])
            elif 800 <= sta_size < 900:
                sta_color = scalarMap.to_rgba(col_values[12])
            elif 900 <= sta_size < 1000:
                sta_color = scalarMap.to_rgba(col_values[13])
            else:
                sta_color = scalarMap.to_rgba(col_values[14])
            m.scatter(x, y, c=sta_color, edgecolor='none', zorder=40,
                      marker='v', s=100)
        except Exception as e:
            print ('[exception] in station %s: %s' % (i, e))
    sta_fio.close()

    # EVENTS
    col_values = range(num_color_grp+4)
    cm_cmap = plt.get_cmap('hot_r')
    cNorm_cmap = colors.Normalize(vmin=0, vmax=col_values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm_cmap, cmap=cm_cmap)

    x_ev, y_ev = m(-360, 0)
    for i in range(num_color_grp-1):
        m.scatter(x_ev, y_ev, 100, color=scalarMap.to_rgba(col_values[i+4]),
                  marker="o", edgecolor="none", zorder=0,
                  label='%s-%s' % (100*i, 100*(i+1)))
    m.scatter(x_ev, y_ev, 100, color=scalarMap.to_rgba(col_values[14]),
              marker="o", edgecolor="none", zorder=0,
              label='>=%s' % ((num_color_grp-1)*100))

    plt.legend(bbox_to_anchor=(1.2, 1.0), loc=2, borderaxespad=0., ncol=2,
               mode='expand', title='   Event        Station', fontsize=18)
    plt.gca().get_legend().get_title().set_fontsize(24)
    plt.gca().get_legend().get_title().set_fontweight('bold')

    ev_fio = open(os.path.join(
        os.path.curdir, 'RESULTS/%s_evloc_%s_%s.txt'
                        % (phase, input_file_name_part, req_band)), 'w')
    num_unique_ev = []
    for i in range(len(ev_info_uniq)):
        try:
            x, y = m(float(ev_info_uniq[i, 1]), float(ev_info_uniq[i, 0]))
            ev_lat_lon_str = '%s_%s_%s' % (round(ev_info_uniq[i, 0], 3),
                                           round(ev_info_uniq[i, 1], 3),
                                           round(ev_info_uniq[i, 2], 3))
            ev_size = len(ev_info_str_arr[ev_info_str_arr[:] ==
                                          ev_lat_lon_str])
            num_unique_ev.append(ev_size)
            ev_fio.writelines('%s   %s   %s\n' % (ev_info_uniq[i, 0],
                                                  ev_info_uniq[i, 1],
                                                  ev_info_uniq[i, 2]))
            if 0 <= ev_size < 100:
                ev_color = scalarMap.to_rgba(col_values[4])
            elif 100 <= ev_size < 200:
                ev_color = scalarMap.to_rgba(col_values[5])
            elif 200 <= ev_size < 300:
                ev_color = scalarMap.to_rgba(col_values[6])
            elif 300 <= ev_size < 400:
                ev_color = scalarMap.to_rgba(col_values[7])
            elif 400 <= ev_size < 500:
                ev_color = scalarMap.to_rgba(col_values[8])
            elif 500 <= ev_size < 600:
                ev_color = scalarMap.to_rgba(col_values[9])
            elif 600 <= ev_size < 700:
                ev_color = scalarMap.to_rgba(col_values[10])
            elif 700 <= ev_size < 800:
                ev_color = scalarMap.to_rgba(col_values[11])
            elif 800 <= ev_size < 900:
                ev_color = scalarMap.to_rgba(col_values[12])
            elif 900 <= ev_size < 1000:
                ev_color = scalarMap.to_rgba(col_values[13])
            else:
                ev_color = scalarMap.to_rgba(col_values[14])
            m.scatter(x, y, c=ev_color, edgecolor='none', zorder=40,
                      marker='o', s=100)
        except Exception as e:
            print ('[exception] in event %s: %s' % (i, e))
    ev_fio.close()

    plt.figure()
    plt.subplot(2, 1, 1)
    num_unique_sta.sort()
    plt.plot(num_unique_sta, 'r', lw=3)
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.xlabel('station number', size='large', weight='bold')
    plt.ylabel('# of repetition', size='large', weight='bold')
    plt.title('Total number of src-rcv pairs: %s' % (sum(num_unique_sta)),
              size='x-large', weight='bold')

    plt.subplot(2, 1, 2)
    num_unique_ev.sort()
    plt.plot(num_unique_ev, 'b', lw=3)
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.xlabel('event number', size='large', weight='bold')
    plt.ylabel('# of repetition', size='large', weight='bold')
    plt.title('Total number of src-rcv pairs: %s' % (sum(num_unique_ev)),
              size='x-large', weight='bold')

    plt.tight_layout()
    plt.savefig('RESULTS/sta_ev_unique.png')

# ###################### write_staev_locfile #############################


def write_staev_locfile(filt_array, inp, req_band):
    """
    Write station and event location files to be used for other calculations
    such as inversion mesh generation
    :param filt_array:
    :param inp:
    :param req_band:
    :return:
    """
    # Keep lat, lon, ele(depth) for stations (events)
    sta_info = np.array([filt_array[:, 2].astype(np.float),
                         filt_array[:, 3].astype(np.float),
                         filt_array[:, 22].astype(np.float)])
    ev_info = np.array([filt_array[:, 27].astype(np.float),
                        filt_array[:, 28].astype(np.float),
                        filt_array[:, 29].astype(np.float)])
    # Each row relates to one event-station pair
    sta_info = np.transpose(sta_info)
    ev_info = np.transpose(ev_info)

    sta_info_str = []
    for i in range(len(sta_info)):
        sta_info_str.append('%s_%s_%s'
                            % (round(sta_info[i, 0], 3),
                               round(sta_info[i, 1], 3),
                               round(sta_info[i, 2], 3)))
    ev_info_str = []
    for i in range(len(ev_info)):
        ev_info_str.append('%s_%s_%s'
                           % (round(ev_info[i, 0], 3),
                              round(ev_info[i, 1], 3),
                              round(ev_info[i, 2], 3)))

    sta_info_str_arr = np.array(sta_info_str)
    ev_info_str_arr = np.array(ev_info_str)

    sta_info_uniq = unique_rows(sta_info)
    ev_info_uniq = unique_rows(ev_info)

    print('\n========================')
    print('#unique stations: %s' % len(sta_info_uniq))
    print('#unique events  : %s' % len(ev_info_uniq))
    print('========================')

    plot_sta_ev_unique(sta_info_uniq, sta_info_str_arr, ev_info_uniq,
                       ev_info_str_arr, inp.phase, inp.input_file_name_part,
                       req_band)


# ###################### axi_kernel_receiver_writer ###########################


def axi_kernel_receiver_writer(evstas, inp):
    """
    Create AXISEM receiver.dat file for further analysis using AXISEM
    """
    print('\n======>> Create AXISEM receiver.dat file for Kernel calculations')
    if not os.path.isdir(os.path.join(inp.out_path, evstas[0, 23, 0])):
        os.mkdir(os.path.join(inp.out_path, evstas[0, 23, 0]))
    else:
        print('Directory already exists: %s' % os.path.join(inp.out_path, evstas[0, 23, 0]))
        return

    band_period = {'band01': 30.0, 'band02': 21.2, 'band03': 15.0,
                   'band04': 10.6, 'band05': 7.5, 'band06': 5.3,
                   'band07': 3.7, 'band08': 2.7}

    line_write = []
    counter = 0
    for j in range(evstas.shape[0]):
        len_freq_avail = len(evstas[j, :, evstas[j, 10, :] != -1])
        if len_freq_avail == 0:
            continue
        third_line = '%s_%s  %s  %s  %s\n' % (evstas[j, 5, 0].split('.')[0],
                                              evstas[j, 5, 0].split('.')[1],
                                              evstas[j, 2, 0], evstas[j, 3, 0],
                                              len_freq_avail)
        line_write.append(third_line)
        for i in range(len(evstas[j, 10, :])):
            if evstas[j, 10, i] == -1:
                continue
            period = band_period[evstas[j, 30, i]]
            forth_line = 'kernel_%s  Gabor_%s  CC  %s  %s\n' \
                         % (str(period).replace('.', '_'), period,
                            evstas[j, 17, i],
                            float(evstas[j, 17, i])+float(evstas[j, 18, i]))
            line_write.append(forth_line)
        counter += 1
    first_line = '%i\n' % counter
    second_line = 'vp  Z\n'
    receiver_fio = open(os.path.join(inp.out_path,
                                     evstas[0, 23, 0], 'receiver.dat'), 'w')
    receiver_fio.writelines(first_line)
    receiver_fio.writelines(second_line)
    for ln in line_write:
        receiver_fio.writelines(ln)
    receiver_fio.close()

    readme_fio = open(
        os.path.join(inp.out_path,
                     evstas[0, 23, 0], 'README.txt'), 'w')
    readme_fio.writelines('number_of_receivers: %i\n' % counter)
    readme_fio.close()

    #plt.ion()
    #plt.figure()
    #plt.subplot(2, 1, 1)
    #plt.plot(filt_array[:, 2].astype(np.float), t_corr, 'r.')
    #plt.xlabel('Latitude', size='large', weight='bold')
    #plt.ylabel('Common Correction Values', size='large', weight='bold')
    #plt.xticks(size='large', weight='bold')
    #plt.yticks(size='large', weight='bold')
    #plt.subplot(2, 1, 2)
    #plt.plot(filt_array[:, 3].astype(np.float), t_corr, 'r.')
    #plt.xlabel('Longitude', size='large', weight='bold')
    #plt.ylabel('Common Correction Values', size='large', weight='bold')
    #plt.xticks(size='large', weight='bold')
    #plt.yticks(size='large', weight='bold')
    #plt.show

# ###################### axi_kernel_receiver_combiner #########################


def axi_kernel_receiver_combiner(measure_par_add, measures):
    """
    For several type of measurements (P, Pdiff, ...), combine the receiver.dat file
    Usage:
    from output_reader import axi_kernel_receiver_combiner
    axi_kernel_receiver_combiner('./RESULTS', ['P', 'Pdiff'])
    """
    ls_all_events = np.array([], dtype='object')
    for i in range(len(measures)):
        ls_all_events = np.append(ls_all_events,
                                  glob.glob(os.path.join(measure_par_add,
                                                         measures[i],
                                                         '*.*.*.*')))
    for i in range(len(ls_all_events)):
        ls_all_events[i] = os.path.basename(ls_all_events[i])
    ls_all_events_unique = np.unique(ls_all_events)
    for i in range(len(ls_all_events_unique)):
        cont_flag = True
        for j in range(len(measures)):
            if not os.path.isfile(os.path.join(measure_par_add, measures[j],
                                               ls_all_events_unique[i],
                                               'receiver.dat')):
                cont_flag = False
        if not cont_flag:
            continue
        all_staevs = []
        len_staev = 0
        for j in range(len(measures)):
            receiver_tmp_fio = open(os.path.join(measure_par_add, measures[j],
                                                 ls_all_events_unique[i],
                                                 'receiver.dat'), 'r')
            receiver_tmp = receiver_tmp_fio.readlines()
            all_staevs.append(receiver_tmp[2:])

            len_staev_fio = open(os.path.join(measure_par_add, measures[j],
                                              ls_all_events_unique[i],
                                              'README.txt'), 'r')
            len_staev += int(len_staev_fio.readlines()[0].split(':')[1])
        receiver_all_fio = open(
            os.path.join(measure_par_add,
                         'receiver_%s.dat'
                         % ls_all_events_unique[i].replace('.', '_')), 'w')
        receiver_all_fio.writelines('%i\n' % len_staev)
        receiver_all_fio.writelines('vp   Z\n')
        for ln_all in all_staevs:
            receiver_all_fio.writelines(ln_all)
        receiver_all_fio.close()

    print('\n=================================')
    print('WARNING: vp and Z are hard coded!')
    print('=================================')

# ###################### compile_raydata_raymatrix ############################


def compile_raydata_raymatrix(inp):
    """
    Compile both raydata and raymatrix for further usage
    :return:
    """
    cur_dir = os.path.abspath(os.curdir)
    os.chdir(os.path.join(os.curdir, 'src', 'raydata_src'))
    
    if inp.crust_corr == 'crust1':
        cprint('utils.py', '[compile_raydata_raymatrix]', bc.bmagenta, 'CRUST1.0')
        os_sys = os.system('./make')
        if not os_sys == 0:
            cprint('utils.py', '[compile_raydata_raymatrix]',
                   bc.red, 'ERROR: make cannot be compiled!')
            sys.exit()
    
    elif inp.crust_corr == 'crust2':
        cprint('utils.py', '[compile_raydata_raymatrix]',
               bc.bold + bc.green, 'CRUST2.0')
        os_sys = os.system('./make2')
        if not os_sys == 0:
            cprint('utils.py', '[compile_raydata_raymatrix]',
                   bc.red, 'ERROR: make2 cannot be compiled!')
            sys.exit()
    
    elif inp.crust_corr == 'crust1io':
        cprint('utility_codes.py', '[compile_raydata_raymatrix]', bc.bmagenta, 'CRUST1.0 + A.M Indian Ocean')
        os_sys = os.system('./makeio')
        if not os_sys == 0:
            cprint('utility_codes.py', '[compile_raydata_raymatrix]',
                   bc.red, 'ERROR: make cannot be compiled!')
            sys.exit()
    
    else:
        cprint('utils.py', '[compile_raydata_raymatrix]',
               bc.red, 'ERROR: this crust correction module is not implemented')

    if inp.run_raymatrix:
        cprint('utils.py', '[compile_raydata_raymatrix]', bc.bmagenta, 'Compiling raymatrix.f')
        os.chdir(os.path.join(cur_dir, 'src', 'raymatrix_src'))
        os_sys = os.system('./make')
        if not os_sys == 0:
            cprint('utils.py', '[compile_raydata_raymatrix]',
                   bc.red, 'ERROR: make cannot be compiled!')
            sys.exit()

    os.chdir(cur_dir)

# ###################### raydata_input_generator #############################


def raydata_input_generator(filt_array, input_file_name, req_band, inp,
                            filt_file='bpf.omega_m'):
    """
    Generate input file compatible with raydata input files
    :param filt_array:
    :param input_file_name:
    :param req_band:
    :param inp:
    :param filt_file:
    :return:
    """
    # filt_file = 'gauss_filter_51_15'
    cprint('utility_codes.py', '[raydata_input_generator]', bc.blue, 'Selected filter: %s' % filt_file)
    import ipdb; ipdb.set_trace()
    if not os.path.isfile(
            os.path.join(os.path.curdir, 'src',
                         'files', 'Pdef_%s' % inp.phase)):
        sys.exit('%s could not be found!'
                 % os.path.join(os.path.curdir, 'src',
                                'files', 'Pdef_%s' % inp.phase))
    phase_def = open(os.path.join(
        os.path.curdir, 'src',
        'files', 'Pdef_%s' % inp.phase), 'r').readlines()

    if not os.path.isfile(
            os.path.join(inp.events_dir, filt_file)):
        sys.exit('%s could not be found!'
                 % os.path.join(inp.events_dir, filt_file))

    filt_def = open(os.path.join(inp.events_dir, filt_file), 'r').readlines()

    inp_lines = []
    inp_lines.append('%s\n' % input_file_name)
    inp_lines.append('%s\n' % inp.twinned)
    inp_lines.append('# Phase : %s\n' % inp.phase)
    inp_lines.append('# Traveltime data, xcorr_min : %s\n' % inp.min_xcorr)
    inp_lines.append('# depth_min : %s, depth_max : %s\n'
                     % (inp.min_depth, inp.max_depth))
    inp_lines.append('# epi_min : %s, epi_max : %s\n'
                     % (inp.min_epi, inp.max_epi))
    inp_lines.append('# Clip : %s\n' % inp.check_clip)
    for ln in phase_def:
        inp_lines.append(ln)
    for ln in filt_def:
        inp_lines.append(ln)

    if inp.phase == 'Pdiff':
        geocen_lat_sta = geocen_array(filt_array[:, 2].astype(np.float))
        geocen_lat_evt = geocen_array(filt_array[:, 27].astype(np.float))
        if not len(geocen_lat_sta) == len(filt_array):
            sys.exit('Something is wrong with the geocentric lats')
        if not len(geocen_lat_evt) == len(filt_array):
            sys.exit('Something is wrong with the geocentric lats')
    counter = 0
    tau_bg = tau.TauPyModel(model='iasp91')
    for sta in filt_array:
        #import ipdb; ipdb.set_trace()
        # progress_bar(counter, len(filt_array))
        counter += 1
        netid = sta[5].split('.')[0]
        if not inp.plot_statistics:
            staid = sta[5].split('.')[1]
        else:
            staid = '%07i' % counter
        if not "x00" in sta:
            chaid = sta[5].split('.')[3]
        else:
            chaid = sta[5].split('.')[4]
        # import ipdb; ipdb.set_trace()
        first_line = '%s %s %s %s     %s   %s %s    %s    %s    %s    %s  %s   %s 1 0 0\n' \
                     % (sta[24], sta[25], sta[26], sta[1], staid, netid, chaid,
                        sta[27], sta[28], sta[29], sta[2], sta[3], sta[22])

        second_line = '1      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0\n'
        if inp.phase == 'Pdiff':
            tB_mfi = travel_time_calc(evla=geocen_lat_evt[counter-1],
                                      evlo=float(sta[28]),
                                      stla=geocen_lat_sta[counter-1],
                                      stlo=float(sta[3]),
                                      evdp=float(sta[29]),
                                      bg_model='iasp91',
                                      tau_bg=tau_bg,
                                      obspytool=True)
        else:
            tB_mfi = sta[17]
        if not tB_mfi:
            # print "Do not include %s.%s.%s -- dist: %s!" \
            #       % (netid, staid, chaid, sta[4])
            continue
        third_line = '    %s  %s  %s  %s   %s   %s\n' \
                     % (sta[7], sta[9], sta[6],
                        int(req_band.split('band')[1]), sta[18], tB_mfi)
        inp_lines.append(first_line)
        inp_lines.append(second_line)
        inp_lines.append(third_line)
    input_file_fio = open(os.path.join(inp.out_path, input_file_name), 'w')
    input_file_fio.writelines(inp_lines)
    input_file_fio.close()

# ###################### geocen_array #############################


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

# ###################### travel_time_calc #############################


def travel_time_calc(evla, evlo, stla, stlo, evdp, bg_model,
                     tau_bg, obspytool=True):
    """
    calculate arrival time of different seismic phases
    :param evla:
    :param evlo:
    :param stla:
    :param stlo:
    :param evdp:
    :param bg_model:
    :param tau_bg:
    :param obspytool:
    :return:
    """
    # --------------- TAUP
    if obspytool:
        dist = locations2degrees(evla, evlo, stla, stlo)
        try:
            tt = tau_bg.get_travel_times(evdp, dist, ('Pdiff',))[0].time
        except Exception as e:
            tt = False
    else:
        # ----------------------- Java method
        taup_process = subprocess.Popen(['taup_time', '-mod', bg_model, '-time',
                                         '-h', str(evdp),
                                         '-ph', 'Pdiff',
                                         '-sta', str(stla), str(stlo),
                                         '-evt', str(evla), str(evlo)],
                                        stdout=subprocess.PIPE)

        tt_raw = taup_process.communicate()[0]
        try:
            tt = tt_raw.split('\n')[0].split()[-1]
            tt = float(tt)
        except Exception as e:
            print('Requested phase: Pdiff ... ERROR: %s' % e)
            print('sta: %s %s' % (stla, stlo))
            print('evt: %s %s' % (evla, evlo))
            tt = False
    # --------------- END TAUP
    return tt

# ###################### raydata_input #############################


def raydata_input(inp, input_file_name):
    """
    make in.input_file_name for raydata
    :param inp:
    :param input_file_name:
    :return:
    """
    in_input_fio = open(os.path.join(inp.out_path, 'in.raydata_%s'
                                                % input_file_name), 'w')
    in_input_fio.write('%s\n' % inp.bg_model)
    in_input_fio.write('%s  %s\n'
                       % (inp.max_num_arrival, inp.delay_wrt_first_arrival))
    in_input_fio.write('%s\n' % input_file_name)
    in_input_fio.write('Pdef_%s' % inp.phase)
    in_input_fio.close()

# ###################### raymatrix_input #############################


def raymatrix_input(inp, input_file_name):
    """
    create in.raymatrix_input_file_name
    :param inp:
    :param input_file_name:
    :return:
    """
    in_input_fio = open(
        os.path.join(inp.out_path, 'in.raymatrix_%s' % input_file_name), 'w')
    in_input_fio.write('%s %s %s\n'
                       % (inp.vp_vs_Qs[0], inp.vp_vs_Qs[1], inp.vp_vs_Qs[2]))
    in_input_fio.write('%s\n' % inp.kernel_quad_km)
    in_input_fio.write('%s\n' % inp.vertex_file)
    in_input_fio.write('%s\n' % inp.facet_file)
    in_input_fio.write('%s' % input_file_name)
    in_input_fio.close()

    in_input_fio = open(
        os.path.join(inp.out_path, 'in.matrixT.%s' % input_file_name), 'w')
    in_input_fio.write('matrixT.%s' % input_file_name)
    in_input_fio.close()

# ###################### prepare_dir #############################


def prepare_dir(input_file_name, inp, vertices, facets):
    """
    prepare directory for one run of raydata and raymatrix
    :param input_file_name:
    :param vertices:
    :param facets:
    :return:
    """
    cprint('utility_codes.py', '[prepare_dir]', bc.blue, 'Prepare output directory at: %s/%s '
           % (inp.out_path, input_file_name))
    cur_dir = os.path.abspath(os.curdir)
    # if os.path.isdir(os.path.join(inp.out_path, '%s_dir' % input_file_name)):
    #     # answer = raw_input("Directory already exists! Remove [R] or Exit [E]?")
    #     answer = 'R'
    #     if answer == 'R':
    #         shutil.rmtree(os.path.join(inp.out_path, '%s_dir' % input_file_name))
    #     else:
    #         sys.exit('Elvis left the building...')

    os.mkdir(os.path.join(inp.out_path, '%s_dir' % input_file_name))

    shutil.move(os.path.join(inp.out_path, input_file_name),
                os.path.join(inp.out_path, '%s_dir' % input_file_name))

    shutil.move(os.path.join(inp.out_path, 'in.raydata_%s' % input_file_name),
                os.path.join(inp.out_path, '%s_dir' % input_file_name))

    shutil.move(os.path.join(inp.out_path, 'in.raymatrix_%s' % input_file_name),
                os.path.join(inp.out_path, '%s_dir' % input_file_name))

    shutil.move(os.path.join(inp.out_path, 'in.matrixT.%s' % input_file_name),
                os.path.join(inp.out_path, '%s_dir' % input_file_name))

    if inp.crust_corr == 'crust1':
        list_files = ['CNelevatio1.txt', 'CNtype1_key.txt', 'CNtype1.txt',
                      'elcordir.tbl', 'IASP91.PREMQ', 'crust1.vp', 'crust1.vs',
                      'crust1.bnds', 'crust1.rho', vertices, facets]

    elif inp.crust_corr == 'crust1io':
        list_files = ['CNelevatio1.txt', 'CNtype1_key.txt', 'CNtype1.txt',
                      'elcordir.tbl', 'IASP91.PREMQ', 'crust1.vp', 'crust1.vs',
                      'crust1.bnds', 'crust1.rho', 'AMcrustIndian.txt',
                      'AMmohoIndian.txt', 'AMsedIndian.txt', vertices, facets]

    elif inp.crust_corr == 'crust2':
        list_files = ['CNelevatio2.txt', 'CNtype2_key.txt', 'CNtype2.txt',
                      'elcordir.tbl', 'IASP91.PREMQ', vertices, facets]

    else:
        cprint('utility_codes.py', '[prepare_dir]',
               bc.red, 'ERROR: this crust correction module is not implemented')

    files_glob = os.path.join(os.path.curdir, 'src', 'files')
    for ls_fi in list_files:
        fi = os.path.join(files_glob, ls_fi)
        shutil.copy(fi, os.path.join(inp.out_path, '%s_dir' % input_file_name))

    if inp.crust_corr == 'crust1':
        shutil.copy(os.path.join(os.curdir, 'src',
                    'raydata_src', 'raydata'), os.path.join(inp.out_path, '%s_dir' % input_file_name))
    elif inp.crust_corr == 'crust1io':
        shutil.copy(os.path.join(os.curdir, 'src',
                    'raydata_src', 'raydataio'), os.path.join(inp.out_path, '%s_dir' % input_file_name))
    elif inp.crust_corr == 'crust2':
        shutil.copy(os.path.join(os.curdir, 'src',
                    'raydata_src', 'raydata2'), os.path.join(inp.out_path, '%s_dir' % input_file_name))
    else:
        cprint('utility_codes.py', '[prepare_dir]',
               bc.red, 'ERROR: this crust correction module is not implemented')

    shutil.copy(os.path.join(os.curdir, 'src',
                             'raymatrix_src', 'raymatrix'),
                os.path.join(inp.out_path, '%s_dir' % input_file_name))
    shutil.copy(os.path.join(os.curdir, 'src',
                             'raymatrix_src', 'mat2asc'),
                os.path.join(inp.out_path, '%s_dir' % input_file_name))

# ###################### run_raydata_raymatrix #############################


def run_raydata_raymatrix(input_file_name, inp, raydata=True, raymatrix=True):
    """
    run both raydata and raymatrix in the directory
    :param input_file_name:
    :param raydata:
    :param raymatrix:
    :return:
    """
    cur_dir = os.path.abspath(os.curdir)

    if raydata:
        if inp.crust_corr == 'crust1':
            cprint('utility_codes.py', '[run_raydata_raymatrix]', bc.bmagenta, 'CRUST1.0')
            cprint ('utility_codes.py', '[run_raydata_raymatrix]', bc.green,
                    'run raydata at %s/%s_dir' % (inp.out_path, input_file_name))
            os.chdir(os.path.join(inp.out_path, '%s_dir' % input_file_name))
            os_sys = os.system('./raydata < in.raydata_%s' % input_file_name)
            if not os_sys == 0:
                cprint('utility_codes.py', '[run_raydata_raymatrix]', bc.red,
                       '[ERROR] raydata was not executed correctly!')

        elif inp.crust_corr == 'crust1io':
            cprint('utility_codes.py', '[run_raydata_raymatrix]', bc.bmagenta, 'CRUST1.0*')
            cprint('utility_codes.py', '[run_raydata_raymatrix]', bc.green,
                   'run raydataio at %s/%s_dir' % (inp.out_path, input_file_name))
            os.chdir(os.path.join(inp.out_path, '%s_dir' % input_file_name))
            os_sys = os.system('./raydataio < in.raydata_%s' % input_file_name)
            if not os_sys == 0:
                cprint('utility_codes.py', '[run_raydata_raymatrix]', bc.red,
                       '[ERROR] raydata was not executed correctly!')

        elif inp.crust_corr == 'crust2':
            cprint('utility_codes.py', '[run_raydata_raymatrix]', bc.magenta + bc.bold, 'CRUST2.0')
            cprint('utility_codes.py', '[run_raydata_raymatrix]', bc.green,
                   'run raydata2 at %s/%s_dir' % (inp.out_path, input_file_name))
            os.chdir(os.path.join(inp.out_path, '%s_dir' % input_file_name))
            os_sys = os.system('./raydata2 < in.raydata_%s' % input_file_name)
            if not os_sys == 0:
                cprint('utility_codes.py', '[run_raydata_raymatrix]', bc.red,
                       '[ERROR] raydata2 was not executed correctly!')
        else:
            cprint('utility_codes.py', '[run_raydata_raymatrix]',
                   bc.red, 'ERROR: this crust correction module is not implemented')

    if raymatrix:
        cprint('utility_codes.py', '[run_raydata_raymatrix]', bc.green,
               'run raymatrix at %s/%s_dir' % (inp.out_path, input_file_name))
        os.chdir(os.path.join(inp.out_path, '%s_dir' % input_file_name))
        os_sys = os.system('./raymatrix < in.raymatrix_%s' % input_file_name)
        if not os_sys == 0:
            cprint('utility_codes.py', '[run_raydata_raymatrix]', bc.red,
                   '[ERROR] raymatrix was not executed correctly!')

    os.chdir(cur_dir)

# ###################### parallel_raydata_raymatrix ###########################


def parallel_raydata_raymatrix(filt_array, input_file_name, req_band, inp):
    """
    To run raydata and raymatrix in parallel
    :param filt_array:
    :param input_file_name:
    :param req_band:
    :param inp:
    :return:
    """


    raydata_input_generator(filt_array=filt_array,
                            input_file_name=input_file_name,
                            req_band=req_band,
                            inp=inp)
    raydata_input(inp=inp, input_file_name=input_file_name)
    raymatrix_input(inp=inp, input_file_name=input_file_name)

    # print '\n======>> prepare output directory at: %s/%s_dir' \
          # % (inp.out_path, input_file_name)
    prepare_dir(input_file_name=input_file_name, inp=inp,
                vertices=inp.vertex_file,
                facets=inp.facet_file)
    run_raydata_raymatrix(input_file_name=input_file_name, inp=inp,
                          raydata=inp.run_raydata,
                          raymatrix=inp.run_raymatrix)

# ###################### raydata_ccorr_reader #########################


def raydata_ccorr_reader(filt_array, input_file_name, corr_io_list,
                         events_dir):
    """
    Read common correction results and rewrite the values in filt_array
    :param filt_array:
    :param input_file_name:
    :param corr_io_list:
    :param events_dir:
    :return:
    """
    # kd (0), station-code (1), network-code (2), point-lat (3),
    # point-lon (4), radius (5), source-lat (6), source-lon (7),
    # ellipticity (8), crust (9), elevation (10), tstar (11),
    # distance (12), trtime_final (13), tobs (14), source-depth (15)
    from datetime import datetime
    t1 = datetime.now()
    ccorr_arr = np.loadtxt(os.path.join(os.path.curdir,
                                        inp.out_path,
                                        '%s_dir' % input_file_name,
                                        'info_collector.%s' % input_file_name),
                           dtype='S', comments='#', delimiter=',')
    # removing those stations that could not be trated in raydata
    grp_ccorr_arr_unique = np.unique(ccorr_arr[:, 1].astype(np.int))
    indx = []
    for sta_idx in np.arange(1, len(filt_array) + 1):
        # progress_bar(sta_idx, (len(filt_array)+1))
        if sta_idx not in grp_ccorr_arr_unique:
            indx.append(False)
        else:
            indx.append(True)
    filt_array = filt_array[np.array(indx)]
    print('Time to load and filter filt_array %s' % (datetime.now() - t1))

    t1 = datetime.now()
    # all the ray segments present in the definition of the ray
    kd_avail = np.unique(ccorr_arr[:, 0].astype(np.int))
    len_ccorr_arr_0 = len(ccorr_arr[ccorr_arr[:, 0].astype(np.int)
                                    == kd_avail[0]])
    t_corr = np.zeros(len_ccorr_arr_0)
    t_ellip = np.zeros(len_ccorr_arr_0)
    t_cc = np.zeros(len_ccorr_arr_0)
    t_elev = np.zeros(len_ccorr_arr_0)

    # loop over ray segments to collect the information
    for i in range(len(kd_avail)):
        progress_bar(i, len(kd_avail))
        tar_ccorr_arr = ccorr_arr[ccorr_arr[:, 0].astype(np.int) ==
                                  kd_avail[i]]
        # check how many repetition we have for the ray segment in ONE Pdef
        if len(tar_ccorr_arr)/len(t_corr) != 1:
            if len(tar_ccorr_arr) % len(t_corr) != 0:
                sys.exit('Number of ray segments is not set correctly for '
                         'all stations!')
            divi = len(tar_ccorr_arr)/len(t_corr)
        else:
            divi = False
        if not divi:
            if kd_avail[i] == 1:
                # collect ellipticity only based on the ray_seg==1
                if corr_io_list[0] == 1:
                    t_ellip += tar_ccorr_arr[:, 8].astype(np.float)
            if corr_io_list[1] == 1:
                t_cc += tar_ccorr_arr[:, 9].astype(np.float)
            if corr_io_list[2] == 1:
                t_elev += tar_ccorr_arr[:, 10].astype(np.float)
        else:
            if kd_avail[i] == 1:
                sys.exit('This case should not happen (common corrections!)')
            if corr_io_list[1] == 1:
                for j in range(0, divi):
                    t_cc += tar_ccorr_arr[j:len(tar_ccorr_arr):divi, 9].\
                        astype(np.float)
            if corr_io_list[2] == 1:
                for j in range(0, divi):
                    t_elev += tar_ccorr_arr[j:len(tar_ccorr_arr):divi, 10]. \
                        astype(np.float)
    print('Time to extract correction values %s' % (datetime.now() - t1))

    sub_ccorr_arr = ccorr_arr[ccorr_arr[:, 0].astype(np.int) == kd_avail[0]]
    t_corr = t_ellip + t_cc + t_elev
    t_final = sub_ccorr_arr[:, 14].astype(np.float) - \
              (sub_ccorr_arr[:,  13].astype(np.float) + t_corr)
    new_col_cr = np.empty([np.shape(filt_array)[0], 7], dtype=object)
    new_col_cr[:, 0] = t_ellip
    new_col_cr[:, 1] = t_cc
    new_col_cr[:, 2] = t_elev
    new_col_cr[:, 3] = t_corr
    new_col_cr[:, 4] = t_final
    new_col_cr[:, 5] = sub_ccorr_arr[:, 13]
    new_col_cr[:, 6] = sub_ccorr_arr[:, 14]
    filt_array = np.append(filt_array, new_col_cr, 1)
    # writing the filt_array
    filt_array_to_write = np.copy(filt_array)
    filt_array_to_write[:, 8] = filt_array_to_write[:, 35]
    raydata_ccorr_writer(filt_array_to_write, events_dir)

    plt.ion()
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.scatter(filt_array[:, 2].astype(np.float), filt_array[:, 34],
                c='r', edgecolors='none', alpha=alpha_plt)
    plt.xlabel('Latitude', size='large', weight='bold')
    plt.ylabel('Common Correction Values', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.subplot(2, 1, 2)
    plt.scatter(filt_array[:, 3].astype(np.float), filt_array[:, 34],
                c='b', edgecolors='none', alpha=alpha_plt)
    plt.xlabel('Longitude', size='large', weight='bold')
    plt.ylabel('Common Correction Values', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.savefig('RESULTS/lat_lon_corr.png')

    plt.ion()
    plt.figure()
    plt.subplot(2, 3, 1)
    plt.scatter(filt_array[:, 2].astype(np.float), filt_array[:, 31], c='r',
                edgecolors='none', alpha=alpha_plt)
    plt.xlabel('Latitude', size='large', weight='bold')
    plt.ylabel('Ellipticity Correction', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.subplot(2, 3, 2)
    plt.scatter(filt_array[:, 2].astype(np.float), filt_array[:, 32], c='r',
                edgecolors='none', alpha=alpha_plt)
    plt.xlabel('Latitude', size='large', weight='bold')
    plt.ylabel('Crustal Correction', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.subplot(2, 3, 3)
    plt.scatter(filt_array[:, 2].astype(np.float), filt_array[:, 33], c='r',
                edgecolors='none', alpha=alpha_plt)
    plt.xlabel('Latitude', size='large', weight='bold')
    plt.ylabel('Topography Correction', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')

    plt.subplot(2, 3, 4)
    plt.scatter(filt_array[:, 3].astype(np.float), filt_array[:, 31], c='b',
                edgecolors='none', alpha=alpha_plt)
    plt.xlabel('Longitude', size='large', weight='bold')
    plt.ylabel('Ellipticity Correction', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.subplot(2, 3, 5)
    plt.scatter(filt_array[:, 3].astype(np.float), filt_array[:, 32], c='b',
                edgecolors='b', alpha=alpha_plt)
    plt.xlabel('Longitude', size='large', weight='bold')
    plt.ylabel('Crustal Correction', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.subplot(2, 3, 6)
    plt.scatter(filt_array[:, 3].astype(np.float), filt_array[:, 33], c='b',
                edgecolors='b', alpha=alpha_plt)
    plt.xlabel('Longitude', size='large', weight='bold')
    plt.ylabel('Topography Correction', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.savefig('RESULTS/lat_lon_each_corr.png')

    plt.figure()
    plt.scatter(filt_array[:, 2].astype(np.float), filt_array[:, 35], c='b',
                edgecolors='none', alpha=alpha_plt)
    plt.xlabel('latitude', size='large', weight='bold')
    plt.ylabel('Corrected Measurements', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.savefig('RESULTS/lat_measurement.png')

    return filt_array

# ###################### raydata_ccorr_writer #############################


def raydata_ccorr_writer(filt_array_corr, events_dir):
    """
    Write similar files as ffproc.ampstt.* but with common correction applied
    :param filt_array_corr:
    :param events_dir:
    :return:
    """
    req_dirs = np.unique(filt_array_corr[:, 23])
    for r_dir in req_dirs:
        if not os.path.isdir(os.path.join(os.path.curdir, inp.out_path, r_dir)):
            os.mkdir(os.path.join(os.path.curdir, inp.out_path, r_dir))
        if not os.path.isdir(os.path.join(os.path.curdir, inp.out_path, r_dir,
                                          'outfiles')):
            os.mkdir(os.path.join(os.path.curdir, inp.out_path, r_dir,
                                  'outfiles'))
        filt_array_this_dir = filt_array_corr[filt_array_corr[:, 23] == r_dir]
        np.savetxt(os.path.join(os.path.curdir, inp.out_path, r_dir, 'outfiles',
                                'ffproc.ampstt.%s'
                                % filt_array_this_dir[0, 30]),
                   filt_array_this_dir[:, 0:22],
                   fmt='%s', delimiter='     ')
        shutil.copy(os.path.join(events_dir, r_dir,
                                 'outfiles', 'ffproc.source'),
                    os.path.join(os.curdir, inp.out_path, r_dir, 'outfiles'))
        shutil.copy(os.path.join(events_dir, r_dir,
                                 'outfiles', 'ffproc.receivers'),
                    os.path.join(os.curdir, inp.out_path, r_dir, 'outfiles'))
        shutil.copy(os.path.join(events_dir, r_dir,
                                 'outfiles', 'ampinv.source'),
                    os.path.join(os.curdir, inp.out_path, r_dir, 'outfiles'))

# ###################### check_par_jobs #############################


def check_par_jobs(jobs, sleep_time=1):
    """
    check whether all the parallel jobs are finished or not
    """
    pp_flag = True
    while pp_flag:
        for proc in jobs:
            if proc.is_alive():
                print('.',)
                # sys.stdout.flush()
                time.sleep(sleep_time)
                pp_flag = True
                break
            else:
                pp_flag = False
    if not pp_flag:
        print ('\n\nAll %s processes are finished...\n' % len(jobs))

# ###################### mat2asc_run #############################


def mat2asc_run(input_file_name, inp):
    """
    run mat2asc in each directory
    :param input_file_name:
    :return:
    """
    # import ipdb; ipdb.set_trace()
    cur_dir = os.path.abspath(os.curdir)
    cprint('utility_codes.py', '[mat2asc_run]', bc.green, 'run mat2asc at %s/%s_dir' % (inp.out_path, input_file_name))
    os.chdir(os.path.join(inp.out_path, '%s_dir' % input_file_name))
    os_sys = os.system('./mat2asc < in.matrixT.%s' % input_file_name)
    if not os_sys == 0:
        cprint('utility_codes.py', '[mat2asc_run]', bc.red, '[ERROR] mat2asc was not executed correctly!')
    os.chdir(cur_dir)

# ###################### vtk_generator #############################


def vtk_generator(inp, req_band, len_dirs):
    """
    VTK file generator out of all the results for one complete run
    INFO:
    - facet file is indexed from 0
    - BUT! ascii.* files are indexed from 1, so we need to subtract them (-1)
    - When describing the cells in terms of point indices, the points must be
    indexed starting at 0.
    """
    # Only plot one kernel:
    # import ipdb; ipdb.set_trace()
    single_kernel = False

    cprint('utility_codes.py', '[VTK GENERATOR]', bc.blue, 'Creating VTK file')
    input_file_name = '%s_%s_%s' % (inp.input_file_name_part, req_band, 1)
    direname = os.path.join(inp.out_path, '%s_dir' % input_file_name)

    cprint('utility_codes.py', '[VTK GENERATOR]', bc.blue, 'Load Vertex file')
    mesh_points = np.loadtxt(os.path.join(direname, inp.vertex_file),
                             skiprows=2, comments='#')
    cprint('utility_codes.py', '[VTK GENERATOR]', bc.blue, 'Load Facet file')
    mesh_facets = np.loadtxt(os.path.join(direname, inp.facet_file),
                             dtype=np.int, skiprows=1, comments='#')

    mat_val_all = [0.0]*len(mesh_points)
    # import ipdb; ipdb.set_trace()
    for nj in range(len_dirs):
        input_file_name = \
            '%s_%s_%s' % (inp.input_file_name_part, req_band, nj+1)
        direname = os.path.join(
            os.path.curdir, inp.out_path, '%s_dir' % input_file_name)
        cprint('utility_codes.py', '[VTK GENERATOR]', bc.blue,
               'Create VTK file at %s/%s_dir' % (inp.out_path, input_file_name))
        ascii_file = 'ascii.matrixT.%s' % input_file_name
        fmatrix = open(os.path.join(direname, ascii_file), 'r')
        fmatrix_r = fmatrix.readlines()
        mat_indx = []
        mat_val = []
        counter = 0
        for j in range(5, len(fmatrix_r)):
            if counter == 3:
                counter = 0
                continue
            if counter == 0:
                # mat_indx: matrix index
                mat_indx_tmp = fmatrix_r[j].split()
                for i in range(len(mat_indx_tmp)):
                    # IF indexing in ASCII file starts from 0,
                    # we do not need -1.
                    # Otherwise we need it!
                    mat_indx.append(int(mat_indx_tmp[i])-1)
            if counter == 1:
                # mat_val: matrix value
                mat_val_tmp = fmatrix_r[j].split()
                for i in range(len(mat_val_tmp)):
                    if inp.abs_vtk:
                        mat_val.append(abs(float(mat_val_tmp[i])))
                    else:
                        mat_val.append(float(mat_val_tmp[i]))
            if single_kernel:
                if j == 7:
                    break
            counter += 1

        if 'p' in inp.phase.lower():
            for i in range(len(mat_indx)):
                mat_val_all[mat_indx[i]] += mat_val[i]
        elif 's' in inp.phase.lower():
            for i in range(len(mat_indx)):
                mat_val_all[mat_indx[i]-len(mat_val_all)] += mat_val[i]
    cprint('utility_codes.py', '[VTK GENERATOR]', bc.blue, 'Sum over all elems: %s' % sum(mat_val_all))

    vtk = pvtk.VtkData(pvtk.UnstructuredGrid(mesh_points, tetra=mesh_facets),
                       pvtk.PointData(pvtk.Scalars(mat_val_all,
                                                   name='kernel_value')),
                       'Inversion Grid')
    print('==================')
    print("WARNING: INDEXING!")
    print('==================')
    vtk.tofile(os.path.join(
        inp.out_path, '%s_%s.vtk'
                        % (input_file_name.split('_')[0], req_band)))

    return sum(mat_val_all)

# ###################### vtk_generator_all #############################


def vtk_generator_all(direname, vertex_file, facet_file):
    """
    VTK file generator out of all the results of all runs
    ATTENTION: ascii.matrix.* should be moved to one dir (direname)
    INFO:
    - facet file is indexed from 0
    - BUT! ascii.* files are indexed from 1, so we need to subtract them (-1)
    - When describing the cells in terms of point indices, the points must be
    indexed starting at 0.
    """
    print ('\n======>> Creating VTK file')
    ascii_files = glob.glob(os.path.join(direname, 'ascii.matrixT.*'))

    print ('------> Load Vertex file')
    mesh_points = np.loadtxt(os.path.join(direname, vertex_file),
                             skiprows=2, comments='#')
    print ('------> Load Facet file')
    mesh_facets = np.loadtxt(os.path.join(direname, facet_file),
                             dtype=np.int, skiprows=1, comments='#')

    mat_val_all = [0.0]*len(mesh_points)
    for nj in range(len(ascii_files)):
        print ('\n------> create VTK file %s' % ascii_files[nj])
        fmatrix = open(os.path.join(ascii_files[nj]), 'r')
        fmatrix_r = fmatrix.readlines()
        mat_indx = []
        mat_val = []
        counter = 0
        for j in range(5, len(fmatrix_r)):
            if counter == 3:
                counter = 0
                continue
            if counter == 0:
                # mat_indx: matrix index
                mat_indx_tmp = fmatrix_r[j].split()
                for i in range(len(mat_indx_tmp)):
                    # IF indexing in ASCII file starts from 0,
                    # we do not need -1.
                    # Otherwise we need it!
                    mat_indx.append(int(mat_indx_tmp[i])-1)
            if counter == 1:
                # mat_val: matrix value
                mat_val_tmp = fmatrix_r[j].split()
                for i in range(len(mat_val_tmp)):
                    mat_val.append(abs(float(mat_val_tmp[i])))
            counter += 1

        for i in range(len(mat_indx)):
            mat_val_all[mat_indx[i]] += mat_val[i]

    vtk = pvtk.VtkData(pvtk.UnstructuredGrid(mesh_points,
                                             tetra=mesh_facets),
                       pvtk.PointData(pvtk.Scalars(mat_val_all,
                                                   name='kernel_value')),
                       'Inversion Grid')
    print ('\n\n==================')
    print ("WARNING: INDEXING!")
    print ('==================')
    vtk.tofile(os.path.join(direname, 'global_all.vtk'))

# ###################### vtk_val_azi #############################


def vtk_val_azi(direname, vertex_file, facet_file):
    """
    VTK file generator out of all the results of all runs
    This version of VTK maker generates a .vtk file that 
    contains the directions (based on azimuth).
    It can be used to see whether a cell is well illuminated
    from all directions.

    ATTENTION: ascii.matrix.* should be moved to one dir (direname)
    INFO:
    - facet file is indexed from 0
    - BUT! ascii.* files are indexed from 1, so we need to subtract them (-1)
    - When describing the cells in terms of point indices, the points must be
    indexed starting at 0.

    Example:
    from output_reader import vtk_val_azi
    direname = './RESULTS/combined'
    vertex_file = 'vertices.ppdiff_staev_200'
    facet_file = 'facets.ppdiff_staev_200'
    vtk_val_azi(direname, vertex_file, facet_file)
    """
    print ('\n======>> Creating VTK file')
    ascii_files = glob.glob(os.path.join(direname, 'ascii.matrixT.*'))

    print ('------> Load Vertex file')
    mesh_points = np.loadtxt(os.path.join(direname, vertex_file),
                             skiprows=2, comments='#')
    print ('------> Load Facet file')
    mesh_facets = np.loadtxt(os.path.join(direname, facet_file),
                             dtype=np.int, skiprows=1, comments='#')
    mesh_points_2 = mesh_points**2
    rad = np.sqrt(mesh_points_2[:,0] + mesh_points_2[:,1] + mesh_points_2[:,2])

    mat_val_all = np.array([0.0]*len(mesh_points))
    mat_azi_qual = np.zeros([len(mesh_points), 4])
    mat_azi_all = {}
    for nj in range(len(ascii_files)):
        print ('\n------> create VTK file %s' % ascii_files[nj])
        fmatrix = open(os.path.join(ascii_files[nj]), 'r')
        fmatrix_r = fmatrix.readlines()
        mat_indx = []
        mat_val = []
        mat_azi = []
        counter = 0
        print ('Length of fmatrix_r: %s' % len(fmatrix_r))
        for j in range(3, len(fmatrix_r)):
            if counter == 0:
                azi = float(fmatrix_r[j].split()[17])*180./np.pi
                #print azi
            if counter == 2:
                # mat_indx: matrix index
                mat_indx_tmp = fmatrix_r[j].split()
                for i in range(len(mat_indx_tmp)):
                    # IF indexing in ASCII file starts from 0,
                    # we do not need -1.
                    # Otherwise we need it!
                    mat_indx.append(int(mat_indx_tmp[i])-1)
            if counter == 3:
                # mat_val: matrix value
                mat_val_tmp = fmatrix_r[j].split()
                for i in range(len(mat_val_tmp)):
                    mat_val.append(abs(float(mat_val_tmp[i])))
                    mat_azi.append(azi)
                counter = -1
            counter += 1

        mat_indx_arr = np.array(mat_indx)
        mat_azi_arr = np.array(mat_azi)

        mat_1 = np.array([0.0]*len(mesh_points))
        mat_1[mat_indx_arr[(0<=mat_azi_arr) & (mat_azi_arr<45)]] = \
                np.bincount(mat_indx_arr[(0<=mat_azi_arr) & (mat_azi_arr<45)])
        mat_azi_qual[:, 0] += mat_1


        mat_1 = np.array([0.0]*len(mesh_points))
        mat_1[mat_indx_arr[(45<=mat_azi_arr) & (mat_azi_arr<90)]] = \
                np.bincount(mat_indx_arr[(45<=mat_azi_arr) & (mat_azi_arr<90)])
        mat_azi_qual[:, 1] += mat_1


        mat_1 = np.array([0.0]*len(mesh_points))
        mat_1[mat_indx_arr[(90<=mat_azi_arr) & (mat_azi_arr<135)]] = \
                np.bincount(mat_indx_arr[(90<=mat_azi_arr) & (mat_azi_arr<135)])
        mat_azi_qual[:, 2] += mat_1


        mat_1 = np.array([0.0]*len(mesh_points))
        mat_1[mat_indx_arr[(135<=mat_azi_arr) & (mat_azi_arr<180)]] = \
                np.bincount(mat_indx_arr[(135<=mat_azi_arr) & (mat_azi_arr<180)])
        mat_azi_qual[:, 3] += mat_1

        mat_1 = np.array([0.0]*len(mesh_points))
        mat_1[mat_indx_arr[(180<=mat_azi_arr) & (mat_azi_arr<225)]] = \
                np.bincount(mat_indx_arr[(180<=mat_azi_arr) & (mat_azi_arr<225)])
        mat_azi_qual[:, 0] += mat_1


        mat_1 = np.array([0.0]*len(mesh_points))
        mat_1[mat_indx_arr[(225<=mat_azi_arr) & (mat_azi_arr<270)]] = \
                np.bincount(mat_indx_arr[(225<=mat_azi_arr) & (mat_azi_arr<270)])
        mat_azi_qual[:, 1] += mat_1


        mat_1 = np.array([0.0]*len(mesh_points))
        mat_1[mat_indx_arr[(270<=mat_azi_arr) & (mat_azi_arr<315)]] = \
                np.bincount(mat_indx_arr[(270<=mat_azi_arr) & (mat_azi_arr<315)])
        mat_azi_qual[:, 2] += mat_1


        mat_1 = np.array([0.0]*len(mesh_points))
        mat_1[mat_indx_arr[(315<=mat_azi_arr) & (mat_azi_arr<360)]] = \
                np.bincount(mat_indx_arr[(315<=mat_azi_arr) & (mat_azi_arr<360)])
        mat_azi_qual[:, 3] += mat_1
        plt.ion()
        plt.figure()
        plt.subplot(2,2,2)
        plt.plot(mat_azi_qual[:,0])
        plt.subplot(2,2,4)
        plt.plot(mat_azi_qual[:,1])
        plt.subplot(2,2,3)
        plt.plot(mat_azi_qual[:,2])
        plt.subplot(2,2,1)
        plt.plot(mat_azi_qual[:,3])
        plt.show()


        #print 'Length of mat_indx: %s' % len(mat_indx)
        mat_val_all[mat_indx] += np.array(mat_val)
        #for i in range(len(mat_indx)):
        #    if i%1000 == 0:
        #        sys.stdout.write('\r')
        #        sys.stdout.write("[%-100s] %d%% ---" % ('='*int(100.*(i+1)/len(mat_indx)),
        #                                                100.*(i+1)/len(mat_indx)))
        #        sys.stdout.flush()

        #    if i%2000 == 0:
        #        sys.stdout.write('\r')
        #        sys.stdout.write("[%-100s] %d%% -|-" % ('='*int(100.*(i+1)/len(mat_indx)),
        #                                                100.*(i+1)/len(mat_indx)))
        #        sys.stdout.flush()
        #    import ipdb; ipdb.set_trace()

        #    if not str(mat_indx[i]) in mat_azi_all.keys():
        #        mat_azi_all[str(mat_indx[i])] = [mat_azi[i]]
        #    else:
        #        mat_azi_all[str(mat_indx[i])].append(mat_azi[i])

    ##mat_azi_qual = [[0.0, 0.0, 0.0, 0.0]]*len(mesh_points)
    #mat_azi_qual = np.zeros([len(mesh_points), 4])
    #for ep in mat_azi_all.keys():
    #    qc_node = np.histogram(mat_azi_all[ep], bins=range(0, 405, 45))
    #    #qc_node = np.histogram(mat_azi_all[ep], bins=[0, 90, 180, 270, 360])
    #    qc_node_combined = np.array([qc_node[0][0] + qc_node[0][4],
    #                                 qc_node[0][1] + qc_node[0][5],
    #                                 qc_node[0][2] + qc_node[0][6],
    #                                 qc_node[0][3] + qc_node[0][7],
    #                                 ])
    #    #if len(np.nonzero(qc_node_combined)[0])
    #    mat_azi_qual[int(ep)] += qc_node_combined

    mat_azi_qual_rank = np.zeros(len(mesh_points))
    for i in range(len(mat_azi_qual)):
        non_zero = np.nonzero(mat_azi_qual[i])[0]
        if len(non_zero) == 4:
            if min(non_zero) >= 10:
                mat_azi_qual_rank[i] = 100
            else:
                mat_azi_qual_rank[i] = 80
        elif len(non_zero) == 3:
            if min(non_zero) >= 10:
                mat_azi_qual_rank[i] = 90
            else:
                mat_azi_qual_rank[i] = 70
        elif len(non_zero) == 2:
            if (non_zero[1] - non_zero[0]) == 2:
                if min(np.nonzero(mat_azi_qual[i])[0]) >= 10:
                    mat_azi_qual_rank[i] = 80
                else:
                    mat_azi_qual_rank[i] = 60
            else:
                if min(np.nonzero(mat_azi_qual[i])[0]) >= 10:
                    mat_azi_qual_rank[i] = 40
                else:
                    mat_azi_qual_rank[i] = 20
        elif len(np.nonzero(mat_azi_qual[i])[0]) == 1:
            if min(np.nonzero(mat_azi_qual[i])[0]) >= 10:
                mat_azi_qual_rank[i] = 20
            else:
                mat_azi_qual_rank[i] = 10
        else:
            mat_azi_qual_rank[i] = 0

    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(mat_val_all)
    plt.subplot(2,1,2)
    plt.plot(mat_azi_qual_rank)
    plt.show()
    
    plt.figure()
    plt.subplot(1,1,1)
    plt.plot(mat_val_all[rad>=3482]*mat_azi_qual_rank[rad>=3482])
    plt.show()

    vtk = pvtk.VtkData(pvtk.UnstructuredGrid(mesh_points,
                                             tetra=mesh_facets),
                       pvtk.PointData(pvtk.Scalars(mat_val_all,
                                                   name='kernel_value')),
                       'Inversion Grid')
    vtk.tofile(os.path.join(direname, 'absolute_value_all.vtk'))

    vtk_azi = pvtk.VtkData(pvtk.UnstructuredGrid(mesh_points,
                                                 tetra=mesh_facets),
                           pvtk.PointData(pvtk.Scalars(mat_azi_qual_rank,
                                                       name='Quality_value')),
                           'Quality check')
    vtk_azi.tofile(os.path.join(direname, 'azi_quality_all.vtk'))
    print('\n\n==================')
    print("WARNING: INDEXING!")
    print('==================')
    raw_input('Press Enter to quit...')

# ###################### make_assemble_dir #############################


def make_assemble_dir(master_add, target_add, vertex_file, facet_file,
                      dlnvp, dlnvs, dlnQs,
                      hypocen_corr,
                      sta_corr_tp, sta_corr_ap,
                      orig_corr, ev_corr_ap,
                      sta_corr_ts, sta_corr_as,
                      corr_switches, out_ident):
    """
    Create assemblematrix directory to run assemblematrix.f
    :param master_add: copy matrices files from this directory
    :param target_add: where everything will be stored
    :param vertex_file:
    :param facet_file:
    :param dlnvp:
    :param dlnvs:
    :param dlnQs:
    :param hypocen_corr:
    :param sta_corr_tp:
    :param sta_corr_ap:
    :param orig_corr:
    :param ev_corr_ap:
    :param sta_corr_ts:
    :param sta_corr_as:
    :param corr_switches:
    :param out_ident:
    :return:
    """
    print ("copying matrixT.* files...",)
    mat_files_glob = glob.glob(os.path.join(master_add, 'matrixT.*'))
    #mat_files_glob = glob.glob(os.path.join(master_add, '*', 'matrixT.*'))
    mat_files_glob.sort()
    #for fi in mat_files_glob:
    #    shutil.copy(fi, os.path.join(target_add))
    #print "DONE"

    print ("copying raymatrix.info.T.* files...",)
    rayT_files_glob = glob.glob(os.path.join(master_add,
                                             'raymatrix.info.T.*'))
    #rayT_files_glob = glob.glob(os.path.join(master_add, '*',
    #                                         'raymatrix.info.T.*'))
    rayT_files_glob.sort()
    #for fi in rayT_files_glob:
    #    shutil.copy(fi, os.path.join(target_add))
    #print "DONE"

    print ("copying vertex and facet files...",)
    shutil.copy(os.path.join(os.curdir, 'src', 'files',
                             vertex_file),
                os.path.join(target_add))
    shutil.copy(os.path.join(os.curdir, 'src', 'files',
                             facet_file),
                os.path.join(target_add))
    print ("DONE")

    print ("creating the input files...",)
    in_assemble_fio = open(os.path.join(target_add, 'in.am'), 'w')
    in_assemble_fio.writelines('%s\n' % vertex_file)
    in_assemble_fio.writelines('%s\n' % facet_file)
    in_assemble_fio.writelines('%s %s\n' % (dlnvp[0], dlnvp[1]))
    in_assemble_fio.writelines('%s %s\n' % (dlnvs[0], dlnvs[1]))
    in_assemble_fio.writelines('%s %s\n' % (dlnQs[0], dlnQs[1]))
    in_assemble_fio.writelines('%s %s\n' % (hypocen_corr[0], hypocen_corr[1]))
    in_assemble_fio.writelines('%s %s\n' % (sta_corr_tp[0], sta_corr_tp[1]))
    in_assemble_fio.writelines('%s %s\n' % (sta_corr_ap[0], sta_corr_ap[1]))
    in_assemble_fio.writelines('%s %s\n' % (orig_corr[0], orig_corr[1]))
    in_assemble_fio.writelines('%s %s\n' % (ev_corr_ap[0], ev_corr_ap[1]))
    in_assemble_fio.writelines('%s %s\n' % (sta_corr_ts[0], sta_corr_ts[1]))
    in_assemble_fio.writelines('%s %s\n' % (sta_corr_as[0], sta_corr_as[1]))
    in_assemble_fio.writelines('%s %s %s %s\n' % (corr_switches[0],
                                                  corr_switches[1],
                                                  corr_switches[2],
                                                  corr_switches[3]))
    in_assemble_fio.writelines('%s\n' % out_ident)
    for fi in mat_files_glob:
        in_assemble_fio.writelines('%s\n' % os.path.basename(fi))
    in_assemble_fio.writelines('stop')
    in_assemble_fio.close()
    print ("DONE")

# ###################### compile_assemblematrix #############################


def compile_assemblematrix(target_add):
    """
    compile assemblematrix.f
    :return:
    """

    cur_dir = os.path.abspath(os.curdir)
    os.chdir(os.path.join(os.curdir, 'src',
                          'assemblematrix'))
    os_sys = os.system('./make')
    if not os_sys == 0:
        os.exit("assemblematrix can not be compiled properly")
    os.chdir(cur_dir)

    shutil.copy(os.path.join(os.curdir, 'src',
                             'assemblematrix', 'assemblematrix'),
                os.path.join(target_add))

# ###################### run_assemblematrix #############################


def run_assemblematrix(target_add):
    """
    run both raydata and raymatrix in the directory
    """
    cur_dir = os.path.abspath(os.curdir)

    print( '\n======>> run assemblematrix at %s' % target_add)
    os.chdir(os.path.join(target_add))
    os_sys = os.system('./assemblematrix < in.am')
    if not os_sys == 0:
        print( 'assemblematrix was not executed correctly!')

    os.chdir(cur_dir)

# ###################### station_filter #############################


def station_filter(ls_stas, min_xcorr=-100, max_xcorr=100, min_epi=0.,
                   max_epi=360., check_clip=True):
    """
    Filters the stations based on the required inputs
    NOT BEING USED AT THIS MOMENT!
    """
    empty_array = np.array([])
    pass_stas_corr_1 = ls_stas[min_xcorr <= ls_stas[:, 6].astype(np.float)]

    if not pass_stas_corr_1.size == 0:
         pass_stas_corr_2 = pass_stas_corr_1[
             max_xcorr > pass_stas_corr_1[:, 6].astype(np.float)]
    else:
        return empty_array

    passed_stas_epi_1 = pass_stas_corr_2[
        min_epi <= pass_stas_corr_2[:, 4].astype(np.float)]
    if not passed_stas_epi_1.size == 0:
        passed_stas_epi_2 = passed_stas_epi_1[
            max_epi > passed_stas_epi_1[:, 4].astype(np.float)]
    else:
        return empty_array
    if passed_stas_epi_2.size == 0:
        return empty_array

    if check_clip:
        passed_stas = passed_stas_epi_2[
            passed_stas_epi_2[:, 19].astype(np.float) < 0.1]

    return passed_stas

# ###################### plot_assemble_dtheor #############################


def plot_assemble_dtheor(assembled_dir, out_ident):
    """
    Plot nonzero, rhs and dtheor created by assemblematrix.f
    :param assembled_dir:
    :param out_ident:
    :return:
    """

    plt.ion()

    fio = open(os.path.join(assembled_dir, 'aux.%s' % out_ident), 'r')
    fi = fio.readlines()
    nonzero = []
    rhs = []
    dtheor = []
    for li in range(0, len(fi)):
        if 'End-of-matrix' in fi[li]:
            break
        nonzero.append(int(fi[li].split()[1]))
        rhs.append(float(fi[li].split()[2]))
        dtheor.append(float(fi[li].split()[14]))
    plt.figure()
    plt.plot(range(1, len(nonzero)+1), nonzero, lw=3)
    plt.xlabel('#measurement', size=24, weight='bold')
    plt.ylabel('#nonzero', size=24, weight='bold')
    plt.xticks(size=18, weight='bold')
    plt.yticks(size=18, weight='bold')

    plt.figure()
    plt.plot(range(1, len(rhs)+1), rhs, lw=3)
    plt.xlabel('#measurement', size=24, weight='bold')
    plt.ylabel('RHS', size=24, weight='bold')
    plt.xticks(size=18, weight='bold')
    plt.yticks(size=18, weight='bold')

    plt.figure()
    plt.plot(range(1, len(dtheor)+1), dtheor, lw=3)
    plt.xlabel('#measurement', size=24, weight='bold')
    plt.ylabel('dtheor', size=24, weight='bold')
    plt.xticks(size=18, weight='bold')
    plt.yticks(size=18, weight='bold')

    plt.show()

# ###################### plot_assemble_stations #############################


def plot_assemble_stations(assembled_dir, out_ident):
    """
    plot two maps: stations and number of usage + dT(av)
    :param assembled_dir:
    :param out_ident:
    :return:
    """

    plt.ion()

    fio = open(os.path.join(assembled_dir,
                            'assemblematrix.stations.%s' % out_ident), 'r')
    fi = fio.readlines()

    stla = []
    stlo = []
    stnn = []
    stdt = []
    for li in range(2, len(fi)):
        stla.append(float(fi[li].split()[6]))
        stlo.append(float(fi[li].split()[7]))
        stnn.append(float(fi[li].split()[2]))
        stdt.append(float(fi[li].split()[4]))

    plt.figure()
    m = Basemap(projection='robin', lon_0=0.0, lat_0=0.0, resolution='c')
    m.fillcontinents()
    m.drawparallels(np.arange(-90., 120., 30.))
    m.drawmeridians(np.arange(0., 420., 60.))
    m.drawmapboundary()

    x, y = m(stlo, stla)
    m.scatter(x, y, s=100, c=stnn, marker="o", edgecolor="none", zorder=10,
            vmax=1000)
    plt.colorbar()
    plt.title('Number of repetition', size=24, weight='bold')

    plt.figure()
    m = Basemap(projection='robin', lon_0=0.0, lat_0=0.0, resolution='c')
    m.fillcontinents()
    m.drawparallels(np.arange(-90., 120., 30.))
    m.drawmeridians(np.arange(0., 420., 60.))
    m.drawmapboundary()

    m.scatter(x, y, s=100, c=stdt, marker="o", edgecolor="none", zorder=10,
            vmin=-7, vmax=7)
    plt.colorbar()
    plt.title('Average dT', size=24, weight='bold')

    plt.show()

# ###################### plot_assemble_columndensity ##########################


def plot_assemble_columndensity(assembled_dir, out_ident,
                                vertex_file, facet_file):
    """
    plot columndensity created by assemblematrix.f
    :param assembled_dir:
    :param out_ident:
    :param vertex_file:
    :param facet_file:
    :return:
    """

    print ('\n======>> Creating VTK file')
    columnden_files = glob.glob(os.path.join(assembled_dir,
                                             'columndensity.*'))
    if not len(columnden_files) == 1:
        sys.exit("Length of columndensity.* files is more than 1!")

    print ('------> Load columndensity file')
    mesh_points = np.loadtxt(os.path.join(columnden_files[0]),
                             skiprows=2, comments='#')
    print ('------> Load facet file')
    mesh_facets = np.loadtxt(os.path.join(assembled_dir, facet_file),
                             dtype=np.int, skiprows=1, comments='#')

    print ('------> Writing to the disk')
    vtk = pvtk.VtkData(pvtk.UnstructuredGrid(mesh_points[:, 0:3],
                                             tetra=mesh_facets),
                       pvtk.PointData(pvtk.Scalars(mesh_points[:, 3],
                                                   name='column_density')),
                       'Inversion Grid')
    print ('\n\n==================')
    print ("WARNING: INDEXING!")
    print ('==================')
    vtk.tofile(os.path.join(assembled_dir, 'columndensity.vtk'))

# ###################### make_mpisolve #############################


def make_mpisolve(assembled_dir, ident1, ident2,
                  chi_rel, ksmooth, epsnorm, epssmooth, ratio, outliers_sigma,
                  max_num_search1, max_num_search2, vertex_file, facet_file):
    """
    create an input file for mpisolvetomo.f
    :param assembled_dir:
    :param ident1:
    :param ident2:
    :param chi_rel:
    :param ksmooth:
    :param epsnorm:
    :param epssmooth:
    :param ratio:
    :param outliers_sigma:
    :param max_num_search1:
    :param max_num_search2:
    :param vertex_file:
    :param facet_file:
    :return:
    """

    print ("creating the input files...",)
    in_mpisolve_fio = open(os.path.join(assembled_dir, 'in.solvetomo'), 'w')
    in_mpisolve_fio.writelines('%s\n' % os.path.abspath(assembled_dir))
    in_mpisolve_fio.writelines('%s\n' % ident1)
    in_mpisolve_fio.writelines('%s\n' % ident2)
    in_mpisolve_fio.writelines('%s %s %s %s %s\n' % (chi_rel,
                                                     ksmooth,
                                                     epsnorm,
                                                     epssmooth,
                                                     ratio))
    in_mpisolve_fio.writelines('%s\n' % outliers_sigma)
    in_mpisolve_fio.writelines('%s %s\n' % (max_num_search1, max_num_search2))
    in_mpisolve_fio.writelines('%s\n' % vertex_file)
    in_mpisolve_fio.writelines('%s' % facet_file)
    in_mpisolve_fio.close()
    print ("DONE")

# ###################### compile_mpisolve #############################


def compile_mpisolve(target_add):
    """
    compile mpisolvetomo.f
    :return:
    """

    cur_dir = os.path.abspath(os.curdir)
    os.chdir(os.path.join(os.curdir, 'src',
                          'mpisolvetomo'))
    os_sys = os.system('./make')
    if not os_sys == 0:
        os.exit("mpisolvetomo.f can not be compiled properly")
    os.chdir(cur_dir)

    shutil.copy(os.path.join(os.curdir, 'src',
                             'mpisolvetomo', 'mpisolvetomo'),
                os.path.join(target_add))

# ###################### create_slurmfile #############################


def create_slurmfile(target_add, partition, sol_name, ntasks, nodes):
    """
    Create slurm file to be used on tethys for running mpisolvetomo
    :param target_add:
    :param partition:
    :param dir_name:
    :param sol_name:
    :param ntasks:
    :param nodes:
    :return:
    """

    print ("creating the slurm file...",)
    """
    #!/bin/bash
    #SBATCH --partition=BIG
    #SBATCH --output=/SCRATCH/hosseini/INVERSION/slurm_output/inversion_%j.out
    #SBATCH --workdir=/SCRATCH/hosseini/INVERSION/00001_P_all_freq
    #SBATCH --job-name="INV-00001"
    #SBATCH --ntasks=264
    #SBATCH --nodes=11
    #SBATCH --mail-type=all
    #SBATCH --time=719:59:59
    # run compute job

    source /home/hosseini/.bashrc
    mpirun.openmpi -np 264 mpisolvetomo < in.st
    """

# ###################### plot_stas_proj #############################


def plot_stas_proj(filt_array_corr, inp, plot_results=True, indx=35):
    """
    Plotting the measurements on the defined projection
    :param filt_array_corr:
    :param inp:
    :return:
    """
    if plot_results:
        plt.figure()
        indx_remove = filt_array_corr[:, 35] == 'k'
        filt_black = filt_array_corr[indx_remove]
        filt_array_corr = filt_array_corr[filt_array_corr[:, 35] != 'k']
        mymap = Basemap(projection=inp.plot_stas_proj_type,
                        lat_0=inp.plot_lat_0,
                        lon_0=inp.plot_lon_0)
        mymap.fillcontinents()
        x, y = mymap(filt_array_corr[:, 3].astype(np.float),
                     filt_array_corr[:, 2].astype(np.float))
        x_k, y_k = mymap(filt_black[:, 3].astype(np.float),
                         filt_black[:, 2].astype(np.float))

    sel_mean = False
    if inp.plot_selected_stas:
        sel_mean = []
        for sel_sta in inp.plot_selected_stas:
            for filt_line in filt_array_corr:
                if sel_sta in filt_line[5]:
                    sel_mean.append(filt_line[indx])
        # print sel_mean
        if len(sel_mean) >= inp.num_sta_thresh_median:
            all_mean = np.mean(sel_mean)
        else:
            all_mean = np.mean(filt_array_corr[:, indx])
    else:
        all_mean = np.mean(filt_array_corr[:, indx])

    if inp.plot_fixed_median:
        all_mean = inp.plot_fixed_median
        print ("Fixed mean is selected to be: %s" % all_mean)
    # 0274.2009.273.a for both P and Pdiff:
    # all_mean = 3.4278005835999998
    # 0224.2009.156.a for both P and Pdiff
    # all_mean = 2.3435461092000001
    if sel_mean:
        print ('mean based on %s GSN stations: %s vs %s (all data)' % \
              (len(sel_mean), all_mean, np.mean(filt_array_corr[:, indx])))

    c_plot = filt_array_corr[:, indx] - all_mean
    if plot_results:
        mymap.scatter(x, y, c=c_plot.astype(np.float),
                      edgecolor='none',
                      zorder=20, s=100, vmin=-2, vmax=2)
        # python dispersion_plot.py npy_output/plot_Pdiff_band01_Pdiff_corr.npy
        # npy_output/plot_Pdiff_band01_Pdiff_corr.npy
        # npy_output/input_plot_Pdiff_Pdiff.pkl
        # plt.title('dT(Pdiff; 30s)', size=46)

        # python dispersion_plot.py npy_output/plot_P_band01_P_corr.npy
        # npy_output/plot_P_band01_P_corr.npy npy_output/input_plot_P_P.pkl
        # plt.title('dT(P; 30s)', size=46)

        # python dispersion_plot.py npy_output/plot_Pdiff_band01_Pdiff_corr.npy
        # npy_output/plot_P_band01_P_corr.npy
        # npy_output/input_plot_Pdiff_Pdiff.pkl
        # plt.title('dT(Pdiff; 30s) - dT(P; 30s)', size=46)

        # python dispersion_plot.py npy_output/plot_Pdiff_band01_Pdiff_corr.npy
        # npy_output/plot_Pdiff_band03_Pdiff_corr.npy
        # npy_output/input_plot_Pdiff_Pdiff.pkl
        # plt.title('dT(Pdiff; 30s) - dT(Pdiff; 15s)', size=46)

        # python dispersion_plot.py npy_output/plot_P_band01_P_corr.npy
        # npy_output/plot_P_band03_P_corr.npy npy_output/input_plot_P_P.pkl
        # plt.title('dT(P; 30s) - dT(P; 15s)', size=46)
        cbar = plt.colorbar()
        # cbar.set_ticks(np.arange(-0.5, 0.75, 0.25))
        cbar.ax.tick_params(labelsize=42)
        mymap.scatter(x_k, y_k, c='k',
                      edgecolor='none',
                      zorder=20, s=100)
        if True:
            # Only for the measurement paper 2015
            mymap.drawgreatcircle(143.4450, 41.8240,
                                  inp.plot_lon_0, inp.plot_lat_0,
                                  ls='-', lw=6, color='r')
            mymap.drawgreatcircle(99.8670, -0.7200,
                                  inp.plot_lon_0, inp.plot_lat_0,
                                  ls='-', lw=6, color='b')
            xmin, ymin = mymap(-135, 20)
            xmax, ymax = mymap(-48, 56)
            plt.axis([xmin, xmax, ymin, ymax])
        plt.show()

        if True:
            # Only for the measurement paper 2015
            from obspy.imaging.beachball import Beach
            plt.figure()
            h = 300000000
            mymap = Basemap(projection='nsper',
                            lat_0=-0.72+40, lon_0=99.867+60,
                            satellite_height=h*1000., resolution='l')
            mymap.fillcontinents()
            mymap.drawparallels(np.arange(-90., 120., 30.))
            mymap.drawmeridians(np.arange(0., 420., 60.))

            x, y = mymap(inp.plot_lon_0, inp.plot_lat_0)
            mymap.scatter(x, y, c='r', s=600, linewidth=5,
                          marker='x', zorder=200)

            x, y = mymap(99.86700, -0.7200000)
            focmecs = [0.163E+21, -0.148E+20, -0.148E+21,
                       0.349E+20, -0.397E+20, -0.157E+21]
            ax = plt.gca()
            b = Beach(focmecs, xy=(x, y), width=7e5,
                      linewidth=1, alpha=0.85, size=400)
            b.set_zorder(10)
            ax.add_collection(b)

            x, y = mymap(143.4450, 41.8240)
            focmecs = [1.810e+18, -1.270e+18, -5.310e+17,
                       8.590e+17, 2.550e+18, -2.860e+17]
            ax = plt.gca()
            b = Beach(focmecs, xy=(x, y), width=7e5,
                      linewidth=1, alpha=0.85, size=400)
            b.set_zorder(10)
            ax.add_collection(b)

            mymap.drawgreatcircle(143.4450, 41.8240,
                                  inp.plot_lon_0, inp.plot_lat_0,
                                  ls='-', lw=6, color='r')
            x_txt, y_txt = mymap(143.4450+18, 41.8240)
            plt.text(x_txt, y_txt, 'EVENT 2', fontsize=24, fontweight='bold',
                     ha='center', va='center', color='k')
            mymap.drawgreatcircle(99.8670, -0.7200,
                                  inp.plot_lon_0, inp.plot_lat_0,
                                  ls='-', lw=6, color='b')
            x_txt, y_txt = mymap(99.8670+10, -0.7200)
            plt.text(x_txt, y_txt, 'EVENT 1', fontsize=24, fontweight='bold',
                     ha='center', va='center', color='k')

            plt.show()

    return all_mean

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


# ------------------- bc color --------------------

class bc:
    yellow = '\033[93m'
    orange = '\033[0;33m'

    bred = '\033[1;31m'
    red = '\033[91m'
    dred = '\033[2;31m'

    bblue = '\033[1;34m'
    blue = '\033[94m'
    dblue = '\033[2;34m'

    bgreen = '\033[1;32m'
    green = '\033[92m'
    dgreen = '\033[2;32m'

    bmagenta = '\033[1;35m'
    magenta = '\033[95m'
    dmagenta = '\033[2;35m'

    bcyan = '\033[1;96m'
    cyan = '\033[96m'
    dcyan = '\033[2;96m'

    bwhite = '\033[1;97m'
    white = '\033[97m'
    black = '\033[0;30m'

    bgrey = '\033[1;90m'
    grey = '\033[90m'
    dgrey = '\033[2;90m'

    end = '\033[0m'
    bold = '\033[1m'
    curs = '\033[3m'
    under = '\033[4m'

# ------------------- get_time --------------------


def get_time():
    time = datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S')
    return time

# ------------------- cprint --------------------


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

    print(bc.green + get_time() + bc.end, \
        bc.curs + bc.magenta + ho_nam + bc.end, \
        bc.cyan + code + bc.end, \
        bc.bwhite + type_info + bc.end, \
        bc_color + text + bc.end)
