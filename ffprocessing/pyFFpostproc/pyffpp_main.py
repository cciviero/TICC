#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  pyffpp_main.py
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

import glob
import os
import sys
import time

from pyffpp_src import input_reader, read_kernel_filter, \
    measurement_info, plot_hists, bc, event_wise_stat_plots
from pyffpp_src import create_vtk_input, make_vtk, plot_amplitude

# -----------------------------------------------------------------------
# ------------------------------ MAIN PART ------------------------------
# -----------------------------------------------------------------------

tic = time.time()

try:
    inp_file = sys.argv[1]
except Exception, e:
    print bc.dred + bc.bold + '[ERROR] No input_pyffpp.ini given! I am out!' + bc.end
    sys.exit(2)

# read input file
inp_ffpp = input_reader(inp_file=inp_file)

path_events = glob.glob(os.path.join(inp_ffpp.input_path, '*.a'))
all_events = []
for i in path_events:
    join = i.split('/')[-1]
    all_events.append(join)
# import ipdb; ipdb.set_trace()
inp_band = read_kernel_filter(inp_ffpp.input_path, all_events[0])
bands_sec = inp_band['dominant_period']

# -----------------------------------------------------------------------

if inp_ffpp.all_statistics:
    print bc.yellow + bc.bold + 'Overall statistic plot of %s ' \
                                 'will be generated!' % inp_ffpp.input_path + bc.end
    final_cc, final_ts, final_amp, final_mag, final_azi, output_dir, count_stations, \
    count_criteria_true, count_criteria_false = measurement_info(inp_ffpp.input_path,
                                                                 all_events, bands_sec, inp_ffpp)

    plot_hists(final_cc, inp_ffpp, bands_sec, count_stations,
               count_criteria_true, count_criteria_false,
               'cross correlation [%]', 'cc', output_dir)

    plot_hists(final_ts, inp_ffpp, bands_sec, count_stations,
               count_criteria_true, count_criteria_false,
               'time shift [s]', 'time_shift', output_dir)

if inp_ffpp.amplitude_plot:
    plot_amplitude(final_cc, final_ts, final_amp, final_mag, final_azi, output_dir)


if inp_ffpp.event_wise_stat:
    print bc.yellow + bc.bold + 'Event wise statistics of %s ' \
                                'will be generated!' % inp_ffpp.input_path + bc.end
    if inp_ffpp.selected_events:
        sel_ev = inp_ffpp.selected_events_list
        event_wise_stat_plots(inp_ffpp, inp_ffpp.input_path, bands_sec, sel_ev)
    else:
        event_wise_stat_plots(inp_ffpp, inp_ffpp.input_path, bands_sec, all_events)

if inp_ffpp.proc_vtk:

    if inp_ffpp.vtk_existing_file:
        print bc.yellow + bc.bold + 'You selected your own txt files for generating vtk files! \n' \
                                    'Lets get started!' + bc.end
        list_txt2vtk_file = []
        for file in inp_ffpp.vtk_input_files:
            list_txt2vtk_file.append(os.path.join(inp_ffpp.vtk_input_path, file))

        make_vtk(list_txt2vtk_file, inp_ffpp.vtk_input_path)

    else:
        print bc.yellow + bc.bold + 'VTK file option //%s//. Creating vtk files for %s!' \
                                 % (inp_ffpp.vtk_file_opt, inp_ffpp.input_path) + bc.end
        list_txt2vtk_file, vtk_path = create_vtk_input(all_events, bands_sec, inp_ffpp)
        make_vtk(list_txt2vtk_file, vtk_path)


# -----------------------------------------------------------------------
# ------------------------------- THE END  ------------------------------
# -----------------------------------------------------------------------

toc = time.time()
print bc.magenta + bc.bold + 'Finished. \n +++ It took %s min. +++ \n' % (round((toc-tic)/60, 2))