#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  parallel_helper.py
#   Purpose:   helping functions for parallelization
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
import numpy as np
import os
import sys
import time

from src.EventClass import read_events, filter_events
from src.InpFFproc import InpFFproc
from src.output_writer import create_dir
import src.utility_codes as util_plot
from src.utility_codes import bc, cprint

# ------------------ rank_0_initialization ---------------------------


def rank_0_initialization(inp_file, size):
    """
    root initializations: inp_ffproc, read_events, filter_events, ...
    :param inp_file:
    :param size:
    :return:
    """
    # ------------------- INPUT -----------------------------
    inp_ffproc_file = os.path.abspath(inp_file)
    inp_ffproc = InpFFproc(inp_ffproc_file=inp_ffproc_file)
    # --------------- Identify all the required events
    all_events = read_events(inp_ffproc=inp_ffproc)
    filter_events(inp_ffproc=inp_ffproc, all_events=all_events)
    create_dir(inp_ffproc.output_dir)
    # Plot all the selected events
    if inp_ffproc.check_evs:
        util_plot.plot_events(all_events=all_events,
                              add_save=inp_ffproc.output_dir)
    # --------------- make sure instaseis is run serially
    all_before_func = np.array(size*[False])
    live_funcs = np.array(size*[True])
    # --------------- Check the filter duration
    if inp_ffproc.filt_mode == 'log-gabor':
        from src.filtering import filter_duration_calc
        filter_duration_calc(inp_ffproc=inp_ffproc)
    return inp_ffproc, all_events, all_before_func, live_funcs

# ------------------ chunk_events ---------------------------


def chunk_events(list2chunk, num_chunks):
    """
    list2chunk into num_chunks sub-lists
    :param list2chunk:
    :param num_chunks:
    :return:
    """
    newseq = []
    splitsize = 1.0 / num_chunks * len(list2chunk)
    for i in range(num_chunks):
        newseq.append(
            list2chunk[int(round(i * splitsize)):int(round((i + 1) * splitsize))])
    return newseq

# ------------------ monitor_par2ser_func ---------------------------


def monitor_par2ser_func(size, comm, all_before_func, free_func, live_funcs):
    """
    monitor the function that should be run in serial
    :param size:
    :param comm:
    :param all_before_func:
    :param free_func:
    :param live_funcs:
    :return:
    """
    # root is running, so np.sum > 1
    while np.sum(live_funcs) > 1:
        time.sleep(0.1)
        for ik in range(1, size):
            if comm.iprobe(source=ik, tag=ik):
                all_before_func[ik] = True
            if comm.iprobe(source=ik, tag=ik+500):
                live_funcs[ik] = False
                dummy_end_msg = comm.recv(source=ik, tag=ik+500)
        num_be4func = all_before_func[all_before_func > 0]
        if len(num_be4func) > 0 and free_func:
            cprint('parallel_helper.py', '[PARALLEL]', bc.dred,
                   '#Processes waiting for instaseis: %s' % len(num_be4func)) 
            rk = np.where(all_before_func > 0)[0][0]
            tmp_before = comm.recv(source=rk, tag=rk)
            comm.send(True, dest=rk, tag=rk+100)
            all_before_func[rk] = False
            free_func = comm.recv(source=rk, tag=rk+1000)
    sys.exit()
