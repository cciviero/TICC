#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  ffprocessing.py
#   Purpose:   Multi-frequency waveform measurement
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")

from mpi4py import MPI
import numpy as np
import os
import sys
import time

import src.utility_codes as util_plot
from src.output_writer import output_writer, kernel_writer
from src.parallel_helper import rank_0_initialization, chunk_events
from src.parallel_helper import monitor_par2ser_func
from src.SNR import cSNR
from src.StationClass import read_syn, read_data_header
from src.StationClass import read_data, filter_stations
from src.STFclass import conv_stf
from src.utility_codes import bc, cprint, exitus_file, check_exitus_file
from src.utility_codes import no_inp_file_exit, print_ffproc_version

# ------------- reading input path from command line
try:
    inp_file = sys.argv[1]
    # ------------- only print the version and exit
    if inp_file.lower() == 'version':
        print_ffproc_version()
except Exception, e:
    no_inp_file_exit()
# ------------- Parallel (initialization)
par_comm = MPI.COMM_WORLD
par_rank = par_comm.rank
par_size = par_comm.size
inp_ffproc = None
all_events = None
list_selected_events = None
# func: instaseis
free_func = True
# ------------- only the root does the first part
if par_rank == 0:
    inp_ffproc, all_events, all_before_func, live_funcs = \
        rank_0_initialization(inp_file, par_size)

    # ------------- Selected Events
    list_selected_events = []
    if inp_ffproc.selected_events:
        dir_selected_events = os.path.join(os.path.curdir, 'selected_events.txt')
        list_selected_events = np.loadtxt(dir_selected_events,
                                          comments='#',
                                          dtype='object')
# import ipdb; ipdb.set_trace()
# ------------- broadcast inp_ffproc and all_events
par_comm.Barrier()
list_selected_events = par_comm.bcast(list_selected_events, root=0)
inp_ffproc = par_comm.bcast(inp_ffproc, root=0)
all_events = par_comm.bcast(all_events, root=0)
if all_events.number_events < 1:
    sys.exit('[EXIT] #events: %s' % all_events.number_events)
par_comm.Barrier()

if par_rank == 0:
    ev_remove_indx = np.zeros(all_events.number_events)
    for ev in range(all_events.number_events):
        if inp_ffproc.selected_events:
            if all_events.name[ev] in list_selected_events:
                ev_remove_indx[ev] = 1
            else:
                ev_remove_indx[ev] = 0
                continue
    
        if not inp_ffproc.force_process:
            complete_event = check_exitus_file(all_events.name[ev], inp_ffproc)
            if not complete_event:
                ev_remove_indx[ev] = 1
            else:
                ev_remove_indx[ev] = 0
        else:
            ev_remove_indx[ev] = 1
    all_events.rmev(ev_remove_indx.astype(np.bool))

par_comm.Barrier()
all_events = par_comm.bcast(all_events, root=0)

# ------------- if more than one processor:
# ffprocessing uses master-slave mpi approach in which:
# rank-0: manager, runs monitor_par2ser_func
# rank-x (x>0): run the rest of the code
if par_rank == 0 and par_size > 1:
    monitor_par2ser_func(par_size, par_comm, all_before_func,
                         free_func, live_funcs)
# ------------- Chunk events for parallel requests
if par_size > 1:
    num_chunk_ev = par_size - 1
else:
    num_chunk_ev = 1

chunk_ev = chunk_events(all_events.name, num_chunk_ev)[par_rank-1]
cprint('ffprocessing.py', '[PARALLEL]', bc.lmagenta,
       'Rank %03i has %s events' % (par_rank, len(chunk_ev)))

# --------------- Loop over events
tic = time.time()
for ev in range(all_events.number_events):
    ev_tic = time.time()
    if not all_events.name[ev] in chunk_ev:
        continue
    if inp_ffproc.selected_events:
        if not all_events.name[ev] in list_selected_events:
            continue
    if not inp_ffproc.force_process:
        complete_event = check_exitus_file(all_events.name[ev], inp_ffproc)
        if complete_event:
            continue
    cprint('ffprocessing.py', '[EVENT]', bc.lblue,
           'Processing: %s' % all_events.name[ev])
    #import ipdb; ipdb.set_trace()
    # --------------- Read data header and choose the relevant ones
    all_stations = read_data_header(inp_ffproc=inp_ffproc,
                                    all_events=all_events,
                                    ev=ev)
    # import ipdb; ipdb.set_trace()
    filter_stations(inp_ffproc=inp_ffproc, all_stations=all_stations)

    # --------------- Fill in the arrival times
    all_stations.add_ttime(inp_ffproc.ph_phase,
                           all_events.inv_dp[ev],
                           inp_ffproc.ph_bg)
    # removing those stations that do not contain the requested phase
    indxs = all_stations.tt_ph != False
    all_stations.rmsta(indxs)
    if len(all_stations.name) == 0:
        cprint('ffprocessing.py', '[WARNING]', bc.orange,
               'No station available for %s -- '
               'Continue with the next event.' % all_events.name[ev])
        continue
    # --------------- Read real data
    # XXX if ?HT or ?HR is selected, first read Z component
    # and then take care of it later in NE2RT function
    real_waveforms = read_data(all_events=all_events,
                               all_stations=all_stations,
                               ev=ev,
                               inp_ffproc=inp_ffproc)
    # import ipdb; ipdb.set_trace()
    #XXXX commented out for now XXX cleaning up all_stations
    all_stations.rmsta(all_stations.real_rdable > 0)
    if len(all_stations.name) == 0:
        cprint('ffprocessing.py', '[WARNING]', bc.red,
               'No station available for %s -- '
               'Continue with the next event.' % all_events.name[ev])
        continue
    # --------------- NE 2 RT
    # XXX based on specified channels, i.e., if we have ?HR or ?HT,
    # it rotates and continues accordingly
    # XXX plot_RT should not be hard coded to True, add a new flag/option in the input file.
    # XXX does this rotate all horizontal components  ??!?!
    if 's' in inp_ffproc.ph_phase.lower():
        from src.rotation_helper import NE_2_RT
        NE_2_RT(all_stations, real_waveforms, inp_ffproc,
                              all_events, ev, plot_RT=True)


    # For additional pre filtering, uncomment:
    # real_waveforms.filter('bandpass', freqmin=1./40, freqmax=1./1)

    # Plotting the theoretical arrival times
    if inp_ffproc.check_ttime:
        util_plot.plot_ttime(all_stations=all_stations,
                             add_save=os.path.join(inp_ffproc.output_dir,
                                                   all_events.name[ev]))
    # Plotting the selected events and stations
    if inp_ffproc.check_ev_stas:
        util_plot.plot_event_stations(all_events=all_events,
                                      all_stations=all_stations,
                                      ev=ev,
                                      inp_ffproc=inp_ffproc)
    # plot the real data
    if inp_ffproc.check_cut_real:
        util_plot.plot_cut_wave(all_stations=all_stations,
                                waveforms=real_waveforms,
                                add_save=os.path.join(inp_ffproc.output_dir,
                                                      all_events.name[ev]),
                                fig_name=
                                'wave_3_reals_%s' % all_events.name[ev],
                                ev_datetime=all_events.datetime[ev])

    # --------------- Convolve STF with Synthetic waveforms
    cprint('ffprocessing.py', '[STF]', bc.lblue, 'Reading STF(s)')
    all_stfs = conv_stf(
        address_stf=os.path.join(inp_ffproc.evproc_path, all_events.name[ev]),
        all_stations=all_stations,
        inp_ffproc=inp_ffproc,
        all_events=all_events,
        ev=ev)

    # --------------- Read synthetic
    if par_size > 1:
        before_insta = True
        par_comm.send(before_insta, dest=0, tag=par_rank)
        allow2cont = par_comm.recv(source=0, tag=par_rank+100)
    syn_waveforms = read_syn(all_events=all_events,
                             all_stations=all_stations,
                             all_stfs=all_stfs,
                             inp_ffproc=inp_ffproc,
                             ev=ev)
    if par_size > 1:
        par_comm.send(True, dest=0, tag=par_rank+1000)

    # cleaning up all_stations based on the readability of synthetics
    # XXX commented out for now, maybe now the best choice... 
    # all_stations.rmsta(all_stations.syn_rdable > 0)

    # removing those synthetic waveforms for which
    # the real data could not be read
    indxs_rm = sorted(np.where(all_stations.real_rdable == 0)[0], reverse=True)
    for ind in indxs_rm:
        # syn_waveforms.remove(syn_waveforms[ind])
        del syn_waveforms[ind]
    if len(all_stations.name) == 0:
        cprint('ffprocessing.py', '[WARNING]', bc.orange,
               'No station available for %s -- '
               'Continue with the next event.' % all_events.name[ev])
        continue
    # print len(real_waveforms[0].data)
    # --------------- If requested integrate real and synthetic waveforms
    if inp_ffproc.int_real:
        cprint('ffprocessing.py', '[WARNING]', bc.orange,
               'Integrate real waveforms!')
        real_waveforms.integrate()
    if inp_ffproc.int_syn:
        cprint('ffprocessing.py', '[WARNING]', bc.orange,
               'Integrate synthetics!')
        syn_waveforms.integrate()

    # plot the synthetic waveforms
    if inp_ffproc.check_cut_syn:
        util_plot.plot_cut_wave(all_stations=all_stations,
                                waveforms=syn_waveforms,
                                add_save=os.path.join(inp_ffproc.output_dir,
                                                      all_events.name[ev]),
                                fig_name=
                                'wave_2_mfis_%s' % all_events.name[ev],
                                ev_datetime=all_events.datetime[ev])
    # print len(real_waveforms[0].data)
    # --------------- save waveforms
    if inp_ffproc.save_smgr:
        from src.output_writer import write_arr2tr
        write_arr2tr(real_waveforms, 'BB_real', all_events.name[ev],
                     inp_ffproc=inp_ffproc, all_stations=all_stations)
        write_arr2tr(syn_waveforms, 'BB_syn', all_events.name[ev],
                     inp_ffproc=inp_ffproc, all_stations=all_stations)
    if inp_ffproc.only_pre_process:
        cprint('ffprocessing.py', '[OUTPUT] [ARR2TR]', bc.cyan,
               'Continue to the next event (only_pre_process option)!')
        continue

    # print len(real_waveforms[0].data)

    # --------------- bandpass real and synthetic data
    # We do not have a main functionality for bandpass since different
    # methods require different inputs/outputs
    if inp_ffproc.filt_mode == 'log-gabor':
        from src.filtering import bandpass_gabor
        realBP, synBP, realREF, synREF = \
            bandpass_gabor(real_waveforms=real_waveforms,
                           syn_waveforms=syn_waveforms,
                           inp_ffproc=inp_ffproc,
                           all_stations=all_stations,
                           all_events=all_events,
                           ev=ev)

    # print len(real_waveforms[0].data), np.shape(realBP)

    # --------------- measurements
    # We do not have a main functionality for measurements since different
    # methods require different inputs/outputs
    if inp_ffproc.mmeant_mode == 'CC-2step':
        from src.measurement import xcorrelation_2step
        cut_realBP, cut_synBP = xcorrelation_2step(
            realBP=realBP,
            synBP=synBP,
            realREF=realREF,
            synREF=synREF,
            inp_ffproc=inp_ffproc,
            all_stations=all_stations,
            real_waveforms=real_waveforms,
            syn_waveforms=syn_waveforms,
            all_events=all_events,
            ev=ev)
    elif inp_ffproc.mmeant_mode == 'CC-2step-env':
        from src.measurement_envelope import xcorrelation_2step_env
        raw_input("CC-2step-env has not tested/cleaned-up!")
        cut_realBP, cut_synBP = xcorrelation_2step_env(
            inp_ffproc=inp_ffproc,
            realBP=realBP,
            synBP=synBP,
            realREF=realREF,
            synREF=synREF,
            all_stations=all_stations,
            real_waveforms=real_waveforms,
            syn_waveforms=syn_waveforms,
            all_events=all_events,
            ev=ev)

    # --------------- OUTPUT
    cprint('ffprocessing.py', '[SNR]', bc.lgreen,
           'Calculating SNR.')
    cSNR(realBP, real_waveforms, inp_ffproc, all_stations, all_events, ev)
    print '\n'
    
    # --------------- Measure particle motion
    # measure particle motion of horizontal component on
    # gabor filtered traces

    if inp_ffproc.particle_motion:
        from src.rotation_helper import calc_part_mot
        calc_part_mot(all_stations, real_waveforms, realBP, synBP,
                      inp_ffproc, all_events, ev)
    print '\n'

    # plot real and synthetic waveforms over each other (the measured travel
    # time is subtracted from the real data)
    if inp_ffproc.check_cut_real_syn:
        util_plot.plot_filt_compare_groups(all_stations=all_stations,
                                           realBP=realBP,
                                           synBP=synBP,
                                           cut_realBP=cut_realBP,
                                           cut_synBP=cut_synBP,
                                           inp_ffproc=inp_ffproc,
                                           add_save=os.path.join(
                                               inp_ffproc.output_dir,
                                               all_events.name[ev]),
                                           fig_name='compare')

    if inp_ffproc.check_cut_real_syn_gabor:
        # XXX
        util_plot.plot_spec_gabor_filter_compare(
            all_stations=all_stations,
            real_waveforms=real_waveforms,
            realBP=realBP,
            synBP=synBP,
            event_time=all_events.datetime[ev],
            cut_realBP=cut_realBP,
            cut_synBP=cut_synBP,
            inp_ffproc=inp_ffproc,
            add_save=os.path.join(inp_ffproc.output_dir, all_events.name[ev]),
            fig_name='and_gabor')

    # projecting the measured travel times on the glob
    if inp_ffproc.check_measured_glob:
        for band in range(inp_ffproc.check_min_bands,
                          inp_ffproc.check_max_bands):
            util_plot.plot_measured_glob(all_stations=all_stations,
                                         all_events=all_events,
                                         ev=ev,
                                         inp_ffproc=inp_ffproc,
                                         band=band)
    # --------------- output_writer
    cprint('ffprocessing.py', '[OUTPUT]', bc.dgreen,
           'Start generating the outputs')
    output_writer(inp_ffproc, all_stations, all_events, ev, all_stfs)

    # --------------- outputs for kernel
    if inp_ffproc.kernel_output:
        kernel_writer(inp_ffproc, all_events, ev, all_stations, all_stfs)

    # --------------- save the smgrs
    # XXX REVIEW
    # if inp_ffproc.save_smgr:
    #     from src.output_writer import write_cut_BB
    #     write_cut_BB(num_sta=all_stations.number_stations,
    #                  inp_ffproc=inp_ffproc,
    #                  cut_realBP=cut_realBP,
    #                  real_waveforms=real_waveforms,
    #                  all_events=all_events,
    #                  cut_synBP=cut_synBP,
    #                  syn_waveforms=syn_waveforms,
    #                  all_stfs=all_stfs,
    #                  ev=ev)

    ev_toc = time.time()
    exitus_file(all_events, ev, all_stfs, inp_ffproc, round((ev_toc-ev_tic)/60., 2))
    print 'Event %s \nTIME: %s min' \
          % (all_events.name[ev], round((ev_toc-ev_tic)/60., 2))

toc = time.time()
print "\n--------------------\n" \
      "TOTAL TIME (rank-%s): %smin" \
      "\n--------------------\n" \
      % (par_rank, round((toc-tic)/60., 2))
if par_size > 1:
    par_comm.send(False, dest=0, tag=par_rank+500)

# The following lines can/may be used later for collecting all the
# relevant variables for reproducibility

# =================================
# import pickle
# fio = open('realBP_%s.pkl' % inp_ffproc.ph_phase, 'w')
# pickle.dump(realBP, fio)
# fio = open('synBP_%s.pkl' % inp_ffproc.ph_phase, 'w')
# pickle.dump(synBP, fio)
# fio = open('inp_ffproc_%s.pkl' % inp_ffproc.ph_phase, 'w')
# pickle.dump(inp_ffproc, fio)
# fio = open('all_stations_%s.pkl' % inp_ffproc.ph_phase, 'w')
# pickle.dump(all_stations, fio)

# =================================
# plot the Green's functions
# if inp_ffproc.check_cut_syn:
#     util_plot.plot_cut_wave(all_stations=all_stations,
#                             waveforms=syn_waveforms,
#                             add_save=os.path.join(inp_ffproc.output_dir,
#                                                   all_events.name[ev]),
#                             fig_name=
#                             'wave_1_grfs_%s' % all_events.name[ev],
#                             ev_datetime=all_events.datetime[ev])

# =================================
# --------------- Convolve STF with REAL seismograms
# if inp_ffproc.conv_real_stf:
#     cprint('ffprocessing.py', '[WARNING]', bc.orange,
#            'Convolve "real" seismograms with STF(s)!')
#     real_waveforms, all_stfs = conv_stf(
#         address_stf=os.path.join(inp_ffproc.evproc_path,
#                                  all_events.name[ev]),
#         all_stations=all_stations,
#         syn_waveforms=real_waveforms,
#         inp_ffproc=inp_ffproc,
#         all_events=all_events,
#         ev=ev)
