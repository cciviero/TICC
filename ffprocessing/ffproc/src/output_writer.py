#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  output_writer.py
#   Purpose:   writing the outputs in the specified format
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

from src.utility_codes import bc, cprint
import src.utility_codes as util_plot

# ------------------ create_dir ---------------------------


def create_dir(output_dir):
    """
    create the top directory, rest of the information/figures will be stored
    in this directory.
    :param output_dir:
    :return:
    """
    #import ipdb; ipdb.set_trace()
    if not os.path.isdir(output_dir):
        os.makedirs(os.path.join(output_dir))
    else:
        cprint('output_writer.py', '[DIRECTORY]', bc.orange,
               '%s already exists. Continue.' % output_dir)

# ------------------ output_writer ---------------------------


def output_writer(inp_ffproc, all_stations, all_events, ev, all_stfs):
    """
    writing the outputs based on the selected format
    :param inp_ffproc:
    :param all_stations:
    :param all_events:
    :param ev:
    :param all_stfs:
    :return:
    """
    if inp_ffproc.output_format == 'original_ffprocessing':
        output_event_based(inp_ffproc, all_stations, all_events, ev, all_stfs)
    else:
        cprint('output_writer.py', '[OUTPUT] [FORMAT]', bc.dred,
               '%s has not implemented! EXIT!' % inp_ffproc.evproc_mode)
        sys.exit()

# ------------------ output_event_based ---------------------------


def output_event_based(inp_ffproc, all_stations, all_events, ev, all_stfs):
    """
    writing the outputs similar to MATLAB (RunFFprocessing)
    :param inp_ffproc:
    :param all_stations:
    :param all_events:
    :param ev:
    :param all_stfs:
    :return:
    """
    if not os.path.isdir(os.path.join(inp_ffproc.output_dir)):
        os.makedirs(os.path.join(inp_ffproc.output_dir))
    ev_name = all_events.name[ev]
    if not os.path.isdir(os.path.join(inp_ffproc.output_dir, ev_name)):
        os.makedirs(os.path.join(inp_ffproc.output_dir, ev_name))

    if not os.path.isdir(os.path.join(inp_ffproc.output_dir, ev_name,
                                      'outfiles')):
        os.makedirs(os.path.join(inp_ffproc.output_dir, ev_name, 'outfiles'))
    else:
        cprint('output_writer.py', '[DIRECTORY]', bc.orange,
               '%s already exists. Continue.'
               % os.path.join(inp_ffproc.output_dir, ev_name, 'outfiles'))

    for band in range(inp_ffproc.check_min_bands, inp_ffproc.check_max_bands):
        outfile_fio = open(os.path.join(inp_ffproc.output_dir,
                                        ev_name, 'outfiles',
                                        'ffproc.ampstt.band%02i'
                                        % (int(band)+1)), 'w')
        outfile_fio.writelines('# ffproc.ampstt\n')
        outfile_fio.writelines('# Event ID: %s\n' % all_events.name[ev])
        outfile_fio.writelines('# Event Time: %s\n' % all_events.datetime[ev])
        outfile_fio.writelines('# idx\tgrp\t'
                               'stla\tstlo\t'
                               'stazie\tstnam\t'
                               'xc_coeff\t'
                               'Tobs\t'
                               'dT\t'
                               'sigma_dT\t'
                               'not_used\tnot_used\t'
                               'A\tnot_used\t'
                               'not_used\tnot_used\t'
                               'tB_smgr\t'
                               'tB_mfi\t'
                               'winlen\t'
                               'clip_taumax\t'
                               'not_used\tSNR\t\t\t'
                               '[ts_step1\tts_step2\tcc_step1\tcc_step2\t'
                               'clip_step1\tclip_step2]\n')
        counter = 1
        cprint('output_writer.py', '[OUTPUT] [UNCERTAINITY]', bc.orange,
               "[OUTPUT] two_sigma is written for the uncertainty!")
        for cha in range(len(all_stations.name)):
            if inp_ffproc.mmeant_mode in ['CC-2step']:
                len_core = inp_ffproc.lgb_filt_duration[band] + \
                           inp_ffproc.ph_static_preset + \
                           inp_ffproc.ph_static_offset
            line_w = '%4d\t%i\t%3.2f\t%3.2f\t%3.2f\t%20s\t%1.2f\t%.2f\t%.2f\t%.2f\t' \
                     '%s\t%s\t%.3f\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%i\t' \
                     '%s\t%.3f\t\t\t%.3f\t%.3f\t%.3f\t%.3f\t%i\t%i\n' \
                     % (counter, all_stfs.grps[cha],
                        all_stations.lat[cha], all_stations.lon[cha],
                        all_stations.dist[cha], all_stations.name[cha],
                        all_stations.cc_step2[band][cha],
                        all_stations.tt_ph[cha] +
                        all_stations.final_time_shift[band][cha],
                        all_stations.final_time_shift[band][cha],
                        all_stations.final_two_sigma[band][cha],
                        '0.00', '0.00',
                        all_stations.final_amp[band][cha], '0.00',
                        '0.00', '0.00',
                        # This is similar to Tobs
                        # Note: the start time + static_preset!
                        all_stations.real_start[band+1][cha] +
                        inp_ffproc.ph_static_preset,
                        # This should be theoretical arrival time
                        all_stations.syn_start[band+1][cha] +
                        inp_ffproc.ph_static_preset,
                        # XXX len_core is used in raydata and raymatrix
                        # in raymatrix, it is being compared with detour time which is calculated
                        # from the arrival time of the seismic phase. 
                        len_core - inp_ffproc.ph_static_preset,
                        all_stations.clip_step2[band][cha],
                        '0.00', all_stations.snr[band][cha],
                        all_stations.ss_shift_step1[cha], all_stations.ss_shift_step2[band][cha],
                        all_stations.cc_step1[cha], all_stations.cc_step2[band][cha],
                        all_stations.clip_step1[cha], all_stations.clip_step2[band][cha])
            outfile_fio.writelines(line_w)
            counter += 1
        outfile_fio.close()

    # ======================= ffproc.source
    outfile_fio = open(os.path.join(inp_ffproc.output_dir, ev_name, 'outfiles',
                                    'ffproc.source'), 'w')
    outfile_fio.writelines('year, julianday, hr, min, sec, msec\n')
    outfile_fio.writelines('%i\t%i\t%i\t%i\t%i\t%s\n'
                           % (all_events.year[ev], all_events.julday[ev],
                              all_events.hr[ev], all_events.minu[ev],
                              all_events.sec[ev],
                              all_events.msec[ev]/1000.))
    outfile_fio.writelines('evlat, evlon, catalog_depth, inverted_depth\n')
    outfile_fio.writelines('%.2f\t%.2f\t%.2f\t%.2f\n'
                           % (all_events.lat[ev], all_events.lon[ev],
                              all_events.cat_dp[ev], all_events.inv_dp[ev]))
    outfile_fio.writelines('scalar_moment, tausrc\n')
    outfile_fio.writelines('%s\t%.2f\n'
                           % (all_events.scmom[ev], all_events.tau[ev]*2.))

    outfile_fio.writelines('xmom_cat(1:6) catalogue moment tensor\n')
    outfile_fio.writelines('%s\t%s\t%s\t%s\t%s\t%s\n'
                           % (all_events.mrr[ev], all_events.mtt[ev],
                              all_events.mpp[ev], all_events.mrt[ev],
                              all_events.mrp[ev], all_events.mtp[ev]))
    outfile_fio.writelines('strike, dip, rake\n')
    outfile_fio.writelines('%s\t%s\t%s\n'
                           % ('0.0', '0.0', '0.0'))
    outfile_fio.writelines('event catalog\n')
    outfile_fio.writelines('pyFFproc\n')
    outfile_fio.writelines('xmom_cat(1:6) inverted moment tensor\n')
    outfile_fio.writelines('%s\t%s\t%s\t%s\t%s\t%s\n'
                           % (all_events.mrr[ev], all_events.mtt[ev],
                              all_events.mpp[ev], all_events.mrt[ev],
                              all_events.mrp[ev], all_events.mtp[ev]))
    outfile_fio.close()

    # ======================= ffproc.receivers
    # XXX scales of elevation and depth should be meter!
    outfile_fio = open(os.path.join(inp_ffproc.output_dir, ev_name,
                                    'outfiles', 'ffproc.receivers'), 'w')
    outfile_fio.writelines('# ffproc.receivers\n')
    counter = 1
    for cha in range(len(all_stations.name)):
        line_w = '%i\t%i\t%s\t%.2f\t%.2f\t%.3f\t%.3f\t%s\n' \
                 % (counter, all_stfs.grps[cha], all_stations.name[cha],
                    all_stations.lat[cha], all_stations.lon[cha],
                    all_stations.el[cha], all_stations.dp[cha],
                    all_stations.name[cha])
        outfile_fio.writelines(line_w)
        counter += 1
    outfile_fio.close()

# ------------------ kernel_writer ---------------------------


def kernel_writer(inp_ffproc, all_events, ev, all_stations, all_stfs):
    """
    Write kernel inputs out of the measurements
    :param inp_ffproc:
    :param all_events:
    :param ev:
    :param all_stations:
    :return:
    """
    if not os.path.isdir(inp_ffproc.output_dir):
        os.makedirs(os.path.join(inp_ffproc.output_dir))
    else:
        cprint('output_writer.py', '[KERNEL]', bc.orange,
               '%s already exists. Continue' % inp_ffproc.output_dir)

    ev_name = all_events.name[ev]
    if not os.path.isdir(os.path.join(inp_ffproc.output_dir, ev_name)):
        os.makedirs(os.path.join(inp_ffproc.output_dir, ev_name))
    else:
        cprint('output_writer.py', '[KERNEL]', bc.orange,
               '%s already exists. Continue!'
               % os.path.join(inp_ffproc.output_dir, ev_name))
    if not os.path.isdir(os.path.join(inp_ffproc.output_dir,
                                      ev_name, 'outfiles')):
        os.makedirs(os.path.join(inp_ffproc.output_dir, ev_name, 'outfiles'))
    else:
        cprint('output_writer.py', '[KERNEL]', bc.orange,
               '%s already exists. Continue!'
               % os.path.join(inp_ffproc.output_dir, ev_name, 'outfiles'))

    kps = inp_ffproc.kernel_ps.upper()
    if kps == 'P':
        model_par = 'vp'
    else:
        cprint('output_writer.py', '[KERNEL]', bc.dred,
               '[KERNEL] %s has not implemented' % kps)
        sys.exit()

    kern_fio = open(os.path.join(inp_ffproc.output_dir, ev_name, 'outfiles',
                                 "kernel_%s.txt" % ev_name), 'w')
    kern_fio.writelines(inp_ffproc.kernel_comp + '\n')
    kernel_min_cc = inp_ffproc.kernel_min_cc

    count_sta = 0
    for i in range(len(all_stations.cc_step2[0])):
        acc_bands = all_stations.cc_step2[:, i] >= kernel_min_cc
        if not len(np.where(acc_bands)[0]) > 0:
            continue
        count_sta += 1
        if all_stations.name[i].split(".")[2].strip() == '':
            sta_name = '%s.--' % all_stations.name[i].split(".")[1]
        else:
            sta_name = '%s.%s' % (all_stations.name[i].split(".")[1],
                                  all_stations.name[i].split(".")[2])

        # XXX geocentric or geographic?
        first_line = '%s  %.2f  %.2f  stf%02i  %i\n'\
                     % (sta_name,
                        all_stations.geocen_lat[i],
                        all_stations.lon[i],
                        all_stfs.grps[i],
                        len(acc_bands[acc_bands > 0]))
        kern_fio.writelines(first_line)
        for j in np.where(acc_bands)[0]:
            if inp_ffproc.mmeant_mode in ['CC-2step-env']:
                len_core = all_stations.length_bp[j][i]
            elif inp_ffproc.mmeant_mode in ['CC-2step']:
                len_core = inp_ffproc.lgb_filt_duration[j] + \
                           inp_ffproc.ph_static_preset + \
                           inp_ffproc.ph_static_offset
            band = '%02i' % (j + 1)
            second_line = '%s_%s band%s %s %f %f %s\n' \
                          % (kps, band, band,
                             inp_ffproc.kernel_misfit,
                             all_stations.syn_start[j+1][i],
                             all_stations.syn_start[j+1][i] + len_core,
                             model_par)
            kern_fio.writelines(second_line)
    kern_fio.close()

    with open(os.path.join(inp_ffproc.output_dir,
                           ev_name, 'outfiles',
                           "kernel_%s.txt" % ev_name), 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(str(count_sta) + '\n' + content)
        f.close()
    
    filt_fio = open(os.path.join(inp_ffproc.output_dir,
                                 ev_name, 'outfiles',
                                 "kernel_filter.txt"), 'w')
    center_period = inp_ffproc.lgb_filt_pmax
    filt_fio.writelines('%s\n' % inp_ffproc.lgb_filt_nscale)
    for i in range(inp_ffproc.lgb_filt_nscale):
        band = '%02i' % (i + 1)
        second_line = 'band%s Gabor  %.2f %s 0.0 0.0\n' \
                      % (band, center_period, inp_ffproc.lgb_filt_sigmaIfc)
        filt_fio.writelines(second_line)
        center_period = center_period/inp_ffproc.lgb_filt_fmult
        filt_fio.writelines('')
    filt_fio.close()

    cmt_lines = []
    cmt_lines.append("XXX\n")
    cmt_lines.append("event name:     %s\n" % all_events.name[ev])
    cmt_lines.append("time shift:     XXX\n")
    cmt_lines.append("half duration:  XXX\n")
    # XXX geocentric or geographic?
    cmt_lines.append("latitude:      %f\n" % all_events.geocen_lat[ev]) 
    cmt_lines.append("longitude:     %f\n" % all_events.lon[ev]) 
    cmt_lines.append("depth:         %f\n" % all_events.inv_dp[ev])
    cmt_lines.append("Mrr:           %.5g\n" % all_events.mrr[ev]) 
    cmt_lines.append("Mtt:           %.5g\n" % all_events.mtt[ev])
    cmt_lines.append("Mpp:           %.5g\n" % all_events.mpp[ev])
    cmt_lines.append("Mrt:           %.5g\n" % all_events.mrt[ev])
    cmt_lines.append("Mrp:           %.5g\n" % all_events.mrp[ev])
    cmt_lines.append("Mtp:           %.5g\n" % all_events.mtp[ev])

    cmt_fio = open(os.path.join(inp_ffproc.output_dir,
                                 ev_name, 'outfiles',
                                 "CMTSOLUTION"), 'w')
    for cmt_l in cmt_lines:
        cmt_fio.writelines(cmt_l)
    cmt_fio.close()

# ------------------ write_arr2tr ---------------------------


def write_arr2tr(waveform, label, event_name,
                 inp_ffproc, all_stations, tr_format='SAC'):
    """
    Save arrays
    :param waveform:
    :param label:
    :param event_name:
    :param inp_ffproc:
    :param tr_format:
    :return:
    """
    cprint('ffprocessing.py', '[OUTPUT] [ARR2TR]', bc.cyan,
            'Saving %s broadband waveforms.' % label)
    save_address = os.path.join(inp_ffproc.output_dir, event_name, 'smgr')
    if not os.path.isdir(save_address):
        os.makedirs(save_address)

    for i, tr in enumerate(waveform):
        address_path = os.path.join(save_address,
                                    '%s.%s_%s.%s.%s.%s'
                                    % (label, event_name,
                                       tr.stats.network,
                                       tr.stats.station,
                                       tr.stats.location,
                                       tr.stats.channel))
        tr.write(address_path, format=tr_format)

# ------------------ write_cut_BB ---------------------------


def write_cut_BB(num_sta, inp_ffproc, cut_realBP, real_waveforms,
                 all_events, cut_synBP, syn_waveforms,
                 all_stfs, ev):
    """
    Write both cut and BB waveforms locally
    :param num_sta:
    :param inp_ffproc:
    :param cut_realBP:
    :param real_waveforms:
    :param all_events:
    :param cut_synBP:
    :param syn_waveforms:
    :param all_stfs:
    :param ev:
    :return:
    """
    cprint('ffprocessing.py', '[OUTPUT] [ARR2TR]', bc.cyan,
           'Writing array to traces.')
    for sta in range(num_sta):
        for band in range(inp_ffproc.check_min_bands,
                          inp_ffproc.check_max_bands):
            # XXX CHECK
            write_arr2tr(cut_realBP[band][sta], 
                         'cut_real_%02i' % band,
                         all_events.name[ev],
                         inp_ffproc=inp_ffproc)
            write_arr2tr(cut_synBP[band][sta][:],
                         'cut_syn_%02i' % band,
                         all_events.name[ev],
                         inp_ffproc=inp_ffproc)

    cprint('ffprocessing.py', '[OUTPUT] [STF]', bc.cyan,
           'Writing STF(s) to stream.')
    stf_save_address = os.path.join(inp_ffproc.output_dir,
                                    all_events.name[ev], 'smgr')
    stf_address_path = os.path.join(stf_save_address,
                                'STF.%s' % all_events.name[ev])
    all_stfs.stfs.write(stf_address_path, format="SAC")
