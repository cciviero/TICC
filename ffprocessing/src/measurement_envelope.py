#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  measurement_envelope.py
#   Purpose:   Collection of utility tools for measurement
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

# Required Python modules will be imported in this part.
from filtering import cut_sample_waveform
import coloredlogs
import logging
import numpy as np
import sys

from obspy.signal.filter import envelope

coloredlogs.install()
"""
Useful commands:
log.debug("A quirky message only developers care about")
log.info("Curious users might want to know this")
log.warn("Something is wrong and any user should be informed")
log.error("Serious stuff, this is red for a reason")
log.critical("OH NO everything is on fire")
"""

# ------------------ xcorrelation_2step ---------------------------


def xcorrelation_2step_env(realBP, synBP, realREF, synREF, inp_ffproc,
                           all_stations, real_waveforms, syn_waveforms,
                           all_events, ev):
    """
    cross-correlation in 2 step: first the waveforms are aligned in
    the lowest frequency band and then CC is applied for all others
    :param realBP:
    :param synBP:
    :param realREF:
    :param synREF:
    :param inp_ffproc:
    :param all_stations:
    :param real_waveforms:
    :param syn_waveforms:
    :param all_events:
    :param ev:
    :return:
    """
    print('\n')
    logging.info('[MEASURE] cross-correlation with 2 steps -- envelope!')
    # STEP 1
    xcorr_refband(realREF, synREF, inp_ffproc, all_stations)
    # STEP 2
    cut_realBP = []
    cut_synBP = []
    for band in range(inp_ffproc.check_min_bands, inp_ffproc.check_max_bands):
        temp_realBP, temp_synBP = \
            xcorr_band(realBP, synBP,
                       real_waveforms, syn_waveforms,
                       inp_ffproc, all_stations, band, all_events, ev)
        cut_realBP.append(temp_realBP)
        cut_synBP.append(temp_synBP)

    # Collecting data and reshaping them
    all_stations.ss_shift_step2 = np.reshape(all_stations.ss_shift_step2,
                                             (inp_ffproc.check_max_bands -
                                              inp_ffproc.check_min_bands, -1))
    all_stations.cc_step2 = np.reshape(all_stations.cc_step2,
                                       (inp_ffproc.check_max_bands -
                                        inp_ffproc.check_min_bands, -1))
    all_stations.clip_step2 = np.reshape(all_stations.clip_step2,
                                         (inp_ffproc.check_max_bands -
                                          inp_ffproc.check_min_bands, -1))
    all_stations.final_amp = np.reshape(all_stations.final_amp,
                                        (inp_ffproc.check_max_bands -
                                         inp_ffproc.check_min_bands, -1))
    all_stations.final_two_sigma = np.reshape(all_stations.final_two_sigma,
                                              (inp_ffproc.check_max_bands -
                                               inp_ffproc.check_min_bands, -1))
    all_stations.cut_sample_real = np.reshape(all_stations.cut_sample_real,
                                              (inp_ffproc.check_max_bands -
                                               inp_ffproc.check_min_bands, -1))
    all_stations.cut_sample_syn = np.reshape(all_stations.cut_sample_syn,
                                             (inp_ffproc.check_max_bands -
                                              inp_ffproc.check_min_bands, -1))
    all_stations.length_bp = np.reshape(all_stations.length_bp,
                                        (inp_ffproc.check_max_bands -
                                         inp_ffproc.check_min_bands, -1))

    if not np.shape(all_stations.ss_shift_step2)[1] == len(all_stations.name):
        sys.exit('[MEASURE] can not reshape the measurement results')

    # calculate the final time based on the samples and
    # the starttime of seismograms:
    logging.info('[MEASURE] Calculating the measured time')

    all_stations.final_time_shift = \
        all_stations.real_start[1:] - all_stations.syn_start[1:]
    # We do not need this anymore as we change the real_start
    # all_stations.ss_shift_step2/all_stations.srate[0]

    return cut_realBP, cut_synBP

# ------------------- xcorr_refband --------------------


def xcorr_refband(realREF, synREF, inp_ffproc, all_stations):
    """
    First step in CC-2step method: aligning the waveforms at the lowest
    frequency band
    :param realREF:
    :param synREF:
    :param inp_ffproc:
    :param all_stations:
    :return:
    """
    for sta in range(len(realREF)):
        sample_shift, cc_factor, corr, two_sigma = \
            xcorr(realREF[sta], synREF[sta])
        # Following case happens if there is something wrong with the
        # cross-correlation
        if cc_factor == False or np.isnan(cc_factor):
            all_stations.ss_shift_step1 = \
                np.append(all_stations.ss_shift_step1, sample_shift)
            all_stations.clip_step1 = np.append(all_stations.clip_step1, 1)
            all_stations.cc_step1 = np.append(all_stations.cc_step1, 0)
            continue
        # Find whether the measured time should be clipped
        max_allowed_shift = \
            inp_ffproc.mmeant_clip_time1*inp_ffproc.ph_sampling_rate + 1
        if sample_shift <= max_allowed_shift:
            all_stations.ss_shift_step1 = \
                np.append(all_stations.ss_shift_step1, sample_shift)
            all_stations.clip_step1 = np.append(all_stations.clip_step1, 0)
        else:
            all_stations.ss_shift_step1 = \
                np.append(all_stations.ss_shift_step1, max_allowed_shift)
            all_stations.clip_step1 = np.append(all_stations.clip_step1, 1)
        all_stations.cc_step1 = np.append(all_stations.cc_step1, cc_factor)

# ------------------- xcorr_band --------------------


def xcorr_band(realBP, synBP, real_waveforms, syn_waveforms, inp_ffproc,
               all_stations, band, all_events, ev):
    """
    Second step in CC-2step: fine alignment of the waveforms which have been
    already shifted based on the lowest frequency
    :param realBP:
    :param synBP:
    :param real_waveforms:
    :param syn_waveforms:
    :param inp_ffproc:
    :param all_stations:
    :param band:
    :param all_events:
    :param ev:
    :return:
    """
    # wave length is defined based on the IR of the filter
    wave_time_IR = inp_ffproc.ph_static_preset +  \
                   inp_ffproc.lgb_filt_duration[0] + \
                   inp_ffproc.ph_static_offset
    wave_length_IR = (inp_ffproc.ph_sampling_rate * wave_time_IR) + 1
    wave_length_IR = int(wave_length_IR)

    cut_realBP = []
    cut_synBP = []
    for sta in range(np.shape(realBP)[1]):
        cont_proc = False
        # ------------------- synthetic data
        # Calculate number of samples from Karin's method
        # This will be used for checking purposes in the next steps
        rel_syn_start_time = syn_waveforms[sta].stats.starttime - \
                             all_events.datetime[ev]
        time_to_slice = all_stations.tt_ph[sta] - \
                        rel_syn_start_time
        time_to_slice += inp_ffproc.tshift[band]
        sample_check = inp_ffproc.ph_sampling_rate * time_to_slice + 1
        sample_check = int(sample_check)

        # Envelope of the synthetic waveform
        # syn_env = envelope(synBP[:, sta, band])
        # start_ttt = sample_check - \
        #             inp_ffproc.ph_static_preset * inp_ffproc.ph_sampling_rate
        # syn_env[0:start_ttt] = 0
        # syn_env[wave_length_IR+start_ttt:] = 0
        # syn_env_norm = syn_env/max(syn_env)
        # import ipdb; ipdb.set_trace()
        syn_env = np.sum(abs(synBP[:, :, band]), axis=1)
        syn_env_norm = syn_env/max(syn_env)

        # XXX INPUT: thres, check_iter, 0.05, 50
        thres = 0.2
        check_iter = 1
        while check_iter > 0 or thres > 0:
            check_iter -= 1
            env_sel = np.where(syn_env_norm > thres)[0]
            # To determine different groups in envelope
            diff_env_sel = env_sel[:-1] - env_sel[1:]
            # First passed group as the beginning time
            btime = env_sel[0]

            group_finder = np.where(diff_env_sel < -2)[0]
            if len(group_finder) > 0:
                etime = btime + group_finder[0]
            else:
                etime = btime + len(env_sel)
            # # Check whether the selected group is reasonable!
            # acceptable_sshift = inp_ffproc.mmeant_clip_time1 * \
            #                     inp_ffproc.ph_sampling_rate
            # if abs(btime - sample_check) > 5*acceptable_sshift:
            #     thres -= 0.05
            #     continue

            # if (etime - btime)/inp_ffproc.ph_sampling_rate > 50:
            #     thres -= 0.05
            #     continue
            cont_proc = True
            break

        all_stations.length_bp = \
            np.append(all_stations.length_bp, etime - btime)
        btime -= inp_ffproc.ph_static_preset*inp_ffproc.ph_sampling_rate
        etime += inp_ffproc.ph_static_offset*inp_ffproc.ph_sampling_rate
        wave_length = etime - btime

        sample_to_slice = btime
        all_stations.cut_sample_syn = \
            np.append(all_stations.cut_sample_syn, sample_to_slice)
        sub_synBP = cut_sample_waveform(waveform=synBP[:, sta, band],
                                        sample_to_slice=sample_to_slice,
                                        wave_length=wave_length)
        all_stations.syn_start[band+1, sta] += \
            (sample_to_slice - 1)/inp_ffproc.ph_sampling_rate
        # make sure that the arrays are in C format before cross-correlation
        sub_synBP = np.require(sub_synBP, dtype=np.float32, requirements=['C'])

        # ------------------- real data
        # the 0 band is the best band for measuring the cc step 1...
        # Add the measured time shift from step-1 to the data
        sample_to_slice = all_stations.ss_shift_step1[sta]
        sample_to_slice = int(sample_to_slice) + btime
        wave_length = etime - btime
        all_stations.cut_sample_real = \
            np.append(all_stations.cut_sample_real, sample_to_slice)
        sub_realBP = cut_sample_waveform(waveform=realBP[:, sta, band],
                                         sample_to_slice=sample_to_slice,
                                         wave_length=wave_length)
        all_stations.real_start[band+1, sta] += \
            (sample_to_slice - 1)/inp_ffproc.ph_sampling_rate
        # make sure that the arrays are in C format before cross-correlation
        sub_realBP = np.require(sub_realBP,
                                dtype=np.float32,
                                requirements=['C'])

        sample_shift, cc_factor, corr, two_sigma = xcorr(sub_realBP, sub_synBP)
        # Following case happens if there is something wrong with the
        # cross-correlation
        if cc_factor == False or np.isnan(cc_factor) or cont_proc == False:
            all_stations.ss_shift_step2 = \
                np.append(all_stations.ss_shift_step2, sample_shift)
            all_stations.clip_step2 = np.append(all_stations.clip_step2, 1)
            all_stations.cc_step2 = np.append(all_stations.cc_step2, 0)
            all_stations.final_two_sigma = \
                np.append(all_stations.final_two_sigma, -12345)
            calc_amp = -12345
            all_stations.final_amp = np.append(all_stations.final_amp, calc_amp)

            cut_realBP.append(sub_realBP)
            cut_synBP.append(sub_synBP)
            continue

        # Find whether the measured time should be clipped
        max_allowed_shift = \
            inp_ffproc.mmeant_clip_time2*inp_ffproc.ph_sampling_rate + 1
        if sample_shift <= max_allowed_shift:
            all_stations.ss_shift_step2 = \
                np.append(all_stations.ss_shift_step2, sample_shift)
            all_stations.clip_step2 = np.append(all_stations.clip_step2, 0)
        else:
            all_stations.ss_shift_step2 = \
                np.append(all_stations.ss_shift_step2, max_allowed_shift)
            all_stations.clip_step2 = np.append(all_stations.clip_step2, 1)
        all_stations.cc_step2 = np.append(all_stations.cc_step2, cc_factor)
        all_stations.final_two_sigma = \
            np.append(all_stations.final_two_sigma, two_sigma)
        all_stations.real_start[band+1, sta] += \
            sample_shift/inp_ffproc.ph_sampling_rate

        # ------------- Calculating the amplitude
        # ------------------- real data
        # the 0 band is the best band for measuring the cc step 1...
        # Add the measured time shift from step-1 to the data
        sample_to_slice = all_stations.ss_shift_step1[sta]
        # For amplitude measurement, shift the trace based on the final CC
        # in addition to the other shifts
        sample_to_slice += all_stations.ss_shift_step2[sta]
        sample_to_slice = int(sample_to_slice) + btime
        wave_length = etime - btime
        sub_realBP_amp = cut_sample_waveform(waveform=realBP[:, sta, band],
                                             sample_to_slice=sample_to_slice,
                                             wave_length=wave_length)
        if (len(sub_realBP_amp) == 0) or (len(sub_synBP) == 0):
            calc_amp = -12345
        else:
            calc_amp = \
                np.sum(sub_synBP*sub_realBP_amp)/np.sum(sub_synBP*sub_synBP)
        all_stations.final_amp = np.append(all_stations.final_amp, calc_amp)

        cut_realBP.append(sub_realBP)
        cut_synBP.append(sub_synBP)

    return cut_realBP, cut_synBP

# ------------------- xcorr --------------------


def xcorr(real_inp, syn_inp):
    """
    Cross correlation function
    :param real_inp:
    :param syn_inp:
    :return:
    """
    if not len(real_inp) or not len(syn_inp):
        return -12345, False, False, False
    real = np.copy(real_inp)
    syn = np.copy(syn_inp)
    # real -= np.mean(real)
    # syn -= np.mean(syn)
    corr = np.correlate(real, syn, 'full')
    sample_shift = np.argmax(corr) - (len(corr)+1)/2 + 1
    cc_factor = np.max(corr)/np.sqrt((real ** 2).sum() * (syn ** 2).sum())
    if abs(cc_factor) > 1.2:
        logging.error('[CC] cc factor is %s, it will be set to zero!'
                      % cc_factor)
        cc_factor = 0
        two_sigma = 1000
    else:
        two_sigma = np.std(corr)/np.sqrt((real ** 2).sum() * (syn ** 2).sum())
    return sample_shift, cc_factor, corr, two_sigma
