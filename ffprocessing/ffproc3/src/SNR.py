#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  SNR.py
#   Purpose:   Calculate SNR (Signal to Noise Ratio)
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
from obspy import Trace
import os

from src.utility_codes import bc, progress_bar, cprint

# ------------------ cSNR ---------------------------


def cSNR(realBP, waveform, inp_ffproc, all_stations, all_events, ev,
         n_tb=0.6, n_ta=0.1, s_tb=0.1, s_ta=0.6, method='squared',
         plot_ph_no=False):
    """
    'c'alculate 'S'ignal to 'N'oise 'R'atio
    calculate SNR, l1 noise level and l2 noise level
    :param realBP:
    :param waveform:
    :param inp_ffproc:
    :param all_stations:
    :param all_events:
    :param ev:
    :param n_tb:
    :param n_ta:
    :param s_tb:
    :param s_ta:
    :param method:
    :param plot_ph_no:
    :return:
    """
    for band in range(inp_ffproc.check_min_bands, inp_ffproc.check_max_bands):
        progress_bar(band, inp_ffproc.check_max_bands)
        for sta in range(all_stations.number_stations):
            if not inp_ffproc.mmeant_calc_snr:
                SNR = 0
            else:
                # trace information
                tr = Trace(realBP[:, sta, band])
                tr.stats.starttime = waveform[sta].stats.starttime
                tr.stats.network = waveform[sta].stats.network
                tr.stats.station = waveform[sta].stats.station
                tr.stats.location = waveform[sta].stats.location
                tr.stats.channel = waveform[sta].stats.channel
                tr.stats.sampling_rate = waveform[sta].stats.sampling_rate

                # finding the noise window
                tt_phase_time = all_events.datetime[ev] + \
                                all_stations.tt_ph[sta]
                smgr_starttime = waveform[sta].stats.starttime

                noise_window = tt_phase_time - smgr_starttime  # in [sec]
                n_before = noise_window * n_tb
                n_after = noise_window * n_ta

                # ====>>>> Slice Noise Window
                noise = tr.copy()
                noise = noise.slice(tt_phase_time - n_before,
                                    tt_phase_time - n_after)

                s_before = noise_window * s_tb
                s_after = noise_window * s_ta
                # ====>>>> Slice Signal Window
                signal = tr.copy()
                signal = signal.slice(tt_phase_time - s_before,
                                      tt_phase_time + s_after)

                if 'squared' in method.lower():
                    if np.sum(np.square(noise)) < 1e-30:
                        sum_square_noise = 1e-30
                    else:
                        sum_square_noise = np.sum(np.square(noise))

                # ====>>>> Calculate SNR
                try:
                    if 'squared' in method.lower():
                        SNR = np.sqrt(np.sum(np.square(signal)) /
                                      sum_square_noise) * \
                              (float(noise.stats.npts)/signal.stats.npts)
                    else:
                        cprint('cSNR.py', '[ERROR] [SNR]', bc.dred,
                               'Method: %s is not implemented.' % method)
                except Exception as e:
                    cprint('cSNR.py', '[ERROR] [SNR]', bc.dred,
                           'ERROR %s: %s \n'
                           'SNR set to "-12345"' % (tr.id, e))
                    SNR = -12345

            # ====>>>> plotting (takes a long time!)
            # XXX REVIEW
            if plot_ph_no:
                plt.clf()
                plt.figure(figsize=(23, 15))

                # plot the trace
                dt = tr.stats.delta
                npts = tr.stats.npts
                t = np.linspace(0., dt*(npts-1), npts)
                c_tr = tr.copy()
                plt.plot(t, c_tr.normalize(norm=max(tr)).data, 'b',
                         lw=3, alpha=0.6, label='real')

                # plot the phase signal
                dt = signal.stats.delta
                npts = signal.stats.npts
                t_start = signal.stats.starttime - tr.stats.starttime
                t = np.linspace(t_start, t_start+dt*(npts-1), npts)
                plt.plot(t, signal.normalize(norm=max(tr)).data,
                         'r', lw=3, alpha=1.0)

                # plot the noise
                dt = noise.stats.delta
                npts = noise.stats.npts
                t_start = noise.stats.starttime - tr.stats.starttime
                t = np.linspace(t_start, t_start+dt*(npts-1), npts)
                plt.plot(t, noise.normalize(norm=max(tr)).data,
                         'g', lw=3, alpha=0.8)

                # plotting the theoretical arrival time
                plt.axvline(noise_window, color='cyan', lw='3',
                            linestyle='dashed',
                            label='theoretical arrival time')

                plt.xlabel('Time [s]', size=18, weight='bold')
                plt.xticks(size=18, weight='bold')
                plt.yticks(size=18, weight='bold')
                plt.ylim(-1, 1)
                plt.legend()
                plt.title('SNR:%s -- %s' % (round(SNR, 4), tr.id),
                          size=25, weight='bold')

                save_address = os.path.join(inp_ffproc.output_dir,
                                            all_events.name[ev], 'snr')
                if not os.path.isdir(save_address):
                    os.makedirs(save_address)
                plt.savefig(os.path.join(save_address,
                                         '%s_%02i.png' % (tr.id, band+1)),
                            bbox_inches='tight')
                plt.clf()
                plt.close()

            # l1_noise = np.sum(np.abs(noise.data))/noise.stats.npts
            # l2_noise = np.sqrt(np.sum(np.square(noise.data))) / \
            #            noise.stats.npts
            all_stations.snr = np.append(all_stations.snr, SNR)

    all_stations.snr = np.reshape(all_stations.snr,
                                  (inp_ffproc.check_max_bands -
                                   inp_ffproc.check_min_bands, -1))

# ------------------ qSNR ---------------------------


def qSNR(waveform, tt_phase, all_events, ev, inp_ffproc,
         s_tb=5,
         s_ta=25,
         n_tb=0.10,
         n_ta=0.10,
         method='squared'):
    """
    'q'uick 'S'ignal to 'N'oise 'R'atio
    This calculator takes the waveform and calculates the SNR depending
    on the phase arrival time
    :param waveform:
    :param tt_phase:
    :param all_events:
    :param ev:
    :param inp_ffproc:
    :param s_tb:
    :param s_ta:
    :param n_tb:
    :param n_ta:
    :param method:
    :return:
    """
    # ====>>>> Define Time Parameters
    event_origin_time = all_events.datetime[ev]
    tt_phase_time = event_origin_time + tt_phase
    smgr_starttime = waveform.stats.starttime

    noise_window = tt_phase_time - smgr_starttime
    n_before = noise_window * n_tb
    n_after = noise_window * n_ta

    # ====>>>> Slice Noise Window
    noise = waveform.copy()
    noise = noise.slice(waveform.stats.starttime + n_before,
                        tt_phase_time - n_after)

    # ====>>>> Slice Signal Window
    signal = waveform.copy()
    signal = signal.slice(tt_phase_time - s_tb, tt_phase_time + s_ta)

    # ====>>>> Calculate SNR
    SNR = np.sqrt(np.sum(np.square(signal)) / np.sum(np.square(noise))) * \
                 (float(noise.stats.npts)/signal.stats.npts)

    return SNR, noise
