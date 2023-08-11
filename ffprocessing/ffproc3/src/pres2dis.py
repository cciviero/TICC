#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  pres2dis.py
#   Purpose:   function for converting pressure smgr to displacement;
#              specifically for using hydrophone data for cross corr
#   Author:    Maria Tsekhmistrenko, Kasra Hosseini
#   Email:     mariat@earth.ox.ac.uk
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
# from filtering import cut_sample_waveform

import matplotlib.pylab as plt
import numpy as np
import glob
import os
import sys

# # --------------------- dis2pres ---------------------
#
# def dis2pres(tr_dis, all_events, ev):
#     """
#     old function that I would use on the synthetic data, but after
#     a conversation with G.Nolet the other way round makes more sense
#
#
#     calculates the pressure smgr from displacement
#     Instaseis cannot generate pressure smgr;
#     recalculate DIS to pressure after SIMONS ET AL 2009/REID ET AL
#     :param tr_dis:
#     :param all_events:
#     :param ev:
#     :return:
#     """
#     # constants according to Simons et al. 2009.
#     density_w = 1000 # kg/m3
#     alpha = 5800 # m/s
#     omega_large = 2 * math.pi  # Hz
#     omega_small = 2 * math.pi * 5  # Hz
#
#     if all_events.mag[ev] > 5:
#         tr_dis.differentiate()
#         tr_dis.data = tr_dis.data * density_w * alpha * omega_large
#     else:
#         tr_dis.differentiate()
#         tr_dis.data = tr_dis.data * density_w * alpha * omega_small
#     return tr_dis


# --------------------- pres2dis ---------------------

def pres2dis(all_stations, real_waveforms, all_events,
                      ev, inp_ffproc, plot_p2d=True):
    """
    calculates the pressure smgr to displacement because Instaseis does have pressure data
    -> recalculate PRES to DIS after Simons et al. 2009/Reid et al.
    :param tr_pres:
    :param all_events:
    :param ev:
    :return:
    """

    # constants according to Simons et al. 2009.
    density_w = 1000 # kg/m3
    alpha = 5800 # m/s
    omega_large = 1 # Hz
    omega_small = 5 # Hz

    for indx, station in enumerate(real_waveforms):
        # import ipdb; ipdb.set_trace()
        # XXX recheck on omega large an small again in the paper, where does this

        original = real_waveforms[indx].copy()
        real_waveforms[indx].integrate()

        if all_events.mag[ev] > 5:
            real_waveforms[indx].data = real_waveforms[indx].data / density_w / alpha / omega_large
        else:
            real_waveforms[indx].data = real_waveforms[indx].data / density_w / alpha / omega_small

        # BDZ to symbolize that we changed it from BDH to displacement
        real_waveforms[indx].stats['channel'] = 'BDZ'
        xaxis = real_waveforms[indx].times()

        if plot_p2d:
            plot_pres2dis(original.data, real_waveforms[indx].data, xaxis,
                          all_stations.name[indx], all_events.name[ev], inp_ffproc)


def plot_pres2dis(original, calculated, xaxis, station_name, event_name, inp_ffproc):

    output_dir = os.path.join(inp_ffproc.output_dir, event_name,
                              'pres2dis')

    file_name = '%s_pres2dis.png' % station_name

    if not os.path.isdir(output_dir):
        os.makedirs(os.path.join(output_dir))

    plt.figure()
    plt.ioff()
    plt.plot(xaxis, original / max(abs(original)), c='r', lw=2, label='pressure')
    plt.plot(xaxis, calculated / max(abs(calculated)), c='b', lw=2, label='calculated displacement')

    plt.legend(fontsize=8)
    plt.ylim(-1, 1)

    plt.title('%s: Pressure to Displacment (normalized individually)' % station_name, weight='bold', fontsize=15)
    plt.savefig(os.path.join(output_dir, file_name), format='png')
    plt.clf()
    plt.close()

