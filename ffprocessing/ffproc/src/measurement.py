#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  measurement.py
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
import matplotlib.pyplot as plt
import numpy as np
import sys

from utility_codes import bc, cprint

# ------------------ xcorrelation_2step ---------------------------


def xcorrelation_2step(realBP, synBP, realREF, synREF, inp_ffproc,
                       all_stations, real_waveforms, syn_waveforms,
                       all_events, ev):
    """
    cross-correlation in 2 step: first the waveforms are aligned in
    the lowest frequency band and then CC is applied for all other bands
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
    print '\n'
    cprint('measurement.py', '[MEASURE]', bc.lgreen,
           'cross-correlation with 2 steps!')
    # STEP 1
    xcorr_refband(realREF, synREF, inp_ffproc, all_stations)
    # STEP 2
    cut_realBP = []
    cut_synBP = []
    for band in range(inp_ffproc.check_min_bands, inp_ffproc.check_max_bands):
        temp_realBP, temp_synBP = xcorr_band(realBP, synBP,
                                             real_waveforms, syn_waveforms,
                                             inp_ffproc, all_stations, band,
                                             all_events, ev)
        cut_realBP.append(temp_realBP)
        cut_synBP.append(temp_synBP)

    # Collecting data and reshaping them
    all_stations.ss_shift_step2 = \
        np.reshape(all_stations.ss_shift_step2,
                   (inp_ffproc.check_max_bands -
                    inp_ffproc.check_min_bands, -1))
    all_stations.cc_step2 = \
        np.reshape(all_stations.cc_step2,
                   (inp_ffproc.check_max_bands -
                    inp_ffproc.check_min_bands, -1))
    all_stations.clip_step2 = \
        np.reshape(all_stations.clip_step2,
                   (inp_ffproc.check_max_bands -
                    inp_ffproc.check_min_bands, -1))
    all_stations.final_amp = \
        np.reshape(all_stations.final_amp,
                   (inp_ffproc.check_max_bands -
                    inp_ffproc.check_min_bands, -1))
    all_stations.final_two_sigma = \
        np.reshape(all_stations.final_two_sigma,
                   (inp_ffproc.check_max_bands -
                    inp_ffproc.check_min_bands, -1))
    all_stations.cut_sample_real = \
        np.reshape(all_stations.cut_sample_real,
                   (inp_ffproc.check_max_bands -
                    inp_ffproc.check_min_bands, -1))
    all_stations.cut_sample_syn = \
        np.reshape(all_stations.cut_sample_syn,
                   (inp_ffproc.check_max_bands -
                    inp_ffproc.check_min_bands, -1))

    if not np.shape(all_stations.ss_shift_step2)[1] == len(all_stations.name):
        cprint('measurement.py', '[MEASURE]', bc.dred,
               'Can not reshape the measurement results! EXIT!')
        sys.exit()

    # calculate the final time based on the samples and
    # the starttime of seismograms:
    cprint('measurement.py', '[MEASURE]', bc.yellow,
           'Calculating the measured time.')
    all_stations.final_time_shift = \
        all_stations.real_start[1:] - all_stations.syn_start[1:]

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
    # import ipdb; ipdb.set_trace()
    for sta in range(len(realREF)):
        # Find whether the measured time should be clipped
        max_allowed_shift = \
            inp_ffproc.mmeant_clip_time1*inp_ffproc.ph_sampling_rate
        sample_shift, cc_factor, two_sigma = \
            xcorr(realREF[sta], synREF[sta], maxlag=int(max_allowed_shift))
        if cc_factor == False or np.isnan(cc_factor):
            all_stations.ss_shift_step1 = \
                np.append(all_stations.ss_shift_step1, sample_shift)
            all_stations.clip_step1 = np.append(all_stations.clip_step1, 1)
            all_stations.cc_step1 = np.append(all_stations.cc_step1, 0)
            continue

        if abs(sample_shift) < abs(max_allowed_shift):
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
    wave_time = inp_ffproc.ph_static_preset +  \
                inp_ffproc.lgb_filt_duration[band] + \
                inp_ffproc.ph_static_offset
    wave_length = (inp_ffproc.ph_sampling_rate * wave_time) + 1
    wave_length = int(wave_length)

    all_stations.length_bp = \
        np.append(all_stations.length_bp, inp_ffproc.lgb_filt_duration[band])

    cut_realBP = []
    cut_synBP = []
    for sta in range(np.shape(realBP)[1]):
        # ------------------- synthetic data
        # slice the waveforms around the theoretical arrival time
        rel_syn_start_time = syn_waveforms[sta].stats.starttime - \
                             all_events.datetime[ev]
        time_to_slice = all_stations.tt_ph[sta] - \
                        rel_syn_start_time - \
                        inp_ffproc.ph_static_preset
        sample_to_slice = (inp_ffproc.ph_sampling_rate * time_to_slice) + 1
        sample_to_slice = int(sample_to_slice)

        all_stations.cut_sample_syn = \
            np.append(all_stations.cut_sample_syn, sample_to_slice)

        sub_synBP = cut_sample_waveform(waveform=synBP[:, sta, band],
                                        sample_to_slice=sample_to_slice,
                                        wave_length=wave_length)

        all_stations.syn_start[band+1, sta] += \
            (sample_to_slice - 1)/inp_ffproc.ph_sampling_rate

        # make sure that the arrays are in C format before cross-correlation
        sub_synBP = np.require(sub_synBP, dtype=np.float32,
                               requirements=['C'])

        # ------------------- real data
        rel_real_start_time = real_waveforms[sta].stats.starttime - \
                              all_events.datetime[ev]
        time_to_slice = all_stations.tt_ph[sta] - \
                        rel_real_start_time - \
                        inp_ffproc.ph_static_preset
        sample_to_slice = (inp_ffproc.ph_sampling_rate * time_to_slice) + 1
        # the 0 band is the best band for measuring the cc step 1...
        # Add the measured time shift from step-1 to the data
        sample_to_slice += all_stations.ss_shift_step1[sta]
        sample_to_slice = int(sample_to_slice)

        all_stations.cut_sample_real = \
            np.append(all_stations.cut_sample_real, sample_to_slice)

        sub_realBP = cut_sample_waveform(waveform=realBP[:, sta, band],
                                         sample_to_slice=sample_to_slice,
                                         wave_length=wave_length)

        all_stations.real_start[band+1, sta] += \
            (sample_to_slice - 1)/inp_ffproc.ph_sampling_rate

        # make sure that the arrays are in C format before cross-correlation
        sub_realBP = np.require(sub_realBP, dtype=np.float32,
                                requirements=['C'])

        # Find whether the measured time should be clipped
        max_allowed_shift = \
            inp_ffproc.mmeant_clip_time2*inp_ffproc.ph_sampling_rate

        sample_shift, cc_factor, two_sigma = \
            xcorr(sub_realBP, sub_synBP, maxlag=int(max_allowed_shift))
        if cc_factor == False or np.isnan(cc_factor):
            all_stations.ss_shift_step2 = \
                np.append(all_stations.ss_shift_step2, sample_shift)
            all_stations.clip_step2 = np.append(all_stations.clip_step2, 1)
            all_stations.cc_step2 = np.append(all_stations.cc_step2, 0)
            all_stations.final_two_sigma = \
                np.append(all_stations.final_two_sigma, -12345)
            all_stations.final_amp = \
                np.append(all_stations.final_amp, -12345)

            cut_realBP.append(sub_realBP)
            cut_synBP.append(sub_synBP)
            continue

        if abs(sample_shift) < abs(max_allowed_shift):
            all_stations.ss_shift_step2 = \
                np.append(all_stations.ss_shift_step2, sample_shift)
            all_stations.clip_step2 = np.append(all_stations.clip_step2, 0)
        else:
            all_stations.ss_shift_step2 = \
                np.append(all_stations.ss_shift_step2, max_allowed_shift)
            all_stations.clip_step2 = np.append(all_stations.clip_step2, 1)
            sample_shift = max_allowed_shift
        all_stations.cc_step2 = np.append(all_stations.cc_step2, cc_factor)
        all_stations.final_two_sigma = \
            np.append(all_stations.final_two_sigma, two_sigma)
        all_stations.real_start[band+1, sta] += \
            sample_shift/inp_ffproc.ph_sampling_rate

        # ------------- Calculating the amplitude
        # ------------------- real data
        rel_real_start_time = real_waveforms[sta].stats.starttime - \
                              all_events.datetime[ev]
        time_to_slice = all_stations.tt_ph[sta] - \
                        rel_real_start_time - \
                        inp_ffproc.ph_static_preset
        sample_to_slice = (inp_ffproc.ph_sampling_rate * time_to_slice) + 1
        # the 0 band is the best band for measuring the cc step 1...
        # Add the measured time shift from step-1 to the data
        sample_to_slice += all_stations.ss_shift_step1[sta]
        # For amplitude measurement, shift the trace based on the final CC
        # in addition to the other shifts
        sample_to_slice += all_stations.ss_shift_step2[sta]
        sample_to_slice = int(sample_to_slice)
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


def xcorr(real_inp, syn_inp, maxlag=10):
    """
    Cross correlation function
    :param real_inp:
    :param syn_inp:
    :param maxlag:
    :return:
    """
    if not len(real_inp) or not len(syn_inp):
        return 0, False, False
    real = np.copy(real_inp)
    syn = np.copy(syn_inp)
    try:
        corr = plt.xcorr(real, syn, maxlags=maxlag)
    except Exception, e:
        # This is because the max allowed time shift is
        # bigger then the trace itself
        if len(real_inp) < maxlag:
            corr = plt.xcorr(real, syn, maxlags=len(real_inp)-1)
        else:
            cprint('measurement.py', '[ERROR] ', bc.dred,
                   'Something is wrong with plt.xcorr. EXIT!')
            sys.exit()
    plt.clf()
    sample_shift = corr[0][np.argmax(corr[1])]
    cc_factor = np.max(corr[1])
    if abs(cc_factor) > 1.2:
        cprint('measurement.py', '[ATTENTION] [CC]', bc.red,
               'cc factor is %s, it will be set to zero!' % cc_factor)
        cc_factor = 0

    if cc_factor < 0.6:
        two_sigma = 10.0
    elif 0.6 <= cc_factor < 0.7:
        two_sigma = 2.0
    elif 0.7 <= cc_factor < 0.8:
        two_sigma = 1.5
    elif 0.8 <= cc_factor < 0.85:
        two_sigma = 1.0
    elif 0.85 <= cc_factor < 0.9:
        two_sigma = 0.6
    elif 0.9 <= cc_factor < 0.95:
        two_sigma = 0.4
    elif 0.95 <= cc_factor:
        two_sigma = 0.3
    else:
        two_sigma = 10.0
    return sample_shift, cc_factor, two_sigma



def xcorrf(data1, data2, shift=None, shift_zero=0, oneside=False,
           demean=True, window=0, ndat1d=0, ndat2d=0, N1=None, N2=None,
           normalize=True,
           freq_domain=False, transform_back=True,
           stdev1=None, stdev2=None):
    """
    FROM: https://github.com/trichter/sito/blob/master/xcorr.py
    Cross-correlation of numpy arrays data1 and data2 in frequency domain.
    
    We define cross-corelation as:
    xcorr[i] = sum_j (tr1[i+j-shift_zero] * tr2[j])
    The data is demeaned before cross-correlation and the result normalized
    after cross-correlation.
    data1, data2: data
    shift:    maximum samples to shift
              (window for i in the above formula)
    shift_zero: shift tr1 before cross-correlation by this amount of samples to
              the right (this means correlation function is shifted to the
              right or better: the window of what you get of the function
              is shifted to the left)
    oneside:  if True only the right/positive side of the correlation function
              is returned. Overrides parameter shift_zero.
    demean:   if True demean data beforehand
    normalize: if True normalize correlation function
              (1 means perfect correlation)
    window:   Use only data in this window for demeaning and normalizing
              0: window = min(ndat1, ndat2)
              >0: window = this parameter
    ndat1d, ndat2d: If >0 use different values for the length of the arrays when
              calculating the mean (defaults to window parameter)
    return:   numpy array with correlation function of length 2*shift+1 for
              oneside=False and of length shift+1 for oneside=True
    """
    if freq_domain and not transform_back:
        return data1 * np.conjugate(data2)
    elif freq_domain:
        min_size = max(2 * shift + 1 + abs(shift_zero),
                       (N1 + N2) // 2 + shift + abs(shift_zero))
        if len(data1) < min_size:
            raise ValueError('NFFT was not large enough to cover the desired '
                          'xcorr!\nnfft: %d, required minimum: %d' %
                          (len(data1), min_size))
        ret = (ifft(data1 * np.conjugate(data2))).real
    else:
        complex_result = (data1.dtype == np.complex or
                          data2.dtype == np.complex)
        N1 = len(data1)
        N2 = len(data2)
        #if isinstance(data1[0], np.integer) or isinstance(data2[0], np.integer):
        data1 = data1.astype('float64')
        data2 = data2.astype('float64')
        #if (N1-N2)%2==1:
        #    raise ValueError('(N1-N2)%2 has to be 0')
        if window == 0:
            window = min(N1, N2)
        if ndat1d == 0:
            ndat1d = window
        if ndat2d == 0:
            ndat2d = window
        # determine indices for demeaning and normalization
        ind1 = max(0, (N1 - window) // 2)
        ind2 = min(N1, (N1 + window) // 2)
        ind3 = max(0, (N2 - window) // 2)
        ind4 = min(N2, (N2 + window) // 2)

        # demean and normalize data
        if demean:
            data1 -= np.sum(data1[ind1:ind2]) / ndat1d
            data2 -= np.sum(data2[ind3:ind4]) / ndat2d
        if normalize:
            data1 /= np.max(data1[ind1:ind2])
            data2 /= np.max(data2[ind3:ind4])

        # Always use 2**n-sized FFT, perform xcorr
        size = max(2 * shift + 1 + abs(shift_zero),
                   (N1 + N2) // 2 + shift + abs(shift_zero))
        nfft = nextpow2(size)
        IN1 = fft(data1, nfft)
        if USE_FFTW3:
            IN1 = IN1.copy()
        IN1 *= np.conjugate(fft(data2, nfft))
        ret = ifft(IN1)
        if not USE_FFTW3:
            del IN1
        if not complex_result:
            ret = ret.real
    # shift data for time lag 0 to index 'shift'

    ret = np.roll(ret, -(N1 - N2) // 2 + shift + shift_zero)[:2 * shift + 1]
    # normalize xcorr
    if normalize:
        if not freq_domain:
            stdev1 = (np.sum(data1[ind1:ind2] ** 2)) ** 0.5
            stdev2 = (np.sum(data2[ind3:ind4] ** 2)) ** 0.5
#            stdev1 = (np.sum(data1 ** 2)) ** 0.5
#            stdev2 = (np.sum(data2 ** 2)) ** 0.5
        if stdev1 == 0 or stdev2 == 0:
            log.warning('Data is zero!!')
            ret[:] = 0.
        else:
            ret /= stdev1 * stdev2
    if oneside:
        ret = ret[shift:]
    return np.copy(ret)


# # # ------------------- xcorr --------------------
# #
# #
# # def xcorr(real_inp, syn_inp):
# #     """
# #     Cross correlation function
# #     :param real_inp:
# #     :param syn_inp:
# #     :return:
# #     """
# #     if not len(real_inp) or not len(syn_inp):
# #         return -12345, False, False, False
# #     real = np.copy(real_inp)
# #     syn = np.copy(syn_inp)
# #     # real -= np.mean(real)
# #     # syn -= np.mean(syn)
# #     corr = np.correlate(real, syn, 'full')
# #     sample_shift = np.argmax(corr) - (len(corr)+1)/2 + 1
# #     cc_factor = np.max(corr)/np.sqrt((real ** 2).sum() * (syn ** 2).sum())
# #     if abs(cc_factor) > 1.2:
# #         cprint('measurement.py', '[ATTENTION][CC]', bc.red,
# #                'cc factor is %s, it will be set to zero!' % cc_factor)
# #         cc_factor = 0
# #     two_sigma = np.std(corr)/np.sqrt((real ** 2).sum() * (syn ** 2).sum())
# #     return sample_shift, cc_factor, corr, two_sigma
