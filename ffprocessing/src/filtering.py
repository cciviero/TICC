#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  filtering.py
#   Purpose:   Collection of tools for filtering
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
from obspy import Trace
import os
import sys

from src.utility_codes import bc, cprint

# ------------------ filter_duration_calc ---------------------------


def filter_duration_calc(inp_ffproc):
    """
    Calculating the Impulse Response of the filters
    (used for defining the time window)
    :param inp_ffproc:
    :return:
    """
    # Dummy data to determine the IR for each filter
    DATdummy = np.array([1, 0, 0])
    trdummy = Trace(DATdummy)
    trdummy.stats.sampling_rate = inp_ffproc.ph_sampling_rate
    ftype = inp_ffproc.filt_mode
    npad = inp_ffproc.lgb_filt_npad
    nscale = inp_ffproc.lgb_filt_nscale
    pmax = inp_ffproc.lgb_filt_pmax
    fmult = inp_ffproc.lgb_filt_fmult
    sigmaIfc = inp_ffproc.lgb_filt_sigmaIfc
    # calculating the Impulse Response
    [tmp, tmp, IR, tshift, lngabor, f] = \
        bandpassfilter(trdummy, ftype, npad, nscale, pmax,
                       fmult, sigmaIfc, calc_IR=True)
    inp_ffproc.tshift = tshift
    inp_ffproc.lgb_filt_IR = IR
    center_period = inp_ffproc.lgb_filt_pmax
    # calculating the filter duration based on IR
    for scal in range(inp_ffproc.lgb_filt_nscale):
        tstmp = np.cumsum(np.power(IR[:, scal], 2))
        tstmp = tstmp/tstmp[-1]
        tstmp = tstmp - inp_ffproc.lgb_filt_energy_frac
        inp_ffproc.lgb_filt_duration.append(
            (np.argmin(abs(tstmp))+1)*(1./inp_ffproc.ph_sampling_rate))
        inp_ffproc.lgb_filt_center_period.append(center_period)
        center_period = round(center_period / inp_ffproc.lgb_filt_fmult, 2)

    inp_ffproc.lgb_filt_duration = np.array(inp_ffproc.lgb_filt_duration)

    # ------------ plotting the filters
    plt.ioff()
    plt.figure(figsize=(10, 10))
    cmap = plt.get_cmap('jet_r')
    colors = cmap(np.linspace(0, 1.0, np.shape(IR)[1]))

    plt.subplot(211)
    for i in range(np.shape(IR)[1]):
        plt.plot(
            np.linspace(0, (len(IR)-1)/inp_ffproc.ph_sampling_rate, len(IR)),
            IR[:, i],
            label=inp_ffproc.lgb_filt_duration[i],
            color=colors[i])
    plt.xlabel('Time [sec]', size=15, weight='bold')
    plt.xlim(0, 100)
    plt.title('Impulse Response of Gabor filters', size=18, weight='bold')
    plt.legend(fontsize=9, title='Filter duration')

    plt.subplot(212)
    for i in range(np.shape(lngabor)[1]):
        plt.plot(f, lngabor[:, i],
                 label=inp_ffproc.lgb_filt_center_period[i],
                 color=colors[i])
        plt.vlines(1/inp_ffproc.lgb_filt_center_period[i],
                   ymin=0, ymax=1,
                   color='k', linestyles='dashed', lw=1)

    plt.xscale('log')
    plt.xlim(0.5E-2, 1E1)
    plt.xlabel('Frequency [1/sec] (logarithmic)', size=15, weight='bold')
    plt.title('Log Gabor Filters', size=18, weight='bold')
    plt.legend(fontsize=9, title='Dominant period')

    plt.tight_layout()
    plt.savefig(os.path.join(inp_ffproc.output_dir, 'gabor_filts.png'),
                format='png', bbox_inches='tight')
    plt.clf()
    plt.ioff()
    plt.close()

    # raydata_raymatrix, resample the gabor filters
    raydata_raymatrix_bpf_input(IR, lngabor, f, inp_ffproc)

    # ------- multiplying the filter duration with nlambda
    inp_ffproc.lgb_filt_duration *= inp_ffproc.lgb_filt_nlambda
    cprint('filtering.py', '[FILTER]', bc.yellow,
           'selected nlambda: %s\nselected window lengths: %s'
           % (inp_ffproc.lgb_filt_nlambda, inp_ffproc.lgb_filt_duration))

# ------------------ bandpass_gabor ---------------------------


def bandpass_gabor(real_waveforms, syn_waveforms, inp_ffproc, all_stations,
                   all_events, ev):
    """
    bandpass real and synthetic data
    :param real_waveforms:
    :param syn_waveforms:
    :param inp_ffproc:
    :param all_stations:
    :param all_events:
    :param ev:
    :return:
    """
    cprint('filtering.py', '[BANDPASS]', bc.yellow, 'log-gabor filter')
    ftype = inp_ffproc.filt_mode
    npad = inp_ffproc.lgb_filt_npad
    nscale = inp_ffproc.lgb_filt_nscale
    pmax = inp_ffproc.lgb_filt_pmax
    fmult = inp_ffproc.lgb_filt_fmult
    sigmaIfc = inp_ffproc.lgb_filt_sigmaIfc

    # -------------- bandpass filter both real and synthetic waveforms
    # accurate duration of the seismogram
    # seismogram_length = preset + filter length (band0) + offset
    acc_time = inp_ffproc.ph_preset + \
               inp_ffproc.lgb_filt_duration[0] + \
               inp_ffproc.ph_offset
    acc_dur = (inp_ffproc.ph_sampling_rate * acc_time) + 1
    # acc_dur = seismogram_length + zero_padding
    acc_dur = acc_dur + inp_ffproc.lgb_filt_npad
    acc_dur = int(acc_dur)

    realBP = np.array([])
    synBP = np.array([])

    # print len(real_waveforms[0].data), np.shape(realBP)

    # syn_start = starttime of the trace - event datetime
    tmp_syn_start = all_stations.syn_start
    tmp_real_start = all_stations.real_start
    # Expand the syn and real start for all bands
    for i in range(inp_ffproc.check_min_bands, inp_ffproc.check_max_bands):
        all_stations.real_start = np.vstack([all_stations.real_start,
                                             tmp_real_start])
        all_stations.syn_start = np.vstack([all_stations.syn_start,
                                            tmp_syn_start])

    for sta in range(len(real_waveforms)):
        # progress_bar(curr_indx=sta, total_number=len(real_waveforms))
        # -------------- real data
        tr_real = real_waveforms[sta]
        tr_real.taper(0.05, 'cosine')
        [tmp_realBP, TAX_real, IR, tshift, lngabor, f] = \
            bandpassfilter(tr_real, ftype, npad, nscale, pmax, fmult, sigmaIfc)
        tmp_realBP = tmp_realBP[:acc_dur]
        if len(realBP) == 0:
            realBP = tmp_realBP
        else:
            realBP = np.append(realBP, tmp_realBP, 1)

        # -------------- synthetic data
        try:
            tr_syn = syn_waveforms[sta]
            tr_syn.taper(0.05, 'cosine')
            [tmp_synBP, TAX_syb, IR, tshift, lngabor, f] = \
                bandpassfilter(tr_syn, ftype, npad, nscale, pmax, fmult, sigmaIfc)
            tmp_synBP = tmp_synBP[:acc_dur]
           #  print np.shape(tmp_synBP)
        except Exception as e2:
            print('e2', e2)

            # XXX HACK for when we don't have horizontal synthetic available 
            tmp_synBP =  np.zeros((acc_dur, 1, 8))
        try:
            if len(synBP) == 0:
                synBP = tmp_synBP
            else:
                synBP = np.append(synBP, tmp_synBP, 1)
        except Exception as e3:
            print('e3', e3)
            import ipdb; ipdb.set_trace()

    # -------------- create reference seismograms
    wave_time = inp_ffproc.ph_static_preset + \
                inp_ffproc.lgb_filt_duration[0] + \
                inp_ffproc.ph_static_offset
    wave_length = (inp_ffproc.ph_sampling_rate * wave_time) + 1
    wave_length = int(wave_length)

    realREF = np.array([])
    for sta in range(len(real_waveforms)):
        real_start_time = real_waveforms[sta].stats.starttime - \
                          all_events.datetime[ev]
        # import ipdb; ipdb.set_trace()
        if not abs(real_start_time - all_stations.real_start[0][sta]) < \
                (1./inp_ffproc.ph_sampling_rate/10.):
            cprint('filtering.py', '[ERROR]', bc.dred,
                   'Start-times Mismatch for real data: %s. EXIT!' % sta)
            sys.exit()

        time_to_slice = all_stations.tt_ph[sta] - \
                        real_start_time - \
                        inp_ffproc.ph_static_preset
        sample_to_slice = (inp_ffproc.ph_sampling_rate * time_to_slice) + 1
        sample_to_slice = int(sample_to_slice)

        tmp_realREF = cut_sample_waveform(waveform=realBP[:, sta, 0],
                                          sample_to_slice=sample_to_slice,
                                          wave_length=wave_length)

        # Adjust the real_start
        all_stations.real_start[0, sta] += \
            (sample_to_slice - 1)/inp_ffproc.ph_sampling_rate
        if len(realREF) == 0:
            realREF = tmp_realREF
        else:
            realREF = np.vstack([realREF, tmp_realREF])
    # print len(real_waveforms[0].data), np.shape(realBP)
    synREF = np.array([])
    for sta in range(len(syn_waveforms)):
        syn_start_time = syn_waveforms[sta].stats.starttime - \
                         all_events.datetime[ev]
        if not abs(syn_start_time - all_stations.syn_start[0][sta]) < \
                (1./inp_ffproc.ph_sampling_rate/10.):
            cprint('filtering.py', '[ERROR]', bc.dred,
                   'Start-times Mismatch for syn data: %s. EXIT!' % sta)
            sys.exit()

        time_to_slice = all_stations.tt_ph[sta] - \
                        syn_start_time - \
                        inp_ffproc.ph_static_preset
        sample_to_slice = (inp_ffproc.ph_sampling_rate * time_to_slice) + 1
        sample_to_slice = int(sample_to_slice)

        tmp_synREF = cut_sample_waveform(waveform=synBP[:, sta, 0],
                                         sample_to_slice=sample_to_slice,
                                         wave_length=wave_length)

        # Adjust the syn_start
        all_stations.syn_start[0, sta] += \
            (sample_to_slice - 1)/inp_ffproc.ph_sampling_rate
        if len(synREF) == 0:
            synREF = tmp_synREF
        else:
            synREF = np.vstack([synREF, tmp_synREF])
    # print len(real_waveforms[0].data), np.shape(realBP)
    if len(np.shape(synREF)) == 1:
        synREF = np.reshape(synREF, [1, len(synREF)])
    if len(np.shape(realREF)) == 1:
        realREF = np.reshape(realREF, [1, len(realREF)])

    return realBP, synBP, realREF, synREF

# ------------------ bandpassfilter ---------------------------


def bandpassfilter(tr, ftype, npad, nscale, pmax, fmult, sigmaIfc,
                   calc_IR=False):
    """
    bandpass filter
    written based on bandpassfilter.m of Karin

    INPUT:
    tr         trace of the seismogram
    plot_resp  (OPTIONAL, default=1) whether or not to plot and save
               filter transfer functions to figures/bpfilter_transfct.jpg

    OUT:
    bpfparms   struct containing filter parameters

    smgrBP    (npad,ndata,nf) bpass-filtered smgrs
              bandpassfilter.m takes care of zeropadding the input waveforms
              DAT and MFI so that there is (almost) no wrap-around energy in
              the outputs smgrBP and mfiBP. npad is hence longer than the
              length of the input tiem series. (TODO: make neater)
    mfiBP     (npad,ndata,nf) bpass-filtered matched filters, does
              include ttime but not amp corrections currently, or empty
              matrix [] if MFI==[]
    TAX       time axis in seconds of filtered time series, starting at t=0

    IR        impulse reponses of the filter, length npad
    """
    DAT = np.reshape(tr.data, (-1, 1))
    dt = 1./tr.stats.sampling_rate

    nptsi = np.shape(DAT)[0]
    ndata = np.shape(DAT)[1]

    # INIT:
    # Length nptso of the filtered time series is max(npad,nptsi),
    # i.e. output has the length of input or of the filter impulse response,
    # depending on which is longer. This is important in order to avoid
    # wrap-around of the filter reponse in the output.

    # In order to avoid wrap-around of the filter, the input arrays DAT and
    # MFI need to get zeropadded: append npad zeros at the end of each, so that
    # the length of each time series passed to the filters will be nptsi+npad,
    # where npad is lthe (sufficiently long) duration of the impulse
    # response(s) determined for the filters.
    # [Length of the output could then be cut arbitrarily short again.]

    # Prepare broadband time series for filtering: append zeros.
    nlong = nptsi + npad
    smgr0 = np.zeros([nlong, ndata])
    smgr0[0:nptsi, 0:ndata] = DAT[:, :]

    # Input time series for impulse response, they are "short" (length npad)
    pulse = np.zeros([npad, 1])
    pulse[0] = 1

    # tshift: Delay each impulse response by 1.1 x center period
    # since by default, acausal filter is max at t=0 and wraps around
    pmid = pmax/(np.power(fmult, np.arange(0, nscale)))
    tshift = 1.1*pmid   # possibly 1.0

    if ftype == 'log-gabor':
        if calc_IR:
            # Compute impulse responses and their durations
            [IR, tmpg1, lngabor, f] = gaborfilter(pulse[:], dt,
                                                  pmax, nscale,
                                                  fmult, sigmaIfc,
                                                  npad, tshift)
            smgrBP = False
            TAX = False
            # Compute the characteristics of the impulse response(s)
            IR = np.squeeze(IR)
        else:
            # Bandpassed time series
            # nf: number of bands = nscale passbands
            # Change 2011/08/10 compared to bpfilter.m:
            # the broadband result smgr0
            # is no longer returned, only the nscale filtered channels.
            nf = nscale
            smgrBP = np.zeros([nlong, ndata, nf])
            # tshift=False to avoid time shift (compatible with kernel code)
            [smgrBP[:, :, 0:nf], fc, lngabor, f] = gaborfilter(smgr0[:, :], dt,
                                                               pmax, nscale,
                                                               fmult, sigmaIfc,
                                                               nlong,
                                                               tshift=False)
            # Make same time axis to go with the output time series,
            # starting at t=0
            TAX = dt*np.arange(0, nlong)
            TAX.transpose()
            IR = False
        return smgrBP, TAX, IR, tshift, lngabor, f
    else:
        cprint('filtering.py', '[ERROR] [FILTERTYPE]', bc.dred,
               '%s has not implemented! EXIT!' % ftype)
        sys.exit()

# ------------------ gaborfilter ---------------------------


def gaborfilter(tsin, dt, pmax, nscale, fmult, sigmaIfc, npad, tshift=False):
    """
    Gabor filter
    written based on gaborfilter.m of Karin

    log-Gabor bandpass filtering (multiplication in the frequency domain by a
    log-Gaussian window)
    G(f,k) = exp( -(ln(f/f_k))^2 / (2*ln(sigmaIfc)^2)  );

    IN:
    tsin    time series to filter (columnwise)
    dt      sampling interval
    pmax    center period in sec of lowest passband filter
    nscale  number of bandpass filters (scales)
    fmult   multiplicative scale factor of filter center freqs
    sigmaIfc  is the (fixed) ratio of sigma divided by fc,
             where sigma is the standard
             deviation of the Gaussian that describes the
             log-Gabor filter in the *time* domain,
             and fc is the filter's center frequency.
             Hence larger sigmaIfc means narrower bandwidth
    npad    desired length of zero-padded time series
             Choose npad>=size(tsin,1)
    tshift  (OPTIONAL) vector of length nscale. Adds a time shift
             to the impulse response of each band (t>0 <--> delay)
             Useful to avoid wrap around on tsin that start at
             or near zero. Choose tshift(k) on the order of
             1/fc(k), where fc is the center frequency of a band.
             Default: tshift(1:nscale) = 0

    OUT:
    tsout   (nn x ndat x nscale) filtered tsin, in all nscale bands
    fc      (nscale x 1)  filter center frequencies
    lngabor (nn x nscale) transfer functions of the log-Gabor filters
    f       (nn x 1) frequency axis for lngabor


    NOTE: The setting of sigmaIfc is empirical given fmult.
    Even coverage of the entire frequency spectrum can be checked
    by plotting the sum of all filter transfer functions.
    E.g., for fmult = sqrt(2), sigmaIfc = .80 is the narrowest-band
    setting that still provides flat coverage.
    sigmaIfc = .5 means the one-sigma interval equals one octave.
    """
    # check whether tshift has been defined or should be continue with
    # its default value (False)
    if np.isscalar(tshift):
        tshift = np.zeros([nscale, 1])

    # time samples
    ns = np.size(tsin, 0)
    # time series
    ndat = np.size(tsin, 1)

    # zero-pad data
    if not npad:
        npad = ns

    # internal (padded length) of each time series
    nn = max(ns, npad)
    tsin = np.append(tsin, np.zeros([npad-ns, ndat]), 0)

    # frequency axis
    # freq resolution
    df = 1./(nn*dt)
    f = df*np.arange(-nn/2., nn/2.)
    f.transpose()

    # *** Temporarily get rid of the zero value in the middle
    # so taking the log does not cause trouble.
    if np.mod(nn, 2) == 0:
        f[int(round(nn/2.+1))-1] = 1

    # filter center frequencies
    fc = (1./pmax)*np.power(fmult, np.arange(0, nscale))
    # in our simple filter gallery:
    # 0.03333333  0.04714     0.06666539  0.09427819  0.13332822  0.18855277
    # 0.26665132  0.3770983

    # make filters
    lngabor = np.zeros([nn, nscale])

    # set
    for k in range(1, nscale+1):
        # f0 = fc[k-1]
        lnf = np.log(abs(f)/fc[k-1])
        lngabor[:, k-1] = np.exp(-np.power(lnf, 2) /
                                 (2.*np.power(np.log(sigmaIfc), 2)))

    # *** Undo the 0-->1 hack
    if np.mod(nn, 2) == 0:
        # bandpass filter has no DC component
        lngabor[int(round(nn/2.)), :] = 0
        f[int(round(nn/2.))] = 0

    tsout = np.zeros([nn, ndat, nscale])
    tsfft = np.zeros([nn, ndat], dtype=complex)
    gb_filter = np.zeros([nn, 1], dtype=complex)
    phas = np.zeros([nn, 1], dtype=complex)   # additional phase shift

    for k in range(0, nscale):
        phas[0:nn, 0] = np.exp(complex(0, -2*np.pi)*f[:]*tshift[k])
        gb_filter[0:nn] = np.fft.fftshift(np.transpose(
            np.column_stack(lngabor[:, k]))*phas[:])
        for idat in range(0, ndat):
            tsfft[:, idat] = np.fft.fft(tsin[:, idat].T)
            tsout[:, idat, k] = np.real(np.fft.ifft(np.transpose(
                np.column_stack(tsfft[:, idat]))*gb_filter[:], axis=0))[:, 0]
            # imag is zero (or rather, less than eps)

    return tsout, fc, lngabor, f

# cutoff frequency at which the filter's power decays to half (-3 dB)
# P/P0 = 0.5 --> 10*log10(P/P0) = -3 db
# A/A0 = sqrt(0.5) --> 20*log10(A/A0) = 10*log10(0.5)= -3 db
#
# if amplitude halves --> -6 db

# ------------------ cut_sample_waveform ---------------------------


def cut_sample_waveform(waveform, sample_to_slice, wave_length):
    """
    cut the trace based on the first sample to start (sample_to_slice) and
    the required length (wave_length).
    This function works with only samples and not time.
    :param waveform:
    :param sample_to_slice:
    :param wave_length:
    :return:
    """
    return waveform[(sample_to_slice-1):(sample_to_slice-1+wave_length)]

# ------------------ raydata_raymatrix_bpf_input ---------------------------


def raydata_raymatrix_bpf_input(IR, lngabor, f, inp_ffproc):
    """
    Generating an input file for raydata_raymatrix ---> resample the gabor filters

    The input file 'bpf_omega_orig.txt' contains the original gabor filters.
    odd columns are time and even columns are gain it goes from 30 sec
    to 2.7 sec dominant period for 8 frequency bands...
    :param IR:
    :param lngabor:
    :param f:
    :param inp_ffproc:
    :return:
    """
    bpf_omega = np.loadtxt(os.path.join('src', 'data', 'bpf_omega_orig.txt'),
                           dtype=float)
    bpf_fio = open(os.path.join(inp_ffproc.output_dir, 'bpf.omega_m'), 'w')
    bpf_fio.writelines('%s # rad. freq/(2*pi*Hz)  abs(gain) (sigma/fc=%s) \n'
                       % (inp_ffproc.lgb_filt_nscale, inp_ffproc.lgb_filt_sigmaIfc))

    cmap = plt.get_cmap('jet_r')
    colors = cmap(np.linspace(0, 1.0, np.shape(IR)[1]))

    plt.ioff()
    plt.figure(figsize=(10, 10))

    counter = 0
    for i in range(np.shape(IR)[1]):
        temp_list = ((x, y) for x, y in zip(f, lngabor[:, i]) if y >= 0.001 and x >= 0)
        temp_x, temp_y = zip(*temp_list)
        x_log = np.logspace(np.log10(min(temp_x)),
                            np.log10(max(temp_x)), 40,
                            endpoint=True, base=10.0)
        f_inter = np.interp(x_log, temp_x, temp_y)

        plt.plot((bpf_omega[:, i+counter] / (2. * np.pi)),
                 bpf_omega[:, i+1+counter],
                 color=colors[i], alpha=0.3,
                 lw=2, label='%s' % inp_ffproc.lgb_filt_center_period[i])
        plt.scatter(x_log, f_inter, alpha=0.5, lw=0, s=7, color=colors[i])
        plt.plot(x_log, f_inter, alpha=0.8, lw=0.5, color=colors[i])

        counter += 1
        line = '%s  # scale=%s, 1/fc= %s s\n' \
               % (len(x_log), counter, inp_ffproc.lgb_filt_center_period[i])
        bpf_fio.writelines(line)
        save_data = np.array([x_log, f_inter]).T
        np.savetxt(bpf_fio, save_data, fmt=['%.4f', '%.4f'])
    bpf_fio.close()

    plt.title('Resampled Gabor filter', weight='bold')
    plt.legend(title='Center Period', fontsize=9)

    plt.xlim(0.8E-2, 10)
    plt.xlabel('Frequency [1/sec]', weight='bold', fontsize=15)
    plt.xscale('log')
    plt.savefig(os.path.join(inp_ffproc.output_dir,
                             'gabor_filts_interp.png'),
                format='png', bbox_inches='tight')

    plt.clf()
    plt.ioff()
    plt.close()
