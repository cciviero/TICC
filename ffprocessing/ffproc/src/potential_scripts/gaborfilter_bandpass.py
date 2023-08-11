"""
gaborfilter and bandpassfilter functions.
They have been tested...
"""

import matplotlib.pyplot as plt
import numpy as np


# //////////////////// gaborfilter ///////////////////////////////////////////


def gaborfilter(tsin, dt, pmax, nscale, fmult, sigmaIfc, npad, tshift=False):
    """
    Gabor filter
    written based on gaborfilter.m of Karin

    Copied from gaborfilter.m
    ==========================
    function [tsout,fc,lngabor,f] = gaborfilter(tsin,dt,pmax,nscale,fmult,...
    sigmaIfc,npad,tshift)
    function [tsout,fc,lngabor,fc] = gaborfilter(tsin,dt,pmax,nscale,fmult,...
    sigmaIfc,npad,tshift)

    log-Gabor bandpass filtering (multiplication in the frequency domain by a
    log-Gaussian window)
    G(f,k) = exp( -(ln(f/f_k))^2 / (2*ln(sigmaIfc)^2)  );

    IN:
    ts      time series to filter (columnwise)
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
             Choose npad>=size(ts,1)
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

    if np.isscalar(tshift):
        tshift = np.zeros([nscale, 1])

    ns = np.size(tsin, 0)  # of time samples
    ndat = np.size(tsin, 1)  # of time series

    # zero-pad data
    if not npad:
        npad = ns

    nn = max(ns, npad)  # internal (padded length) of each time series
    tsin = np.append(tsin, np.zeros([npad-ns, ndat]), 0)

    # frequency axis
    df = 1./(nn*dt)  # freq resolution
    f = df*np.arange(-nn/2, nn/2)
    f.transpose()

    # *** Temporarily get rid of the zero value in the middle
    # so taking the log does not cause trouble.
    if np.mod(nn, 2) == 0:
        f[int(round(nn/2+1))-1] = 1

    # filter center frequencies
    fc = (1./pmax)*np.power(fmult, np.arange(0, nscale))
    # in our simple filter gallery:
    # 0.03333333  0.04714     0.06666539  0.09427819  0.13332822  0.18855277
    # 0.26665132  0.3770983

    # make filters
    lngabor = np.zeros([nn, nscale])

    # set
    for k in range(1, nscale+1):
        f0 = fc[k-1]
        lnf = np.log(abs(f)/fc[k-1])
        lngabor[:, k-1] = np.exp(-np.power(lnf, 2) /
                                 (2.*np.power(np.log(sigmaIfc), 2)))

    # *** Undo the 0-->1 hack
    if np.mod(nn, 2) == 0:
        # bandpass filter has no DC component
        lngabor[int(round(nn/2.)), :] = 0
        f[int(round(nn/2.))] = 0

    tsout = np.zeros([nn, ndat, nscale])
    tsfft = np.zeros([nn, ndat], dtype=np.complex)
    filter = np.zeros([nn, 1], dtype=np.complex)
    phas = np.zeros([nn, 1], dtype=np.complex)   # additional phase shift

    for k in range(0, nscale):
        phas[0:nn, 0] = np.exp(complex(0, -2*np.pi)*f[:]*tshift[k])
        filter[0:nn] = np.fft.fftshift(np.transpose(
            np.column_stack(lngabor[:, k]))*phas[:])
        for idat in range(0, ndat):
            tsfft[:, idat] = np.fft.fft(tsin[:, idat].T)
            tsout[:, idat, k] = np.real(np.fft.ifft(np.transpose(
                np.column_stack(tsfft[:, idat]))*filter[:], axis=0))[:, 0]
            # imag is zero (or rather, less than eps)

    return tsout, fc, lngabor, f

# cutoff frequency at which the filter's power decays to half (-3 dB)
# P/P0 = 0.5 --> 10*log10(P/P0) = -3 db
# A/A0 = sqrt(0.5) --> 20*log10(A/A0) = 10*log10(0.5)= -3 db
#
# if amplitude halves --> -6 db

# //////////////////// bandpassfilter ////////////////////////////////////////


def bandpassfilter(DAT, dt):
    """
    bandpass filter
    written based on bandpassfilter.m of Karin

    function [bpfparms,smgrBP,mfiBP,TAX,IR] = bandpassfilter(DAT,MFI,dt,...
    filterdir,plot_resp);
    function [params,smgrBP,mfiBP,TAX,IR] = bandpassfilter(DAT,MFI,dt,...
    filterdir,plot_resp);

    Compute finite-frequency amplitudes and time delays
    using log-gabor filters (or other kinds of filters, like
    Butterworth, but that's not completely implemented at the moment).

    INPUT:
    DAT        (nptsi,ndata) column matrix of broadband seismograms
    dt         sampling interval
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

    nptsi = np.shape(DAT)[0]
    ndata = np.shape(DAT)[1]

    # HARD CODED based on a specific filter that I have used for all my
    # measurements, this should be more generic, similar to Karin's code to
    # read a text file in order to set the following parameters
    ftype = 'log-Gabor'
    npad = 1024
    nscale = 8
    pmax = 30
    fmult = 1.4142
    sigmaIfc = 0.5
    nout = 40

    #===================================================
    #   START FILTERING SECTION
    #===================================================

    # INIT:

    # Length nptso of the filtered time series is max(npad,nptsi),
    # i.e. output has the length of input or of the filter impulse response,
    # depending on which is longer. This is important in order to avoid
    # wrap-around of the filter reponse in the output.

    # In order to avoid wrap-around of the filter, the input arrays DAT and
    ## MFI
    # need to get zeropadded: append npad zeros at the end of each, so that
    ## the
    # length of each time series passed to the filters will be nptsi+npad,
    # where npad is lthe (sufficiently long) duration of the impulse
    ## response(s)
    # determined for the filters.
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

    # Compute impulse responses and their durations
    [IR, tmpg1, tmpg2, tmpg3] = gaborfilter(pulse[:], dt, pmax, nscale, fmult,
                                      sigmaIfc, npad, tshift)

    # Bandpassed time series
    # nf: number of bands = nscale passbands
    # Change 2011/08/10 compared to bpfilter.m: the broadband result smgr0
    ## is no longer returned, only the nscale filtered channels.
    nf = nscale
    smgrBP = np.zeros([nlong, ndata, nf])
    [smgrBP[:, :, 0:nf], fc, lngabor, f] = \
        gaborfilter(smgr0[:, :], dt, pmax, nscale, fmult, sigmaIfc, nlong,
                    tshift)

    # Make sime time axis to go with the output time series, starting at t=0
    TAX = dt*np.arange(0, nlong)
    TAX.transpose()

    # Compute the characteristics of the impulse response(s)
    IR  = np.squeeze(IR)

    return smgrBP, TAX, IR



# ====================== TRASH ===============================================
# ====================== BANDPASS ============================================
    ### duration = time until fraction "thres" of energy has been spent
    ##thres = 0.90   # <~ 1 (e.g. 0.90 )
    ##len_n = np.zeros([nscale, 1])
    ##tmp = np.cumsum(np.power(IR, 2), axis=0)
    ##for k in range(0, nscale):
    ##    ii = min(np.where(tmp[:, k]/tmp[-1, k] >= thres))
    ##    if not np.isscalar(ii):
    ##        import ipdb; ipdb.set_trace()
    ##        len_n[k] = ii
    ##    else:
    ##        len_n[k] = npad

    # In case bpfparms was filled manually:
    # Save center frequencies in Hz
    # Save response duration in sec
    #new_parms_struct = 0;
    #if(not(isfield(bpfparms,'response_duration')))
    #  bpfparms.response_duration = dt*len_n(:);
    #  new_parms_struct = 1;
    #end
    #if(not(isfield(bpfparms,'fc')))
    #  bpfparms.fc = fc;
    #  new_parms_struct = 1;
    #end
    #if(not(isfield(bpfparms,'flag')))
    #  bpfparms.flag = 0; % default exit flag
    #  new_parms_struct = 1;
    #end

    #% struct bpfparms is complete -- write-protect it
    #% Make sure this file is READ-ONLY (so that earlier versions
    #% cannot get overwritten. Write them to different filterdir instead!)
    #if(new_parms_struct)
    #  [succs msg] = fileattrib(fnam,'-w','a');
    #end