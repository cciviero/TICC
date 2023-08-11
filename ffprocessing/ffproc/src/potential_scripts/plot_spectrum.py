#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  plot_spectrum.py
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GNU Lesser General Public License, Version 3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
import numpy as np
try:
    from obspy.signal.util import nextpow2
except Exception, e:
    from obspy.signal.util import next_pow_2 as nextpow2

# ##################### spectrum_calc ##################################


def spectrum_calc(tr):
    """
    Simple code to calculate the spectrum of a trace
    :param tr: obspy Trace
    :return:
    freqs, spec_tr
    which are frequencies and spectrum of the trace

    To plot the spectrum:

    import matplotlib.pyplot as plt
    plt.loglog(freqs, spec_tr)
    """

    nfft = nextpow2(tr.stats.npts)
    freqs = np.fft.rfftfreq(nfft) / tr.stats.delta

    spec_tr = np.abs(np.fft.rfft(tr.data, n=nfft))

    return freqs, spec_tr
