#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
# Filename:  test_InpFFproc.py
#   Purpose:   testing InpFFproc class
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

# Required Python and Obspy modules will be imported in this part.
from ffprocessing.src.InpFFproc import InpFFproc

# ##################### test_InpFFproc ###############################


def test_inpffproc():
    """
    Testing the generated input class from InpFFproc
    :return:
    """
    inp_ffproc = InpFFproc(
        inp_ffproc_file="./tests/files/input_ffprocessing.ini")
    assert inp_ffproc.ph_phase == 'P'
    assert inp_ffproc.ph_channel == ['BHZ']
    assert inp_ffproc.ph_min_epi == 49.
    assert inp_ffproc.ph_max_epi == 50.
    assert inp_ffproc.ph_bg == 'iasp91'
    assert inp_ffproc.ph_sampling_rate == 10.
    assert inp_ffproc.ph_preset == 20.
    assert inp_ffproc.ph_offset == 20.
    assert inp_ffproc.check_evs is True
    assert inp_ffproc.check_ttime is True
    assert inp_ffproc.check_ev_stas is True
    assert inp_ffproc.check_cut_syn is True
    assert inp_ffproc.check_cut_real is True
    assert inp_ffproc.check_cut_real_syn is True
    assert inp_ffproc.check_bands == '1-8'
    assert inp_ffproc.check_min_bands == 0
    assert inp_ffproc.check_max_bands == 8
    assert inp_ffproc.check_min_cc == 0.8
    assert inp_ffproc.check_measured_glob is True
    assert inp_ffproc.real_mode == 'read'
    assert inp_ffproc.real_path == '/home/hosseini/Work/Scripts/gitHUB/' \
                                   'pyffproc/data/real_event'
    assert inp_ffproc.real_name_format == 'BH'
    assert inp_ffproc.syn_mode == 'read'
    assert inp_ffproc.syn_path == '/home/hosseini/Work/Scripts/gitHUB/' \
                                  'pyffproc/data/syn_event'
    assert inp_ffproc.syn_name_format == 'SAC_realName'
    assert inp_ffproc.evproc_mode == 'read'
    assert inp_ffproc.evproc_path == '/home/hosseini/Work/Scripts/gitHUB/' \
                                     'pyffproc/data/processed_event'
    assert inp_ffproc.evproc_name_format == '*.*.*.*'
    assert inp_ffproc.evproc_min_mag == 5.
    assert inp_ffproc.evproc_max_mag == 8.5
    assert inp_ffproc.evproc_min_depth == 60.
    assert inp_ffproc.evproc_max_depth == 1000.
    assert inp_ffproc.evproc_min_year == 2009
    assert inp_ffproc.evproc_max_year == 2010
    assert inp_ffproc.stf_mode == 'read'
    assert inp_ffproc.mmeant_mode == 'CC-2step'
    assert inp_ffproc.mmeant_clip_time1 == 7.
    assert inp_ffproc.mmeant_clip_time2 == 2.
    assert inp_ffproc.filt_mode == 'log-gabor'
    assert inp_ffproc.lgb_filt_pmax == 30.
    assert inp_ffproc.lgb_filt_nscale == 8
    assert inp_ffproc.lgb_filt_fmult == 1.4142
    assert inp_ffproc.lgb_filt_sigmaIfc == 0.5
    assert inp_ffproc.lgb_filt_npad == 1024
    assert inp_ffproc.lgb_filt_duration == []
    assert inp_ffproc.lgb_filt_energy_frac == 0.95
    assert inp_ffproc.lgb_filt_nlambda == 1.5