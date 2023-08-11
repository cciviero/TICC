#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  test_maker.py
#   Purpose:   Create output files similar to RunFFprocessing
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

# Required Python modules will be imported in this part.
try:
    from obspy.geodetics import locations2degrees
except Exception, e:
    from obspy.core.util import locations2degrees
import os
import shutil


# ############### INPUT
# Since this code is meant to be fully compatible with RunFFProcessing outputs,
# it just accepts ONE event and several stations!
evdp = 0.0
evlat = 90.
evlon = 0.

stel = [0.0]
# stlat = range(-10, -60, -10)
stlat = [-30]
stlon = [0]

bands = ['band01']

# Other inputs:
event_name = '6666.6666.666.a'
# ############### END INPUT

# ############################# FUNCTIONS #####################################


def ampinv_source_writer(evlat, evlon, evdp, add):
    """
    Write ampinv.source with the following format:
     year, julianday, hr, min, sec, msec
      2009   273    10    16     9   249
     evlat, evlon, catalog_depth, inverted_depth
     -0.7200000       99.86700       84.00000       82.00000
     scalar_moment, tausrc
      2.2752803E+20   19.44691
     xmom_oest(1:6)
       0.175E+21  -0.260E+20  -0.149E+21   0.760E+20  -0.510E+20  -0.129E+21
     strike, dip, rake
       187.0000       52.00000       53.00000
     event catalog
     NEIC
     inverted moment tensor
       0.163E+21  -0.148E+20  -0.148E+21   0.349E+20  -0.397E+20  -0.157E+21
    """
    line01 = 'year, julianday, hr, min, sec, msec\n'
    line02 = ' 2666   66    1    1     1   1\n'
    line03 = 'evlat, evlon, catalog_depth, inverted_depth\n'
    line04 = '%s    %s    %s    %s\n' % (evlat, evlon, evdp, evdp)
    line05 = 'scalar_moment, tausrc\n'
    line06 = ' 1.0E+20   20.0\n'
    line07 = 'xmom_oest(1:6)\n'
    line08 = '  1.E+20  1.0E+20  1.0E+20   1.0E+20  1.0E+20  1.0E+20\n'
    line09 = 'strike, dip, rake\n'
    line10 = '  45.0000       45.00000       45.00000\n'
    line11 = 'event catalog\n'
    line12 = 'NEIC\n'
    line13 = 'inverted moment tensor\n'
    line14 = '  1.E+20  1.0E+20  1.0E+20   1.0E+20  1.0E+20  1.0E+20'

    fio = open(os.path.join(add, 'ampinv.source'), 'w')
    fio.writelines(line01)
    fio.writelines(line02)
    fio.writelines(line03)
    fio.writelines(line04)
    fio.writelines(line05)
    fio.writelines(line06)
    fio.writelines(line07)
    fio.writelines(line08)
    fio.writelines(line09)
    fio.writelines(line10)
    fio.writelines(line11)
    fio.writelines(line12)
    fio.writelines(line13)
    fio.writelines(line14)
    fio.close()

    shutil.copy(os.path.join(add, 'ampinv.source'), os.path.join(add, 'ffproc.source'))


def ffproc_receivers_writer(stlat, stlon, stel, add):
    """
    # 0274.2009.273.a  idx_ai;  grpaff; BASENAME_UNIQUE; stla; stlo; elevation(meters); burial(meters); sacfilename
      1    1                TA.109C..x00.BHZ     32.8889   -117.1051       150.0         0.0   dis.109C..BHZ
      2    1                TA.113A..x00.BHZ     32.7683   -113.7667       118.0         0.0   dis.113A..BHZ
      3    1                TA.120A..x00.BHZ     32.5466   -108.6330      1528.0         0.0   dis.120A..BHZ
    """
    fio = open(os.path.join(add, 'ffproc.receivers'), 'w')
    line01 = '# 6666.6666.666.a  idx_ai;  grpaff; BASENAME_UNIQUE; stla; stlo; elevation(meters); burial(meters); ' \
             'sacfilename\n'
    fio.writelines(line01)
    for i in range(len(stlat)):
        sta_name = '%3s' % i
        sta_name = sta_name.replace(' ', '0')
        fio.writelines('%s    1         TS.TEST%s..x00.BHZ  %s  %s  %s  0.0  dis.%s..BHZ\n'
                       % ((i+1), sta_name, stlat[i], stlon[i], stel[i], sta_name))
    fio.close()


def ffproc_ampstt_bands_writer(stlat, stlon, evlat, evlon, evdp, band):
    """
    # idx grp    stla      stlo       stazie                stnam            xc_coeff    Tobs         dT      sigma_dT
    not_used   not_used      A          sigma_A      not_used     not_used     tB_smgr  tB_mfi    winlen clip_taumax
    not_used  not_used
    1   1    32.8889  -117.1051   132.6603              TA.109C..x00.BHZ    0.948     966.00       5.40       0.40
    0.00       0.00     1.1734e+08   0.0000e+00   0.0000e+00   0.0000e+00    966.00   960.60    65.10        0
    0.00       0.00
    """
    winlen_dic = {
        'band01': 65.10,
        'band02': 46.20,
        'band03': 32.70,
        'band04': 23.10,
        'band05': 16.40,
        'band06': 11.60,
        'band07': 8.30,
        'band08': 5.90
    }
    fio = open(os.path.join(add, 'ffproc.ampstt.%s' % band), 'w')
    for i in range(len(stlat)):
        sta_name = '%3s' % i
        sta_name = sta_name.replace(' ', '0')
        epicen = locations2degrees(evlat, evlon, stlat[i], stlon[i])
        fio.writelines('%s  1  %s  %s  %s  TS.TEST%s..x00.BHZ  1.0  1.0  1.0  0.2  0.00  0.00  1.0  0.0  0.0  0.0  1'
                       '.0  1.0  %s  0 0.0  0.0\n' % (i+1, stlat[i], stlon[i], epicen, sta_name, winlen_dic[band]))
    fio.close()

# ############################ MAIN PROGRAM ###################################


if not os.path.isdir(os.path.join(os.path.curdir, event_name)):
    os.mkdir(os.path.join(os.path.curdir, event_name))
else:
    print 'WARNING: %s exists!' % os.path.join(os.path.curdir, event_name)
if not os.path.isdir(os.path.join(os.path.curdir, event_name, 'outfiles')):
    os.mkdir(os.path.join(os.path.curdir, event_name, 'outfiles'))
else:
    print 'WARNING: %s exists!' % os.path.join(os.path.curdir, event_name, 'outfiles')
add = os.path.join(os.path.curdir, event_name, 'outfiles')
ampinv_source_writer(evlat=evlat, evlon=evlon, evdp=evdp, add=add)
ffproc_receivers_writer(stlat=stlat, stlon=stlon, stel=stel, add=add)

for band in bands:
    ffproc_ampstt_bands_writer(stlat=stlat, stlon=stlon, evlat=evlat, evlon=evlon, evdp=evdp, band=band)
