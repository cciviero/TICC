#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  rotate_RR_statioins.py
#   Purpose:   edit metadata file of RHUM RUM stations to fill with
#              azimuth value derived from JRS measurements
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
import glob
import matplotlib.pyplot as plt
import numpy as np
import os

from obspy import read_inventory

# =====>>>>> read text file with information

info_rot = np.loadtxt('RR_rotations2.txt', dtype=object, delimiter='\t', usecols=(0,6))

# =====>>>>> read stationxml files

path_2_xml = '/net/siglochnas1/volume1/data1/mariat/RHUM-RUM/data/rhum_rum_2011_2015'
read_all_events = glob.glob(os.path.join(path_2_xml, '*.*.*.*'))

for event in read_all_events:
    print '================= %s =================' % os.path.basename(event)
    for station in info_rot:
        stat = station[0]
        azimuth = station[1]
        print '................. %s .................' % stat
        try:
            path_xml_bh1 = glob.glob(os.path.join(event, 'resp', 'STXML.YV.' + stat + '.00.BH1'))[0]
            if os.path.exists(path_xml_bh1):
                inv_bh1 = read_inventory(path_xml_bh1)
                # =====>>>>> rename stationxml files to '*orig*' files
                inv_bh1.write(os.path.join(event, 'resp', 'orig_STXML.YV.' + stat + '.00.BH1'), format='STATIONXML')
                print 'azimuth_orig BH1 = ', inv_bh1[0][0][0].azimuth
                inv_bh1[0][0][0].azimuth = azimuth
                print 'azimuth_new BH1 = ', inv_bh1[0][0][0].azimuth
                # =====>>>>> write new stationxml files with aizmuth values
                inv_bh1.write(os.path.join(event, 'resp', 'STXML.YV.' + stat + '.00.BH1'), format='STATIONXML')

        except Exception, exp:
            print exp
            print '%s - %s - %s does not exist...' % (os.path.basename(event), stat, 'BH1')

        print ' - - - - - - -'

        try:
            path_xml_bh2 = glob.glob(os.path.join(event, 'resp', 'STXML.YV.' + stat + '.00.BH2'))[0]
            if os.path.exists(path_xml_bh2):
                inv_bh2 = read_inventory(path_xml_bh2)
                # =====>>>>> rename stationxml files to '*orig*' files
                inv_bh2.write(os.path.join(event, 'resp', 'orig_STXML.YV.' + stat + '.00.BH2'), format='STATIONXML')
                print 'azimuth_orig BH2 = ', inv_bh2[0][0][0].azimuth
                bh2_az = eval(azimuth) + 90
                if bh2_az > 360:
                    bh2_az = bh2_az - 360
                inv_bh2[0][0][0].azimuth = bh2_az
                print 'azimuth_new BH1= ', inv_bh2[0][0][0].azimuth
                # =====>>>>> write new stationxml files with aizmuth values
                inv_bh2.write(os.path.join(event, 'resp', 'STXML.YV.' + stat + '.00.BH2'), format='STATIONXML')

        except Exception, exp:
            print exp
            print '%s - %s - %s does not exist...' % (os.path.basename(event), stat, 'BH2')






