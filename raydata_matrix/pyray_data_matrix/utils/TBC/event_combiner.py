#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  event_combiner.py
#   Purpose:   Combine the events (for different frequency bands) to be used with other tools
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import glob
import os
import shutil
import sys


source_dir = './RESULTS'
target_dir = './event_comparison'

events = glob.glob(os.path.join(source_dir, '*.*.*.*'))
for i in range(len(events)):
    source_cp_all = glob.glob(os.path.join(events[i], 'outfiles', 'ffproc.ampstt.*'))
    if len(source_cp_all) != 1:
        sys.exit('There is something wrong! Believe me! :)')
    source_cp = source_cp_all[0]
    target_dir_cp = os.path.join(target_dir, os.path.basename(events[i]), 'outfiles')
    shutil.copy(source_cp, target_dir_cp)
