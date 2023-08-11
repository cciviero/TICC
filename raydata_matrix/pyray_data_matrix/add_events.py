#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  add_events.py
#   Purpose:   Generate event list to be read by raydata_raymatrix.py
# -------------------------------------------------------------------

import glob
import os

# -------------------------------------------------------------------

req_add = '/mnt/seismodata/MT/SEA-SEIS_TOMO/ffproc/SEA-SEIS_skewed_15.10.2021b'
fmt_dir = '*_*.a'
# req_add = '/mnt/seismodata/MT/GLOBAL-TOMO/Pdiff_95_180_OX'
# req_add = '/mnt/seismodata/MT/ANTARCTICA-TOMO/pyffproc_ant_18022021'
# fmt_dir = '*_*.*'
# fmt_dir = '*.*.*.*'
starting_id = 6000

# -------------------------------------------------------------------

print('ATTENTION old event list will be overwritten!')
sel_evs = open('./selected_events_indexed_SEA-SEIS_skewed_15.10.2021b.txt', 'w+')
glob_dirs = glob.glob(os.path.join(req_add, fmt_dir))

# -------------------------------------------------------------------

# dummy event
dummy_event = '/home/mariat/Codes/TICC/raydata_matrix/pyray_data_matrix/utils/6666.6666.666.a'
sel_evs.writelines('%s,%s,??.??,??.??,?,?.?,?,?,?,?,?,?,\n' % (starting_id + 0, os.path.basename(dummy_event)))

# -------------------------------------------------------------------

# events you selected
counter = 1
for g_dir in glob_dirs:
    print ("%s %s" % (counter+starting_id, os.path.basename(g_dir)))
    sel_evs.writelines('%s,%s,??.??,??.??,?,?.?,?,?,?,?,?,?,\n' % (counter+starting_id, os.path.basename(g_dir)))
    counter += 1
sel_evs.close()

print('FINISHED!')
