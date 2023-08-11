#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  add_events.py
#   Purpose:   Generate event list to be read by raydata_raymatrix.py
# -------------------------------------------------------------------

import glob
import os

# -------------------------------------------------------------------

# req_add = '/data/maria/UPFLOW/TICC/measurement_check/pyffproc_2022-12-01_scardec_iasp'
# req_add = '/data/maria/global_tomography/measurements/Pdiff_95_180_OX'
# req_add = '/data/maria/tomography/measurements/UPFLOW_2021-2022'
req_add = '/Volumes/EUROPE_DATA/eu-data-outfiles2012'
fmt_dir = '*_*.a'
# fmt_dir = '*_*.*'
#fmt_dir = '*.*.*.*'
starting_id = 0

# -------------------------------------------------------------------

print('ATTENTION old event list will be overwritten!')
sel_evs = open('./selected_events_indexed_eu_tomo2012.txt', 'w+')
glob_dirs = glob.glob(os.path.join(req_add, fmt_dir))

# -------------------------------------------------------------------

# dummy event
dummy_event = '/Users/cciviero/Documents/SO_PROJECT/Global_tomography/TICC/TICC-light_chiara/raydata_matrix/pyray_data_matrix/utils/6666.6666.666.a'
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
