#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  run_assemblematrix.py
#   Purpose:   Compile and run assemblematrix.f (step03 of inversion)
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")

import datetime
import os
import sys

import utility_codes as outread


# ------------------- INPUT -----------------------------
# vertices file name with grid nodes x, y, z
vertex_file = 'vertices.SEA-SEIS'
# facets file name
facet_file = 'facets.SEA-SEIS'
# 1) invert for dlnVp (y/n), estimated parameter sigma (0.01=1% vp anomalies)
dlnvp = [1, 0.02]
# 2) dlnVs
dlnvs = [0, 0.02]
# 3) dlnQs
dlnQs = [0, 0.3]
# 4) hypoctr corr (in km)
hypocen_corr = [1, 20.0]
# 5) station corr t-P (sec)
sta_corr_tp = [0, 1.0]
# 6) station corr dlnA-P (dimless)
sta_corr_ap = [0, 0.15]
# 7) origin time correction (sec)
# 1: one per event
# 2: one per cluster
orig_corr = [2, 2.5]
# 8) event correction dlnA-P
ev_corr_ap = [0, 0.5]
# 9) station corr t-S
sta_corr_ts = [0, 1.0]
# 10) station corr dlnA-S
sta_corr_as = [0, 1.0]
# 11) EXTENDED hypoctr corr (in km) - better if this switch mirror the one from 4)
hypocen_corr_depth = [1, 5.0]
# corrections: ellipticity, crust, station elevation, Q-dispersion
# corr_switches = [1, 1, 1, 0]
corr_switches = [1, 1, 1, 0]
# "ident" for output files (mat.ident, out.ident, etc) this is also the ident that you get from 
# the ray-data-ray-matrix codes
out_ident = 'ssp'

# master directory
add_master_mat = '/mnt/seismodata/MT/SEA-SEIS_TOMO/rdrm/SEA-SEIS_skewed_15.10.2021b'
add_list_mat = False

#assembled_dir where to copy the files
# assembled_dir = os.path.join(add_master_mat, 'assemble_dir')
assembled_dir = '/mnt/seismodata/MT/SEA-SEIS_TOMO/assmat/SEA-SEIS_skewed_15.10.2021b'

run_assemblematrix = True
do_analysis_plots = True
# ------------------- END INPUT -----------------------------

print "================ ASSEMBLE MATRIX ===================="
t_start = datetime.datetime.now()

if run_assemblematrix:
    if not os.path.isdir(assembled_dir):
       os.mkdir(assembled_dir)

    if add_master_mat:
        outread.make_assemble_dir(add_master_mat,
                                  assembled_dir,
                                  vertex_file,
                                  facet_file,
                                  dlnvp, dlnvs, dlnQs,
                                  hypocen_corr,
                                  sta_corr_tp, sta_corr_ap,
                                  orig_corr, ev_corr_ap,
                                  sta_corr_ts, sta_corr_as,
                                  hypocen_corr_depth,
                                  corr_switches,
                                  out_ident)
    else:
        sys.exit('This mode has not implemented yet!')

    t_cr = datetime.datetime.now()

    outread.compile_assemblematrix(target_add=assembled_dir)
    outread.run_assemblematrix(target_add=assembled_dir)

    t_run = datetime.datetime.now()

    print '\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print "Time to create/copy the required files: %s" % (t_cr - t_start)
    print "Time to run: %s" % (t_run - t_cr)
    print '\nREGULAR END --- %s sec' % (datetime.datetime.now() - t_start)
    print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

if do_analysis_plots:
    outread.plot_assemble_dtheor(assembled_dir, out_ident)
    outread.plot_assemble_stations(assembled_dir, out_ident)
    # outread.plot_assemble_columndensity(assembled_dir, out_ident,
    #                                     vertex_file, facet_file)


raw_input('\nPlease press enter to finish the program!')
