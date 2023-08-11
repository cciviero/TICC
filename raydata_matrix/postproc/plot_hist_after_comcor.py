#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  plot_hist_trti.py
#   Purpose:   plot bandwise histogram of travel times with common
#              corrections already applied
#   Author:    Mary
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import sys
import os
import glob

# -----------------------------------------------------------------------
# ------------------------------- INPUT ---------------------------------
# -----------------------------------------------------------------------

# ====>>>> input file: 'synthetic_calculated. ...'
# input_path = sys.argv[1]
#path = '/Users/maria/PhD/Codes/Measurements/FFINVERSION/INVERSION/Programs/' \
#             'pyfiles/pyray_data_matrix/RESULTS_BACKUPs/RESULTS_CRUST1'

# path = '/Users/maria/PhD/Codes/Measurements/FFINVERSION/INVERSION/Programs/pyfiles/pyray_data_matrix/RESULTS_22.08.2016'
# path = '/Users/maria/PhD/Codes/Inversion/2016_inversion/homo_400_rhum_rum_P_v0103/raymatrix_results'
path =  '/data/maria/global_tomography/tomo_P_Pdiff_NA_UCL/rdrm_NA_P_Pdiff_UCL_2016-2022'
bands = range(8)
cores = range(20)

# -----------------------------------------------------------------------
# ----------------------------- Functions -------------------------------
# -----------------------------------------------------------------------


def read_info_file(input_path):

    ## column information of those info collector files
    # 0 kd
    # 1 stationcode
    # 2 netw
    # 3 ptlat
    # 4 ptlon
    # 5 rtarget
    # 6 slat
    # 7 slon

    #*** 8 ecorr - ell
    #*** 9 tau - crustal correction
    #*** 10 telev - eleve

    # 11 tstar - attenuation
    # 12 del - epicentral distance?
    #*** 13 trtime_final - beginning of match filter - theoretical arrival time on geocentric
    #*** 14 tobs - tobs
    # 15 sdep - source depth

    dt = np.dtype([('kd', int), ('sc', 'S10'), ('net', 'S10'), ('ptlat', float), ('ptlon', float), ('rtarget', float),
                   ('slat', float), ('slon', float), ('ecorr', float), ('tau', float), ('telev', float),
                   ('tstar', float), ('del', float), ('trtime', float), ('tobs', float), ('sdep', float)])

    info_file = glob.glob(input_path[0])

    print ('Working on %s' % info_file[0])

    crust_info = np.loadtxt(info_file[0], dtype=dt, delimiter=',')
    up = crust_info[1::2]
    down = crust_info[::2]

    return crust_info, up, down

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# -----------------------------------------------------------------------
# ------------------------------ Main part ------------------------------
# -----------------------------------------------------------------------
# plotting information
fig = plt.figure(frameon=True)
fig.set_size_inches(30, 20)
# for the specific band different colors and linewidths
cmap = plt.get_cmap('jet_r')
colors = cmap(np.linspace(0, 1.0, len(bands)))
line_width = [x + 1 for x in bands[::-1]]
band_freq = [30.0, 21.21, 15.0, 10.6, 7.50, 5.30, 3.70, 2.65]
bins_bin = 20
plt.axvline(x=0, color='r', lw=3, linestyle='--')
plt.axvline(x=0.5, color='r', lw=2, linestyle='--')
plt.axvline(x=-0.5, color='r', lw=2, linestyle='--')
plt.axvline(x=1, color='r', lw=1, linestyle='--')
plt.axvline(x=-1, color='r', lw=1, linestyle='--')
plt.axvline(x=-2, color='orange', lw=1, linestyle='--')
plt.axvline(x=2, color='orange', lw=1, linestyle='--')

counter = 0


for band in bands:
    trti_corr = np.array([])

    for core in cores:

        band_dir = '/*band0%s_%s_*dir' % (band+1, core+1)
        path2info = glob.glob(os.path.join(path + band_dir +'/info_collector.*'))
        # print(core, cores, band)
        # import ipdb; ipdb.set_trace()
        crust_info, up, down = read_info_file(path2info)

        for i in range(int(len(crust_info)/2)):
            # *** 8 ecorr - ell
            # *** 9 tau - crustal correction
            # *** 10 telev - eleve
            # *** 13 trtime_final
            # *** 14 tobs - tobs
            syn_corr = up[i]['tobs'] + down[i]['ecorr'] +\
                       up[i]['tau'] + down[i]['tau'] +\
                       up[i]['telev'] + down[i]['telev']

            final_corr = up[i]['trtime'] - syn_corr
            trti_corr = np.append(trti_corr, -final_corr)
            counter +=1

    bins_range = [np.min(trti_corr), np.max(trti_corr)]
    y, binEdges = np.histogram(trti_corr, bins=bins_bin, range=bins_range)
    bin_centers = 0.5 * (binEdges[1:] + binEdges[:-1])

    # import ipdb; ipdb.set_trace()
    plt.plot(bin_centers, y, color=colors[band], lw=line_width[band],
             label='B0%s - %2.2fs: %s data points' % (band + 1, band_freq[band], len(trti_corr)))

plt.xlabel('travel time', size=30, weight='bold')
plt.ylabel('number of measurements', size=30, weight='bold')
plt.title('Common corrected travel time measurements', size=40, weight='bold')
# \nCriteria: cc>=%s -- ts <=%s -- clip=%s'
#           % (cc_factor, ts_factor, clip_factor), size=40, weight='bold')
plt.legend(fontsize=25, loc=1)

plt.grid(True)
# to but the text relative to the x/y axis

plt.annotate('Number measurements (all bands): %s\n' % counter, xy=(0.02, 0.90), size=23, xycoords='axes fraction',
              bbox=dict(boxstyle="round", facecolor='red', alpha=0.8))

plt.xticks(weight='bold', size=18)
plt.yticks(weight='bold', size=18)

# plt.show()

output_dir = path

plt.savefig(os.path.join(output_dir + '/common_corrected_measurements.png'))


# import ipdb; ipdb.set_trace()
















