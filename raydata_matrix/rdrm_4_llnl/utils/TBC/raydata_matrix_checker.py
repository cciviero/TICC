#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  raydata_matrix_checker.py
#   Purpose:   Check the raydata and raymatrix outputs
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

"""
Folder structure:
main_dir = RESULTS
                       RESULTS
                         |
                  ---------------
                  |              |
                  P_....         ...
                  |
           -------------------------
           |                      |
  synthetic_calculated.xxxx      ...
"""

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
import sys

# ------------------- INPUT -----------------------------
main_dir = '../RESULTS'
req_phase = 'P,Pdiff'
# ------------------- END INPUT -------------------------

ls_out = glob.glob(os.path.join(main_dir, '*', 'synthetic_calculated.*'))
print '=========='
print '#files: %s' % len(ls_out)
print '=========='

res = np.loadtxt(ls_out[0], delimiter=',', dtype='S')
for i in range(1, len(ls_out)):
    sys.stdout.flush()
    res = np.vstack([res, np.loadtxt(ls_out[i], delimiter=',', dtype='S')])

# Initializing a figure:
plt.ion()

if len(res) < 100:
    alpha_plt = 1.0
elif 100 <= len(res) < 10000:
    alpha_plt = 0.2
elif 10000 <= len(res) < 50000:
    alpha_plt = 0.1
elif 50000 <= len(res) < 200000:
    alpha_plt = 0.01
else:
    alpha_plt = 0.002

plt.figure()
plt.subplot(2, 1, 1)
plt.scatter(res[:, 3].astype(np.float),
            res[:, 8].astype(np.float), c='r', edgecolors='none',
            alpha=alpha_plt)
plt.xlabel('Station latitude', size='large', weight='bold')
plt.ylabel('Ellipticity Correction', size='large', weight='bold')
plt.xticks(size='large', weight='bold')
plt.yticks(size='large', weight='bold')

plt.vlines(50.0, min(res[:, 8].astype(np.float)),
           max(res[:, 8].astype(np.float)),
           'k', linestyles='-')
plt.vlines(26.0, min(res[:, 8].astype(np.float)),
           max(res[:, 8].astype(np.float)),
           'k', linestyles='-')
plt.vlines(60.0, min(res[:, 8].astype(np.float)),
           max(res[:, 8].astype(np.float)),
           'k', linestyles='-')
plt.vlines(34.0, min(res[:, 8].astype(np.float)),
           max(res[:, 8].astype(np.float)),
           'k', linestyles='-')

plt.subplot(2, 1, 2)
plt.scatter(res[:, 4].astype(np.float),
            res[:, 8].astype(np.float), c='r', edgecolors='none',
            alpha=alpha_plt)
plt.xlabel('Station longitude', size='large', weight='bold')
plt.ylabel('Ellipticity Correction', size='large', weight='bold')
plt.xticks(size='large', weight='bold')
plt.yticks(size='large', weight='bold')

plt.vlines(-127.0, min(res[:, 8].astype(np.float)),
           max(res[:, 8].astype(np.float)), 'k', linestyles='-')
plt.vlines(-60.0, min(res[:, 8].astype(np.float)),
           max(res[:, 8].astype(np.float)), 'k', linestyles='-')
plt.vlines(-12.0, min(res[:, 8].astype(np.float)),
           max(res[:, 8].astype(np.float)), 'k', linestyles='dotted')
plt.vlines(41.0, min(res[:, 8].astype(np.float)),
           max(res[:, 8].astype(np.float)), 'k', linestyles='dotted')

taup_time = []
for val in range(np.shape(res)[0]):
    taup_process = subprocess.Popen(['taup_time', '-mod', 'iasp91', '-time',
                                     '-h', str(float(res[val, 17])),
                                     '-ph', req_phase,
                                     '-deg', str(float(res[val, 12]))],
                                    stdout=subprocess.PIPE)
    tt_raw = taup_process.communicate()[0]

    try:
        tt = tt_raw.split('\n')[0].split()[-1]
    except Exception, e:
        print 'Can not calculate the time with taup!'
        tt = -9999

    taup_time.append(tt)


plt.figure()
plt.subplot(2, 1, 1)
plt.scatter(res[:, 12].astype(np.float),
            res[:, 14].astype(np.float) - res[:, 13].astype(np.float),
            c='r', edgecolors='none',
            alpha=alpha_plt)
plt.xlabel('Epicentral distance', size='large', weight='bold')
plt.ylabel('Calculated - MFI', size='large', weight='bold')
plt.xticks(size='large', weight='bold')
plt.yticks(size='large', weight='bold')

plt.subplot(2, 1, 2)
plt.scatter(res[:, 12].astype(np.float),
            res[:, 15].astype(np.float),
            c='r', edgecolors='none',
            alpha=alpha_plt)
plt.xlabel('Epicentral distance', size='large', weight='bold')
plt.ylabel('%(Calculated - MFI)/MFI', size='large', weight='bold')
plt.xticks(size='large', weight='bold')
plt.yticks(size='large', weight='bold')

plt.figure()
plt.scatter(res[:, 12].astype(np.float),
            res[:, 16].astype(np.float) - res[:, 13].astype(np.float),
            c='r', edgecolors='none',
            alpha=alpha_plt)
plt.xlabel('Epicentral distance', size='large', weight='bold')
plt.ylabel('tobs - tmfi', size='large', weight='bold')
plt.xticks(size='large', weight='bold')
plt.yticks(size='large', weight='bold')

plt.figure()
plt.scatter(res[:, 12].astype(np.float),
            res[:, 14].astype(np.float) - np.array(taup_time).astype(np.float),
            c='r', edgecolors='none',
            alpha=alpha_plt)
plt.xlabel('Epicentral distance', size='large', weight='bold')
plt.ylabel('calculated - taup', size='large', weight='bold')
plt.xticks(size='large', weight='bold')
plt.yticks(size='large', weight='bold')


plt.figure()
plt.scatter(res[:, 12].astype(np.float),
            abs(res[:, 14].astype(np.float) -
                np.array(taup_time).astype(np.float)) /
            np.array(taup_time).astype(np.float) * 100,
            c='r', edgecolors='none',
            alpha=alpha_plt)
plt.xlabel('Epicentral distance', size='large', weight='bold')
plt.ylabel('%(calculated - taup)/taup', size='large', weight='bold')
plt.xticks(size='large', weight='bold')
plt.yticks(size='large', weight='bold')

raw_input('Press enter...')
