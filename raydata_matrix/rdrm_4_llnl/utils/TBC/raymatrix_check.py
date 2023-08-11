#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  raymatrix_check.py
#   Purpose:   Check the raymatrix output (raymatrix.out.*)
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
           ------------------
           |                |
    raymatrix.out.xxxx      ...
"""

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# ------------------- INPUT -----------------------------
main_dir = '../RESULTS'
# ------------------- END INPUT -------------------------

# Initializing a figure:
plt.figure()
plt.ion()
plt.show()

ls_matout = glob.glob(os.path.join(main_dir, '*', 'raymatrix.out.*'))
print '========================='
print '#raymatrix.out. files: %s' % len(ls_matout)
print '========================='

counter = 1
#rayout_array = np.array([])
for i in range(len(ls_matout)):
    sys.stdout.flush()
    fio_rayout = open(ls_matout[i], 'r')
    fi_rayout = fio_rayout.readlines()
    rayout_array_tmp = np.array([])
    for j in range(10, len(fi_rayout)):
        if fi_rayout[j].strip() == '':
            break
        if np.size(rayout_array_tmp) == 0:
            rayout_array_tmp = np.array(fi_rayout[j].split())
        else:
            rayout_array_tmp = np.vstack([rayout_array_tmp,
                                          np.array(fi_rayout[j].split())])

    #if np.size(rayout_array) == 0:
    #    rayout_array = np.array(rayout_array_tmp)
    #else:
    #    rayout_array = np.vstack([rayout_array, rayout_array_tmp])

    if np.size(rayout_array_tmp) == 14:
        rayout_array_tmp = np.reshape(rayout_array_tmp, [1, 14])

    if np.size(rayout_array_tmp[(rayout_array_tmp[:, 12].astype(np.float) > 105.0), 12]) != 0:
        color = 'r'
        print '\n%s/%s' % (i+1, len(ls_matout))
        print '%s' % ls_matout[i]
        print rayout_array_tmp[rayout_array_tmp[:, 12].astype(np.float) > 105.0, 1]
    elif np.size(rayout_array_tmp[(rayout_array_tmp[:, 12].astype(np.float) < 95.0), 12]) != 0:
        color = 'r'
        print '\n%s/%s' % (i+1, len(ls_matout))
        print '%s' % ls_matout[i]
        print rayout_array_tmp[rayout_array_tmp[:, 12].astype(np.float) < 95.0, 1]
    else:
        print '%s/%s' % (i+1, len(ls_matout)),
        color = 'b'

    # Accuray:
    plt.plot(np.linspace(counter, counter+1,
                         len(rayout_array_tmp[:, 12].astype(np.float))),
             rayout_array_tmp[:, 12].astype(np.float) - 100.0, '.', c=color)
    plt.xlabel('Groups', size='xx-large', weight='bold')
    plt.ylabel('Accuracy (Sparse), Deviation from 100%', size='xx-large',
               weight='bold')
    plt.xticks(size='xx-large', weight='bold')
    plt.yticks(size='xx-large', weight='bold')
    counter += 1
    plt.draw()
raw_input('\n\nPress enter to close the program!')
