"""
Simple script to combine all filt_array's (*.npy)
in order to plot figures for all frequencies or phases or...
"""

import glob
import numpy as np
import os
import sys

sys.path.append('..')

from output_reader import check_selection
from output_reader import write_staev_locfile

# ================= INPUT
npy_dir = '../OUTPUTS/P_Pdiff/npy_P_Pdiff'
# Local variables:
phase = 'band1_8'
fnam = 'P_Pdiff_'
# Criteria:
min_xcorr = 0.8
max_xcorr = 1.01
min_epi = 32.0
max_epi = 160.01
# =================

npy_add = glob.glob(os.path.join(npy_dir, '*.npy'))
npy_add.sort()

print 'Appending: %s' % npy_add[0]
filt_array = np.load(npy_add[0])
for i in range(1, len(npy_add)):
    print 'Appending: %s' % npy_add[i]
    filt_array = np.append(filt_array, np.load(npy_add[i]), 0)

print 'Check the selected event-station pairs!'
check_selection(filt_array, min_xcorr, max_xcorr, min_epi, max_epi)
print 'Write event-station location files!'
write_staev_locfile(filt_array, phase, fnam, phase)

