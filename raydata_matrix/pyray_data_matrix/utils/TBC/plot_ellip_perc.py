"""
Simple plotting tool to show the difference between
ellipticity correction methods (ray-tracer, BLNK)
"""

import glob
import matplotlib.pyplot as plt
import os
import sys

dirname = os.path.abspath(sys.argv[1])
ellip_files = glob.glob(os.path.join(dirname, '*', 'ellipticity_comparison.*'))

diff = []
for filename in ellip_files:
    fio = open(filename, 'r')
    fi = fio.readlines()

    for i in range(len(fi)):
        diff.append(float(fi[i].split(',')[2]))

plt.hlines(0.1, 0, len(diff), lw=1, linestyle='dashed')
plt.plot(diff, 'k', lw=3)
plt.xticks(size='large', weight='bold')
plt.yticks(size='large', weight='bold')
plt.xlabel('Measurements', size='x-large', weight='bold')
plt.ylabel('abs(BLNK - ray_tracer)', size='x-large', weight='bold')
plt.show()
