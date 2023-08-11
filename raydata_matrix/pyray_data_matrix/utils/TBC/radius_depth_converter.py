"""
Change the format of station files in makemymesh from radius to depth
Usually the station file comes from taup:
    taup_path -ph P -h 0 -sta 40 0 -evt 0 0
and first it should ordered as:
    lat   lon    radius
"""

import numpy as np
import sys

# ============= INPUT
radius = 6371.
# =============

filename = sys.argv[1]

filename_np = np.loadtxt(filename, skiprows=1)
fio = open('./corrected.txt', 'w')

filename_np[:, 1] = radius - filename_np[:, 1]
for f_np in filename_np:
    fio.writelines('%s   %s   %s\n' % (f_np[2], f_np[3], f_np[1]))
fio.close()
