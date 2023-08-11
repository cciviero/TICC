#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  grouping_interval.py
#   Purpose:   Generating groups of nodes based on Sensitivity Kernel
#               values by using intervals
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import matplotlib.pyplot as plt
import numpy as np
import sys

# ----------------- INPUT
# Absolute (abs) or relative (rel) value?
abs_rel = 'abs'
# Earth's radius in km
eradius = 6371.
# -----------------

# read the columndensity file
colden = np.loadtxt(sys.argv[1], skiprows=2)
colden = np.c_[colden, np.arange(0, len(colden))]

# Change the x,y,z to lat,lon,dep
rxy = np.sqrt(colden[:, 0]**2 + colden[:, 1]**2)
lon = np.arctan2(colden[:, 1], colden[:, 0])*180./np.pi
lat = np.arctan2(colden[:, 2], rxy)*180./np.pi
dep = eradius - np.sqrt(colden[:, 0]**2 + colden[:, 1]**2 + colden[:, 2]**2)

# Choose the relative or absolute value
if abs_rel == 'abs':
    val = colden[:, 4]
elif abs_rel == 'rel':
    val = colden[:, 3]
else:
    sys.exit('%s has not implemented!' % abs_rel)

# Just as starting values
bin_num = 10
bin_num_last = int(bin_num)
while bin_num:
    plt.subplot(2, 1, 1)
    plt.plot(dep, val, '.', lw=3)
    plt.xlabel("Depth (km)", size=18, weight='bold')
    plt.ylabel("Sensitivity Kernel", size=18, weight='bold')

    bins = np.linspace(0, np.max(val), bin_num)
    plt.axvline(2889, color='k', ls='--')
    for i in bins:
        plt.axhline(i, color='r', ls='--')

    val[val < 0.000001] = 0.000001

    vallog = np.log10(val)

    plt.subplot(2, 1, 2)
    plt.plot(dep, vallog, '.', lw=3)
    plt.xlabel("Depth (km)", size=18, weight='bold')
    plt.ylabel("Sensitivity Kernel (log10)", size=18, weight='bold')

    bins = np.linspace(0, np.max(vallog), bin_num)
    plt.axvline(2889, color='k', ls='--')
    for i in bins:
        plt.axhline(i, color='r', ls='--')
    plt.show()

    # Digitize the logarithmic values
    valdig = np.digitize(vallog, bins, right=False)

    bin_num = raw_input('Enter the number of bins... (f: finalize)\n')
    if bin_num.lower() == 'f':
        bin_num = False
    else:
        bin_num = int(bin_num)
        bin_num_last = int(bin_num)

print "All nodes with depth >= 3000 are in group: %s" % (bin_num_last + 1)
valdig[dep >= 3000] = bin_num_last + 1
print "Create a histogram map..."
hist, bins_hist = np.histogram(valdig, bins=range(0, bin_num_last+3))
width = 0.7 * (bins_hist[1] - bins_hist[0])
center = (bins_hist[:-1] + bins_hist[1:]) / 2
plt.subplot(2, 1, 1)
plt.bar(center, hist, align='center', width=width)
plt.xlabel("Digitized values", size=18, weight='bold')
plt.ylabel("#Vertices", size=18, weight='bold')
plt.xlim(0, bin_num_last+2)
hi = []
hval = []
for i in range(min(valdig), max(valdig)+1):
    hi.append(i)
    hval.append(np.sum(vallog[valdig == i])/len(vallog[valdig == i]))

hval_arr = np.array(hval)
hval_arr = hval_arr[-2] - hval_arr + 1.
hval_arr[-1] = hval_arr[0]

plt.subplot(2, 1, 2)
plt.plot(hi, hval_arr, 'bo')
plt.xlim(0, bin_num_last+2)
plt.show()

# ================= group_density
wf_array = np.c_[colden[:, 0], colden[:, 1], colden[:, 2], valdig]
np.savetxt('group_density', wf_array)

header_fio = open(sys.argv[1], 'r')
header_fi = header_fio.readlines()

with open('group_density', 'r+') as f:
    content = f.read()
    f.seek(0, 0)
    f.write(header_fi[0] + header_fi[1] + content)

# ================= group_value
wf_array = np.c_[colden[:, 0], colden[:, 1], colden[:, 2], val]
np.savetxt('group_value', wf_array)

with open('group_value', 'r+') as f:
    content = f.read()
    f.seek(0, 0)
    f.write(header_fi[0] + header_fi[1] + content)

damp_val = np.copy(valdig)
damp_val = damp_val.astype(float)
for i in hi:
    damp_val[valdig == i] = hval_arr[i-hi[0]]

import ipdb; ipdb.set_trace()
# ================= damp_density
wf_array = np.c_[colden[:, 0], colden[:, 1], colden[:, 2], damp_val]
np.savetxt('damp_density', wf_array)

with open('damp_density', 'r+') as f:
    content = f.read()
    f.seek(0, 0)
    f.write(header_fi[0] + header_fi[1] + content)
