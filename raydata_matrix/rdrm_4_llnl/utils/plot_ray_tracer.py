"""
Simple script to plot the ray created by raydata

INPUT
-----
filename, e.g.: raydata.mytest_band01_1
line number: in the following example: 348
 ...
 347  1  467    1    609.42     -0.53      0.00   0.770E+03   0.393E+03
 348   467
 349   6371.000   2.77539   0.00000   5.80000  0.001667  0.999000E+03  0.999000E+03
 350   6356.997   2.77454   0.00084   5.80000  0.001667  0.115165E-01  0.115252E-01
 351   6351.000   2.77418   0.00121   5.80000  0.001667  0.806943E-02  0.807808E-02
 352   6351.000   2.72728   0.00121   6.50000  0.001667  0.838745E-02  0.807807E-02
 353   6337.272   2.72633   0.00216   6.50000  0.001667  0.462769E-02  0.453807E-02
 354   6336.000   2.72624   0.00225   6.50000  0.001667  0.444348E-02  0.436147E-02
 ...
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import sys

filename = os.path.abspath(sys.argv[1])
ICB = 1217.9998
CMB = 3482.0
Rad = 6371.0
print '\n=================='
print 'ICB: %s' % ICB
print 'CMB: %s' % CMB
print 'Rad: %s' % Rad
print '=================='

# Opening the input file
fio = open(filename, 'r')
fi = fio.readlines()

line_to_start = int(raw_input('Enter the line number: (starting from 1)\n'))
number_lines_2_read = int(fi[line_to_start-1].split()[-1])

# r: radius, delta: epicentral distance (radians)
r = []
delta = []
for i in range(line_to_start, line_to_start+number_lines_2_read):
    fi_tmp = fi[i].split()
    if len(fi_tmp) == 7:
        r.append(float(fi_tmp[0]))
        #delta.append(float(fi_tmp[2])*180./np.pi)
        delta.append(float(fi_tmp[2]))

fname = raw_input('Enter a taup_path file: (taup_path -mod iasp91 -h 0 -ph Pdiff -deg 129.81 -gmt)')

plt.ion()
ax = plt.subplot(111, polar=True)

if fname:
    fi = np.loadtxt(fname)
    ax.plot(fi[:, 0]/180.*np.pi, fi[:, 1], 'k--', linewidth=4, zorder=10)

ax.plot(delta, r, color='r', linewidth=3)

# DIRTY! plot ICB and CMB
ax.plot(np.linspace(0, 360, 10000), np.ones(10000)*ICB, lw=1.0, color='k')
ax.plot(np.linspace(0, 360, 10000), np.ones(10000)*CMB, lw=1.0, color='k')

# change the theta direction to clockwise!
ax.set_theta_direction(-1)
ax.set_theta_zero_location('N')

ax.set_thetagrids([0, 45, 90, 135, 180, 225, 270, 315], size='x-large', weight='bold')
# ax.set_rgrids([ICB, CMB], labels=['ICB', 'CMB'], size='x-large', weight='bold')
ax.set_rgrids([ICB, CMB], labels=[' ', ' '], size='x-large', weight='bold')
ax.set_rmax(Rad)
plt.show()

# =========== Non-polar plot

plt.figure()
ax = plt.subplot(111)
ax.plot(np.array(delta)*180./np.pi, r, color='r', linewidth=3, label='Tracer')
if fname:
    ax.plot(fi[:, 0], fi[:, 1], 'k--', linewidth=4, zorder=10, label='TauP')

plt.hlines(CMB, min(np.array(delta)*180./np.pi), max(np.array(delta)*180./np.pi), 
        lw=3, color='k', linestyle='dotted')
plt.xlim(min(np.array(delta)*180./np.pi)-5, max(np.array(delta)*180./np.pi)+5)
plt.xlabel('Distance (deg)', size=24, weight='bold')
plt.ylabel('Depth (km)', size=24, weight='bold')
plt.xticks(size=18, weight='bold')
plt.yticks(size=18, weight='bold')
plt.legend(loc='upper center', prop={'size': 24, 'weight': 'bold'})
plt.show()

raw_input('Press enter to exit...')
