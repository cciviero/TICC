"""
Function to calculate the backazimuth and
R, T components of horizontal components
"""

import matplotlib.pyplot as plt
import numpy as np

from obspy import read
from obspy.core.util import gps2DistAzimuth
from obspy.signal import rotate

import subprocess

# ---------------- INPUT
smgr_add = './grf.II.AAK.10.x00.BH*'
bg_model = 'iasp91'
req_phase = 'S'
# ---------------- END INPUT


st = read(smgr_add)
st.taper()
st.filter('lowpass', freq=0.5)
st.filter('highpass', freq=0.012)


tr_N = st.select(channel='BHN')[0]
tr_E = st.select(channel='BHE')[0]

(dist, azi, bazi) = gps2DistAzimuth(tr_N.stats.sac.evla,
                                    tr_N.stats.sac.evlo,
                                    tr_N.stats.sac.stla,
                                    tr_N.stats.sac.stlo)

tt = False

# --------------- TAUP
taup_process = subprocess.Popen(['taup_time', '-mod', bg_model, '-time',
                                 '-h', str(tr_N.stats.sac.evdp),
                                 '-ph', req_phase,
                                 '-deg', str(dist)],
                                stdout=subprocess.PIPE)

tt_raw = taup_process.communicate()[0]
try:
    tt = tt_raw.split('\n')[0].split()[-1]
    tt = float(tt)
except Exception, e:
    print 'Requested phase: %s ... ERROR: %s' % (req_phase, e)
    tt = False
# --------------- END TAUP

(tr_data_R, tr_data_T) = rotate.rotate_NE_RT(tr_N.data, tr_E.data, bazi)

tr_R = tr_N.copy()
tr_T = tr_N.copy()

tr_R.data = tr_data_R
tr_T.data = tr_data_T

# --------------- Particle movement
tr_N_plt = tr_N.copy()
tr_E_plt = tr_E.copy()

tr_N_plt.trim(tr_N_plt.stats.starttime, tr_N_plt.stats.starttime+925)
tr_E_plt.trim(tr_E_plt.stats.starttime, tr_E_plt.stats.starttime+925)
plt.figure()
plt.plot(tr_E_plt.data, tr_N_plt.data, '.')
plt.xticks(size=18, weight='bold')
plt.yticks(size=18, weight='bold')
plt.xlabel('E component', size=22, weight='bold')
plt.ylabel('N component', size=22, weight='bold')
plt.show()
# --------------- END Particle movement

plt.ion()

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(np.linspace(0, (tr_N.stats.npts-1)/tr_N.stats.sampling_rate,
                     tr_N.stats.npts),
         tr_N.data, lw=3)
if tt:
    plt.vlines(tt, min(tr_N.data), max(tr_N.data), color='k')
plt.xticks(size=18, weight='bold')
plt.yticks(size=18, weight='bold')
plt.title('N component', size=22, weight='bold')

plt.subplot(2, 1, 2)
plt.plot(np.linspace(0, (tr_E.stats.npts-1)/tr_E.stats.sampling_rate,
                     tr_E.stats.npts),
         tr_E.data, lw=3, color='r')
if tt:
    plt.vlines(tt, min(tr_E.data), max(tr_E.data), color='k')
plt.xticks(size=18, weight='bold')
plt.yticks(size=18, weight='bold')
plt.title('E component', size=22, weight='bold')

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(np.linspace(0, (tr_R.stats.npts-1)/tr_R.stats.sampling_rate,
                     tr_R.stats.npts),
         tr_R.data, lw=3)
if tt:
    plt.vlines(tt, min(tr_R.data), max(tr_R.data), color='k')
plt.xticks(size=18, weight='bold')
plt.yticks(size=18, weight='bold')
plt.title('R component', size=22, weight='bold')

plt.subplot(2, 1, 2)
plt.plot(np.linspace(0, (tr_T.stats.npts-1)/tr_T.stats.sampling_rate,
                     tr_T.stats.npts),
         tr_T.data, lw=3, color='r')
if tt:
    plt.vlines(tt, min(tr_T.data), max(tr_T.data), color='k')
plt.xticks(size=18, weight='bold')
plt.yticks(size=18, weight='bold')
plt.title('T component', size=22, weight='bold')
plt.show()

raw_input('Press enter to quite...')
