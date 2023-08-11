"""
Plot filters presented in a filter file like: bpf.omega_m
"""
import matplotlib.pyplot as plt
import numpy as np

# ========== INPUT
filename = '../src_raydata_raymatrix/files/bpf.omega_m'
#filename = '../src_raydata_raymatrix/files/gauss_filter_51_15'
# Define the frequencies in each block:
freqs = [1./30, 1/21.2, 1/15., 1/10.6, 1/7.5, 1/5.3, 1/3.7, 1/2.7]
#freqs = [1./51, 1./30., 1./22.5, 1./15]
# ========== END INPUT

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

fio = open(filename, 'r')
fi = fio.readlines()
num_filts = int(fi[0].split()[0])

line_num = 0
freq_filt_all = []
gain_filt_all = []
for i in range(num_filts):
    line_num += 1
    print line_num
    filt_freqs_num = int(fi[line_num].split()[0])
    freq_filt = []
    gain_filt = []
    for j in range(filt_freqs_num):
        line_num += 1
        freq_filt.append(float(fi[line_num].split()[0]))
        gain_filt.append(float(fi[line_num].split()[1]))
    print line_num
    ax.plot(freq_filt, gain_filt, 'k', lw=4)

for fr in freqs:
    ax.plot([2*np.pi*fr, 2*np.pi*fr], [0, 1], 'r--', lw=2)
ax.set_xscale('log')
plt.xlabel('rad. freq (2*pi*Hz)', size='xx-large', weight='bold')
plt.ylabel('Absolute gain', size='xx-large', weight='bold')
plt.xticks(size='x-large', weight='bold')
plt.yticks(size='x-large', weight='bold')
plt.show()
