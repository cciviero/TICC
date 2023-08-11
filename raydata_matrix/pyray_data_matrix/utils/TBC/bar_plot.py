"""
Simple script to bins the measurements by depth and 
generate some statistics and figures for that
"""

import glob
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

############### INPUT
add_files = '../OUTPUTS/P_Pdiff/npy_P_Pdiff/Pd_gl*'
# depth intervals
dp_interval = 6
# min/maximum depth
min_depth = 0
max_depth = 600
# Plotting
plot_hist = False
plot_hist_final = True

# Does not really matter if you dont want to change the colors
colors = ['r', 'g', 'b', 'y', 'r', 'g', 'b', 'y']
zord = [1, 2, 3, 4, 5, 6, 7, 8]
############### END INPUT

bins = np.arange(0, 600+dp_interval, dp_interval)

ls_npy = glob.glob(add_files)
ls_npy.sort()
# ===================== TO SELECT NOT ALL THE FILES
ls_npy = ls_npy[0:-2]

hists = []
rans = []
for i in range(len(ls_npy)):
    print ls_npy[i]
    arr_npy = np.load(ls_npy[i])
    histo, ran = np.histogram(arr_npy[:, 29].astype(np.float), bins=bins)
    hists.append(histo)
    rans.append(ran)
    if plot_hist:
        plt.ion()
        fig = plt.figure()
        cs = [colors[i]]*len(ran[:-1])
        ax = fig.add_subplot(111)
        ax.bar(ran[:-1], histo, color=cs, alpha=0.8)
        plt.title('band0%s' % (i+1))
        plt.show()

arr_hists = np.array(hists)
print 'Total number of measurements: %s' % np.sum(arr_hists) 
print 'Each group: \n%s' % np.sum(arr_hists, 1)

if plot_hist_final:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(ran[:-1], np.sum(arr_hists, 0), 'bo-')
    plt.show()

all_num = np.sum(arr_hists, 0)
all_num_depth = np.vstack([all_num, ran[:-1]]).transpose()
all_num_depth_sorted = np.copy(all_num_depth)
all_num_depth_sorted = sorted(all_num_depth_sorted,
                              key=lambda all_num_depth_sorted:
                              float(all_num_depth_sorted[0]))

output_fio = open('depth_number.txt', 'w')
for i in range(len(all_num_depth_sorted)-1, -1, -1):
    select_dp_1 = arr_npy[float(all_num_depth_sorted[i][1]) <=
                          arr_npy[:, 29].astype(np.float)]
    select_dp_2 = select_dp_1[select_dp_1[:, 29].astype(np.float) <
                              float(all_num_depth_sorted[i][1])+dp_interval]
    ev_names = []
    for j in range(len(select_dp_2)):
        ev_names.append(select_dp_2[j, 23])
    str_ev_names = str(ev_names)[1:-1]
    output_fio.writelines('%s,%s,%s\n' % (all_num_depth_sorted[i][1],
                                          all_num_depth_sorted[i][0],
                                          str_ev_names))
output_fio.close()

#fig = plt.figure()
#ax = fig.add_subplot(111,projection='3d')
#for i in range(len(ls_npy)):
#    print ls_npy[i]
#    arr_npy = np.load(ls_npy[i])
#    hist, ran = np.histogram(arr_npy[:, 29].astype(np.float),
#                             bins=np.arange(0, 600))
#    cs = [colors[i]]*len(ran[:-1])
#    ax.bar(ran[:-1], hist, zs=zord[i], zdir='y', color=cs, alpha=0.8)
