#!/usr/bin/env python
# -*- coding: utf-8 -*-


# ################### WARNING
# USE WITH CAUTION
# ################### WARNING


# -------------------------------------------------------------------
#   Filename:  grouping_kdtree.py
#   Purpose:   Generating groups of nodes based on Sensitivity Kernel
#               values by using KDTREE
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

from kh_utils import SphericalNearestNeighbour

# ----------------- INPUT
# Absolute (abs) or relative (rel) value?
abs_rel = 'abs'
# some parameters for searching algorithm
num_neigh = 5000
den_thre = 10.
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

oindx = colden[:, 5].astype(int)

# initialization for the while loop
all_acc_indx = []
all_acc_grp = []
reg_dict = {}
read_indx = 0
icount = 0
contin = True
exist_dict = False
while contin:
    kd = SphericalNearestNeighbour(lat, lon, dep, eradius=eradius)
    qindx = kd.query(np.array([lat[read_indx]]),
                     np.array([lon[read_indx]]),
                     np.array([dep[read_indx]]), num_neigh)

    # extract only column densities of the neighbors
    q_colden = val[qindx[1][0]]
    # calculate the mean
    # q_mean = np.mean(q_colden)
    # deviation from the mean
    # q_dev = abs(q_colden) - q_mean
    compare_val = q_colden[read_indx]
    q_dev = q_colden - compare_val

    # determine the accepted neighbors
    acc_indx = qindx[1][0, abs(q_dev) < den_thre]
    oacc_indx = oindx[qindx[1][0, abs(q_dev) < den_thre]]

    if len(reg_dict.keys()) < 1:
        add_grp = icount
        reg_dict[icount] = [q_colden[read_indx], oacc_indx]
    else:
        for i in reg_dict.keys():
            if abs(compare_val-reg_dict[i][0]) < den_thre:
                print abs(compare_val-reg_dict[i][0])
                add_grp = i
                reg_dict[add_grp][1] = np.append(reg_dict[add_grp][1], oacc_indx)
                exist_dict = True
                break
            else:
                add_grp = icount
                reg_dict[icount] = [q_colden[read_indx], oacc_indx]
                exist_dict = False
    # region dictionary
    # reg_dict[icount] = [q_mean, oacc_indx]
    all_acc_indx.extend(list(oacc_indx))
    grps = list(np.zeros(len(oacc_indx)).astype(int)+add_grp)
    all_acc_grp.extend(grps)

    len_lat_1 = len(lat)
    mask = np.ones(len(lat), dtype=bool)
    mask[acc_indx] = False
    # removing the grouped lat, lon, dep, val and index
    lat = lat[mask]
    lon = lon[mask]
    dep = dep[mask]
    val = val[mask]
    oindx = oindx[mask]

    len_lat_2 = len(lat)

    if exist_dict:
        icount -= 1

    icount += 1

    # if 5000 < len_lat_2 < num_neigh:
    #     num_neigh -= 5000
    #     print(len(lat))
    #     continue
    if len(lat) < num_neigh:
        reg_dict[icount] = [abs(np.mean(val)), oindx]
        all_acc_indx.extend(list(oindx))
        grps = list(np.zeros(len(oindx)).astype(int)+icount)
        all_acc_grp.extend(grps)
        contin = False
    if len_lat_1 == len_lat_2:
        # if two consecutive searches resulted in the same number of
        # latitude and longitude, we go to the next node...otherwise,
        # back to the first node for consistency
        read_indx += 1
    else:
        read_indx = 0
    print len(lat)

grp_sorted = [x for (y, x) in sorted(zip(all_acc_indx, all_acc_grp))]
wf_array = np.c_[colden[:, 0], colden[:, 1], colden[:, 2], grp_sorted]
np.savetxt('group_den', wf_array)

# from mpl_toolkits.mplot3d import Axes3D
# for i in reg_dict.keys():
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     colden_sel = colden[reg_dict[i][1], :]
#     xs = colden_sel[:, 0]
#     ys = colden_sel[:, 1]
#     zs = colden_sel[:, 2]
#     c = np.ones(len(xs)) * reg_dict[i][0]
#     ax.scatter(xs, ys, zs, c=c)
#     plt.show()
