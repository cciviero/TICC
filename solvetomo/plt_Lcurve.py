# /usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  plot_tomo_vtk.py
#   Purpose:   plot lcurve and some analysis plots
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

import matplotlib

matplotlib.use("Agg")

import glob
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import os
import itertools

# ------------------- INPUT -----------------------------

min_trade = 1
max_trade = 11

# path where your sol, solx, cor etc files are...
# path = ['/net/siglochnas1/volume1/data1/mariat/RHUM-RUM/measurement_redo/P_evolution/resolutiontest_IPdRa_v2/inversion_dir_ISC_Pdiff_RRall']


# path = glob.glob('/Users/maria/PhD/Codes/__TEST_DATA/OUTPUT/'
#              'P_32_85_STFnc_IASP91_2s_03.06.2017_the17/'
#              'rdrm_cc0.8_crust*/inversion_dir/sd_0.3')

# path = glob.glob('/Volumes/AB_EXTRA/2017_inversions/grid_iter4/PPdiPP_OX_LMU_ISC_RRP_L_0_3/rsync_me')

path = glob.glob('/mnt/seismodata2/MT/RAYTRACER/TEST_OUTPUT/TESTME_modelfiles_12.08.2021/inv_12.08.2021/sd_0.3')
path.sort()
print path

tradeoffs_name = ['0.3']
# tradeoffs_name  = ['crust1', 'crust1+', 'crust2']
# tradeoffs_name = ['0.3', '0.5', '0.7', '0.9', '1.2']
# tradeoffs_name = ['0.3 - no corr', '0.3 - *2', '0.3 - *0.1', '1.0 - *0.1', '0.3 - *0.3', '1.0 - *0.3', '1.0 - *0.5']

# DEFAULT
output_dir = path[0]
# import ipdb; ipdb.set_trace()
# output_dir = '/net/siglochnas1/volume1/data1/mariat/RHUM-RUM/measurement_redo/S_19.04.2018'
# output_dir = '/Users/maria/PhD/Codes/__TEST_DATA/OUTPUT/P_32_85_STFnc_IASP91_2s_03.06.2017_the17'

# tradeoffs_name = [
#         '0.001',
#         '0.01',
#         '0.1',
#         '0.2',
#         '0.3',
#         '0.4',
#         '0.5',
#         '1.0',
#         'sn\_0',
#         'sn\_3',
#         'sn\_6',
#         's\_6\_20']


colors = ["r", "orange", "g", "b", "indigo", "violet",
          "r", "orange", "g", "b", "indigo", "violet",
          "r", "orange", "g", "b", "indigo", "violet"]
line_styles = ['-', '-', '-', '-', '-', '-',
               '--', '--', '--', '--', '--', '--',
               ':', ':', ':', ':', ':', ':']

# >>>>> Analysis plots settings
trade_range = range(min_trade, max_trade)

# 3rd plot: dVp/Vp
dpar = 0.01  # [%] uncertainty that we assign to dlnvp (assemblematrix)
vp_binsbin = 61  # amount of bins
vp_binsrange = [-0.03, 0.03]  # boundaries for the bins

# 4th plot: various corrections


# 5th plot: residuals
nrow = 4440585  # make sure you check this value in e.g diagnostics.so02.08.201
# nrow = 11666570
res_binsbin = 61
res_binsrange = [-5, 5]

# ------------------- 1st Plot -----------------------------
plt.ioff()

plt.figure(figsize=(15, 10))

for j, tradeoff in enumerate(path):
    print "%s -- %s" % (j, tradeoff)
    xytradeoff = glob.glob(os.path.join(tradeoff, 'xy.tradeoff.*'))
    # import ipdb; ipdb.set_trace()
    xytradeoff_name = xytradeoff[0].split('.')[-1]
    troff_arr = np.loadtxt(xytradeoff[0], comments='#')

    plt.plot(troff_arr[min_trade:max_trade, 1], troff_arr[min_trade:max_trade, 0],
             lw=3, color=colors[j], ls=line_styles[j],
             label=tradeoffs_name[j])

    plt.scatter(troff_arr[min_trade:max_trade, 1], troff_arr[min_trade:max_trade, 0],
                lw=3, color=colors[j])

    plt.axvline(1, 0, 1, c='k', ls='--', lw=3)

    for i in range(min_trade, min_trade + len(troff_arr[min_trade:max_trade, 1])):
        plt.text(troff_arr[i, 1], troff_arr[i, 0],
                 '%s' % (i + 1), color=colors[j], size='x-large')

plt.ylabel(r'$\|m\|_{2}$', size=42, weight='bold')
plt.xlabel(r'$\chi^2$', size=42, weight='bold')
plt.xticks(size=32, weight='bold')
plt.yticks(size=32, weight='bold')
plt.grid()
plt.title('L-curve', size=52, weight='bold')
plt.tight_layout()
plt.legend(prop={'size': '24'}, bbox_to_anchor=(1, 1))
plt.xlim(-1)
plt.ylim(-1)
# import ipdb; ipdb.set_trace()
plt.savefig(os.path.join(output_dir, 'L_curve.png'))
# import ipdb; ipdb.set_trace()


# ------------------- 2nd Plot -----------------------------
plt.figure(figsize=(15, 10))
for j, tradeoff in enumerate(path):
    print "%s -- %s" % (j, tradeoff)
    texsolve = glob.glob(os.path.join(tradeoff, 'tex.solve.*'))
    texsolve_name = texsolve[0].split('.')[-1]
    troff_arr = np.loadtxt(os.path.join(tradeoff,
                                        'tex.solve.%s' % texsolve_name),
                           skiprows=7, usecols=(0, 1, 2, 3, 4, 5, 6), dtype='object')
    plt.plot(troff_arr[min_trade:max_trade, 2].astype(np.float),
             troff_arr[min_trade:max_trade, 4].astype(np.float) * dpar * 100,
             lw=3, color=colors[j], ls=line_styles[j],
             label=tradeoffs_name[j])
    plt.scatter(troff_arr[min_trade:max_trade, 2].astype(np.float),
                troff_arr[min_trade:max_trade, 4].astype(np.float) * dpar * 100,
                lw=3, color=colors[j])
    plt.axvline(1, 0, 1, c='k', ls='--', lw=1)
    for i in range(min_trade, min_trade + len(troff_arr[min_trade:max_trade, 1])):
        plt.text(float(troff_arr[i, 2]),
                 float(troff_arr[i, 4]) * dpar * 100,
                 '%s' % (i + 1), color=colors[j], size='x-large')
plt.ylabel(r'$\|m\|_{\inf}$ (%)', size=42, weight='bold')
plt.xlabel(r'$\chi^2$', size=42, weight='bold')
plt.title('L-curve (maximum model-parameter values)', size=20, weight='bold')
plt.xticks(size=32, weight='bold')
plt.yticks(size=32, weight='bold')
plt.grid()
plt.tight_layout()
plt.legend(prop={'size': '24'}, bbox_to_anchor=(1, 1))
# plt.xlim(1, 3.5)
plt.savefig(os.path.join(output_dir, 'm_inf_curve'))
plt.show()

# ------------------- 3rd Plot -----------------------------

bins_bin = vp_binsbin
sol_path = glob.glob(os.path.join(path[0], 'solx.*'))

for i, spath in enumerate(sol_path[::-1]):
    print 'Working on dVp/Vp'
    plt.figure(figsize=(15, 10))
    dvpvp = np.loadtxt(spath, usecols=(3,), dtype=float, skiprows=2)

    # bins_range = [np.min(dvpvp), np.max(dvpvp)]
    bins_range = vp_binsrange
    hist, bins = np.histogram(dvpvp, bins=bins_bin, range=bins_range)
    # width = 0.7 * (bins[1] - bins[0])
    # center = (bins[:-1] + bins[1:]) / 2
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    # import ipdb; ipdb.set_trace()

    # import ipdb; ipdb.set_trace()
    plt.axvline(dpar, 0, 1, c='r', ls='--', lw=2)
    plt.axvline(-dpar, 0, 1, c='r', ls='--', lw=2)

    # import ipdb; ipdb.set_trace()
    # plt.xlim(-dpar-2*np.std(dvpvp), dpar+2*np.std(dvpvp))
    # plt.xlim(-dpar-0.3*dpar, dpar+0.3*dpar)
    # plt.xlim(np.min(dvpvp), np.max(dvpvp))
    # plt.xlim(-0.3, 0.3)
    if i == 0:
        # import ipdb; ipdb.set_trace()
        ylim_value = 1.05 * np.max(hist)
    plt.ylim(0, ylim_value)
    plt.title('dVp/Vp [%]', weight='bold', size=12)
    plt.xlabel('prior uncertainity: %s' % dpar, weight='bold')
    trade_nr = trade_range[::-1][i]
    plt.bar(center, hist, align='center', width=width)
    plt.savefig(os.path.join(output_dir, 'dVpVp_solx_%s' % trade_nr))

# ------------------- 4th Plot -----------------------------


cor_path = glob.glob(os.path.join(path[0], 'cor.*'))

plot_names = ['Hypocentre correction [km]', 'Station correction t-P [sec]',
              'Station correction dlnA-P [dimless]', 'Origin time correction [sec]',
              'Event correction dlnA-P [dimless]', 'Station correction t-S [sec]',
              'Station correction dlnA-S [dimless]']

units = ['[km]', '[sec]', '[dimless]', '[sec]', '[dimless]', '[sec]', '[dimless]']

save_name = ['hypocent_corr', 'sta_corr_t-P', 'sta_corr_dlnA-P', 'origtime_corr',
             'ev_corr_dlnA-P', 'sta_corr_t-S', 'sta_corr_dlnA-P']

bins_size = 51
bin_range_multiplier = [2, 3, 2, 3, 3, 3, 2]

ylimits = []
for i, cpath in enumerate(cor_path[::-1]):
    f = open(cpath)
    lines = f.readlines()
    l2 = np.array(lines[1].split(), dtype=int)  # switch on or not
    l3 = np.array(lines[2].split(), dtype=float)  # which values
    l4 = np.array(lines[3].split(), dtype=int)  # how many values are down the list

    values = np.loadtxt(cpath, usecols=(0,), dtype=float, skiprows=6)

    for j, jval in enumerate(l2):
        plt.figure(figsize=(15, 10))

        if jval > 0:

            print 'Working on: %s' % plot_names[j]
            corr_list = values[:l4[j]]
            corr_value = l3[j]

            bins_range = [-l3[j] * bin_range_multiplier[j], l3[j] * bin_range_multiplier[j]]
            hist, bins = np.histogram(corr_list, bins=bins_size, range=bins_range)
            width = 0.7 * (bins[1] - bins[0])
            center = (bins[:-1] + bins[1:]) / 2

            plt.axvline(corr_value, 0, 1, c='r', ls='--', lw=2)
            plt.axvline(-corr_value, 0, 1, c='r', ls='--', lw=2)

            trade_nr = trade_range[::-1][i]
            plt.bar(center, hist, align='center', width=width)
            plt.title(plot_names[j], weight='bold', size=12)
            plt.xlabel('prior uncertainity: %s %s' % (corr_value, units[j]), weight='bold')

            if i == 0:
                ylimits.append(np.max(hist))

            plt.ylim(0, ylimits[j])
            # plt.xlim(-corr_value - 2*np.std(corr_list), corr_value + 2*np.std(corr_list))
            plt.savefig(os.path.join(output_dir, '%s_cor%s' % (save_name[j], trade_nr)))
        else:
            if i == 0:
                ylimits.append(0)

print ylimits
print 'DONE!'

# ------------------- 5th Plot -----------------------------

bins_bin = res_binsbin

res_path = glob.glob(os.path.join(path[0], 'res.*'))

cur_dir = os.getcwd()

os.chdir(glob.glob(os.path.join(path[0]))[0])
os.chdir('..')
path2inverion = os.getcwd()
aux_path = glob.glob(os.path.join(path2inverion, 'aux.*'))
# aux_path = glob.glob('aux.*')
# import ipdb; ipdb.set_trace()
with open(aux_path[0]) as f_input:
    aux_bands = np.loadtxt(itertools.islice(f_input, 0, nrow), usecols=(9,), dtype=int)

os.chdir(cur_dir)
# import ipdb; ipdb.set_trace()
set_aux_bands = set(aux_bands)
nr_bands = len(set_aux_bands)
cmap = plt.get_cmap('jet_r')
colors = cmap(np.linspace(0, 1.0, nr_bands))

max_list = []
for i, rpath in enumerate(res_path):
    print 'working on:', rpath
    plt.figure(figsize=(15, 10))
    for j, band in enumerate(set_aux_bands):
        # plot band0/ISC data separately
        if band == 0:
            plt.figure(figsize=(15, 10))
            filter_list = aux_bands == band
            values = np.loadtxt(rpath, usecols=(0,), dtype=float)

            try:
                filtered_values = values[filter_list]
                bins_range = res_binsrange
            except Exception, exp:
                print 'Ups:', exp
                import ipdb;

                ipdb.set_trace()

            hist, bins = np.histogram(filtered_values, bins=bins_bin, range=bins_range)
            width = 0.7 * (bins[1] - bins[0])
            center = (bins[:-1] + bins[1:]) / 2
            plt.axvline(dpar, 0, 1, c='r', ls='--', lw=2)
            plt.axvline(dpar, 0, 1, c='r', ls='--', lw=2)

            trade_nr = trade_range[i] + 1

            plt.bar(center, hist, align='center', width=width, color=colors[j],
                    label='band0%s' % (band), alpha=0.6)
            plt.title('data resuduals d-A*m', weight='bold', size=12)
            plt.xlabel('res [sec]', weight='bold')
            plt.legend()

            # plt.xlim(-dpar - 0.8*dpar, dpar + 0.8*dpar)
            # max_list.append(np.max(hist))
            # print max_list
            if j == 0 and i == 0:
                ylim_value_0 = np.max(hist)
            plt.ylim(0, ylim_value_0)
            plt.savefig(os.path.join(output_dir, 'res_B0_%s' % trade_nr))
            plt.clf()

        else:
            filter_list = aux_bands == band
            values = np.loadtxt(rpath, usecols=(0,), dtype=float)

            try:
                filtered_values = values[filter_list]
                bins_range = res_binsrange
            except Exception, exp:
                print 'Ups:', exp
                import ipdb;

                ipdb.set_trace()

            hist, bins = np.histogram(filtered_values, bins=bins_bin, range=bins_range)
            width = 0.7 * (bins[1] - bins[0])
            center = (bins[:-1] + bins[1:]) / 2
            plt.axvline(dpar, 0, 1, c='r', ls='--', lw=2)
            plt.axvline(dpar, 0, 1, c='r', ls='--', lw=2)

            trade_nr = trade_range[i] + 1

            plt.bar(center, hist, align='center', width=width, color=colors[j],
                    label='band0%s' % (band), alpha=0.6)
            plt.title('data resuduals d-A*m', weight='bold', size=12)
            plt.xlabel('res [sec]', weight='bold')
            plt.legend()

            # plt.xlim(-dpar - 0.8*dpar, dpar + 0.8*dpar)
            # max_list.append(np.max(hist))
            # print max_list
            if j == 1 and i == 0:
                ylim_value = np.max(hist)
            plt.ylim(0, ylim_value)
    plt.savefig(os.path.join(output_dir, 'res_B_%s' % trade_nr))

print 'Done...'

