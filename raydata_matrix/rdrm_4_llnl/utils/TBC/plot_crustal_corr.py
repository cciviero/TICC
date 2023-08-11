#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  plot_crustal_corr.py
#   Purpose:   simple plot for overview of common corrections for one
#              event, all stations, one band
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import sys
import os

class bc:
    grey = '\033[90m'     # broing information
    yellow = '\033[93m'   # FYI
    orange = '\033[0;33m'  # Warning

    lred = '\033[1;31m'  # there is smoke
    red = '\033[91m'    # fire!
    dred = '\033[2;31m'  # Everything is on fire

    lblue = '\033[1;34m'
    blue = '\033[94m'
    dblue = '\033[2;34m'

    lgreen = '\033[1;32m'  # all is normal
    green = '\033[92m'    # something else
    dgreen = '\033[2;32m'  # even more interesting

    lmagenta = '\033[1;35m'
    magenta = '\033[95m'  # for title
    dmagenta = '\033[2;35m'

    cyan = '\033[96m'    # system time
    white = '\033[97m'   # final time

    black = '\033[0;30m'

    end = '\033[0m'
    bold = '\033[1m'
    under = '\033[4m'

# -----------------------------------------------------------------------
# ---------------------------- Input part -------------------------------
# -----------------------------------------------------------------------

# ====>>>>
input_file = sys.argv[1]

# ====>>>> read file; column information:
#  0  1 segment
#  1  2 ievt
#  2  3 idate
#  3  4 sta_code
#  4  5 net
#  5  6 sta_lat
#  6  7 sta_lon
#  7  8 depth in radius
#  8  9 source lat
#  9 10 source lon
# 10 11 ellipticity corr
# 11 12 crust thickness (NEW!!!)
# 12 13 crust corr
# 13 14 elevation      (NEW!!!)
# 14 15 elevation corr
# 15 16 attenuation
# 16 17 dist [deg]
# 17 18 tbmfi
# 18 19 trtime(iar)
# 19 20 ???
# 20 21 tobs
# 21 22 source depth
dt = np.dtype([('segment', int), ('ievt', int), ('idate', int), ('station code', int), ('net', '|S3'), ('sta_lat', float),
               ('sta_lon', float), ('depth_radius', float), ('source_lat', float), ('source_lon', float),
               ('ellipticity_corr', float), ('crust_thick', float), ('crust_corr', float), ('elevation', float),
               ('elevation_corr', float), ('attenuation', float),
               ('dist_deg', float), ('tbmfi', float), ('trtimeiar', float), ('???', float),
               ('tobs', float), ('source depth', float)])
ccorr_arr = np.genfromtxt(input_file, delimiter=',', dtype=dt)

# -----------------------------------------------------------------------
# ------------------------------ Main part ------------------------------
# -----------------------------------------------------------------------

# get unique events
uni_ev = list(set(ccorr_arr['ievt']))
# list_uni_ev = [int(x) for x in uni_ev]
print 'Amount of events:', len(uni_ev)

# ====>>>> This part sorts for the different events for one band
master = []
for ev in uni_ev:
    temp = []
    for station in ccorr_arr:
        if ev == station[1]:
            # print ev, station[1]
            temp.append(station)
            # print 'here'
            # print temp
    master.append(temp)

for event in master:
    plt.ioff()
    plt.figure(figsize=(15, 10))
    ev = np.array(event)
    plt.suptitle('ievt: %s' % ev[0][1], weight='bold', size=15)

    print bc.green + bc.bold + 'Working on Event(ID): %s' % ev[0][1] + bc.end

    # ====>>>> calculating ddif:
    dtheor = ev['trtimeiar'] + ev['ellipticity_corr']
    dtheor = dtheor + ev['crust_corr']
    dtheor = dtheor + ev['elevation_corr']

    ddif = ev['tobs'] - dtheor

    # ====>>>> 1st plot: total time -> ddif = dobs - dtheor
    plt.subplot(3, 2, 1)
    mymap = Basemap(projection='robin', lon_0=0, lat_0=0)

    mymap.drawcoastlines(color='black', linewidth=1, zorder=15)
    mymap.drawmeridians(np.arange(0, 360, 30), zorder=20)
    mymap.drawparallels(np.arange(-90, 90, 30), zorder=20)

    # only define once the source coordinates
    sx, sy = mymap(ev[0][8], ev[0][7])
    mymap.scatter(sx, sy, s=80, marker='*', c='y', lw=1, zorder=30, edgecolors='g')

    # get the boundaries for the colorbar (same for the second subplot - only defined once)
    col_abs = np.mean(ddif) + np.std(ddif)

    x, y = mymap(ev['sta_lon'], ev['sta_lat'])
    mymap.scatter(x, y, c=ddif, cmap=plt.cm.jet, s=20, lw=0, zorder=10,
                  edgecolors=None, vmin=-col_abs, vmax=col_abs, alpha=0.8)

    cbar = plt.colorbar(orientation='horizontal', shrink=0.5)
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.set_label('time [s]')
    plt.title('Total Time: ddif = dobs - dtheor', weight='bold')

    # ====>>>> calculating dt = tobs - tpred
    dt = ev['tobs'] - ev['trtimeiar']

    # ====>>>> 2nd plot: dt = tobs - tpred
    plt.subplot(3, 2, 2)
    mymap = Basemap(projection='robin', lon_0=0, lat_0=0)
    mymap.drawcoastlines(color='black', linewidth=1, zorder=15)
    mymap.drawmeridians(np.arange(0, 360, 30), zorder=20)
    mymap.drawparallels(np.arange(-90, 90, 30), zorder=20)

    mymap.scatter(sx, sy, s=80, marker='*', c='y', lw=1, zorder=30, edgecolors='g')

    x, y = mymap(ev['sta_lon'], ev['sta_lat'])
    mymap.scatter(x, y, c=dt, cmap=plt.cm.jet, s=20, lw=0, zorder=10, edgecolors=None,
                  vmin=-col_abs, vmax=col_abs, alpha=0.8)

    cbar = plt.colorbar(orientation='horizontal', shrink=0.5)
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.set_label('time [s]')
    plt.title('dt = tobs - tpred', weight='bold')

    # ====>>>> 3rd plot: ellipticity correction
    plt.subplot(3, 2, 3)
    mymap = Basemap(projection='robin', lon_0=0, lat_0=0)
    mymap.drawcoastlines(color='black', linewidth=1, zorder=15)
    mymap.drawmeridians(np.arange(0, 360, 30), zorder=20)
    mymap.drawparallels(np.arange(-90, 90, 30), zorder=20)

    mymap.scatter(sx, sy, s=80, marker='*', c='y', lw=1, zorder=30, edgecolors='g')

    # get the boundaries for the colorbar
    col_abs = np.mean(ev['ellipticity_corr']) + np.std(ev['ellipticity_corr'])

    x, y = mymap(ev['sta_lon'], ev['sta_lat'])
    mymap.scatter(x, y, c=ev['ellipticity_corr'], cmap=plt.cm.jet, s=20, lw=0, zorder=10, edgecolors=None,
                  vmin=-col_abs, vmax=col_abs, alpha=0.8)

    cbar = plt.colorbar(orientation='horizontal', shrink=0.5)
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.set_label('time [s]')
    plt.title('Ellipticity Correction', weight='bold')

    # ====>>>> 4th plot: Crustal correction
    plt.subplot(3, 2, 4)
    mymap = Basemap(projection='robin', lon_0=0, lat_0=0)
    mymap.drawcoastlines(color='black', linewidth=1, zorder=15)
    mymap.drawmeridians(np.arange(0, 360, 30), zorder=20)
    mymap.drawparallels(np.arange(-90, 90, 30), zorder=20)

    mymap.scatter(sx, sy, s=80, marker='*', c='y', lw=1, zorder=30, edgecolors='g')
    # import ipdb; ipdb.set_trace()
    # get the boundaries for the colorbar
    # col_abs = np.mean(ev['crust_corr']) + np.std(ev['crust_corr'])
    col_abs = np.mean(ev['crust_corr'])

    x, y = mymap(ev['sta_lon'], ev['sta_lat'])
    mymap.scatter(x, y, c=ev['crust_corr'], cmap=plt.cm.jet, s=20, lw=0, zorder=10, edgecolors=None,
                  vmin=-col_abs, vmax=col_abs, alpha=0.8)

    cbar = plt.colorbar(orientation='horizontal', shrink=0.5)
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.set_label('time [s]')
    plt.title('Crustal Correction', weight='bold')

    # ====>>>> 5th plot: Attenuation
    plt.subplot(3, 2, 5)
    mymap = Basemap(projection='robin', lon_0=0, lat_0=0)
    mymap.drawcoastlines(color='black', linewidth=1, zorder=15)
    mymap.drawmeridians(np.arange(0, 360, 30), zorder=20)
    mymap.drawparallels(np.arange(-90, 90, 30), zorder=20)

    mymap.scatter(sx, sy, s=80, marker='*', c='y', lw=1, zorder=30, edgecolors='g')

    # get the boundaries for the colorbar
    col_abs = np.mean(ev['attenuation']) + np.std(ev['attenuation'])

    x, y = mymap(ev['sta_lon'], ev['sta_lat'])
    mymap.scatter(x, y, c=ev['attenuation'], cmap=plt.cm.jet, s=20, lw=0, zorder=10, edgecolors=None,
                  vmin=0, vmax=col_abs, alpha=0.8)

    cbar = plt.colorbar(orientation='horizontal', shrink=0.5)
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.set_label('time [s]')
    plt.title('Attenuation', weight='bold')

    # ====>>>> 6th plot: elevation Correction
    plt.subplot(3, 2, 6)
    mymap = Basemap(projection='robin', lon_0=0, lat_0=0)
    mymap.drawcoastlines(color='black', linewidth=1, zorder=15)
    mymap.drawmeridians(np.arange(0, 360, 30), zorder=20)
    mymap.drawparallels(np.arange(-90, 90, 30), zorder=20)

    mymap.scatter(sx, sy, s=80, marker='*', c='y', lw=1, zorder=30, edgecolors='g')

    # get the boundaries for the colorbar



    x, y = mymap(ev['sta_lon'], ev['sta_lat'])
    mymap.scatter(x, y, c=ev['elevation_corr'], cmap=plt.cm.jet, s=20, lw=0, zorder=10, edgecolors=None,
                  vmin=-col_abs, vmax=col_abs, alpha=0.8)

    cbar = plt.colorbar(orientation='horizontal', shrink=0.5)
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.set_label('time [s]')
    plt.title('Elevation Correction', weight='bold')

    plt.tight_layout()

    # ====>>>> saving the plots
    save_path = input_file[0:-len(input_file.split('/')[-1])]
    save_folder = 'saved_corr_plots'
    save_dir = os.path.join(save_path, save_folder)
    if not os.path.isdir(save_dir):
        os.makedirs(os.path.join(save_dir))
    plt.savefig(save_dir + '/' +'common_corr_ievt_%s.png' % ev[0][1], dpi=300)

    plt.close()
    plt.clf()
    print bc.yellow + 'Saved figure in folder: %s' % save_dir + bc.end

    # SEPERATE PLOTS
    # ====>>>> 1st plot: crustal thickness
    plt.ioff()
    plt.figure(figsize=(15, 10))
    plt.suptitle('ievt: %s' % ev[0][1], weight='bold', size=15)

    plt.subplot(1, 2, 1)
    mymap = Basemap(projection='robin', lon_0=0, lat_0=0)
    mymap.drawcoastlines(color='black', linewidth=1, zorder=15)
    mymap.drawmeridians(np.arange(0, 360, 30), zorder=20)
    mymap.drawparallels(np.arange(-90, 90, 30), zorder=20)

    mymap.scatter(sx, sy, s=80, marker='*', c='y', lw=1, zorder=30, edgecolors='g')

    # get the boundaries for the colorbar
    col_abs = np.mean(ev['crust_thick']) + np.std(ev['crust_thick'])

    x, y = mymap(ev['sta_lon'], ev['sta_lat'])
    mymap.scatter(x, y, c=ev['crust_thick'], cmap=plt.cm.jet, s=20, lw=0, zorder=10, edgecolors=None,
                  vmin=-col_abs, vmax=col_abs, alpha=0.8)

    cbar = plt.colorbar(orientation='horizontal', shrink=0.5)
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.set_label('[km]')
    plt.title('Cruatl Thickness', weight='bold')

    # ====>>>> 2nd plot: Elevation
    plt.subplot(1, 2, 2)
    mymap = Basemap(projection='robin', lon_0=0, lat_0=0)
    mymap.drawcoastlines(color='black', linewidth=1, zorder=15)
    mymap.drawmeridians(np.arange(0, 360, 30), zorder=20)
    mymap.drawparallels(np.arange(-90, 90, 30), zorder=20)

    mymap.scatter(sx, sy, s=80, marker='*', c='y', lw=1, zorder=30, edgecolors='g')

    col_abs = np.mean(ev['elevation']) + np.std(ev['elevation'])

    x, y = mymap(ev['sta_lon'], ev['sta_lat'])
    mymap.scatter(x, y, c=ev['elevation'], cmap=plt.cm.jet, s=20, lw=0, zorder=10, edgecolors=None,
                  vmin=-col_abs, vmax=col_abs, alpha=0.8)

    cbar = plt.colorbar(orientation='horizontal', shrink=0.5)
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.set_label('[m]')
    plt.title('Elevation', weight='bold')

    # ====>>>> saving the plots
    save_path = input_file[0:-len(input_file.split('/')[-1])]
    save_folder = 'saved_corr_plots'
    save_dir = os.path.join(save_path, save_folder)
    if not os.path.isdir(save_dir):
        os.makedirs(os.path.join(save_dir))
    plt.savefig(save_dir + '/' + 'crthick_elevation_%s.png' % ev[0][1], dpi=300)

    plt.close()
    plt.clf()
    print bc.yellow + 'Saved figure in folder: %s' % save_dir + bc.end
    # import ipdb;ipdb.set_trace()



'''
OLD PART
# ====>>>> Plotting
plt.ion()
plt.figure(figsize=(15, 10))

plt.subplot(2, 2, 1)
mymap = Basemap(projection='robin', lon_0=0, lat_0=0)
mymap.drawcoastlines(color='black', linewidth=1, zorder=15)
mymap.drawmeridians(np.arange(0, 360, 30), zorder=20)
mymap.drawparallels(np.arange(-90, 90, 30), zorder=20)

x, y = mymap(ccorr_arr[:, 4], ccorr_arr[:, 3])
mymap.scatter(x, y, c=ccorr_arr[:, 6], cmap=plt.cm.jet, s=20, lw=1, zorder=10, edgecolors='k')

cbar = plt.colorbar(orientation='horizontal', shrink=0.5)
cbar.ax.tick_params(labelsize=

)
plt.title('Ellipticity Correction', weight='bold')

plt.subplot(2, 2, 2)
mymap = Basemap(projection='robin', lon_0=0, lat_0=0)
mymap.drawcoastlines(color='black', linewidth=1, zorder=15)
mymap.drawmeridians(np.arange(0, 360, 30), zorder=20)
mymap.drawparallels(np.arange(-90, 90, 30), zorder=20)

x, y = mymap(ccorr_arr[:, 4], ccorr_arr[:, 3])
mymap.scatter(x, y, c=ccorr_arr[:, 7], cmap=plt.cm.jet, s=20, lw=1, zorder=10, edgecolors='k')

cbar = plt.colorbar(orientation='horizontal', shrink=0.5)
cbar.ax.tick_params(labelsize=10)
plt.title('Crustal Correction', weight='bold')

plt.subplot(2, 2, 3)
mymap = Basemap(projection='robin', lon_0=0, lat_0=0)
mymap.drawcoastlines(color='black', linewidth=1, zorder=15)
mymap.drawmeridians(np.arange(0, 360, 30), zorder=20)
mymap.drawparallels(np.arange(-90, 90, 30), zorder=20)

x, y = mymap(ccorr_arr[:, 4], ccorr_arr[:, 3])
mymap.scatter(x, y, c=ccorr_arr[:, 8], cmap=plt.cm.jet, s=20, lw=1, zorder=10, edgecolors='k')

cbar = plt.colorbar(orientation='horizontal', shrink=0.5)
cbar.ax.tick_params(labelsize=10)
plt.title('Elevation Correction', weight='bold')

plt.subplot(2, 2, 4)
mymap = Basemap(projection='robin', lon_0=0, lat_0=0)
mymap.drawcoastlines(color='black', linewidth=1, zorder=15)
mymap.drawmeridians(np.arange(0, 360, 30), zorder=20)
mymap.drawparallels(np.arange(-90, 90, 30), zorder=20)

x, y = mymap(ccorr_arr[:, 4], ccorr_arr[:, 3])
mymap.scatter(x, y, c=ccorr_arr[:, 5], cmap=plt.cm.jet, s=20, lw=1, zorder=10, edgecolors='k')

cbar = plt.colorbar(orientation='horizontal', shrink=0.5)
cbar.ax.tick_params(labelsize=10)
plt.title('Depth [rad]', weight='bold')

plt.suptitle('%s' % ell_ccor_fio.split('/')[-1], weight='bold', size=15)

plt.show()
'''
