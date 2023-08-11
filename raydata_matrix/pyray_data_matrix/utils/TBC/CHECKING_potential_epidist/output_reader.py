import glob
import math as m_math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.basemap import Basemap
import numpy as np
from obspy import UTCDateTime
import os
import subprocess
import shutil
import sys
import time
import pyvtk as pvtk

####################### GENERAL INPUT ##################################
# Transparency of the scatterer plots
global alpha_plt
alpha_plt = 0.01

####################### output_file_reader #############################


def output_file_reader(evinfo, req_band='band01'):
    """
    This function reads one output file and pass all the values
    ATTENTION: there is not filtering at this stage, it reads all...
    """
    add_output = evinfo[0]
    empty_array = np.array([])
    if not os.path.isfile(os.path.join(add_output, 'outfiles',
                                       'ffproc.ampstt.%s' % req_band)):
        print "%s is not found!" % os.path.join(add_output, 'outfiles',
                                                'ffproc.ampstt.%s' % req_band)
        return empty_array
    if not os.path.isfile(os.path.join(add_output, 'outfiles',
                                       'ffproc.receivers')):
        print "%s is not found!" % os.path.join(add_output, 'outfiles',
                                                'ffproc.receivers')
        return empty_array
    rd_output = np.loadtxt(os.path.join(add_output, 'infiles',
                                        'yspec', 'in.matlab.ttime'),
                           dtype='S', comments='#',
                           delimiter=',')
    ### sta_output = np.loadtxt(os.path.join(add_output, 'outfiles',
    ###                                      'ffproc.receivers'),
    ###                         dtype='S', comments='#')
    ### # To avoid any problem in case that there is only one src-rcv available:
    ### if np.size(rd_output) == 22:
    ###     rd_output = np.reshape(rd_output, [1, 22])
    ### if np.size(sta_output) == 8:
    ###     sta_output = np.reshape(sta_output, [1, 8])
    ### # Generate an empty array with size of: rd_output_rows + 9
    ### new_col_cr = np.empty([np.shape(rd_output)[0], 9], dtype=object)
    ### new_col_cr[:, 0] = (sta_output[:, 5].astype(float) -
    ###                     sta_output[:, 6].astype(float))/1000.
    ### new_col_cr[:, 1] = os.path.basename(add_output)
    ### new_col_cr[:, 2] = evinfo[1]
    ### new_col_cr[:, 3] = evinfo[2]
    ### new_col_cr[:, 4] = evinfo[3]
    ### new_col_cr[:, 5] = evinfo[4]
    ### new_col_cr[:, 6] = evinfo[5]
    ### new_col_cr[:, 7] = evinfo[6]
    ### new_col_cr[:, 8] = req_band
    ### # Now append the rd_output_rows X 9 array to the original rd_output
    ### output_sta_evname = np.append(rd_output, new_col_cr, 1)
    return rd_output

####################### event_filter #############################


def event_filter(par_add, selected_events_add, all_events=True,
                 min_dp=-10, max_dp=1000):
    """
    Filters the events in one par_dir based on the required inputs
    """
    if not os.path.isdir(par_add):
        print "%s is not a valid directory!" % par_add
        return False, False

    selected_events = np.loadtxt(fname=selected_events_add, dtype='S',
                                 comments='#', delimiter=',')
    event_adds = []
    for i in range(len(selected_events[:, 1])):
        if all_events != True:
            if not selected_events[:, 1][i] in all_events: continue
        event_adds.append([os.path.join(par_add, selected_events[:, 1][i]),
                           selected_events[:, 0][i]])
    passed_event_adds = []
    for i in range(len(event_adds)):
        add_flag = False
        try:
            fio_source = open(os.path.join(event_adds[i][0], 'outfiles',
                                           'ffproc.source'), 'r')
            f_source = fio_source.readlines()
            ev_year, ev_julianday, ev_hr, ev_min, ev_sec, ev_msec = \
                f_source[1].split()
            evlat, evlon, catalog_depth, inverted_depth = f_source[3].split()
            try:
                mrr, mtt, mpp, mrt, mrp, mtp = f_source[13].split()
            except Exception, e:
                print 'WARNING: inverted moment tensor was not found: %s' \
                      % event_adds[i][0]
                mrr, mtt, mpp, mrt, mrp, mtp = f_source[7].split()
            ev_date = UTCDateTime(year=int(ev_year), julday=int(ev_julianday))
            ev_date_str = '%4s%3s' % (ev_date.year, ev_date.julday)
            ev_date_str = ev_date_str.replace(' ', '0')
            ev_time = '%2s%2s%2s' % (ev_hr, ev_min, ev_sec)
            ev_time = ev_time.replace(' ', '0')
            ev_id = event_adds[i][1]

            # Check for depth
            if min_dp <= float(inverted_depth) < max_dp:
                add_flag = True
            if add_flag:
                passed_event_adds.append([event_adds[i][0], ev_date_str,
                                          ev_time, ev_id, evlat, evlon,
                                          inverted_depth])
        except Exception, e:
            print 'ERROR: %s' % e

    # Sort the events by depth and then pass it
    return sorted(passed_event_adds,
                  key=lambda passed_event_adds: float(passed_event_adds[6]))

####################### array_station_filter #############################


def array_station_filter(passed_array, min_xcorr=-100, max_xcorr=100,
                         min_epi=0., max_epi=360., check_clip=True):
    """
    Filters the stations in an array based on the required inputs
    """
    empty_array = np.array([])

    # --------------- XCORRELATION ---------------
    pass_stas_corr_1 = passed_array[
        passed_array[:, 4].astype(float) >= min_epi]
    if not pass_stas_corr_1.size == 0:
        pass_stas_final = pass_stas_corr_1[
            pass_stas_corr_1[:, 4].astype(float) < max_epi]
        return pass_stas_final
    else:
        return empty_array

    # --------------- EPICENTRAL ---------------
    #passed_stas_epi_1 = pass_stas_final[
    #    pass_stas_final[:, 4].astype(float) >= min_epi]
    #if not passed_stas_epi_1.size == 0:
    #    passed_stas_final = passed_stas_epi_1[
    #        passed_stas_epi_1[:, 4].astype(float) < max_epi]
    #    return passed_stas_final
    #else:
    #    return empty_array

    # --------------- CHECK CLIPS ---------------
    ##if check_clip:
    ##    passed_stas_final = passed_stas_final[
    ##        passed_stas_final[:, 19].astype(float) < 0.1]

    ##if not passed_stas_final.size == 0:
    ##    return passed_stas_final
    ##else:
    ##    return empty_array

####################### array_station_filter_mark #############################


def array_station_filter_mark(all_output_files, min_xcorr=-100, max_xcorr=100,
                              min_epi=0., max_epi=360., check_clip=True):
    """
    Filters the stations in an array based on the required inputs
    """
    # --------------- XCORRELATION ---------------
    all_output_files[:, 10, :][all_output_files[:, 6, :].astype(float) <
                               min_xcorr] = -1
    all_output_files[:, 10, :][all_output_files[:, 6, :].astype(float) >=
                               max_xcorr] = -1

    # --------------- EPICENTRAL ---------------
    all_output_files[:, 10, :][all_output_files[:, 4, :].astype(float) <
                               min_epi] = -1
    all_output_files[:, 10, :][all_output_files[:, 4, :].astype(float) >=
                               max_epi] = -1

    # --------------- CHECK CLIPS ---------------
    if check_clip:
        all_output_files[:, 10, :] \
            [all_output_files[:, 19, :].astype(float) > 0.1] = -1

    return all_output_files

####################### check_selection #############################


def check_selection(filt_array, min_xcorr, max_xcorr, min_epi, max_epi,
                    check_clip=True):
    """
    Check whether the selection procedure worked well
    """
    global alpha_plt

    if len(filt_array) < 100:
        alpha_plt = 1.0
    elif 100 <= len(filt_array) < 10000:
        alpha_plt = 0.2
    elif 10000 <= len(filt_array) < 50000:
        alpha_plt = 0.1
    elif 50000 <= len(filt_array) < 200000:
        alpha_plt = 0.01
    else:
        alpha_plt = 0.002

    plt.ion()
    plt.figure()
    plt.subplot(2, 2, 1)
    plt.scatter(filt_array[:, 2].astype(float),
                filt_array[:, 6].astype(float), c='r', edgecolors='none',
                alpha=alpha_plt)
    plt.xlabel('Station latitude', size='large', weight='bold')
    plt.ylabel('xcorr factor', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.hlines(min_xcorr, np.min(filt_array[:, 2].astype(float)),
               np.max(filt_array[:, 2].astype(float)), 'k',
               linestyles='--')
    plt.hlines(max_xcorr, np.min(filt_array[:, 2].astype(float)),
               np.max(filt_array[:, 2].astype(float)), 'k',
               linestyles='--')
    plt.vlines(50.0, min_xcorr, max_xcorr, 'k', linestyles='-')
    plt.vlines(26.0, min_xcorr, max_xcorr, 'k', linestyles='-')
    plt.vlines(60.0, min_xcorr, max_xcorr, 'k', linestyles='dotted')
    plt.vlines(34.0, min_xcorr, max_xcorr, 'k', linestyles='dotted')

    plt.subplot(2, 2, 2)
    plt.scatter(filt_array[:, 3].astype(float),
                filt_array[:, 6].astype(float), c='r', edgecolors='none',
                alpha=alpha_plt)
    plt.xlabel('Station longitude', size='large', weight='bold')
    plt.ylabel('xcorr factor', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.hlines(min_xcorr, np.min(filt_array[:, 3].astype(float)),
               np.max(filt_array[:, 3].astype(float)), 'k',
               linestyles='--')
    plt.hlines(max_xcorr, np.min(filt_array[:, 3].astype(float)),
               np.max(filt_array[:, 3].astype(float)), 'k',
               linestyles='--')
    plt.vlines(-127.0, min_xcorr, max_xcorr, 'k', linestyles='-')
    plt.vlines(-60.0, min_xcorr, max_xcorr, 'k', linestyles='-')
    plt.vlines(-12.0, min_xcorr, max_xcorr, 'k', linestyles='dotted')
    plt.vlines(41.0, min_xcorr, max_xcorr, 'k', linestyles='dotted')

    plt.subplot(2, 2, 3)
    plt.scatter(filt_array[:, 4].astype(float),
                filt_array[:, 6].astype(float), c='b', edgecolors='none',
                alpha=alpha_plt)
    plt.xlabel('Epicentral distance', size='large', weight='bold')
    plt.ylabel('xcorr factor', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.vlines(min_epi, min_xcorr, max_xcorr, 'k', linestyles='--')
    plt.vlines(max_epi, min_xcorr, max_xcorr, 'k', linestyles='--')

    plt.subplot(2, 2, 4)
    plt.scatter(filt_array[:, 19].astype(float),
                filt_array[:, 6].astype(float), c='r', edgecolors='none',
                alpha=alpha_plt)
    plt.xlabel('Clip value', size='large', weight='bold')
    plt.ylabel('xcorr factor', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')

    # Importing as a MATLAB object for further processing
    import py2mat_mod
    py2mat_mod.py2mat(filt_array[:, 2].astype(float), 'stlat', 'stlat')
    py2mat_mod.py2mat(filt_array[:, 3].astype(float), 'stlon', 'stlon')
    py2mat_mod.py2mat(filt_array[:, 4].astype(float), 'epi_dist',
                      'epi_dist')
    py2mat_mod.py2mat(filt_array[:, 6].astype(float), 'xcorr', 'xcorr')

    """
    ----------------------
    Useful matlab commands
    ----------------------
    load('stlat.mat')
    load('stlon.mat')
    load('epi_dist.mat')
    load('xcorr.mat')

    dat = [epi_dist; xcorr];
    dat = dat';
    %n = log(hist3(dat, [123, 20])); % default is to 10x10 bins
    n = hist3(dat, [123, 20]); % default is to 10x10 bins
    n1 = n';
    n1(size(n,2) + 1, size(n,1) + 1) = 0;
    xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
    yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,2)+1);
    h = pcolor(xb,yb,n1);
    caxis([0, 1500])
    colormap(jet(15));

    shading flat
    """

    plt.tight_layout()
    plt.show()

####################### unique_rows #############################


def unique_rows(data):
    """
    Make a unique array based on the array's row
    """
    data_mid =np.ascontiguousarray(data).view(
        np.dtype((np.void, data.dtype.itemsize * data.shape[1])))
    _, idx = np.unique(data_mid, return_index=True)
    unique_data = data[idx]
    return unique_data

####################### plot_sta_ev_unique #############################


def plot_sta_ev_unique(sta_info_uniq, sta_info_str_arr, ev_info_uniq,
                       ev_info_str_arr, phase, input_file_name_part,
                       req_band, num_color_grp=11):
    """
    Plot all stations and events which have been passed the criteria
    """
    print 'Plotting events and stations...'

    plt.figure()

    plt.subplot2grid((4, 4), (0, 0), colspan=3, rowspan=3)

    m = Basemap(projection='robin', lon_0=0.0, lat_0=0.0, resolution='c')
    m.fillcontinents()
    m.drawparallels(np.arange(-90., 120., 30.))
    m.drawmeridians(np.arange(0., 420., 60.))
    m.drawmapboundary()

    # STATIONS
    col_values = range(num_color_grp+4)
    cm_cmap = plt.get_cmap('Greens_r')
    cNorm_cmap = colors.Normalize(vmin=0, vmax=col_values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm_cmap, cmap=cm_cmap)

    x_ev, y_ev = m(-360, 0)
    for i in range(num_color_grp-1):
        m.scatter(x_ev, y_ev, 100, color=scalarMap.to_rgba(col_values[i]),
                  marker="o", edgecolor="none", zorder=0,
                  label='%s-%s' % (100*i, 100*(i+1)))
    m.scatter(x_ev, y_ev, 100, color=scalarMap.to_rgba(col_values[-5]),
              marker="o", edgecolor="none", zorder=0,
              label='>=%s' % ((num_color_grp-1)*100))

    sta_fio = open(os.path.join(
        os.path.curdir, '%s_staloc_%s_%s.txt'
                        % (phase, input_file_name_part, req_band)), 'w')
    num_unique_sta = []
    for i in range(len(sta_info_uniq)):
        try:
            x, y = m(float(sta_info_uniq[i, 1]), float(sta_info_uniq[i, 0]))
            sta_lat_lon_str = '%s_%s_%s' % (sta_info_uniq[i, 0],
                                            sta_info_uniq[i, 1],
                                            sta_info_uniq[i, 2])
            sta_size = len(sta_info_str_arr[sta_info_str_arr[:] ==
                                            sta_lat_lon_str])
            num_unique_sta.append(sta_size)
            #sta_fio.writelines('%s   %s   %s  %s\n' % (sta_info_uniq[i, 0], sta_info_uniq[i, 1],
            #                                            sta_info_uniq[i, 2], sta_size))
            sta_fio.writelines('%s   %s   %s \n' % (sta_info_uniq[i, 0],
                                                    sta_info_uniq[i, 1],
                                                    sta_info_uniq[i, 2]))
            #sta_size = np.log(sta_size)
            if 0 <= sta_size < 100:
                sta_color = scalarMap.to_rgba(col_values[0])
            elif 100 <= sta_size < 200:
                sta_color = scalarMap.to_rgba(col_values[1])
            elif 200 <= sta_size < 300:
                sta_color = scalarMap.to_rgba(col_values[2])
            elif 300 <= sta_size < 400:
                sta_color = scalarMap.to_rgba(col_values[3])
            elif 400 <= sta_size < 500:
                sta_color = scalarMap.to_rgba(col_values[4])
            elif 500 <= sta_size < 600:
                sta_color = scalarMap.to_rgba(col_values[5])
            elif 600 <= sta_size < 700:
                sta_color = scalarMap.to_rgba(col_values[6])
            elif 700 <= sta_size < 800:
                sta_color = scalarMap.to_rgba(col_values[7])
            elif 800 <= sta_size < 900:
                sta_color = scalarMap.to_rgba(col_values[8])
            elif 900 <= sta_size < 1000:
                sta_color = scalarMap.to_rgba(col_values[9])
            else:
                sta_color = scalarMap.to_rgba(col_values[10])
            m.scatter(x, y, c=sta_color, edgecolor='none', zorder=40,
                      marker='v', s=100)
        except Exception, e:
            print '[exception] in station %s: %s' % (i, e)
    sta_fio.close()

    # EVENTS
    col_values = range(num_color_grp+4)
    cm_cmap = plt.get_cmap('hot')
    cNorm_cmap = colors.Normalize(vmin=0, vmax=col_values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm_cmap, cmap=cm_cmap)

    x_ev, y_ev = m(-360, 0)
    for i in range(num_color_grp-1):
        m.scatter(x_ev, y_ev, 100, color=scalarMap.to_rgba(col_values[i]),
                  marker="o", edgecolor="none", zorder=0,
                  label='%s-%s' % (100*i, 100*(i+1)))
    m.scatter(x_ev, y_ev, 100, color=scalarMap.to_rgba(col_values[-5]),
              marker="o", edgecolor="none", zorder=0,
              label='>=%s' % ((num_color_grp-1)*100))

    plt.legend(bbox_to_anchor=(1.2, 1.0), loc=2, borderaxespad=0., ncol=2,
               mode='expand', title='   Event        Station', fontsize=18)
    plt.gca().get_legend().get_title().set_fontsize(24)
    plt.gca().get_legend().get_title().set_fontweight('bold')

    ev_fio = open(os.path.join(
        os.path.curdir, '%s_evloc_%s_%s.txt'
                        % (phase, input_file_name_part, req_band)), 'w')
    num_unique_ev = []
    for i in range(len(ev_info_uniq)):
        try:
            x, y = m(float(ev_info_uniq[i, 1]), float(ev_info_uniq[i, 0]))
            ev_lat_lon_str = '%s_%s_%s' % (ev_info_uniq[i, 0],
                                           ev_info_uniq[i, 1],
                                           ev_info_uniq[i, 2])
            ev_size = len(ev_info_str_arr[ev_info_str_arr[:] ==
                                          ev_lat_lon_str])
            num_unique_ev.append(ev_size)
            #ev_fio.writelines('%s   %s   %s   %s\n' % (ev_info_uniq[i, 0], ev_info_uniq[i, 1],
            #                                           ev_info_uniq[i, 2], ev_size))
            ev_fio.writelines('%s   %s   %s\n' % (ev_info_uniq[i, 0],
                                                  ev_info_uniq[i, 1],
                                                  ev_info_uniq[i, 2]))
            #ev_size = np.log(ev_size)
            if 0 <= ev_size < 100:
                ev_color = scalarMap.to_rgba(col_values[0])
            elif 100 <= ev_size < 200:
                ev_color = scalarMap.to_rgba(col_values[1])
            elif 200 <= ev_size < 300:
                ev_color = scalarMap.to_rgba(col_values[2])
            elif 300 <= ev_size < 400:
                ev_color = scalarMap.to_rgba(col_values[3])
            elif 400 <= ev_size < 500:
                ev_color = scalarMap.to_rgba(col_values[4])
            elif 500 <= ev_size < 600:
                ev_color = scalarMap.to_rgba(col_values[5])
            elif 600 <= ev_size < 700:
                ev_color = scalarMap.to_rgba(col_values[6])
            elif 700 <= ev_size < 800:
                ev_color = scalarMap.to_rgba(col_values[7])
            elif 800 <= ev_size < 900:
                ev_color = scalarMap.to_rgba(col_values[8])
            elif 900 <= ev_size < 1000:
                ev_color = scalarMap.to_rgba(col_values[9])
            else:
                ev_color = scalarMap.to_rgba(col_values[10])
            m.scatter(x, y, c=ev_color, edgecolor='none', zorder=40, marker='o', s=100)
        except Exception, e:
            print '[exception] in event %s: %s' % (i, e)
    ev_fio.close()

    plt.figure()
    plt.subplot(2, 1, 1)
    num_unique_sta.sort()
    plt.plot(num_unique_sta, 'r', lw=3)
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.xlabel('station number', size='large', weight='bold')
    plt.ylabel('# of repetition', size='large', weight='bold')
    plt.title('Total number of src-rcv pairs: %s' % (sum(num_unique_sta)),
              size='x-large', weight='bold')

    plt.subplot(2, 1, 2)
    num_unique_ev.sort()
    plt.plot(num_unique_ev, 'b', lw=3)
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.xlabel('event number', size='large', weight='bold')
    plt.ylabel('# of repetition', size='large', weight='bold')
    plt.title('Total number of src-rcv pairs: %s' % (sum(num_unique_ev)),
              size='x-large', weight='bold')

    plt.tight_layout()
    plt.show()
    raw_input(' ')

####################### write_staev_locfile #############################


def write_staev_locfile(filt_array, phase, input_file_name_part, req_band):
    """
    Write station and event location files to be used for other calculations
    such as inversion mesh generation
    """
    # Keep lat, lon, ele(depth) for stations (events)
    sta_info = np.array([filt_array[:, 2].astype(float),
                         filt_array[:, 3].astype(float),
                         filt_array[:, 22].astype(float)])
    ev_info = np.array([filt_array[:, 27].astype(float),
                        filt_array[:, 28].astype(float),
                        filt_array[:, 29].astype(float)])
    # Each row relates to one event-station pair
    sta_info = np.transpose(sta_info)
    ev_info = np.transpose(ev_info)

    sta_info_str = []
    for i in range(len(sta_info)):
        sta_info_str.append('%s_%s_%s'
                            % (sta_info[i, 0], sta_info[i, 1], sta_info[i, 2]))

    ev_info_str = []
    for i in range(len(ev_info)):
        ev_info_str.append('%s_%s_%s'
                           % (ev_info[i, 0], ev_info[i, 1], ev_info[i, 2]))

    sta_info_str_arr = np.array(sta_info_str)
    ev_info_str_arr = np.array(ev_info_str)

    sta_info_uniq = unique_rows(sta_info)
    ev_info_uniq = unique_rows(ev_info)

    print '\n========================'
    print '#unique stations: %s' % len(sta_info_uniq)
    print '#unique events  : %s' % len(ev_info_uniq)
    print '========================'

    plot_sta_ev_unique(sta_info_uniq, sta_info_str_arr, ev_info_uniq,
                       ev_info_str_arr, phase, input_file_name_part, req_band)


####################### axi_kernel_receiver_writer ############################


def axi_kernel_receiver_writer(evstas):
    """
    Create AXISEM receiver.dat file for further analysis using AXISEM
    """
    print '\n======>> Create AXISEM receiver.dat file for Kernel calculations'
    if not os.path.isdir(os.path.join(os.path.curdir,
                                      'RESULTS', evstas[0, 23, 0])):
        os.mkdir(os.path.join(os.path.curdir, 'RESULTS', evstas[0, 23, 0]))
    else:
        print 'Directory already exists: %s' % os.path.join(
            os.path.curdir, 'RESULTS', evstas[0, 23, 0])
        return

    band_period = {'band01': 30.0, 'band02': 21.2, 'band03': 15.0,
                   'band04': 10.6, 'band05': 7.5, 'band06': 5.3,
                   'band07': 3.7, 'band08': 2.7}

    line_write = []
    counter = 0
    for j in range(evstas.shape[0]):
        len_freq_avail = len(evstas[j, :, evstas[j, 10, :] != -1])
        if len_freq_avail == 0:
            continue
        third_line = '%s_%s  %s  %s  %s\n' % (evstas[j, 5, 0].split('.')[0],
                                              evstas[j, 5, 0].split('.')[1],
                                              evstas[j, 2, 0], evstas[j, 3, 0],
                                              len_freq_avail)
        line_write.append(third_line)
        for i in range(len(evstas[j, 10, :])):
            if evstas[j, 10, i] == -1:
                continue
            period = band_period[evstas[j, 30, i]]
            forth_line = 'kernel_%s  Gabor_%s  CC  %s  %s\n' \
                         % (str(period).replace('.', '_'), period,
                            evstas[j, 17, i],
                            float(evstas[j, 17, i])+float(evstas[j, 18, i]))
            line_write.append(forth_line)
        counter += 1
    first_line = '%i\n' % counter
    second_line = 'vp  Z\n'
    receiver_fio = open(os.path.join(os.path.curdir, 'RESULTS',
                                     evstas[0, 23, 0], 'receiver.dat'), 'w')
    receiver_fio.writelines(first_line)
    receiver_fio.writelines(second_line)
    for ln in line_write:
        receiver_fio.writelines(ln)
    receiver_fio.close()

    readme_fio = open(
        os.path.join(os.path.curdir, 'RESULTS',
                     evstas[0, 23, 0], 'README.txt'), 'w')
    readme_fio.writelines('number_of_receivers: %i\n' % counter)
    readme_fio.close()

    #plt.ion()
    #plt.figure()
    #plt.subplot(2, 1, 1)
    #plt.plot(filt_array[:, 2].astype(float), t_corr, 'r.')
    #plt.xlabel('Latitude', size='large', weight='bold')
    #plt.ylabel('Common Correction Values', size='large', weight='bold')
    #plt.xticks(size='large', weight='bold')
    #plt.yticks(size='large', weight='bold')
    #plt.subplot(2, 1, 2)
    #plt.plot(filt_array[:, 3].astype(float), t_corr, 'r.')
    #plt.xlabel('Longitude', size='large', weight='bold')
    #plt.ylabel('Common Correction Values', size='large', weight='bold')
    #plt.xticks(size='large', weight='bold')
    #plt.yticks(size='large', weight='bold')
    #plt.show

####################### axi_kernel_receiver_combiner ##########################


def axi_kernel_receiver_combiner(measure_par_add, measures):
    """
    For several type of measurements (P, Pdiff, ...), combine the receiver.dat file
    Usage:
    from output_reader import axi_kernel_receiver_combiner
    axi_kernel_receiver_combiner('./RESULTS', ['P', 'Pdiff'])
    """
    ls_all_events = np.array([], dtype='object')
    for i in range(len(measures)):
        ls_all_events = np.append(ls_all_events,
                                  glob.glob(os.path.join(measure_par_add,
                                                         measures[i],
                                                         '*.*.*.*')))
    for i in range(len(ls_all_events)):
        ls_all_events[i] = os.path.basename(ls_all_events[i])
    ls_all_events_unique = np.unique(ls_all_events)
    for i in range(len(ls_all_events_unique)):
        cont_flag = True
        for j in range(len(measures)):
            if not os.path.isfile(os.path.join(measure_par_add, measures[j],
                                               ls_all_events_unique[i],
                                               'receiver.dat')):
                cont_flag = False
        if not cont_flag:
            continue
        all_staevs = []
        len_staev = 0
        for j in range(len(measures)):
            receiver_tmp_fio = open(os.path.join(measure_par_add, measures[j],
                                                 ls_all_events_unique[i],
                                                 'receiver.dat'), 'r')
            receiver_tmp = receiver_tmp_fio.readlines()
            all_staevs.append(receiver_tmp[2:])

            len_staev_fio = open(os.path.join(measure_par_add, measures[j],
                                              ls_all_events_unique[i],
                                              'README.txt'), 'r')
            len_staev += int(len_staev_fio.readlines()[0].split(':')[1])
        receiver_all_fio = open(
            os.path.join(measure_par_add,
                         'receiver_%s.dat'
                         % ls_all_events_unique[i].replace('.', '_')), 'w')
        receiver_all_fio.writelines('%i\n' % len_staev)
        receiver_all_fio.writelines('vp   Z\n')
        for ln_all in all_staevs:
            receiver_all_fio.writelines(ln_all)
        receiver_all_fio.close()

    print '\n================================='
    print 'WARNING: vp and Z are hard coded!'
    print '================================='

####################### compile_raydata_raymatrix #############################


def compile_raydata_raymatrix():
    """
    Compile both raydata and raymatrix for further usage
    """
    cur_dir = os.path.abspath(os.curdir)
    os.chdir(os.path.join(os.curdir, 'src_raydata_raymatrix', 'raydata_src'))
    os_sys = os.system('./make')
    if not os_sys == 0:
        print "raydata can not be compiled properly"
    os.chdir(cur_dir)

    os.chdir(os.path.join(os.curdir, 'src_raydata_raymatrix', 'raymatrix_src'))
    os_sys = os.system('./make')
    if not os_sys == 0:
        print "raydata can not be compiled properly"
    os.chdir(cur_dir)

####################### raydata_input_generator #############################


def raydata_input_generator(filt_array, input_file_name, twinned, phase,
                            min_xcorr, min_depth, max_depth,
                            min_epi, max_epi, check_clip, req_band):
    """
    Generate input file compatible with raydata input files
    """
    filt_file = 'bpf.omega_m'
    #filt_file = 'gauss_filter_51_15'
    print '\n\n==================='
    print 'Selected filter: %s' % filt_file
    print '==================='
    if not os.path.isfile(os.path.join(os.path.curdir, 'src_raydata_raymatrix',
                                       'files', 'Pdef_%s' % phase)):
        sys.exit('%s could not be found!'
                 % os.path.join(os.path.curdir, 'src_raydata_raymatrix',
                                'files', 'Pdef_%s' % phase))
    phase_def = open(os.path.join(
        os.path.curdir, 'src_raydata_raymatrix',
        'files', 'Pdef_%s' % phase), 'r').readlines()

    if not os.path.isfile(os.path.join(os.path.curdir, 'src_raydata_raymatrix',
                                       'files', filt_file)):
        sys.exit('%s could not be found!'
                 % os.path.join(os.path.curdir, 'src_raydata_raymatrix',
                                'files', filt_file))
    filt_def = open(os.path.join(
        os.path.curdir, 'src_raydata_raymatrix',
        'files', filt_file), 'r').readlines()

    inp_lines = []
    inp_lines.append('%s\n' % input_file_name)
    inp_lines.append('%s\n' % twinned)
    inp_lines.append('# Phase : %s\n' % phase)
    inp_lines.append('# Traveltime data, xcorr_min : %s\n' % min_xcorr)
    inp_lines.append('# depth_min : %s, depth_max : %s\n'
                     % (min_depth, max_depth))
    inp_lines.append('# epi_min : %s, epi_max : %s\n' % (min_epi, max_epi))
    inp_lines.append('# Clip : %s\n' % check_clip)
    for ln in phase_def:
        inp_lines.append(ln)
    for ln in filt_def:
        inp_lines.append(ln)
    for sta in filt_array:
        netid = sta[5].split('.')[0]
        staid = sta[5].split('.')[1]
        chaid = sta[5].split('.')[4]
        first_line = '%s %s %s %s     %s   %s %s    %s    %s    %s    %s  %s   %s 1 0 0\n' \
                     % (sta[24], sta[25], sta[26], sta[1], staid, netid, chaid,
                        sta[27], sta[28], sta[29], sta[2], sta[3], sta[22])
        second_line = '1      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0      0.0\n'
        if phase == 'Pdiff':
            tB_mfi = geocen_dist(evla=float(sta[27]), evlo=float(sta[28]),
                                 stla=float(sta[2]), stlo=float(sta[3]),
                                 evdp=float(sta[29]), bg_model='iasp91')
        else:
            tB_mfi = sta[17]
        if not tB_mfi:
            print "Do not include this one!"
            continue
        third_line = '    %s  %s  %s  %s   %s   %s\n' \
                     % (sta[7], sta[9], sta[6],
                        int(req_band.split('band')[1]), sta[18], tB_mfi)
        inp_lines.append(first_line)
        inp_lines.append(second_line)
        inp_lines.append(third_line)
    input_file_fio = open(os.path.join('RESULTS', input_file_name), 'w')
    input_file_fio.writelines(inp_lines)

####################### geocen_dist #############################


def geocen_dist(evla, evlo, stla, stlo, evdp, bg_model):
    """
    calculate distance between event and station using geocentric latitudes
    """
    fac = 0.993305621334896

    # ----------------- First for station:
    colat = 90.0 - stla
    if colat == 0.:
        colat = 1.0e-5

    # Formula:
    # arg = colat*rpd
    # geocen=pi2-atan(fac*cos(arg)/(max(1.0e-30,sin(arg))))
    geocen_lat = m_math.pi/2. - \
                 m_math.atan(fac*m_math.cos(colat*m_math.pi/180.)/
                             max(1.0e-30, m_math.sin(colat*m_math.pi/180.)))
    stcola_geocen = geocen_lat*180./m_math.pi
    stla_geocen = 90.0 - stcola_geocen

    # ----------------- Second for event:
    colat = 90.0 - evla
    if colat == 0.:
        colat = 1.0e-5

    geocen_lat = m_math.pi/2. - \
                 m_math.atan(fac*m_math.cos(colat*m_math.pi/180.)/
                             max(1.0e-30, m_math.sin(colat*m_math.pi/180.)))
    evcola_geocen = geocen_lat*180./m_math.pi
    evla_geocen = 90.0 - evcola_geocen

    # --------------- TAUP
    #print str(evdp),
    #print str(stla_geocen),
    #print str(stlo),
    #print str(evla_geocen),
    #print str(evlo),
    taup_process = subprocess.Popen(['taup_time', '-mod', bg_model, '-time',
                                     '-h', str(evdp),
                                     '-ph', 'Pdiff',
                                     '-sta', str(stla_geocen), str(stlo),
                                     '-evt', str(evla_geocen), str(evlo)],
                                    stdout=subprocess.PIPE)

    tt_raw = taup_process.communicate()[0]
    try:
        tt = tt_raw.split('\n')[0].split()[-1]
        tt = float(tt)
    except Exception, e:
        print 'Requested phase: Pdiff ... ERROR: %s' % e
        print 'sta: %s %s' % (stla_geocen, stlo) 
        print 'evt: %s %s' % (evla_geocen, evlo) 
        tt = False
    # --------------- END TAUP
    #print tt

    return tt

####################### raydata_input #############################


def raydata_input(bg_model, input_file_name, phase, max_num_arrival,
                  delay_wrt_first_arrival):
    """
    make in.input_file_name for raydata
    """
    in_input_fio = open(os.path.join('RESULTS', 'in.raydata_%s'
                                                % input_file_name), 'w')
    in_input_fio.write('%s\n' % bg_model)
    in_input_fio.write('%s  %s\n' % (max_num_arrival, delay_wrt_first_arrival))
    in_input_fio.write('%s\n' % input_file_name)
    in_input_fio.write('Pdef_%s' % phase)
    in_input_fio.close()

####################### raymatrix_input #############################


def raymatrix_input(vp_vs_Qs, kernel_quad_km, vertex_file, facet_file,
                    input_file_name):
    """
    create in.raymatrix_input_file_name
    """
    in_input_fio = open(
        os.path.join('RESULTS', 'in.raymatrix_%s' % input_file_name), 'w')
    in_input_fio.write('%s %s %s\n' % (vp_vs_Qs[0], vp_vs_Qs[1], vp_vs_Qs[2]))
    in_input_fio.write('%s\n' % kernel_quad_km)
    in_input_fio.write('%s\n' % vertex_file)
    in_input_fio.write('%s\n' % facet_file)
    in_input_fio.write('%s' % input_file_name)
    in_input_fio.close()

    in_input_fio = open(
        os.path.join('RESULTS', 'in.matrixT.%s' % input_file_name), 'w')
    in_input_fio.write('matrixT.%s' % input_file_name)
    in_input_fio.close()

####################### prepare_dir #############################


def prepare_dir(input_file_name):
    """
    prepare directory for one run of raydata and raymatrix
    """
    cur_dir = os.path.abspath(os.curdir)
    if os.path.isdir(os.path.join(
            os.path.curdir, 'RESULTS', '%s_dir' % input_file_name)):
        sys.exit("Directory already exists!")
    os.mkdir(os.path.join(
        os.path.curdir, 'RESULTS', '%s_dir' % input_file_name))
    shutil.move(os.path.join('RESULTS', input_file_name),
                os.path.join(os.path.curdir,
                             'RESULTS', '%s_dir' % input_file_name))
    shutil.move(os.path.join('RESULTS', 'in.raydata_%s' % input_file_name),
                os.path.join(os.path.curdir, 'RESULTS',
                             '%s_dir' % input_file_name))
    shutil.move(os.path.join('RESULTS', 'in.raymatrix_%s' % input_file_name),
                os.path.join(os.path.curdir,
                             'RESULTS', '%s_dir' % input_file_name))
    shutil.move(os.path.join('RESULTS', 'in.matrixT.%s' % input_file_name),
                os.path.join(os.path.curdir,
                             'RESULTS', '%s_dir' % input_file_name))

    files_glob = glob.glob(os.path.join(os.path.curdir,
                                        'src_raydata_raymatrix', 'files', '*'))
    for fi in files_glob:
        shutil.copy(fi, os.path.join('RESULTS', '%s_dir' % input_file_name))


    shutil.copy(os.path.join(os.curdir, 'src_raydata_raymatrix',
                             'raydata_src', 'raydata'),
                os.path.join(cur_dir, 'RESULTS', '%s_dir' % input_file_name))
    shutil.copy(os.path.join(os.curdir, 'src_raydata_raymatrix',
                             'raymatrix_src', 'raymatrix'),
                os.path.join(cur_dir, 'RESULTS', '%s_dir' % input_file_name))
    shutil.copy(os.path.join(os.curdir, 'src_raydata_raymatrix',
                             'raymatrix_src', 'mat2asc'),
                os.path.join(cur_dir, 'RESULTS', '%s_dir' % input_file_name))

####################### run_raydata_raymatrix #############################


def run_raydata_raymatrix(input_file_name, raydata=True, raymatrix=True):
    """
    run both raydata and raymatrix in the directory
    """
    cur_dir = os.path.abspath(os.curdir)
    if raydata:
        print '\n======>> run raydata at ./RESULTS/%s_dir' % input_file_name
        os.chdir(os.path.join(cur_dir, 'RESULTS', '%s_dir' % input_file_name))
        os_sys = os.system('./raydata < in.raydata_%s' % input_file_name)
        if not os_sys == 0:
            print 'raydata was not executed correctly!'

    if raymatrix:
        print '\n======>> run raymatrix at ./RESULTS/%s_dir' % input_file_name
        os.chdir(os.path.join(cur_dir, 'RESULTS', '%s_dir' % input_file_name))
        os_sys = os.system('./raymatrix < in.raymatrix_%s' % input_file_name)
        if not os_sys == 0:
            print 'raymatrix was not executed correctly!'

    os.chdir(cur_dir)

####################### parallel_raydata_raymatrix ############################


def parallel_raydata_raymatrix(filt_array, input_file_name, twinned, phase,
                               min_xcorr, min_depth, max_depth,
                               min_epi, max_epi, check_clip, bg_model,
                               vp_vs_Qs, kernel_quad_km, vertex_file,
                               facet_file, max_num_arrival,
                               delay_wrt_first_arrival,
                               run_raydata, run_raymatrix,
                               req_band):
    """
    To run raydata and raymatrix in parallel
    """
    raydata_input_generator(filt_array=filt_array,
                            input_file_name=input_file_name,
                            twinned=twinned,
                            phase=phase,
                            min_xcorr=min_xcorr,
                            min_depth=min_depth,
                            max_depth=max_depth,
                            min_epi=min_epi,
                            max_epi=max_epi,
                            check_clip=check_clip,
                            req_band=req_band)
    raydata_input(bg_model=bg_model,
                  input_file_name=input_file_name,
                  phase=phase,
                  max_num_arrival=max_num_arrival,
                  delay_wrt_first_arrival=delay_wrt_first_arrival)
    raymatrix_input(vp_vs_Qs=vp_vs_Qs,
                    kernel_quad_km=kernel_quad_km,
                    vertex_file=vertex_file,
                    facet_file=facet_file,
                    input_file_name=input_file_name)

    print '\n======>> prepare output directory at: ./RESULTS/%s_dir' \
          % input_file_name
    prepare_dir(input_file_name=input_file_name)
    run_raydata_raymatrix(input_file_name=input_file_name,
                          raydata=run_raydata,
                          raymatrix=run_raymatrix)

####################### raydata_ccorr_reader #############################


def raydata_ccorr_reader(filt_array, input_file_name, corr_io_list):
    """
    Read common correction results and rewrite the values in filt_array
    """
    ccorr_arr = np.loadtxt(os.path.join(os.path.curdir,
                                        'RESULTS', '%s_dir' % input_file_name,
                                        'ell_ccor.%s' % input_file_name),
                           dtype='S', comments='#', delimiter=',')
    kd_avail = np.unique(ccorr_arr[:, 0].astype(int))
    t_corr = np.zeros(
        len(ccorr_arr[ccorr_arr[:, 0].astype(int) == kd_avail[0]]))
    t_ellip = np.zeros(
        len(ccorr_arr[ccorr_arr[:, 0].astype(int) == kd_avail[0]]))
    t_cc = np.zeros(
        len(ccorr_arr[ccorr_arr[:, 0].astype(int) == kd_avail[0]]))
    t_elev = np.zeros(
        len(ccorr_arr[ccorr_arr[:, 0].astype(int) == kd_avail[0]]))

    for i in range(len(kd_avail)):
        tar_ccorr_arr = ccorr_arr[
            ccorr_arr[:, 0].astype(int) == kd_avail[i]]
        if len(tar_ccorr_arr)/len(t_corr) != 1:
            if len(tar_ccorr_arr) % len(t_corr) != 0:
                sys.exit('Number of ray segments is not set correctly for '
                         'all stations!')
            divi = len(tar_ccorr_arr)/len(t_corr)
        else:
            divi = False
        if not divi:
            if kd_avail[i] == 1:
                if corr_io_list[0] == 1:
                    t_ellip += tar_ccorr_arr[:, 6].astype(float)
            if corr_io_list[1] == 1:
                t_cc += tar_ccorr_arr[:, 7].astype(float)
            if corr_io_list[2] == 1:
                t_elev += tar_ccorr_arr[:, 8].astype(float)
        else:
            if kd_avail[i] == 1:
                print 'This case should not happen (common corrections!)'
                if corr_io_list[0] == 1:
                    e_corrector = np.array([])
                    for j in range(0, len(tar_ccorr_arr), divi):
                        e_corrector = np.append(
                            e_corrector,
                            np.sum(tar_ccorr_arr[j:j+divi, 6].astype(float)))
                    t_ellip += e_corrector
            if corr_io_list[1] == 1:
                cc_corrector = np.array([])
                for j in range(0, len(tar_ccorr_arr), divi):
                    cc_corrector = np.append(
                        cc_corrector,
                        np.sum(tar_ccorr_arr[j:j+divi, 7].astype(float)))
                t_cc += cc_corrector
            if corr_io_list[2] == 1:
                el_corrector = np.array([])
                for j in range(0, len(tar_ccorr_arr), divi):
                    el_corrector = np.append(
                        el_corrector,
                        np.sum(tar_ccorr_arr[j:j+divi, 8].astype(float)))
                t_elev += el_corrector

    t_corr = t_ellip + t_cc + t_elev
    plt.ion()
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.scatter(filt_array[:, 2].astype(float), t_corr,
                c='r', edgecolors='none', alpha=alpha_plt)
    plt.xlabel('Latitude', size='large', weight='bold')
    plt.ylabel('Common Correction Values', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.subplot(2, 1, 2)
    plt.scatter(filt_array[:, 3].astype(float), t_corr,
                c='b', edgecolors='none', alpha=alpha_plt)
    plt.xlabel('Longitude', size='large', weight='bold')
    plt.ylabel('Common Correction Values', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.show()

    plt.ion()
    plt.figure()
    plt.subplot(2, 3, 1)
    plt.scatter(filt_array[:, 2].astype(float), t_ellip, c='r',
                edgecolors='none', alpha=alpha_plt)
    plt.xlabel('Latitude', size='large', weight='bold')
    plt.ylabel('Ellipticity Correction', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.subplot(2, 3, 2)
    plt.scatter(filt_array[:, 2].astype(float), t_cc, c='r',
                edgecolors='none', alpha=alpha_plt)
    plt.xlabel('Latitude', size='large', weight='bold')
    plt.ylabel('Crustal Correction', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.subplot(2, 3, 3)
    plt.scatter(filt_array[:, 2].astype(float), t_elev, c='r',
                edgecolors='none', alpha=alpha_plt)
    plt.xlabel('Latitude', size='large', weight='bold')
    plt.ylabel('Topography Correction', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')

    plt.subplot(2, 3, 4)
    plt.scatter(filt_array[:, 3].astype(float), t_ellip, c='b',
                edgecolors='none', alpha=alpha_plt)
    plt.xlabel('Longitude', size='large', weight='bold')
    plt.ylabel('Ellipticity Correction', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.subplot(2, 3, 5)
    plt.scatter(filt_array[:, 3].astype(float), t_cc, c='b',
                edgecolors='b', alpha=alpha_plt)
    plt.xlabel('Longitude', size='large', weight='bold')
    plt.ylabel('Crustal Correction', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.subplot(2, 3, 6)
    plt.scatter(filt_array[:, 3].astype(float), t_elev, c='b',
                edgecolors='b', alpha=alpha_plt)
    plt.xlabel('Longitude', size='large', weight='bold')
    plt.ylabel('Topography Correction', size='large', weight='bold')
    plt.xticks(size='large', weight='bold')
    plt.yticks(size='large', weight='bold')
    plt.show()

    filt_array[:, 8] = filt_array[:, 8].astype(float) - t_corr
    return filt_array

####################### raydata_ccorr_writer #############################


def raydata_ccorr_writer(filt_array_corr, events_dir):
    """
    Write similar files as ffproc.ampstt.* but with common correction applied
    """
    req_dirs = np.unique(filt_array_corr[:, 23])
    for r_dir in req_dirs:
        if not os.path.isdir(os.path.join(os.path.curdir, 'RESULTS', r_dir)):
            os.mkdir(os.path.join(os.path.curdir, 'RESULTS', r_dir))
        if not os.path.isdir(os.path.join(os.path.curdir, 'RESULTS', r_dir,
                                          'outfiles')):
            os.mkdir(os.path.join(os.path.curdir, 'RESULTS', r_dir,
                                  'outfiles'))
        filt_array_this_dir = filt_array_corr[filt_array_corr[:, 23] == r_dir]
        np.savetxt(os.path.join(os.path.curdir, 'RESULTS', r_dir, 'outfiles',
                                'ffproc.ampstt.%s'
                                % filt_array_this_dir[0, 30]),
                   filt_array_this_dir[:, 0:22],
                   fmt='%s', delimiter='     ')
        shutil.copy(os.path.join(events_dir, r_dir,
                                 'outfiles', 'ffproc.source'),
                    os.path.join(os.curdir, 'RESULTS', r_dir, 'outfiles'))
        shutil.copy(os.path.join(events_dir, r_dir,
                                 'outfiles', 'ffproc.receivers'),
                    os.path.join(os.curdir, 'RESULTS', r_dir, 'outfiles'))
        shutil.copy(os.path.join(events_dir, r_dir,
                                 'outfiles', 'ampinv.source'),
                    os.path.join(os.curdir, 'RESULTS', r_dir, 'outfiles'))

####################### check_par_jobs #############################


def check_par_jobs(jobs, sleep_time=1):
    """
    check whether all the parallel jobs are finished or not
    """
    pp_flag = True
    while pp_flag:
        for proc in jobs:
            if proc.is_alive():
                print '.',
                sys.stdout.flush()
                time.sleep(sleep_time)
                pp_flag = True
                break
            else:
                pp_flag = False
    if not pp_flag:
        print '\n\nAll %s processes are finished...\n' % len(jobs)

####################### mat2asc_run #############################


def mat2asc_run(input_file_name):
    """
    run mat2asc in each directory
    """
    cur_dir = os.path.abspath(os.curdir)
    print '\n======>> run mat2asc at ./RESULTS/%s_dir' % input_file_name
    os.chdir(os.path.join(cur_dir, 'RESULTS', '%s_dir' % input_file_name))
    os_sys = os.system('./mat2asc < in.matrixT.%s' % input_file_name)
    if not os_sys == 0:
        print 'mat2asc was not executed correctly!'
    os.chdir(cur_dir)

####################### vtk_generator #############################


def vtk_generator(input_file_name_part, req_band, vertex_file, facet_file,
                  parallel_exec, len_dirs, abs_vtk):
    """
    VTK file generator out of all the results for one complete run
    INFO:
    - facet file is indexed from 0
    - BUT! ascii.* files are indexed from 1, so we need to subtract them (-1)
    - When describing the cells in terms of point indices, the points must be
    indexed starting at 0.
    """
    # Only plot one kernel:
    single_kernel = False

    print '\n======>> Creating VTK file'
    input_file_name = '%s_%s_%s' % (input_file_name_part, req_band, 1)
    direname = os.path.join(os.path.curdir, 'RESULTS', '%s_dir'
                                                       % input_file_name)

    print '------> Load Vertex file'
    mesh_points = np.loadtxt(os.path.join(direname, vertex_file),
                             skiprows=2, comments='#')
    print '------> Load Facet file'
    mesh_facets = np.loadtxt(os.path.join(direname, facet_file),
                             dtype=int, skiprows=1, comments='#')

    mat_val_all = [0.0]*len(mesh_points)
    for nj in range(len_dirs):
        input_file_name = '%s_%s_%s' % (input_file_name_part, req_band, nj+1)
        direname = os.path.join(os.path.curdir,
                                'RESULTS', '%s_dir' % input_file_name)
        print '------> create VTK file at ./RESULTS/%s_dir' % input_file_name
        ascii_file = 'ascii.matrixT.%s' % input_file_name
        fmatrix = open(os.path.join(direname, ascii_file), 'r')
        fmatrix_r = fmatrix.readlines()
        mat_indx = []
        mat_val = []
        counter = 0
        for j in range(5, len(fmatrix_r)):
            if counter == 3:
                counter = 0
                continue
            if counter == 0:
                # mat_indx: matrix index
                mat_indx_tmp = fmatrix_r[j].split()
                for i in range(len(mat_indx_tmp)):
                    # IF indexing in ASCII file starts from 0,
                    # we do not need -1.
                    # Otherwise we need it!
                    mat_indx.append(int(mat_indx_tmp[i])-1)
            if counter == 1:
                # mat_val: matrix value
                mat_val_tmp = fmatrix_r[j].split()
                for i in range(len(mat_val_tmp)):
                    if abs_vtk:
                        mat_val.append(abs(float(mat_val_tmp[i])))
                    else:
                        mat_val.append(float(mat_val_tmp[i]))
            if single_kernel:
                if j == 7:
                    break
            counter += 1

        for i in range(len(mat_indx)):
            mat_val_all[mat_indx[i]] += mat_val[i]

    print '\n\n=================='
    print 'Sum over all elems: %s' % sum(mat_val_all)

    vtk = pvtk.VtkData(pvtk.UnstructuredGrid(mesh_points,
                                             tetra=mesh_facets),
                       pvtk.PointData(pvtk.Scalars(mat_val_all,
                                                   name='kernel_value')),
                       'Inversion Grid')
    print '=================='
    print "WARNING: INDEXING!"
    print '=================='
    vtk.tofile(os.path.join(
        os.path.curdir, '%s_%s.vtk'
                        % (input_file_name.split('_')[0], req_band)))

####################### vtk_generator_all #############################


def vtk_generator_all(direname, vertex_file, facet_file):
    """
    VTK file generator out of all the results of all runs
    ATTENTION: ascii.matrix.* should be moved to one dir (direname)
    INFO:
    - facet file is indexed from 0
    - BUT! ascii.* files are indexed from 1, so we need to subtract them (-1)
    - When describing the cells in terms of point indices, the points must be
    indexed starting at 0.
    """
    print '\n======>> Creating VTK file'
    ascii_files = glob.glob(os.path.join(direname, 'ascii.matrixT.*'))

    print '------> Load Vertex file'
    mesh_points = np.loadtxt(os.path.join(direname, vertex_file),
                             skiprows=2, comments='#')
    print '------> Load Facet file'
    mesh_facets = np.loadtxt(os.path.join(direname, facet_file),
                             dtype=int, skiprows=1, comments='#')

    mat_val_all = [0.0]*len(mesh_points)
    for nj in range(len(ascii_files)):
        print '\n------> create VTK file %s' % ascii_files[nj]
        fmatrix = open(os.path.join(ascii_files[nj]), 'r')
        fmatrix_r = fmatrix.readlines()
        mat_indx = []
        mat_val = []
        counter = 0
        for j in range(5, len(fmatrix_r)):
            if counter == 3:
                counter = 0
                continue
            if counter == 0:
                # mat_indx: matrix index
                mat_indx_tmp = fmatrix_r[j].split()
                for i in range(len(mat_indx_tmp)):
                    # IF indexing in ASCII file starts from 0,
                    # we do not need -1.
                    # Otherwise we need it!
                    mat_indx.append(int(mat_indx_tmp[i])-1)
            if counter == 1:
                # mat_val: matrix value
                mat_val_tmp = fmatrix_r[j].split()
                for i in range(len(mat_val_tmp)):
                    mat_val.append(abs(float(mat_val_tmp[i])))
            counter += 1

        for i in range(len(mat_indx)):
            mat_val_all[mat_indx[i]] += mat_val[i]

    vtk = pvtk.VtkData(pvtk.UnstructuredGrid(mesh_points,
                                             tetra=mesh_facets),
                       pvtk.PointData(pvtk.Scalars(mat_val_all,
                                                   name='kernel_value')),
                       'Inversion Grid')
    print '\n\n=================='
    print "WARNING: INDEXING!"
    print '=================='
    vtk.tofile(os.path.join(direname, 'global_all.vtk'))

####################### vtk_val_azi #############################


def vtk_val_azi(direname, vertex_file, facet_file):
    """
    VTK file generator out of all the results of all runs
    This version of VTK maker generates a .vtk file that 
    contains the directions (based on azimuth).
    It can be used to see whether a cell is well illuminated
    from all directions.

    ATTENTION: ascii.matrix.* should be moved to one dir (direname)
    INFO:
    - facet file is indexed from 0
    - BUT! ascii.* files are indexed from 1, so we need to subtract them (-1)
    - When describing the cells in terms of point indices, the points must be
    indexed starting at 0.

    Example:
    from output_reader import vtk_val_azi
    direname = './RESULTS/combined'
    vertex_file = 'vertices.ppdiff_staev_200'
    facet_file = 'facets.ppdiff_staev_200'
    vtk_val_azi(direname, vertex_file, facet_file)
    """
    print '\n======>> Creating VTK file'
    ascii_files = glob.glob(os.path.join(direname, 'ascii.matrixT.*'))

    print '------> Load Vertex file'
    mesh_points = np.loadtxt(os.path.join(direname, vertex_file),
                             skiprows=2, comments='#')
    print '------> Load Facet file'
    mesh_facets = np.loadtxt(os.path.join(direname, facet_file),
                             dtype=int, skiprows=1, comments='#')
    mesh_points_2 = mesh_points**2
    rad = np.sqrt(mesh_points_2[:,0] + mesh_points_2[:,1] + mesh_points_2[:,2])

    mat_val_all = np.array([0.0]*len(mesh_points))
    mat_azi_qual = np.zeros([len(mesh_points), 4])
    mat_azi_all = {}
    for nj in range(len(ascii_files)):
        print '\n------> create VTK file %s' % ascii_files[nj]
        fmatrix = open(os.path.join(ascii_files[nj]), 'r')
        fmatrix_r = fmatrix.readlines()
        mat_indx = []
        mat_val = []
        mat_azi = []
        counter = 0
        print 'Length of fmatrix_r: %s' % len(fmatrix_r)
        for j in range(3, len(fmatrix_r)):
            if counter == 0:
                azi = float(fmatrix_r[j].split()[17])*180./np.pi
                #print azi
            if counter == 2:
                # mat_indx: matrix index
                mat_indx_tmp = fmatrix_r[j].split()
                for i in range(len(mat_indx_tmp)):
                    # IF indexing in ASCII file starts from 0,
                    # we do not need -1.
                    # Otherwise we need it!
                    mat_indx.append(int(mat_indx_tmp[i])-1)
            if counter == 3:
                # mat_val: matrix value
                mat_val_tmp = fmatrix_r[j].split()
                for i in range(len(mat_val_tmp)):
                    mat_val.append(abs(float(mat_val_tmp[i])))
                    mat_azi.append(azi)
                counter = -1
            counter += 1

        mat_indx_arr = np.array(mat_indx)
        mat_azi_arr = np.array(mat_azi)

        mat_1 = np.array([0.0]*len(mesh_points))
        mat_1[mat_indx_arr[(0<=mat_azi_arr) & (mat_azi_arr<45)]] = \
                np.bincount(mat_indx_arr[(0<=mat_azi_arr) & (mat_azi_arr<45)])
        mat_azi_qual[:, 0] += mat_1


        mat_1 = np.array([0.0]*len(mesh_points))
        mat_1[mat_indx_arr[(45<=mat_azi_arr) & (mat_azi_arr<90)]] = \
                np.bincount(mat_indx_arr[(45<=mat_azi_arr) & (mat_azi_arr<90)])
        mat_azi_qual[:, 1] += mat_1


        mat_1 = np.array([0.0]*len(mesh_points))
        mat_1[mat_indx_arr[(90<=mat_azi_arr) & (mat_azi_arr<135)]] = \
                np.bincount(mat_indx_arr[(90<=mat_azi_arr) & (mat_azi_arr<135)])
        mat_azi_qual[:, 2] += mat_1


        mat_1 = np.array([0.0]*len(mesh_points))
        mat_1[mat_indx_arr[(135<=mat_azi_arr) & (mat_azi_arr<180)]] = \
                np.bincount(mat_indx_arr[(135<=mat_azi_arr) & (mat_azi_arr<180)])
        mat_azi_qual[:, 3] += mat_1

        mat_1 = np.array([0.0]*len(mesh_points))
        mat_1[mat_indx_arr[(180<=mat_azi_arr) & (mat_azi_arr<225)]] = \
                np.bincount(mat_indx_arr[(180<=mat_azi_arr) & (mat_azi_arr<225)])
        mat_azi_qual[:, 0] += mat_1


        mat_1 = np.array([0.0]*len(mesh_points))
        mat_1[mat_indx_arr[(225<=mat_azi_arr) & (mat_azi_arr<270)]] = \
                np.bincount(mat_indx_arr[(225<=mat_azi_arr) & (mat_azi_arr<270)])
        mat_azi_qual[:, 1] += mat_1


        mat_1 = np.array([0.0]*len(mesh_points))
        mat_1[mat_indx_arr[(270<=mat_azi_arr) & (mat_azi_arr<315)]] = \
                np.bincount(mat_indx_arr[(270<=mat_azi_arr) & (mat_azi_arr<315)])
        mat_azi_qual[:, 2] += mat_1


        mat_1 = np.array([0.0]*len(mesh_points))
        mat_1[mat_indx_arr[(315<=mat_azi_arr) & (mat_azi_arr<360)]] = \
                np.bincount(mat_indx_arr[(315<=mat_azi_arr) & (mat_azi_arr<360)])
        mat_azi_qual[:, 3] += mat_1
        plt.ion()
        plt.figure()
        plt.subplot(2,2,2)
        plt.plot(mat_azi_qual[:,0])
        plt.subplot(2,2,4)
        plt.plot(mat_azi_qual[:,1])
        plt.subplot(2,2,3)
        plt.plot(mat_azi_qual[:,2])
        plt.subplot(2,2,1)
        plt.plot(mat_azi_qual[:,3])
        plt.show()


        #print 'Length of mat_indx: %s' % len(mat_indx)
        mat_val_all[mat_indx] += np.array(mat_val)
        #for i in range(len(mat_indx)):
        #    if i%1000 == 0:
        #        sys.stdout.write('\r')
        #        sys.stdout.write("[%-100s] %d%% ---" % ('='*int(100.*(i+1)/len(mat_indx)),
        #                                                100.*(i+1)/len(mat_indx)))
        #        sys.stdout.flush()

        #    if i%2000 == 0:
        #        sys.stdout.write('\r')
        #        sys.stdout.write("[%-100s] %d%% -|-" % ('='*int(100.*(i+1)/len(mat_indx)),
        #                                                100.*(i+1)/len(mat_indx)))
        #        sys.stdout.flush()
        #    import ipdb; ipdb.set_trace()

        #    if not str(mat_indx[i]) in mat_azi_all.keys():
        #        mat_azi_all[str(mat_indx[i])] = [mat_azi[i]]
        #    else:
        #        mat_azi_all[str(mat_indx[i])].append(mat_azi[i])

    ##mat_azi_qual = [[0.0, 0.0, 0.0, 0.0]]*len(mesh_points)
    #mat_azi_qual = np.zeros([len(mesh_points), 4])
    #for ep in mat_azi_all.keys():
    #    qc_node = np.histogram(mat_azi_all[ep], bins=range(0, 405, 45))
    #    #qc_node = np.histogram(mat_azi_all[ep], bins=[0, 90, 180, 270, 360])
    #    qc_node_combined = np.array([qc_node[0][0] + qc_node[0][4],
    #                                 qc_node[0][1] + qc_node[0][5],
    #                                 qc_node[0][2] + qc_node[0][6],
    #                                 qc_node[0][3] + qc_node[0][7],
    #                                 ])
    #    #if len(np.nonzero(qc_node_combined)[0])
    #    mat_azi_qual[int(ep)] += qc_node_combined

    mat_azi_qual_rank = np.zeros(len(mesh_points))
    for i in range(len(mat_azi_qual)):
        non_zero = np.nonzero(mat_azi_qual[i])[0]
        if len(non_zero) == 4:
            if min(non_zero) >= 10:
                mat_azi_qual_rank[i] = 100
            else:
                mat_azi_qual_rank[i] = 80
        elif len(non_zero) == 3:
            if min(non_zero) >= 10:
                mat_azi_qual_rank[i] = 90
            else:
                mat_azi_qual_rank[i] = 70
        elif len(non_zero) == 2:
            if (non_zero[1] - non_zero[0]) == 2:
                if min(np.nonzero(mat_azi_qual[i])[0]) >= 10:
                    mat_azi_qual_rank[i] = 80
                else:
                    mat_azi_qual_rank[i] = 60
            else:
                if min(np.nonzero(mat_azi_qual[i])[0]) >= 10:
                    mat_azi_qual_rank[i] = 40
                else:
                    mat_azi_qual_rank[i] = 20
        elif len(np.nonzero(mat_azi_qual[i])[0]) == 1:
            if min(np.nonzero(mat_azi_qual[i])[0]) >= 10:
                mat_azi_qual_rank[i] = 20
            else:
                mat_azi_qual_rank[i] = 10
        else:
            mat_azi_qual_rank[i] = 0

    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(mat_val_all)
    plt.subplot(2,1,2)
    plt.plot(mat_azi_qual_rank)
    plt.show()
    
    plt.figure()
    plt.subplot(1,1,1)
    plt.plot(mat_val_all[rad>=3482]*mat_azi_qual_rank[rad>=3482])
    plt.show()

    vtk = pvtk.VtkData(pvtk.UnstructuredGrid(mesh_points,
                                             tetra=mesh_facets),
                       pvtk.PointData(pvtk.Scalars(mat_val_all,
                                                   name='kernel_value')),
                       'Inversion Grid')
    vtk.tofile(os.path.join(direname, 'absolute_value_all.vtk'))

    vtk_azi = pvtk.VtkData(pvtk.UnstructuredGrid(mesh_points,
                                                 tetra=mesh_facets),
                           pvtk.PointData(pvtk.Scalars(mat_azi_qual_rank,
                                                       name='Quality_value')),
                           'Quality check')
    vtk_azi.tofile(os.path.join(direname, 'azi_quality_all.vtk'))
    print '\n\n=================='
    print "WARNING: INDEXING!"
    print '=================='
    raw_input('Press Enter to quit...')

####################### make_assemble_dir #############################


def make_assemble_dir(master_add, target_add, vertex_file, facet_file,
                      dlnvp, dlnvs, dlnQs,
                      hypocen_corr,
                      sta_corr_tp, sta_corr_ap,
                      orig_corr, ev_corr_ap,
                      sta_corr_ts, sta_corr_as,
                      corr_switches, out_ident):
    """
    Create assemblematrix directory to run assemblematrix.f
    :param master_add: copy matrices files from this directory
    :param target_add: where everything will be stored
    :param vertex_file:
    :param facet_file:
    :param dlnvp:
    :param dlnvs:
    :param dlnQs:
    :param hypocen_corr:
    :param sta_corr_tp:
    :param sta_corr_ap:
    :param orig_corr:
    :param ev_corr_ap:
    :param sta_corr_ts:
    :param sta_corr_as:
    :param corr_switches:
    :param out_ident:
    :return:
    """
    print "copying matrixT.* files...",
    mat_files_glob = glob.glob(os.path.join(master_add, '*', 'matrixT.*'))
    mat_files_glob.sort()
    for fi in mat_files_glob:
        shutil.copy(fi, os.path.join(target_add))
    print "DONE"

    print "copying raymatrix.info.T.* files...",
    rayT_files_glob = glob.glob(os.path.join(master_add, '*',
                                             'raymatrix.info.T.*'))
    rayT_files_glob.sort()
    for fi in rayT_files_glob:
        shutil.copy(fi, os.path.join(target_add))
    print "DONE"

    print "copying vertex and facet files...",
    shutil.copy(os.path.join(os.curdir, 'src_raydata_raymatrix', 'files',
                             vertex_file),
                os.path.join(target_add))
    shutil.copy(os.path.join(os.curdir, 'src_raydata_raymatrix', 'files',
                             facet_file),
                os.path.join(target_add))
    print "DONE"

    print "creating the input files...",
    in_assemble_fio = open(os.path.join(target_add, 'in.am'), 'w')
    in_assemble_fio.writelines('%s\n' % vertex_file)
    in_assemble_fio.writelines('%s\n' % facet_file)
    in_assemble_fio.writelines('%s %s\n' % (dlnvp[0], dlnvp[1]))
    in_assemble_fio.writelines('%s %s\n' % (dlnvs[0], dlnvs[1]))
    in_assemble_fio.writelines('%s %s\n' % (dlnQs[0], dlnQs[1]))
    in_assemble_fio.writelines('%s %s\n' % (hypocen_corr[0], hypocen_corr[1]))
    in_assemble_fio.writelines('%s %s\n' % (sta_corr_tp[0], sta_corr_tp[1]))
    in_assemble_fio.writelines('%s %s\n' % (sta_corr_ap[0], sta_corr_ap[1]))
    in_assemble_fio.writelines('%s %s\n' % (orig_corr[0], orig_corr[1]))
    in_assemble_fio.writelines('%s %s\n' % (ev_corr_ap[0], ev_corr_ap[1]))
    in_assemble_fio.writelines('%s %s\n' % (sta_corr_ts[0], sta_corr_ts[1]))
    in_assemble_fio.writelines('%s %s\n' % (sta_corr_as[0], sta_corr_as[1]))
    in_assemble_fio.writelines('%s %s %s %s\n' % (corr_switches[0],
                                                  corr_switches[1],
                                                  corr_switches[2],
                                                  corr_switches[3]))
    in_assemble_fio.writelines('%s\n' % out_ident)
    for fi in mat_files_glob:
        in_assemble_fio.writelines('%s\n' % os.path.basename(fi))
    in_assemble_fio.writelines('stop')
    in_assemble_fio.close()
    print "DONE"

####################### compile_assemblematrix #############################


def compile_assemblematrix(target_add):
    """
    compile assemblematrix.f
    :return:
    """

    cur_dir = os.path.abspath(os.curdir)
    os.chdir(os.path.join(os.curdir, 'src_raydata_raymatrix',
                          'assemblematrix'))
    os_sys = os.system('./make')
    if not os_sys == 0:
        os.exit("assemblematrix can not be compiled properly")
    os.chdir(cur_dir)

    shutil.copy(os.path.join(os.curdir, 'src_raydata_raymatrix',
                             'assemblematrix', 'assemblematrix'),
                os.path.join(target_add))

####################### run_assemblematrix #############################


def run_assemblematrix(target_add):
    """
    run both raydata and raymatrix in the directory
    """
    cur_dir = os.path.abspath(os.curdir)

    print '\n======>> run assemblematrix at %s' % target_add
    os.chdir(os.path.join(target_add))
    os_sys = os.system('./assemblematrix < in.am')
    if not os_sys == 0:
        print 'assemblematrix was not executed correctly!'

    os.chdir(cur_dir)

####################### station_filter #############################


def station_filter(ls_stas, min_xcorr=-100, max_xcorr=100, min_epi=0.,
                   max_epi=360., check_clip=True):
    """
    Filters the stations based on the required inputs
    NOT BEING USED AT THIS MOMENT!
    """
    empty_array = np.array([])
    pass_stas_corr_1 = ls_stas[min_xcorr <= ls_stas[:, 6].astype(float)]

    if not pass_stas_corr_1.size == 0:
         pass_stas_corr_2 = pass_stas_corr_1[
             max_xcorr > pass_stas_corr_1[:, 6].astype(float)]
    else:
        return empty_array

    passed_stas_epi_1 = pass_stas_corr_2[
        min_epi <= pass_stas_corr_2[:, 4].astype(float)]
    if not passed_stas_epi_1.size == 0:
        passed_stas_epi_2 = passed_stas_epi_1[
            max_epi > passed_stas_epi_1[:, 4].astype(float)]
    else:
        return empty_array
    if passed_stas_epi_2.size == 0:
        return empty_array

    if check_clip:
        passed_stas = passed_stas_epi_2[
            passed_stas_epi_2[:, 19].astype(float) < 0.1]

    return passed_stas

####################### plot_assemble_dtheor #############################


def plot_assemble_dtheor(assembled_dir, out_ident):
    """
    Plot nonzero, rhs and dtheor created by assemblematrix.f
    :param assembled_dir:
    :param out_ident:
    :return:
    """

    plt.ion()

    fio = open(os.path.join(assembled_dir, 'aux.%s' % out_ident), 'r')
    fi = fio.readlines()
    nonzero = []
    rhs = []
    dtheor = []
    for li in range(0, len(fi)):
        if 'End-of-matrix' in fi[li]:
            break
        nonzero.append(int(fi[li].split()[1]))
        rhs.append(float(fi[li].split()[2]))
        dtheor.append(float(fi[li].split()[14]))
    plt.figure()
    plt.plot(range(1, len(nonzero)+1), nonzero, lw=3)
    plt.xlabel('#measurement', size=24, weight='bold')
    plt.ylabel('#nonzero', size=24, weight='bold')
    plt.xticks(size=18, weight='bold')
    plt.yticks(size=18, weight='bold')

    plt.figure()
    plt.plot(range(1, len(rhs)+1), rhs, lw=3)
    plt.xlabel('#measurement', size=24, weight='bold')
    plt.ylabel('RHS', size=24, weight='bold')
    plt.xticks(size=18, weight='bold')
    plt.yticks(size=18, weight='bold')

    plt.figure()
    plt.plot(range(1, len(dtheor)+1), dtheor, lw=3)
    plt.xlabel('#measurement', size=24, weight='bold')
    plt.ylabel('dtheor', size=24, weight='bold')
    plt.xticks(size=18, weight='bold')
    plt.yticks(size=18, weight='bold')

    plt.show()

####################### plot_assemble_stations #############################


def plot_assemble_stations(assembled_dir, out_ident):
    """
    plot two maps: stations and number of usage + dT(av)
    :param assembled_dir:
    :param out_ident:
    :return:
    """

    plt.ion()

    fio = open(os.path.join(assembled_dir,
                            'assemblematrix.stations.%s' % out_ident), 'r')
    fi = fio.readlines()

    stla = []
    stlo = []
    stnn = []
    stdt = []
    for li in range(2, len(fi)):
        stla.append(float(fi[li].split()[6]))
        stlo.append(float(fi[li].split()[7]))
        stnn.append(float(fi[li].split()[2]))
        stdt.append(float(fi[li].split()[4]))

    plt.figure()
    m = Basemap(projection='robin', lon_0=0.0, lat_0=0.0, resolution='c')
    m.fillcontinents()
    m.drawparallels(np.arange(-90., 120., 30.))
    m.drawmeridians(np.arange(0., 420., 60.))
    m.drawmapboundary()

    x, y = m(stlo, stla)
    m.scatter(x, y, s=100, c=stnn, marker="o", edgecolor="none", zorder=10,
            vmax=1000)
    plt.colorbar()
    plt.title('Number of repetition', size=24, weight='bold')

    plt.figure()
    m = Basemap(projection='robin', lon_0=0.0, lat_0=0.0, resolution='c')
    m.fillcontinents()
    m.drawparallels(np.arange(-90., 120., 30.))
    m.drawmeridians(np.arange(0., 420., 60.))
    m.drawmapboundary()

    m.scatter(x, y, s=100, c=stdt, marker="o", edgecolor="none", zorder=10,
            vmin=-7, vmax=7)
    plt.colorbar()
    plt.title('Average dT', size=24, weight='bold')

    plt.show()

####################### plot_assemble_columndensity ##########################


def plot_assemble_columndensity(assembled_dir, out_ident,
                                vertex_file, facet_file):
    """
    plot columndensity created by assemblematrix.f
    :param assembled_dir:
    :param out_ident:
    :param vertex_file:
    :param facet_file:
    :return:
    """

    print '\n======>> Creating VTK file'
    columnden_files = glob.glob(os.path.join(assembled_dir,
                                             'columndensity.*'))
    if not len(columnden_files) == 1:
        sys.exit("Length of columndensity.* files is more than 1!")

    print '------> Load columndensity file'
    mesh_points = np.loadtxt(os.path.join(columnden_files[0]),
                             skiprows=2, comments='#')
    print '------> Load facet file'
    mesh_facets = np.loadtxt(os.path.join(assembled_dir, facet_file),
                             dtype=int, skiprows=1, comments='#')

    print '------> Writing to the disk'
    vtk = pvtk.VtkData(pvtk.UnstructuredGrid(mesh_points[:, 0:3],
                                             tetra=mesh_facets),
                       pvtk.PointData(pvtk.Scalars(mesh_points[:, 3],
                                                   name='column_density')),
                       'Inversion Grid')
    print '\n\n=================='
    print "WARNING: INDEXING!"
    print '=================='
    vtk.tofile(os.path.join(assembled_dir, 'columndensity.vtk'))
