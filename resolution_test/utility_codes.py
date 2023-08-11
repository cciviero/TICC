#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  utility_codes.py
#   Author:    Kasra Hosseini
#   Email:     kasra.hosseinizad@earth.ox.ac.uk
#   License:   GNU Lesser General Public License, Version 3
# -------------------------------------------------------------------

# -----------------------------------------------------
# ----------------Import required Modules -------------
# -----------------------------------------------------
import copy
import configparser
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import os
from pylab import cm
import pickle
from scipy.interpolate import interp1d
import sys

import warnings
warnings.filterwarnings("ignore")
plt.ioff()

# ------------------ InpSubTomo ---------------------------


class InpSubTomo:
    def __init__(self, inp_tomo_depth):
        config_ini = ConfigParser.ConfigParser()
        config_ini.read(inp_tomo_depth)
        # ---------------plot_params
        self.path_save = eval(config_ini.get('plot_params', 'path_save'))
        self.remove_mean = eval(config_ini.get('plot_params', 'remove_mean'))
        self.ref_bg_model = eval(config_ini.get('plot_params', 'ref_bg_model'))
        self.map_proj = eval(config_ini.get('plot_params', 'map_proj'))
        self.lon_0 = float(config_ini.get('plot_params', 'lon_0'))
        self.lat_0 = float(config_ini.get('plot_params', 'lat_0'))
        # ---------------colorbar
        self.vmin = float(config_ini.get('colorbar', 'vmin'))
        self.vmax = float(config_ini.get('colorbar', 'vmax'))
        self.sel_cmap_name = eval(config_ini.get('colorbar', 'sel_cmap_name'))
        self.sel_cmap_numb = int(config_ini.get('colorbar', 'sel_cmap_numb'))
        # ---------------plate_hotspot
        self.plot_coastlines = eval(config_ini.get('plate_hotspot', 'plot_coastlines'))
        self.plot_plate_bound = eval(config_ini.get('plate_hotspot', 'plot_plate_bound'))
        self.plot_hotspot = eval(config_ini.get('plate_hotspot', 'plot_hotspot'))

        self.plot_recons_plate = eval(config_ini.get('plate_hotspot', 'plot_recons_plate'))
        self.plot_subduction = eval(config_ini.get('plate_hotspot', 'plot_subduction'))
        self.plot_ridge = eval(config_ini.get('plate_hotspot', 'plot_ridge'))
        self.plot_recon_coastline = eval(config_ini.get('plate_hotspot', 'plot_recon_coastline'))

        self.option_recons_time = eval(config_ini.get('plate_hotspot', 'option_recons_time'))
        self.user_recons_time = float(config_ini.get('plate_hotspot', 'user_recons_time'))

        self.plate_sinking_rate_um_low = float(config_ini.get('plate_hotspot', 'plate_sinking_rate_um_low'))
        self.plate_sinking_rate_lm_low = float(config_ini.get('plate_hotspot', 'plate_sinking_rate_lm_low'))
        self.plate_sinking_rate_um_high = float(config_ini.get('plate_hotspot', 'plate_sinking_rate_um_high'))
        self.plate_sinking_rate_lm_high = float(config_ini.get('plate_hotspot', 'plate_sinking_rate_lm_high'))

        self.age_low_plate_color = eval(config_ini.get('plate_hotspot', 'age_low_plate_color'))
        self.age_high_plate_color = eval(config_ini.get('plate_hotspot', 'age_high_plate_color'))
        self.age_low_coastline_color = eval(config_ini.get('plate_hotspot', 'age_low_coastline_color'))
        self.age_high_coastline_color = eval(config_ini.get('plate_hotspot', 'age_high_coastline_color'))
        self.recon_plate_linwidth = float(config_ini.get('plate_hotspot', 'recon_plate_linwidth'))
        self.recon_coast_linwidth = float(config_ini.get('plate_hotspot', 'recon_coast_linwidth'))
        # ---------------interpolation
        self.tar_dirname = eval(config_ini.get('interpolation', 'tar_dirname'))
        self.thresh_max_diff_depth = int(config_ini.get('interpolation', 'thresh_max_diff_depth'))
        self.thresh_min_diff_depth = int(config_ini.get('interpolation', 'thresh_min_diff_depth'))
        # ---------------tomo_models
        self.GyPSuM_P = eval(config_ini.get('tomo_models', 'GyPSuM_P'))
        self.GyPSuM_S = eval(config_ini.get('tomo_models', 'GyPSuM_S'))
        self.PRI_P05 = eval(config_ini.get('tomo_models', 'PRI_P05'))
        self.PRI_S05 = eval(config_ini.get('tomo_models', 'PRI_S05'))
        self.GAP_P4 = eval(config_ini.get('tomo_models', 'GAP_P4'))
        self.savani = eval(config_ini.get('tomo_models', 'savani'))
        self.LLNL_G3Dv3 = eval(config_ini.get('tomo_models', 'LLNL_G3Dv3'))
        self.SEMUCB_WM1 = eval(config_ini.get('tomo_models', 'SEMUCB_WM1'))
        self.MITP_2011 = eval(config_ini.get('tomo_models', 'MITP_2011'))
        self.S40RTS = eval(config_ini.get('tomo_models', 'S40RTS'))
        self.UU_P07 = eval(config_ini.get('tomo_models', 'UU_P07'))
        self.TX2011 = eval(config_ini.get('tomo_models', 'TX2011'))
        self.HMSL_P06 = eval(config_ini.get('tomo_models', 'HMSL_P06'))
        self.HMSL_S06 = eval(config_ini.get('tomo_models', 'HMSL_S06'))
        self.MITP08 = eval(config_ini.get('tomo_models', 'MITP08'))
        self.s10mean = eval(config_ini.get('tomo_models', 's10mean'))
        self.S20RTS = eval(config_ini.get('tomo_models', 'S20RTS'))
        self.DETOX_P01 = eval(config_ini.get('tomo_models', 'DETOX_P01'))
        self.SL2013sv = eval(config_ini.get('tomo_models', 'SL2013sv'))
        self.SP12RTS_P = eval(config_ini.get('tomo_models', 'SP12RTS_P'))
        self.SP12RTS_S = eval(config_ini.get('tomo_models', 'SP12RTS_S'))
        self.SPani_P = eval(config_ini.get('tomo_models', 'SPani_P'))
        self.SPani_S = eval(config_ini.get('tomo_models', 'SPani_S'))
        self.S362ANI_M = eval(config_ini.get('tomo_models', 'S362ANI_M'))
        self.SEMum = eval(config_ini.get('tomo_models', 'SEMum'))
        self.SAW642ANb = eval(config_ini.get('tomo_models', 'SAW642ANb'))
        self.MITP_USA_2016MAY = eval(config_ini.get('tomo_models', 'MITP_USA_2016MAY'))
        self.TX2015 = eval(config_ini.get('tomo_models', 'TX2015'))
        self.threeD2016_09Sv = eval(config_ini.get('tomo_models', 'threeD2016_09Sv'))
        self.Sigloch_NAm_2011 = eval(config_ini.get('tomo_models', 'Sigloch_NAm_2011'))
        self.Zaroli2016 = eval(config_ini.get('tomo_models', 'Zaroli2016'))
        self.SGLOBE_rani = eval(config_ini.get('tomo_models', 'SGLOBE_rani'))
        self.SEISGLOB1 = eval(config_ini.get('tomo_models', 'SEISGLOB1'))
        self.SEISGLOB2 = eval(config_ini.get('tomo_models', 'SEISGLOB2'))
        self.BT_MODEL001 = eval(config_ini.get('tomo_models', 'BT_MODEL001'))
        self.TX2019slab_P = eval(config_ini.get('tomo_models', 'TX2019slab_P'))
        self.TX2019slab_S = eval(config_ini.get('tomo_models', 'TX2019slab_S'))
        self.DETOXP1 = eval(config_ini.get('tomo_models', 'DETOX-P1'))
        self.DETOXP2 = eval(config_ini.get('tomo_models', 'DETOX-P2'))
        self.DETOXP3 = eval(config_ini.get('tomo_models', 'DETOX-P3'))

# ------------------ InpSubVote ---------------------------


class InpSubVote:
    def __init__(self, inp_vote_depth):
        config_ini = ConfigParser.ConfigParser()
        config_ini.read(inp_vote_depth)
        # ---------------plot_params
        self.path_save = eval(config_ini.get('plot_params', 'path_save'))
        self.remove_mean = eval(config_ini.get('plot_params', 'remove_mean'))
        self.ref_bg_model = eval(config_ini.get('plot_params', 'ref_bg_model'))
        self.map_proj = eval(config_ini.get('plot_params', 'map_proj'))
        self.lon_0 = float(config_ini.get('plot_params', 'lon_0'))
        self.lat_0 = float(config_ini.get('plot_params', 'lat_0'))
        # ---------------vote_maps
        self.vote_method = eval(config_ini.get('vote_maps', 'vote_method'))
        self.high_low = eval(config_ini.get('vote_maps', 'high_low'))
        self.vmin_target = int(config_ini.get('vote_maps', 'vmin_target'))
        # ---------------plot_reconstructions
        self.plate_sinking_rate_um_low = eval(config_ini.get('plot_reconstructions',
                                                             'plate_sinking_rate_um_low'))
        self.plate_sinking_rate_um_high = eval(config_ini.get('plot_reconstructions',
                                                              'plate_sinking_rate_um_high'))
        self.plate_sinking_rate_lm_low = eval(config_ini.get('plot_reconstructions',
                                                             'plate_sinking_rate_lm_low'))
        self.plate_sinking_rate_lm_high = eval(config_ini.get('plot_reconstructions',
                                                              'plate_sinking_rate_lm_high'))

        self.plot_recons_plate = eval(config_ini.get('plot_reconstructions', 'plot_recons_plate'))
        self.plot_subduction = eval(config_ini.get('plot_reconstructions', 'plot_subduction'))
        self.plot_ridge = eval(config_ini.get('plot_reconstructions', 'plot_ridge'))
        self.plot_recon_coastline = eval(config_ini.get('plot_reconstructions', 'plot_recon_coastline'))

        self.age_low_plate_color = eval(config_ini.get('plot_reconstructions', 'age_low_plate_color'))
        self.age_high_plate_color = eval(config_ini.get('plot_reconstructions', 'age_high_plate_color'))
        self.age_low_coastline_color = eval(config_ini.get('plot_reconstructions', 'age_low_coastline_color'))
        self.age_high_coastline_color = eval(config_ini.get('plot_reconstructions', 'age_high_coastline_color'))
        # ---------------plate_hotspot
        self.plot_coastline = eval(config_ini.get('plate_hotspot', 'plot_coastline'))
        self.plot_plate_bound = eval(config_ini.get('plate_hotspot', 'plot_plate_bound'))
        self.plot_hotspot = eval(config_ini.get('plate_hotspot', 'plot_hotspot'))
        self.recon_plate_linwidth = eval(config_ini.get('plate_hotspot', 'recon_plate_linwidth'))
        self.recon_coast_linwidth = eval(config_ini.get('plate_hotspot', 'recon_coast_linwidth'))
        # ---------------interpolation
        self.tar_dirname = eval(config_ini.get('interpolation', 'tar_dirname'))
        self.thresh_max_diff_depth = int(config_ini.get('interpolation', 'thresh_max_diff_depth'))
        self.thresh_min_diff_depth = int(config_ini.get('interpolation', 'thresh_min_diff_depth'))
        # ---------------tomo_models
        self.GyPSuM_P = eval(config_ini.get('tomo_models', 'GyPSuM_P'))
        self.GyPSuM_S = eval(config_ini.get('tomo_models', 'GyPSuM_S'))
        self.PRI_P05 = eval(config_ini.get('tomo_models', 'PRI_P05'))
        self.PRI_S05 = eval(config_ini.get('tomo_models', 'PRI_S05'))
        self.GAP_P4 = eval(config_ini.get('tomo_models', 'GAP_P4'))
        self.savani = eval(config_ini.get('tomo_models', 'savani'))
        self.LLNL_G3Dv3 = eval(config_ini.get('tomo_models', 'LLNL_G3Dv3'))
        self.SEMUCB_WM1 = eval(config_ini.get('tomo_models', 'SEMUCB_WM1'))
        self.MITP_2011 = eval(config_ini.get('tomo_models', 'MITP_2011'))
        self.S40RTS = eval(config_ini.get('tomo_models', 'S40RTS'))
        self.UU_P07 = eval(config_ini.get('tomo_models', 'UU_P07'))
        self.TX2011 = eval(config_ini.get('tomo_models', 'TX2011'))
        self.HMSL_P06 = eval(config_ini.get('tomo_models', 'HMSL_P06'))
        self.HMSL_S06 = eval(config_ini.get('tomo_models', 'HMSL_S06'))
        self.MITP08 = eval(config_ini.get('tomo_models', 'MITP08'))
        self.s10mean = eval(config_ini.get('tomo_models', 's10mean'))
        self.S20RTS = eval(config_ini.get('tomo_models', 'S20RTS'))
        self.DETOX_P01 = eval(config_ini.get('tomo_models', 'DETOX_P01'))
        self.SL2013sv = eval(config_ini.get('tomo_models', 'SL2013sv'))
        self.SP12RTS_P = eval(config_ini.get('tomo_models', 'SP12RTS_P'))
        self.SP12RTS_S = eval(config_ini.get('tomo_models', 'SP12RTS_S'))
        self.SPani_P = eval(config_ini.get('tomo_models', 'SPani_P'))
        self.SPani_S = eval(config_ini.get('tomo_models', 'SPani_S'))
        self.S362ANI_M = eval(config_ini.get('tomo_models', 'S362ANI_M'))
        self.SEMum = eval(config_ini.get('tomo_models', 'SEMum'))
        self.SAW642ANb = eval(config_ini.get('tomo_models', 'SAW642ANb'))
        self.MITP_USA_2016MAY = eval(config_ini.get('tomo_models', 'MITP_USA_2016MAY'))
        self.TX2015 = eval(config_ini.get('tomo_models', 'TX2015'))
        self.threeD2016_09Sv = eval(config_ini.get('tomo_models', 'threeD2016_09Sv'))
        self.Sigloch_NAm_2011 = eval(config_ini.get('tomo_models', 'Sigloch_NAm_2011'))
        self.Zaroli2016 = eval(config_ini.get('tomo_models', 'Zaroli2016'))
        self.SGLOBE_rani = eval(config_ini.get('tomo_models', 'SGLOBE_rani'))
        self.SEISGLOB1 = eval(config_ini.get('tomo_models', 'SEISGLOB1'))
        self.SEISGLOB2 = eval(config_ini.get('tomo_models', 'SEISGLOB2'))
        self.BT_MODEL001 = eval(config_ini.get('tomo_models', 'BT_MODEL001'))
        self.TX2019slab_P = eval(config_ini.get('tomo_models', 'TX2019slab_P'))
        self.TX2019slab_S = eval(config_ini.get('tomo_models', 'TX2019slab_S'))
        self.DETOXP1 = eval(config_ini.get('tomo_models', 'DETOX-P1'))
        self.DETOXP2 = eval(config_ini.get('tomo_models', 'DETOX-P2'))
        self.DETOXP3 = eval(config_ini.get('tomo_models', 'DETOX-P3'))

# ------------------ InpSubCross ---------------------------


class InpSubCross:
    def __init__(self, inp_tomo_cross):
        config_ini = ConfigParser.ConfigParser()
        config_ini.read(inp_tomo_cross)
        # ---------------plot_params
        self.path_save = eval(config_ini.get('plot_params', 'path_save'))
        self.remove_mean = eval(config_ini.get('plot_params', 'remove_mean'))
        self.ref_bg_model = eval(config_ini.get('plot_params', 'ref_bg_model'))
        self.vtk_var_name = eval(config_ini.get('plot_params', 'vtk_var_name'))
        # ---------------colorbar
        self.vmin = float(config_ini.get('colorbar', 'vmin'))
        self.vmax = float(config_ini.get('colorbar', 'vmax'))
        self.sel_cmap_name = eval(config_ini.get('colorbar', 'sel_cmap_name'))
        self.sel_cmap_numb = int(config_ini.get('colorbar', 'sel_cmap_numb'))
        # ---------------tomo_models
        self.GyPSuM_P = eval(config_ini.get('tomo_models', 'GyPSuM_P'))
        self.GyPSuM_S = eval(config_ini.get('tomo_models', 'GyPSuM_S'))
        self.PRI_P05 = eval(config_ini.get('tomo_models', 'PRI_P05'))
        self.PRI_S05 = eval(config_ini.get('tomo_models', 'PRI_S05'))
        self.GAP_P4 = eval(config_ini.get('tomo_models', 'GAP_P4'))
        self.savani = eval(config_ini.get('tomo_models', 'savani'))
        self.LLNL_G3Dv3 = eval(config_ini.get('tomo_models', 'LLNL_G3Dv3'))
        self.SEMUCB_WM1 = eval(config_ini.get('tomo_models', 'SEMUCB_WM1'))
        self.MITP_2011 = eval(config_ini.get('tomo_models', 'MITP_2011'))
        self.S40RTS = eval(config_ini.get('tomo_models', 'S40RTS'))
        self.UU_P07 = eval(config_ini.get('tomo_models', 'UU_P07'))
        self.TX2011 = eval(config_ini.get('tomo_models', 'TX2011'))
        self.HMSL_P06 = eval(config_ini.get('tomo_models', 'HMSL_P06'))
        self.HMSL_S06 = eval(config_ini.get('tomo_models', 'HMSL_S06'))
        self.MITP08 = eval(config_ini.get('tomo_models', 'MITP08'))
        self.s10mean = eval(config_ini.get('tomo_models', 's10mean'))
        self.S20RTS = eval(config_ini.get('tomo_models', 'S20RTS'))
        self.DETOX_P01 = eval(config_ini.get('tomo_models', 'DETOX_P01'))
        self.SL2013sv = eval(config_ini.get('tomo_models', 'SL2013sv'))
        self.SP12RTS_P = eval(config_ini.get('tomo_models', 'SP12RTS_P'))
        self.SP12RTS_S = eval(config_ini.get('tomo_models', 'SP12RTS_S'))
        self.SPani_P = eval(config_ini.get('tomo_models', 'SPani_P'))
        self.SPani_S = eval(config_ini.get('tomo_models', 'SPani_S'))
        self.S362ANI_M = eval(config_ini.get('tomo_models', 'S362ANI_M'))
        self.SEMum = eval(config_ini.get('tomo_models', 'SEMum'))
        self.SAW642ANb = eval(config_ini.get('tomo_models', 'SAW642ANb'))
        self.MITP_USA_2016MAY = eval(config_ini.get('tomo_models', 'MITP_USA_2016MAY'))
        self.TX2015 = eval(config_ini.get('tomo_models', 'TX2015'))
        self.threeD2016_09Sv = eval(config_ini.get('tomo_models', 'threeD2016_09Sv'))
        self.Sigloch_NAm_2011 = eval(config_ini.get('tomo_models', 'Sigloch_NAm_2011'))
        self.Zaroli2016 = eval(config_ini.get('tomo_models', 'Zaroli2016'))
        self.SGLOBE_rani = eval(config_ini.get('tomo_models', 'SGLOBE_rani'))
        self.SEISGLOB1 = eval(config_ini.get('tomo_models', 'SEISGLOB1'))
        self.SEISGLOB2 = eval(config_ini.get('tomo_models', 'SEISGLOB2'))
        self.BT_MODEL001 = eval(config_ini.get('tomo_models', 'BT_MODEL001'))
        self.TX2019slab_P = eval(config_ini.get('tomo_models', 'TX2019slab_P'))
        self.TX2019slab_S = eval(config_ini.get('tomo_models', 'TX2019slab_S'))
        self.DETOXP1 = eval(config_ini.get('tomo_models', 'DETOX-P1'))
        self.DETOXP2 = eval(config_ini.get('tomo_models', 'DETOX-P2'))
        self.DETOXP3 = eval(config_ini.get('tomo_models', 'DETOX-P3'))

# -----------plot_sinking_plates--------------------


def plot_sinking_plates(m, base_path='/Volumes/HAL_Storage/backups/20181027_SubMachineModels/models/platerecon/plates',
                        model='Matthews_2016_subduction',
                        tar_depth=2800,
                        sinking_rate_um=1, sinking_rate_lm=1,
                        path_save=False, plate_color='m',
                        linwidth=4, map_proj='cea', option_recons_time='option_rt_2',
                        req_age_inp=1000., return_lon_lat=False):
    """
    Plot subduction/ridges/coastlines
    :param m:
    :param base_path:
    :param model:
    :param tar_depth:
    :param sinking_rate_um:
    :param sinking_rate_lm:
    :param path_save:
    :param plate_color:
    :param linwidth:
    :param map_proj:
    :param option_recons_time:
    :param req_age_inp:
    :param return_lon_lat:
    :return:
    """
    if option_recons_time == 'option_rt_2':
        if float(tar_depth) <= 660:
            req_age = float(tar_depth)/sinking_rate_um
        else:
            req_age = float(660)/sinking_rate_um + \
                      float(float(tar_depth)-660)/sinking_rate_lm
    else:
        req_age = req_age_inp
    if 'Matthews_2016' in model:
        if req_age > 410:
            plt.figure()
            ax = plt.gca()
            ax.axis('off')
            ax.set_title("Matthews_2016_subduction model is selected.\n\n"
                         "Supported age: < %s Ma\nCalculated age: %s Ma"
                         % (410, req_age), loc='left')
            plt.tight_layout()
            plt.savefig(path_save, format='JPEG', bbox_inches='tight')
            sys.exit()
    elif 'Zahirovic_2016' in model:
        if req_age > 230:
            plt.figure()
            ax = plt.gca()
            ax.axis('off')
            ax.set_title("Zahirovic_2016_subduction model is selected.\n\n"
                         "Supported age: < %s Ma\nCalculated age: %s Ma"
                         % (230, req_age), loc='left')
            plt.tight_layout()
            plt.savefig(path_save, format='JPEG', bbox_inches='tight')
            sys.exit()
    elif 'Seton_2012' in model:
        if req_age > 200:
            plt.figure()
            ax = plt.gca()
            ax.axis('off')
            ax.set_title("Seton_2012 model is selected.\n\n"
                         "Supported age: < %s Ma\nCalculated age: %s Ma"
                         % (200, req_age), loc='left')
            plt.tight_layout()
            plt.savefig(path_save, format='JPEG', bbox_inches='tight')
            sys.exit()

    inv_file = np.loadtxt(os.path.join(base_path, '%s_pkl' % model, 'inv_%s.txt' % model),
                          delimiter=',', dtype='object')
    inv_file_cp = copy.deepcopy(inv_file)
    inv_file_cp[:, 0] = inv_file_cp[:, 0].astype(np.float) - req_age
    tar_index = np.argmin(inv_file_cp[inv_file_cp[:, 0] >= 0][:, 0])
    tar_model_path = inv_file_cp[inv_file_cp[:, 0] >= 0][tar_index][1]
    lonlat_pkl = pickle.load(open(os.path.join(base_path, '%s_pkl' % model, tar_model_path), 'r'))
    x, y = m(lonlat_pkl[:, 0], lonlat_pkl[:, 1])

    x[np.isnan(lonlat_pkl[:, 0])] = np.nan
    y[np.isnan(lonlat_pkl[:, 0])] = np.nan
    dx = np.abs(x[1:] - x[:-1])

    if map_proj in ['cea', 'hammer', 'moll', 'robin']:
        dx_thresh = 1e7
    elif map_proj in ['cyl']:
        dx_thresh = 50
    else:
        dx_thresh = 1e7
    ind_jump, = np.nonzero(dx > dx_thresh)
    x[ind_jump] = np.nan
    m.plot(x, y, '-', c=plate_color, lw=linwidth, zorder=100)
    if return_lon_lat:
        return req_age, lonlat_pkl
    else:
        return req_age

# ---------------- plot_hotspots


def plot_hotspots(m, base_path='/Volumes/HAL_Storage/backups/20181027_SubMachineModels/models//platerecon/plates'):
    """
    Plot hotspots
    :param m:
    :param base_path:
    :return:
    """
    global_hotspots = np.load(os.path.join(base_path, 'Whittaker_etal_2013_global_hotspots.npy'))
    deep_hotspots = np.load(os.path.join(base_path, 'Whittaker_etal_2015_deep_plumes.npy'))

    #x, y = m(global_hotspots[:, 0], global_hotspots[:, 1])
    # m.scatter(x, y, zorder=100, c='#d24dff', s=20, linewidth='0.7', edgecolor='k')
    #m.scatter(x, y, zorder=100, c='magenta', s=80, linewidth='0.7', edgecolor='k')

    x, y = m(deep_hotspots[:, 0], deep_hotspots[:, 1])
    # m.scatter(x, y, zorder=100, c='#d24dff', s=80, linewidth='0.7', edgecolor='k')
    m.scatter(x, y, zorder=100, c='magenta', s=200, linewidth='0.7', edgecolor='k', marker='^')

    x, y = m(-111, 44)
    # m.scatter(x, y, zorder=100, c='#d24dff', s=80, linewidth='0.7', edgecolor='k')
    m.scatter(x, y, zorder=100, c='magenta', s=200, linewidth='0.7', edgecolor='k', marker='^')


# -----------_get_colormap--------------------


def _get_colormap(colors, colormap_name, num_colors, divide_by=255.):
    """
    A simple helper function facilitating linear colormap creation.
    :param colors:
    :param colormap_name:
    :param num_colors:
    :return:
    """
    # Sort and normalize from 0 to 1.
    indices = np.array(sorted(colors.iterkeys()))
    normalized_indices = (indices - indices.min()) / indices.ptp()

    # Create the colormap dictionary and return the colormap.
    cmap_dict = {"red": [], "green": [], "blue": []}
    for _i, index in enumerate(indices):
        color = colors[index]
        cmap_dict["red"].append((normalized_indices[_i],
                                 color[0]/divide_by,
                                 color[0]/divide_by))
        cmap_dict["green"].append((normalized_indices[_i],
                                   color[1]/divide_by,
                                   color[1]/divide_by))
        cmap_dict["blue"].append((normalized_indices[_i],
                                  color[2]/divide_by,
                                  color[2]/divide_by))
    return LinearSegmentedColormap(colormap_name, cmap_dict, num_colors)

# -----------mean_weighted--------------------


def mean_weighted(tar_arr, tar_lats):
    """
    Calculate weighted mean (latitude)
    :param tar_arr:
    :param tar_lats:
    :return:
    """
    weights_lats = np.cos(np.pi/180.*tar_lats)
    return np.sum(tar_arr*weights_lats)/np.sum(weights_lats)

# -----------bg_model_interpolate--------------------


def bg_model_interpolate(tar_bg, tar_depth,
                         add_bg_models='../models/tomography/bg_models'):
    """
    inteprolate over the background model
    :param tar_bg:
    :param tar_depth:
    :param add_bg_models:
    :return:
    """
    if tar_bg == 'prem':
        prem_model = np.loadtxt(
            os.path.join(add_bg_models, "PREM_1s.csv"), delimiter=',')
        inter_p = interp1d(prem_model[:, 1], prem_model[:, 3])(tar_depth)
        inter_s = interp1d(prem_model[:, 1], prem_model[:, 5])(tar_depth)
    elif tar_bg == 'prem500':
        prem_500 = np.loadtxt(
            os.path.join(add_bg_models, "PREM500.csv"), delimiter=',')
        inter_p = 0
        inter_s = interp1d(6371. - prem_500[:, 0]/1000., prem_500[:, 3]/1000.)(tar_depth)
    elif tar_bg == 'iasp91':
        iasp91 = np.loadtxt(
            os.path.join(add_bg_models, "IASP91.csv"), delimiter=',')
        inter_p = interp1d(iasp91[:, 0], iasp91[:, 2])(tar_depth)
        inter_s = interp1d(iasp91[:, 0], iasp91[:, 3])(tar_depth)
    elif tar_bg == 'ak135':
        ak135 = np.loadtxt(
            os.path.join(add_bg_models, "AK135F_AVG.csv"), delimiter=',')
        inter_p = interp1d(ak135[:, 0], ak135[:, 2])(tar_depth)
        inter_s = interp1d(ak135[:, 0], ak135[:, 3])(tar_depth)
    elif tar_bg == 'tna-sna':
        tna_sna = np.loadtxt(
            os.path.join(add_bg_models,
                         "StartingVsModel_TNA-SNA-average.txt"))
        inter_p = 0
        inter_s = interp1d(6371. - tna_sna[:, 0], tna_sna[:, 1])(tar_depth)
    elif tar_bg == 'gyp_prem':
        gyp_prem = np.loadtxt(
            os.path.join(add_bg_models, "GyPSumP_PREM.txt"))
        inter_p = interp1d(6371. - gyp_prem[:, 0],
                           gyp_prem[:, 1])(tar_depth)
        inter_s = 0
    elif tar_bg == 'gap':
        gap_model = np.loadtxt(
            os.path.join(add_bg_models, "GAP_bg_model.txt"))
        inter_p = interp1d(gap_model[:, 0], gap_model[:, 1])(tar_depth)
        inter_s = 0
    elif tar_bg == 'semucb':
        gap_model = np.loadtxt(os.path.join(add_bg_models, "SEMUCB.ref"))
        inter_p = interp1d(6371. - gap_model[:, 0]/1000.,
                           gap_model[:, 2]/1000.)(tar_depth)
        inter_s = interp1d(6371. - gap_model[:, 0]/1000.,
                           gap_model[:, 3]/1000.)(tar_depth)
    elif tar_bg == 'TX2011':
        tx_model = np.loadtxt(os.path.join(add_bg_models, "TX2011.ref"))
        inter_p = 0
        inter_s = interp1d(tx_model[:, 0], tx_model[:, 1])(tar_depth)
    elif tar_bg == 'STW105':
        stw_model = np.loadtxt(os.path.join(add_bg_models, "STW105.ref"))
        inter_p = 0
        inter_s = interp1d(6371. - stw_model[:, 0]/1000., stw_model[:, 3]/1000.)(tar_depth)
    elif tar_bg == 'vsref_3D2016_03Sv':
        vsref_3d2016 = np.loadtxt(os.path.join(add_bg_models, "vsref_vsglob_eye_3D2016_03Sv.ref"))
        inter_p = 0
        inter_s = interp1d(vsref_3d2016[:, 0], vsref_3d2016[:, 1])(tar_depth)
    elif tar_bg == 'llnl_g3dv3':
        llnl_g3dv3_ref = np.loadtxt(os.path.join(add_bg_models, "LLNL_G3Dv3.LayerAverages.txt"))
        inter_p = interp1d(llnl_g3dv3_ref[:, 2], llnl_g3dv3_ref[:, 3])(tar_depth)
        inter_s = 0
    return [inter_p, inter_s]

# -----------cmap_gallery--------------------


def cmap_gallery(sel_cmap_name, sel_cmap_numb, plushalf=False):
    """
    gallery of colormaps
    :param sel_cmap_name:
    :param sel_cmap_numb:
    :param plushalf:
    :return:
    """
    if not plushalf:
        move_colorbar = 0
    else:
        move_colorbar = 0.5
    if sel_cmap_name == 'custom-001':
        cmap = _get_colormap({
            -1.0: [127, 0, 0],
            -0.833333333333: [195, 0, 0],
            -0.833333333333: [195, 0, 0],
            -0.666666666667: [255, 8, 0],
            -0.666666666667: [255, 8, 0],
            -0.5: [255, 77, 0],
            -0.5: [255, 77, 0],
            -0.333333333333: [255, 145, 0],
            -0.333333333333: [255, 145, 0],
            -0.166666666667: [255, 255, 140],
            -0.166666666667: [255, 255, 140],
            0.0: [255, 255, 255],
            0.0: [255, 255, 255],
            0.166666666667: [147, 213, 255],
            0.166666666667: [147, 213, 255],
            0.333333333333: [54, 175, 255],
            0.333333333333: [54, 175, 255],
            0.5: [26, 127, 255],
            0.5: [26, 127, 255],
            0.666666666667: [10, 60, 155],
            0.666666666667: [10, 60, 155],
            0.833333333333: [0, 39, 130],
            0.833333333333: [0, 39, 130],
            1.0: [0, 0, 100],
          },
          "colormap_custom", sel_cmap_numb + move_colorbar)
    elif sel_cmap_name == 'custom-002':
        cmap = create_log_based_custom_002(sel_cmap_numb)
    elif sel_cmap_name in ['arctic', 'ETOPO1', 'GMT_haxby', 'GMT_panoply', 'meyers', 'muted-d-09',
                           'cpt_seismic', 'bilbao', 'lajolla', 'oslo']:
        cmap = create_cpt_city_colorbar(sel_cmap_name, sel_cmap_numb + move_colorbar)
        if sel_cmap_name == 'oslo':
            cmap = reverse_colourmap(cmap, sel_cmap_numb + move_colorbar)
    elif sel_cmap_name == 'custom-003':
        cmap = _get_colormap({
        0.0:   [ 0.0157473,      0.00332647 ,      0              ],
     0.0625:   [ 0.299413 ,      0.000366217,      0.000549325    ],
     0.125:    [ 0.508812 ,      0          ,      0              ],
     0.1875:   [ 0.672862 ,      0.139086   ,      0.00270085     ],
     0.25:     [ 0.777096 ,      0.330175   ,      0.000885023    ],
     0.3125:   [ 0.854063 ,      0.510857   ,      0              ],
     0.375:    [ 0.904479 ,      0.690486   ,      0              ],
     0.4375:   [ 0.9514   ,      0.835615   ,      0.449271       ],
     0.5:      [ 0.881376 ,      0.912184   ,      0.818097       ],
     0.5625:   [ 0.60824  ,      0.892164   ,      0.935546       ],
     0.625:    [ 0.327672 ,      0.784939   ,      0.873426       ],
     0.6875:   [ 0.188952 ,      0.641306   ,      0.792096       ],
     0.75:     [ 0.128054 ,      0.492592   ,      0.720287       ],
     0.8125:   [ 0.0552529,      0.345022   ,      0.659495       ],
     0.875:    [ 0        ,      0.216587   ,      0.524575       ],
     0.9375:   [ 0        ,      0.120394   ,      0.302678       ],
     1.0:      [ 0        ,      0          ,      0              ]
     },
          "colormap_custom", sel_cmap_numb + move_colorbar, divide_by=1.)

    elif sel_cmap_name == 'custom-004':
        cmap = _get_colormap({
        0.0:   [ 0.299413 ,      0.000366217,      0.000549325    ],
     0.125:    [ 0.508812 ,      0          ,      0              ],
     0.1875:   [ 0.672862 ,      0.139086   ,      0.00270085     ],
     0.25:     [ 0.777096 ,      0.330175   ,      0.000885023    ],
     0.3125:   [ 0.854063 ,      0.510857   ,      0              ],
     0.375:    [ 0.904479 ,      0.690486   ,      0              ],
     0.4375:   [ 0.9514   ,      0.835615   ,      0.449271       ],
     0.5:      [ 0.881376 ,      0.912184   ,      0.818097       ],
     0.5625:   [ 0.60824  ,      0.892164   ,      0.935546       ],
     0.625:    [ 0.327672 ,      0.784939   ,      0.873426       ],
     0.6875:   [ 0.188952 ,      0.641306   ,      0.792096       ],
     0.75:     [ 0.128054 ,      0.492592   ,      0.720287       ],
     0.8125:   [ 0.0552529,      0.345022   ,      0.659495       ],
     0.875:    [ 0        ,      0.216587   ,      0.524575       ],
     1.0:      [ 0        ,      0.120394   ,      0.302678       ],
     },
          "colormap_custom", sel_cmap_numb + move_colorbar, divide_by=1.)

    elif sel_cmap_name == 'custom-005':
        cmap = _get_colormap({
        0.0:   [ 0.508812 ,      0          ,      0              ],
     0.1875:   [ 0.672862 ,      0.139086   ,      0.00270085     ],
     0.25:     [ 0.777096 ,      0.330175   ,      0.000885023    ],
     0.3125:   [ 0.854063 ,      0.510857   ,      0              ],
     0.375:    [ 0.904479 ,      0.690486   ,      0              ],
     0.4375:   [ 0.9514   ,      0.835615   ,      0.449271       ],
     0.5:      [ 0.881376 ,      0.912184   ,      0.818097       ],
     0.5625:   [ 0.60824  ,      0.892164   ,      0.935546       ],
     0.625:    [ 0.327672 ,      0.784939   ,      0.873426       ],
     0.6875:   [ 0.188952 ,      0.641306   ,      0.792096       ],
     0.75:     [ 0.128054 ,      0.492592   ,      0.720287       ],
     0.8125:   [ 0.0552529,      0.345022   ,      0.659495       ],
     1.0:      [ 0        ,      0.216587   ,      0.524575       ],
     },
          "colormap_custom", sel_cmap_numb + move_colorbar, divide_by=1.)



    elif sel_cmap_name == 'VM-001':
        cmap = _get_colormap({
          # 0.00000: [100, 100, 100],
          0.07143: [4, 45, 104],
          0.14286: [31, 92, 165],
          0.21429: [12, 146, 187],
          0.28571: [56, 182, 196],
          0.35714: [123, 206, 187],
          0.42857: [197, 235, 182],
          0.50000: [236, 250, 180],
          0.57143: [255, 239, 163],
          0.64286: [255, 219, 123],
          0.71429: [255, 180, 82],
          0.78571: [255, 142, 65],
          0.85714: [255, 80, 46],
          0.92857: [193, 7, 39],
          1.00000: [0, 0, 0]
          },
          "colormap_custom", sel_cmap_numb + move_colorbar)
    elif sel_cmap_name == 'VM-002':
        cmap = _get_colormap({
          # 0.00000: [100, 100, 100],
          0.07143: [4, 45, 104],
          0.14286: [31, 92, 165],
          0.21429: [12, 146, 187],
          0.28571: [56, 182, 196],
          0.35714: [123, 206, 187],
          0.42857: [197, 235, 182],
          0.50000: [236, 250, 180],
          0.57143: [255, 239, 163],
          0.64286: [255, 219, 123],
          0.71429: [255, 180, 82],
          0.78571: [255, 142, 65],
          0.85714: [255, 80, 46],
          0.92857: [193, 7, 39],
          # 1.00000: [0, 0, 0]
          },
          "colormap_custom", sel_cmap_numb + move_colorbar)
    else:
        cmap = cm.get_cmap(sel_cmap_name, sel_cmap_numb + move_colorbar)
    return cmap

# -----------create_log_based_custom_002--------------------


def create_log_based_custom_002(sel_cmap_numb=81):
    """
    create log-based custom colormaps
    :param sel_cmap_numb:
    :return:
    """

    rgb_arr = np.array(
        [[0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],
        [0.698 ,  0.0941,  0.1686],

        [0.8392,  0.3765,  0.302 ],
        [0.8392,  0.3765,  0.302 ],
        [0.8392,  0.3765,  0.302 ],
        [0.8392,  0.3765,  0.302 ],
        [0.8392,  0.3765,  0.302 ],
        [0.8392,  0.3765,  0.302 ],
        [0.8392,  0.3765,  0.302 ],
        [0.8392,  0.3765,  0.302 ],
        [0.8392,  0.3765,  0.302 ],
        [0.8392,  0.3765,  0.302 ],

        [0.9569,  0.6471,  0.5098],
        [0.9569,  0.6471,  0.5098],
        [0.9569,  0.6471,  0.5098],
        [0.9569,  0.6471,  0.5098],
        [0.9569,  0.6471,  0.5098],

        [0.9922,  0.8588,  0.7804],
        [0.9922,  0.8588,  0.7804],
        [0.9922,  0.8588,  0.7804],

        [0.9686,  0.9686,  0.9686],
        [0.9686,  0.9686,  0.9686],
        [0.9686,  0.9686,  0.9686],
        [0.9686,  0.9686,  0.9686],
        [0.9686,  0.9686,  0.9686],

        [0.8196,  0.898 ,  0.9412],
        [0.8196,  0.898 ,  0.9412],
        [0.8196,  0.898 ,  0.9412],

        [0.5725,  0.7725,  0.8706],
        [0.5725,  0.7725,  0.8706],
        [0.5725,  0.7725,  0.8706],
        [0.5725,  0.7725,  0.8706],
        [0.5725,  0.7725,  0.8706],

        [0.2627,  0.5765,  0.7647],
        [0.2627,  0.5765,  0.7647],
        [0.2627,  0.5765,  0.7647],
        [0.2627,  0.5765,  0.7647],
        [0.2627,  0.5765,  0.7647],
        [0.2627,  0.5765,  0.7647],
        [0.2627,  0.5765,  0.7647],
        [0.2627,  0.5765,  0.7647],
        [0.2627,  0.5765,  0.7647],
        [0.2627,  0.5765,  0.7647],

        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745],
        [0.1294,  0.4   ,  0.6745]])

    cmap_dict = {"red": [], "green": [], "blue": []}
    for i, ind in enumerate(np.linspace(0, 1, len(rgb_arr))):
        cmap_dict["red"].append((ind, rgb_arr[i, 0], rgb_arr[i, 0]))
        cmap_dict["green"].append((ind, rgb_arr[i, 1], rgb_arr[i, 1]))
        cmap_dict["blue"].append((ind, rgb_arr[i, 2], rgb_arr[i, 2]))
    cmap = LinearSegmentedColormap("colormap_custom", cmap_dict, sel_cmap_numb)
    return cmap

# -----------create_cpt_city_colorbar--------------------


def create_cpt_city_colorbar(sel_cmap_name, sel_cmap_numb):
    """
    create color bar from cpt_city
    :param sel_cmap_name:
    :param sel_cmap_numb:
    :return:
    """
    cmap_dict = gmtColormap(sel_cmap_name)
    cmap = LinearSegmentedColormap("colormap_custom", cmap_dict, sel_cmap_numb)
    return cmap

# -----------gen_list_of_models--------------------


def gen_list_of_models():
    """
    List of tomography models, according to the dirs, e.g. KH
    :return:
    """
    list_of_models = ['GyPSuM-P', 'GyPSuM-S', 'PRI-P05',
                      'PRI-S05', 'GAP-P4', 'savani',
                      'LLNL_G3Dv3', 'SEMUCB-WM1', 'MITP_USA_2011MAR',
                      'S40RTS', 'UU-P07', 'TX2011',
                      'HMSL-P06', 'HMSL-S06', 'MITP08',
                      's10mean', 'S20RTS', 'KH',
                      'SL2013sv', 'SP12RTS-P', 'SP12RTS-S',
                      'SPani-P', 'SPani-S',
                      'S362ANI+M',
                      'SEMum', 'SAW642ANb', 'MIT2016',
                      'TX2015', '3D2016_09Sv', 'Sigloch_NAM',
                      'Zaroli2016', 'SGLOBE_rani',
                      'SEISGLOB1', 'SEISGLOB2',
                      'BT_MODEL001',
                      'TX2019slab-P', 'TX2019slab-S',
                      'DETOX-P1', 'DETOX-P2', 'DETOX-P3'
                      ]
    return list_of_models

# -----------gen_model_names--------------------


def gen_model_names():
    """
    List of tomography model names, e.g. DETOX-P01
    :return:
    """
    model_names = ['GyPSuM-P', 'GyPSuM-S', 'PRI-P05',
                   'PRI-S05', 'GAP-P4', 'savani',
                   'LLNL_G3Dv3', 'SEMUCB-WM1', 'MITP_2011',
                   'S40RTS', 'UU-P07', 'TX2011',
                   'HMSL-P06', 'HMSL-S06', 'MITP08',
                   's10mean', 'S20RTS', 'DETOX-P01',
                   'SL2013sv', 'SP12RTS-P', 'SP12RTS-S',
                   'SPani-P', 'SPani-S',
                   'S362ANI+M',
                   'SEMum', 'SAW642ANb', 'MITP_USA_2016MAY',
                   'TX2015', '3D2016_09Sv', 'Sigloch_NAm_2011',
                   'Zaroli2016', 'SGLOBE-rani',
                   'SEISGLOB1', 'SEISGLOB2',
                   'BT_MDOEL001',
                   'TX2019slab-P', 'TX2019slab-S',
                   'DETOX-P1', 'DETOX-P2', 'DETOX-P3'
                   ]
    return model_names

# -----------gen_list_of_bg--------------------


def gen_list_of_bg():
    """
    List of background models and phase type (0: P, 1: S)
    :return:
    """
    list_of_bg = ['gyp_prem', 'tna-sna', 'iasp91',
                  'iasp91', 'gap', 'prem',
                  'llnl_g3dv3', 'semucb', 'ak135',
                  'prem', 'ak135', 'TX2011',
                  'ak135', 'ak135', 'ak135',
                  'prem', 'prem', 'iasp91',
                  'ak135', 'prem', 'prem',
                  'prem', 'prem',
                  'STW105',
                  'prem', 'prem500', 'ak135',
                  'TX2011', 'vsref_3D2016_03Sv', 'iasp91',
                  'iasp91', 'prem',
                  'prem', 'prem',
                  'prem',
                  'ak135', 'tna-sna',
                  'iasp91', 'iasp91', 'iasp91'
                  ]
    list_of_ph = [0, 1, 0,
                  1, 0, 1,
                  0, 1, 0,
                  1, 0, 1,
                  0, 1, 0,
                  1, 1, 0,
                  1, 0, 1,
                  0, 1,
                  1,
                  1, 1, 0,
                  1, 1, 0,
                  1, 1,
                  1, 1,
                  0,
                  0, 1,
                  0, 0, 0
                  ]
    return list_of_bg, list_of_ph

# -----------shoot--------------------


def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    """
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.

    EPS= 0.00000000005
    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
        alert("Only N-S courses are meaningful, starting at a pole!")

    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2. * np.arctan2 (tu, cf)

    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs (y - c) > EPS):
        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
              d / 4. - cz) * sy * d + tu

    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi

    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)

    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi

    return (glon2, glat2, baz)

# -----------rotate_vector--------------------


def rotate_vector(normal_plane, target_axis):
    """
    rotate a vector
    :param normal_plane:
    :param target_axis:
    :return:
    """
    G = np.array([[np.dot(normal_plane, target_axis),
                   -1*np.linalg.norm(np.cross(normal_plane, target_axis)),
                   0],
                  [np.linalg.norm(np.cross(normal_plane, target_axis)),
                   np.dot(normal_plane, target_axis),
                   0],
                  [0, 0, 1]])
    F = np.array(
        [normal_plane,
         (target_axis -
          np.dot(normal_plane, target_axis)*normal_plane) /
         np.linalg.norm(target_axis -
                        np.dot(normal_plane, target_axis)*normal_plane),
         np.cross(target_axis, normal_plane)]).T
    U = np.dot(F, np.dot(G, np.linalg.inv(F)))
    return U

# -----------convertLatLon2xyz--------------------


def convertLatLon2xyz(lat, lon, eradius=6371.009):
    """
    convert latitude/longitude to xyz
    :param lat:
    :param lon:
    :param eradius:
    :return:
    """
    r = eradius
    # Convert data from lat/lng to x/y/z.
    colat = 90.0 - lat
    x = r * np.sin(np.deg2rad(colat)) * np.cos(np.deg2rad(lon))
    y = r * np.sin(np.deg2rad(colat)) * np.sin(np.deg2rad(lon))
    z = r * np.cos(np.deg2rad(colat))
    return np.array([x, y, z]).T

# -----------gmtColormap--------------------


def gmtColormap(fileName, GMTPath='../helpers/cpt_city'):

      filePath = os.path.join(GMTPath, fileName + '.cpt')

      try:
          f = open(filePath)
      except:
          print("file ",filePath, "not found")
          return None

      lines = f.readlines()
      f.close()

      x = []
      r = []
      g = []
      b = []
      colorModel = "RGB"
      for l in lines:
          ls = l.split()
          if l[0] == "#":
             if ls[-1] == "HSV":
                 colorModel = "HSV"
                 continue
             else:
                 continue
          if ls[0] == "B" or ls[0] == "F" or ls[0] == "N":
             pass
          else:
              x.append(float(ls[0]))
              r.append(float(ls[1]))
              g.append(float(ls[2]))
              b.append(float(ls[3]))
              xtemp = float(ls[4])
              rtemp = float(ls[5])
              gtemp = float(ls[6])
              btemp = float(ls[7])

      x.append(xtemp)
      r.append(rtemp)
      g.append(gtemp)
      b.append(btemp)

      nTable = len(r)
      x = np.array( x , np.float)
      r = np.array( r , np.float)
      g = np.array( g , np.float)
      b = np.array( b , np.float)
      if colorModel == "HSV":
         for i in range(r.shape[0]):
             rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
             r[i] = rr ; g[i] = gg ; b[i] = bb
      if colorModel == "HSV":
         for i in range(r.shape[0]):
             rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
             r[i] = rr ; g[i] = gg ; b[i] = bb
      if colorModel == "RGB":
          r = r/255.
          g = g/255.
          b = b/255.
      xNorm = (x - x[0])/(x[-1] - x[0])

      red = []
      blue = []
      green = []
      for i in range(len(x)):
          red.append([xNorm[i],r[i],r[i]])
          green.append([xNorm[i],g[i],g[i]])
          blue.append([xNorm[i],b[i],b[i]])
      colorDict = {"red":red, "green":green, "blue":blue}
      return (colorDict)

# -----------extract_plates_info--------------------


def extract_plates_info(base_path='../models/platerecon/plates',
                        model='Matthews_2016_subduction',
                        req_age=1000.):
    """
    :param base_path:
    :param model:
    :param req_age_inp:
    :return:
    """
    if 'Matthews_2016' in model:
        if req_age > 410:
            sys.exit()
    elif 'Zahirovic_2016' in model:
        if req_age > 230:
            sys.exit()
    elif 'Seton_2012' in model:
        if req_age > 200:
            sys.exit()

    inv_file = np.loadtxt(os.path.join(base_path, '%s_pkl' % model, 'inv_%s.txt' % model),
                          delimiter=',', dtype='object')
    inv_file_cp = copy.deepcopy(inv_file)
    inv_file_cp[:, 0] = inv_file_cp[:, 0].astype(np.float) - req_age
    tar_index = np.argmin(inv_file_cp[inv_file_cp[:, 0] >= 0][:, 0])
    tar_model_path = inv_file_cp[inv_file_cp[:, 0] >= 0][tar_index][1]
    lonlat_pkl = pickle.load(open(os.path.join(base_path, '%s_pkl' % model, tar_model_path), 'r'))
    return lonlat_pkl

# -----------reverse_colourmap--------------------


def reverse_colourmap(cmap, sel_cmap_numb, name = 'my_cmap_r'):
    """
    In:
    cmap, name
    Out:
    my_cmap_r

    Explanation:
    t[0] goes from 0 to 1
    row i:   x  y0  y1 -> t[0] t[1] t[2]
                   /
                  /
    row i+1: x  y0  y1 -> t[n] t[1] t[2]

    so the inverse should do the same:
    row i+1: x  y1  y0 -> 1-t[0] t[2] t[1]
                   /
                  /
    row i:   x  y1  y0 -> 1-t[n] t[2] t[1]
    """
    reverse = []
    k = []

    for key in cmap._segmentdata:
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:
            data.append((1-t[0],t[2],t[1]))
        reverse.append(sorted(data))

    LinearL = dict(zip(k,reverse))
    my_cmap_r = LinearSegmentedColormap(name, LinearL, sel_cmap_numb)
    return my_cmap_r
