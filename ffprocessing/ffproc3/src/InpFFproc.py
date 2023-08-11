#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  InpFFproc.py
#   Purpose:   Input generator for ffprocessing
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
import configparser
import sys
from src.utility_codes import bc, cprint

# ------------------ input_reader ---------------------------


class InpFFproc:
    def __init__(self, inp_ffproc_file):
        cprint('InpFFproc.py', '[INPUT]', bc.yellow, 'Creating input class!')
        config_ini = configparser.ConfigParser()
        config_ini.read(inp_ffproc_file)
        # ---------------section_phase
        self.ph_phase = eval(config_ini.get('section_phase', 'ph_phase'))
        self.ph_channel = eval(config_ini.get('section_phase', 'ph_channel'))
        self.ph_horizontal = eval(config_ini.get('section_phase', 'ph_horizontal'))
        self.proc_hydrophone = eval(config_ini.get('section_phase', 'proc_hydrophone'))
        self.ph_min_epi = float(config_ini.get('section_phase', 'ph_min_epi'))
        self.ph_max_epi = float(config_ini.get('section_phase', 'ph_max_epi'))
        self.ph_bg = eval(config_ini.get('section_phase', 'ph_bg'))
        self.ph_sampling_rate = float(config_ini.get('section_phase',
                                                     'ph_sampling_rate'))
        self.ph_preset = float(config_ini.get('section_phase', 'ph_preset'))
        self.ph_offset = float(config_ini.get('section_phase', 'ph_offset'))

        self.ph_static_preset = float(config_ini.get('section_phase',
                                                     'ph_static_preset'))
        self.ph_static_offset = float(config_ini.get('section_phase',
                                                     'ph_static_offset'))
        # ---------------section_check
        self.force_process = eval(config_ini.get('section_check',
                                                 'force_process'))
        self.check_evs = eval(config_ini.get('section_check', 'check_evs'))
        self.check_ttime = eval(config_ini.get('section_check',
                                               'check_ttime'))
        self.check_ev_stas = eval(config_ini.get('section_check',
                                                 'check_ev_stas'))
        self.check_cut_syn = eval(config_ini.get('section_check',
                                                 'check_cut_syn'))
        self.check_cut_real = eval(config_ini.get('section_check',
                                                  'check_cut_real'))
        self.check_cut_real_syn = eval(config_ini.get('section_check',
                                                      'check_cut_real_syn'))
        self.check_cut_real_syn_gabor = \
            eval(config_ini.get('section_check', 'check_cut_real_syn_gabor'))
        self.check_bands = config_ini.get('section_check', 'check_bands')
        if len(self.check_bands.split('-')) == 1:
            self.check_min_bands = int(self.check_bands) - 1
            self.check_max_bands = int(self.check_bands) + 1
        elif len(self.check_bands.split('-')) == 2:
            self.check_min_bands = int(self.check_bands.split('-')[0]) - 1
            self.check_max_bands = int(self.check_bands.split('-')[1])
        else:
            cprint('InpFFproc.py', '[ERROR]', bc.dred,
                   'bands: %s are not defined correctly. Format: 1-8'
                   % self.check_bands)
        self.check_min_cc = float(config_ini.get('section_check',
                                                 'check_min_cc'))
        self.check_measured_glob = eval(config_ini.get('section_check',
                                                       'check_measured_glob'))
        # ---------------section_real
        self.real_mode = eval(config_ini.get('section_real', 'real_mode'))
        self.real_path = eval(config_ini.get('section_real', 'real_path'))
        self.real_name_format = eval(config_ini.get('section_real',
                                                    'real_name_format'))
        self.real_id_state = eval(config_ini.get('section_real', 'real_id_state'))
        self.real_id_list = eval(config_ini.get('section_real', 'real_id_list'))
        self.real_id = eval(config_ini.get('section_real', 'real_id'))
        self.real_rect = config_ini.get('section_real', 'real_rect')
        self.real_min_azi = float(config_ini.get('section_real',
                                                 'real_min_azi'))
        self.real_max_azi = float(config_ini.get('section_real',
                                                 'real_max_azi'))
        self.int_real = eval(config_ini.get('section_real',
                                            'int_real'))        
        # ---------------section_syn
        self.syn_mode = eval(config_ini.get('section_syn', 'syn_mode'))
        self.syn_path = eval(config_ini.get('section_syn', 'syn_path'))
        self.syn_name_format = eval(config_ini.get('section_syn',
                                                   'syn_name_format'))
        self.int_syn = eval(config_ini.get('section_syn',
                                           'int_syn'))        
        # ---------------section_event
        self.evproc_mode = eval(config_ini.get('section_event', 'evproc_mode'))
        self.evproc_path = eval(config_ini.get('section_event', 'evproc_path'))
        self.evproc_name_format = eval(config_ini.get('section_event',
                                                      'evproc_name_format'))
        self.evproc_min_mag = float(config_ini.get('section_event',
                                                   'evproc_min_mag'))
        self.evproc_max_mag = float(config_ini.get('section_event',
                                                   'evproc_max_mag'))
        self.evproc_min_depth = float(config_ini.get('section_event',
                                                     'evproc_min_depth'))
        self.evproc_max_depth = float(config_ini.get('section_event',
                                                     'evproc_max_depth'))
        self.evproc_min_year = int(config_ini.get('section_event',
                                                  'evproc_min_year'))
        self.evproc_max_year = int(config_ini.get('section_event',
                                                  'evproc_max_year'))
        self.evproc_year_list = eval(config_ini.get('section_event',
                                                  'evproc_year_list'))
        self.evproc_rect = config_ini.get('section_event',
                                           'evproc_rect')
        self.selected_events = eval(config_ini.get('section_event',
                                                   'selected_events'))
        # ---------------section_stf
        self.stf_mode = eval(config_ini.get('section_stf', 'stf_mode'))
        self.stf_scardec_db = eval(config_ini.get('section_stf', 'stf_scardec_db'))
        self.stf_min_qual = float(config_ini.get('section_stf',
                                                 'stf_min_qual'))
        # ---------------section_measurement
        self.mmeant_mode = eval(config_ini.get('section_measurement',
                                               'mmeant_mode'))
        self.mmeant_clip_time1 = float(config_ini.get('section_measurement',
                                                      'mmeant_clip_time1'))
        self.mmeant_clip_time2 = float(config_ini.get('section_measurement',
                                                      'mmeant_clip_time2'))
        self.mmeant_calc_snr = eval(config_ini.get('section_measurement',
                                                   'mmeant_calc_snr'))
        # ---------------section_output
        self.output_format = eval(config_ini.get('section_output',
                                                 'output_format'))
        self.output_dir = eval(config_ini.get('section_output', 'output_dir'))
        self.save_smgr = eval(config_ini.get('section_output', 'save_smgr'))
        self.only_pre_process = eval(config_ini.get('section_output',
                                                    'only_pre_process'))
        # ---------------section_kernel
        self.kernel_output = eval(config_ini.get('section_kernel',
                                                 'kernel_output'))
        self.kernel_ps = eval(config_ini.get('section_kernel', 'kernel_ps'))
        self.kernel_comp = eval(config_ini.get('section_kernel',
                                               'kernel_comp'))
        self.kernel_misfit = eval(config_ini.get('section_kernel',
                                                 'kernel_misfit'))
        self.kernel_min_cc = float(config_ini.get('section_kernel',
                                                  'kernel_min_cc'))
        self.kernel_stf_gauss = eval(config_ini.get('section_kernel',
                                                    'kernel_stf_gauss'))
        # ---------------section_filter
        self.filt_mode = eval(config_ini.get('section_filter', 'filt_mode'))
        cprint('InpFFproc.py', '[FILTER]', bc.yellow,
               'Mode: %s' % self.filt_mode)
        # ---------------section log_gabor
        if self.filt_mode == 'log-gabor':
            self.lgb_filt_pmax = float(config_ini.get('section_log_gabor',
                                                      'lgb_filt_pmax'))
            self.lgb_filt_nscale = int(config_ini.get('section_log_gabor',
                                                      'lgb_filt_nscale'))
            self.lgb_filt_fmult = float(config_ini.get('section_log_gabor',
                                                       'lgb_filt_fmult'))
            self.lgb_filt_sigmaIfc = float(config_ini.get('section_log_gabor',
                                                          'lgb_filt_sigmaIfc'))
            self.lgb_filt_npad = int(config_ini.get('section_log_gabor',
                                                    'lgb_filt_npad'))
            self.lgb_filt_duration = []
            self.lgb_filt_center_period = []
            self.lgb_filt_energy_frac = \
                float(config_ini.get('section_log_gabor',
                                     'lgb_filt_energy_frac'))
            self.lgb_filt_nlambda = \
                float(config_ini.get('section_log_gabor',
                                     'lgb_filt_nlambda'))
            self.lgb_filt_IR = []
        else:
            cprint('InpFFproc.py', '[ERROR]', bc.dred,
                   '%s has not implemented!' % self.filt_mode)
            sys.exit()
        # ---------------section ami
        self.data_structure_ami = eval(config_ini.get('section_ami', 'data_structure_ami'))

        # ---------------section obs
        self.correction_hydrophone = eval(config_ini.get('section_obs', 'correction_hydrophone'))
        self.correction_seismometer = eval(config_ini.get('section_obs', 'correction_seismometer'))

        self.trace_pre_filt = eval(config_ini.get('section_obs', 'trace_pre_filt'))
        self.trace_unit = eval(config_ini.get('section_obs', 'trace_unit'))
        self.trace_waterlevel = eval(config_ini.get('section_obs', 'trace_waterlevel'))
        self.trace_taper = eval(config_ini.get('section_obs', 'trace_taper'))

        self.particle_motion = eval(config_ini.get('section_obs', 'particle_motion'))
