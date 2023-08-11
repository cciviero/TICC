#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  rdrm4llnl.py
#   Purpose:   
#   Author:    Kasra Hosseini + Maria Tsekhmistrenko
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

# Required Python modules will be imported in this part.

# next two lines are necessary for running it on nas
# import matplotlib
# matplotlib.use("Agg")

import datetime
import glob
import multiprocessing
import numpy as np
import os
import pickle
import shutil
import sys

import utils as util
from utils import bc

# ------------------- INPUT -----------------------------
if not sys.argv[1]:
    sys.exit('No input file is entered!')
else:
    input_raydata_matrix = sys.argv[1]
# ------------------- END INPUT -------------------------

# ------------------- READING INPUT
inp = util.InpClass(input_raydata_matrix)
# ------------------- Loop over the bands
read_band = 0
req_band_prev = False
isc_iteration = 1

outfile_fio = open(os.path.join(inp.out_path, 'exitus_info.txt'), 'w+')
# import ipdb; ipdb.set_trace()
for req_band in inp.req_bands:
    if read_band % 2 == 1:
        fill_in_mode = True
    else:
        fill_in_mode = False
        req_band_prev = req_band
    # (F1) if fill_in_mode should be diabled:
    fill_in_mode = False
    req_band_prev = req_band
    # ----------- print input
    util.print_inp(inp, req_band)
    util.cprint('raydata_raymatrix.py', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', bc.white, '')
    util.cprint('raydata_raymatrix.py', '[PROGRAM START]', bc.bblue, '')
    t_start = datetime.datetime.now()

    # ################### ####################

    input_file_names = glob.glob(os.path.join(inp.events_dir,f'{inp.input_file_name_part}_{inp.req_bands[0]}_1'))
    # input_file_names.sort(reverse=True)
    input_file_names.sort()
    isc_iteration = len(input_file_names)
    #import ipdb; ipdb.set_trace()

    for isc_item in range(0, isc_iteration):
        util.cprint('raydata_raymatrix.py', 'input files (raydata and raymatrix) generator', bc.white, '')
        if not inp.parallel_exec:
            util.compile_raydata_raymatrix(inp=inp)
            if req_band == 'band00':
                #import ipdb; ipdb.set_trace()
                
                shutil.copy(input_file_names[isc_item], inp.out_path)
                input_file_name = os.path.basename(input_file_names[isc_item])
            
            # elif req_band == 'band00':
            #     pass
            else:
                pass

            util.raydata_input(inp=inp, input_file_name=input_file_name)
            util.raymatrix_input(inp=inp, input_file_name=input_file_name)
            # ------------ prepare output directory
            util.prepare_dir(input_file_name=input_file_name, inp=inp,
                             vertices=inp.vertex_file,
                             facets=inp.facet_file)
            util.run_raydata_raymatrix(input_file_name=input_file_name, inp=inp,
                                       raydata=inp.run_raydata,
                                       raymatrix=inp.run_raymatrix)
            util.cprint('raydata_raymatrix.py', '[FINISHED]', bc.bgreen, 'raydata and raymatrix finished')

            if inp.run_raymatrix:
                util.mat2asc_run(input_file_name=input_file_name, inp=inp)
            # if not inp.plot_statistics:
            #     mat_val_all = util.vtk_generator(inp=inp, req_band=req_band, len_dirs=1)

            # import ipdb; ipdb.set_trace()
        
        
        
        else:
            util.compile_raydata_raymatrix(inp=inp)
            if not np.mod(len(filt_array), inp.np_req) == 0:
                sub_filt_array = np.split(
                    filt_array[0:-np.mod(len(filt_array), inp.np_req)], inp.np_req)
                sub_filt_array.append(
                    filt_array[-np.mod(len(filt_array), inp.np_req):])
            else:
                sub_filt_array = np.split(filt_array, float(inp.np_req))
            par_job = []
            for nj in range(len(sub_filt_array)):
                input_file_name = '%s_%s_%s' % (inp.input_file_name_part,
                                                req_band,
                                                nj+1)
                par_job.append(
                    multiprocessing.Process(
                        target=util.parallel_raydata_raymatrix,
                        args=(sub_filt_array[nj], input_file_name,
                              req_band, inp)))
            for nj in range(len(par_job)):
                par_job[nj].start()
            for nj in range(len(par_job)):
                util.check_par_jobs(par_job, 0.2)
            mat2asc_par_job = []
            for nj in range(len(sub_filt_array)):
                input_file_name = '%s_%s_%s' % (inp.input_file_name_part,
                                                req_band,
                                                nj+1)
                mat2asc_par_job.append(
                    multiprocessing.Process(target=util.mat2asc_run,
                                            args=(input_file_name, inp)))
                # outread.mat2asc_run(input_file_name=input_file_name)

            for nj in range(len(mat2asc_par_job)):
                mat2asc_par_job[nj].start()
            util.check_par_jobs(mat2asc_par_job, 0.2)

            # mat_val_all = util.vtk_generator(inp=inp,
            #                    req_band=req_band,
            #                    len_dirs=len(sub_filt_array))

    if not req_band == 'band00':
        print ('\n')
        util.cprint('raydata_raymatrix.py', '[#EVENT-STATION PAIRS]', bc.bmagenta, '%s' % len(filt_array))

    util.cprint('raydata_raymatrix.py', '[REGULAR END]', bc.white, '%s sec' % (datetime.datetime.now() - t_start))
    util.cprint('raydata_raymatrix.py', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', bc.white, '')
    print ('\n\n')

    outfile_fio.writelines('>>> %s <<<\n' % req_band)
    if not req_band == 'band00':
        print ('\n')
        outfile_fio.writelines('EVENT-STATION PAIRS: %s \n' % len(filt_array))
    # if inp.parallel_exec:
    #     outfile_fio.writelines('Sum over all elems: %s \n' % mat_val_all)
    outfile_fio.writelines('========================= \n')
    outfile_fio.writelines('All events   : %s \n' % len(passed_event_adds))
    outfile_fio.writelines('Used events  : %s \n' % succ_event_count)
    outfile_fio.writelines('Missed events: %s \n' % (len(passed_event_adds) - succ_event_count))
    outfile_fio.writelines('========================= \n\n')


outfile_fio.writelines('========================= \n')
outfile_fio.writelines('INPUT CHECK \n')
outfile_fio.writelines('========================= \n')

outfile_fio.writelines('plot_statistics = %s \n' % inp.plot_statistics)
outfile_fio.writelines('run_raydata = %s \n' % inp.run_raydata)
outfile_fio.writelines('run_raymatrix = %s \n' % inp.run_raymatrix)
outfile_fio.writelines('parallel_exec = %s \n' % inp.parallel_exec)
outfile_fio.writelines('np_req = %s \n' % inp.np_req)
outfile_fio.writelines('phase = %s \n' % inp.phase)
outfile_fio.writelines('events_dir = %s \n' % inp.events_dir)
outfile_fio.writelines('min_epi = %s \n' % inp.min_epi)
outfile_fio.writelines('max_epi = %s \n' % inp.max_epi)
outfile_fio.writelines('req_bands = %s \n' % inp.req_bands)
outfile_fio.writelines('all_events = %s \n' % inp.all_events)
outfile_fio.writelines('input_file_name_part = %s \n' % inp.input_file_name_part)
outfile_fio.writelines('twinned = %s \n' % inp.twinned)
outfile_fio.writelines('max_num_arrival = %s \n' % inp.max_num_arrival)
outfile_fio.writelines('delay_wrt_first_arrival = %s \n' % inp.delay_wrt_first_arrival)
outfile_fio.writelines('corr_io_list = %s \n' % inp.corr_io_list)
outfile_fio.writelines('crust1 = %s \n' % inp.crust_corr)
outfile_fio.writelines('vp_vs_Qs = %s \n' % inp.vp_vs_Qs)
outfile_fio.writelines('kernel_quad_km = %s \n' % inp.kernel_quad_km)
outfile_fio.writelines('vertex_file = %s \n' % inp.vertex_file)
outfile_fio.writelines('facet_file = %s \n' % inp.facet_file)
outfile_fio.writelines('min_depth = %s \n' % inp.min_depth)
outfile_fio.writelines('max_depth = %s \n' % inp.max_depth)
outfile_fio.writelines('min_xcorr = %s \n' % inp.min_xcorr)
outfile_fio.writelines('max_xcorr = %s \n' % inp.max_xcorr)
outfile_fio.writelines('bg_model = %s \n' % inp.bg_model)
outfile_fio.writelines('check_clip = %s \n' % inp.check_clip)
outfile_fio.writelines('out_path = %s \n' % inp.out_path)

outfile_fio.close()

# raw_input('\nPlease press enter to finish the program!')