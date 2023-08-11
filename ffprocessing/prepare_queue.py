#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  prepare_queue.py
#   Purpose:   Prepare kernel outputs for a queue system
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
import glob
import numpy as np
import os
from shutil import copyfile

# ------------------ INPUT
output_ffproc_dir = '/net/siglochnas1/volume1/data1/kasrash/2017_measurements/P_Pdiff_10_180_LMU'
evs_name_format = '*.*.*.*'
output_dir = './kernel_dir'
tar_depth = 0

kernel_supermuc = '/gss/scratch/pr63qo/di29kub/kernel_dir'
kernel_bin = '/gpfs/work/pr63qo/di29kub/codes/mc_kernel/bin/mc_kernel'
vertice_add = '/gpfs/work/pr63qo/di29kub/codes/mc_kernel/Meshes/vertices.iter4'
facet_add = '/gpfs/work/pr63qo/di29kub/codes/mc_kernel/Meshes/facets.iter4'
fwd_dir = '/gpfs/work/pr63qo/di29kub/codes/axisem/SOLVER/kernel_5s_fwd_0km'
bwd_dir = '/gpfs/work/pr63qo/di29kub/codes/axisem/SOLVER/iasp91_5s_bwd'
# ------------------

ev_dirs = glob.glob(os.path.join(output_ffproc_dir, evs_name_format))
ev_names = []
for ev_d in ev_dirs:
    ev_names.append(os.path.basename(ev_d))

if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

for ev_n in ev_names:
    kernel_dir = os.path.join(output_dir, 'kernel_%s' % ev_n)

    ev_depth = float(np.loadtxt(os.path.join(
                           output_ffproc_dir, ev_n,
                          'outfiles', 'CMTSOLUTION'), dtype='object', skiprows=6)[0][1])
    if not abs(tar_depth - ev_depth) < 1e-5:
        continue
    print(ev_n)

    if not os.path.isdir(kernel_dir):
        os.mkdir(kernel_dir)

    copyfile(os.path.join(output_ffproc_dir, ev_n,
                          'outfiles', 'CMTSOLUTION'),
             os.path.join(kernel_dir, 'CMTSOLUTION'))

    copyfile(os.path.join(output_ffproc_dir, ev_n,
                          'outfiles', 'kernel_%s.txt' % ev_n),
             os.path.join(kernel_dir, 'kernel_%s.txt' % ev_n))

    copyfile(os.path.join(output_ffproc_dir, ev_n,
                          'outfiles', 'kernel_filter.txt'),
             os.path.join(kernel_dir, 'kernel_filter.txt'))

    copyfile(os.path.join(output_ffproc_dir, ev_n,
                          'outfiles', 'stf.dat'),
             os.path.join(kernel_dir, 'stf.dat'))

    job_cmd = '#@ output=job_$(jobid).out\n'
    job_cmd += '#@ error=job_$(jobid).err\n'
    job_cmd += '#@ job_type=parallel\n'
    job_cmd += '#@ network.MPI=sn_all,not_shared,us\n'
    job_cmd += '#@ notification=always\n'
    job_cmd += '#@ notify_user=hosseini@geophysik.uni-muenchen.de\n'
    job_cmd += '#@ energy_policy_tag=MCKernel\n'
    job_cmd += '#@ minimize_time_to_solution=yes\n'
    # ------------------------------
    job_cmd += '#@ class=micro\n'
    job_cmd += '#@ tasks_per_node=25\n'
    job_cmd += '#@ first_node_tasks=1\n'
    job_cmd += '#@ node=21\n'
    job_cmd += '#@ wall_clock_limit=19:59:00\n'
    # ------------------------------
    initialdir = '%s/kernel_%s' % (kernel_supermuc, ev_n)
    job_cmd += '#@ job_name=kernel_%s\n' % ev_n
    job_cmd += '#@ initialdir=%s\n' % initialdir
    job_cmd += '#@ queue\n'
    job_cmd += '. /etc/profile\n'
    job_cmd += '. /etc/profile.d/modules.sh\n'
    job_cmd += '\n'
    job_cmd += 'module load netcdf/mpi\n'
    job_cmd += 'module load fftw\n'
    job_cmd += '\n'
    job_cmd += 'mkdir %s/Filters\n' % initialdir
    job_cmd += 'mkdir %s/Seismograms\n' % initialdir
    job_cmd += 'cp %s %s\n' % (kernel_bin, initialdir)
    job_cmd += '\n'
    job_cmd += 'poe %s  %s/inparam_%s 2>&1  > %s/OUTPUT_0000\n' \
               % (kernel_bin, initialdir, ev_n, initialdir)

    job_cmd_fio = open(os.path.join(kernel_dir, 'job.cmd'), 'w')
    job_cmd_fio.writelines(job_cmd)
    job_cmd_fio.close()

    inparam_kerner = "FWD_DIR    '%s'\n" % fwd_dir
    inparam_kerner += "BWD_DIR    '%s'\n" % bwd_dir
    inparam_kerner += "PARALLEL_READING    true\n"
    inparam_kerner += "SRC_FILE    '%s/CMTSOLUTION'\n" % initialdir
    inparam_kerner += "REC_FILE    '%s/kernel_%s.txt'\n" % (initialdir, ev_n)
    inparam_kerner += "FILT_FILE    'kernel_filter.txt'\n"
    inparam_kerner += "STF_FILE    'stf.dat'\n"
    inparam_kerner += "MESH_FILE_TYPE    'tetrahedral'\n"
    inparam_kerner += "MESH_FILE_VERTICES        '%s'\n" % vertice_add
    inparam_kerner += "MESH_FILE_FACETS        '%s'\n" % facet_add
    inparam_kerner += "KERNEL_FOR_ABSOLUTE_PERTURBATIONS    false\n"
    inparam_kerner += "NO_INT_OVER_VOLUME    false\n"
    inparam_kerner += "INT_OVER_BACKGROUND_MODEL    false\n"
    inparam_kerner += "INT_OVER_3D_HETEROGENEITIES    false\n"
    inparam_kerner += "HET_FILE    'tests/savani.rtpv'\n"
    inparam_kerner += "OUT_PREFIX    'kerner'\n"
    inparam_kerner += "DUMP_TYPE    'xdmf'\n"
    inparam_kerner += "WRITE_SEISMOGRAMS    true\n"
    inparam_kerner += "ALLOWED_ERROR    1e-17\n"
    inparam_kerner += "ALLOWED_RELATIVE_ERROR    1e-2\n"
    inparam_kerner += "POINTS_PER_MC_STEP    4\n"
    inparam_kerner += "MAXIMUM_ITERATIONS    1000\n"
    inparam_kerner += "WRITE_DETAILED_CONVERGENCE    false\n"
    # XXX BUFFER is defined by the geometry of the mesh and the parameters of the simulation
    # I got these numbers from Simon 
    inparam_kerner += "STRAIN_BUFFER_SIZE    100\n"
    inparam_kerner += "DISPL_BUFFER_SIZE    10\n"
    inparam_kerner += "ELEMENTS_PER_TASK    100\n"
    inparam_kerner += "NO_SORT_MESH_ELEMENTS    false\n"
    inparam_kerner += "USE_PSEUDORANDOM_NUMBERS    false\n"
    inparam_kerner += "MASK_SOURCE_RECEIVER    false\n"
    inparam_kerner += "DAMP_RADIUS_SOURCE_RECEIVER    100.d3\n"
    inparam_kerner += "NO_DECONVOLVE_STF    false\n"
    inparam_kerner += "INTEGRATION_SCHEME    parseval\n"
    inparam_kerner += "FFTW_PLAN    MEASURE\n"
    inparam_kerner += "PLOT_WAVEFIELDS    false\n"
    inparam_kerner += "INT_TYPE    'onvertices'\n"
    inparam_kerner += "CREATE_INTERMEDIATE    true\n"

    inparam_kerner_fio = open(os.path.join(kernel_dir,
                                           'inparam_%s' % ev_n), 'w')
    inparam_kerner_fio.writelines(inparam_kerner)
    inparam_kerner_fio.close()
