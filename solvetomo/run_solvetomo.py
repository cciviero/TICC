#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  run_solvetomo.py
#   Purpose:   Compile and run mpisolvetomo.f
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
# Required Python modules will be imported in this part.
import datetime
import os
import sys
import os
import glob
import shutil

# -----------------------------------------------------------------------
# ----------------------------- INPUT -----------------------------------
# -----------------------------------------------------------------------

# number of cores
nc = 2
# which mpi to use - previous problem with compiler and mpi... eg. you can chose between
# compiler: 'mpif90.openmpi'
# execute: 'mpirun.openmpi'
# OR
#compiler: 'mpifort'
#exectue: 'mpiexec'

#compiler = 'mpif90.openmpi'
#execute = 'mpirun.openmpi'
compiler = 'mpifort'
execute = 'mpiexec'

# input path; usually the assemble directory
# in_path = '/data/maria/global_tomography/tomo_P_Pdiff_NA_UCL/ass_PortugueseNationalSeismicNetwork_2009-2013'
# in_path = '/data/maria/global_tomography/tomo_P_Pdiff_NA_UCL/_ass_cominded_all'
# in_path = '/data/maria/tomography/tomo_ASAP/ass_2023-02-14'
# in_path = '/data/maria/tomography/P_evolution_ASAP/ass_ISC'
in_path = '/Volumes/EUROPE_DATA/ass_eutomo2012'
# in_path = '/data/maria/tomography/P_evolution_ASAP/ass_ISC_PdOX'
# in_path = '/data/maria/tomography/P_evolution_ASAP/ass_ISC_PdOX1-8'

# where the files created by assemblematrix are moved for inversion
# attention: depending on your smoothing/damping, output will be sorted in a subfolder
# out_path = '/mnt/seismodata2/MT/SEA-SEIS/resolution/ISC_29.05.2020/res_op2_n0'
# out_path = '/data/maria/global_tomography/tomo_P_Pdiff_NA_UCL/inv_PortugueseNationalSeismicNetwork_2009-2013'
# out_path = '/data/maria/global_tomography/tomo_P_Pdiff_NA_UCL/_inv_combined_all'
# out_path = '/data/maria/tomography/tomo_ASAP/inv_2023-02-14'
out_path = '/Volumes/EUROPE_DATA/inv_eutomo2012'
# out_path = '/data/maria/tomography/P_evolution_ASAP/inv_ASAP'
# out_path = '/data/maria/tomography/P_evolution_ASAP/inv_ISC_PdOX1-8'

# ytry for L curve in mpisolvetome.f ~Line 106:
# ytry_lcurve must have exactly 11 elements and the first two are
# not actually read (the program automatically sets them
# to 1.0 and 0.03, respectivelyy
# original:
ytry = [1.0, 0.0, 0.0005, 0.001, 0.0025, 0.003, 0.0035, 0.005, 0.008, 0.05, 0.1]
#ytry = [1.0, 0.0, 0.0001, 0.0003, 0.0005, 0.001, 0.003, 0.005, 0.01, 0.1, 0.3]

# ytry = [1.,0.,.0001,.001,.003,.005,.007,.01,.025,.05,.1]
# ytry = [1., 0., 0.01, .025, .05, .1, .2, .5, 1., 2, 5]
# ytry = [1.0, 0.0, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.05, 0.1, 0.5, 0.9]

# ytry = [1.0, 0.0, 0.0001, 0.0003, 0.0005, 0.001, 0.003, 0.005, 0.1, 0.3, 0.5]
# ytry = [1.0, 0.0, 0.005, 0.007, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]


# vertices and facet name:
vertex_file = 'vertices.h4'
facet_file = 'facets.h4'

# indent the one that you define in assemblematrix
indent = 'P_eutomo'

# EXTRA for naming dlnV?
phase = 'dlnVp'

# smoothing damping parameters: chi^2/N - ksmooth (method of smoothing) -
#                               epsnorm (damping) - epssmooth (smoothness) - ratio
sd = [0, 1, 3.0, 1.0, 0.3]
# sd = [0, 1, 7.0, 1.0, 0.7]

# max mag of outliers
max_outliers = 3.0

# max nr of iterations during search, rootfinding
max_iterations = [5000, 5000]

# 9th line must be zero [XXX WHY??]
line9 = 0

# -----------------------------------------------------------------------
# ------------------------------- RUN -----------------------------------
# -----------------------------------------------------------------------

# =====>>>> 1st: move necessary files to out_path form in_path
if not os.path.isdir(out_path):
    os.makedirs(out_path)

print('Following files are moved...')

# AUX.<INDENT>
print('aux.%s' % indent, '(auxiliary or data file)')
path_2_aux = os.path.join(in_path, 'aux.%s' % indent)
if os.path.isfile(path_2_aux):
    os.rename(path_2_aux, os.path.join(out_path, 'aux.%s' % indent))
elif os.path.isfile(os.path.join(out_path, 'aux.%s' % indent)):
    print ('aux.%s' % indent, 'is already in destination folder')
else:
    print ('aux.%s' % indent, 'cannot be found. EXIT!')
    sys.exit('ERROR')

# MAT.<INDENT>
print ('mat.%s' % indent, '(matrix file)')
path_2_mat = os.path.join(in_path, 'mat.%s' % indent)
if os.path.isfile(path_2_mat):
    os.rename(path_2_mat, os.path.join(out_path, 'mat.%s' % indent))
elif os.path.isfile(os.path.join(out_path, 'mat.%s' % indent)):
    print ('mat.%s' % indent, 'is already in destination folder')
else:
    print ('mat.%s' % indent, 'cannot be found. EXIT!')
    sys.exit('ERROR')

# ATAMX.<INDENT>
print ('atamx.%s' % indent)
path_2_atamx = os.path.join(in_path, 'atamx.%s' % indent)
if os.path.isfile(path_2_atamx):
    os.rename(path_2_atamx, os.path.join(out_path, 'atamx.%s' % indent))
elif os.path.isfile(os.path.join(out_path, 'atamx.%s' % indent)):
    print( 'atamx.%s' % indent, 'is already in destination folder')
else:
    print ('atamx.%s' % indent, 'cannot be found. EXIT!')
    sys.exit('ERROR')


print ('\nFollwoing files are copied...')

try: 
    # VERTEX.<INDENT>
    print ('%s' % vertex_file)
    path_2_vertex = os.path.join(in_path, '%s' % vertex_file)
    if os.path.isfile(path_2_vertex):
        shutil.copy(path_2_vertex, os.path.join(out_path, '%s' % vertex_file))
    elif os.path.isfile(os.path.join(out_path, '%s' % vertex_file)):
        print ('%s' % vertex_file, 'is already in destination folder')
    else:
        print ('%s' % vertex_file, 'cannot be found. EXIT!')
        sys.exit('ERROR')
    
    # FACET.<INDENT>
    print ('%s' % facet_file)
    path_2_facet = os.path.join(in_path, '%s' % facet_file)
    if os.path.isfile(path_2_facet):
        shutil.copy(path_2_facet, os.path.join(out_path, '%s' % facet_file))
    elif os.path.isfile(os.path.join(out_path, '%s' % facet_file)):
        print ('%s' % facet_file, 'is already in destination folder')
    else:
        print ('%s' % facet_file, 'cannot be found. EXIT!')
        sys.exit('ERROR')
except Exception as exp:
    print(exp)
    print('Just continue...')


# =====>>>> 2nd: create in.st file
# information for in.st
# 1. out_path
# 2. indent
# 3. indent
# 4. list of smoothing damping parameters
# 5. max magnitued of outliers allowed
# 6. max nr of iterations during search, rootfinding
# 7. vertice file name
# 8. facet file name
# 9. 0??

print ('Creating in.st file!\n')
outfile_fio = open(os.path.join(out_path, 'in.st'), 'w+')
outfile_fio.writelines('%s\n' % out_path)
outfile_fio.writelines('%s\n' % indent)
outfile_fio.writelines('%s\n' % indent)
outfile_fio.writelines('%s  %s  %s  %s  %s\n' % (sd[0], sd[1], sd[2], sd[3], sd[4]))
outfile_fio.writelines('%s\n' % max_outliers)
outfile_fio.writelines('%s  %s\n' % (max_iterations[0], max_iterations[1]))
outfile_fio.writelines('%s\n' % vertex_file)
outfile_fio.writelines('%s\n' % facet_file)
outfile_fio.writelines('%s\n' % line9)
outfile_fio.writelines('%s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s\n' % (ytry[0], ytry[1], ytry[2], ytry[3], ytry[4],
                                                                ytry[5], ytry[6], ytry[7], ytry[8], ytry[9],
                                                                ytry[10]))
outfile_fio.close()

# =====>>>> 3rd: compile mpisolvetomo.f
if compiler == 'mpifort':
    os_sys = os.system('./makev2')
elif compiler == 'mpif90.openmpi':
    os_sys = os.system('./make_serverv2')
else:
    sys.exit('This compiler is not implemented...EXIT!')
if not os_sys == 0:
    os.exit("mpisolvetomo_v2.f can not be compiled properly")
try:
    os.rename('mpisolvetomov2', os.path.join(out_path, 'mpisolvetomov2'))
except Exception as exp:
    print (exp)
    shutil.move('mpisolvetomov2', os.path.join(out_path, 'mpisolvetomov2'))

# =====>>>> 4th: execute mpisolvetomo
cur_dir = os.path.abspath(os.curdir)
print ('\n======>> run mpisolvetomov2 at %s' % out_path)

os.chdir(os.path.join(out_path))
if execute == 'mpiexec':
    os_sys = os.system('mpiexec -np %s mpisolvetomov2 < in.st' % nc)
elif execute == 'mpirun.openmpi':
    os_sys = os.system('mpirun.openmpi -np %s mpisolvetomov2 < in.st' % nc)
else:
    sys.exit('This mpi is not impolemented...EXIT!')
if not os_sys == 0:
    print ('mpisolvetomov2 was not executed correctly!')
os.chdir(cur_dir)

# =====>>>> 5th: OUTPUT
# What we get:
print('OUTPUT in files in folder: %s' % out_path)
print('solx.%s (solution)' % indent)
print( 'sol.%s (full solution, raw chis)' % indent)
print( 'res.%s (data residuals file)' % indent)
print( 'diagnostics.solve.%s' % indent)
print( 'outl.%s outliers' % indent)
print( 'xy.tradeoff.%s (GMT file for Lcurve)' % indent)
print( 'out.solve.%s' % indent)
print( 'tex.solve.%s \n' % indent)


# =====>>>> 6th: move output files to a designated folder

res_files = glob.glob(os.path.join(out_path, 'res.*.%s' % indent))
sd_folder_name = 'sd_%s' % sd[-1]
sd_folder = os.path.join(out_path, sd_folder_name)
if not os.path.isdir(sd_folder):
    os.mkdir(sd_folder)

os.chdir(os.path.join(out_path))
# all outputs are numbered from 2 to 11
for file_nr in range(2, 12):
    os.rename('res.%02d.%s' % (file_nr, indent), os.path.join(sd_folder, 'res.%02d.%s' % (file_nr, indent)))
    os.rename('sol.%02d.%s' % (file_nr, indent), os.path.join(sd_folder, 'sol.%02d.%s' % (file_nr, indent)))
    os.rename('solx.%02d.%s.%s' % (file_nr, phase, indent), os.path.join(sd_folder, 'solx.%02d.%s.%s' % (file_nr, phase, indent)))
    os.rename('cor.%02d.%s' % (file_nr, indent), os.path.join(sd_folder, 'cor.%02d.%s' % (file_nr, indent)))

os.rename('diagnostics.solve.%s' % indent, os.path.join(sd_folder, 'diagnostics.solve.%s' % indent))
os.rename('outl.%s' % indent, os.path.join(sd_folder, 'outl.%s' % indent))
os.rename('xy.tradeoff.%s' % indent, os.path.join(sd_folder, 'xy.tradeoff.%s' % indent))
os.rename('out.solve.%s' % indent, os.path.join(sd_folder, 'out.solve.%s' % indent))
os.rename('tex.solve.%s' % indent, os.path.join(sd_folder, 'tex.solve.%s' % indent))
os.rename('in.st', os.path.join(sd_folder, 'in.st'))

print ('DONE!')
