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
nc = 30

# which mpi to use - previous problem with compiler and mpi... eg. you can chose between
# compiler: 'mpif90.openmpi'
# execute: 'mpirun.openmpi'
# OR
# compiler: 'mpifort'
# exectue: 'mpiexec'

compiler = 'mpif90.openmpi'
execute = 'mpirun.openmpi'
# compiler = 'mpifort'
# execute = 'mpiexec'

# input path; usually the assemble directory
in_path = '/mnt/seismodata2/MT/RAYTRACER/DETOX_OUTPUT/4TICC/ass_04.08.2021'

# where the files created by assemblematrix are moved for inversion
# attention: depending on your smoothing/damping, output will be sorted in a subfolder
out_path = '/mnt/seismodata2/MT/RAYTRACER/DETOX_OUTPUT/4TICC/inv_04.08.2021'

# vertices and facet name:
vertex_file = 'vertices.h4'
facet_file = 'facets.h4'

# indent
# indent = 'icep'
indent = 'e3dP'

# smoothing damping parameters: chi^2/N - ksmooth (method of smoothing) -
#                               epsnorm (damping) - epssmooth (smoothness) - ratio
sd = [0, 1, 3.0, 10.0, 0.3]

# max mag of outliers
max_outliers = 3.0

# max nr of iterations during search, rootfinding
max_iterations = [6000,6000]

# 9th line must be zero [XXX WHY??]
line9 = 0

# -----------------------------------------------------------------------
# ------------------------------- RUN -----------------------------------
# -----------------------------------------------------------------------

# =====>>>> 1st: move necessary files to out_path form in_path
if not os.path.isdir(out_path):
    os.mkdir(out_path)

print 'Following files are moved...'

# AUX.<INDENT>
print 'aux.%s' % indent, '(auxiliary or data file)'
path_2_aux = os.path.join(in_path, 'aux.%s' % indent)
if os.path.isfile(path_2_aux):
    os.rename(path_2_aux, os.path.join(out_path, 'aux.%s' % indent))
elif os.path.isfile(os.path.join(out_path, 'aux.%s' % indent)):
    print 'aux.%s' % indent, 'is already in destination folder'
else:
    print 'aux.%s' % indent, 'cannot be found. EXIT!'
    sys.exit('ERROR')

# MAT.<INDENT>
print 'mat.%s' % indent, '(matrix file)'
path_2_mat = os.path.join(in_path, 'mat.%s' % indent)
if os.path.isfile(path_2_mat):
    os.rename(path_2_mat, os.path.join(out_path, 'mat.%s' % indent))
elif os.path.isfile(os.path.join(out_path, 'mat.%s' % indent)):
    print 'mat.%s' % indent, 'is already in destination folder'
else:
    print 'mat.%s' % indent, 'cannot be found. EXIT!'
    sys.exit('ERROR')

# ATAMX.<INDENT>
print 'atamx.%s' % indent
path_2_atamx = os.path.join(in_path, 'atamx.%s' % indent)
if os.path.isfile(path_2_atamx):
    os.rename(path_2_atamx, os.path.join(out_path, 'atamx.%s' % indent))
elif os.path.isfile(os.path.join(out_path, 'atamx.%s' % indent)):
    print 'atamx.%s' % indent, 'is already in destination folder'
else:
    print 'atamx.%s' % indent, 'cannot be found. EXIT!'
    sys.exit('ERROR')


print '\nFollwoing files are copied...'

# VERTEX.<INDENT>
print '%s' % vertex_file
path_2_vertex = os.path.join(in_path, '%s' % vertex_file)
if os.path.isfile(path_2_vertex):
    shutil.copy(path_2_vertex, os.path.join(out_path, '%s' % vertex_file))
elif os.path.isfile(os.path.join(out_path, '%s' % vertex_file)):
    print '%s' % vertex_file, 'is already in destination folder'
else:
    print '%s' % vertex_file, 'cannot be found. EXIT!'
    sys.exit('ERROR')

# FACET.<INDENT>
print '%s' % facet_file
path_2_facet = os.path.join(in_path, '%s' % facet_file)
if os.path.isfile(path_2_facet):
    shutil.copy(path_2_facet, os.path.join(out_path, '%s' % facet_file))
elif os.path.isfile(os.path.join(out_path, '%s' % facet_file)):
    print '%s' % facet_file, 'is already in destination folder'
else:
    print '%s' % facet_file, 'cannot be found. EXIT!'
    sys.exit('ERROR')

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

print 'Creating in.st file!\n'
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
outfile_fio.close()

# =====>>>> 3rd: compile mpisolvetomo.f
if compiler == 'mpifort':
    os_sys = os.system('./make')
elif compiler == 'mpif90.openmpi':
    os_sys = os.system('./make_server')
else:
    sys.exit('This compiler is not implemented...EXIT!')
if not os_sys == 0:
    os.exit("mpisolvetomo.f can not be compiled properly")
try:
    os.rename('mpisolvetomo', os.path.join(out_path, 'mpisolvetomo'))
except Exception, exp:
    print exp
    shutil.move('mpisolvetomo', os.path.join(out_path, 'mpisolvetomo'))

# =====>>>> 4th: execute mpisolvetomo
cur_dir = os.path.abspath(os.curdir)
print '\n======>> run mpisolvetomo at %s' % out_path
print out_path
os.chdir(os.path.join(out_path))
if execute == 'mpiexec':
    os_sys = os.system('mpiexec -np %s mpisolvetomo < in.st' % nc)
elif execute == 'mpirun.openmpi':
    os_sys = os.system('mpirun.openmpi -np %s mpisolvetomo < in.st' % nc)
else:
    sys.exit('This mpi is not impolemented...EXIT!')
if not os_sys == 0:
    print 'mpisolvetomo was not executed correctly!'
os.chdir(cur_dir)

# =====>>>> 5th: OUTPUT
# What we get:
print 'OUTPUT in files in folder: %s' % out_path
print 'solx.%s (solution)' % indent
print 'sol.%s (full solution, raw chis)' % indent
print 'res.%s (data residuals file)' % indent
print 'diagnostics.solve.%s' % indent
print 'outl.%s outliers' % indent
print 'xy.tradeoff.%s (GMT file for Lcurve)' % indent
print 'out.solve.%s' % indent
print 'tex.solve.%s \n' % indent


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
    os.rename('solx.%02d.dlnVp.%s' % (file_nr, indent), os.path.join(sd_folder, 'solx.%02d.dlnVp.%s' % (file_nr, indent)))
    os.rename('cor.%02d.%s' % (file_nr, indent), os.path.join(sd_folder, 'cor.%02d.%s' % (file_nr, indent)))

os.rename('diagnostics.solve.%s' % indent, os.path.join(sd_folder, 'diagnostics.solve.%s' % indent))
os.rename('outl.%s' % indent, os.path.join(sd_folder, 'outl.%s' % indent))
os.rename('xy.tradeoff.%s' % indent, os.path.join(sd_folder, 'xy.tradeoff.%s' % indent))
os.rename('out.solve.%s' % indent, os.path.join(sd_folder, 'out.solve.%s' % indent))
os.rename('tex.solve.%s' % indent, os.path.join(sd_folder, 'tex.solve.%s' % indent))

print 'DONE!'