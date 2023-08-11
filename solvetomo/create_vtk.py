#/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  plot_tomo_vtk.py
#   Purpose:   Plot tomographic results in VTK format
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
import glob
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import os
import pyvtk as pvtk
import sys

# ------------------- INPUT -----------------------------
compare_tomo = False
average_tomo = False
damp_smooth_ratio = False

column_density = True

# tomo_files = glob.glob('/Users/maria/PhD/RR_RESULTS/P_18.08.2017_evolution/9_ISC_Pdiff_RRall/columndensity.1.ISC_Pdiff_RRall')
# tomo_files = glob.glob('/Users/maria/PhD/RR_RESULTS/P_32_85_Apr_01/columndensity.1.P_RRx2')

#tomo_files = glob.glob('/Users/maria/PhD/RR_RESULTS/P_18.08.2017_evolution/10_RR_P_PP_Pdiff/columndensity.1.RR_4Pd')
tomo_files = glob.glob('/mnt/seismodata2/MT/SEA-SEIS/resolution/SEA-ISC_29.05.2020/res_op2_n0/columndensity.1.SEA-ISC')
tomo_files.sort()
#tomo_files = glob.glob('/Users/maria/PhD/Codes/__TEST_DATA/OUTPUT/rdrm_rrx2/inversion_dir/volume4vtk.txt')
dir_facet = '/mnt/home_geo/mariat/Codes/TICC/raydata_matrix/pyray_data_matrix/src_raydata_raymatrix/files'

facet_file = 'facets.h4'
vertex_file = 'vertices.h4'
# -------------------------------------------------------

plt.ion()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

tomo_counter = 0
for tomo_file in tomo_files:
    if tomo_file:
        # this is because the column 'hitcount' is saved in the 5th column and not in the 3rd
        if column_density:
            column_index = 4
        else:
            column_index = 3

        print 'Load vertices and facets files...'
        output_file = 'vel_%02i_%s.vtk' % (tomo_counter,
                                           os.path.basename(tomo_file))
        mesh_points = np.loadtxt(os.path.join(dir_facet, vertex_file),
                                 skiprows=2, comments='#')
        mesh_facets = np.loadtxt(os.path.join(dir_facet, facet_file),
                                 dtype=np.int, skiprows=1, comments='#')

        fio_tomo = open(tomo_file, 'r')
        fi_tomo = fio_tomo.readlines()

        print 'Reading %s and extract the velocity structure...' % tomo_file


        mat_val_all = []
        for fi_line in range(2, len(fi_tomo)):
            mat_val_all.append(float(fi_tomo[fi_line].split()[column_index]))

        # import ipdb;
        #
        # ipdb.set_trace()

        if compare_tomo:
            print 'compare_tomo == True'
            fio_tomo_sub = open(tomo_files[1], 'r')
            fi_tomo_sub = fio_tomo_sub.readlines()

            fio_compare = open('text_compare_' + output_file, 'w')
            fio_compare.writelines(fi_tomo_sub[0])
            fio_compare.writelines(fi_tomo_sub[1])
            print 'Reading %s and extract the velocity structure...' \
                  % tomo_files[1]

            mat_val_all_sub = []
            radiuses = []
            for fi_line in range(2, len(fi_tomo_sub)):
                mat_val_all_sub.append(
                    float(fi_tomo_sub[fi_line].split()[3]) -
                    mat_val_all[fi_line - 2])
                radiuses.append(
                    np.sqrt(float(fi_tomo_sub[fi_line].split()[0])**2 +
                            float(fi_tomo_sub[fi_line].split()[1])**2 +
                            float(fi_tomo_sub[fi_line].split()[2])**2))
                fio_compare.writelines('%s    %s    %s    %s\n'
                                       % (fi_tomo_sub[fi_line].split()[0],
                                          fi_tomo_sub[fi_line].split()[1],
                                          fi_tomo_sub[fi_line].split()[2],
                                          mat_val_all_sub[fi_line-2]))
            fio_compare.close()

            mat_val_all_sub = np.array(mat_val_all_sub)
            radiuses = np.array(radiuses)
            mat_selected = mat_val_all_sub[radiuses >= 3482]

            print 'converting the data to VTK format....'
            vtk = pvtk.VtkData(
                pvtk.UnstructuredGrid(mesh_points, tetra=mesh_facets),
                pvtk.PointData(
                    pvtk.Scalars(mat_val_all_sub,
                                 name='%s' % output_file.split('.vtk')[0])),
                'Velocity Structure')
            print 'write the VTK file at:\n%s' \
                  % os.path.join(os.curdir, output_file)
            vtk.tofile(os.path.join(os.curdir, 'compare_' + output_file))
            sys.exit('Compare mode is selected...exit!')

        if average_tomo:
            print 'average_tomo == True'
            mat_val_all_sub = []
            for av_file in range(1, len(tomo_files)):
                fio_tomo_sub = open(tomo_files[av_file], 'r')
                fi_tomo_sub = fio_tomo_sub.readlines()
                print 'Reading %s and extract the velocity structure...' \
                      % tomo_files[1]
                for fi_line in range(2, len(fi_tomo_sub)):
                    if av_file == 1:
                        mat_val_all_sub.append(
                            float(fi_tomo_sub[fi_line].split()[3]) +
                            mat_val_all[fi_line - 2])
                    else:
                        mat_val_all_sub[fi_line - 2] += \
                            float(fi_tomo_sub[fi_line].split()[3])
            mat_val_all_sub = np.array(mat_val_all_sub)/len(tomo_files)
            print 'converting the data to VTK format....'
            vtk = pvtk.VtkData(
                pvtk.UnstructuredGrid(mesh_points, tetra=mesh_facets),
                pvtk.PointData(
                    pvtk.Scalars(mat_val_all_sub,
                                 name='%s' % output_file.split('.vtk')[0])),
                'Velocity Structure')
            print 'write the VTK file at:\n%s' \
                  % os.path.join(os.curdir, output_file)
            vtk.tofile(os.path.join(os.curdir, 'compare_' + output_file))
            sys.exit('Average mode is selected...exit!')

        print 'converting the data to VTK format....'
        vtk = pvtk.VtkData(
            pvtk.UnstructuredGrid(mesh_points, tetra=mesh_facets),
            pvtk.PointData(pvtk.Scalars(
                mat_val_all,
                name='%s' % output_file.split('.vtk')[0])),
            'Velocity Structure')
        print 'write the VTK file at:\n%s' \
              % os.path.join(os.path.dirname(tomo_file), output_file)
        vtk.tofile(os.path.join(os.path.dirname(tomo_file), output_file))
        tomo_counter += 1

print 'DONE!'
