#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  vtk_generator.py
#   Purpose:   Generates vtk files (sum of all bands)
#   Author:    Kasra Hosseini (edit:MT)
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

# Required Python modules will be imported in this part.

import glob
import numpy as np
import os
import shutil
import pyvtk as pvtk

# -----------------------------------------------------------------------
# ---------------------------- INPUT ------------------------------------
# -----------------------------------------------------------------------

# where are the band dir from raydata_matrix
path_dir = '/Users/maria/PhD/Codes/__TEST_DATA/OUTPUT/S_32_85_08.03.2017_TEST_the17/OUTPUT_crust2'

# which phase are you using
phase = 'P'

# name of vertex and facet file
vertex_file = 'vertices.h4'
facet_file = 'facets.h4'

# -----------------------------------------------------------------------
# ---------------------------- FUNCTIONS --------------------------------
# -----------------------------------------------------------------------

# ###################### vtk_generator_all #############################


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
                             dtype=np.int, skiprows=1, comments='#')

    mat_val_all = [0.0] * len(mesh_points)
    mat_val_dt_all = [0.0] * len(mesh_points)
    for nj in range(len(ascii_files)):
        print '\n------> create VTK file %s' % ascii_files[nj]
        fmatrix = open(os.path.join(ascii_files[nj]), 'r')
        fmatrix_r = fmatrix.readlines()
        mat_indx = []
        mat_val = []
        mat_val_dt = []
        counter = 0
        for j in range(5, len(fmatrix_r)):
            if np.mod(j, 1000) == 0:
                print '%s/%s' % (j, len(fmatrix_r))
            if counter == 3:
                counter = 0
                continue
            if counter == 0:
                # mat_indx: matrix index
                mat_indx_tmp = np.array(fmatrix_r[j].split()).astype(np.int)
                # IF indexing in ASCII file starts from 0,
                # we do not need -1.
                # Otherwise we need it!
                mat_indx = np.append(mat_indx, mat_indx_tmp - 1)
            if counter == 1:
                # mat_val: matrix value
                mat_val_tmp = np.array(fmatrix_r[j].split()).astype(np.float)

                ddif = (float(fmatrix_r[j - 3].split()[3]) -
                        (float(fmatrix_r[j - 3].split()[4]) +
                         float(fmatrix_r[j - 3].split()[6]) +
                         float(fmatrix_r[j - 3].split()[7]) +
                         float(fmatrix_r[j - 3].split()[8]))) / float(fmatrix_r[j - 3].split()[5])

                mat_val_dt_tmp = np.array(mat_val_tmp).astype(np.float) / np.sum(
                    np.array(mat_val_tmp).astype(np.float)) * ddif
                mat_val = np.append(mat_val, np.abs(mat_val_tmp))
                mat_val_dt = np.append(mat_val_dt, mat_val_dt_tmp)
            counter += 1

        if 'p' in phase.lower():
            for i in range(len(mat_indx)):
                mat_val_all[int(mat_indx[i])] += mat_val[i]
                mat_val_dt_all[int(mat_indx[i])] += mat_val_dt[i]
        elif 's' in phase.lower():
            for i in range(len(mat_indx)):
                mat_val_all[int(mat_indx[i]-len(mat_val_all))] += mat_val[i]
                mat_val_dt_all[int(mat_indx[i]-len(mat_val_all))] += mat_val_dt[i]
        for i in range(len(mesh_points)):
            num_repetition = np.count_nonzero(mat_indx == np.array(mesh_points[i]))
            if num_repetition < 1:
                continue
            mat_val_dt_all /= num_repetition

    vtk = pvtk.VtkData(pvtk.UnstructuredGrid(mesh_points,
                                             tetra=mesh_facets),
                       pvtk.PointData(pvtk.Scalars(mat_val_all,
                                                   name='kernel_value')),
                       'Inversion Grid')
    print '\n\n=================='
    print "WARNING: INDEXING!"
    print '=================='
    vtk.tofile(os.path.join(direname, 'global_all.vtk'))

    vtk = pvtk.VtkData(pvtk.UnstructuredGrid(mesh_points,
                                             tetra=mesh_facets),
                       pvtk.PointData(pvtk.Scalars(mat_val_dt_all,
                                                   name='dt_projected')),
                       'Inversion Grid')
    vtk.tofile(os.path.join(direname, 'global_all_dt.vtk'))


# ###################### vtk_val_azi #############################


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
                             dtype=np.int, skiprows=1, comments='#')
    mesh_points_2 = mesh_points ** 2
    rad = np.sqrt(mesh_points_2[:, 0] + mesh_points_2[:, 1] + mesh_points_2[:, 2])

    mat_val_all = np.array([0.0] * len(mesh_points))
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
                azi = float(fmatrix_r[j].split()[17]) * 180. / np.pi
                # print azi
            if counter == 2:
                # mat_indx: matrix index
                mat_indx_tmp = fmatrix_r[j].split()
                for i in range(len(mat_indx_tmp)):
                    # IF indexing in ASCII file starts from 0,
                    # we do not need -1.
                    # Otherwise we need it!
                    mat_indx.append(int(mat_indx_tmp[i]) - 1)
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

        mat_1 = np.array([0.0] * len(mesh_points))
        mat_1[mat_indx_arr[(0 <= mat_azi_arr) & (mat_azi_arr < 45)]] = \
            np.bincount(mat_indx_arr[(0 <= mat_azi_arr) & (mat_azi_arr < 45)])
        mat_azi_qual[:, 0] += mat_1

        mat_1 = np.array([0.0] * len(mesh_points))
        mat_1[mat_indx_arr[(45 <= mat_azi_arr) & (mat_azi_arr < 90)]] = \
            np.bincount(mat_indx_arr[(45 <= mat_azi_arr) & (mat_azi_arr < 90)])
        mat_azi_qual[:, 1] += mat_1

        mat_1 = np.array([0.0] * len(mesh_points))
        mat_1[mat_indx_arr[(90 <= mat_azi_arr) & (mat_azi_arr < 135)]] = \
            np.bincount(mat_indx_arr[(90 <= mat_azi_arr) & (mat_azi_arr < 135)])
        mat_azi_qual[:, 2] += mat_1

        mat_1 = np.array([0.0] * len(mesh_points))
        mat_1[mat_indx_arr[(135 <= mat_azi_arr) & (mat_azi_arr < 180)]] = \
            np.bincount(mat_indx_arr[(135 <= mat_azi_arr) & (mat_azi_arr < 180)])
        mat_azi_qual[:, 3] += mat_1

        mat_1 = np.array([0.0] * len(mesh_points))
        mat_1[mat_indx_arr[(180 <= mat_azi_arr) & (mat_azi_arr < 225)]] = \
            np.bincount(mat_indx_arr[(180 <= mat_azi_arr) & (mat_azi_arr < 225)])
        mat_azi_qual[:, 0] += mat_1

        mat_1 = np.array([0.0] * len(mesh_points))
        mat_1[mat_indx_arr[(225 <= mat_azi_arr) & (mat_azi_arr < 270)]] = \
            np.bincount(mat_indx_arr[(225 <= mat_azi_arr) & (mat_azi_arr < 270)])
        mat_azi_qual[:, 1] += mat_1

        mat_1 = np.array([0.0] * len(mesh_points))
        mat_1[mat_indx_arr[(270 <= mat_azi_arr) & (mat_azi_arr < 315)]] = \
            np.bincount(mat_indx_arr[(270 <= mat_azi_arr) & (mat_azi_arr < 315)])
        mat_azi_qual[:, 2] += mat_1

        mat_1 = np.array([0.0] * len(mesh_points))
        mat_1[mat_indx_arr[(315 <= mat_azi_arr) & (mat_azi_arr < 360)]] = \
            np.bincount(mat_indx_arr[(315 <= mat_azi_arr) & (mat_azi_arr < 360)])
        mat_azi_qual[:, 3] += mat_1
        plt.ion()
        plt.figure()
        plt.subplot(2, 2, 2)
        plt.plot(mat_azi_qual[:, 0])
        plt.subplot(2, 2, 4)
        plt.plot(mat_azi_qual[:, 1])
        plt.subplot(2, 2, 3)
        plt.plot(mat_azi_qual[:, 2])
        plt.subplot(2, 2, 1)
        plt.plot(mat_azi_qual[:, 3])
        plt.show()

        # print 'Length of mat_indx: %s' % len(mat_indx)
        mat_val_all[mat_indx] += np.array(mat_val)
        # for i in range(len(mat_indx)):
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
    # mat_azi_qual = np.zeros([len(mesh_points), 4])
    # for ep in mat_azi_all.keys():
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
    plt.subplot(2, 1, 1)
    plt.plot(mat_val_all)
    plt.subplot(2, 1, 2)
    plt.plot(mat_azi_qual_rank)
    plt.show()

    plt.figure()
    plt.subplot(1, 1, 1)
    plt.plot(mat_val_all[rad >= 3482] * mat_azi_qual_rank[rad >= 3482])
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



# -----------------------------------------------------------------------
# ------------------------------ RUN ------------------------------------
# -----------------------------------------------------------------------


# =====>>>>> 1st: put all ascii files in one folder

# create an ascii_all folder
ascii_all = os.path.join(path_dir, 'ascii_all')
if not os.path.exists(ascii_all):
    os.makedirs(ascii_all)

# mv from all dirs the ascii files in one folder
all_dirs = glob.glob(os.path.join(path_dir, '*_dir'))
for dir in all_dirs:
    file_cur_dir = glob.glob(os.path.join(dir, 'ascii.*'))
    try:
        file_name = os.path.basename(file_cur_dir[0])
    except Exception, exp:
        continue
    file_new_dir = os.path.join(ascii_all, file_name)
    os.rename(file_cur_dir[0], file_new_dir)

# copy just once the vertices and facet file there
cur_vertex_dir = os.path.join(all_dirs[0], vertex_file)
new_vertex_dir = os.path.join(ascii_all, vertex_file)

cur_facet_dir = os.path.join(all_dirs[0], facet_file)
new_facet_dir = os.path.join(ascii_all, facet_file)

shutil.copy(cur_vertex_dir, new_vertex_dir)
shutil.copy(cur_facet_dir, new_facet_dir)

# run vtk_generator_all
vtk_generator_all(ascii_all, vertex_file, facet_file)


