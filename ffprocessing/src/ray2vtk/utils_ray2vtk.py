#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  utils_ray2vtk.py
#   Purpose:   Utility codes for ray2vtk.py
#   Author:    Kasra Hosseini, Maria Tsekhmistrenko, Simon Staehler
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import numpy as np
from obspy.geodetics.base import gps2dist_azimuth
from obspy.taup import TauPyModel
import os

# ------------------ azi_dist_r_to_xyz ---------------------------


def azi_dist_r_to_xyz(azi, dist, radius):
    """
    Convert point coordinates given in azimuth, distance and radius (with
    respect to the north pole) to cartesian coordinates
    :param azi:
    :param dist:
    :param radius:
    :return:
    """
    lat1 = (np.pi / 2 - dist)
    lon1 = azi

    xyz = np.ndarray(3)
    xyz[0] = np.cos(lon1)*np.cos(lat1)
    xyz[1] = np.sin(lon1)*np.cos(lat1)
    xyz[2] = np.sin(lat1)
    xyz *= radius

    return xyz

# ------------------ create_rot_mat ---------------------------


def create_rot_mat(srclat, srclon):
    """
    Create a rotation matrix to rotate a source at coordinates srclat, srclon
    to the north pole.
    Uses the rotation matrix as it is defined in AxiSEM.
    :param srclat:
    :param srclon:
    :return:
    """
    srccolat = (90.0 - srclat) * np.pi / 180.0
    srclon *= np.pi / 180.0

    rot_mat = np.ndarray((3, 3))
    rot_mat[0, 0] = np.cos(srccolat) * np.cos(srclon)
    rot_mat[1, 1] = np.cos(srclon)
    rot_mat[2, 2] = np.cos(srccolat)
    rot_mat[1, 0] = np.cos(srccolat) * np.sin(srclon)
    rot_mat[2, 0] = -np.sin(srccolat)
    rot_mat[2, 1] = 0.0
    rot_mat[0, 1] = -np.sin(srclon)
    rot_mat[0, 2] = np.sin(srccolat) * np.cos(srclon)
    rot_mat[1, 2] = np.sin(srccolat) * np.sin(srclon)

    return rot_mat

# ------------------ path_to_vtk ---------------------------


def path_to_vtk(evlat, evlon, reclat, reclon, evdepth_in_km, phase,
                bg_model="iasp91"):
    """
    Create a VTK file with the path from (evlat, evlon, evdepth_in_km) to
    (reclat, reclon).
    :param evlat:
    :param evlon:
    :param reclat:
    :param reclon:
    :param evdepth_in_km:
    :param phase:
    :param bg_model:
    :return:
    """
    # Calculate event/station distance and azimuth
    distance, azi, bazi = gps2dist_azimuth(evlat, evlon, reclat, reclon)
    distance /= 111000.0

    # Azimuth is defined clockwise, longitude anti-clockwise (multiply by -1)
    azi *= -1 * np.pi / 180.

    # Azimuth has zero in north, colat in south
    azi += np.pi

    # Calculate path with TauPy
    tp = TauPyModel(bg_model)
    path = tp.get_ray_paths(source_depth_in_km=evdepth_in_km,
                            distance_in_degree=distance,
                            phase_list=[phase])

    # Select first path
    try:
        path_P = path[0]
    except Exception as e:
        return False, False, False

    x = []
    y = []
    z = []

    # Create rotation matrix to rotate source from (evlat, evlon) to the
    # northpole and use the inverse. (We later want to rotate a source from the
    # north pole to (evlat, evlon).
    rot_mat = np.linalg.inv(create_rot_mat(srclat=evlat,
                                           srclon=evlon))

    if "diff" in path_P.name:
        for pi in range(len(path_P.path)):
            if abs(path_P.path[pi][3] - 2889.0) < 1.:
                d_new = np.linspace(path_P.path[pi]['dist'],
                                    path_P.path[pi+1]['dist'], 10)
                t_new = np.linspace(path_P.path[pi]['time'],
                                    path_P.path[pi+1]['time'], 10)
                p_inp = path_P.path[pi]['p']
                depth_inp = path_P.path[pi]['depth']
                for ni in range(len(d_new)):
                    path_P.path = np.insert(path_P.path, pi+1+ni,
                                            (p_inp, t_new[ni], d_new[ni],
                                             depth_inp))
                break

    # Go along the path
    for step in path_P.path:
        # Convert the values for distance and depth to Cartesian coordinates,
        # assuming a source at the north pole.
        v = azi_dist_r_to_xyz(azi=azi,
                              radius=6371.0-step[3],
                              dist=step[2])

        # Rotate the path, such that the source is at the correct location.
        v_rot = np.dot(v, rot_mat)

        x.append(v_rot[0])
        y.append(v_rot[1])
        z.append(v_rot[2])

    return x, y, z

# ------------------ write_vtk ---------------------------


def write_vtk(X, Y, Z, B, L, sum_L, range_L, output_dir, output_name):
    """
    :param X:
    :param Y:
    :param Z:
    :param B:
    :param L:
    :param sum_L:
    :param range_L:
    :param output_dir:
    :param output_name:
    :return:
    """
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    file_name = os.path.join(output_dir, output_name)

    with open(file_name, 'w') as f:
        f.write('# vtk DataFile Version 2.0\n')
        f.write('VTK file test\n')
        f.write('ASCII\n')
        f.write('DATASET POLYDATA\n')
        f.write('POINTS %s float\n' % int(sum_L))
    f.close()

    with open(file_name, 'a') as f:
        np.savetxt(f, np.c_[X, Y, Z], fmt='%10.3f')
    f.close()

    with open(file_name, 'a') as f:
        f.write('LINES %s %s\n' % (len(L), len(L)+sum_L))
    f.close()

    with open(file_name, 'a') as f:
        len_tmp = 0
        for i in range(len(L)):
            save_array = np.array([])
            save_array = np.append(save_array, (int(L[i])))
            save_array = np.append(save_array,
                                   range_L[len_tmp:int(L[i]) + len_tmp])
            np.savetxt(f, save_array[None], fmt='%10d')
            len_tmp += int(L[i])
    f.close()

    with open(file_name, 'a') as f:
        f.write('CELL_DATA %d\n' % len(L))
        f.write('SCALARS cell_scalars float\n')
        f.write('LOOKUP_TABLE default\n')
    f.close()

    with open(file_name, 'a') as f:
        np.savetxt(f, B, fmt='%10.5f')
    f.close()

# ------------------ event2vtk ---------------------------


def event2vtk(evlats, evlons, evdepths_in_km):
    """
    create a vtk file for events
    :param evlats:
    :param evlons:
    :param evdepths_in_km:
    :return:
    """
    counter = 0
    xyz = []
    for i in range(len(evlats)):
        elat = evlats[i]*np.pi/180.
        elon = evlons[i]*np.pi/180.
        edp = evdepths_in_km[i]
        x = (6371-edp) * np.cos(elat) * np.cos(elon)
        y = (6371-edp) * np.cos(elat) * np.sin(elon)
        z = (6371-edp) * np.sin(elat)
        xyz.append("%s %s %s" % (x, y, z))
        counter += 1
    fio = open('events.vtk', 'w')
    fio.writelines('# vtk DataFile Version 3.0\n')
    fio.writelines('vtk output\n')
    fio.writelines('ASCII\n')
    fio.writelines('DATASET UNSTRUCTURED_GRID\n')
    fio.writelines('POINTS %i float\n' % counter)
    for i in range(len(xyz)):
        fio.writelines('%s\n' % xyz[i])

    fio.writelines('\n')
    fio.writelines('CELLS %i %i\n' % (counter, counter*2))
    for i in range(len(xyz)):
        fio.writelines('1 %s\n' % i)

    fio.writelines('\n')
    fio.writelines('CELL_TYPES %i\n' % counter)
    for i in range(len(xyz)):
        fio.writelines('1\n')
    fio.close()
