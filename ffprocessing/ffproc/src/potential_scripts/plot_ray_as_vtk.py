#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Scripts to create a VTK file with a taup-generated path between an
earthquake and a receiver.
Copyright: Simon Stähler, 2016
Author: Simon Stähler, LMU München

Usage:
In [2]: from plot_ray_as_vtk import path_to_vtk
In [3]: path_to_vtk(evlat=54.0, evlon=12.0, reclat=24., reclon=59.0,
                    evdepth_in_km=100, filename='raypath.vtp')
"""

from obspy.taup import TauPyModel
from obspy.geodetics import gps2dist_azimuth
import numpy as np
import vtk


def azi_dist_r_to_xyz(azi, dist, radius):
    # Convert point coordinates given in azimuth, distance and radius (with
    # respect to the north pole) to cartesian coordinates
    lat1 = (np.pi / 2 - dist)
    lon1 = azi

    xyz = np.ndarray(3)
    xyz[0] = np.cos(lon1)*np.cos(lat1)
    xyz[1] = np.sin(lon1)*np.cos(lat1)
    xyz[2] = np.sin(lat1)
    xyz *= radius

    return xyz


def create_rot_mat(srclat, srclon):
    # Create a rotation matrix to rotate a source at coordinates srclat, srclon
    # to the north pole.
    # Uses the rotation matrix as it is defined in AxiSEM.
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


def path_to_vtk(evlat, evlon, reclat, reclon, evdepth_in_km, filename):
    # Create a VTK file with the path from (evlat, evlon, evdepth_in_km) to
    # (reclat, reclon).
    # The file is saved in 'filename'.

    # Calculate event/station distance and azimuth
    distance, azi, bazi = gps2dist_azimuth(evlat, evlon,
                                           reclat, reclon)
    distance /= 111000.0
    print 'Distance: ', distance, ', Azimuth:', azi, ', Backazimuth:', bazi

    # Azimuth is defined clockwise, longitude anti-clockwise (multiply by -1)
    azi *= -1 * np.pi / 180.

    # Azimuth has zero in north, colat in south
    azi += np.pi

    # Calculate path with TauPy
    tp = TauPyModel()
    path = tp.get_ray_paths(source_depth_in_km=evdepth_in_km,
                            distance_in_degree=distance,
                            phase_list='P')

    # Select first path
    path_P = path[0]
    x = []
    y = []
    z = []

    # Create rotation matrix to rotate source from (evlat, evlon) to the
    # northpole and use the inverse. (We later want to rotate a source from the
    # north pole to (evlat, evlon).
    rot_mat = np.linalg.inv(create_rot_mat(srclat=evlat,
                                           srclon=evlon))

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

    # Define VTK Points along the path
    npoints = len(x)
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(npoints)
    for i in range(0, npoints):
        points.SetPoint(i, x[i], y[i], z[i])

    # Create line segments along the path
    lines = vtk.vtkCellArray()
    lines.InsertNextCell(npoints)
    for i in range(0, npoints):
        lines.InsertCellPoint(i)

    # Create Polydata from points and lines
    polygon = vtk.vtkPolyData()
    polygon.SetPoints(points)
    polygon.SetLines(lines)

    # Save PolyData object to VTK file
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(filename)
    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput(polygon)
    else:
        writer.SetInputData(polygon)
    writer.Write()
