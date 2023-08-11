#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  ray2vtk.py
#   Purpose:   Generate a VTK file for seismic rays
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import numpy as np
import os

from utils_ray2vtk import path_to_vtk, write_vtk

# ========== INPUT
ls_input_file = ["./inputs/event_0001.txt",
                 "./inputs/event_0002.txt"]

"""
Format of the text file:
event_lat, event_lon, event_depth(km), station_lat, station_lon, phase, value

1. (recommended) it is better to have one file per event. This way, you can
load them separately in Paraview and etc.
2. value is a float number that can be used to color each ray.
"""
# ==========

for inp_file in ls_input_file:
    if not os.path.exists(inp_file):
        print("file: %s does not exist! -- continue with the next event"
              % inp_file)
    print("file: %s" % inp_file)

    L = np.array([])  # len of the ray calculated
    B = np.array([])
    X = np.array([])
    Y = np.array([])
    Z = np.array([])
    all_rays = np.loadtxt(inp_file, dtype="object",
                          comments="#", delimiter=",")

    for one_ray in all_rays:
        evla = float(one_ray[0])
        evlo = float(one_ray[1])
        evdp = float(one_ray[2])
        stla = float(one_ray[3])
        stlo = float(one_ray[4])
        phase = one_ray[5].strip()
        val = float(one_ray[6])

        x, y, z = path_to_vtk(evla, evlo, stla, stlo, evdp, phase)
        if not x: continue

        L = np.append(L, len(x))
        X = np.append(X, x)
        Y = np.append(Y, y)
        Z = np.append(Z, z)
        B = np.append(B, val)

    sum_L = int(np.sum(L))
    range_L = range(sum_L)

    file_base_name = os.path.basename(os.path.abspath(inp_file))
    vtk_file_name = "rays_%s.vtk" % os.path.splitext(file_base_name)[0]
    write_vtk(X, Y, Z, B, L, sum_L, range_L, "vtk", vtk_file_name)
