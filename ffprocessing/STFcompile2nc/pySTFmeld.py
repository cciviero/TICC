#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  pySTFmeld.py
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# --------------- Import required Modules (Python and Obspy) ------------
# -----------------------------------------------------------------------

# import matplotlib
# matplotlib.use("Agg")

from mpi4py import MPI
import numpy as np
import os
import sys
import time

# %load_ext autoreload
# %autoreload 2
from STFutils import *
from netCDF4 import Dataset

import matplotlib.pylab as plt


# -----------------------------------------------------------------------
# -------------------------- Initialise Code ----------------------------
# -----------------------------------------------------------------------

# =====>>>>> INPUT

input_file = 'input_stf.ini'
stf_input = STFreadit(inp_file=input_file)


# XXX during testing
try:
    os.remove(os.path.join(stf_input.path_out, '%s.nc' % stf_input.nc_name))
except Exception, exp:
    print exp


# -----------------------------------------------------------------------
# ----------------------------- Main Part -------------------------------
# -----------------------------------------------------------------------

cprint('pySTFmeld.py', '[PART I]', bc.lblue,
       'Initialise netCDF with %s' % stf_input.in_ref)

if stf_input.in_ref == 'database':
    stf_database = read_events_dmt(stf_input)
else:
    cprint('pySTFmeld.py', '[PART I] [INPUT ERROR]', bc.red,
           'Not implemented yet. Work harder!')


cprint('pySTFmeld.py', '[PART II]', bc.lblue, 'Update netCDF')

# output_nc = os.path.join(stf_input.path_out)
file_name = os.path.join(stf_input.path_out, '%s.nc' % stf_input.nc_name)
# volcgrp = Dataset(file_name, 'r', format="NETCDF4")

# print '%s' % volcgrp.groups['2015.11.13-20:51:31-0234'].groups['DMT'].variables['catalog'][0]

update(stf_input, file_name)

# output_nc = os.path.join(stf_input.path_out)
# file_name = os.path.join(output_nc, '%s.nc' % stf_input.nc_name)
# volcgrp = Dataset(file_name, 'r', format="NETCDF4")

# -----------------------------------------------------------------------
# ------------------------------- Plot it -------------------------------
# -----------------------------------------------------------------------

# plot_it(volcgrp)

