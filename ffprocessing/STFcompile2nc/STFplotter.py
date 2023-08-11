#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  pySTFmeld.py
#   License:   GPLv3 
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# --------------- Import required Modules (Python and Obspy) ------------
# -----------------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")

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

input_file = 'input_stf.ini'
stf_input = STFreadit(inp_file=input_file)

output_nc = os.path.join(stf_input.path_out)
file_name = os.path.join(output_nc, '%s.nc' % stf_input.nc_name)
volcgrp = Dataset(file_name, 'r', format="NETCDF4")

# -----------------------------------------------------------------------
# ------------------------------- Plot it -------------------------------
# -----------------------------------------------------------------------

plot_it(volcgrp)
