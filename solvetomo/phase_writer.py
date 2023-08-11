#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""phase_writer: write 'phase codes' as defined by mapping_dict for each element in d vector"""
__version__ = "0.0.1"
__author__ = "Kasra Hosseini (kasra.hosseinizad@earth.ox.ac.uk)"
__copyright__ = "(C) 2018 Kasra Hosseini. GNU GPL 3."

# TODO:
#
import numpy as np

# ---------------- INPUT
inp_file_assemble_out = '/net/siglochnas1/volume1/data1/mariat/RHUM-RUM/measurement_redo/P_evolution/inversion_dir_ISC_Pdiff_RRall/assemblematrix.out.ISC_Pdiff_RRall'
mapping_dict = \
    {'iscP': 1,
     'P_RR_temp': 2,
     'P_RR_perm': 2,
     'KH_Pdiff': 3
     }

 # 'PM_LMU': 1,
 # 'PM_OX': 1,
 # 'P_RR': 5,
 # 'PPM_LMU': 2,
 # 'PPM_OX': 2,
 # 'PdiM_LMU': 3,
 # 'PdiM_OX': 3,
 # 'iscP': 4,
 # 'iscPn1': 4,
 # 'iscPP': 6,
 # 'iscPKPbcn1': 7

measured_phases = mapping_dict.keys()
# ---------------- END INPUT

print "Reading %s file..." % inp_file_assemble_out
assemble_out_fio = open(inp_file_assemble_out, 'r')
assemble_out_fi = assemble_out_fio.readlines()
phase_arr = []
for i in range(1, len(assemble_out_fi)):
    ass_line = assemble_out_fi[i].split()
    # don't read the lines at the bottom of assemblematrix
    if len(ass_line) < 4:
        if 'demean' in assemble_out_fi[i].lower():
            continue
        else:
            break
    mapped_value = mapping_dict[ass_line[3].split('.')[1].split('_band')[0]]
    phase_arr.extend(int(ass_line[1])*[mapped_value])
print "#measurements (all): %s" % len(phase_arr)
print "Writing phase_type.txt..."
phase_arr = np.array(phase_arr)
np.savetxt('phase_type.txt', phase_arr, fmt='%d')
