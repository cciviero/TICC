#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  abaq2inv.py
#   Purpose:   convert ABAQUS format to Karin's one
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import sys

# ######################### INPUT
basename = 'mc_kernel_test'
# #########################

# reading the ABAQUS file
abaq_fi = sys.argv[1]

vertices = []
facets = []
abaq_fio = open(abaq_fi, 'r')
for line in abaq_fio:
    if '*' in line:
        continue
    line = line.split('\n')[0]
    line_split = line.split(',')
    if len(line_split) == 4:
        vertices.append([line_split[1], line_split[2], line_split[3]])
    elif len(line_split) == 5:
        facets.append([int(line_split[1])-1,
                       int(line_split[2])-1,
                       int(line_split[3])-1,
                       int(line_split[4])-1])
    else:
        print 'WARNING: strange length: %s' % line_split
        print 'Dont forget one of the last lines!'

facets_fio = open('facets.%s' % basename, 'w')
facets_fio.writelines('%i\n' % len(facets))
for fac in facets:
    facets_fio.writelines('%s %s %s %s\n' % (fac[0], fac[1], fac[2], fac[3]))
facets_fio.close()

vertices_fio = open('vertices.%s' % basename, 'w')
vertices_fio.writelines('3\n')
vertices_fio.writelines('%i\n' % len(vertices))
for vert in vertices:
    vertices_fio.writelines('%s   %s   %s\n' % (vert[0], vert[1], vert[2]))
vertices_fio.close()
