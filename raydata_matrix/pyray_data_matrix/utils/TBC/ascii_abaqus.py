"""
Simple script to change the format from ASCII tetrahedron (can be created by
Karin's code) into ABAQUS format which seems to work better for Kernel code
"""

import numpy as np
import os

# ------------ INPUT
par_add = '/home/hosseini/FFINVERSION/INVERSION/Programs/pyfiles/' \
          'pyray_data_matrix/src_raydata_raymatrix/files'
tet_name = 'ppdiff_homo_500'
# ------------ END INPUT

vertices_fi = os.path.join(par_add, 'vertices.%s' % tet_name)
facets_fi = os.path.join(par_add, 'facets.%s' % tet_name)

vertices = np.loadtxt(vertices_fi, skiprows=2)
facets = np.loadtxt(facets_fi, skiprows=1)
facets += 1

abaq_writer = []

abaq = '*HEADING\n'
abaq += 'cubit(t/neptun-dunkles/axisem/kerner/unit_tests/tetrahedron.inp):02' \
        '/17/2014: 17\n'
abaq += 'version: 13.0\n'
abaq += '**\n'
abaq += '********************************** P A R T S ' \
        '**********************************\n'
abaq += '*PART, NAME=Part-Default\n'
abaq += '**\n'
abaq += '********************************** N O D E S ' \
        '**********************************\n'
abaq += '*NODE, NSET=ALLNODES\n'

abaq_writer.append(abaq)

for i in range(np.shape(vertices)[0]):
    abaq_writer.append('%s, %s, %s, %s\n' % (i+1,
                                             vertices[i, 0],
                                             vertices[i, 1],
                                             vertices[i, 2]))
abaq = '**\n'
abaq += '********************************** E L E M E N T S ' \
        '****************************\n'
abaq += '*ELEMENT, TYPE=C3D4, ELSET=EB2\n'

abaq_writer.append(abaq)

for i in range(np.shape(facets)[0]):
    abaq_writer.append('%s, %s, %s, %s, %s\n' % (i+1,
                                                 int(facets[i, 0]),
                                                 int(facets[i, 1]),
                                                 int(facets[i, 2]),
                                                 int(facets[i, 3])))

abaq = '**\n'
abaq += '********************************** P R O P E R T I E S ' \
        '************************\n'
abaq += '*SOLID SECTION, ELSET=EB2, MATERIAL=Default-Steel\n'
abaq += '**\n'
abaq += '*END PART\n'

abaq_writer.append(abaq)

writer_fio = open(os.path.join(par_add, '%s.inp' % tet_name), 'w')
for line in abaq_writer:
    writer_fio.writelines(line)
writer_fio.close()

