#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  run_resolutiontest.py
#   Purpose:   Compile and run resolutiontest.f
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------
# Required Python modules will be imported in this part.
import datetime
import os
import sys
import os
import glob
import shutil
import time

# -----------------------------------------------------------------------
# ----------------------------- INPUT -----------------------------------
# -----------------------------------------------------------------------

# path to assemblematrix folder
path_assemble = '/mnt/seismodata2/MT/SEA-SEIS/assmat/ISC_29.05.2020'
# path to solvetome
path_solvetomo = '/mnt/seismodata2/MT/SEA-SEIS/solvetomo/ISC_29.05.2020'
# path to resolution test folder
path_resolution = '/mnt/seismodata2/MT/SEA-SEIS/resolution/ISC_29.05.2020/res_op2_n0'

# vertices/facets names, assuming they are in assemble_dir
vertices = 'vertices.h4'
facets = 'facets.h4'

# indent
indent = 'iscP'
# new indent for test file
new_indent = 'res_iscP'

# XXX  create in.resolution ??
input_res = True
input_name = 'in.res_op2_n0'

# -----------------------------------------------------------------------
# ------------------------------- COPY ----------------------------------
# -----------------------------------------------------------------------

tic = time.time()
cur_dir = os.path.abspath(os.curdir)

# create resolution folder
if not os.path.isdir(path_resolution):
    os.mkdir(path_resolution)

# XXX for now copy all necessary files to resolution folder saver that way, even if
# it takes around five to ten minutes

# assemblematrix.events.<INDENT>
path_2_assemblematrix = glob.glob(os.path.join(path_assemble, 'assemblematrix.*.%s' % indent))
for assemble in path_2_assemblematrix:
    bsname = os.path.basename(assemble)
    if os.path.isfile(os.path.join(path_resolution, bsname)):
        print bsname, 'is already in destination folder'
    elif os.path.isfile(assemble):
        shutil.copyfile(assemble, os.path.join(path_resolution, bsname))
        print 'coping assemblematrix.*.%s \n' % indent
    else:
        print bsname, 'cannot be found. EXIT!'
        import ipdb;
        ipdb.set_trace()
        sys.exit('ERROR')

# ata.1.<INDENT>

path_2_ata = os.path.join(path_assemble, 'ata.1.%s' % indent)
ata_bsname = os.path.basename(path_2_ata)
if os.path.isfile(os.path.join(path_resolution, ata_bsname)):
    print ata_bsname, 'is already in destination folder'
elif os.path.isfile(path_2_ata):
    shutil.copyfile(path_2_ata, os.path.join(path_resolution, ata_bsname))
    print 'coping ata.1.%s \n' % indent
else:
    print ata_bsname, 'cannot be found. EXIT!'
    sys.exit('ERROR')

# columndensity.1.<INDENT>

path_2_columndensity = os.path.join(path_assemble, 'columndensity.1.%s' % indent)
columndensity_bsname = os.path.basename(path_2_columndensity)
if os.path.isfile(os.path.join(path_resolution, columndensity_bsname)):
    print columndensity_bsname, 'is already in destination folder'
elif os.path.isfile(path_2_columndensity):
    shutil.copyfile(path_2_columndensity, os.path.join(path_resolution, columndensity_bsname))
    print 'coping columndensity.1.%s \n' % indent
else:
    print columndensity_bsname, 'cannot be found. EXIT!'
    import ipdb;
    ipdb.set_trace()
    sys.exit('ERROR')
    
# colindx.<INDENT>

path_2_colindx = os.path.join(path_assemble, 'colindx.%s' % indent)
colindx_bsname = os.path.basename(path_2_colindx)
if os.path.isfile(os.path.join(path_resolution, colindx_bsname)):
    print colindx_bsname, 'is already in destination folder'
elif os.path.isfile(path_2_colindx):
    shutil.copyfile(path_2_colindx, os.path.join(path_resolution, colindx_bsname))
    print 'coping colindx.%s \n' % indent
else:
    print colindx_bsname, 'cannot be found. EXIT!'
    import ipdb;
    ipdb.set_trace()
    sys.exit('ERROR')

# data.<INDENT>

path_2_data = os.path.join(path_assemble, 'data.%s' % indent)
data_bsname = os.path.basename(path_2_data)
if os.path.isfile(os.path.join(path_resolution, data_bsname)):
    print data_bsname, 'is already in destination folder'
elif os.path.isfile(path_2_data):
    shutil.copyfile(path_2_data, os.path.join(path_resolution, data_bsname))
    print 'coping data.%s \n' % indent
else:
    print data_bsname, 'cannot be found. EXIT!'
    import ipdb;
    ipdb.set_trace()
    sys.exit('ERROR')

# vertices

path_2_vertices = os.path.join(path_assemble, vertices)
if os.path.isfile(os.path.join(path_resolution, vertices)):
    print vertices, 'is already in destination folder'
elif os.path.isfile(path_2_vertices):
    shutil.copyfile(path_2_vertices, os.path.join(path_resolution, vertices))
    print 'coping', vertices
else:
    print vertices, 'cannot be found. EXIT!'
    import ipdb;

    ipdb.set_trace()
    sys.exit('ERROR')

# facets

path_2_facets = os.path.join(path_assemble, facets)
if os.path.isfile(os.path.join(path_resolution, facets)):
    print facets, 'is already in destination folder'
elif os.path.isfile(path_2_facets):
    shutil.copyfile(path_2_facets, os.path.join(path_resolution, facets))
    print 'coping', facets
else:
    print facets, 'cannot be found. EXIT!'
    import ipdb;

    ipdb.set_trace()
    sys.exit('ERROR')

# mat.<INDENT>

path_2_mat= os.path.join(path_solvetomo, 'mat.%s' % indent)
matbsname = os.path.basename(path_2_mat)
if os.path.isfile(os.path.join(path_resolution, matbsname)):
    print matbsname, 'is already in destination folder'
elif os.path.isfile(path_2_mat):
    shutil.copyfile(path_2_mat, os.path.join(path_resolution, matbsname))
    print 'coping mat.%s \n' % indent
else:
    print matbsname, 'cannot be found. EXIT!'
    import ipdb;

    ipdb.set_trace()
    sys.exit('ERROR')

# aux.<INDENT>

path_2_aux = os.path.join(path_solvetomo, 'aux.%s' % indent)
aux_bsname = os.path.basename(path_2_aux)
if os.path.isfile(os.path.join(path_resolution, aux_bsname)):
    print aux_bsname, 'is already in destination folder'
elif os.path.isfile(path_2_aux):
    shutil.copyfile(path_2_aux, os.path.join(path_resolution, aux_bsname))
    print 'coping aux.%s \n' % indent
else:
    print aux_bsname, 'cannot be found. EXIT!'
    import ipdb;

    ipdb.set_trace()
    sys.exit('ERROR')

# atamx.<INDENT>

path_2_atamx = os.path.join(path_solvetomo, 'atamx.%s' % indent)
atamx_bsname = os.path.basename(path_2_atamx)
if os.path.isfile(os.path.join(path_resolution, atamx_bsname)):
    print atamx_bsname, 'is already in destination folder'
elif os.path.isfile(path_2_atamx):
    shutil.copyfile(path_2_atamx, os.path.join(path_resolution, atamx_bsname))
    print 'coping atamx.%s \n' % indent
else:
    print atamx_bsname, 'cannot be found. EXIT!'
    import ipdb;

    ipdb.set_trace()
    sys.exit('ERROR')

toc = time.time()

print 'Time lapsed: %s min' % (round((toc-tic)/60., 2))


# -----------------------------------------------------------------------
# ------------------------------- MAIN ----------------------------------
# -----------------------------------------------------------------------

outfile_fio = open(os.path.join(path_resolution, 'in.link'), 'w+')
outfile_fio.writelines('%s\n' % indent)
outfile_fio.writelines('%s' % new_indent)
outfile_fio.close()

try:
    os_sys = os.system('./make')
    shutil.move('resolutiontest', os.path.join(path_resolution, 'resolutiontest'))
    shutil.copyfile('run_resolution.sh', os.path.join(path_resolution, 'run_resolution.sh'))
    os.chdir(os.path.join(path_resolution))
    print '\n------------------------'
    print '- Run RESOLUTIONTEST.F -'
    print '------------------------\n'
    if input_res:
        os.system('./resolutiontest < %s' % input_name)
    else:
        os.system('./resolutiontest')
    print '\n-----------------'
    print '- Linking FILES -'
    print '-----------------\n'

    try:
        os.system('chmod +x run_resolution.sh')
    except:
        pass
    os.system('./run_resolution.sh < in.link')

except Exception, exp:
    print exp
    import ipdb; ipdb.set_trace()
    sys.exit('Automatic EXIT!')


