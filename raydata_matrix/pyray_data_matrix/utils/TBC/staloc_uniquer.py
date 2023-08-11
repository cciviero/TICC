"""
Simple script to generate one single file out of all location files
"""
import glob
import numpy as np
import os

dir_add = '../OUTPUTS/USA10/station_ev'

####################### unique_rows #############################


def unique_rows(data):
    """
    Make a unique array based on the array's row
    """
    data_mid = np.ascontiguousarray(data).view(np.dtype((np.void, data.dtype.itemsize * data.shape[1])))
    _, idx = np.unique(data_mid, return_index=True)
    unique_data = data[idx]
    return unique_data


####################### MAIN PROGRAM #############################

sta_add = glob.glob(os.path.join(dir_add, '*'))
for i in range(len(sta_add)):
    st_tmp = np.loadtxt(sta_add[i])
    if i == 0:
        st_all = st_tmp
    else:
        st_all = np.append(st_all, st_tmp, 0)
st_unique = unique_rows(st_all)
sta_fio = open('./ALL_staloc.txt', 'w')
for i in range(len(st_unique)):
    sta_fio.writelines('%s   %s   %s\n' % (st_unique[i, 0], st_unique[i, 1], st_unique[i, 2]))
sta_fio.close
