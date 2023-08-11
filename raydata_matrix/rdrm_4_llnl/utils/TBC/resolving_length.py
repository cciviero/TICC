"""
Plot the resolving length of each cell of inversion grids
"""

import matplotlib.pyplot as plt
import numpy as np
from pylab import cm as pltcm
import sys

######################## INPUT
# case=1: 200_660
# case=2: 900_300
# case=3: staev_200
case = [1, 2, 3]
man_vmax = 3.2
gridsize = 100
earth_radius = 6371
########################

plt.ion()

for icase in range(len(case)):
    plt.figure()

    if case[icase] == 1:
        address = '../OUTPUTS/P_Pdiff/INVERSION/200_660/' \
                  'mesh.resolvinglength.ppdiff_200_660'
        plt_title = 'hmax: 400km\n' \
                    '900   200km\n' \
                    '2889  400km\n' \
                    'r:1400km 66 --> 200km'
        vl1 = earth_radius - 900
        vl2 = earth_radius - 2889
        vl3 = earth_radius - 1400
    elif case[icase] == 2:
        address = '../OUTPUTS/P_Pdiff/INVERSION/900_300/' \
                  'mesh.resolvinglength.ppdiff_900_300'
        plt_title = 'hmax: 300km\n' \
                    '900   200km\n' \
                    '2889  300km\n' \
                    'r:1400km 66 --> 200km'
        vl1 = earth_radius - 900
        vl2 = earth_radius - 2889
        vl3 = earth_radius - 1400
    elif case[icase] == 3:
        address = '../OUTPUTS/P_Pdiff/INVERSION/staev_200/' \
                  'mesh.resolvinglength.ppdiff_staev_200'
        plt_title = 'hmax: 300km\n' \
                    '2889 200\n' \
                    '4000 250\n' \
                    'r:1400km 66 --> 200km'
        vl1 = earth_radius - 2889
        vl2 = earth_radius - 4000
        vl3 = earth_radius - 1400
    else:
        sys.exit('This case %s has not implemented!' % icase)


    cmap = pltcm.get_cmap('afmhot_r', 18)
    mesh_res = np.loadtxt(address, skiprows=2)
    radius_mesh = np.sqrt(np.power(mesh_res[:, 0], 2) +
                          np.power(mesh_res[:, 1], 2) +
                          np.power(mesh_res[:, 2], 2))

    plt.hexbin(radius_mesh, mesh_res[:, 3], bins='log',  gridsize=gridsize,
               cmap=cmap, vmax=man_vmax)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=24)

    plt.axvline(x=vl1, color='k', linestyle='--', lw=3)
    plt.text(vl1+30, 410, '%skm' % (earth_radius-vl1), size=18, weight='bold')
    plt.axvline(x=vl2, color='k', linestyle='--', lw=3)
    plt.text(vl2+30, 420, '%skm' % (earth_radius-vl2), size=18, weight='bold')
    plt.axvline(x=vl3, color='k', linestyle='--', lw=3)
    plt.text(vl3+30, 430, '%skm' % (earth_radius-vl3), size=18, weight='bold')
    plt.axhline(y=66, color='k', linestyle='--', lw=3)
    plt.ylim(50, 450)
    plt.xticks(size=24, weight='bold')
    plt.yticks(size=24, weight='bold')
    plt.xlabel('Radius (km)', size=24, weight='bold')
    plt.ylabel('Edge length (km)', size=24, weight='bold')
    plt.title(plt_title)
    #xy = np.vstack([radius_mesh, mesh_res[:,3]])
    #z_dense = gaussian_kde(xy)(xy)
    #plt.scatter(radius_mesh, mesh_res, c=z_dense, cmap=cmap, edgecolor='none')
    #plt.plot(radius_mesh, mesh_res[:,3], '.')
plt.show()
raw_input('Press enter to exit the program...')
