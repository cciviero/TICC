"""
This program tried to plot a figure with 4 subplots in which
the intensity or in other words the values should be plotted
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.basemap import Basemap
import numpy as np
import pickle

fio_pkl = open('final_pickled.pkl', 'r')
pkl_obj = pickle.load(fio_pkl)

plt.ion()
plt.figure(figsize=(20, 10))
plt.clf()

plt.subplot(2, 2, 1)
m = Basemap(projection='cyl', lon_0=180.0, lat_0=0.0, resolution='c')
m.drawcoastlines()
m.drawmapboundary()
m.imshow(pkl_obj[:, :, 0],
         interpolation='none',
         norm=LogNorm(vmin=1, vmax=np.max(pkl_obj[:, :, 0])))
#         cmap='hot_r')
cbar = plt.colorbar(orientation='horizontal', shrink=0.8)
cbar.ax.tick_params(labelsize=18)
plt.title('0-90 azimuth\n#%i' % np.sum(pkl_obj[:, :, 0]),
          size=18, weight='bold')

print 'Plot first figure: 2'
plt.subplot(2, 2, 2)
m = Basemap(projection='cyl', lon_0=180.0, lat_0=0.0, resolution='c')
m.drawcoastlines()
m.drawmapboundary()
m.imshow(pkl_obj[:, :, 1],
         interpolation='none',
         norm=LogNorm(vmin=1, vmax=np.max(pkl_obj[:, :, 1])))
cbar = plt.colorbar(orientation='horizontal', shrink=0.8)
cbar.ax.tick_params(labelsize=18)
plt.title('90-180 azimuth\n#%i' % np.sum(pkl_obj[:, :, 1]),
          size=18, weight='bold')

print 'Plot first figure: 3'
plt.subplot(2, 2, 3)
m = Basemap(projection='cyl', lon_0=180.0, lat_0=0.0, resolution='c')
m.drawcoastlines()
m.drawmapboundary()
m.imshow(pkl_obj[:, :, 2],
         interpolation='none',
         norm=LogNorm(vmin=1, vmax=np.max(pkl_obj[:, :, 2])))
#         cmap='Reds')
cbar = plt.colorbar(orientation='horizontal', shrink=0.8)
cbar.ax.tick_params(labelsize=18)
plt.title('180-270 azimuth\n#%i' % np.sum(pkl_obj[:, :, 2]),
          size=18, weight='bold')

print 'Plot first figure: 4'
plt.subplot(2, 2, 4)
m = Basemap(projection='cyl', lon_0=180.0, lat_0=0.0, resolution='c')
m.drawcoastlines()
m.drawmapboundary()
m.imshow(pkl_obj[:, :, 3],
         interpolation='none',
         norm=LogNorm(vmin=1, vmax=np.max(pkl_obj[:, :, 3])))
#         cmap='Reds')
cbar = plt.colorbar(orientation='horizontal', shrink=0.8)
cbar.ax.tick_params(labelsize=18)
plt.title('270-360 azimuth\n#%i' % np.sum(pkl_obj[:, :, 3]),
          size=18, weight='bold')


raw_input('press enter...')

# ######################### 3D...not sure whether it is useful
# from mpl_toolkits.mplot3d.axes3d import Axes3D
# import numpy as np
# fig = plt.figure()
# plt.ion()
# ax = Axes3D(fig)
# ys = np.linspace(0, 360, np.shape(pkl_obj)[0])
# xs = np.linspace(0, 180, np.shape(pkl_obj)[1])
# xx, yy = np.meshgrid(xs, ys)
#
# fig = plt.figure()
# ax = Axes3D(fig)
# ax.plot_surface(xx, yy, pkl_obj[:, :, 0], rstride=1, cstride=1, cmap='hot')
# plt.show()
#
# fig = plt.figure()
# ax = Axes3D(fig)
# ax.plot_surface(xx, yy, pkl_obj[:, :, 1], rstride=1, cstride=1, cmap='hot')
# plt.show()
#
# fig = plt.figure()
# ax = Axes3D(fig)
# ax.plot_surface(xx, yy, pkl_obj[:, :, 2], rstride=1, cstride=1, cmap='hot')
# plt.show()
#
# fig = plt.figure()
# ax = Axes3D(fig)
# ax.plot_surface(xx, yy, pkl_obj[:, :, 3], rstride=1, cstride=1, cmap='hot')
# plt.show()
