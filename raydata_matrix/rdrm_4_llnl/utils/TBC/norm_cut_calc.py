#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  norm_cut_calc.py
#   Purpose:   Return norm vectors (to be used in Paraview)
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

############################### shoot #########################################


def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq
    """
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.
 
    EPS= 0.00000000005
    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
        alert("Only N-S courses are meaningful, starting at a pole!")
 
    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2. * np.arctan2 (tu, cf)
 
    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs (y - c) > EPS):
 
        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
              d / 4. - cz) * sy * d + tu
 
    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi    
 
    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)
 
    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi
 
    return (glon2, glat2, baz)
 
############################### great #########################################


def great(m, startlon, startlat, azimuth,*args, **kwargs):
    glon1 = startlon
    glat1 = startlat
    glon2 = glon1
    glat2 = glat1
 
    step = 50
 
    glon2, glat2, baz = shoot(glon1, glat1, azimuth, step)
    if azimuth-180 >= 0:
        while glon2 <= startlon:
            m.drawgreatcircle(glon1, glat1, glon2, glat2,del_s=50,**kwargs)
            azimuth = baz + 180.
            glat1, glon1 = (glat2, glon2)
 
            glon2, glat2, baz = shoot(glon1, glat1, azimuth, step)
    elif azimuth-180 < 0:
        while glon2 >= startlon:
            m.drawgreatcircle(glon1, glat1, glon2, glat2,del_s=50,**kwargs)
            azimuth = baz + 180.
            glat1, glon1 = (glat2, glon2)
 
            glon2, glat2, baz = shoot(glon1, glat1, azimuth, step)
 
############################### vect_gen #########################################


def vect_gen(lat1, lon1, lat2, lon2):
    """
    INPUT: latitudes and longitudes of two points
    OUTPUT: two vectors to be used for np.cross
    """
    rho = 1.0
    lats = np.array([lat1, lat2])
    lons = np.array([lon1, lon2])

    theta = (90.-lats)*(np.pi/180)
    phi = lons*(np.pi/180)

    zs = rho*np.cos(theta)
    xs = rho*np.sin(theta)*np.cos(phi)
    ys = rho*np.sin(theta)*np.sin(phi)
    vec1 = [xs[0], ys[0], zs[0]]
    vec2 = [xs[1], ys[1], zs[1]]
    return vec1, vec2

############################### norm_cut_calc #########################################


def norm_cut_calc(lat1, lon1, lat2, lon2):
    """
    Use google earth to have latitudes and longitudes of the required line.
    """
    plt.ion()
    m = Basemap(projection='ortho', lon_0=(lon1+lon2)/2., lat_0=(lat1+lat2)/2, resolution='c')
    m.fillcontinents()
    m.drawparallels(np.arange(-90., 120., 30.))
    m.drawmeridians(np.arange(0., 420., 60.))
    m.drawmapboundary()

    m.drawgreatcircle(lon1, lat1, lon2, lat2, 1000, linewidth=3, c='r')

    vec1, vec2 = vect_gen(lat1, lon1, lat2, lon2)
    vec_cross = np.cross(vec1, vec2)
    print '------------------------'
    print 'Requested cross-section:'
    print vec_cross/np.sqrt(np.sum(np.power(vec_cross, 2)))

    lat2_2 = lat1 + 0.01
    vec1, vec2 = vect_gen(lat1, lon1, lat2_2, lon1)
    vec_cross = np.cross(vec1, vec2)
    print '------------------------'
    print 'Same longitude as the first point:'
    print vec_cross/np.sqrt(np.sum(np.power(vec_cross, 2)))

    lon2_3 = lon1 + 0.01
    vec1, vec2 = vect_gen(lat1, lon1, lat1, lon2_3)
    vec_cross = np.cross(vec1, vec2)
    print '------------------------'
    print 'Same latitude as the first point:'
    print vec_cross/np.sqrt(np.sum(np.power(vec_cross, 2)))

    raw_input('Press Enter to finish the program...')


if True:
    latlons = raw_input('Enter lat1,lon1,lat2,lon2\n')
    lat1, lon1, lat2, lon2 = latlons.split(',')
    if float(lon1) < 0:
        lon1 = 360 - abs(float(lon1))
    if float(lon2) < 0:
        lon2 = 360 - abs(float(lon2))
    norm_cut_calc(float(lat1), float(lon1), float(lat2), float(lon2))

else:
    glon1 = 0.
    glat1 = 0.
    azimuth = 30.
    maxdist = 9000

    glon2, glat2, baz = shoot(glon1, glat1, azimuth, maxdist)
    print '%s,%s,%s,%s' % (glat1, glon1, glat2, glon2)

