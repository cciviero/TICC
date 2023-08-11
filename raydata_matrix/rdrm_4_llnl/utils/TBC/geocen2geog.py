import math as m_math
import matplotlib.pyplot as plt
import numpy as np

############## geocen2geog_single ###################


def geocen2geog_single(lat_gc):
    fac = 0.993305621334896

    colat_gc = 90.0 - lat_gc
    if colat_gc == 0.:
        colat_gc = 1.0e-5

    colat_geo = m_math.atan(
        fac*m_math.cos(m_math.pi/2. - colat_gc*m_math.pi/180.)/
        (m_math.copysign(1, m_math.sin(m_math.pi/2. -
                                       colat_gc*m_math.pi/180.))*
         max(1.0e-30,
             abs(m_math.sin(m_math.pi/2. - colat_gc*m_math.pi/180.)))))
    colat_geo = colat_geo*180./m_math.pi
    lat_geo_final = m_math.copysign(1, colat_geo)*(90 - abs(colat_geo))
    return lat_geo_final

############## geog2geocen_single ###################


def geog2geocen_single(lat_geo):
    fac = 0.993305621334896

    colat_geo = 90.0 - lat_geo
    if colat_geo == 0.:
        colat_geo = 1.0e-5

    colat_gc = m_math.pi/2. - \
                 m_math.atan(fac*m_math.cos(colat_geo*m_math.pi/180.)/
                             max(1.0e-30,
                                 m_math.sin(colat_geo*m_math.pi/180.)))
    colat_gc = colat_gc*180./m_math.pi
    lat_geo_final = 90.0 - colat_gc

    return lat_geo_final

############## geocen2geog ###################


def geocen2geog(evlagc, evlogc, stlagc, stlogc):
    fac = 0.993305621334896

    # ----------------- First for station:
    colatgc = 90.0 - stlagc
    if colatgc == 0.:
        colatgc = 1.0e-5

    # Formula:
    # arg = colat*rpd
    # geocen=pi2-atan(fac*cos(arg)/(max(1.0e-30,sin(arg))))
    geog_lat = m_math.atan(
        fac*m_math.cos(m_math.pi/2. - colatgc*m_math.pi/180.)/
        (m_math.copysign(1, m_math.sin(m_math.pi/2. -
                                       colatgc*m_math.pi/180.))*
         max(1.0e-30, abs(m_math.sin(m_math.pi/2. - colatgc*m_math.pi/180.)))))
    stcola_geog = geog_lat*180./m_math.pi
    stla_geog = m_math.copysign(1, stcola_geog)*(90 - abs(stcola_geog))

    # ----------------- Second for event:
    colatgc = 90.0 - evlagc
    if colatgc == 0.:
        colatgc = 1.0e-5

    geog_lat = m_math.atan(
        fac*m_math.cos(m_math.pi/2. - colatgc*m_math.pi/180.)/
        (m_math.copysign(1, m_math.sin(m_math.pi/2. -
                                       colatgc*m_math.pi/180.))*
         max(1.0e-30, abs(m_math.sin(m_math.pi/2. - colatgc*m_math.pi/180.)))))
    evcola_geog = geog_lat*180./m_math.pi
    evla_geog = m_math.copysign(1, evcola_geog)*(90 - abs(evcola_geog))

    return evla_geog, stla_geog

############## geog2geocen ###################


def geog2geocen(evla, evlo, stla, stlo):
    fac = 0.993305621334896

    # ----------------- First for station:
    colat = 90.0 - stla
    if colat == 0.:
        colat = 1.0e-5

    # Formula:
    # arg = colat*rpd
    # geocen=pi2-atan(fac*cos(arg)/(max(1.0e-30,sin(arg))))
    geocen_lat = m_math.pi/2. - \
                 m_math.atan(fac*m_math.cos(colat*m_math.pi/180.)/
                             max(1.0e-30, m_math.sin(colat*m_math.pi/180.)))
    stcola_geocen = geocen_lat*180./m_math.pi
    stla_geocen = 90.0 - stcola_geocen

    # ----------------- Second for event:
    colat = 90.0 - evla
    if colat == 0.:
        colat = 1.0e-5

    geocen_lat = m_math.pi/2. - \
                 m_math.atan(fac*m_math.cos(colat*m_math.pi/180.)/
                             max(1.0e-30, m_math.sin(colat*m_math.pi/180.)))
    evcola_geocen = geocen_lat*180./m_math.pi
    evla_geocen = 90.0 - evcola_geocen

    return evla_geocen, stla_geocen




############## TESTING ###################

import matplotlib.pyplot as plt

for i in np.arange(-90, 90, 0.1):
    gc = geog2geocen_single(i)
    plt.subplot(2, 1, 1)
    plt.plot(i, gc - i, 'r.')
    plt.title('difference between geocentric - geographic')
    gg = geocen2geog_single(gc)
    plt.subplot(2, 1, 2)
    plt.plot(i, gg - i, 'b.')
    plt.title('difference between geographics (retrieved - input)')
plt.tight_layout()
plt.show()

