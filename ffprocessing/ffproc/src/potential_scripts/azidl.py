#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------
#   Filename:  azidl.py
#   Purpose:   Azimuth and Distance calculation on an ellipsoidal Earth
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
#
#   Copyright (C) 2014 Kasra Hosseini
#-------------------------------------------------------------------

#-----------------------------------------------------------------------
#----------------Import required Modules (Python and Obspy)-------------
#-----------------------------------------------------------------------

# Required Python modules will be imported in this part.

# Added this line for python 2.5 compatibility
from __future__ import with_statement
import numpy as np

# INFO:
# Header of the original FORTRAN code: (wkbj_solo.f90 in our archive)
#   SUBROUTINE AZIDL(SLAT,SLON,ELAT,ELON,DEL,DELS,AZIS,AZIE)
#   AZIS IS AZIMUTH UIT STATION, AZIE UIT EPICENTRUM
#   DELS IS THE DISTANCE OVER THE SURFACE OF AN ELLIPSOIDAL EARTH
#   CALCULATED WITH FORM 29 OF THOMAS, JGR 70, PP 3331-3340, 1965
#   ALLES IN RADIALEN
#   rclat, rclon, slat, slon must be given in radians, e.g.
#   rclat=stla(i)*pi2/90.

# ----------------- FUNCTIONS ---------------------
def COLAT(X):
    return 1.57079632-np.arctan([.993277*np.tan(X)])

def AZIDL(SLAT, SLON, ELAT, ELON):
    SLAT = float(SLAT)*np.pi/180.
    SLON = float(SLON)*np.pi/180.
    ELAT = float(ELAT)*np.pi/180.
    ELON = float(ELON)*np.pi/180.
    SCOLAT=COLAT(SLAT)
    ECOLAT=COLAT(ELAT)
    SC1=np.sin(SCOLAT)
    SC2=np.cos(SCOLAT)
    SC3=np.sin(SLON)
    SC4=np.cos(SLON)
    EC1=np.sin(ECOLAT)
    EC2=np.cos(ECOLAT)
    EC3=np.sin(ELON)
    EC4=np.cos(ELON)
    AE=EC1*EC4
    BE=EC1*EC3
    AZI1=(AE-SC3)**2+(BE+SC4)**2+EC2*EC2-2.
    AZI2=(AE-SC2*SC4)**2+(BE-SC2*SC3)**2+(EC2+SC1)**2-2.

    if AZI2 == 0.:
        AZIS=np.pi-np.sign([AZI1])*np.pi/2.
    else:
        AZIS=np.arctan(complex(AZI2, AZI1))
        CODELB=SC1*(AE*SC4+BE*SC3)+SC2*EC2

    DEL=np.arccos(CODELB)
    AS=SC1*SC4
    BS=SC1*SC3
    AZI1=(AS-EC3)**2+(BS+EC4)**2+SC2*SC2-2.
    AZI2=(AS-EC2*EC4)**2+(BS-EC2*EC3)**2+(SC2+EC1)**2-2.

    if AZI2 == 0.:
        AZIE=np.pi-np.sign([AZI1])*np.pi/2.
    else:
        AZIE=np.arctan(complex(AZI2, AZI1))

    SCOLAT=np.pi/2.-SLAT
    ECOLAT=np.pi/2.-ELAT
    SC1=np.sin(SCOLAT)
    SC2=np.cos(SCOLAT)
    EC1=np.sin(ECOLAT)
    EC2=np.cos(ECOLAT)
    AE=EC1*EC4
    BE=EC1*EC3
    DELBOD=np.arccos(SC1*(AE*SC4+BE*SC3)+SC2*EC2)
    COSDEL=np.cos(DELBOD)
    X1=((SC2+EC2)**2)/(1.+COSDEL)
    X2=((SC2-EC2)**2)/(1.-COSDEL)
    X=X1+X2
    Y=X1-X2
    COTDEL=1./np.tan(DELBOD)
    SINDEL=np.sin(DELBOD)
    DEL2=DELBOD*DELBOD
    A=64.*DELBOD+16.*DEL2*COTDEL
    D=48.*SINDEL+8.*DEL2/COSDEL
    B=-(D+D)
    E=30.*np.sin(DELBOD+DELBOD)
    C=-30.*DELBOD-8.*DEL2*COTDEL-E/2.
    DELS=6378.2064*(DELBOD-.000847518825*(X*DELBOD-3.*Y*SINDEL)+.0897860195E-6*(X*(A+C*X+D*Y)+Y*(B+E*Y)))

    DEL_degree = DEL[0]*180./np.pi
    return DEL_degree
