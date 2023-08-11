#!/usr/bin/env python
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------
#   Filename:  paraview_solx_vtk.py
#   Purpose:   Convenient tool for working with and visualizing VTK files
#              This code should be run from Paraview and not separately
#   Author:    Kasra Hosseini
#   Email:     hosseini@geophysik.uni-muenchen.de
#   License:   GPLv3
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ----------------Import required Modules (Python and Obspy)-------------
# -----------------------------------------------------------------------

# Required Python modules will be imported in this part.
import numpy as np
import os

try: 
    paraview.simple
except:
    from paraview.simple import *

# ------------------- INPUT -----------------------------
# VTK file to be analyzed
solx1_file = '/Users/maria/PhD/Codes/Inversion/2016_inversion/homo_400_rhum_rum_P_v0103/inversion_L_0_3/vel_00_solx.05.dlnVp.homo_400.vtk'
continent_file = '/Users/maria/PhD/ParaView_Files/continents.vtk'
plates_file = '/Users/maria/PhD/ParaView_Files/plates.vtk'
hotspot_file = '/Users/maria/PhD/ParaView_Files/hotspots.vtk'

# TDC
latmin = -60
lonmin = -20
latmax = -26
lonmax = 10
# Africa
latmin = -60
lonmin = -20
latmax = 45
lonmax = 70

# Pitcairn
#latmin = -35
#lonmin = -140
#latmax = -15
#lonmax = -100
# NA
#latmin = 16
#lonmin = -123
#latmax = 55
#lonmax = -47

# La Reunion
latmin = -40
lonmin = 35
latmax = 0
lonmax = 70

# Louisville
#latmin = -70
#lonmin = -160
#latmax = -10
#lonmax = -70
# Indian Ocean
#latmin = -50
#lonmin = 80
#latmax = -15
#lonmax = 120

# Unknown!
lat_cut1 = -37.11
lon_cut1 = -12.19
lat_cut2 = -17.07
lon_cut2 = 13.11
# ------------------- END INPUT -----------------------------


# ############################## vect_gen #####################################


def vect_gen(lat1, lon1, lat2, lon2):
    """
    INPUT: latitudes and longitudes of two points
    OUTPUT: two vectors to be used for np.cross
    :param lat1:
    :param lon1:
    :param lat2:
    :param lon2:
    :return:
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

# ############################## norm_cut_calc ################################


def norm_cut_calc(lat1, lon1, lat2, lon2):
    """
    Calculate the vector norms.
    Suggestion:
    Use google earth to have latitudes and longitudes of the required line.
    """
    vec1, vec2 = vect_gen(lat1, lon1, lat2, lon2)
    vec_cross = np.cross(vec1, vec2)
    print '------------------------'
    print 'Requested cross-section:'
    print vec_cross/np.sqrt(np.sum(np.power(vec_cross, 2)))
    return vec_cross/np.sqrt(np.sum(np.power(vec_cross, 2)))

print '\n\n================================================================='
print '!!!!!!!!!!!!!!!!!!!! PROGRAM START !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
print '================================================================='

value_vtk = os.path.basename(solx1_file)
# Calculate the norm vector of the slices
(x1, y1, z1) = norm_cut_calc(latmin, lonmin, latmin, lonmax)
(x2, y2, z2) = norm_cut_calc(latmin, lonmin, latmax, lonmin)
(x3, y3, z3) = norm_cut_calc(latmax, lonmin, latmax, lonmax)
(x4, y4, z4) = norm_cut_calc(latmin, lonmax, latmax, lonmax)
(xs, ys, zs) = norm_cut_calc(lat_cut1, lon_cut1, lat_cut2, lon_cut2)

# Create the colors:
paraview.simple._DisableFirstRenderCameraReset()
solx1_PiecewiseFunction = \
    CreatePiecewiseFunction(Points=[-0.00999999977648258, 1.0, 0.5,
                                    0.0, 0.01, 1.0, 0.5, 0.0])
solx1_PVLookupTable = CreateLookupTable(
    RGBPoints=[-0.01, 0.0157473, 0.00332647, 0.0,
               -0.00882355, 0.299413, 0.000366217, 0.000549325,
               -0.00756864, 0.508812, 0.0, 0.0,
               -0.00631374, 0.672862, 0.139086, 0.00270085,
               -0.0050588, 0.777096, 0.330175, 0.000885023,
               -0.00380392, 0.854063, 0.510857, 0.0,
               -0.00254902, 0.904479, 0.690486, 0.0,
               -0.00129412, 0.9514, 0.835615, 0.449271,
               -3.92199999999995e-05, 0.881376, 0.912184, 0.818097,
               0.0012157, 0.60824, 0.892164, 0.935546,
               0.00247058, 0.327672, 0.784939, 0.873426,
               0.0037255, 0.188952, 0.641306, 0.792096,
               0.00498038, 0.128054, 0.492592, 0.720287,
               0.0062353, 0.0552529, 0.345022, 0.659495,
               0.0074902, 0.0, 0.216587, 0.524575,
               0.0087451, 0.0, 0.120394, 0.302678,
               0.01, 0.0, 0.0, 0.0],
    VectorMode='Magnitude',
    NanColor=[0.498039, 0.0, 0.0],
    ScalarOpacityFunction=solx1_PiecewiseFunction,
    ColorSpace='Lab',
    LockScalarRange=1)

solx1_PVLookupTable = CreateLookupTable(
    VectorMode='Magnitude',
    NanColor=[0.498039, 0.0, 0.0],
    ScalarOpacityFunction=solx1_PiecewiseFunction,
    ColorSpace='Lab',
    LockScalarRange=1,
    RGBPoints=
    [-4.0          ,   0.6980,    0.0941,    0.1686,
     -3.89873417722,   0.6980,    0.0941,    0.1686,
     -3.79746835443,   0.6980,    0.0941,    0.1686,
     -3.69620253165,   0.6980,    0.0941,    0.1686,
     -3.59493670886,   0.6980,    0.0941,    0.1686,
     -3.49367088608,   0.6980,    0.0941,    0.1686,
     -3.39240506329,   0.6980,    0.0941,    0.1686,
     -3.29113924051,   0.6980,    0.0941,    0.1686,
     -3.18987341772,   0.6980,    0.0941,    0.1686,
     -3.08860759494,   0.6980,    0.0941,    0.1686,
     -2.98734177215,   0.6980,    0.0941,    0.1686,
     -2.88607594937,   0.6980,    0.0941,    0.1686,
     -2.78481012658,   0.6980,    0.0941,    0.1686,
     -2.6835443038 ,   0.6980,    0.0941,    0.1686,
     -2.58227848101,   0.6980,    0.0941,    0.1686,
     -2.48101265823,   0.6980,    0.0941,    0.1686,
     -2.37974683544,   0.6980,    0.0941,    0.1686,
     -2.27848101266,   0.6980,    0.0941,    0.1686,
     -2.17721518987,   0.6980,    0.0941,    0.1686,
     -2.07594936709,   0.6980,    0.0941,    0.1686,
     -1.9746835443 ,   0.8392,    0.3765,    0.3020,
     -1.87341772152,   0.8392,    0.3765,    0.3020,
     -1.77215189873,   0.8392,    0.3765,    0.3020,
     -1.67088607595,   0.8392,    0.3765,    0.3020,
     -1.56962025316,   0.8392,    0.3765,    0.3020,
     -1.46835443038,   0.8392,    0.3765,    0.3020,
     -1.36708860759,   0.8392,    0.3765,    0.3020,
     -1.26582278481,   0.8392,    0.3765,    0.3020,
     -1.16455696203,   0.8392,    0.3765,    0.3020,
     -1.06329113924,   0.8392,    0.3765,    0.3020,
     -0.96202531645,   0.9569,    0.6471,    0.5098,
     -0.86075949367,   0.9569,    0.6471,    0.5098,
     -0.75949367088,   0.9569,    0.6471,    0.5098,
     -0.65822784810,   0.9569,    0.6471,    0.5098,
     -0.55696202531,   0.9569,    0.6471,    0.5098,
     -0.45569620253,   0.9922,    0.8588,    0.7804,
     -0.35443037974,   0.9922,    0.8588,    0.7804,
     -0.25316455696,   0.9922,    0.8588,    0.7804,
     -0.15189873417,   0.9686,    0.9686,    0.9686,
     -0.05063291139,   0.9686,    0.9686,    0.9686,
     0.050632911392,   0.9686,    0.9686,    0.9686,
     0.151898734177,   0.9686,    0.9686,    0.9686,
     0.253164556962,   0.9686,    0.9686,    0.9686,
     0.354430379747,   0.8196,    0.8980,    0.9412,
     0.455696202532,   0.8196,    0.8980,    0.9412,
     0.556962025316,   0.5725,    0.7725,    0.8706,
     0.658227848101,   0.5725,    0.7725,    0.8706,
     0.759493670886,   0.5725,    0.7725,    0.8706,
     0.860759493671,   0.5725,    0.7725,    0.8706,
     0.962025316456,   0.5725,    0.7725,    0.8706,
     1.06329113924 ,   0.2627,    0.5765,    0.7647,
     1.16455696203 ,   0.2627,    0.5765,    0.7647,
     1.26582278481 ,   0.2627,    0.5765,    0.7647,
     1.36708860759 ,   0.2627,    0.5765,    0.7647,
     1.46835443038 ,   0.2627,    0.5765,    0.7647,
     1.56962025316 ,   0.2627,    0.5765,    0.7647,
     1.67088607595 ,   0.2627,    0.5765,    0.7647,
     1.77215189873 ,   0.2627,    0.5765,    0.7647,
     1.87341772152 ,   0.2627,    0.5765,    0.7647,
     1.9746835443  ,   0.2627,    0.5765,    0.7647,
     2.07594936709 ,   0.1294,    0.4000,    0.6745,
     2.17721518987 ,   0.1294,    0.4000,    0.6745,
     2.27848101266 ,   0.1294,    0.4000,    0.6745,
     2.37974683544 ,   0.1294,    0.4000,    0.6745,
     2.48101265823 ,   0.1294,    0.4000,    0.6745,
     2.58227848101 ,   0.1294,    0.4000,    0.6745,
     2.6835443038  ,   0.1294,    0.4000,    0.6745,
     2.78481012658 ,   0.1294,    0.4000,    0.6745,
     2.88607594937 ,   0.1294,    0.4000,    0.6745,
     2.98734177215 ,   0.1294,    0.4000,    0.6745,
     3.08860759494 ,   0.1294,    0.4000,    0.6745,
     3.18987341772 ,   0.1294,    0.4000,    0.6745,
     3.29113924051 ,   0.1294,    0.4000,    0.6745,
     3.39240506329 ,   0.1294,    0.4000,    0.6745,
     3.49367088608 ,   0.1294,    0.4000,    0.6745,
     3.59493670886 ,   0.1294,    0.4000,    0.6745,
     3.69620253165 ,   0.1294,    0.4000,    0.6745,
     3.79746835443 ,   0.1294,    0.4000,    0.6745,
     3.89873417722 ,   0.1294,    0.4000,    0.6745,
     4.0           ,   0.1294,    0.4000,    0.6745])

# Do the actual clipping
solx1_vtk = LegacyVTKReader(guiName="solx1", FileNames=[solx1_file])
Clip1 = Clip(guiName="c1", Scalars=['POINTS', value_vtk],
             ClipType="Plane", Value=0.0004882216453552246)
Clip1.ClipType.Normal = [x1, y1, z1]

Clip2 = Clip(guiName="c2", InsideOut=1, Scalars=['POINTS', 'c1'],
             ClipType="Plane", Value=1)
Clip2.ClipType.Normal = [x2, y2, z2]

Clip3 = Clip(guiName="c3", InsideOut=1, Scalars=['POINTS', 'c2'],
             ClipType="Plane", Value=1)
Clip3.ClipType.Normal = [x3, y3, z3]

Clip4 = Clip(guiName="c4", Scalars=['POINTS', 'c3'],
             ClipType="Plane", Value=0.0007011685520410538)
Clip4.ClipType.Normal = [x4, y4, z4]

SetActiveSource(solx1_vtk)
Clip5 = Clip(guiName="cmb", InsideOut=1, Scalars=['POINTS', value_vtk],
             Value=0.0004882216453552246, ClipType="Sphere")
Clip5.ClipType.Center = [0.0, 0.0, 0.0]
Clip5.ClipType.Radius = 3482.0

# For the clipped volume, create some horizontal slices too
for dept in [0, 100, 200, 410, 660, 900, 1100, 1300, 1500, 1800,
             2000, 2200, 2400, 2600, 2800]:
    SetActiveSource(Clip4)
    Slice1 = Slice(guiName="S-%s" % dept, SliceOffsetValues=[0.0],
                   SliceType="Sphere")
    Slice1.SliceType.Radius = 6371 - dept
    Slice1.SliceType.Center = [0.0, 0.0, 0.0]
    SetActiveSource(Slice1)
    DataRepresentation3 = Show()

# Slice (based on lat_cut1, lon_cut1, lat_cut2, lon_cut2)
SetActiveSource(solx1_vtk)
Slice2 = Slice(guiName="Slice2", SliceOffsetValues=[0.0], SliceType="Plane")
Slice2.SliceType.Origin = [0, 0, 0]
Slice2.SliceType.Normal = [xs, ys, zs]

# Add the continents, plates and hotspots
continents_vtk = LegacyVTKReader(guiName="continents.vtk",
                                 FileNames=[continent_file])
plates_vtk = LegacyVTKReader(guiName="plates.vtk", FileNames=[plates_file])
hotspot_vtk = LegacyVTKReader(guiName="hotspot.vtk", FileNames=[hotspot_file])

# only show the hotspots for the selected region, to avoid confusion!
SetActiveSource(hotspot_vtk)
hs1 = Clip(guiName="hs1", Scalars=['POINTS', value_vtk],
           ClipType="Plane", Value=0.0004882216453552246)
hs1.ClipType.Normal = [x1, y1, z1]
hs2 = Clip(guiName="hs2", InsideOut=1, Scalars=['POINTS', 'hs1'],
           ClipType="Plane", Value=1)
hs2.ClipType.Normal = [x2, y2, z2]
hs3 = Clip(guiName="hs3", InsideOut=1, Scalars=['POINTS', 'hs2'],
           ClipType="Plane", Value=1)
hs3.ClipType.Normal = [x3, y3, z3]
hs4 = Clip(guiName="hs4", Scalars=['POINTS', 'hs3'],
           ClipType="Plane", Value=0.0007011685520410538)
hs4.ClipType.Normal = [x4, y4, z4]

SetActiveSource(solx1_vtk)
RenderView1 = CreateRenderView()
RenderView1.CacheKey = 0.0
RenderView1.StereoType = 0
RenderView1.UseLight = 1
RenderView1.StereoRender = 0
RenderView1.CameraPosition = \
    [-0.11324999999987995, -0.23999999999978172, 42967.9512208794]
RenderView1.StereoCapableWindow = 0
RenderView1.OrientationAxesVisibility = 0
RenderView1.CameraClippingRange = [27132.33286101435, 63005.5539315754]
RenderView1.RemoteRenderThreshold = 3.0
RenderView1.Background = \
    [0.31999694819562063, 0.3400015259021897, 0.4299992370489052]
RenderView1.LightIntensity = 0.51
RenderView1.CameraFocalPoint = \
    [-0.11324999999987995, -0.23999999999978172, -0.07980000000043219]
RenderView1.CenterAxesVisibility = 0
RenderView1.CameraParallelScale = 11120.944758759495
RenderView1.CenterOfRotation = \
    [-0.11324999999987995, -0.23999999999978172, -0.07980000000043219]
RenderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]

# DataRepresentation1 = Show()
# DataRepresentation1.EdgeColor = [0.0, 0.0, 0.5000076295109483]
# DataRepresentation1.SelectionPointFieldDataArrayName = 'solx1_vtk'
# DataRepresentation1.ScalarOpacityFunction = solx1_PiecewiseFunction
# DataRepresentation1.ColorArrayName = ('POINT_DATA', 'solx1_vtk')
# DataRepresentation1.ScalarOpacityUnitDistance = 173.83040694044624
# DataRepresentation1.Visibility = 0
# DataRepresentation1.LookupTable = solx1_PVLookupTable
# DataRepresentation1.ScaleFactor = 1284.1448

SetActiveSource(Clip4)
DataRepresentation3 = Show()
DataRepresentation3.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation3.SelectionPointFieldDataArrayName = value_vtk
DataRepresentation3.ScalarOpacityFunction = solx1_PiecewiseFunction
DataRepresentation3.ColorArrayName = ('POINT_DATA', value_vtk)
DataRepresentation3.ScalarOpacityUnitDistance = 201.10756532082027
DataRepresentation3.LookupTable = solx1_PVLookupTable
DataRepresentation3.ScaleFactor = 1284.12626

SetActiveSource(Clip5)
DataRepresentation3 = Show()
DataRepresentation3.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation3.SelectionPointFieldDataArrayName = value_vtk
DataRepresentation3.ScalarOpacityFunction = solx1_PiecewiseFunction
DataRepresentation3.ColorArrayName = ('POINT_DATA', value_vtk)
DataRepresentation3.ScalarOpacityUnitDistance = 201.10756532082027
DataRepresentation3.LookupTable = solx1_PVLookupTable
DataRepresentation3.ScaleFactor = 1284.12626

SetActiveSource(Slice2)
DataRepresentation39 = Show()
DataRepresentation39.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation39.SelectionPointFieldDataArrayName = value_vtk
DataRepresentation39.ColorArrayName = ('POINT_DATA', value_vtk)
DataRepresentation39.LookupTable = solx1_PVLookupTable
DataRepresentation39.ScaleFactor = 1284.0803055356419

SetActiveSource(continents_vtk)
DataRepresentation10 = Show()
DataRepresentation10.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation10.Scale = [6371.0, 6371.0, 6371.0]
DataRepresentation10.SelectionPointFieldDataArrayName = ' Zcoord'
DataRepresentation10.SelectionCellFieldDataArrayName = 'Material Id'
DataRepresentation10.ScalarOpacityFunction = solx1_PiecewiseFunction
DataRepresentation10.ColorArrayName = ('POINT_DATA', ' Zcoord')
DataRepresentation10.ScalarOpacityUnitDistance = 0.20179812480658527
DataRepresentation10.LookupTable = solx1_PVLookupTable
DataRepresentation10.LineWidth = 3.0
DataRepresentation10.ScaleFactor = 0.20090600252151491

SetActiveSource(plates_vtk)
DataRepresentation15 = Show()
DataRepresentation15.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation15.Scale = [6471.0, 6471.0, 6471.0]
DataRepresentation15.DiffuseColor = [0.0, 0.6666666666666666, 0.0]
DataRepresentation15.ScalarOpacityUnitDistance = 0.1881526719997256
DataRepresentation15.LineWidth = 2.0
DataRepresentation15.ScaleFactor = 0.19997270107269288

SetActiveSource(hs4)
DataRepresentation15 = Show()
DataRepresentation15.EdgeColor = [0.0, 1.0, 0.0]
DataRepresentation15.Scale = [6471.0, 6471.0, 6471.0]
DataRepresentation15.DiffuseColor = [0.0, 0.6666666666666666, 0.0]
DataRepresentation15.ScalarOpacityUnitDistance = 0.1881526719997256
DataRepresentation15.LineWidth = 2.0
DataRepresentation15.PointSize = 10.0
DataRepresentation15.ScaleFactor = 0.19997270107269288

Render()
