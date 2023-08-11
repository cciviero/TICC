#  how to run the resolution test

create new folder for the resolution test of this model version

./make

move output to the folder you want to do the resolution test
and copy the ruh_resolution.sh too

files needed:

vertices
facets

- mat.indent
ata.indent
- atamax.indent
columndensity.1.indent
- aux.indent

options 7
radius 300
anomaly -3
noise no


## Input file for resolution test needs following information

Example for ./resoltuion input: bold part is your input

 Give vertex file name (eg vertices.all):
**vertices.RRx2**

 Read vertices.RRx2
 total nr of vertices np=       56709
 Give ident of the matrix file:
**P_RRx2**

 Assembled matrix file is: mat.P_RRx2
 Auxiliary (data) file is: aux.P_RRx2

 Matrix nrow,ncol=       76791       62964
 Last row in aux file:       76791
 nr of stations:       17541
 nr of clusters:        4380
 Finished reading        76791  data from        17541  stations
 and        4380  events/clusters

 Cannot read kssw from aux file - ignored

 Test model options for resolutiontest:

 -1=no structural features, scramble data
 0=no structural features, noise only
 1=Gaussian anomalies at user specified locations
 2=Gaussian anomalies at regular intervals
 3=Point anomalies at user specified locations
 4=Point anomalies at regular intervals
 5=Global spherical harmonics
 6=Checkerboard
 7=Vertical plumes at specified locations

 Other:
 11=generate dT values from user-provided model
 Give test option:
**7**

 Specify locations, depth extent, and radius of plumes
 (radius = 1/e width of Gaussian shaped cross-section)
 Give 1/e radius of plume conduit in km:
**300**

 Give amplitude of maximum anomaly in %:
**-3**

 After last plume, type -1 -1 -1 -1
 Give lon (deg), lat (deg), min_depth, max_depth (km):
**55 -21 0 600
35.5 -27 600 1200
29.5 -25 1200 2900
24.5 -27.5 2400 2900
-1 -1 -1 -1**

 Give ident for synthetic model output file syn.<ident>:
**plume_300_obliqII**

Writing model to file syn.plume_300_obliqII  
 Model min/max values are:  -3.94608490E-02   0.00000000  
 Nonzero matrix elements:        2001  out of        56709  
 Statistics:  
 Largest element of matrix =    11.0933704  
 Largest predicted datum =    8.07435417  
 Before adding noise:  
 Data min/max values are:  -3.32151838E-02   8.07435417  
 Nonzero data:       42371  out of        76791  

 Noise options:  
 0: no noise  
 1: Gaussian noise  
 Give noise option:  
**0**

 Give ident for new aux (data) output file:
**plume_300_obliqII**

 New auxiliary (data) file is: aux.plume_300_obliqII


 You can also create an **in.resolution** file and give input! Eg:

**vertices.RRx2  
  P_RRx2  
  7  
  300  
  -3  
  55 -21 0 600  
  35.5 -27 600 1200  
  29.5 -25 1200 2900  
  24.5 -27.5 2400 2900  
  -1 -1 -1 -1  
  plume_300_obli qII  
  0  
  plume_300_obliqII**

