# raydata raymatrix

## how to run
1. run add_eventx.py
1. get your input/input...ini file right
1. run raydata_raymatrix.py:
    * python raydata_raymatrix.py input/input_raydata_raymatrix.ini

## To dos

* fix issue with raymatrix not compiling/executing right!
* look at vtk_generator, if matrix file is empty

## Important notes on the codes

### raymatrix.f

Where to find:
[raymatrix](src_raydata_raymatrix/raymatrix_src/raymatrix.f)

##### Kernel information for P and S
Kernel information is stored separately for P and S in the matrix

L928: j0=(jtype-1)*np

for P jtype=1
for S jtype=2

Example: if we have a matrix of np=100 (np=number of points), then:
j0=0 for P and
j0=100 for S

L929: do j=1,np is the same for P and S but in

L959: ja(kount)=j0+j we get ja=j for P but ja=np+j for S

Hence our S matrix is bigger then our defined grid, which can cause problems
for generating vtk files later from the ascii files if not taken care of.

(23.03.2017 MT: taken care of for vtk generating functions...)



