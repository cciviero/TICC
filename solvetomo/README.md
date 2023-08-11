# solvetomo

*last step: inversion*

## run it

open run_solvetomo.py and edit the input part of the code, save
and execute make sure that you chose the right mpi that runs on your
computer...

mpirun.openmpi -np 16 mpisolvetomo < in.st

or

mpiexec -np 16 mpisolvetomo < in.st

## Encountered Problems

If your anaconda environment is on AND mpi does not work, try:

*soure deactivate*

This might happen on your local machine and/or on a server system...


## INPUT

## OUTPUT


* solx.< ident >: solution - velocity model resulting from iteration 
* sol.< ident >: full solution, raw chis)
* res.< ident >: data residuals d-A*m
* diagnostics.solve.< ident >
* outl.< ident > outliers
* xy.tradeoff.< ident >: GMT file for Lcurve
* out.solve.< ident >
* tex.solve.< ident >: latex tables


