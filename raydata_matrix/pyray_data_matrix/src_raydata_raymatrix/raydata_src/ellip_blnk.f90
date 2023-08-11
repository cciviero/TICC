SUBROUTINE ellip_blnk(phase, del, sdep, scolat, az, tcor)
! =======================================
! phase:
! del: epicentral distance to station (in degrees) 
! sdep: depth of event (KM)
! scolat: co-latitude of the source (in degree)
! az: azimuth from source to station (in degrees)
! tcor: time correction for path (ellipticity)
! =======================================


IMPLICIT NONE

character(8), INTENT(IN) :: phase
real, INTENT(IN) :: del, sdep, scolat, az
real, INTENT(OUT) :: tcor
real :: degrad, edist, azim, edepth, elat, ecolat, zs
real, dimension(2) :: usrc
logical :: abrt

!degree ---> radian
degrad = 45.0/atan(1.0)

! open ellipticity correction file
open(21,file='elcordir.tbl',access='direct',form='formatted',recl=80)

edist = del
azim = az/degrad
zs = sdep

call depset(zs,usrc)
edepth = zs

ecolat = scolat/degrad
call ellref(ecolat)

!phase  : a  string specifying the PHASE,   -e.g P, ScP etc.  
!edist  :  epicentral distance to station (in degrees)     
!edepth :  depth of event         
!ecolat :  epicentral co-latitude of source (in radians) 
!azim   :  azimuth from source to station (in radians)
!tcor   :  time correction for path to allow for ellipticity
call ellcor(phase, edist, edepth, ecolat, azim, tcor, abrt)
!WRITE(*,*) "corrected time: ", tcor

END SUBROUTINE ellip_blnk
