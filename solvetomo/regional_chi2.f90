program regional_chi2
  !-----------------------
  ! regional_chi2
  !-----------------------
  include "../includes/solvet.h"
  character(NFN) :: ident, res_file
  character(256) :: path_dir
  integer :: lastrow, krow, outl, counter, kount
  integer :: ibandkrow, kstatkrow, kklustkrow, noutl
  dimension :: rhs(NMAX),sigma(NMAX),ievt(NMAX)
  dimension :: kklust(NMAX),kstat(NMAX)
  dimension :: res(NMAX)
  character(20) :: stationcodekrow,stationcode(NMAX), networkcode(NMAX)
  character(20) :: networkcodekrow
  character(10) :: phase_mapper(10)
  dimension :: iband(NMAX),knt(NMAX), grpw(NMAX)
  dimension :: igrp(NMAX),dtheor(NMAX)
  dimension :: statlat(NESMAX),statlon(NESMAX),evlat(NESMAX),evlon(NESMAX),evdep(NESMAX)
  real :: statlat_tmp, statlon_tmp, evlat_tmp, evlon_tmp, evdep_tmp
  real :: avdt_tmp, avda_tmp
  integer :: num_def_regions=10, num_def_phases=10
  real(8) :: chi2_phase(10), chi2_regional(10), chi2_total
  integer :: ph_type, reg_type, ph_count(10), reg_count(10)
  real :: lon1, lon2, lat1, lat2
  integer :: j, ldr

  ! ====================== PROGRAM STARTS ======================================
  do j=1,num_def_phases
    ph_count(j) = 0
    chi2_phase(j) = 0
  end do
  do j=1,num_def_regions
    reg_count(j) = 0
    chi2_regional(j) = 0
  end do

  ! INPUTS:
  ! path_dir: of the inversion
  ! ident, you know what that is...if not, refer to all the inversion codes
  print *, "read path_dir..."
  read(5, fmt='(a)') path_dir
  ldr=len(trim(path_dir))
  if(path_dir(ldr:ldr).ne.'/') then
    ldr=ldr+1
    path_dir(ldr:ldr)='/'
  endif
  print *, "read ident and res_file..."
  read(5, fmt='(a)') ident
  read(5, fmt='(a)') res_file
  print *, 'path:', path_dir(1:ldr)
  print *, 'ident:', adjustl(trim(ident))
  print *, 'res file:', adjustl(trim(res_file))

  ! ----- reading outl.<ident> file
  print *, 'Opening outl.<ident>'
  open(15, file=path_dir(1:ldr)//'outl.'//ident, status='old', FORM='FORMATTED')

  ! ----- reading the assemblematrix.stations
  print *, 'Opening assemblematrix.stations file...'
  open(22, iostat=ios, file=path_dir(1:ldr)//'assemblematrix.stations.'//ident, form='formatted')
  read(22, *)
  read(22, *)
  i = 0
  do
    read(22, 340, iostat=ios) statlist_tmp,netwlist_tmp,mstat_1,mstat_2, &
                avdt_tmp,avda_tmp,statlat_tmp,statlon_tmp,statelev_tmp
    340   format(a16,1x,a8,2i8,2f10.3,3f10.3)
    if (ios.lt.0) exit
    i = i+1
    statlat(i) = statlat_tmp
    statlon(i) = statlon_tmp
  end do
  num_stats = i
  print *, 'Number of stations: ', num_stats
  close(22)
  write(*,*) "-------------------"

  ! ----- reading the assemblematrix.events
  print *, 'Opening assemblematrix.events file...'
  open(33,iostat=ios, file=path_dir(1:ldr)//'assemblematrix.events.'//ident, form='formatted')
  read(33, *)
  read(33, *)
  i = 0
  do
    read(33, *, iostat=ios) evlist_tmp,klustlist_tmp,mevt_1,mevt_2, &
                        avdt_tmp,avda_tmp,evlat_tmp,evlon_tmp,evdep_tmp
    !342   format(i8,1x,i6,2i8,2f10.3,3f10.3)
    !print *, ios
    if (ios.ne.0) exit
    i = i+1
    evlat(i) = evlat_tmp
    evlon(i) = evlon_tmp
    evdep(i) = evdep_tmp
  end do
  num_events = i
  print *, 'number of events: ', num_events
  close(33)
  write(*,*) "-------------------"

  ! ----- reading aux.<ident> file
  print *, 'Opening aux.<ident>'
  open(2,file=path_dir(1:ldr)//'aux.'//ident, status='old', form='formatted')

  print *, 'Opening res file'
  open(33,file=path_dir(1:ldr)//adjustl(trim(res_file))//'.'//adjustl(trim(ident)), status='old', form='formatted')

  print *, 'Opening phase_type file'
  open(44,file='./phase_type.txt', status='old', form='formatted')

  print *, 'Opening phase_region.txt file'
  open(55,file='./phase_region.txt', form='formatted')

  print *, 'Start reading aux rows...'
  print *, 'x will be printed every 5K'
  print *, 'info will be printed every 100K'
  counter = 0
  chi2_total = 0
  noutl = 0
  do
    counter = counter + 1
    read(2,*,iostat=ios) krow,kount,rhskrow,sigmakrow,ievtkrow, &
              kklustkrow,stationcodekrow,networkcodekrow,kstatkrow, &
              ibandkrow,krtyp,klust1,igroup,grpwkrow,dthkrow
    if (krow.eq.0) exit
    if (mod(counter, 5000) == 0) write(6,'(a1,$)') 'x'
    if (mod(counter, 100000) == 0) write(*,*) counter, krow

    lastrow=krow
    rhs(krow)=rhskrow
    sigma(krow)=sigmakrow
    ievt(krow)=ievtkrow
    kklust(krow)=kklustkrow
    kstat(krow)=kstatkrow
    stationcode(krow)=stationcodekrow
    networkcode(krow)=networkcodekrow
    iband(krow)=ibandkrow
    knt(krow)=kount
    grpw(krow)=grpwkrow
    igrp(krow)=igroup
    dtheor(krow)=dthkrow

    ! This row is an outlier or not?
    read(15,*,iostat=ios2) outl
    read(33, *, iostat=ios) res(krow)
    read(44, *, iostat=ios) ph_type

    colat1 = 90. - evlat(kkluskrow)
    lat1 = evlat(kkluskrow)
    lon1 = evlon(kkluskrow)
    colat2 = 90. - statlat(kstatkrow)
    lon2 = statlon(kstatkrow)
    lat2 = statlat(kstatkrow)

    ! Calculate distance and azimuth
    call stadis(colat1,lon1,colat2,lon2,del,az,t0,p0,t2,p1)
    reg_type = 0

    if (outl < 0.5) then
      chi2_total = chi2_total + res(krow)**2
      chi2_phase(ph_type) = chi2_phase(ph_type) + res(krow)**2
      ph_count(ph_type) = ph_count(ph_type) + 1

      ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ADD REGIONS
      if ((lon2.ge.-135).and.(lon2.le.-50).and.(lat2.ge.20).and.(lat2.le.50).and.(ph_type==1.or.ph_type==5)) then
        chi2_regional(1) = chi2_regional(1) + res(krow) ** 2
        reg_count(1) = reg_count(1) + 1
      endif
      if ((lon2.ge.-135).and.(lon2.le.-50).and.(lat2.ge.-4).and.(lat2.le.20).and.(ph_type==1.or.ph_type==5)) then
        chi2_regional(2) = chi2_regional(2) + res(krow) ** 2
        reg_count(2) = reg_count(2) + 1
      endif
      if ((lon2.ge.-82).and.(lon2.le.-32).and.(lat2.ge.-57).and.(lat2.le.12).and.(ph_type==1.or.ph_type==5)) then
        chi2_regional(3) = chi2_regional(3) + res(krow) ** 2
        reg_count(3) = reg_count(3) + 1
      endif
      if ((lon2.ge.-13).and.(lon2.le.50).and.(lat2.ge.26).and.(lat2.le.58).and.(ph_type==1.or.ph_type==5)) then
        chi2_regional(4) = chi2_regional(4) + res(krow) ** 2
        reg_count(4) = reg_count(4) + 1
      endif
      if ((lon2.ge.90).and.(lon2.le.154).and.(lat2.ge.-14).and.(lat2.le.23).and.(ph_type==1.or.ph_type==5)) then
        chi2_regional(5) = chi2_regional(5) + res(krow) ** 2
        reg_count(5) = reg_count(5) + 1
      endif
      if ((lon2.ge.41).and.(lon2.le.72).and.(lat2.ge.-38).and.(lat2.le.-8).and.(ph_type==1.or.ph_type==5)) then
        chi2_regional(6) = chi2_regional(6) + res(krow) ** 2
        reg_count(6) = reg_count(6) + 1
      endif
      if ((lon2.ge.-174).and.(lon2.le.-138).and.(lat2.ge.-6).and.(lat2.le.34).and.(ph_type==1.or.ph_type==5)) then
        chi2_regional(7) = chi2_regional(7) + res(krow) ** 2
        reg_count(7) = reg_count(7) + 1
      endif
      if ((lon2.ge.-25).and.(lon2.le.-10).and.(lat2.ge.60).and.(lat2.le.68).and.(ph_type==1.or.ph_type==5)) then
        chi2_regional(8) = chi2_regional(8) + res(krow) ** 2
        reg_count(8) = reg_count(8) + 1
      endif
      if ((lon2.ge.-120).and.(lon2.le.-32.5).and.(lat2.ge.-60).and.(lat2.le.20.0).and.(ph_type==1.or.ph_type==5)) then
        chi2_regional(9) = chi2_regional(9) + res(krow) ** 2
        reg_count(9) = reg_count(9) + 1
      endif
      if ((lon2.ge.-7.1).and.(lon2.le.-1.5).and.(lat2.ge.56.1).and.(lat2.le.59.2).and.(ph_type==1.or.ph_type==5)) then
        chi2_regional(10) = chi2_regional(10) + res(krow) ** 2
        reg_count(10) = reg_count(10) + 1
      endif
      ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> END ADD REGIONS

    else
      noutl = noutl + 1
      reg_type = -12345
    endif
    write(55, *) reg_type, ph_type
  end do
  print *, 'Last row:', lastrow
  nrow=lastrow
  print *, res_file
  print *, 'chi2 (total) = ', chi2_total/(nrow - noutl)
  print *, '------------Phase'
  phase_mapper(1) = 'PM'
  phase_mapper(2) = 'PPM'
  phase_mapper(3) = 'PdiM'
  phase_mapper(4) = 'iscP'
  phase_mapper(5) = 'P_RR'
  phase_mapper(6) = 'iscPP'
  phase_mapper(7) = 'iscPKPbcn1'
  do j=1,10
    print *, j, phase_mapper(j), chi2_phase(j), ph_count(j), chi2_phase(j)/ph_count(j)
  end do
  print *, '------------Regions'
  do j=1,10
    print *, j, chi2_regional(j), reg_count(j), chi2_regional(j)/reg_count(j)
  end do
end program regional_chi2


subroutine stadis(colat,colon,scolat,scolon,del,az,t0,p0,t1,p1)

! Computes the epicentral distance and azimuth from source to receiver.
! Latitudes are converted to geocentric latitudes prior to performing
! the computations (it is assumed that input latitudes are geographic).
!  input:
!    colat = source colatitude, (degrees)
!    colon =   "    colongitude,  "
!    scolat    = station colatidue,
!    scolon    =    "    colongitude,
!  output:
!    del   = epicentral distance (degrees)
!    az    = azimuth from source to receiver, measure from North (degrees)
      data rpd/1.745329252e-2/,dpr/57.29577951/

!  first do eq coords.
      if(colat.eq.0.) colat=1.0e-5
      t0=colat*rpd
      p0=colon*rpd
      c0= cos(t0)
      s0= sin(t0)
!     print *,t0,p0,c0,s0
!  now do station coords.
      t1=scolat*rpd
      c1= cos(t1)
      s1= sin(t1)
      p1=scolon*rpd
!     print *,t1,c1,s1,p1
!  now calculate distance
      dp=p1-p0
      co=c0*c1+s0*s1*cos(dp)
      si=sqrt(1.-co*co)
!     print *,dp,co,si
      del=atan2(si,co)*dpr
!  now calculate azimuth
      caz=(c1-c0*co)/max(1.0e-30,si*s0)
      dp2=-dp
      saz=-s1*sin(dp2)/si
      az= atan2(saz,caz)*dpr
      if(az.lt.0.0) az=360.0 + az
!     print *,'stadis returns',del,az
      return
      end

! convert a point from cartesian to spherical coordinates

subroutine xyz2rtp(x,y,z,r,t,p)

! input: cartesian coordinates x,y,z
! output: spherical coordinates r,t (theta),p (phi)

data halfpi/1.570796327/

r=sqrt(x*x+y*y+z*z)       !radius
t=acos(z/r)               ! colatitude
if (x.eq.0.0) then        ! make sure return is +pi/2
   p=halfpi
   return
endif
p=atan2(y,x)              ! longitude

return
end

subroutine roundme(x, ndigits, x_round)
real(8) :: x, x_round
integer :: ndigits
! input: x, a real number
! input: ndigits, number of digits
! output: x_round, rounded number
B=FLOAT (INT(A * 1000.0 + 0.5)) / 1000.0
x_round = float(int(x*10.*ndigits + 0.5)) / (10.*ndigits)
return
end
