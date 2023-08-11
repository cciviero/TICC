      program resolutiontest

      ! matrix rows  for at most MMAX nonzero elements per row, 
      ! NDIM unknowns (including corrections) and NMAX data
      ! NOTE: Longitude input must range from -180 to 180 degrees 
      !
      ! Compile:
      ! g77  -o resolutiontest resolutiontest.f
      ! Note: static complation needed on alhazen -static
      ! g77 -g -fbounds-check -fno-automatic -W -o resolutiontest resolutiontest.f
      ! 2010/11/08: Compilation at LMU: 
      ! gfortran -std=legacy -o resolutiontest resolutiontest.f
      !
      ! Options for resolution test patterns are:
      ! 1: isolated Gaussians
      ! 2: regularly spaced Gaussians
      ! 3: point anomalies at user-specified locations
      ! 4: point anomalies at regular intervals
      ! 5: global spherical harmonics
      ! 6: checkerboard (i.e., hard thresholding) 
      ! 7: vertical plumes at specified locations (2-D Gaussian columns)
      ! 11: resolution test from user-provided model
      ! 12:generate dT values from user-provided model
      !
      ! Input:
      ! (1) vertices.<mesh_ident> : tetrahedral mesh file
      ! (2) aux.<drun_ident>  (real) data file, rhs, sigma normalized 
      ! (3) mat.<drun_ident>  kernel matrix file, binary
      !   where drun_ident is the ident of a previous (real data) 
      !   run, to which the resolution test is to correpond.
      !   aux.<drun_ident> is read for dimensional info only.
      ! Input files must be present locally (in runspace).
      !
      ! Output:
      ! (4) resolutiontest.<rt_ident> : log file, ascii
      ! (1) syn.<rt_ident> : synthetic resolution test
      !       model file (not normalized by prior 
      !       parameter uncertainties)
      ! (2) aux.<rt_ident> : the synthetic "data"
      !       corresponding to the syn file (with noise,
      !       if desired, sigma-normalized)
      !
      ! Example of command-line input (prompted)
      ! (cat in.resolutiontest.rt400_3.USA101_RT)
      !
      ! vertices.USA10         # mesh file
      ! USA101_RT              # drun_ident
      ! 2                      # Gaussians
      ! 200                    # half-width of Gaussian peaks 40
      ! 3                      # amplitude in %             1000
      ! -160 -40    8          # lonmin, lonmax, dlon
      ! 18   58      8         # latmin, latmax, dlat
      ! 3300 6300  600         # rmin,   rmax,   dr
      ! rt400_3.USA101_RT      # rt_ident (for syn)
      ! 1                      # noise y/n
      ! 0.6                    # noise std (i.e., 1/SNR)
      ! rt400_3.USA101_RT      # rt_ident2 (for aux. file)
      !                        # usually identical to rt_ident
 

      ! history:
      ! june 2010 added option to compute theoretical travel times
      ! 201/08/26 K.S.: Line 
      ! parameter(MMAX=200000,NMAX=2000000,NDIM=350000
      ! was replaced by 
      ! include "../includes/solvet.h" 

      include "/includes/solvet.h"
!      parameter(MMAX=200000,NMAX=2000000,NDIM=350000)
!      parameter(MMAX=200000,NMAX=2000000,NDIM=200000)
      dimension a(MMAX),ja(MMAX),knt(NMAX)

      ! data
      dimension rhs(NMAX),sigma(NMAX),ievt(NMAX),dtheor(NMAX)
      dimension iband(NMAX),kklust(NMAX),sol(NDIM),kstat(NMAX)
      dimension counter(NDIM)
      ! plates needs to be loaded as cartesian coordinates
      dimension plates(6048)
      ! gerror = Given ERROR
      dimension gerror(NMAX)
      integer evlist(500000),klustlist(500000)

      ! header files for Delauney grid
      include '/includes/nn.param'
      common/nndim/np,nt
      include '/includes/setdel2.h' ! common for points,centres,vertices etc

      ! other
      dimension jssw(10),jlast(11),nssw(3),kssw(4),kountpar(3),dpar(10)
      dimension grpw(NMAX),igrp(NMAX)

      ! various character variables 
      character*256 fname, processing
      character*256 ident,ident2,ident3, remval
      character*16 stationcodekrow,stationcode(NMAX),statlist(500000)
      character*8 phase,networkcode,netwlist(500000)
      character*256 directory


      ! Normally jdebug=0. Set -1 if input and output files are ascii i.o.
      ! binary. Set 1 for some debugging output. Set 2 if avpu and atupv
      ! were changed so it just does the dot-product tests, then stops.

      jdebug=0    
      pi=4.0*atan(1.0)
      d2r=pi/180.0


c read vertices      
      print *,'Give vertex file name (eg vertices.all):'
      read(5,fmt='(a)') fname
      open(1,file=fname)
      read(1,*) nd
      if(nd.ne.3) stop 'nd not equal to 3'
      read(1,*) np
      if(np.gt.np_max) stop 'np>np_max'
      do i=1,np  
        read(1,*) (points(j,i),j=1,3)
      enddo  
      write(6,*) 'Read ',fname,' total nr of vertices np=',np
      close(1)

c get matrix
! open the assembled matrix file (binary format only)
      print *,'Give ident of the matrix file:'
      read(5,fmt='(a)') ident
      print *,'Assembled matrix file is: mat.'//ident
      print *,'Auxiliary (data) file is: aux.'//ident
      open(3,file='mat.'//ident,form='unformatted')
      read(3) nrow,ncol
      read(3) jlast
      read(3) jssw
      read(3) dpar
      npar=jlast(3)
      open(4,file='resolutiontest.out.'//ident)
      write(4,*) 'Resolution test for grid: ',fname
      write(4,*) 'and matrix: mat.'//ident
      print *,'Matrix nrow,ncol=',nrow,ncol
      write(4,*) ' nrow,ncol=',nrow,ncol
      nrow0=nrow
      if(ncol.gt.NDIM) stop 'ncol>NDIM'
      open(2,file='aux.'//ident)
      if(jdebug.gt.0) write(13,*) 'Now reading aux file'

      ! read aux file with data information

      ! The data in the aux file (rhs) are scaled to variance=1
      ! but sigma stores the original standard deviation
      ! (they are also scaled by grpwkrow=chiw if data are grouped)
100   read(2,*) krow,kount,rhskrow,sigmakrow,ievtkrow,
     &      kklustkrow,stationcodekrow,networkcode,kstatkrow,
     &      ibandkrow,krtyp,klust1,igroup,grpwkrow,dthkrow
      if(jdebug.gt.0) write(13,*) 'krow=',krow,' kount=',kount,ievtkrow,
     &  stationcodekrow
      if(krow.eq.0) goto 140
      lastrow=krow
      rhs(krow)=rhskrow
      sigma(krow)=sigmakrow
      ievt(krow)=ievtkrow
      kklust(krow)=kklustkrow
      kstat(krow)=kstatkrow
      stationcode(krow)=stationcodekrow
      iband(krow)=ibandkrow
      knt(krow)=kount
      grpw(krow)=grpwkrow
      igrp(krow)=igroup
      dtheor(krow)=dthkrow

      goto 100
140   nrow=lastrow              ! number of data
      print *,'Last row in aux file:',lastrow
      ! read rest of aux file
      read(2,*) nstat
      print *,'nr of stations:',nstat
      read(2,*) (nn,statlist(i),netwlist(i),i=1,nstat)
150   format(i5,/,(i5,1x,a16,1x,a8))
      read(2,*) nklust
      print *,'nr of clusters:',nklust
      read(2,*) (nn,evlist(i),klustlist(i),i=1,nklust)
160   format(i5,/,(i5,2i10))
      print *,'Finished reading ',nrow,' data from ',nstat,' stations'
      print *,'and',nklust,' events/clusters'
      write(4,*) ' reading ',nrow,' data from ',nstat,' stations'
      write(4,*) ' and',nklust,' events/clusters'
      do i=1,ngroup
        read(2,*,iostat=ios) ii,n1,c1
        if(ios.ne.0) then
          print *,'Cannot read group info from aux file - ignored'
          goto 161
        endif  
      enddo
      ! older aux files do not have kwws at the end, so we use iostat
161   read(2,*,iostat=ios) kssw
      if(ios.ne.0) then
        print *,'Cannot read kssw from aux file - ignored'
        do i=1,4
          kssw(i)=-99
        enddo  
      endif  

      if(nrow.ne.nrow0) stop 'nr of data incompatible with matrix size'

      ! zero the test model (incluing corrections)
      do i=1,ncol
        sol(i)=0.
      enddo  

      print *,'Test model options for resolutiontest:'
      print *,'-1=no structural features, scramble data'
      print *,'0=no structural features, noise only'
      print *,'1=Gaussian anomalies at user specified locations'
      print *,'2=Gaussian anomalies at regular intervals'
      print *,'3=Point anomalies at user specified locations'
      print *,'32=Point anomalies at user specified'
      print *, 'locations and amplitude'
      print *,'4=Point anomalies at regular intervals'
      print *,'5=Global spherical harmonics'
      print *,'6=Checkerboard'
      print *,'7=Vertical plumes at specified locations'
      print *,'8=Vertical plumes - specify xlon,xlat,&
     &depthmin,depthmax and width'
      print *,'9=Vertical plumes - specify xlon,xlat,&
     &depthmin,depthmax, width and max amplitude'
      print *,'10=the snake: mid-ocean ridge horizontal columns'
      print *,'11=resolution test from user-provided model'
      print *,'12=resolution test from user-provided model WITH'
      print *,'zero valued Gaussian anomalies at'
      print *,'user specific loactions'
      print *,'13=resolution test from user-provided model WITH'
      print *,'a band of specified by rmin/rmax set to zero'
      print *,'14=resolution test from user-provided model WITHOUT'
      print *,'either high or low velocities'
      print *
      print *,'Other:'
      print *,'15=generate dT values from user-provided model'
      print *
      print *


      print *,'Give test option:'
      read *,koption
      print *, 'Test option used: ',koption
      print *

      ntst=0

      if(koption.eq.-1) then
        noption=2
        goto 600

      ! Noise only
      else if(koption.eq.0) then
        do i=1,nrow
          rhs(i)=0.
        enddo  
        noption=1
        goto 600

      !-----------------------------------
      ! Which input pattern to generate:
      !-----------------------------------

      !-----------------------------------
      ! 1: Isolated Gaussians
      !-----------------------------------
      else if(koption.eq.1) then
        print *,'Give 1/e half-width of Gaussian peak in km:'
        read *,width
        w2=width**2
        dmax=3.0*width
        dmax2=dmax**2
        print *,'Give amplitude of maximum in %:'
        read *,amp
        write(4,*) 'Half-width: ',width,', amp:',amp,'%'
        amp=0.01*amp    ! convert from %
        print *,'After last peak type 0,0,0-----------------------'
        print *,'Give radius (km), longitude (deg), latitude (deg):'
        write(4,'("Radius    lon     lat")')
200     read *,r,xlon,xlat
        if(r.le.0.) goto 400
        write(4,'(f6.0,2f8.2)') r,xlon,xlat
        rlon=d2r*xlon
        rcolat=d2r*(90.-xlat)
        x0=r*cos(rlon)*sin(rcolat)
        y0=r*sin(rlon)*sin(rcolat)
        z0=r*cos(rcolat)
        do i=1,np       ! loop over all nodes
          xp=points(1,i) 
          yp=points(2,i) 
          zp=points(3,i) 
          d2=(xp-x0)**2+(yp-y0)**2+(zp-z0)**2
          if(d2.lt.dmax2) then
            sol(i)=sol(i)+amp*exp(-d2/w2)
          endif
        enddo
        goto 200

      !-----------------------------------
      ! 2: Regularly spaced Gaussians
      ! This yields a faithful rendering only if 
      ! the pulsewidth is significantly 
      ! smaller than the pulse spacing; else cancellation of several pulses
      !-----------------------------------
      else if(koption.eq.2) then
        print *,'Give 1/e width of Gaussian peak in km:'
        read *,width
        w2=width**2
        dmax=3.0*width ! farther away than that, assume anomaly is negligible
        dmax2=dmax**2
        print *,'Give amplitude of maximum in %:'
        read *,amp
        write(4,*) 'Width: ',width,', amp:',amp,'%'
        amp=0.01*amp    ! convert from %
        print *,'Give min, max longitude and spacing (deg):'
        read *,xlon1,xlon2,dlon
        print *,'Give min,max latitude and spacing :'
        read *,xlat1,xlat2,dlat
        print *,'Give min,max radius and spacing (km):'
        read *,r1,r2,dr
        write(4,*) 'Longitudes:',xlon1,xlon2,dlon
        write(4,*) 'Latitudes:',xlat1,xlat2,dlat
        write(4,*) 'Radii:',r1,r2,dr

        rlon1=d2r*xlon1
        rlon2=d2r*xlon2
        drlon=dlon*d2r
        drlat=dlat*d2r
        rcolat1=d2r*(90.-xlat2)
        rcolat2=d2r*(90.-xlat1)

        rsign=-1
        r=r2+dr
        do while(r.ge.r1)
          r=r-dr
          rsign=-rsign
          ysign= rsign
          rcolat=rcolat1-drlat
          do while(rcolat.le.rcolat2)
            rcolat=rcolat+drlat
            ysign=-ysign
            xsign=ysign
            rlon=rlon1-drlon
            do while(rlon.le.rlon2)
              rlon=rlon+drlon
              xsign=-xsign
              x0=r*cos(rlon)*sin(rcolat)
              y0=r*sin(rlon)*sin(rcolat)
              z0=r*cos(rcolat)
              do i=1,np       ! loop over all nodes
                xp=points(1,i) 
                yp=points(2,i) 
                zp=points(3,i) 
                d2=(xp-x0)**2+(yp-y0)**2+(zp-z0)**2
                if(d2.lt.dmax2) then
                  sol(i)=sol(i)+xsign*amp*exp(-d2/w2)
                endif
              enddo
            enddo
          enddo
        enddo  


      !-----------------------------------
      ! 3: Isolated spheres
      !-----------------------------------
      else if(koption.eq.3) then
        print *,'Give amplitude of spike in %:'
        read *,amp
        write(4,*) ' Spike Amp:',amp,'%'
        amp=0.01*amp    ! convert from %
        print *,'After last peak type 0,0,0-----------------------'
        print *,'Give radius (km), longitude (deg), latitude (deg):'
        write(4,'("Radius    lon     lat")')
210     read *,r,xlon,xlat
        if(r.le.0.) goto 400
        write(4,'(f6.0,2f8.2)') r,xlon,xlat
        ! find closest node in grid
        rlon=d2r*xlon
        rcolat=d2r*(90.-xlat)
        x0=r*cos(rlon)*sin(rcolat)
        y0=r*sin(rlon)*sin(rcolat)
        z0=r*cos(rcolat)
        dmax2=1.0e30
        iclosest=1
        do i=1,np       ! loop over all nodes
          xp=points(1,i) 
          yp=points(2,i) 
          zp=points(3,i) 
          d2=(xp-x0)**2+(yp-y0)**2+(zp-z0)**2
          if(d2.lt.dmax2) then
            iclosest=i
            dmax2=d2
          endif
        enddo
        sol(iclosest)=amp
        goto 210

      !-----------------------------------
      ! 32: Isolated spheres
      !-----------------------------------
      else if(koption.eq.32) then

        print *,'After last peak type 0,0,0,0------'
        print *,'Give radius (km), longitude
     &(deg), latitude (deg), amp:'
        write(4,'("Radius    lon     lat")')
211     read *,r,xlon,xlat, amp
        if(r.le.0.) goto 400
        amp=0.01*amp    ! convert from %
        write(4,'(f6.0,2f8.2)') r,xlon,xlat
        ! find closest node in grid
        rlon=d2r*xlon
        rcolat=d2r*(90.-xlat)
        x0=r*cos(rlon)*sin(rcolat)
        y0=r*sin(rlon)*sin(rcolat)
        z0=r*cos(rcolat)
        dmax2=1.0e30
        iclosest=1
        do i=1,np       ! loop over all nodes
          xp=points(1,i)
          yp=points(2,i)
          zp=points(3,i)
          d2=(xp-x0)**2+(yp-y0)**2+(zp-z0)**2
          if(d2.lt.dmax2) then
            iclosest=i
            dmax2=d2
          endif
        enddo
        sol(iclosest)=amp
        goto 211

      !-----------------------------------
      ! 33: Isolated spheres
      !-----------------------------------
      ! ATTENTION this functions makes sense
      ! as a hitcount visualisation, less
      ! as a resolution test...
      else if(koption.eq.33) then
        print *,'>>hit<<:hitcount or >>average<<:average dt'
        read *, processing
        print *, 'You choose:', processing

        print *,'Give path to x, y, z, ts file:'
        read(5,fmt='(a)') fname
        print *, 'You are using:', fname
        ! open(8,file=fname)
        open (unit=888, file=fname, STATUS = 'OLD')
212     read(888,*) x0,y0,z0,amp
        ! print *, x0, y0, z0, amp
        if((x0.eq.0).and.(y0.eq.0).and.(z0.eq.0)) then
            if (processing.eq.'average') then
                do i=1,np
                    ! print *, counter(i)
                    if (counter(i).gt.0) then
                        !print *, sol(i), counter(i)
                        sol(i) = sol(i)/counter(i)*0.01
                        !print *, sol(i)
                    end if
                end do
            end if
            goto 400
        end if
        ! find closest node in grid
        dmax2=1.0e30
        iclosest=1
        ! for finding the 4 closest points
        dmax3=1.0e30
        dmax4=1.0e30
        dmax5=1.0e30

        iclosest3=1
        iclosest4=1
        iclosest5=1

        do i=1,np       ! loop over all nodes
          xp=points(1,i)
          yp=points(2,i)
          zp=points(3,i)
          d2=(xp-x0)**2+(yp-y0)**2+(zp-z0)**2
          if(d2.lt.dmax5) then
            iclosest5=i
            dmax5=d2
            endif
          if(d2.lt.dmax4) then
            iclosest4=i
            dmax5=dmax4
            dmax4=d2
            endif
          if(d2.lt.dmax3) then
            iclosest3=i
            dmax4=dmax3
            dmax3=d2
            endif
          if(d2.lt.dmax2) then
            iclosest=i
            dmax3=dmax2
            dmax2=d2
            endif
        enddo
        ! now the four closest nodes are set to one
        ! and if the same nodes are affected it will
        ! count up - hit count map

        if(processing.eq.'hit') then
            sol(iclosest) = sol(iclosest) + 1
            sol(iclosest3) = sol(iclosest3) + 1
            sol(iclosest4) = sol(iclosest4) + 1
            sol(iclosest5) = sol(iclosest5) + 1

        else if(processing.eq.'average') then
            sol(iclosest) = sol(iclosest) + amp/4
            sol(iclosest3) = sol(iclosest3) + amp/4
            sol(iclosest4) = sol(iclosest4) + amp/4
            sol(iclosest5) = sol(iclosest5) + amp/4

            counter(iclosest) = counter(iclosest) + 1
            counter(iclosest3) = counter(iclosest3) + 1
            counter(iclosest4) = counter(iclosest4) + 1
            counter(iclosest5) = counter(iclosest5) + 1
            ! print *, sol(iclosest), amp, counter(iclosest)
        else
            stop 'Not implented. EXIT.'
        end if

        goto 212


      !-----------------------------------
      ! 4: Regularly spaced spheres
      !-----------------------------------
      else if(koption.eq.4) then
        print *,'Give amplitude of maximum in %:'
        read *,amp
        write(4,*) 'Spike amp:',amp,'%'
        amp=0.01*amp    ! convert from %
        print *,'Give min, max longitude and spacing (deg):'
        read *,xlon1,xlon2,dlon
        print *,'Give min,max latitude and spacing :'
        read *,xlat1,xlat2,dlat
        print *,'Give min,max radius and spacing (km):'
        read *,r1,r2,dr
        write(4,*) 'Longitudes:',xlon1,xlon2,dlon
        write(4,*) 'Latitudes:',xlat1,xlat2,dlat
        write(4,*) 'Radii:',r1,r2,dr

        rlon1=d2r*xlon1
        rlon2=d2r*xlon2
        drlon=dlon*d2r
        drlat=dlat*d2r
        rcolat1=d2r*(90.-xlat2)
        rcolat2=d2r*(90.-xlat1)

        rsign=-1
        r=r1-dr
        do while(r.le.r2)
          r=r+dr
          rsign=-rsign
          ysign=rsign
          rcolat=rcolat1-drlat
          do while(rcolat.le.rcolat2)
            rcolat=rcolat+drlat
            ysign=-ysign
            xsign=ysign
            rlon=rlon1
            do while(rlon.le.rlon2)
              rlon=rlon+drlon
              xsign=-xsign
              x0=r*cos(rlon)*sin(rcolat)
              y0=r*sin(rlon)*sin(rcolat)
              z0=r*cos(rcolat)
              ! Find closest node
              dmax2=1.0e30
              iclosest=1
              do i=1,np       ! loop over all nodes
                xp=points(1,i) 
                yp=points(2,i) 
                zp=points(3,i) 
                d2=(xp-x0)**2+(yp-y0)**2+(zp-z0)**2
                if(d2.lt.dmax2) then
                  iclosest=i
                  dmax2=d2
                endif
              enddo
              sol(iclosest)=xsign*amp
              print *,iclosest,rlon*57.3,90-rcolat*57.3,xsign,amp
            enddo
          enddo
        enddo  

      !-----------------------------------
      ! 5: Spherical harmonic
      !-----------------------------------
      else if(koption.eq.5) then
        print *,'Give orders l and m:'
        read *,l,m
        print *,'Give minimum and maximum radius (km):'
        read *,r1,r2
        print *,'Give wavelength for radial dependence:'
        read *,rwl
        print *,'Give amplitude (%):'
        read *,amp
        amp=0.01*amp
        write(4,*) 'l,m=',l,m
        write(4,*) 'Radii:',r1,r2,' wavelength:',rwl
        do i=1,np       ! loop over all nodes
          xp=points(1,i) 
          yp=points(2,i) 
          zp=points(3,i) 
          call xyz2rtp(xp,yp,zp,r,theta,phi)
          if(r.gt.r2.or.r.lt.r1) goto 300
          fr=sin(2*pi*(r-r1)/rwl)
          plm=plgndr(l,m,cos(theta))
          sol(i)=sol(i)+amp*fr*sin(m*phi)*plm
300       continue
        enddo

      !-----------------------------------
      ! 6: Checkerboard
      ! i.e. hard thresholding instead of smoothly tapered Gaussians
      !-----------------------------------
      else if(koption.eq.6) then
        print *,'Give 1/e width of Gaussian peak in km:'
        read *,width
        w2=width**2
        dmax=3.0*width
        dmax2=dmax**2
        print *,'Give amplitude of maximum in %:'
        read *,amp
        write(4,*) 'Width: ',width,', amp:',amp,'%'
        amp=0.01*amp    ! convert from %
        print *,'Give min, max longitude and spacing (deg):'
        read *,xlon1,xlon2,dlon
        print *,'Give min,max latitude and spacing :'
        read *,xlat1,xlat2,dlat
        print *,'Give min,max radius and spacing (km):'
        read *,r1,r2,dr
        write(4,*) 'Longitudes:',xlon1,xlon2,dlon
        write(4,*) 'Latitudes:',xlat1,xlat2,dlat
        write(4,*) 'Radii:',r1,r2,dr

        rlon1=d2r*xlon1
        rlon2=d2r*xlon2
        drlon=dlon*d2r
        drlat=dlat*d2r
        rcolat1=d2r*(90.-xlat2)
        rcolat2=d2r*(90.-xlat1)

        rsign=-1
        r=r2+dr ! from surface downward
        do while(r.ge.r1)
          r=r-dr
          rsign=-rsign
          ysign= rsign
          rcolat=rcolat1-drlat
          do while(rcolat.le.rcolat2)
            rcolat=rcolat+drlat
            ysign=-ysign
            xsign=ysign
            rlon=rlon1-drlon
            do while(rlon.le.rlon2)
              rlon=rlon+drlon
              xsign=-xsign
              x0=r*cos(rlon)*sin(rcolat)
              y0=r*sin(rlon)*sin(rcolat)
              z0=r*cos(rcolat)

c              write(6,1111) xsign, 6371.-r, 90.-(rcolat/d2r), rlon/d2r   !!
c 1111         format(f3.0,2x,f12.2,2x,f12.2,2x,f12.2)                    !!

              do i=1,np       ! loop over all nodes
                xp=points(1,i)
                yp=points(2,i)
                zp=points(3,i)

                call xyz2rtp(xp,yp,zp,rp,theta,phi) ! in radians, -pi,...,pi
cc                if(phi.lt.0) phi = phi + 2*pi       !! if input lon [0,360]

                d1 = abs(r-rp)
                d2 = abs(rcolat-theta)
                d3 = abs(rlon-phi) !!

                if(d1.le.dr/2.and.d2.le.drlat/2.and.d3.le.drlon/2)then
                  sol(i) = xsign*amp
                endif

              enddo
            enddo
          enddo
        enddo

      !-----------------------------------
      ! 7: Isolated Plumes, having a 2-D Gaussian as cross-section.
      !-----------------------------------
      else if(koption.eq.7) then
        print *,'Specify locations, depth extent, and radius of plumes'
        print *,'(radius = 1/e width of Gaussian shaped cross-section)'
        print *,'Give 1/e radius of plume conduit in km:'
        read *,width
        w2=width**2
        dmax=3.0*width
        dmax2=dmax**2
        print *,'Give amplitude of maximum anomaly in %:'
        read *,amp
        write(4,*) 'Plume radius: ',width,', amp:',amp,'%'
        amp=0.01*amp    ! convert from %
        print *,'After last plume, type -1 -1 -1 -1 '
        print *,'Give lon (deg), lat (deg), min_depth, max_depth (km):'
        write(4,'("lon     lat  min_depth  max_depth")')
320     read *,xlon,xlat,depthmin,depthmax
        if(depthmin.le.-1.) goto 400 ! done
        write(4,'(2f8.2,2f6.0)') xlon,xlat,depthmin,depthmax
        r_dmin = 6371.-depthmin
        r_dmax = 6371.-depthmax

        ! convert plume coordinates to cartesian
        rlon=d2r*xlon
        rcolat=d2r*(90.-xlat)
        do i=1,np       ! loop over all mesh vertices
          xp=points(1,i) 
          yp=points(2,i) 
          zp=points(3,i) 
          
          ! find closest point on plume axis (x0,y0,z0)
          rp = sqrt(xp**2+yp**2+zp**2) ! dist from earth center
          r0 = rp ! closest point on axis has same depth
          ! if within the desired plume depth range and distance, then add
          if(r0.le.r_dmin.and.r0.ge.r_dmax)then
            x0=r0*cos(rlon)*sin(rcolat)
            y0=r0*sin(rlon)*sin(rcolat)
            z0=r0*cos(rcolat)

            d2=(xp-x0)**2+(yp-y0)**2+(zp-z0)**2
            if(d2.lt.dmax2) then
              sol(i)=sol(i)+amp*exp(-d2/w2)
            endif
          endif
        enddo
        goto 320

      !-----------------------------------
      ! 8: Isolated Plumes, having a 2-D Gaussian as cross-section.
      !-----------------------------------
      else if(koption.eq.8) then
        print *,'Specify locations, depth extent, and radius of plumes'
        print *,'(radius = 1/e width of Gaussian shaped cross-section)'          
        print *,'Give amplitude of maximum anomaly in %:'
        read *,amp
        write(4,*) 'Plume radius: ',width,', amp:',amp,'%'
        amp=0.01*amp    ! convert from %
        print *,'After last plume, type -1 -1 -1 -1 -1'
        print *,'Give lon (deg), lat (deg), min_depth,'
        print *,'max_depth(km), width (km):'
        write(4,'("lon     lat  min_depth  max_depth  width")')
321     read *,xlon,xlat,depthmin,depthmax, width
        w2=width**2
        dmax=3.0*width
        dmax2=dmax**2
        if(depthmin.le.-1.) goto 400 ! done
        write(4,'(2f8.2,2f6.0)') xlon,xlat,depthmin,depthmax
        r_dmin = 6371.-depthmin
        r_dmax = 6371.-depthmax

        ! convert plume coordinates to cartesian
        rlon=d2r*xlon
        rcolat=d2r*(90.-xlat)
        do i=1,np       ! loop over all mesh vertices
          xp=points(1,i) 
          yp=points(2,i) 
          zp=points(3,i) 
        
          ! find closest point on plume axis (x0,y0,z0)
          rp = sqrt(xp**2+yp**2+zp**2) ! dist from earth center
          r0 = rp ! closest point on axis has same depth
          ! if within the desired plume depth range and distance, then add
          if(r0.le.r_dmin.and.r0.ge.r_dmax)then
            x0=r0*cos(rlon)*sin(rcolat)
            y0=r0*sin(rlon)*sin(rcolat)
            z0=r0*cos(rcolat)

            d2=(xp-x0)**2+(yp-y0)**2+(zp-z0)**2
            if(d2.lt.dmax2) then
              sol(i)=sol(i)+amp*exp(-d2/w2)
            endif
          endif
        enddo
        goto 321
        
      !-----------------------------------
      ! 9: Isolated Plumes, having a 2-D Gaussian as cross-section.
      !-----------------------------------
      else if(koption.eq.9) then
        print *,'Specify locations, depth extent, and radius of plumes'
        print *,'After last plume, type -1 -1 -1 -1 -1'
        print *,'Give lon (deg), lat (deg), min_depth,'
        print *,'max_depth(km), width (km):'
        write(4,'("lon     lat  min_depth  max_depth  width amp")')
322     read *,xlon,xlat,depthmin, depthmax, width, amp
        write(4,*) 'Plume radius: ',width,', amp:',amp,'%'
        amp=0.01*amp    ! convert from %
        w2=width**2
        dmax=3.0*width
        dmax2=dmax**2
        if(depthmin.le.-1.) goto 400 ! done
        write(4,'(2f8.2,2f6.0)') xlon,xlat,depthmin,depthmax
        r_dmin = 6371.-depthmin
        r_dmax = 6371.-depthmax

        ! convert plume coordinates to cartesian
        rlon=d2r*xlon
        rcolat=d2r*(90.-xlat)
        do i=1,np       ! loop over all mesh vertices
          xp=points(1,i) 
          yp=points(2,i) 
          zp=points(3,i) 
        
          ! find closest point on plume axis (x0,y0,z0)
          rp = sqrt(xp**2+yp**2+zp**2) ! dist from earth center
          r0 = rp ! closest point on axis has same depth
          ! if within the desired plume depth range and distance, then add
          if(r0.le.r_dmin.and.r0.ge.r_dmax)then
            x0=r0*cos(rlon)*sin(rcolat)
            y0=r0*sin(rlon)*sin(rcolat)
            z0=r0*cos(rcolat)

            d2=(xp-x0)**2+(yp-y0)**2+(zp-z0)**2
            if(d2.lt.dmax2) then
              sol(i)=sol(i)+amp*exp(-d2/w2)
            endif
          endif
        enddo
        goto 322

      !-----------------------------------
      ! 10: The MOR snake
      !-----------------------------------
      else if(koption.eq.10) then
        print *,'Specify depth extent and radius of MOR columns'
        print *,'(radius = 1/e width of Gaussian shaped cross-section)'
        print *,'Give 1/e radius of MOR column in km:'
        read *,width
        w2=width**2
        dmax=3.0*width
        dmax2=dmax**2
        print *,'Give amplitude of maximum anomaly in %:'
        read *,amp
        print *,'Give depth extend (depthmin & depthmax) in km:'
        read *, depthmin, depthmax
        write(4,*) 'MOR radius: ',width,', amp:',amp,'%'
        amp=0.01*amp    ! convert from %
        
        print *,'Read x_mor, y_mor, z_mor from file'
        open (unit=444, file='plates.csv', STATUS = 'OLD')
        
        do j=1,6049 ! length of plates file
            read(444,*) x_mor, y_mor, z_mor, something
            ! print *, j, x_mor, y_mor, z_mor
            write(4,'("x_mor    y_mor    z_mor")')
            if(something.le.-1.) goto 400 ! done
            write(4,'(2f8.2,2f6.0)') xlon,xlat,depthmin,depthmax
            r_dmin = 6371.-depthmin
            r_dmax = 6371.-depthmax

            do i=1,np       ! loop over all mesh vertices
              xp=points(1,i) 
              yp=points(2,i) 
              zp=points(3,i) 
            
              ! find closest point on plume axis (x0,y0,z0)
              rp = sqrt(xp**2+yp**2+zp**2) ! dist from earth center
              r0 = rp ! closest point on axis has same depth
              ! if within the desired plume depth range and distance, then add
              if(r0.le.r_dmin.and.r0.ge.r_dmax)then
                x0=x_mor*6371.
                y0=y_mor*6371.
                z0=z_mor*6371.

                d2=(xp-x0)**2+(yp-y0)**2+(zp-z0)**2
                if(d2.lt.dmax2) then
                  sol(i)=sol(i)+amp*exp(-d2/2/w2)
                  !sol(i)=sol(i)+amp
                endif
              endif
            enddo
        enddo 

      !-----------------------------------
      ! 11: resolution test from user provided input model
      !-----------------------------------
      else if(koption.eq.11) then
        print *,'Give model file name (e.g. solx.PRI-P05):'
        read(5,fmt='(a)') fname
        open(8,file=fname)
        read(8,*)
        read(8,*) npm
        if(npm.ne.np) stop 'model dimension incompatible with matrix'
        smax=0.
        do i=1,np
          read(8,*) xp,yp,zp,sol(i)
        enddo


      !-----------------------------------
      ! 111: resolution test from user provided input model
      !-----------------------------------
      else if(koption.eq.111) then
        print *,'Give model file name (e.g. 2500_2900_3.txt):'
        read(5,fmt='(a)') fname
        open(8,file=fname)
        read(8,*) npm
        smax=0.

222     read(8,*) r, xlon, xlat, ampl
        if(r.le.0.) goto 400
        write(4,'(f6.0,2f8.2)') r, xlon, xlat
        ! find closest node in grid
        rlon=d2r*xlon
        rcolat=d2r*(90.-xlat)

        x0=r*cos(rlon)*sin(rcolat)
        y0=r*sin(rlon)*sin(rcolat)
        z0=r*cos(rcolat)

        dmax2=1.0e30
        iclosest=1

        do i=1,np       ! loop over all nodes
          xp=points(1,i)
          yp=points(2,i)
          zp=points(3,i)
          d2=(xp-x0)**2+(yp-y0)**2+(zp-z0)**2
          if(d2.lt.dmax2) then
            iclosest=i
            dmax2=d2
          endif
        enddo
        sol(iclosest)=ampl
        goto 222

      !-----------------------------------
      ! 12: resolution test from user provided input model
      !     MINUS zero valued gaussian anomaly
      !-----------------------------------
      else if(koption.eq.12) then
        print *,'Give model file name (e.g. solx.PRI-P05):'
        read(5,fmt='(a)') fname
        open(8,file=fname)
        read(8,*)
        read(8,*) npm
        if(npm.ne.np) stop 'model dimension incompatible with matrix'
        smax=0.
        do i=1,np
          read(8,*) xp,yp,zp,sol(i)
        enddo

        print *,'Give width, amp, radius(km), lon, lat(deg):'
        print *,'After last peak type 0,0,0,0,0'
201     read *, width, amp, r, xlon, xlat
        if(r.le.0.) goto 400
        w2=width**2
        dmax=3.0*width
        dmax2=dmax**2
        rlon=d2r*xlon
        rcolat=d2r*(90.-xlat)
        x0=r*cos(rlon)*sin(rcolat)
        y0=r*sin(rlon)*sin(rcolat)
        z0=r*cos(rcolat)
        do i=1,np       ! loop over all nodes
          xp=points(1,i)
          yp=points(2,i)
          zp=points(3,i)
          d2=(xp-x0)**2+(yp-y0)**2+(zp-z0)**2
          if(d2.lt.dmax2) then
            sol(i)= amp*exp(-d2/w2)
          endif
        enddo
        goto 201

      !-----------------------------------
      ! 13: resolution test from user provided input model
      !      MINUS band defined by rmin/rmax valued zero
      !-----------------------------------
      else if(koption.eq.13) then
        print *,'Give model file name (e.g. solx.PRI-P05):'
        read(5,fmt='(a)') fname
        open(8,file=fname)
        read(8,*)
        read(8,*) npm
        if(npm.ne.np) stop 'model dimension incompatible with matrix'
        smax=0.
        do i=1,np
          read(8,*) xp,yp,zp,sol(i)
        enddo

        ! at this moment we only read ONE entry
        print *, 'Give rmin and rmax for zero out band around'
        print *, 'input model: (eg: 5961, 6371 )'
        read *, rmin, rmax
        print *, rmin, rmax
        ! amp is set to zero! hard coded
        amp = 0
        dmax2=1.0e30

        do i=1,np       ! loop over all mesh vertices
          xp=points(1,i)
          yp=points(2,i)
          zp=points(3,i)
          ! find closest point on plume axis (x0,y0,z0)
          rp = sqrt(xp**2+yp**2+zp**2) ! dist from earth center
          r0 = rp ! closest point on axis has same depth
          ! if within the desired depth range
          ! then put to zero
          if(r0.ge.rmin.and.r0.le.rmax)then
              sol(i)=amp
            endif
        enddo

      !-----------------------------------
      ! 14: resolution test from user provided input model
      !    MINUS all low/high velocities
      !-----------------------------------
      else if(koption.eq.14) then
        print *,'Give model file name (e.g. solx.PRI-P05):'
        read(5,fmt='(a)') fname
        print *, fname
        open(8,file=fname)
        read(8,*)
        read(8,*) npm
        if(npm.ne.np) stop 'model dimension incompatible with matrix'
        smax=0.
        do i=1,np
          read(8,*) xp,yp,zp,sol(i)
        enddo

        ! at this moment we only read ONE entry
        print *, 'What velocities shall I remove? >high< or >low<'
        read *, remval
        print *, remval
        ! amp is set to zero! hard coded
        amp = 0
        dmax2=1.0e30

        do i=1,np       ! loop over all mesh vertices
          xp=points(1,i)
          yp=points(2,i)
          zp=points(3,i)
          ! find closest point on plume axis (x0,y0,z0)
          rp = sqrt(xp**2+yp**2+zp**2) ! dist from earth center
          r0 = rp ! closest point on axis has same depth
          ! if within the desired depth range
          ! then put to zero
          if (remval.eq.'low') then
              if(r0.ge.rmin.and.sol(i).le.0)then
                  sol(i)=amp
                endif
          else if (remval.eq.'high') then
              if(r0.ge.rmin.and.sol(i).ge.0)then
                  sol(i)=amp
                endif
          else
              stop 'I did not implement this. EXIT.'
          endif
        enddo


      !-----------------------------------
      ! 15: get dt from a user provided model
      !-----------------------------------
      else if(koption.eq.15) then
        print *,'Give model file name (e.g. solx.PRI-P05):'
        read(5,fmt='(a)') fname
        open(8,file=fname)
        read(8,*)
        read(8,*) npm
        if(npm.ne.np) stop 'model dimension incompatible with matrix'
        smax=0.
        do i=1,np
          read(8,*) xp,yp,zp,sol(i)
          smax=max(smax,abs(sol(i)))
        enddo  
      
        ! Now scale the model with a priori uncertainties
        write(4,*) 'Model scaling'
        write(4,*) '    i1    i2       dpar'
        do k=1,3
          write(4,'(2i6,f10.3)')jlast(k)+1,jlast(k+1),dpar(k)
          do i=jlast(k)+1,jlast(k+1)
            sol(i)=sol(i)/dpar(k)
          enddo
        enddo  
  
        call asol(sol,rhs,nrow,ncol)      ! compute dT, store in rhs
        print *,'Largest model element =',smax

        open(9,file='times.'//ident)      ! open output file
        print *,'Output is in file: times.'//ident
        if(kssw(1).eq.0) then
          write(9,*) 'no ellipticity correction'
        else if(kssw(1).eq.1) then
          write(9,*) 'ellipticity correction applied'
        else
          write(9,*) 'ellipticity correction unknown'
        endif  
        if(kssw(2).eq.0) then
          write(9,*) 'No crustal correction'
        else if(kssw(1).eq.1) then
          write(9,*) 'crustal correction applied'
        else
          write(9,*) 'crustal correction unknown'
        endif  
        if(kssw(3).eq.0) then
          write(9,*) 'No station elevation correction'
        else if(kssw(1).eq.1) then
          write(9,*) 'station elevation correction applied'
        else
          write(9,*) 'station elevation correction unknown'
        endif  
        if(kssw(4).eq.0) then
          write(9,*) 'No dispersion correction'
        else if(kssw(1).ge.1) then
          write(9,*) 'dispersion correction applied'
        else
          write(9,*) 'dispersion correction unknown'
        endif  
        write(9,'(a,a)') '       i    Tbckgr        dT       Sum',
     &     '  Band   Evt Klust  Stat'

        ! the matrix elements were scaled by grpw*dpar/sigma in assemblematrix
        ! we have already divided sol by dpar, now we multiply the resulting
        ! right hand side by sigma to get an unscaled delay, and output:

        do i=1,nrow
          rhs(i)=rhs(i)*sigma(i)/grpw(i)
          write(9,310) i,dtheor(i),rhs(i),dtheor(i)+rhs(i),iband(i),
     &        ievt(i),kklust(i),kstat(i),stationcode(i)
310       format(i8,3f10.2,4i6,1x,a16)
        enddo
        stop 'Normal end of program for koption=11'
      else
        stop 'koption not allowed'
      endif  


400   continue

      ! write out the synthetic model before scaling with dpar
      print *,'Give ident for synthetic model output file syn.<ident>:'
      read(*,fmt='(a)') ident3
      open(1,file='syn.'//ident3)
      write(6,fmt='(a)') 'Writing model to file syn.'//ident3
      write(1,*)3
      write(1,*) np
      solmax=0.
      solmin=0.
      nsol=0
      do i=1,np  
        write(1,*) (points(j,i),j=1,3),sol(i)
        solmax=max(sol(i),solmax)
        solmin=min(sol(i),solmin)
        if(sol(i).ne.0.) nsol=nsol+1
      enddo  
      close(1)
      print *,'Model min/max values are:',solmin,solmax
      print *,'Nonzero matrix elements:',nsol,' out of ',np
      write(4,*) ' Model min/max values are:',solmin,solmax
      write(4,*) ' Nonzero matrix elements:',nsol,' out of ',np
      
      ! Now scale the model with a priori uncertainties
      write(4,*) 'Model scaling'
      write(4,*) '    i1    i2       dpar'
      do k=1,3
        write(4,'(2i6,f10.3)')jlast(k)+1,jlast(k+1),dpar(k)
        do i=jlast(k)+1,jlast(k+1)
          sol(i)=sol(i)/dpar(k)
        enddo
      enddo  

      call asol(sol,rhs,nrow,ncol)      ! replace rhs by synthetic data
      ! Note K.S.: rhs = kernel_matrix*sol, where model vector sol is
      ! already scaled by prior parameter uncertainties. 
      ! rhs (in the case of real data) would correspond to the contents
      ! of the aux.<ident> file, i.e. data normalized by their 
      ! measurement uncertainties.

      rhsmax=0.
      rhsmin=0.
      nrhs=0
      do i=1,nrow
        rhsmax=max(rhs(i),rhsmax)
        rhsmin=min(rhs(i),rhsmin)
        if(rhs(i).ne.0.) nrhs=nrhs+1
      enddo  
      print *,'Before adding noise:'
      print *,'Data min/max values are:',rhsmin,rhsmax
      print *,'Nonzero data:',nrhs,' out of ',nrow
      write(4,*) ' Before adding noise:'
      write(4,*) ' Data min/max values are:',rhsmin,rhsmax
      write(4,*) ' Nonzero data:',nrhs,' out of ',nrow


      !-----------------------------
      !   ADD NOISE
      !-----------------------------

      print *,'Noise options:'
      print *,'0: no noise'
      print *,'1: Gaussian noise'
      print *,'3: noise from standard deviation given in aux.<indent>'
      if(koption.le.0) 
     &    print *,'2: scrambled data acts as noise (for very low S/N)'
      print *,'Give noise option:'
      read *, noption
      print *, 'Noise option ', noption

      if(koption.eq.0.and.noption.eq.0) stop 'Synthetic data zero!'

600   idum=-1   ! initializes ran1 or gasdev

      if(noption.eq.0) then
        goto 700

      ! Add Gaussian noise
      else if(noption.eq.1) then
        ! the program will deal with the group weighting, so
        ! the correct input for noise should normally be 1
        print *,'Give scaled standard deviation (normally 1):'
        read *,stdev
        write(4,*) 'Noise stand.dev.:',stdev
        ermax=0.
        ermin=0.
        rms=0.
        open(123, file='check_noise_ouput.txt')
        do i=1,nrow
          error=stdev*gasdev(idum)*grpw(i)
          write(123, *) stdev, gasdev(idum), idum, grpw(i), 
     & error , '|', rhs(i)
          rhs(i)=rhs(i)+ error
          write(123, *) 'new rhs(i) = ', rhs(i)
          ermax=max(ermax,error)
          ermin=min(ermin,error)
          rms=rms+error**2
        enddo  
        rms=sqrt(rms/nrow)
        print *,'Min/max error:',ermin,ermax,', RMS=',rms

        rhsmax=0.
        rhsmin=0.
        nrhs=0
        do i=1,nrow
          rhsmax=max(rhs(i),rhsmax)
          rhsmin=min(rhs(i),rhsmin)
          if(rhs(i).ne.0.) nrhs=nrhs+1
        enddo  
        print *,'After adding noise:'
        print *,'Data min/max values are:',rhsmin,rhsmax
        print *,'Nonzero:,',nrhs,' out of ',nrow
        write(4,*) ' After adding noise:'
        write(4,*) ' Data min/max values are:',rhsmin,rhsmax
        write(4,*) ' Nonzero:,',nrhs,' out of ',nrow
        
      ! Scramble data
      else if(noption.eq.2) then
        if(koption.gt.0) stop 'noise and model option not allowed'
        print *,'We assume a fraction of the data is noise'
        print *,'Give fraction (e.g. 1):'
        read *,fraction
        write(4,*) 'N/S ratio:',fraction
        ! exchange each datum with a random other (this assumes the
        ! probability distribution of the data is that of the noise)
        do i=1,nrow
          new=max(1,min(nrow,(nint(nrow*ran1(idum)+1.0))))
          r=rhs(new)
          rhs(new)=rhs(i)
          rhs(i)=r
        enddo
        
      ! stdev/sigma from aux file 
      else if(noption.eq.3) then
          ! the program will deal with the group weighting, so
          ! the correct input for noise should normally be 1
          print *,'Reading stdev from aux file!'
          ermax=0.
          ermin=0.
          rms=0.
          open(123, file='check_noise_ouput.txt')
          do i=1,nrow
            error=sigma(i)*gasdev(idum)*grpw(i)
            write(123, *) sigma(i), gasdev(idum), idum, grpw(i), 
     & error , '|', rhs(i)
            rhs(i)=rhs(i)+ error
            write(123, *) 'new rhs(i) = ', rhs(i)
            ermax=max(ermax,error)
            ermin=min(ermin,error)
            rms=rms+error**2
          enddo  
          rms=sqrt(rms/nrow)
          print *,'Min/max error:',ermin,ermax,', RMS=',rms

          rhsmax=0.
          rhsmin=0.
          nrhs=0
          do i=1,nrow
            rhsmax=max(rhs(i),rhsmax)
            rhsmin=min(rhs(i),rhsmin)
            if(rhs(i).ne.0.) nrhs=nrhs+1
          enddo  
          print *,'After adding noise:'
          print *,'Data min/max values are:',rhsmin,rhsmax
          print *,'Nonzero:,',nrhs,' out of ',nrow
          write(4,*) ' After adding noise:'
          write(4,*) ' Data min/max values are:',rhsmin,rhsmax
          write(4,*) ' Nonzero:,',nrhs,' out of ',nrow
      
      ! reading an arbitary txt file which has the same amount
      ! of rows as the aux file!
      else if(noption.eq.4) then
        ! the program will deal with the group weighting, so
        ! the correct input for noise should normally be 1
        print *,'Reading error from text file!'
        print *,'Give text file name/path:'
        read(5,fmt='(a)') directory
        print *, directory
        open(321,file=directory)
        do i = 1, nrow
            read(321,*) gerror(i)
        end do
        ermax=0.
        ermin=0.
        rms=0.
        open(123, file='check_noise_ouput.txt')
        do i=1,nrow
            error=gerror(i)
            write(123, *) error , '|', rhs(i)
            rhs(i)=rhs(i)+ error
            write(123, *) 'new rhs(i) = ', rhs(i)
            ermax=max(ermax,error)
            ermin=min(ermin,error)
            rms=rms+error**2
        enddo  
        rms=sqrt(rms/nrow)
        print *,'Min/max error:',ermin,ermax,', RMS=',rms
                
        rhsmax=0.
        rhsmin=0.
        nrhs=0
        do i=1,nrow
          rhsmax=max(rhs(i),rhsmax)
          rhsmin=min(rhs(i),rhsmin)
          if(rhs(i).ne.0.) nrhs=nrhs+1
        enddo  
        print *,'After adding noise:'
        print *,'Data min/max values are:',rhsmin,rhsmax
        print *,'Nonzero:,',nrhs,' out of ',nrow
        write(4,*) ' After adding noise:'
        write(4,*) ' Data min/max values are:',rhsmin,rhsmax
        write(4,*) ' Nonzero:,',nrhs,' out of ',nrow
  
      else
        stop 'noption not allowed'
      endif  

700   continue
      print *,'Give ident for new aux (data) output file:'
      read(*,fmt='(a)') ident2

c matrix file
! open the synthetic data auxiliary file (binary format only -- no, ascii)
      print *,'New auxiliary (data) file is: ','aux.'//ident2
      write(4,*) ' New auxiliary (data) file is: ','aux.'//ident2
      open(2,file='aux.'//ident2)

      krtyp=0           ! we don't bother to store type & kluster
      klust1=0
      networkcode='QQ'  ! not stored, give nonsense

      do i=1,nrow
        write(2,120) i,knt(i),rhs(i),sigma(i),ievt(i),
     &      kklust(i),stationcode(i),networkcode,kstat(i),
     &      iband(i),krtyp,klust1,igrp(i),grpw(i)
      enddo
      write(2,120) 0,0,0.,0.,0,0,'X','X',0,0,0,0,0,0.
120   format(2i8,e14.5,e12.3,2i8,1x,a16,1x,a8,2i8,3i7,e14.5)  ! or 2e14.5??
!      write(2,120) 0,0,0.,0.,0,0,'X','X',0,0,0,0,0,0.
!120   format(2i8,e14.5,e12.3,2i8,1x,a16,1x,a8,2i8,3i4,e14.5)
      write(2,810) nstat,(i,statlist(i),netwlist(i),i=1,nstat)
810   format(i5,/,(i5,1x,a16,1x,a8))
      write(2,820) nklust,(i,evlist(i),klustlist(i),i=1,nklust)
820   format(i5,/,(i5,2i10))

      end

      ! Note: Format 120 (for aux.<syn_ident> file) changed in assemblematrix and 
      !    mpisolvetomo on 2010/08/26. Changed it here accordingly,
      !    but the second entry of 2e14.5 does not actually get written
      !    in the case of synthetic data??

! ----------------------------------------------------------------

      subroutine asol(sol,b,m,n)		
      
      ! Forward prediction (matrix multiplication)
      ! computes b=a*sol
      ! where vector b has dimension m*1 and sol n*1, but where
      ! the data prediction b is only computed for b(1)-b(nrow)

      parameter(MMAX=200000)
      dimension a(MMAX),ja(MMAX)
      dimension sol(n),b(m)

      jdebug=0

      ! binary format only
      ! this is the mat.<ident> assembled kernel matrix file,
      ! whose header contains some dimensional info (ncol,krow, etc.)
      rewind 3
      read(3) nrow,ncol
      read(3)                 ! not used
      read(3)                 ! not used
      read(3)
      if(ncol.ne.n) stop 'ncol on file not equal to n'

      amax=0.
      bmax=0.
     

      do irow=1,nrow 				! work row by row

        b(irow)=0.
        
        read(3) krow,kount,iband,ddif,s,kklust,kstat,groupweight
        read(3) (ja(i),i=1,kount)
        read(3) (a(i),i=1,kount)
        if(jdebug.gt.0) write(13,*) 'irow=',irow,' kount=',kount
        
        sum=0
        do j=1,kount
          b(irow)=b(irow)+a(j)*sol(ja(j))  
          sum=sum+a(j)
          if(jdebug.gt.0) write(13,*) ja(j),a(j),sol(ja(j))
          amax=max(amax,abs(a(j)))
        end do
        bmax=max(bmax,abs(b(irow)))

        if(jdebug.gt.0) write(13,*) 'b(',irow,')=',b(irow),
     &     ' sum Aij=',sum
        if(irow.gt.1000) jdebug=0       ! don't print too much

      end do

      print *,'Statistics:'
      print *,'Largest element of matrix = ',amax
      print *,'Largest predicted datum = ',bmax

      return 
      end
        

      FUNCTION ran1(idum)

c random numbers between 0 and 1 (exclusive of the endpoints)
c start with idum=-1

      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END

      FUNCTION GASDEV(IDUM)
      DATA ISET/0/
      IF (ISET.EQ.0) THEN
1       V1=2.*RAN1(IDUM)-1.
        V2=2.*RAN1(IDUM)-1.
        R=V1**2+V2**2
        IF(R.GE.1.)GO TO 1
        FAC=SQRT(-2.*LOG(R)/R)
        GSET=V1*FAC
        GASDEV=V2*FAC
        ISET=1
      ELSE
        GASDEV=GSET
        ISET=0
      ENDIF
      RETURN
      END

      FUNCTION PLGNDR(L,M,X)
      IF(M.LT.0.OR.M.GT.L.OR.ABS(X).GT.1.)stop 'plgndr: bad arguments'
      PMM=1.
      IF(M.GT.0) THEN
        SOMX2=SQRT((1.-X)*(1.+X))
        FACT=1.
        DO 11 I=1,M
          PMM=-PMM*FACT*SOMX2
          FACT=FACT+2.
11      CONTINUE
      ENDIF
      IF(L.EQ.M) THEN
        PLGNDR=PMM
      ELSE
        PMMP1=X*(2*M+1)*PMM
        IF(L.EQ.M+1) THEN
          PLGNDR=PMMP1
        ELSE
          DO 12 LL=M+2,L
            PLL=(X*(2*LL-1)*PMMP1-(LL+M-1)*PMM)/(LL-M)
            PMM=PMMP1
            PMMP1=PLL
12        CONTINUE
          PLGNDR=PLL
        ENDIF
      ENDIF
      RETURN
      END

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
