        program raymatrix

      ! compile g77 -o raymatrix raymatrix.f  

      ! Changed from Raffaella Montelli's globmat.f & bldmat.f to
      ! accommodate general ray definitions and user-specified
      ! background models (as in raydata.f), as well as multiple
      ! frequency bands. Added amplitude data and changed all variables
      ! to single precision except for Delauney routines.
      ! Made it more efficient where possible.

      ! INPUT:
      ! Interactive or piped parameter input; Example:
      !
      ! 1 0 0    ! on/off switches for model parameters Vp, Vs, Qs
      ! 20.0     ! typical step size (km) for quadrature of kernels
      ! newvertices.xyz     ! vertex file name
      ! newfacets.all       ! facet file name
      ! P001_WUSA.P         ! ident (identifier for input AND output
      !                       file names. infile = raydata.<ident>)

      ! OUTPUT:
      ! 
      ! raymatrix.out.<ident>     log of kernel info written to binaries.
      ! raymatrix.info.T.<ident>  log of data info written to binaries.
      ! raymatrix.info.A.<ident>  Files get read by assemblematrix masking scheme.
      !                         
      ! (11) matrixT.<ident>         binary, ttime kernels for each datum
      ! (12) matrixA.<ident>         binary, amplitude kernels for each datum
      ! histo.xy.<ident>        histogram of epicentral distances
      ! matx:                   debug plot output
      
      ! Changes
      ! Sept 2006: In subroutine dintw, xclm is pre-read and frequency dependent 
      !  Make dtmax and spatial integral (drho, dtheta) frequency dependent
      !  (subroutines getk and filla)
      ! Oct 2006:  Deal with negative tdetour (subroutine getk)
      ! Treat the region near source/receiver & in crust as a homogeneous 
      ! sphere and use ray theory within the sphere (after "Handle source first")
      ! Output warning message of imprecision of kernel integral test if 
      ! jdebug>0 (after rho loop and iband loop).
      ! Compute period() in subroutine rdband2 and output period() after np,jssw

      ! *********** Experimental version under testing ************
      ! Look for TODO to find programming incompletes
      ! Does not yet handle:
      !   accurate time windowing
      !   backscatter for pP,PP
      !   later arrivals (?)
      !   linked data
      ! GN 2006.

      ! raymatrix is the third step in the tomographic inversion:
      ! [1] run raydata to construct a finite-frequency data file
      ! [2] run springs3d.f to parameterize the model
      ! [3] run raymatrix.f to build the raw matrix for ray data.
      ! [4] run assemblematrix to combine all data & raw matrices
      ! [5] run solvetomo to get a tomographic solution

      ! If datafile is linked (for differential data such as PP-P)
      ! the matrix file needs to be de-linked with diffmatrix.f
      ! after running raymatrix.f; this creates a matrix for the
      ! arrival time *differences*.

      ! If no other data are to be added, raymatrix will be followed
      ! by a run of tomoinversion.f

      ! CHANGE LOG
      ! 2010/08/12 K.S. change name of common block from /mod/ to /modl/
      !  since old name was a fortran standard violation and ifort would
      ! not compile. Also changed in raydata.f and raytracesubs.f
      ! raymatrix.f(367): error #8038: A referenced intrinsic procedure 
      ! can not have the same name as a common block.   [MOD]


      ! header files for Delauney triangulation
      include '../includes/nn.param'
      include '../includes/setdel3.h'       ! common for np,nt (Delauney)
      include '../includes/setdel2.h'       ! common for points,centres,vertices etc

c     common/debug/asum

! ray structure information (segments, reflection points etc)

      dimension rseg(20),ktseg(20),kdwn(20)             ! ray segments
      dimension ptlat(20),ptlon(20)                     ! bounce pts
      dimension xsgl(20),ysgl(20),hsgl(20)              ! singular points
      dimension tau(20),telev(20)                       ! corrections
      dimension rotm(3,3)                               ! transformation

! filter bank information (for now assumed same for all rays
! in this raydata file)

      parameter (NFR=20)        ! total (time+ampl) nr of freq bands
      dimension nfreq(NFR),bndomega(200,NFR),dotm(200,NFR) ! spectral bands
      dimension detourlim(NFR),winlen(NFR)  ! effective limit for kernels
      dimension xclm(NFR)       ! length of xcorr window
      dimension wbar(NFR)       ! average frequency of a particular source spectrum
      dimension dtmax(NFR),dstmax(NFR)   ! boundary of spatial integral
      dimension jssw(3)         ! parameter type sense switches
      dimension kountpar(3)     ! individual kounts for Vp,Vs,Qs
      character*2 cpar(3)       ! parameter codes
      data cpar/'Vp','Vs','Qs'/

! We calculate the frequency integrals for time and amplitude
! only once and interpolate over at most 250 detour time points.
! tknl(i,j,k) has the kernel for detour time i and frequency
! band k, with j=1 giving sin(w*dT) and j=2 the cos(w*dT) kernel.

      parameter (NTM=1500)       ! must be conform dintw.f
      dimension tknl(NTM,2,NFR),aknl(NTM,2,NFR),dtbl(NTM),ntbl(NTM)

! matrices for the Delauney tetrahedra
      dimension ab(4,4),abinv(4,4,nt_max)       ! transformations

! coordinate transformation matrices (TODO: remove unused g2e)
      dimension e2g(3,3),g2e(3,3)

! matrix rows aa for total of NFR frequency bands and column number ja
      dimension aa(3*np_max,NFR),asparse(3*np_max),ja(3*np_max)  
      dimension dknl(2,NFR)                    ! temporary kernel storage

! data
! We allow for NFR/2 frequency bands in each of the amplitude and
! time data (or NFR total if not equal). Note that the data for each
! source-receiver path are numbered differently from the numbering
! of the frequency bands. Thus, tobs(i) is observed travel time
! (NOT delay...) in band nbt(i), with error tsig(i) and
! x-correlation coefficient corcoeft(i). Similar for amplitude
! anomalies d lnA in aobs(i). time is the theoretical arrival
! time not yet corrected for crust or ellipticity.
! rmsb(k) is the average noise level in band k. 
! xl and xla are x-correlation window length for time, amplitude 
! (usually equal)

      ! delaytime data
      dimension tobs(NFR),nbt(NFR),tsig(NFR),corcoeft(NFR)
      ! amplitude data
      dimension aobs(NFR),nba(NFR),asig(NFR),corcoefa(NFR)  
      dimension rmsb(NFR),xl(NFR),xla(NFR),nflag(NFR),period(NFR)
      integer histodistance(361)

! ray information
! yp(i,j,k) contains k'th arrival, node j; coordinates i=1:r,
! i=2: angle i, i=3: arc distance delta. h11 and h22 contain
! d^2T/dq_1 and d^2T/dq_2, resp., and rayvel(j) is the velocity
! in the background model at node j.
! we allow for at most three arrivals within the correlation time
! window (e.g. a triplication or S+ScS).

      parameter (NDR=9000,NDL=190)
      dimension yp(3,NDR)            ! ray vector
      dimension h11(NDR),h22(NDR)    ! Hessians M11,M22
      dimension rayvel(NDR),rayq(NDR) ! velocity and 1/Q at ray nodes
      dimension legend(NDL),ktype(NDR)

! various character variables 
      character*72 fname,fname1
      character*16 stationcode,chray
      character*8 phase,netw
      character*30 dataf,dataf1,dataf2
      character*3 comp
      character*1 yesno
      character*1 dtyp
      
      logical ruthere

! background model (qvp,qvs and qqs are ignored)
      common /modl/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &      Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

! homogeneous sphere radius around singular points
      data r0/20./ ! radius of homogeneous sphere at source and receiver
                   ! the same as in raydata
      data r1/20./ ! radius of homogeneous sphere at caustics
      data r2/300./ ! radius of sphere with smaller quadrature stepsize(drho) near caustics

      ! The following allows for debugging. Normally jdebug is 0.
      ! Or set jdebug=-1 if you wish files to be ascii rather than binary.
      jdebug=0      ! set to 1 for limited, 2 for extensive output, 
                     ! to n for lots of output at node n. 

      ! set universal constants at machine precision
      pi=4.0*atan(1.0)
      twopi=2.0*pi
      halfpi=0.5*pi
      r2d=180.0/pi
      do i=1,361
        histodistance(i)=0
      enddo  

      print *,'Read parameters sense switches (Vp,Vs,Qs) 0=no,1=yes:'
      print *,'(eg 1 0 0 for Vp only, or 0 1 1 for Vs & Qs)'
      read *,jssw

      ! set some computational constants; lower step0 for high frequency
      print *,'Give typical step size for quadrature (eg 20 km):'
      read *,step0
      dstlimit=1000.      ! limit radius of Frechet kernel 

! input file with model parametrization
      print *,'Input 3D model grid...'
      call inpgrd
      print *,'3D model has ',np,' grid points'
      np3=np*3
      if(np.gt.np_max) stop 'np>np_max; raise dimensions in nn.param'

! compute inverse matrices for interpolation coefficients
      print *,'Computing coefficient matrices...'
      do it=1,nt        ! loop over nt tetrahedra (nt in setdel3.h)
	x1=points(1,vertices(1,it))    ! x for first vertex
	y1=points(2,vertices(1,it))    ! y for first vertex
	z1=points(3,vertices(1,it))    ! z for first vertex
	do jj=1,4       ! points and vertices are in setdel2.h
	  ab(1,jj)=points(1,vertices(jj,it))-x1    ! dx for point jj
	  ab(2,jj)=points(2,vertices(jj,it))-y1    ! dy for point jj
	  ab(3,jj)=points(3,vertices(jj,it))-z1    ! dz for point jj
	  ab(4,jj)=1
	enddo
        call invrt(ab,abinv(1,1,it),4)
      enddo  

! open the processed data file for input (END OF USER INPUT)
      print *,'Give ident of data file name (eg PP for raydata.PP):'
      read(5,fmt='(a)') dataf

      ! Do not inadvertently overwrite an earlier matrix file
      fname='matrixT.'//dataf
      inquire(file=fname,exist=ruthere)
      if(ruthere) then
        print *,'CAREFUL! an earlier matrix file with this ident'
        print *,'exists. Can I overwrite it (y/n) ?'
        read(*,'(a)') yesno
        if(yesno.ne.'y') stop 'User bailed out'
      endif  
      fname='raydata.'//dataf
      
      print *,'Opening processed data file ',fname
      open(2,file=fname,status='old')

      ! read first lines of data file
      read(2,fmt='(a)') fname
      read(2,*) maxarr,tdifmax
      read(2,fmt='(a)') dataf1
      read(2,fmt='(a)') dataf2
      linked=0
      if(dataf2.ne.'None') linked=1

! open output files
      print *,'Opening file raymatrix.info.T.'//dataf1
      open(17,file='raymatrix.info.T.'//dataf1)
      write(17,89)   
      print *,'Opening file raymatrix.info.A.'//dataf1
      open(18,file='raymatrix.info.A.'//dataf1)
      write(18,89)
89    format( 'idat ievt kluster type band stnam netwk rlat
     &  rlon    relev  delta    azi  xcorr rms_noi    kunit comp idate
     &iotime slat   slon   sdep ')
      
      
      print *,'Opening file raymatrix.out.'//dataf1
      open(1,file='raymatrix.out.'//dataf1)
      itm1=time()        ! nonstandard Fortran
      write(1,fmt='(a)') ctime(itm1)
      write(1,fmt='(2a)') ' Data file: ',dataf1
      write(1,fmt='(2a)') ' Linked as differential data to: ',dataf2
      write(1,fmt='(a,i2,1x,a,f7.1)') ' Data file has maxarr=',maxarr,
     &      ' and tdifmax=',tdifmax
      if(maxarr.gt.0) write(1,*) 'WARNING: later arrivals ignored'

      print *,'Reading 1D background model file: ',fname
      write(1,fmt='(2a)') ' Model file: ',fname
      call model1(fname)
      write(6,*) 'Radius of Earth model:',r(n)
      write(6,*) 'Radius of CMB:',r(noc)
      write(6,*) 'Radius of ICB:',r(nic)
      write(6,*) 'Radius of Moho:',r(nmoh)
      rn=r(n)
      rcmb=r(noc)

! open matrix files
      fname='matrixT.'//dataf
      fname1='matrixA.'//dataf
      if(jdebug.eq.0) then      ! unformatted binary
        print *,'Matrix files are unformatted'
        print *,'Opening matrix file ',fname
        open(11,file=fname,form='unformatted')
        print *,'Opening matrix file ',fname1
        open(12,file=fname1,form='unformatted')
        write(11) linked,dataf2
        write(11) np,jssw
        write(12) linked,dataf2
        write(12) np,jssw
      else  
        print *,'Matrix files are formatted, jdebug=',jdebug
        print *,'Opening matrix file ',fname
        open(11,file=fname)
        print *,'Opening matrix file ',fname1
        open(12,file=fname1)
        write(11,fmt='(i2,1x,a)') linked,dataf2
        write(11,*) np,jssw
        write(12,fmt='(i2,1x,a)') linked,dataf2
        write(12,*) np,jssw
      endif

! read ray info from the data file
      print *,'Reading ray info from processed data file'
      call rdray(2,0,chray,phase,rseg,ktseg,kdwn,nseg)
      print *,'Ray info read for ',chray,' with ',nseg,' segments'

! read band info from the data file 
      print *,'Reading band info from processed data file'
      call rdband2(2,0,bndomega,dotm,nfreq,nband,period)
      print *,nband,' frequency bands read'
      if(nband.gt.NFR/2) stop 'nband exceeds NFR'
      if(jdebug.eq.0) then
        write(11) nband,(period(i),i=1,nband)
        write(12) nband,(period(i),i=1,nband)
      else  
        write(11,fmt='(i3,(20f7.1))') nband,(period(i),i=1,nband)
        write(12,fmt='(i3,(20f7.1))') nband,(period(i),i=1,nband)
      endif  

! read max. xcorr length from the xclm file
      open(8,file='xclm.'//dataf,status='old')
      read(8,fmt='(a72)') line  ! skip the 1st line
      do ib=1,nband
        read(8,fmt='(i4,f10.3)') jb,xclm(ib)
        if(jb.ne.ib) then
          print *, 'Missing band in xclm.'//dataf
          stop
        endif
      enddo
      close(8)

! construct integral tables tknl and aknl and find sensitivity region
      print *,'Computing ',nband,' frequency integral tables'
      call dintw(tknl,aknl,dtbl,ntbl,dotm,bndomega,nfreq,nband,xclm)
      ! find detour limits, adjust ntbl
      call findtlim(tknl,nband,dtbl,ntbl,detourlim)

! average frequency wbar(j) for band j
      call meanw(wbar,dotm,bndomega,nfreq,nband,xclm)  


      print *,'Reading data, every x is 100 data processed'
      ndata=0
      nrow =0     ! number of output matrix rows 
      nrowt=0     ! number of output matrix rows for ttimes
      nrowa=0     ! number of output matrix rows for ttimes
      rmsb(1)=0.  ! safeguard in case of nband=0 (ray theory)
      write(1,90) 
90    format(/,'Sparse matrix winnowing',//,
     &  'Date     Station  B  P D kountp  non0  Theor/A+    sum/A-',
     &  ' sumsparse sparse(%) reduct(%)  accursparse(%) accur(%)')

!     LOOP: read next source-receiver path from raydata.<ident>
10    read(2,20,end=900,iostat=ios) idate,iotime,ievt,kluster,
     &      stationcode,netw,comp,slat,slon,sdep,
     &      rlat,rlon,relev,rdel,raz,del,az,narr
20    format(3i8,i4,1x,a16,1x,a8,1x,a3,2f9.3,f7.1,2f9.3,f7.3,
     &       2f9.5,2f9.3,i3)
21    format(3f9.2,i3,f9.1)
22    format(i2,16f10.1)
      ndata=ndata+1
      if(ios.ne.0.or.idate.le.0) then
        print *,'Error in data file header at datum #',ndata
        print *,'idate,iotime,ievt,stat=',idate,iotime,ievt,
     &    stationcode
        stop
      endif  
      nhist=max(1,min(nint(del),361))
      histodistance(nhist)=histodistance(nhist)+1
c     write(6,25) ndata,ievt,stationcode,comp
      if(mod(ndata,5000).eq.0) write(6,*)'-',ndata
      if(mod(ndata,100).eq.0) write(6,'(a1,$)') 'x'
25    format(i6," Event nr",i8,1x,"in ",a16,1x,a3)
      read(2,22,iostat=ios) kunit,rms0,(rmsb(i),i=1,nband)
      if(ios.ne.0) then
        print *,'Error in data file rms line at datum #',ndata
        print *,'idate,iotime,ievt,stat=',idate,iotime,ievt,
     &    stationcode
        stop
      endif  

      ! set nflag (1=datum for band ib,0=no datum) to zero
      ! winlen(ib) is xcorr window length in band ib for this path
      do ib=1,nband
        nflag(ib)=0
        winlen(ib) = xclm(ib)
      enddo

      ! read travel time in different bands for this path
      read(2,*,iostat=ios) nobst
      if(ios.ne.0) then
        print *,'Error in data file nobst line at datum #',ndata
        print *,'idate,iotime,ievt,stat=',idate,iotime,ievt,
     &    stationcode
        stop
      endif  
      do i=1,nobst
        read(2,21,iostat=ios) tobs(i),tsig(i),corcoeft(i),nbt(i),xl(i)
        if(ios.ne.0) then
          print *,'Error in data file tobs line at datum #',ndata
          print *,'idate,iotime,ievt,stat=',idate,iotime,ievt,
     &      stationcode
          stop
        endif

        if(nbt(i).gt.0) then      ! only for ff data --  nbt(i)=0  means ray theory
           nflag(nbt(i))=1          
           winlen(nbt(i)) = xl(i) ! winlen(ib) is xcorr window length in band ib for this path
        endif
      enddo  

      ! find detour time limit for kernel integraion for this path
      do ib=1,nband
        dtmax(ib) = min(detourlim(ib),winlen(ib))
      enddo

      ! read amplitude anomalies for this path
      read(2,*,iostat=ios) nobsa
      if(ios.ne.0) then
        print *,'Error in data file nobsa line at datum #',ndata
        print *,'idate,iotime,ievt,stat=',idate,iotime,ievt,
     &    stationcode
        stop
      endif  
      if(nobst.gt.0.and.nbt(1).le.0.and.nobsa.gt.0) then
        print *,'nobsa>0 while processing ray theoretical datum'
        print *,'Error occurs in data file at datum #',ndata
        print *,'idate,iotime,ievt,stat=',idate,iotime,ievt,
     &    stationcode
        stop
      endif
      do i=1,nobsa
        read(2,21,iostat=ios) aobs(i),asig(i),corcoefa(i),nba(i),xla(i)
        if(ios.ne.0) then
          print *,'Error in data file aobs line at datum #',ndata
          print *,'idate,iotime,ievt,stat=',idate,iotime,ievt,
     &      stationcode
          stop
        endif  
        if(nba(i).gt.0) nflag(nba(i))=1
      enddo  

      if(narr.le.0) goto 10     ! skip if no arrival

      nbtot=nobst+nobsa         ! total # of frequency bands for A+T
c     asum=0.                   ! debug on sum of Aij

      ! debug of reflections
c     do ib=1,nband
c       dtmin=0.
c       if(linked.gt.0) dtmin=-dtmax(ib)      ! symmetric windowing for PP,pP etc
c       if(jdebug.gt.0.and.ib.eq.1) write(13,*) 'band 1, dtmin,max=',
c    &     dtmin,dtmax(1)
c     enddo


      ! Read reflection points, sum corrections
      if(jdebug.ge.1) write(13,1305)
1305  format('   k  kd    plat    plon      rseg     tau    elev')
      corcrust=0.
      corelev=0.
      do k=1,nseg  
        read(2,35,iostat=ios) kseg,kdwn(k),ptlat(k),ptlon(k),rtarget,
     &        tau(k),telev(k)
35      format(2i3,2f8.2,f8.1,2f8.3)
        if(ios.ne.0) then
          print *,'Error in data file ptlat/lon line at datum #',ndata
          print *,'idate,iotime,ievt,stat=',idate,iotime,ievt,
     &      stationcode
          print *,'Possible problem: source and rcv too closely spaced.'
          print *,'Quitting.'
          stop
        endif  
        ! take care of one bogus data line and one empty line written by
        ! raydata in case of no convergence (could not be handled more elegantly)
        if(kseg.le.0) then
          read(2,*)
          read(2,*)
          goto 10
        endif  
        if(jdebug.ge.1) write(13,1306) k,kdwn(k),ptlat(k),ptlon(k),
     &        rtarget,tau(k),telev(k)
1306    format(2i4,2f8.2,f10.1,2f8.2)
        corcrust=corcrust+tau(k)
        corelev=corelev+telev(k)
      enddo         ! end of loop over ray segments 
      if(jdebug.ge.1) write(13,*) 'Total: ',corcrust,corelev

      ! set matrix row elements to 0
      do j=1,nbtot
        do i=1,np3
          aa(i,j)=0
        enddo
      enddo

      ! Read first arrival for this path
      read(2,45,iostat=ios) jar,nray,nlegs,trtime,ecorr,tstar,qray,
     &      slowness
45    format(i2,2i5,3f10.2,3e12.3)
      if(ios.ne.0) then
        print *,'Error in data file arrival# line at datum #',ndata
        print *,'idate,iotime,ievt,stat=',idate,iotime,ievt,
     &    stationcode
        stop
      endif  
      
      if(nray.gt.NDR) stop 'Too many ray nodes, increase NDR'
      if(nlegs.gt.NDL) stop 'Too many ray legs, increase NDL'
      if(ecorr.lt.-90.) ecorr=0.         ! ecorr=-99 if unavailable
      read(2,*,iostat=ios) (legend(i),i=1,nlegs)
      if(ios.ne.0) then
        print *,'Error in data file leg end line at datum #',ndata
        print *,'idate,iotime,ievt,stat=',idate,iotime,ievt,
     &    stationcode
        stop
      endif  
      if(nray.lt.1) goto 10    ! it should never get here...
      iq=n

      ! read r(km), i(rad), phi(rad), c0(km/s), h11 and h22 (s/km^2)
      ! meanwhile find caustics
      nsgl = 2    ! # of singular points (first 2 are source and receiver)
      h11p = 0.  ! h11 at previous ray node
      do i=1,nray
        read(2,46,iostat=ios) (yp(j,i),j=1,3),rayvel(i),rayq(i),
     &        h11(i),h22(i)
        if(ios.ne.0) then
          print *,'Error in data file h11/h22 line at datum #',ndata
          print *,'idate,iotime,ievt,stat=',idate,iotime,ievt,
     &      stationcode
          stop
        endif  
        ktype(i)=1
        if(rayvel(i).lt.0.) ktype(i)=2  ! Vs is identified by <0
        rayvel(i)=abs(rayvel(i))
46      format(f10.3,3f10.5,f10.6,2e14.6)
        ! find ray node and coordinates of caustics (where h11 changes sign)
        if(h11p*h11(i).lt.0)  then
          nsgl = nsgl+1    ! # of singular point
          h11c = h11(i)
          rsgl = yp(1,i)+h11c*(yp(1,i-1)-yp(1,i))/(h11c-h11p) ! caustic r
          psgl = yp(3,i)+h11c*(yp(3,i-1)-yp(3,i))/(h11c-h11p) ! caustic phi
          xsgl(nsgl) = rsgl*cos(psgl)           ! cartesian coordinates of caustic
          ysgl(nsgl) = rsgl*sin(psgl)
          if(jdebug.gt.0) write(13,*) 'caustic ',nsgl-2,' near ray node'
     &      ,i,', r,phi=',rsgl,psgl
        endif
        h11p = h11(i)
      enddo
      angle0=yp(2,1)    ! take-off angle with vertical (radians, pi=up)
      vsrc=rayvel(1)    ! velocity at source
      krtyp=ktype(nray) ! 1=P or 2=S at the receiver
      ! first two singular points (source and receiver)
      xsgl(1) = yp(1,1)*cos(yp(3,1))  
      ysgl(1) = yp(1,1)*sin(yp(3,1))  
      xsgl(2) = yp(1,nray)*cos(yp(3,nray))
      ysgl(2) = yp(1,nray)*sin(yp(3,nray))

      ! Skip ray information for later arrivals, if any
      do iar=2,narr
        read(2,45,iostat=ios) jar,nray2,nleg2
        read(2,*,iostat=ios2) (leg2,i=1,nleg2)
        if(ios.ne.0.or.ios2.ne.0) then
          print *,'Error while skipping later arrivals, for datum #',
     &          ndata
          print *,'idate,iotime,ievt,stat=',idate,iotime,ievt,
     &    stationcode
          stop
        endif  
        iq=n
        do i=1,nray2
          read(2,46,iostat=ios) 
        if(ios.ne.0.or.ios2.ne.0) then
          print *,'Error while skipping h11/h22 lines, for datum #',
     &          ndata
          print *,'idate,iotime,ievt,stat=',idate,iotime,ievt,
     &    stationcode
          stop
        endif  
        enddo
      enddo
      
      if(jdebug.ge.1) write(13,*) narr,' arrivals, ellipcor=',ecorr

      scolat=(90.-slat)/r2d
      sph=slon/r2d
      rcolat=(90.0-rlat)/r2d
      rph=rlon/r2d
      sth=geocen(scolat)       ! ellipticity corrected latitude
      rth=geocen(rcolat)

      ! compute the coordinate transformation matrices:
      ! e2g to map from equatorial (x,y,0) coordinates to geocentric
      ! g2e to map from geocentric to equatorial cartesian
      call euler_rot(sth,sph,rth,rph,g2e,e2g)
      if(jdebug.gt.0) write(13,*) 'e2g=',e2g

c------------------------------------------------------------------
c construction of the nonzero elements of the matrix A
c------------------------------------------------------------------

      loc=0             ! Delauney node initialization
      ntot=nobst+nobsa
      ileg=1    ! legs separated by jumps in quadrature plane orientation
      nquadtot=0         ! counts number of integration volumes
      if(jdebug.ge.2)
     &write(13,fmt="(4a)")'  k0      diskkt         1/c0   errkt(%)',
     &'    diskpstv     diskngtv   errka(%)',
     &'    diskka       theor      errkq(%)',
     &'      n      alpha    rhomax     dlmin etc'
      x0=yp(1,1)*cos(yp(3,1))    ! Cartesian coordinates of source
      y0=yp(1,1)*sin(yp(3,1))    
      tstarq=0.         ! for debug KQ

      do k0=2,nray     ! loop over ray nodes 
        
        diskint=0.       ! surface integral to check eq (114),dT
        diskpstv=0.      ! surface integral of positive kernel, dA
        diskngtv=0.      ! surface integral of negative kernel, dA
        diskq=0.         ! surface integral of Q kernel, dA

        x1=x0            ! x1,y1 is previous node
        y1=y0
        k1=k0-1
        x0=yp(1,k0)*cos(yp(3,k0))    ! x=r*cos phi
        y0=yp(1,k0)*sin(yp(3,k0))    ! y=r*sin phi
        xc=0.5*(x0+x1)               ! midpoint on ray step   
        yc=0.5*(y0+y1)

        if(k1.ge.legend(ileg)) then  ! reflection
          ileg=ileg+1
          goto 80                    ! skip 
        endif  

        c0=rayvel(k0)              ! Vp or Vs at ray node
        cav=0.5*(rayvel(k1)+c0)    ! average velocity over segment
        q0=rayq(k0)                ! 1/Q at ray node
        qav=0.5*(rayq(k1)+q0)      ! average 1/Q over segment
        ! Use cartesian coordinates in equatorial plane
        h=sqrt((x0-x1)**2+(y0-y1)**2)     ! distance from last node
        dlmin=h
        dlmax=h
c       tsum=tsum+h/cav       ! travel time (for debug)   
        tstarq=tstarq+h*qav/cav

        ! RAY THEORY? Then do this and jump to end of ray node loop:
        if(nobst.gt.0.and.nbt(1).eq.0) then
          dV=h      ! (line integral, this simple trick fools filla)
          dknl(1,1)=-1.0/cav               ! ray kernel for dv/v (only one datum!)
          maslov=0
c         xsg=e2g(1,1)*x0+e2g(1,2)*y0      ! coordinate of ray segment center point
c         ysg=e2g(2,1)*x0+e2g(2,2)*y0
c         zsg=e2g(3,1)*x0+e2g(3,2)*y0
          call xyzscat(x0,y0,0.,0.,rotm,e2g,xsg,ysg,zsg,rsg,tsg,psg)
          if(jdebug.ge.2) write(13,*) k0,'  ray theory: lat,lon,r=',
     &      90-tsg*r2d,psg*r2d,rsg
          call filla(xsg,ysg,zsg,nobst,nobsa,aa,abinv,loc,dV,dknl,
     &      maslov,q0,ktype(k0),jssw,-1,-1)
          goto 85                  ! jump to end of k0 loop
        endif
 
        ! Near singular points (source, receiver, caustics)?
        ! Then use ray theory for T and Q kernels, use constant velocity perturbation
        ! for A kernel, and jump to end of ray node loop
        do js=1,nsgl           ! loop through singular points
          hsgl(js)=sqrt((x0-xsgl(js))**2+(y0-ysgl(js))**2) ! distance to singular point
          ! near singular points (r0 for source/receiver, r1 for caustics)
          if( (js.le.2 .and. hsgl(js).lt.r0) .or. 
     &        (js.gt.2 .and. hsgl(js).lt.r1) ) then  
            do ll=1,nobst
              dknl(1,ll)=-1.0/cav   ! ray kernel for traveltime
              dknl(2,ll)=-1.0/cav  
            enddo
            do ll=1,nobsa
              la=nobst+ll     ! amplitude datum number
              dknl(1,la)=0.   ! constant velocity perturbation --> dA=0
              dknl(2,la)=0.   
              if(jssw(3).eq.1) then      ! Q kernel
                ! find dominant frequency for datum la
                ib=nba(ll)               ! band number of amplitude datum
                qk=-wbar(ib)*qav/2./cav     ! ray kernel for Q
                dknl(1,la)=qk            ! add on top of zero A kernel
                dknl(2,la)=qk
              endif
            enddo
            dV=h      ! (line integral, this simple trick fools filla)
            maslov=0
            call xyzscat(x0,y0,0.,0.,rotm,e2g,xsg,ysg,zsg,rsg,tsg,psg)
            if(jdebug.ge.2) write(13,*) k0,'  ray theory: lat,lon,r=',
     &        90-tsg*r2d,psg*r2d,rsg
            call filla(xsg,ysg,zsg,nobst,nobsa,aa,abinv,loc,dV,dknl,
     &        maslov,q0,ktype(k0),jssw,-1,-1)
            goto 85                  ! jump to end of k0 loop
          endif
        enddo
 
        ! alpha is the angle between subsequent planes
        alpha=yp(3,k1)-yp(3,k0)+yp(2,k1)-yp(2,k0)
        if(alpha.eq.0.) then
          rr=1.0e10
        else  
          rr=h/alpha          ! curvature radius of raypath
        endif  

        if(h11(k0).eq.0.) stop 'h11=0 in ray file'
        hmin=min(abs(h11(k0)),abs(h22(k0)))
        hmax=max(abs(h11(k0)),abs(h22(k0)))
        hdiff=max(0.0001,abs(h11(k0)-h22(k0)))
        ! The following assumes Omega=1 (i.e. forward P->P or S->S
        ! scattering only). For mode conversion or large angle
        ! scattering implement eq. (57) into defenition of spr
        spr=sqrt(abs(h11(k0)*h22(k0)))/c0    ! sqrt(|M'+M"|)/c, eq 104
        call findmas(h11(k0),h22(k0),maslov) ! find Maslov index

        ! find spacial quadrature limit
        do ib=1,nband
          dstmax(ib)=sqrt(2*dtmax(ib)/hmin)      ! max distance from ray
          dstmax(ib)=min(dstmax(ib),dstlimit)    ! but don't go too far
        enddo
        call coordl(yp(2,k0),yp(3,k0),rotm) ! q1,q2 directions

        if(jdebug.gt.2.and.jdebug.eq.k0) then
          write(13,*) '------------------------------- node ',k0
          write(13,*) 'h11,h22=',h11(k0),h22(k0)
          write(13,*) 'phi,i=',yp(3,k0),yp(3,k1),yp(2,k0),yp(2,k1)
          write(13,*) 'h,rr, alpha(deg)=',h,rr,alpha*r2d
          write(13,*) 'hmin,dstlimit=',hmin,dstlimit
          write(13,*) 'c0,q0,maslov=',c0,q0,maslov
          write(13,*) 'dstmax(i)=',(dstmax(i),i=1,nband)
          write(13,*) 'step0,drho(i)=',step0,(0.05*dstmax(i),i=1,nband)
          write(13,*) 'rotm=',rotm
          write(13,1307)
1307      format('bd1',3x,'rhoc',6x,'drho',4x,'dtheta',9x,'n',4x,
     &           'detour',5x,'diskint',4x,'diskpstv',4x,'diskngtv',
     &           4x,'diskq')
        endif

        ! spatial quadrature step size is frequency dependent
        do ibd=1,nband
          if(nflag(ibd).eq.0) goto 65  ! go to next band

          ! do quadrature in polar coordinates to keep it simple
          rho=0.                          ! polar radius
          nquad=0
          drho=min(0.05*dstmax(ibd),step0)

          ! Near caustics (within r2 and out of r1)?
          ! Then use smaller drho to avoid large error
          do js=3,nsgl   ! loop through caustics
            if(hsgl(js).ge.r1 .and. hsgl(js).lt.r2) 
     &        drho=max(drho/4.,5.)
          enddo

          ! Get idt,ida: datum number for this path
          id=0
          idt=0
          ida=0
          do ib=1,nobst
            id=id+1                 ! datum index (=ib for time data)
            jb=nbt(ib)              ! time datum ib is in band jb
            if(jb.eq.ibd) idt=id
          enddo
          do ib=1,nobsa
            id=id+1                 ! datum index (=nobst+ib for amp data)
            jb=nba(ib)              ! amplitude datum ib is in band jb
            if(jb.eq.ibd) ida=id
          enddo
  
50        continue   
          rho=rho+drho
          rhoc=rho-0.5*drho       ! rho at center of dS
          hlfrhoc2=0.5*rhoc*rhoc
          rh=hlfrhoc2*hdiff
          tlargest=0.
  
          nrok=0          ! checks if tdetour still OK for this rho
  
          ! step size in theta 
          dtheta=drho/rhoc
          ntheta=twopi/dtheta+1
          dtheta=twopi/ntheta
          dS=rhoc*drho*dtheta     ! surface element

          theta=0.0       ! q=(cos theta, sin theta), ie q1~SV, q2~SH
55        theta=theta+dtheta
          thetac=theta-0.5*dtheta         ! theta at center of dS
          q1=rhoc*cos(thetac)
          q2=rhoc*sin(thetac)
  
          ! frequency integral interpolation for dT
          tdetour=0.5*(h11(k0)*q1*q1+h22(k0)*q2*q2)
          tlargest=max(tlargest,tdetour)
          if(tdetour.gt.dtmax(ibd)) goto 58    ! skip if outside max ellipse
          nrok=1

          ! interpolate kernels for band ibd and detour time tdetour
          call getk(dknl,nobst,nobsa,spr,tdetour,tknl,aknl,dtbl,ntbl,
     &              nbt,nba,ibd)

          dl=(rr-q1)*alpha       ! distance to last plane
          dlmin=min(dl,dlmin)
          dlmax=max(dl,dlmax)
          dV=dl*dS               ! volume element
          ! note that dl could occasionally be negative, in which case it is
          ! correcting the previous one. Next ray node will add again.

          call xyzscat(x0,y0,q1,q2,rotm,e2g,xsg,ysg,zsg,rsg,tsg,psg)
c         if(jdebug.gt.2.and.jdebug.eq.k0) 
c    &        write(13,*) 'Scatterer at',90-tsg*r2d,psg*r2d,rsg

          ! fold back at boundary if total reflection
          ! TODO: check du sign for folded scatterer
          if(rsg.gt.rn) call fold(rsg,tsg,psg,xsg,ysg,zsg,rn)
          if(rsg.lt.rcmb.and.ktype(k0).eq.2)
     &      call fold(rsg,tsg,psg,xsg,ysg,zsg,rcmb)
            
          ! TODO check early backscatter for pP, PP etc

          ! add to matrix elements
          call filla(xsg,ysg,zsg,nobst,nobsa,aa,abinv,loc,dV,dknl,
     &          maslov,q0,ktype(k0),jssw,idt,ida)
 
          ! diagnostic for first dT, dA datum
          im=mod(abs(maslov),2)+1   !im=1 for sin, 2 for cos dependence
          ima=3-im                  !change sin <--> cos for amplitudes
          if(ibd.eq.nbt(1)) then    ! diagnostic for first dT datum
            nquad=nquad+1
            dskt=dknl(im,1)         
            diskint=diskint+dS*dskt   
          endif
          if(ibd.eq.nba(1)) then    ! diagnostic for first dA datum
            dska=dknl(ima,nobst+1)  
            if(dska.gt.0) diskpstv=diskpstv+dS*dska
            if(dska.lt.0) diskngtv=diskngtv+dS*dska
            if(jssw(3).eq.1) then     ! sensitivity to dQ
              dskq=dknl(im,nobst+1)/2.*q0
              diskq=diskq+dS*dskq
            endif
          endif

58        if(theta.lt.twopi-0.0001) goto 55      ! end of theta loop
60        continue

          ! diagnostic for first dT,dA datum
          if( jdebug.gt.2.and.jdebug.eq.k0.and.
     &        (ibd.eq.nbt(1).or.ibd.eq.nba(1)) ) 
     &      write(13,1309) rhoc,drho,dtheta,ntheta,tlargest,diskint, 
     &      diskpstv,diskngtv,diskq 
1309      format(f10.0,f10.1,f10.3,i10,f10.1,4e12.3)

          if(rho.lt.dstmax(ibd).and.nrok.gt.0) goto 50     ! end of rho loop

          ! Check disk quadrature errors for first dT,dA datum
          if(ibd.eq.nbt(1)) then    ! for dT
            diskerrt=(-diskint*c0-1.0)*100.  
!            if(abs(diskerrt).gt.9.9) print *,'WARNING: KT quadrature ',
!     &        'imprecision for band',ibd,' ray node',k0,
!     &        ' diskerrt=',nint(diskerrt),'%'
          endif
          if(ibd.eq.nba(1)) then    ! for A kernel
            diskerra=(1-abs(diskpstv/diskngtv))*100.  
            if(abs(diskerra).gt.18.0) print *,'WARNING: KA quadrature ',
     &        'imprecision for band',ibd,' ray node',k0,
     &        ' diskerra=',nint(diskerra),'%'
          endif
          if(jssw(3).eq.1 .and. ibd.eq.nba(1)) then    ! for Q kernel
            diskerrq=(-diskq/wbar(ibd)/q0*2.*c0-1.)*100.
            if(abs(diskerrq).gt.20.0) print *,'WARNING: KQ quadrature ',
     &        'imprecision for band',ibd,' ray node',k0,
     &        ' diskerrq=',nint(diskerrq),'%'
          endif

          if(jdebug.ge.2.and.ibd.eq.nbt(1)) write(13,1310) k0,diskint, ! for dT
c         if(jdebug.ge.2.and.ibd.eq.nba(1)) write(13,1310) k0,diskint, ! for dA
     &      -1.0/c0,diskerrt,diskpstv,diskngtv,diskerra,
     &      diskq,-wbar(ibd)*q0/2./c0,diskerrq,
     &      nquad,alpha*r2d,rhoc,dlmin,h,dlmax
1310      format(i4,2e13.4,f10.3,2e13.4,f10.3,2e13.4,f10.3,i8,7f10.2)

65      enddo  ! end of ibd=1,nband loop

80      continue
        nquadtot=nquadtot+nquad
c       if(jdebug.gt.1) jdebug=1         ! limit debug output
85    enddo     ! end of ray node loop


!-------------------------------------------------------------

!     remove very small elements from A and write out

!-------------------------------------------------------------


      ! Note on the column order in matrix Aij
      ! Vp: columns 1-np, Vs: columns np+1 to 2*np,
      ! 1/Qs: columns 2*np+1 to 3*np

      alimfraction=0.01
      if(jdebug.gt.0) write(13,*) 'Sparse matrix winnowing for lim ',
     &   alimfraction
c1319  format(/,'Sparse matrix winnowing for lim',f8.5,//,'  Date',
c     &  '  Station      nb  kount   non0      time         sum',
c     &  ' sumsparse sparsity(%) reduct(%)  accur(%)')
      
      do nb=1,ntot  ! datum
      
        ! do statistics and sparse matrix format. Note that we number columns
        ! for each parameter (Vp,Vs,Qs) even if jssw=0. E.g.,
        ! the columns for Vs always start at column np+1.
        ! c.f. line:  kountpar(jtype)=kountp
        kount=0 ! nr of sparse elements of the complete row (Vp+Vs+Qs)

        do jtype=1,3
          kountp=0  ! individual count of sparse elements for one parameter (Vp or Vs or Qs)
          if(jssw(jtype).eq.1) then
            absum=0.
            sum0=0.
            sum0pstv=0.
            sum0ngtv=0.
            sumq0=0.
            amx=0.
            non0=0  !nr of nonzero elements for one parameter (Vp or Vs or Qs)
            j0=(jtype-1)*np
            do j=1,np
              aa1=aa(j0+j,nb)
              aabs=abs(aa1)
              absum=absum+aabs
              if(aabs.gt.0.) non0=non0+1
              amx=max(aabs,amx)
              if(nb.le.nobst.and.jtype.ne.3) then ! KT 
                sum0=sum0+aa1
              elseif(nb.gt.nobst.and.jtype.ne.3) then   ! KA 
                if(aa1.gt.0.) then ! KA positive 
                  sum0pstv=sum0pstv+aa1
                else                   ! KA negative 
                  sum0ngtv=sum0ngtv+aa1 
                endif
              elseif(nb.gt.nobst.and.jtype.eq.3) then  !KQ
                sumq0=sumq0+aa1
              endif
            enddo
            average=absum/non0
            alim=alimfraction*average
    
            sumsparse=0.
            sumpstv=0.
            sumngtv=0.
            sumq=0.
            do j=1,np
              aa1=aa(j0+j,nb)
              if(abs(aa1).gt.alim) then
                kount=kount+1
                kountp=kountp+1
                ja(kount)=j0+j
                asparse(kount)=aa1
                if(nb.le.nobst.and.jtype.ne.3) then ! KT 
                  sumsparse=sumsparse+aa1  
                elseif(nb.gt.nobst.and.jtype.ne.3) then   !KA 
                  if(aa1.gt.0) then ! KA positive 
                    sumpstv=sumpstv+aa1
                  else                  ! KA negative 
                    sumngtv=sumngtv+aa1
                  endif
                elseif(nb.gt.nobst.and.jtype.eq.3) then        !KQ
                  sumq=sumq+aa1  
                endif
              endif
            enddo
            
            ! write raymatrix.out.<ident> and raymatrix.info.<ident>
            sparsity=(100.*kountp)/np
            reduction=(100.*kountp)/non0
            if(nb.le.nobst.and.jtype.ne.3) then ! KT
              acc=100.*abs(sumsparse/trtime)
              acc0=100.*abs(sum0/trtime)
              iband=nbt(nb)
              xco  =corcoeft(nb)
              nrowt=nrowt+1 
              nrow = nrow+1   ! total number of matrix rows written (nT+nA)
              ! line in raymatrix.out.<ident>
              write(1,100) idate,stationcode(1:8),iband,cpar(jtype),
     &          'T',kountp,non0,trtime,sum0,sumsparse,sparsity,
     &          reduction,acc,acc0
!              if(abs(acc-100.).gt.5.0) print *,'WARNING: KT quadrature
!     &          ',' imprecision for band',iband,' acc=',nint(acc),'%'
            elseif(nb.gt.nobst.and.jtype.ne.3) then    ! KA 
              iband=nba(nb-nobst)
              xco  =corcoefa(nb-nobst)
              nrowa = nrowa+1
              nrow = nrow+1   ! total number of matrix rows written (nT+nA)
              ! line in raymatrix.out.<ident>
              acc =100.*abs(sumpstv/sumngtv)
              acc0=100.*abs(sum0pstv/sum0ngtv)
              write(1,100) idate,stationcode(1:8),iband,
     &          cpar(jtype),'A',kountp,non0,sumpstv,sumngtv,0.,
     &          sparsity,reduction,acc,acc0
              if(abs(acc-100.).gt.10.) print *,'WARNING: KA quadrature
     &          ',' imprecision for band',iband,' acc=',nint(acc),'%'
            elseif(nb.gt.nobst.and.jtype.eq.3) then      !KQ
              iband=nba(nb-nobst)
              xco  =corcoefa(nb-nobst)
              ! line in raymatrix.out.<ident>
              sumtheo=wbar(iband)*tstarq/2.
              acc=100.*abs(sumq/sumtheo)
              acc0=100.*abs(sumq0/sumtheo)
              write(1,100) idate,stationcode(1:8),iband,
     &          cpar(jtype),'A',kountp,non0,sumtheo*100.,sumq0*100.,
     &          sumq*100.,sparsity,reduction,acc,acc0
              if(abs(acc-100.).gt.10.) print *,'WARNING: KQ quadrature
     &          ',' imprecision for band',iband,' acc=',nint(acc),'%'
            endif   
          endif ! if(jssw(jtype).eq.1)
          kountpar(jtype)=kountp
        enddo  ! endof jtype loop

100     format(i8,1x,a8,i2,1x,a,1x,a1,i6,i7,7f10.1)
104     format(i7,1x,i6,1x,i2,1x,a1,1x,i2,1x,a8,1x,a4,1x,
     & 2f9.3,1x,f6.3,1x,2f8.3,1x,f6.3,f10.1,1x,i1,1x,a4,
     & i10,i10,2f9.3,1x,f6.1)
        
        if(nb.gt.nobst) then    ! data are dA/A
          iout=12
          jout=18
          dtyp='A'
          iband=nba(nb-nobst)
          dobs=aobs(nb-nobst)
          dtheor=1.
          s=asig(nb-nobst)
        else                    ! data are travel times
          iout=11
          jout=17
          dtyp='T'
          iband=nbt(nb)
          ! TODO handle linked data correctly here
          dobs=tobs(nb)
          dtheor=trtime       
          s=tsig(nb)
        endif

        ! write one matrix row, rhs element & corrections
        ! the three 0's at the end are placeholders for a linked ray
        if(jdebug.eq.0) then                    ! unformatted  
          write(iout) ndata,kount,iband,dobs,dtheor,s,ecorr,corcrust,
     &      corelev,tstar,qray,ievt,kluster,stationcode,netw,comp,
     &      slowness,raz,rdel,angle0,vsrc,krtyp,0.,0.,0
          write(iout) kountpar
          write(iout) (ja(i),i=1,kount)
          write(iout) (asparse(i),i=1,kount)
        else  
          write(iout,110)ndata,kount,iband,dobs,dtheor,s,ecorr,corcrust,
     &      corelev,tstar,qray,ievt,kluster,stationcode,netw,comp,
     &      slowness,raz,rdel,angle0,vsrc,krtyp,0.,0.,0
110       format(2i8,i4,2f10.2,5f7.2,f10.1,i8,i3,1x,a16,1x,a8,1x,a3,
     &           f10.1,3f10.5,f8.3,i2,2f10.2,i8)
          write(iout,*) kountpar
          write(iout,*) (ja(i),i=1,kount)
          write(iout,*) (asparse(i),i=1,kount)
        endif  

c-------New location: write raymatrix.info files
        ! line in raymatrix.info.<dtyp>.<ident>
        rmsbiband=0.
        if(iband.gt.0) rmsbiband=rmsb(iband)
        write(jout,104) ndata,ievt,kluster,dtyp,iband,
     &  stationcode(1:8),netw,rlat,rlon,relev,
     &  del,az,xco,rmsbiband,kunit,comp,      
     &  idate,iotime, slat,slon,sdep
c-----  END New location: write raymatrix.info files
        
      enddo  ! end of nb(datum) loop

      goto 10   ! read data for next ray

c     END of reading/writing data rows
900   continue

c      close(17)
      
      
! debug plot output of the 1st source-receiver path
      if(jdebug.gt.0 .and. ndata.eq.1 ) then   
        open(15,file='matx')    ! matx has same format as 3D models
        write(15,*) 3
        write(15,*) np
        do j=1,np
          x=points(1,j)
          y=points(2,j)
          z=points(3,j)
          write(15,*) x,y,z,aa(j,1)  ! only for the 1st datum
        enddo
        close(15)
      endif

!     output file: stats of epicentral distances
      open(4,file='histo.xy.'//dataf)
      do i=360,1,-1
        j=i
        if(histodistance(i).gt.0) goto 910
      enddo
910   do i=1,j
        write(4,*) i,histodistance(i)
      enddo
      
      write(1,fmt='(//,"Histogram of epicentral distances")')
      do i=1,360,10
        ksum=0
        do j=0,9
          ksum=ksum+histodistance(i+j)
        enddo
        write(1,fmt='(i3,"-",i3, " deg:",i8)') i,i+9,ksum
      enddo
      write(1,*) '>360 deg:',histodistance(361)
      itm2=time()
      write(1,fmt='(//,a)') ctime(itm2)
      write(1,*) 'Computation time:',itm2-itm1,' sec'
      print *,'End of raymatrix, computation time:',itm2-itm1,' sec'
      print *,'Number of src/rcv pairs read (ndata)   :',ndata
      print *,'Number of   matrix rows written (nrow) :',nrow
      print *,'Number of T matrix rows written (nrowt):',nrowt
      print *,'Number of A matrix rows written (nrowa):',nrowa
      
      end

!----------------------------------------------------------------

!  SUBROUTINES

!----------------------------------------------------------------


      subroutine cart(theta,phi,x,y,z)
      ! used in euler_rot
      ! colatatitude theta and longitude phi -> Cartesian x,y,z
      s=sin(theta)
      x=s*cos(phi)
      y=s*sin(phi)
      z=cos(theta)
      return
      end
      
!-----------------------------------------------------------

      subroutine coordl(ai,ph,rotm)

      DIMENSION rotm(3,3)

c computes the matrix rotm to find ray-centered coordinate directions
c  x'=rotm * x; where x' is the coordinates in the x original frame
c  rotm= | p1 sv1 sh1 |
c        | p2 sv2 sh2 |
c        | p3 sv3 sh3 |
c here p=vector tangential to ray, sv=vector normal to ray in equatorial
c plane, sh= p x sv (cross product).

c cf. Dahlen et al. (GJI, 2000): q1 is in direction of sv, q2 of sh.

c input: ray incidence angle ai, epicentral distance ph (radians)
c output: rotm

      cp=cos(ph)
      sp=sin(ph)
      ci=cos(ai)
      si=sin(ai)
      p1=cp*ci-sp*si
      p2=sp*ci+cp*si
      rotm(1,1)=p1
      rotm(2,1)=p2
      rotm(3,1)=0.
      rotm(1,2)=p2
      rotm(2,2)=-p1
      rotm(3,2)=0.
      rotm(1,3)=0.
      rotm(2,3)=0.
      rotm(3,3)=-p1*p1-p2*p2


      return
      end

c Bug fix from Guust, 2009/04/25:
c I just discovered a bug in subroutine coordl that can be serious, 
c depending on the behaviour of the Fortran compiler. The last lines:
c 
c 1193       rotm(3,1)=0.
c 1194       rotm(3,2)=0.
c 1195       rotm(3,3)=-p1*p1-p2*p2
c
c make no sense because they fill the last *row* of rotm, rather than 
c the column. It overwrites zeroes in (3,1) and (3,2) that were already 
c there; (3,3) is written in the correct place, but the first two elements 
c of the third column are never filled in. I suspect the g77 compiler set 
c them to 0 anyway, but that is not official Fortran. The code should be:

c 1193       rotm(1,3)=0.
c 1194       rotm(2,3)=0.
c 1195       rotm(3,3)=-p1*p1-p2*p2
c 
c The routine is used in raymatrix.f but I do not think anywhere else.
c

!--------------------------------------------------


      subroutine cross(sx,sy,sz,rx,ry,rz,px,py,pz)

      ! cross product of (sx,sy,sz) with (rx,ry,rz) gives (px,py,pz)
     
      px = sy*rz - sz*ry
      py = sz*rx - sx*rz
      pz = sx*ry - sy*rx
      return
      end

!--------------------------------------------------

        subroutine dintw(tknl,aknl,dtbl,ntbl,dotm,bndomega,nfreq,nband,
     &                   xclm)

        ! computes integrals with sin(w*dt) and cos(w*dt) terms 
        ! for fast interpolation - see Dahlen, eq (78), eg:
        ! tknl = -(1/2*pi) int [w^3 dotm^2 sin Phi]/ int [w^2 dotm^2]
        ! and
        ! aknl = -(1/2*pi) int [w^2 dotm^2 cos Phi]/ int [dotm^2]

        ! input: dotm,bndomega,nfreq,nband and xclm
        !        dotm(i,j) is amplitude spectrum; i ranges from 1 to
        !        nfreq(j), j from 1 to nband. bndomega(i,j) is the
        !        circle frequency (rad/s) for dotm(i,j).
        !        xclm(j) is the max. xcorr length for band j
        ! output: tknl,aknl,dtbl and ntbl
        !        tknl(i,j,k) is the value of the quotient of the
        !        frequency integrals in Dahlen eq (78) for detour time
        !        (i-1)*dtbl(k) and frequency band k. 
        !        j=1 for the sin(w*dT)
        !        integral (Maslow=0), j=2 for the cos one (Maslow -1).
        !        aknl, as tknl but for amplitudes, j=1 for sin(w*dT)
        !        (Maslow -1), j=2 for cos (Maslow 0).
        !        Both j=1,2 are for positive sign in front of the integral.
        !        For band k, i ranges from 1 to ntbl(k) for both T and A.

        ! tested on analytical solution, 4/25/06 by GN using tstdintw.f

        parameter(NTM=1500, NFR=20)      ! keep same as in main

        dimension tknl(NTM,2,NFR),aknl(NTM,2,NFR),nfreq(NFR)
        dimension bndomega(200,NFR),dotm(200,NFR),dtbl(NFR),ntbl(NFR),
     &            xclm(NFR)
        dimension gt1(NTM),gt2(NTM),ga1(NTM),ga2(NTM)
        dimension gt1old(NTM),gt2old(NTM),ga1old(NTM),ga2old(NTM)

        data twopi/6.2831853/
c       data dtmax/50./         ! max detour time 50 seconds

        jdebug=0

        if(nband.gt.NFR/2) stop 'dintw.f: too many frequency bands'

        do kb=1,nband           ! we have data for nband freq bands
          
c         if(xclm(kb).eq.0.) goto 105  ! go to next band
          if(jdebug.gt.0) write(13,*) 'dintw band:',kb

          do j=1,2
            do i=1,NTM
              tknl(i,j,kb)=0.   ! kernel numerators for time i, band kb
              aknl(i,j,kb)=0.   ! j=1 for sin(w*dT), j=2 for cos(w*dT)
            end do
          end do  
          suma=0.               ! amplitude kernel denominator (mdot^2)
          sumt=0.               ! time kernel denominator (t^2*mdot^2)

          nwb=nfreq(kb)
          if(nwb.gt.200) stop 'dintw.f: too many frequencies in dotm'
          wmax=bndomega(nwb,kb)
          dwmin=0.01*twopi/xclm(kb)  ! min step size in omega*delay is pi/xclm(kb)
          print *,'dwmin=',dwmin,' for wmax=',wmax

          dtbl(kb)=0.1*twopi/wmax      ! step for dT interpolation
          ntbl(kb)=xclm(kb)/dtbl(kb)+1
          if(ntbl(kb).gt.NTM) then
            print *,'kb,xclm,wmax=',kb,xclm(kb),wmax
            print *,'dtbl,ntbl=',dtbl(kb),ntbl(kb)
            stop 'dintw: xclm too large (nt>NTM)'
          endif  
          w=bndomega(1,kb)
          w2=w*w
          dm2=dotm(1,kb)**2
          ft=w2*dm2             ! denominator for time kernels
          fa=dm2                ! and for amplitude kernels
          do i=1,ntbl(kb)
            tdelay=(i-1)*dtbl(kb)
            gt1(i)=w*w2*dm2*sin(w*tdelay)
            gt2(i)=w*w2*dm2*cos(w*tdelay)
            ga1(i)=w2*dm2*sin(w*tdelay)
            ga2(i)=w2*dm2*cos(w*tdelay)
          enddo

          do ib=2,nwb           ! dotm is specified at nwb nodes

            wold=bndomega(ib-1,kb)
            dwbnd=bndomega(ib,kb)-wold
            nw=dwbnd/dwmin+1
            dw=dwbnd/nw
            ddm=(dotm(ib,kb)-dotm(ib-1,kb))*dw/dwbnd
            dm0=dotm(ib-1,kb)

            ! trapezoidal quadrature, factor 0.5 omitted
            do j=1,nw           ! divide into nw subintervals
              ftold=ft          ! integrands at left side of trapezium
              faold=fa
              do i=1,ntbl(kb)
                gt1old(i)=gt1(i)
                gt2old(i)=gt2(i)
                ga1old(i)=ga1(i)
                ga2old(i)=ga2(i)
              enddo
              w=wold+j*dw
              w2=w*w
              dm=dm0+j*ddm
              dm2=dm*dm
              ft=w2*dm2         ! integrands on right side
              fa=dm2
              suma=suma+dw*(faold+fa)
              sumt=sumt+dw*(ftold+ft)
              do i=1,ntbl(kb)
                tdelay=(i-1)*dtbl(kb)
                gt1(i)=w*w2*dm2*sin(w*tdelay)
                gt2(i)=w*w2*dm2*cos(w*tdelay)
                ga1(i)=w2*dm2*sin(w*tdelay)
                ga2(i)=w2*dm2*cos(w*tdelay)
                tknl(i,1,kb)=tknl(i,1,kb)+dw*(gt1(i)+gt1old(i))
                tknl(i,2,kb)=tknl(i,2,kb)+dw*(gt2(i)+gt2old(i))
                aknl(i,1,kb)=aknl(i,1,kb)+dw*(ga1(i)+ga1old(i))
                aknl(i,2,kb)=aknl(i,2,kb)+dw*(ga2(i)+ga2old(i))
              enddo     ! end loop over detour times dT
            end do      ! end w loop over subintervals
          end do        ! end loop over nodes for dotm

          ! all integrals are known now for band kb. Divide numerator
          ! and denominator (missing factor of 0.5 divides out at this 
          ! point)
          do j=1,2
            do i=1,ntbl(kb)
              tknl(i,j,kb)=-tknl(i,j,kb)/(twopi*sumt)
              aknl(i,j,kb)=-aknl(i,j,kb)/(twopi*suma)
              if(jdebug.gt.0) write(13,*) i,j,(i-1)*dtbl(kb),
     &              tknl(i,j,kb),aknl(i,j,kb)
            end do
          end do  
  
105     end do          ! end loop for frequency bands

        if(jdebug.eq.1) then
          open(31,file='T.xy')
          do k=1,nband
            do i=1,ntbl(k)
              write(31,*) (i-1)*dtbl(k),-twopi*tknl(i,1,k)
            enddo
            write(31,fmt='(">")')
          enddo
        endif  

        return
        end


!--------------------------------------------------

      subroutine euler_rot(ts,ps,tr,pr,g2e,e2g)

! calculates the Euler matrices e2g and its inverse g2e to go from
! equatorial coordinates cartesian along source-receiver path to 
! geocentric (again cartesian) and back

! input: spherical coordinates (ts,ps): theta and phi for vector s and
!        (tr,pr) for vector r, both in radians
! output: matrices g2e and e2g

! (ts,ps) and (tr,pr) define a great circle n the sphere that acts
! as the equator in a new coordinate system, where (ts,ps) is the
! origin and the longitude increases in the direction of (tr,pr)

! Usage of g2e: compute (x,y,z) geocentric, g2e*(x,y,z) gives
! cartesian coordinates equatorial.  Eg, (ts,pr) -> (0,0) and
! (tr,ps) -> (0,Delta) when expressed back in angles.

! e2g maps new (cartesian) coordinate back to the old system.

! validated with tsteuler.f

      dimension g2e(3,3),e2g(3,3)

c     jdebug=0

      call cart(ts,ps,xi,yi,zi)
      call cart(tr,pr,xr,yr,zr)

c     if(jdebug.gt.0) then
c       write(13,*) 'euler called with',ts,ps,tr,pr
c       write(13,*) 'xi,yi,zi=',xi,yi,zi
c       write(13,*) 'xr,yr,zr=',xr,yr,zr
c     endif  

      xk=yi*zr-yr*zi
      yk=xr*zi-xi*zr
      zk=xi*yr-xr*yi

      dk=sqrt(xk*xk+yk*yk+zk*zk)
      xk=xk/dk
      yk=yk/dk
      zk=zk/dk

      xj=yk*zi-yi*zk
      yj=xi*zk-xk*zi
      zj=xk*yi-xi*yk

      dj=sqrt(xj*xj+yj*yj+zj*zj)
      xj=xj/dj
      yj=yj/dj
      zj=zj/dj

c rotation matrix cos(x',x)

      g2e(1,1)=xi
      g2e(2,1)=xj
      g2e(3,1)=xk
      g2e(1,2)=yi
      g2e(2,2)=yj
      g2e(3,2)=yk
      g2e(1,3)=zi
      g2e(2,3)=zj
      g2e(3,3)=zk

      e2g(1,1)=xi
      e2g(2,1)=yi
      e2g(3,1)=zi
      e2g(1,2)=xj
      e2g(2,2)=yj
      e2g(3,2)=zj
      e2g(1,3)=xk
      e2g(2,3)=yk
      e2g(3,3)=zk

      return
      end

      !--------------------------------------------------

      subroutine filla(x,y,z,nobst,nobsa,aa,abinv,loc,dV,dknl,maslov,
     &      q0,ktype,jssw,idt,ida)

      ! fills matrix aa at the columns for location x,y,z
      ! aa(i,j) has element for node i, datum j, for one source-receiver path
      ! Quadrature is done as a Riemann sum over volume elements dV

      ! input: x,y,z=location of the quadrature element (center of dV)
      !        nobst,nobsa=nr data (freq bands for dT & dA)
      !        abinv=inverse interpolation matrix (from subroutine invrt)
      !        loc=estimate for tetrahedron index, or 0
      !        dV=volume of quadrature element h*dq1*dq2
      !        dknl(m,id)=integral kernel values at x,y,z for datum id
      !           the first nobst are for travel time kernels, followed
      !           by nobsa amplitude kernels, both sin and cos, see getk
      !        maslov= negative maslov index difference (# of -pi/2 phase changes)
      !        q0=1/Q at ray node, where Q is Qp or Qs depending ktype
      !        aa=tomography matrix, for one source-receiver path
      !        ktype=1 for P, 2 for S wave at ray node
      !        jssw=sense switches (0/1) for parameters Vp,Vs,Qs
      !        idt,ida=datum nr of dT,dA in current band; 
      !                0 if no datum in current band; -1 if all bands together
      ! output:aa, adjusted for contribution from x,y,z
      !        dskt,dska,dskq=integral kernel value for band 1 (for debug purposes)
      !        loc is index of the tetrahedron

      ! header files for Delauney triangulation
      include '../includes/nn.param'
      include '../includes/setdel3.h'       ! common for np,nt (Delauney)
      include '../includes/setdel2.h'       ! common for points,centres,vertices etc

      real*8 xd(3)              ! real*8 conforms to Sambridge's code
      
      parameter (NFR=20)        ! total (time+ampl) nr of freq bands
      dimension abinv(4,4,nt_max),bb(4)
      dimension aa(3*np_max,NFR),jssw(3)
      dimension dknl(2,NFR),qfac(3)
      data qfac/1.,1.,0.5/

c     common/debug/asum

      jdebug=0

      ! find tetrahedron index loc and node interpolation coefficients bb
      xd(1)=x         ! move to double precision for find_tet
      xd(2)=y
      xd(3)=z
      call find_tet(xd,loc)     ! loc = tetrahedron number 
      ii=vertices(1,loc)        ! first vertex in this tetrahedron
      dx=x-points(1,ii)         ! location of (x,y,z) wrt first vertex
      dy=y-points(2,ii)
      dz=z-points(3,ii)
      ! compute interpolation coefficient for each vertex
      bb(1)=abinv(1,1,loc)*dx+abinv(1,2,loc)*dy+abinv(1,3,loc)*dz+
     &      abinv(1,4,loc)
      bb(2)=abinv(2,1,loc)*dx+abinv(2,2,loc)*dy+abinv(2,3,loc)*dz+
     &      abinv(2,4,loc)
      bb(3)=abinv(3,1,loc)*dx+abinv(3,2,loc)*dy+abinv(3,3,loc)*dz+
     &      abinv(3,4,loc)
      bb(4)=abinv(4,1,loc)*dx+abinv(4,2,loc)*dy+abinv(4,3,loc)*dz+
     &      abinv(4,4,loc)

c     bsum=0.           ! debug: check on b(kk) summing to 1
c     do k=1,4  
c       bsum=bsum+bb(k)
c     enddo  
c     if(abs(bsum-1.0).gt.0.01) then
c       print *,bb
c       print *,'bsum=',bsum
c     endif  

      ! see if kernel has +/- sin/cos (omega*dT)
      im=mod(abs(maslov),2)+1   !im=1 for sin, 2 for cos dependence
      ima=3-im                  !change sin <--> cos for amplitudes
      sdV=dV                    ! apply sign to volume element
      sdVa=dV
      if(maslov.lt.0) sdVa=-dV  
      if(maslov.eq.-2) sdV=-dV

      dskt=dknl(im,1)         ! compute diagnostic for first dT datum 
      dska=dknl(ima,nobst+1)  ! and diagnostic for first dA datum 
      dskq=dknl(im,nobst+1)/2.*q0
      
      ! add to matrix elements

      do jtype=1,3

        ! Make sure data and parameter match:
        if(jssw(jtype).eq.1.and.(jtype.eq.ktype.or.jtype.eq.3)) then

          i0=(jtype-1)*np       ! starting column in matrix -1
          f=qfac(jtype)         ! dknl has already the factor -1/2pi
          if(jtype.eq.3) then   ! if parameter is Q
            f=f*q0
            ima=im
            sdVa=sdV
          endif  

          do k=1,4                  ! loop over the four vertices
            ii=vertices(k,loc)      ! get node numbers for vertices

            ! First handle delay time data
            if(jtype.lt.3) then     ! skip if parameter is Q
              if(idt.eq.-1) then    ! if ray theory, all bands the same
                do id=1,nobst
                  aa(i0+ii,id)=aa(i0+ii,id)+sdV*bb(k)*dknl(im,id)
                enddo
c               asum=asum+sdV*bb(k)*dknl(im,1)  !  for debug only
              else if(idt.gt.0) then            ! if we have delay data for this band
                aa(i0+ii,idt)=aa(i0+ii,idt)+sdV*bb(k)*dknl(im,idt)
c               if(idt.eq.1) asum=asum+sdV*bb(k)*dknl(im,1)  !  for debug only
              endif
            endif  

            if(jdebug.gt.0) write(13,1301) ii,x-points(1,ii),
     &        y-points(2,ii),z-points(3,ii),bb(k),loc
1301        format(i8,3f10.1,f10.4,i8)

            ! Now handle amplitude data, if any
            if(ida.eq.-1) then          ! ray theory (at source,receiver,caustics)
              do id=nobst+1,nobst+nobsa
                aa(i0+ii,id)=aa(i0+ii,id)+f*sdVa*bb(k)*dknl(ima,id)
                if(jdebug.gt.0) write(13,*) 'ii,aa=',i0+ii,aa(i0+ii,id)
              enddo
            else if(ida.gt.0) then      ! if we have amp data for this band
              aa(i0+ii,ida)=aa(i0+ii,ida)+f*sdVa*bb(k)*dknl(ima,ida)
              if(jdebug.gt.0) write(13,*) 'ii,aa=',i0+ii,aa(i0+ii,ida)
            endif
          enddo  ! do loop over the four vertices

        endif 
      enddo  ! do jtype=1,3
      
      return
      end

! --------------------------------------------------------------
        subroutine findtlim(tknl,nband,dtbl,ntbl,detourlim)

        ! find out for which detour time the kernels become negligible

        parameter (NFR=20,NTM=1500)       ! must be conform dintw.f
        dimension tknl(NTM,2,NFR),nfreq(NFR)
        dimension dtbl(NFR),ntbl(NFR),detourlim(NFR)

        ! input: tknl,nband,dtbl,ntbl
        !        tknl(i,j,k) is the value of the quotient of the
        !        frequency integrals in Dahlen eq (78) for detour time
        !        (i-1)*dtbl(k) and frequency band k, k=1,nband.
        !        We only use j=1 for the sin(w*dT) integral (Maslow=0)
        !        For band k, i ranges from 1 to ntbl(k) for both T and A.
        ! output: largest relevant detour time detourlim(k),k=1,nband
        !        ntbl(k) is adjusted to this time.

        do k=1,nband
          tkmax=0.
          do i=1,ntbl(k)
            tkmax=max(abs(tknl(i,1,k)),tkmax)
          enddo
          ii=ntbl(k)
          tklim=0.001*tkmax
          do while (abs(tknl(ii,1,k)).lt.tklim)
            ii=ii-1
            if(ii.lt.2) stop 'Error in findtlim'        ! debug only
          enddo  
          detourlim(k)=(ii-1)*dtbl(k)
          ntbl(k)=ii
        enddo

        return
        end

! --------------------------------------------------------------
        subroutine findmas(h11,h22,maslov)

        ! Dahlen et al (2004) eq. 105:
        ! maslov = Mxs+Mxr-Mrs = sig(Mxs+Mxr)/2 - 1
        ! where sig()= # positive - negative Hii
        ! thus maslov=0 if both h11 and h22 >0,
        ! =-1 if only one < 0, =-2 if both <0.

        ! For travel times:
        ! K ~ sin[omega*dT - dM*pi/2]

        maslov=0
        if(h11.gt.0.0.and.h22.gt.0.) return
        maslov=-1
        if(h11*h22.lt.0.) return
        maslov=-2
        return

        end

! --------------------------------------------------------------

        subroutine fold(rsg,tsg,psg,xsg,ysg,zsg,rfold)

        ! folds back at surface if rsg>rn, or CMB if rsg<rcmb 

        ! input: rsg, tsg, psg spherical coordinates of point
        !        xsg,ysg,zsg idem, cartesian
        !        rfold: folding radius (surface r or cmb)
        ! output: xsg,ysg,zsg and rsg, revised for image source
        !        at the other side of the boundary

        rsg=rsg+2.0*(rfold-rsg)
        rt=rsg*sin(tsg)
        xsg=rt*cos(psg)
        ysg=rt*sin(psg)
        zsg=rsg*cos(tsg)
        return
        end

!---------------------------------------------------------------

      function geocen(arg)
! input:  arg    = geographic colatitude (radians)
! output: geocen = geocentric colatitude (radians)
! fac=(1-f)**2
      data halfpi,fac/1.570796326794895,0.993305621334896/
      geocen=halfpi-atan(fac*cos(arg)/(max(1.0e-30,sin(arg))))
      return
      end

!---------------------------------------------------------------

      subroutine getk(dknl,nobst,nobsa,spr,tdetour,tknl,aknl,dtbl,ntbl,
     &      nbt,nba,ibd)

      ! interpolates kernel at time tdetour after multiplying tknl and
      ! aknl with location-dependent factor (in spr)

      ! input: tdetour is extra time needed to hit scatterer (in sec)
      !        spr is the factor sqrt(det[M'+M"]) in Dahlen et al (104)
      !           divided by ray node velocity (and * Omega if the
      !           scattering coefficient is not 1 - eg for conversions,
      !           and * N if the product of refl/transm coefficients
      !           is not equal to 1).
      !        nobst,nobsa= number of observations for time, ampl
      !        nbt,nba=frequency band indices
      !        tknl,aknl= time, amplitude kernels in table (see dintw)
      !        dtbl=time spacing in table (first is for time 0)
      !        ntbl=number of entries in table
      !           (see dintw for a more precise definition of these)
      !        ibd=band nr of the current call of getk
      ! output: dknl(i,j) is interpolated value of datum j
      !        sin (i=1) or cos (i=2) dependence on dT. 
      ! Amplitude data follow time data in dknl, and the second
      ! index for dknl is NOT the same as the second index in
      ! the parent arrays tknl and aknl!

      parameter (NTM=1500)       ! parameters must be conform main program
      parameter (NFR=20)        ! total (time+ampl) nr of freq bands

      dimension tknl(NTM,2,NFR),aknl(NTM,2,NFR),dtbl(NTM),ntbl(NTM)
      dimension dknl(2,NFR),nbt(NFR),nba(NFR) 

      jdebug=0

      id=0
      if(jdebug.gt.0) write(13,*)'getk for detour time',tdetour
      if(jdebug.gt.0) 
     &  write(13,*)'   k   id   jb   dtbl    tdetour    dknlt  dknla'

      ! deal with tdetour < 0
      sg=1.0
      if(tdetour.lt.0) sg=-1.0
      tdetour=abs(tdetour)

      ! kernels for delay time data
      do ib=1,nobst
        id=id+1                 ! datum index (=ib for time data)
        jb=nbt(ib)              ! time datum ib is in band jb
        if(jb.eq.ibd) then      ! in current band 
          k=tdetour/dtbl(jb)+1
          if(k.ge.ntbl(jb)) then
            dknl(1,id)=0.
            dknl(2,id)=0.
          else                    ! linear interpolation
            tdif=tdetour-(k-1)*dtbl(jb)
            dknl(1,id)=(tknl(k,1,jb)+tdif*(tknl(k+1,1,jb)-tknl(k,1,jb))/
     &      dtbl(jb))*spr
            dknl(1,id)=sg*dknl(1,id) 
            dknl(2,id)=(tknl(k,2,jb)+tdif*(tknl(k+1,2,jb)-tknl(k,2,jb))/
     &      dtbl(jb))*spr
          endif
          if(jdebug.gt.0) write(13,1301) k,id,jb,dtbl(jb),tdetour,
     &                    dknl(1,id),dknl(2,id)
1301      format(3i5,2f10.3,2e12.3)
        endif
      enddo  

      ! now the same for amplitude data. These simply follow in dknl
      ! (note that id is equal to nobst at this point)
      do ib=1,nobsa
        id=id+1                 ! datum index (=nobst+ib for amp data)
        jb=nba(ib)              ! amplitude datum ib is in band jb
        if(jb.eq.ibd) then      ! in current band
          k=tdetour/dtbl(jb)+1
          if(k.ge.ntbl(jb)) then
            dknl(1,id)=0.
            dknl(2,id)=0.
          else
            tdif=tdetour-(k-1)*dtbl(jb)
            dknl(1,id)=(aknl(k,1,jb)+tdif*(aknl(k+1,1,jb)-aknl(k,1,jb))/
     &      dtbl(jb))*spr
            dknl(1,id)=sg*dknl(1,id) 
            dknl(2,id)=(aknl(k,2,jb)+tdif*(aknl(k+1,2,jb)-aknl(k,2,jb))/
     &      dtbl(jb))*spr
          endif
        endif
      enddo  

      return
      end
      
      !--------------------------------------------------

      subroutine inpgrd

      include '../includes/nn.param'        ! parameter statements
      include '../includes/setdel3.h'       ! common for np,nt
      include '../includes/setdel2.h'       ! common for points, centres, data etc

      idebug=0
      iwrite=0
      iextend_outside_hull=0
      icall_find_node=1

      call nn_init(idebug,iwrite,icall_find_node,
     *iextend_outside_hull,nt,np)

      print *,'Read model parameterization, nt,np=',nt,np

      return
      end

      !-----------------------------------------------------

      subroutine invrt(a,ainv,n)

      ! adaptation from NR subroutine gaussj
      ! input: a(n,n), matrix of size and dimension n
      ! output: ainv, the inverse of a

      dimension a(n,n),ainv(n,n),ipiv(n),indxr(n),indxc(n)

      jdebug=0

      do j=1,n
        do i=1,n
          ainv(i,j)=a(i,j)
        enddo
      enddo  
      do j=1,n
        ipiv(j)=0
      enddo

      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(ainv(j,k)).ge.big)then
                  big=abs(ainv(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                stop 'invrt: singular matrix'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=ainv(irow,l)
            ainv(irow,l)=ainv(icol,l)
            ainv(icol,l)=dum
14        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (ainv(icol,icol).eq.0.) pause 'singular matrix.'
        pivinv=1./ainv(icol,icol)
        ainv(icol,icol)=1.
        do 16 l=1,n
          ainv(icol,l)=ainv(icol,l)*pivinv
16      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=ainv(ll,icol)
            ainv(ll,icol)=0.
            do 18 l=1,n
              ainv(ll,l)=ainv(ll,l)-ainv(icol,l)*dum
18          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=ainv(k,indxr(l))
            ainv(k,indxr(l))=ainv(k,indxc(l))
            ainv(k,indxc(l))=dum
23        continue
        endif
24    continue

      if(jdebug.gt.0) write(13,*) ainv

      return
      end

!--------------------------------------------------

        subroutine meanw(wbar,dotm,bndomega,nfreq,nband,xclm)

        ! computes average frequency:
        ! wbar = int [w dotm^2] / int [dotm^2]

        ! input: dotm,bndomega,nfreq,nband,xclm
        !        dotm(i,j) is amplitude spectrum; i ranges from 1 to
        !        nfreq(j), j from 1 to nband. bndomega(i,j) is the
        !        circle frequency (rad/s) for dotm(i,j).
        !        xclm(j) is the max. xcorr length for band j
        ! output: wbar(j) is average frequency for band j

        parameter(NFR=20)      ! keep same as in main

        dimension nfreq(NFR)
        dimension bndomega(200,NFR),dotm(200,NFR),xclm(NFR),wbar(NFR)

        data twopi/6.2831853/

        jdebug=0

        if(nband.gt.NFR/2) stop 'meanw.f: too many frequency bands'

        do kb=1,nband           ! we have data for nband freq bands

          if(jdebug.gt.0) write(13,*) 'meanw band:',kb
          nwb=nfreq(kb)
          if(nwb.gt.200) stop 'meanw.f: too many frequencies in dotm'

          ! initial value
          wbar(kb)=0.   ! numerator (int [w dotm^2])
          suma=0.    ! denominator (int [dotm^2])
          wmax=bndomega(nwb,kb)
          dwmin=0.01*twopi/xclm(kb)  ! min step size in omega*delay is pi/xclm(kb)
          w=bndomega(1,kb)
          dm2=dotm(1,kb)**2
          fa=dm2                ! integrand of denominator
          ga=w*dm2              ! integrand of numerator

          do ib=2,nwb           ! dotm is specified at nwb nodes
            ! step size of freq. integral
            wold=bndomega(ib-1,kb)
            dwbnd=bndomega(ib,kb)-wold
            nw=dwbnd/dwmin+1
            dw=dwbnd/nw
            ddm=(dotm(ib,kb)-dotm(ib-1,kb))*dw/dwbnd
            dm0=dotm(ib-1,kb)

            ! trapezoidal quadrature, factor 0.5 omitted
            do j=1,nw           ! divide into nw subintervals
              faold=fa          ! integrands at left side of trapezium
              gaold=ga
              w=wold+j*dw
              w2=w*w 
              dm=dm0+j*ddm
              dm2=dm*dm
              fa=dm2         ! integrands on right side
              suma=suma+dw*(faold+fa)
              ga=w*dm2
              wbar(kb)=wbar(kb)+dw*(ga+gaold)
            end do      ! end w loop over subintervals
          end do        ! end loop over nodes for dotm

          ! all integrals are known now for band kb. Divide numerator
          ! and denominator (missing factor of 0.5 divides out at this point)
          if(jdebug.gt.0) write(13,*) 'int[w dotm^2], int[dotm^2]=',wbar(kb),suma
          wbar(kb)=wbar(kb)/suma
          if(jdebug.gt.0) write(13,*) 'wbar=',wbar(kb)

        end do          ! end loop for frequency bands

        return
        end

!---------------------------------------------------------------

      subroutine model1(fname)

      ! as model.f but without spline initialization

      ! Reads a model from file <fname>
      ! Model units are assumed to be MKS, unless the first density
      ! read is less than 100, in which case km and km/s units are
      ! assumed
      ! This is the same file format as used for love.f and rayleigh.f
      ! (with moho node added in line 3)

      ! Format of the model file [unit 1, opened/closed by the routine]:
      ! Line 1: ignored
      ! Line 2: ignored
      ! Line 3: n (# of nodes), nic (top node of inner core), noc (top
      !         node of outer core), nmoh (top node of moho)
      ! Lines 4ff: radius, density, Vp, Vs

      ! The nodes are ordered with INCREASING radius
      ! Discontinuities are indicated by two subsequent nodes
      ! with the same radius. Nodes closer than 0.01 km are collapsed
      ! into a discontinuity.

      ! input: model file name fname (character*72)
      ! output: common block /modl/ is filled with model information
      !         including spline coefficients.

      character*72 fname
      
      ! background model
      common /modl/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &      Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh
      
      jdebug=0

      nmoh=0            ! safeguard against old model files
      open(10,file=fname)
      read(10,*) 
      read(10,*) 
      read(10,*) n,nic,noc,nmoh
      if(nmoh.le.noc.or.nmoh.gt.n) stop 'model: Error in Moho index'
      if(n.gt.500) stop 'Cannot handle models with n>500'
      mks=1
      read(10,*) r(1),rho,vp(1),vs(1)
      if(rho.lt.100.) mks=0
      backspace(10)
      do i=1,n
        read(10,*) r(i),rho,vp(i),vs(i),Qs(i)
        if(mks.eq.1) then
          r(i)=0.001*r(i)
          vp(i)=0.001*vp(i)
          vs(i)=0.001*vs(i)
        endif
        if(jdebug.gt.0) write(13,*) i,r(i),vp(i),vs(i),Qs(i)
      end do  
      close(10)
      
   60 nsl=n   
      if(vs(nsl).gt.0.) go to 70
   65 nsl=nsl-1
      if(vs(nsl).le.0.) go to 65
   70 nicp1=nic+1
      nocp1=noc+1
      nslp1=nsl+1

      if(nsl.ne.n) stop 'Oceanic reference model not allowed'

      rn=r(n)
      do 45 i=2,n
      if(abs(r(i)-r(i-1)).lt.0.01) r(i)=r(i-1)
   45 continue

      return
      end        

!-----------------------------------------------------

      subroutine rdray(kin,kout,chray,phase,rseg,ktseg,kdwn,nseg)

      ! read ray definition

      ! input:
      ! kin - I/O unit number for input (5 for screen) [assumed open]
      ! kout - I/O unit number for output (0 for no outp) [assumed open]
      ! output:
      ! chray - ray identification (e.g. PKiKP) from first line of kin
      ! phase - ray ident for ellipticity corrections (max 8 char)
      ! rseg(i) - start radius of segment i
      ! ktseg(i) - wave type of segment i (1=P,2=S,last point of ray 
      !            may be 0, but will be set to 1 or 2 on output)
      ! kdwn(i) - ray direction (1=down,0=up,2=turning point,3=refl,
      !           4=transmission, 5=last point of ray)
      ! nseg - number of segments

      ! the ray is defined by segments and their starting radii
      ! the first (source) one specifies: radius, up/down, P/S
      ! for each subsequent segment: radius,refl/trans/turn,P/S
      ! and for the last point of the ray: (radius,5,0))
      ! note that a turning point separates two different segments
      ! (the radius is then interpreted as the MINIMUM radius for 
      ! the turning point).

      ! Examples (segments only):
      ! P wave, source at 10 km depth:
      ! 6361 1 1
      ! 3480 2 1
      ! 6371 5 0
      ! PcS wave, source 400 km depth:
      ! 5971 1 1
      ! 3480 3 2
      ! 6371 5 0
      ! pPKP wave, 100 km depth:
      ! 6271 0 1
      ! 6371 3 1
      ! 3480 4 1
      ! 1221.5 2 1
      ! 3480 4 1
      ! 6371 5 0


      dimension rseg(20),ktseg(20),kdwn(20)
      character*16 chray
      character*1 wtyp(0:2)
      character*4 ktyp(0:5)
      character*8 phase

      data wtyp/' ','P','S'/ 
      data ktyp/'up','down','turn','refl','tran','end'/
      jdebug=0

      if(kin.eq.5) then
        print *,'Give ray description, then ellipticity ident, then'
        print *,'Input of ray segments (one per line):'
        print *,'Give radius (km), kdwn, wavetype (P=1,S=2),'
        print *,'where kdwn=1 for downward start, 0 for upward start,'
        print *,'2 for max turning depth, 3 for reflection, 4 for'
        print *,'transmission, 5 for end of ray'
      endif

      i=1
      read(kin,fmt='(a)') chray
      print *,'Ray ident: ',chray
      if(kout.ne.0) write(kout,fmt='(a)') chray
      read(kin,fmt='(a)') phase
      if(kout.ne.0) write(kout,fmt='(a)') phase
      read(kin,*) rseg(i),kdwn(i),ktseg(i)
      if(kout.ne.0) 
     &      write(kout,40) rseg(i),kdwn(i),ktseg(i),ktyp(kdwn(i)),
     &      wtyp(ktseg(i))
40    format(f8.1,2i5,2x,a4,2x,a1)
      if(kdwn(1).gt.1) stop 'First ray point should have kdwn 1 or 0'
      rmin=rseg(1)
      rmax=rmin
      do while(kdwn(i).ne.5)
        i=i+1
        if(i.gt.20) stop 'Max segments is 20'
        read(kin,*) rseg(i),kdwn(i),ktseg(i)
        if(ktseg(i).eq.0) ktseg(i)=ktseg(i-1)
        if(kout.ne.0) 
     &        write(kout,40) rseg(i),kdwn(i),ktseg(i),ktyp(kdwn(i)),
     &        wtyp(ktseg(i))
        if(kdwn(i).lt.2) stop 'kdwn < 2 only allowed for first segment'
        rmin=min(rmin,rseg(i))
        rmax=max(rmax,rseg(i))
        if(jdebug.eq.1) print *,rseg(i),kdwn(i),ktseg(i),ktyp(kdwn(i)),
     &        wtyp(ktseg(i))
      end do
      if(jdebug.eq.1) print *,'nseg=',i
      nseg=i
      print *,'Rays between r=',rmin,rmax

      return
      end

!---------------------------------------------------------------

      subroutine rdband2(kin,kout,bndomega,dotm,nfreq,nband,period)

      ! reads spectral band info
      ! input: I/O unit numbers for input and output kin,kout
      !        (assumed open, no output if kout=0).
      ! output: bndomega=frequencies (rad/s), dotm=mdot(omega), nfreq=# of frequencs

      parameter (NFR=20)
      dimension bndomega(200,NFR),dotm(200,NFR),nfreq(NFR),period(NFR)

      pi=4.0*atan(1.0)
      twopi=2.0*pi

      read(kin,*) nband
      if(nband.gt.NFR) then
        print *, 'nband >', NFR, ': increase dimensions'
        stop
      endif
      if(kout.ne.0) write(kout,fmt='(i5)') nband
      do ib=1,nband
        read(kin,*) nfreq(ib)
        if(nfreq(ib).gt.200) stop 'nfreq>200, increase dimensions'
        if(kout.ne.0) write(kout,fmt='(i5)') nfreq(ib)
        dmmax=0.
        do j=1,nfreq(ib)
          read(kin,*) bndomega(j,ib),dotm(j,ib)
          if(dotm(j,ib).gt.dmmax) then
            dmmax=dotm(j,ib)
            pmax=bndomega(j,ib)
          endif
          if(kout.ne.0) 
     &          write(kout,fmt='(2f8.4)') bndomega(j,ib),dotm(j,ib)
        end do
        period(ib)=twopi/pmax   ! assume dominant period = max dotm
      end do

      return
      end

!---------------------------------------------------------------

      subroutine rtp2xyz(r,t,p,x,y,z)

      ! convert spherical coordinates to cartesian

      ! input: r,t=theta,p=phi
      ! output: cartesian coordinates x,y,z

      rst=r*sin(t)
      x=rst*cos(p)
      y=rst*sin(p)
      z=r*cos(t)

      return
      end

!-------------------------------------

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
!-------------------------------------

      subroutine xyzscat(xl,yl,q1,q2,rotm,e2g,xsg,ysg,zsg,rsg,
     &      tsg,psg)
      DIMENSION rotm(3,3), e2g(3,3)
c     data jdebug/0/

! find scatterer's position in quadrature plane

! input:
!     xl,yl: ray cartesian coord in equatorial plane (source-rec plane)
!     q1,q2: scatterer location in ray coordinates 
!     rotm: rotation matrix computed by coordl
!     e2g: equatorial-to-geocentric tranformation matrix
! output:
!     xsg,ysg,zsg: scatterer cartesian coordinates, geocentric
!     rsg,tsg,psg: scatterer radius, theta, phi, geocentric

      jdebug=0
      if(jdebug.gt.0) write(13,*) 'xyzscat for xl,yl=',xl,yl,
     &   ', q1,q2=',q1,q2

      ! find equatorial cartesian coordinates of scatterer
      xs=xl+rotm(1,2)*q1+rotm(1,3)*q2
      ys=yl+rotm(2,2)*q1+rotm(2,3)*q2
      zs=rotm(3,2)*q1+rotm(3,3)*q2
      if(jdebug.gt.0) write(13,*) 'xs,ys,zs=',xs,ys,zs
     

      ! convert to spherical
      ! call xyz2rtp(xs,ys,zs,rs,ts,ps)
      if(jdebug.gt.0) write(13,*) 'rs,ts,ps=',rs,ts,ps

      ! convert to geocentric cartesian
      xsg=e2g(1,1)*xs+e2g(1,2)*ys+e2g(1,3)*zs
      ysg=e2g(2,1)*xs+e2g(2,2)*ys+e2g(2,3)*zs
      zsg=e2g(3,1)*xs+e2g(3,2)*ys+e2g(3,3)*zs
      if(jdebug.gt.0) write(13,*) 'xsg,ysg,zsg=',xsg,ysg,zsg

      ! convert to spherical
      call xyz2rtp(xsg,ysg,zsg,rsg,tsg,psg)
      if(jdebug.gt.0) write(13,*) 'rsg,tsg,psg=',rsg,tsg,psg

      return
      end
   
c-------------------------------------------------------------

c MALCOLM SAMBRIDGE'S DELAUNEY SUBROUTINES

c-------------------------------------------------------------

      ! the following routines are from nnset.f, written by Malcolm
      ! Sambridge, but somewhat adapted for raymatrix.f (GN March 2006).

        Subroutine nn_init
     &            (idebug,iwrite,
     &             icall_find_node,
     &             iextend_outside_hull,nt,np)
        include	'../includes/nn.param'
        include '../includes/setdel2.h'

        common/nnbuildlists/nnlist(nnsum_max),
     &                      ntlist(ntsum_max),
     &                      nn_pointer(np_max2+1),
     &                      nt_pointer(np_max2+1),
     &                      node(np_max2)

	common/nnswitches/nnwrite,lud
        common/find_nodec/usefindnode
        common/find_tetc/extendtetra,
     &                   hullinfo(2,nf_max),
     &                   hullfaces(3,nf_max),
     &                   hullneighbour(3,nf_max)
        common/find_tetc2/hullnorm(4,nf_max)
        common/build_lists3dc/writeout

	real*8		hullnorm
	integer		hullinfo
	integer		hullfaces
	integer		hullneighbour
        integer         nhwork1(nh_max+1)
        integer         nhwork2(nwork2d)
        integer         nhwork3(nwork2d)
	logical		node

	logical		nnwrite

        logical         consistent
        real*8		xd(3)
        real*8		ax,ay,az,bx,by,bz,dist,a,b,c,d
        logical         out
        logical         debug,writeout,usefindnode,extendtetra
        integer		cyc(4)
        data		cyc/2,3,4,1/

        character*60 nodesfile,vertexfile

      print *,'Give vertex file name (eg vertices.all.xyz):'
      read(5,fmt='(a)') vertexfile
      open(51,file=vertexfile,status='old')
      read(51,*) nd
      read(51,*) np
      if(nd.ne.3)then
         write(*,*)' nd is not equal to 3 dimensions'
         stop 'error in nn_init'
      else if(np.gt.np_max)then
         write(*,*)' Too many points in input file'
         stop 'error in nn_init'
      end if
      do i=1,np
        read(51,*)(points(j,i),j=1,nd)
      end do
      close(51)

      print *,'Give Delauney file name (eg facets.all):'
      read(5,fmt='(a)') nodesfile

      nnwrite = .false.
      debug = .false.
      writeout = .false.
      usefindnode = .false.
      extendtetra = .false.
      if(idebug.eq.1)debug = .true.
      if(iwrite.eq.1)writeout = .true.
      if(icall_find_node.eq.1)usefindnode = .true.
      if(iextend_outside_hull.eq.1)extendtetra = .true.
      if(debug)write(*,*)' Total number of points read in:',np
      if(debug)write(*,*)' calling nn3d_setup'
      
      open(50,file=nodesfile,status='old')
      lu_in=50

      call nn3d_setup
     &     (np,nt_max,np_max,nnpn_max,nwork3d,
     &      points,lu_in,data,nt,vertices,
     &      centres,neighbour,loc,nnn,nflist,
     &      ntrilist)
      close(50)
      if(debug)write(13,*)' done nn3d_setup'
 
      if(debug)write(13,fmt='(/"  Number of points = ",i7,
     &              " Tetrahedra = ",i7/)')np,nt
      if(extendtetra)then

        do i=1,np
           nnn(i) = 0
        end do
        nf = 0
        nhc = 0
	do i=1,nt
           do j=1,4
             if(neighbour(j,i).eq.0)then
               j1 = cyc(j)
               j2 = cyc(j1)
               j3 = cyc(j2)
               nf = nf + 1
               if(nf.gt.nf_max)then
                  write(*,200)
                  stop
               end if
               node1 = vertices(j1,i) 
               node2 = vertices(j2,i) 
               node3 = vertices(j3,i) 
	       if(nnn(node1).eq.0)then
                  nhc = nhc + 1
                  nodel = nhc
                  nnn(node1) = nodel
               else
                  nodel = nnn(node1)
               end if
               hullfaces(1,nf) = nodel
	       if(nnn(node2).eq.0)then
                  nhc = nhc + 1
                  nodel = nhc
                  nnn(node2) = nodel
               else
                  nodel = nnn(node2)
               end if
               hullfaces(2,nf) = nodel
	       if(nnn(node3).eq.0)then
                  nhc = nhc + 1
                  nodel = nhc
                  nnn(node3) = nodel
               else
                  nodel = nnn(node3)
               end if
               hullfaces(3,nf) = nodel
	       neighbour(j,i) = -nf
               hullinfo(1,nf) = i
               hullinfo(2,nf) = j
               ax =  points(1,node2) - points(1,node1)
               ay =  points(2,node2) - points(2,node1)
               az =  points(3,node2) - points(3,node1)
               bx =  points(1,node3) - points(1,node1)
               by =  points(2,node3) - points(2,node1)
               bz =  points(3,node3) - points(3,node1)
               a = ay*bz-az*by
               b = az*bx-ax*bz
               c = ax*by-bx*ay
               d = sqrt(a*a + b*b + c*c)
               dist = points(1,node1)*a +
     &                points(2,node1)*b +
     &                points(3,node1)*c
               if(dist.lt.0.0)then
                 d = -d
               end if
               hullnorm(1,nf) = a/d
               hullnorm(2,nf) = b/d
               hullnorm(3,nf) = c/d
               hullnorm(4,nf) = dist/d
             end if
           end do
        end do
	nh = (4+nf)/2
	if(debug)write(13,*) 'Number of faces on the hull =',nf
	if(debug)write(13,*) 'Number of nodes on the hull =',nh
	if(debug)write(13,*) 'Counted number of nodes on the hull ='
     &           ,nhc
 	if(writeout)then
	   write(*,*)' vertices of triangles and normals on hull '
	   write(*,*)' (indices are local to hull)'
        end if
        call build_hullneighbour
     &       (nh,hullfaces,nf,nwork2d,
     &        hullneighbour,nhwork1,nhwork2,nhwork3)
	if(debug)then
           write(13,*)' checking hullneighbour array'
           call check_neighbour
     &          (hullneighbour,2,hullfaces,nf,consistent)
           write(13,*)' done hullneighbour array'
	end if
        if(consistent)then
           write(*,*)' '
           write(*,*)' Result of hullneighbour check: ',
     &               'hullneighbour matrix is consistent'
           write(*,*)' '
        end if
      end if
      if(usefindnode)then
        if(writeout)write(*,*)' building nnlist and ntlist'
        call build_lists3D
     &             (np,vertices,nt,np_max,
     &              ntsum_max,nnsum_max,node,
     &              nn_pointer,nnlist,nnsum,
     &              nt_pointer,ntlist,ntsum)

        if(writeout)write(*,*)' built nnlist and ntlist'
      end if
      if(.not.debug)return
      write(13,*)' checking neighbour array'
      call check_neighbour(neighbour,3,vertices,nt,consistent)
      write(13,*)' done neighbour array'
      if(consistent)then
         write(*,*)' '
         write(*,*)' Result of neighbour check: ',
     &             'neighbour matrix is consistent'
         write(*,*)' '
      end if
      iok = 0
      ic = 0
      write(*,*)' debug check: locating centres of every tetrahedron'
      do i=1,nt
         xd(1) = 0.d0
         xd(2) = 0.d0
         xd(3) = 0.d0

         do j=1,4
            xd(1) = xd(1) + points(1,vertices(j,i))
            xd(2) = xd(2) + points(2,vertices(j,i))
            xd(3) = xd(3) + points(3,vertices(j,i))
         end do

         xd(1) = 0.25d0*xd(1)
         xd(2) = 0.25d0*xd(2)
         xd(3) = 0.25d0*xd(3)
         call tetloc(xd,points,vertices,neighbour,loc,out)
         if(out)write(13,*)' center of tet',i,' is outside of 
     &the convex hull'
         if(i.ne.loc)then
             write(13,*)' center of',i,' in tetrahedra',loc
             write(13,*)' x:',(xd(j),j=1,3)
             iok = 1
         end if
      end do
      if(iok.eq.0)then
         write(*,*)' '
         write(*,*)' Result of tetloc test:',
     &             ' All points located successfully'
         write(*,*)' '
      end if
 100  format(i2,' : ',25(i2,1x))
 101  format(i3,' : ',20(i3,1x))
 102  format(i4,' : ',15(i4,1x))
 103  format(i6,' : ',10(i6,1x))
 200  format(1x,//' Error - maximum number of faces',
     &             ' on convex hull exceeded'//,
     &             ' increase size of parameter nf_max and recompile'/)
      return
      end

!---------------------------------------------------------

        subroutine find_tet (xd,loc)
 
        include '../includes/nn.param'
        include '../includes/setdel2.h'

        common/find_tetc/extendtetra,
     &                   hullinfo(2,nf_max),
     &                   hullfaces(3,nf_max),
     &                   hullneighbour(3,nf_max)
        common/find_tetc2/hullnorm(4,nf_max)
	common/tet_loc_face/i

        real*8 		xd(3)
	real*8		hullnorm
	real*8		a,d,dnew
	integer		hullinfo
	integer		hullfaces
	integer		hullneighbour
        logical         out
        logical         extendtetra

 
        ! qhull numbering of tets starts at 0 ,we start at 1
	if(loc.eq.0) loc = 1 

        ! in case the first try comes back with "out of hull message"
        locold = loc
        ktriedalready = 0

 5      call tetloc(xd,points,vertices,neighbour,loc,out)
	if(out)then
           if(extendtetra)then  ! allow extrapolation outside hull
             int = -neighbour(i,loc)
             a = hullnorm(1,int)*xd(1)+
     &           hullnorm(2,int)*xd(2)+
     &           hullnorm(3,int)*xd(3)
             d = hullnorm(4,int)/a 
 10          continue
             do k=1,3
                int_new = hullneighbour(k,int)
                a = hullnorm(1,int_new)*xd(1)+
     &              hullnorm(2,int_new)*xd(2)+
     &              hullnorm(3,int_new)*xd(3)
                dnew = hullnorm(4,int_new)/a 
                if(dnew.gt.0.and.dnew.lt.d)then
                   d = dnew
                   int = int_new
                   go to 10
                end if
             end do
             loc = hullinfo(1,int)
           else  ! don't accept points outside hull
             if(ktriedalready.eq.0)then ! give it a second try (culprit might be roundoff error)
               ktriedalready = 1
               xd(1)=xd(1)+0.05   ! try for slightly different point
               xd(2)=xd(2)-0.1
               xd(3)=xd(3)+0.05
               loc = locold
               goto 5 
             else     ! give up
               write(*,*)'points outside hull'
               write(*,*)'xd(1),xd(2),xd(3)',xd(1),xd(2),xd(3)
	       loc = 0 
               stop
             endif

	   end if
        end if
        return
        end
c------------------------------------------------------------------------
	Subroutine nn3d_setup
     &             (np,nt_max,np_max,nnpn_max,nwork3d,
     &              points,lu_in,data,nt,vertices,
     &              centres,neighbour,loc,nnn,nflist,
     &              ntrilist)

	real*8		points(3,*)
	real*8		centres(4,*)
        real*8		data(*)
	integer		vertices(4,*)
	integer		neighbour(4,*)
	integer		nnn(*)
	integer		nflist(2,*)
	integer		ntrilist(*)
	logical		nnwrite

	common/nnswitches/nnwrite,lud
           nt = 0
           read(lu_in,*)

  1        read(lu_in,*,end=3,err=100)
     &     vertices(1,nt+1),vertices(2,nt+1),
     &     vertices(3,nt+1),vertices(4,nt+1)

           nt = nt + 1

           if(nt.ge.nt_max) go to 200
           go to 1
  3        continue


        if(nt.ge.nt_max) go to 200
	do 5 i=1,nt
           do 5 j=1,4
	   vertices(j,i) = vertices(j,i) + 1
 5      continue
	write(*,*)'ccentres3d nt',nt
        call ccentres3d
     &             (points,vertices,nt,centres)
        write(*,*)'build_nv3d np,nt,np_max,nwork3d',np,nt,np_max,nwork3d
        call build_nv3d
     &  (np,vertices,nt,np_max,nwork3d,
     &   neighbour,nnn,nflist,ntrilist)
        loc = nt/2
	return
  100   write(*,*)
     &  'Error in nn_setup: read error in Delaunay input file'
        stop
  200   write(*,*) 
        write(*,*) 'Error in nn_setup: too many tetrahedra'
        write(*,*) '      current number of tetrahedra :',nt
        write(*,*) '      current value of nt_max      :',nt_max
        write(*,*) 'Remedy: increase size of parameter nt_max'
        write(*,*) '        in calling program (see file nn.param).'
        write(*,*) 
        stop 
	end
c------------------------------------------------------------------------
	Subroutine ccentres3d(points,vertices,nt,centres)
	real*8		points(3,*)
	real*8		centres(4,*)
	real*8		p(3,4)
        real*8          x,y,z
	integer		vertices(4,*)
	logical		nnwrite

	common/nnswitches/nnwrite,lud
	write(*,*)'ccentres nt',nt
	do 5 i= 1,nt

           do 6 j=1,3
              do 7 k=1,4
	      p(j,k) = points(j,vertices(k,i))
              ! Kasra
              !WRITE(*,*) 'i, nt, p: ', i, nt, p(j,k)
              ! End Kasra
  7           continue
  6        continue

	   call circumsphere(p,x,y,z)
	   centres(1,i) = x
	   centres(2,i) = y
	   centres(3,i) = z
	   centres(4,i) = (x-p(1,1))*(x-p(1,1)) 
     &                  + (y-p(2,1))*(y-p(2,1)) 
     &                  + (z-p(3,1))*(z-p(3,1))
 5	continue
	return
	end
c------------------------------------------------------------------------
	Subroutine circumsphere(p,x,y,z)

	real*8		p(3,4)
        real*8          A(3,3),b(3),d,x,y,z,p0,xs
	integer		indx(3)
	logical		nnwrite
	common/nnswitches/nnwrite,lud
        do 5 i=1,3
           p0 = p(i,1)
           do 10 j=1,3
              ! KASRA
              !WRITE(*,*) 'p, p0', p(i,j+1), p0
              ! END KASRA
              A(j,i) = p(i,j+1) - p0
 10        continue
 5      continue
        do 15 i=1,3
           b(i) = 0.d0
 15     continue
        do 20 i=1,3
           p0 = p(i,1)
           xs = p0*p0
           do 25 j=1,3
              p0 = p(i,j+1)
              b(j) = b(j) + p0*p0 - xs 
 25        continue
 20     continue
        do 30 i=1,3
           b(i) = b(i)/2.d0
 30     continue
        call ludcmp(A,3,3,indx,d)
        call lubksb(A,3,3,indx,b)
	x = b(1)
	y = b(2)
	z = b(3)
	return
	end
c------------------------------------------------------------------------
	Subroutine build_nv3d
     &             (np,vertices,nt,np_max,nwork3d,
     &              neighbour,nnn,nflist,ntrilist)
	integer		vertices(4,*)
	integer		neighbour(4,*)
	integer		nnn(*)
	integer		nflist(2,*)
	integer		ntrilist(*)
	logical		nnwrite
	common/nnswitches/nnwrite,lud
	if(nnwrite)write(*,*)' Building neighbour v ...'
	do 5 i = 1,4
	   do 5 j = 1,nt
	      neighbour(i,j) = 0
 5      continue
        do 6 i = 1,nwork3d
	   nflist(1,i) = 0
	   nflist(2,i) = 0
	   ntrilist(i) = 0
 6      continue
        do 7 i = 1,np
	   nnn(i) = 0
 7      continue
	do 10 it = 1,nt
	   i1 = vertices(1,it)
	   i2 = vertices(2,it)
	   i3 = vertices(3,it)
	   i4 = vertices(4,it)
	   nnn(i1) = nnn(i1) + 1
	   nnn(i2) = nnn(i2) + 1
	   nnn(i3) = nnn(i3) + 1
	   nnn(i4) = nnn(i4) + 1
 10     continue
	itemp = nnn(1)+1
	nnn(1) = 1
	do 20 j = 2,np+1
	   itemp2  = itemp 
	   itemp   = itemp + nnn(j)+1
	   nnn(j) = itemp2 + 1
 20     continue
	if(nnn(np+1).ge.nwork3d)then
           write(*,*)' '
           write(*,*)' Error: array sizes too small in subroutine '
     &               ,'build_nv3d'
           write(*,*)'        maximum number of attached tetrahedra'
           write(*,*)'        is too small: current value =',nwork3d
           write(*,*)'                   must be at least =',itemp2+1
           write(*,*)' '
           write(*,*)' Remedy: Increase size of parameter nwork3d'
           write(*,*)'         in main program and recompile'
	   stop
	end if
	do 25 it = 1,nt
           do 30 j = 1,4
	      j2 = mod(j,4)+1
	      j3 = mod(j2,4)+1
	      j4 = mod(j3,4)+1
	      i1 = vertices(j,it) 
	      i2 = vertices(j2,it) 
	      i3 = vertices(j3,it) 
	      i4 = vertices(j4,it) 
              kstart = nnn(i2)
              kstop =  nnn(i2+1) - 1
              imax = max(i3,i4)
              imin = min(i3,i4)
	      jt = 0
              do 35 k = kstart,kstop
	         if(nflist(1,k).eq.0)then
                    nflist(1,k) = imax
                    nflist(2,k) = imin
                    ntrilist(k) = it
                    go to 36
	         else if(nflist(1,k).eq.imax.and.
     &	                 nflist(2,k).eq.imin)then
                    jt = ntrilist(k)
		    call update_neighbour
     &              (neighbour,vertices,it,jt,j,i2,i3,i4)
                    go to 30
                 end if
 35           continue
 36           continue
              kstart = nnn(i3)
              kstop =  nnn(i3+1) - 1
              imax = max(i2,i4)
              imin = min(i2,i4)
	      jt = 0
              do 135 k = kstart,kstop
	         if(nflist(1,k).eq.0)then
                    go to 136
	         else if(nflist(1,k).eq.imax.and.
     &	                 nflist(2,k).eq.imin)then
                    jt = ntrilist(k)
		    call update_neighbour
     &              (neighbour,vertices,it,jt,j,i2,i3,i4)
                    go to 30
                 end if
 135           continue
 136           continue
              kstart = nnn(i4)
              kstop =  nnn(i4+1) - 1
              imax = max(i2,i3)
              imin = min(i2,i3)
	      jt = 0
              do 235 k = kstart,kstop
	         if(nflist(1,k).eq.0)then
                    go to 236
	         else if(nflist(1,k).eq.imax.and.
     &	                 nflist(2,k).eq.imin)then
                    jt = ntrilist(k)
		    call update_neighbour
     &              (neighbour,vertices,it,jt,j,i2,i3,i4)
                    go to 30
                 end if
 235           continue
 236           continue
 30        continue
 25     continue
	if(nnwrite)write(*,*)' built neighbour v'
	return
	end
c------------------------------------------------------------------------
        Subroutine update_neighbour
     &  (neighbour,vertices,it,jt,j,i2,i3,i4)

	integer		neighbour(4,*)
	integer		vertices(4,*)
	integer		v1,v2,v3,v4
        neighbour(j,it) = jt
        v1 = vertices(1,jt)
        v2 = vertices(2,jt)
        v3 = vertices(3,jt)
        v4 = vertices(4,jt)
        if(v1.ne.i2.and.v1.ne.i3.and.v1.ne.i4)then
            neighbour(1,jt) = it
        else if(v2.ne.i2.and.v2.ne.i3.and.v2.ne.i4)then
            neighbour(2,jt) = it
        else if(v3.ne.i2.and.v3.ne.i3.and.v3.ne.i4)then
            neighbour(3,jt) = it
        else if(v4.ne.i2.and.v4.ne.i3.and.v4.ne.i4)then
            neighbour(4,jt) = it
        end if
	return
	end
c------------------------------------------------------------------------
	Subroutine check_neighbour(neighbour,nd,vertices,nt,consistent)

	integer		neighbour(nd+1,*)
	integer		vertices(nd+1,*)
	logical		consistent

        consistent = .true.
        if(nd.eq.2)then

	do 10 i=1,nt
	   do 20 j=1,3
              k = neighbour(j,i)
              j1 = mod(j,3)+1
              j2 = mod(j1,3)+1
              i1 = vertices(j1,i)
              i2 = vertices(j2,i)
              if(k.le.0)go to 20
	      do 30 m=1,3
                 m1 = mod(m,3)+1
                 m2 = mod(m1,3)+1
                 k1 = vertices(m1,k) 
                 k2 = vertices(m2,k) 
                 if(neighbour(m,k).eq.i)then
                    if(k1.eq.i1.and.k2.eq.i2.or.
     &                 k1.eq.i2.and.k2.eq.i1)go to 20
                 else if(m.eq.3)then
                 write(*,*)' 2D neighbour matrix inconsistent'
                 write(*,*)' triangle',i,' has neighbour',k,' index',j
                 write(*,*)' triangle',k,' has no neighbour',i
                 consistent = .false.
                 end if
 30           continue
              consistent = .false.
              write(*,*)
     &        ' neighbour found but vertices inconsistent'
              write(*,*)
     &        ' neighbours triangle',i,' and',k
              write(*,*)
     &        ' vertices of triangle',i,' :',(vertices(kk,i),kk=1,3)
              write(*,*)
     &        ' vertices of triangle',k,' :',(vertices(mm,k),mm=1,3)
 20        continue
 10     continue
	else if(nd.eq.3)then

	do 110 i=1,nt
	   do 120 j=1,4
              k = neighbour(j,i)
              j1 = mod(j,4)+1
              j2 = mod(j1,4)+1
              j3 = mod(j2,4)+1
              i1 = vertices(j1,i)
              i2 = vertices(j2,i)
              i3 = vertices(j3,i)
              if(k.le.0)go to 120
	      do 130 m=1,4
                 m1 = mod(m,4)+1
                 m2 = mod(m1,4)+1
                 m3 = mod(m2,4)+1
                 k1 = vertices(m1,k) 
                 k2 = vertices(m2,k) 
                 k3 = vertices(m3,k) 
                 if(neighbour(m,k).eq.i)then
                    if((k1.eq.i1.and.k2.eq.i2.and.k3.eq.i3).or.
     &                 (k1.eq.i1.and.k2.eq.i3.and.k3.eq.i2).or.
     &                 (k1.eq.i2.and.k2.eq.i1.and.k3.eq.i3).or.
     &                 (k1.eq.i2.and.k2.eq.i3.and.k3.eq.i1).or.
     &                 (k1.eq.i3.and.k2.eq.i1.and.k3.eq.i2).or.
     &                 (k1.eq.i3.and.k2.eq.i2.and.k3.eq.i1))go to 120
                 else if(m.eq.4)then
                 write(*,*)' 3D neighbour matrix inconsistent'
                 write(*,*)' tetrahedron',i,' has neighbour',k,
     &                     ' index',j
                 write(*,*)' tetrahedron',k,' has no neighbour',i
                 consistent = .false.
                 end if
 130          continue
              consistent = .false.
              write(*,*)
     &        ' neighbour found but vertices inconsistent'
              write(*,*)
     &        ' neighbours tetrahedron',i,' and',k
              write(*,*)
     &        ' vertices of tetrahedron',i,' :',i1,i2,i3,vertices(j,i)
              write(*,*)
     &        ' vertices of tetrahedron',k,' :',k1,k2,k3,vertices(m,k)
 120       continue
 110    continue
	end if
	return
	end
c------------------------------------------------------------------------

	Subroutine tetloc(x,points,vertices,neighbour,loc,out)
	real*8		points(3,*)
	real*8		x(3),del1,del2
	real*8		a1,a2,a3,b1,b2,b3
	real*8		perp(3)
	integer		vertices(4,*)
	integer		neighbour(4,*)
	integer		v1,v2,v3,v4
        logical		out
        integer		cyc(4)
        data		cyc/2,3,4,1/
        parameter(nitmx=1000)

	common/tet_loc_face/i
        common/raffaella/IERDO10

        IERDO10=0
	out = .false.
        icdo10 =0
 10     continue
        icdo10=icdo10+1
        if(out)return
        if(icdo10.gt.nitmx) then
          IERDO10=999
          return 
        endif

        do 20 i=1,4
	   i2 = cyc(i)
	   i3 = cyc(i2)
	   i4 = cyc(i3)
           v1 = vertices(i,loc)
           v2 = vertices(i2,loc)
           v3 = vertices(i3,loc)
           v4 = vertices(i4,loc)
           a1 = points(1,v3)-points(1,v2)    ! vector a=v3-v2
           a2 = points(2,v3)-points(2,v2) 
           a3 = points(3,v3)-points(3,v2) 
           b1 = points(1,v4)-points(1,v2)    ! vector b=v4-v2
           b2 = points(2,v4)-points(2,v2) 
           b3 = points(3,v4)-points(3,v2) 
           perp(1) = a2*b3-b2*a3             ! p = a x b   
           perp(2) = a3*b1-b3*a1
           perp(3) = a1*b2-b1*a2
           del1 = 0.d0
           do 11 j=1,3          ! inner product (p,v1-v2)
              del1 = del1 + perp(j)*(points(j,v1)-points(j,v2))
  11       continue 
           del2 = 0.d0
           do 12 j=1,3          ! inner product (p,x-v2)
              del2 = del2 + perp(j)*(x(j)-points(j,v2))
  12       continue 
	   if((del1.lt.0.d0.and.del2.gt.0.d0).or.
     &        (del1.gt.0.d0.and.del2.lt.0.d0))then

	      if(neighbour(i,loc).le.0)then
                 out = .true.
	      else
	         loc = neighbour(i,loc)
	      end if
	      go to 10
	   end if
 20     continue
	return
	end

      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      Implicit real*8 (A-H,O-Z)
      PARAMETER (NMAX=100,TINY=1.0E-20)
      DIMENSION A(NP,NP),INDX(N),VV(NMAX)
      D=1.d0
      DO 12 I=1,N
        AAMAX=0.d0
        DO 11 J=1,N
          IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
11      CONTINUE
        ! Kasra
        !WRITE(*,*) 'AAMAX', AAMAX
        ! End Kasra
        IF (AAMAX.EQ.0.d0) PAUSE 'Singular matrix.'
        VV(I)=1.d0/AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.d0
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*DABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1.d0/A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.d0)A(N,N)=TINY
      RETURN
      END
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      Implicit real*8 (A-H,O-Z)
      DIMENSION A(NP,NP),INDX(N),B(N)
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
c------------------------------------------------------------------------
        subroutine find_node (alat,alon,depth,nodein,ioff,walk)
 
        include '../includes/nn.param'
        include '../includes/setdel2.h'

        common/nnbuildlists/nnlist(nnsum_max),
     &                      ntlist(ntsum_max),
     &                      nn_pointer(np_max+1),
     &                      nt_pointer(np_max+1),
     &                      node(np_max)

	common/xyzswitch/xyzsw

        common/find_nodec/usefindnode

        real*8 		xd(3)
	real*4		alat
	real*4		alon
        logical		node
        logical		usefindnode
        integer		walk
        data		degtorad/0.017453292/
        data		rad/6371.0/
	if(.not.usefindnode)then
           write(*,*)' '
	   write(*,*)' Error: findnode arrays have not been built'
	   write(*,*)'        findnode can only be called if'
	   write(*,*)'        initialiation is performed in nn3d_init'
	   write(*,*)'        '
	   write(*,*)' Remedy: input parameter icall_find_node must'
	   write(*,*)'         be changed to 1 in nn3d_init.'
	   write(*,*)'         Also make sure that array size '
	   write(*,*)'         parameters ntsum_max and nnsum_max'
	   write(*,*)'         are set big enough in nn.param'
	   write(*,*)'         before recompiling'
	   write(*,*)'        '
	   stop
	end if
        walk = 1
	xyzsw = 1
        if(xyzsw.eq.1)then
          xd(1) = alat
          xd(2) = alon
          xd(3) = depth
        else
          alat2 = (90-alat)*degtorad
          alon2 = alon*degtorad
 	  xd(1) = (rad-depth)*dble(sin(alat2)*cos(alon2))
 	  xd(2) = (rad-depth)*dble(sin(alat2)*sin(alon2))
	  xd(3) = (rad-depth)*dble(cos(alat2))
        end if
        nodein = noden-ioff
        if(noden.le.0)noden = 1
        a = points(1,noden+ioff)-xd(1)
        b = points(2,noden+ioff)-xd(2)
        c = points(3,noden+ioff)-xd(3)
        dmin = a*a + b*b + c*c
 20     continue
        do i=nn_pointer(noden),nn_pointer(noden+1)-1
           n = nnlist(i)
           a = points(1,n+ioff)-xd(1)
           b = points(2,n+ioff)-xd(2)
           c = points(3,n+ioff)-xd(3)
           dist = a*a + b*b + c*c
           if(dist.lt.dmin)then
              noden=n
              dmin = dist
              walk = walk + 1
              go to 20
           end if
        end do
        nodein = noden+ioff
        return
        end
c------------------------------------------------------------------------
	Subroutine build_lists3D
     &             (np,vertices,nt,np_max,
     &              ntsum_max,nnsum_max,node,
     &              nn_pointer,nnlist,nnsum,
     &              nt_pointer,ntlist,ntsum)
	integer		vertices(4,*)
	integer		nn_pointer(*)
	integer		nnlist(*)
	integer		nt_pointer(*)
	integer		ntlist(*)
	logical		node(*)
	logical		writeout
	integer		pairs(2,6)
        common/build_lists3dc/writeout
        maxsum = 0
	pairs(1,1) = 1
	pairs(2,1) = 2
	pairs(1,2) = 1
	pairs(2,2) = 3
	pairs(1,3) = 1
	pairs(2,3) = 4
	pairs(1,4) = 2
	pairs(2,4) = 3
	pairs(1,5) = 2
	pairs(2,5) = 4
	pairs(1,6) = 3
	pairs(2,6) = 4

        do i = 1,np+1
           nt_pointer(i) = 0
           nn_pointer(i) = 0
        end do
        do i = 1,np
           node(i) = .false.
        end do
	do i=1,nt
           do j=1,4
              k=vertices(j,i)
              nt_pointer(k) = nt_pointer(k) + 1
           end do
        end do
        itemp = nt_pointer(1)
        nt_pointer(1) = 1
        do j = 2,np+1
           itemp2 = itemp
           itemp = itemp + nt_pointer(j)
           nt_pointer(j) = itemp2 + 1
        end do
        ntsum = itemp2
        if(writeout)write(*,*)' Size of ntlist packed array =',ntsum
        if(ntsum.gt.ntsum_max)then
           write(*,*)
           write(*,*)'Error: array sizes too small in '
     &     ,'subroutine build_lists3D'
           write(*,*)
           write(*,*)'       maximum size of ntlist array =',ntsum_max
           write(*,*)'       value required               =',ntsum
           write(*,*)
           write(*,*)'Increase size of parameter ntsum_max',
     &               ' and recompile'
           write(*,*)
           stop
        end if
	do i=1,nt
           do j=1,4
              k=vertices(j,i)
              ntlist(nt_pointer(k))=i
              nt_pointer(k) = nt_pointer(k)+1
           end do
        end do
        do i = np,2,-1
           nt_pointer(i) = nt_pointer(i-1)
        end do
        nt_pointer(1) = 1
        nnb = 0
        nn_pointer(1) = 1
        do i = 1,np
           nna = nnb + 1
           do j = nt_pointer(i),nt_pointer(i+1)-1
              it = ntlist(j)
              do k=1,6
                 i1 = pairs(1,k)
                 i2 = pairs(2,k)
                 n1 = vertices(i1,it)
                 n2 = vertices(i2,it)
                 if(n1.eq.i.and..not.node(n2))then
                    node(n2) = .true.
                    nnb = nnb + 1
                    nnlist(nnb) = n2 
                 else if(n2.eq.i.and..not.node(n1))then
                    node(n1) = .true.
                    nnb = nnb + 1
                    nnlist(nnb) = n1 
                 end if
              end do
           end do
           nn_pointer(i+1) = nnb+1
           do k=nna,nnb
              node(nnlist(k)) = .false.
           end do
           if(nnb.gt.nnsum_max)then
              write(*,*)'Error: array sizes too small in '
     &        ,'subroutine build_lists3D'
              write(*,*)'       maximum size of nnlist array =',
     &        nnsum_max
              write(*,*)'       value required is at least   =',nnb
              write(*,*)'       reached while processing node ',i
              write(*,*)'       Total number of nodes        =',np
              write(*,*)' Increase size of parameter nnsum_max'
              stop
           end if
        end do
        nnsum = nnb

        if(writeout)write(*,*)' Size of nnlist packed array =',nnsum
	return
	end
c------------------------------------------------------------------------
	Subroutine build_hullneighbour
     &             (nh,vertices,nt,nwork2d,
     &              neighbour,nnn,nnlist,ntrilist)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		nnn(*)
	integer		nnlist(*)
	integer		ntrilist(*)
	logical		nnwrite

	common/nnswitches/nnwrite,lud

	if(nnwrite)write(*,*)' Building neighbour v ...'
	do 5 i = 1,3
	   do 5 j = 1,nt
	      neighbour(i,j) = 0
 5      continue
        do 6 i = 1,nwork2d
	   nnlist(i) = 0
	   ntrilist(i) = 0
 6      continue
        do 7 i = 1,nh
	   nnn(i) = 0
 7      continue
	do 10 it = 1,nt
	   i1 = vertices(1,it)
	   i2 = vertices(2,it)
	   i3 = vertices(3,it)
	   nnn(i1) = nnn(i1) + 1
	   nnn(i2) = nnn(i2) + 1
	   nnn(i3) = nnn(i3) + 1
 10     continue
	itemp = nnn(1)+1
	nnn(1) = 1
	do 20 j = 2,nh+1
	   itemp2  = itemp 
	   itemp   = itemp + nnn(j)+1
	   nnn(j) = itemp2 + 1
 20     continue
	if(nnn(nh+1).ge.nwork2d)then
           write(*,*)' '
           write(*,*)' Error: array sizes too small in subroutine'
     &               ,'build_hullneighbour'
           write(*,*)'       maximum number of neighbours for all nodes'
           write(*,*)'       on the hull is too small: current value ='
     &     ,nwork2d
           write(*,*)'       Increase size of parameter nwork2d'
           write(*,*)'       to at least',nnn(nh+1)
           write(*,*)'       An upper limit for nwork2d is 3 times'
           write(*,*)'       the maximum number of faces on the hull' 
           write(*,*)'       nwork2d = 3*nf_max + nh_max' 
           write(*,*)' '
	   stop
	end if
	do 25 it = 1,nt
	   i1 = vertices(1,it) 
	   i2 = vertices(2,it) 
	   i3 = vertices(3,it) 
	   j1 = nnn(i1)
	   j2 = nnn(i1+1) - 1 
	   jt = 0
	   do 30 j = j1,j2
	      if(nnlist(j).eq.0)then
	         nnlist(j) = i2
	         ntrilist(j) = it
	         go to 31
	      else if(nnlist(j).eq.i2.and.ntrilist(j).ne.it)then
                 jt = ntrilist(j)
	         go to 31
	      end if
  30       continue
  31       continue
	   if(jt.eq.0)then
	      j1 = nnn(i2)
	      j2 = nnn(i2+1) - 1 
	      do 32 j = j1,j2
	         if(nnlist(j).eq.0)then
	            nnlist(j) = i1
	            ntrilist(j) = it
	            go to 33
	         end if
  32          continue
	   end if
  33       continue

	   if(jt.ne.0)then
	      neighbour(3,it) = jt
	      k1 = vertices(1,jt)
	      k2 = vertices(2,jt)
	      k3 = vertices(3,jt)
	      if(k1.ne.i1.and.k1.ne.i2)then 
	         neighbour(1,jt) = it
	      else if(k2.ne.i1.and.k2.ne.i2)then 
	         neighbour(2,jt) = it
	      else
	         neighbour(3,jt) = it
	      end if
           end if
	   jt = 0
	   j1 = nnn(i1)
	   j2 = nnn(i1+1) - 1 
	   do 130 j = j1,j2
	      if(nnlist(j).eq.0)then
	         nnlist(j) = i3
	         ntrilist(j) = it
	         go to 131
	      else if(nnlist(j).eq.i3.and.ntrilist(j).ne.it)then
                 jt = ntrilist(j)
	         go to 131
	      end if
  130      continue
  131      continue
	   if(jt.eq.0)then
	      j1 = nnn(i3)
	      j2 = nnn(i3+1) - 1 
	      do 132 j = j1,j2
	         if(nnlist(j).eq.0)then
	            nnlist(j) = i1
	            ntrilist(j) = it
	            go to 133
	         end if
  132         continue
	      end if
  133      continue
	   if(jt.ne.0)then
	     neighbour(2,it) = jt
	     k1 = vertices(1,jt) 
	     k2 = vertices(2,jt)
	     k3 = vertices(3,jt)
	     if(k1.ne.i1.and.k1.ne.i3)then 
	        neighbour(1,jt) = it
	     else if(k2.ne.i1.and.k2.ne.i3)then 
	        neighbour(2,jt) = it
	     else
	        neighbour(3,jt) = it
	     end if
           end if
	   jt = 0
	   j1 = nnn(i2)
	   j2 = nnn(i2+1) - 1 
	   do 230 j = j1,j2
	      if(nnlist(j).eq.0)then
	         nnlist(j) = i3
	         ntrilist(j) = it
	         go to 231
 	      else if(nnlist(j).eq.i3.and.ntrilist(j).ne.it)then
                 jt = ntrilist(j)
 	         go to 231
	      end if
  230      continue
  231      continue
	   if(jt.eq.0)then
	      j1 = nnn(i3)
	      j2 = nnn(i3+1) - 1 
	      do 232 j = j1,j2
	         if(nnlist(j).eq.0)then
	            nnlist(j) = i2
	            ntrilist(j) = it
	            go to 233
	         end if
  232         continue
	      end if
  233      continue
	   if(jt.ne.0)then
	     neighbour(1,it) = jt
	     k1 = vertices(1,jt) 
	     k2 = vertices(2,jt)
	     k3 = vertices(3,jt)
	     if(k1.ne.i2.and.k1.ne.i3)then 
	        neighbour(1,jt) = it
	     else if(k2.ne.i2.and.k2.ne.i3)then 
	        neighbour(2,jt) = it
	     else
	        neighbour(3,jt) = it
	     end if
           end if

 25     continue

	if(nnwrite)write(*,*)' built hullneighbour v'

	return
	end

!     --------------- END OF DELAUNEY STUFF -----------------------

