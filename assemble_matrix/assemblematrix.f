      program assemblematrix

      ! Compile: ifort -o $target assemblematrix.f
      ! Or the cautious way for very large matrices:
      ! ifort -mcmodel=large -shared-intel -o $target assemblematrix.f
      ! gfortran compiler produces erratic output, 
      ! see below (2010/08/05)!!!
      ! 
      ! GN 2006.
      ! Revision log:
      ! Oct 06: Output column sum instead of log(A'A/max(A'A)) to file columndensity. 
      !         Introduced kssw(4) switch to force
      !         correction for Q-dispersion for period(iband).
      !         Output data vector for diagnostic
      !         Add output of data error to matrix file
      !         Solve the demean problem by adding a scratch aux. file
      !         Change in cfevt to get correct nevent when discriminate
      !         among clusters

      ! *********** Experimental version under testing ************
      ! Does not yet handle:
      !  linked matrices for differential time data
      !  split output for parallel version of solvetomo
      ! 
      ! 2006/11/17 K.S.: added option to mask out (reject) data
      ! based on criteria specified in a file in.mask
      ! 2007/2/4 G.N: Added data groups, weighted equally as a group in 
      ! the total chi square misfit, or with weights specified by user.
      !
      ! 2010/04/08 Doubled the size of some character buffers 
      ! (fname2, ident_in). They must be twice as long as
      ! fname and ident, else runtime error for both gfortran
      ! and ifort in Munich. dataf2 must be no longer than 30 chars
      ! (=its size in raymatrix), else runtime error when reading first
      ! line of matrixT file in the unformatted case.
      ! 2010/04/08 K.S. changed format 120, 3i4 --> 3i7 to prevent overflow.
      ! CAUTION: Needed to make the same adjustment 
      ! in mpisolvetomo.f and solvetomo.f
      ! (but not resolutiontest.f, since free format read there).
      ! 2010/04/08 K.S. renamed "common /mod/" to "common /modl/", since
      ! ifort gives a compile error:
      ! 97> ifort -o $target assemblematrix.f
      ! assemblematrix.f(551): error #8038: A referenced intrinsic procedure 
      ! can not have the same name as a common block.   [MOD]
      !       if(mod(nrow,100).eq.1) write(6,fmt='(a1,$)') 'x' !XXX  
      ! ---------^
      ! compilation aborted for assemblematrix.f (code 1)
      ! Change nowhere else: no other routine uses this common block,
      ! all other "mod" occurrences refer to "modulo" function.
      ! * 2010/04/08 K.S. Changed "fname2" to "dataf2" 
      !    write(6,'(2a)') 'Appears to be linked to: ',dataf2
      !    write(6,'(2a)') 'Appears to be linked to: ',dataf2
      ! Must have been a bug (fname2 not initialized) that 
      ! would only have materialized for linked data files.
      ! * 2010/04/08 K.S. increased max allowable matrix dimensions to
      ! NMAX=2000000,NDIM=250000 -- this is consistenst with
      ! the current values of includes/solvet.h in mpisolvetomo.f

      ! * 2010/08/03 K.S. increased max allowable unknowns NDIM
      ! to 350,000 since new global grid G30 hast 294,000 unknowns alone.
      ! Also changed NDIM in includes/solvet.h
      ! Note that this will still not suffice for joint T&A inversions.
      ! Also not that gfortran compiler for unknown reasons did
      ! NOT print the error message "ncol>NDIM. Increase NDIM." but
      ! bailed out silently! (ifort gave the message.)

      ! 2010/08/26 K.S. gfortran works ONLY with option -finit-local-zero
      ! gfortran -std=legacy -o $target assemblematrix.f ! BUGGY!!!
      ! gfortran -finit-local-zero -std=legacy -o $target assemblematrix.f !!OK
      !
      ! ifort works:
      ! ifort -o $target assemblematrix.f  ! works
      ! Precautionry compile, Jens advice, in case large common 
      ! blocks/static are used.
      ! ifort -mcmodel=large -shared-intel -o $target assemblematrix.f
      !
      ! 2010/08/26 K.S. now using include '../includes/solvet.h'
      ! instead of parameter(MMAX=200000,NMAX=2000000,NDIM=350000)
      ! (Also changed that in resolutiontest.f)

      ! TROUBLE TICKET K.S. 2010/08/05:
      ! gfortran produces erratic runtime results, there seems to be
      !  a memory issue.
      ! gfortran -std=legacy -o $target assemblematrix.f
      ! test case: newident = test.G10T.BB.USA10 (full BB global 
      ! data set on old grid.)
      ! gfortran produces a different matrix in every run, best seen 
      ! in the final output line "Largest element of A'A="
      ! Sometime output is so bad that program bails out.
      ! Output of a few subsequent runs:
      ! Largest element of A'A=   18.111338     , sqrt=   4.2557416    
      ! Largest element of A'A=  1.20876095E+25 , sqrt=  3.47672386E+12
      ! Largest element of A'A=  5.12861225E-41 , sqrt=  7.16143328E-21
      ! Largest element of A'A=   10806.03     , sqrt=   103.9520   
      ! Large or infinite ata in column       74784
      ! Largest element of A'A=   0.0000000     , sqrt=   0.0000000  

      ! ifort (with or without options -mcmodel=large -shared-intel)
      ! seems to not have this problem and produces consistent results:
      ! Largest element of A'A=   10806.03     , sqrt=   103.9520  


      !================================================================
      ! 
      ! assemblematrix is the fourth step in the tomographic inversion:
      ! [1] run raydata to construct a finite-frequency data file
      ! [2] run springs3d.f to parameterize the model
      ! [3] run raymatrix.f to build the raw matrix for ray data.
      ! [3a] run diffmatrix.f for differential data, if any
      ! [4] run assemblematrix to combine all data & raw matrices
      ! [5] run solvetomo to get a tomographic solution

      ! assemblematrix reads one or more matrix files, combines them
      ! into one large matrix  and adds to that
      ! the columns for various station- and event corrections

      ! adding rows to regularize the inversion is left for
      ! solvetomo.f, to retain full flexibility

      ! The matrix is stored on disk, one row at a time, with the
      ! rhs (the datum) and the row divided by the standard deviation

      ! matrix scaling:
      ! ALL data and rows are divided by the standard deviation <s> for the datum
      ! If data groups are defined, data and matrix rows are multiplied by
      !    <chiw(igroup)>= N/(Ni*Ng), where N=total nr of data, Ni number of data
      !    in group i, Ng number of data groups.
      ! If the weighting option is chosen, chiw(i) is specified by the user.
      ! Giving weigh=1 for every group mean each *datum* (not group) will
      ! receive the same weight.

      ! INPUT:
      ! interactive input of parameters and file names (in.assemblematrix)
      ! assemblematrix.mask.<ident>   (optional) 2-column file provided by
      !           user. Contains criteria to mask (reject)
      !           data. These are used to screen the info read from files
      !           raymatrix.info.T.<ident> or raymatrix.info.T.<ident>
      !           Example format:
      !           IGNORE_BAND 2
      !           IGNORE_BAND 4
      !           XCORR_MIN   0.90
      !
      ! OUTPUT:
      ! mat.<ident>      assembled matrix file (binary), input for solvetomo.f
      ! aux.<ident>      auxiliary file (ascii),         input for solvetomo.f 
      ! ata.<k>.<ident>  diagonal of  A'A matrix         input for solvetomo.f 
      !                  k=1 for Vp, 2 for Vs, 3 for Q
      !
      ! assemblematrix.out.<ident>      main diagnostics file, ascii
      ! data.<ident>     ascii, two columns, ndata rows:
      !                  absolute value (ddif), normalized value (rhs)
      !
      ! dum.<ident>       scratch file for matrix output, written before mat.<ident>
      !                   (temporary binary output)
      ! dumaux.<ident>    scratch file for auxiliary output, written before aux.<ident>
      !                  (also temporary ascii output for debugging)
      !
      ! inmapc.<k>:        input for k'th kernel density map
      ! columndensity.<k>.<ident> input file for map.f, scaled ata.1.P 
      !                   (nparms,4) entries:
      !             node x,y,z,sum(column(A))/max(A),
      !             i.e. relative sensitivity at x,y,z
      ! runmapc:           shell script to plot A column sum (spatial kernel density)
      !
      ! tex.<ident>       tex file contains some station/src statistics in tables
      !
      ! File units:
      ! 1 : dum.<ident>
      ! 2 : aux.<ident>
      ! 3 : mat.<ident>
      ! 4 : assemblematrix.out.<ident>
      ! 7 : tex.<ident>; columndensity.[1,2,3].<ident>
      ! 8 : inmapc.[1,2,3]; ata.[1,2,3].<ident>; atamx.<ident>
      ! 9 : runmapc
      ! 10: data.<ident> 
      ! 11: matrixT, matrixA input files (binary)
      ! 12: dumaux.<ident>
      ! 19: assemblematrix.mask.<ident>
      ! 20: assemblematrix.stations.<ident>  NEW: station statistics
      ! 23: assemblematrix.mask.<ident>  IN: optional, user-provided file

      ! header files for Delauney triangulation
      include 'includes/nn.param'
      include 'includes/setdel3.h'       ! common for np,nt (Delauney)
      include 'includes/setdel2.h'       ! common for points,centres,vertices etc
      include 'includes/solvet.h'        ! for matrix dims MMAX, NMAX, NDIM
      ! matrix rows: allow for 
      ! MMAX nonzero elements in one row, 
      ! MMAX data (rows)
      ! NDIM unknowns (models parameters + corrections)
      ! (should be consistent with mpisolvetomos solvet.h)
      ! TODO: line below should be replaced by include "../includes/solvet.h"
c      parameter(MMAX=200000,NMAX=2000000,NDIM=350000)   ! 2010/08/04 K.S.
c      parameter(MMAX=200000,NMAX=2000000,NDIM=250000)  ! 2010/04/09 K.S.
c      parameter(MMAX=200000,NMAX=1000000,NDIM=200000)  ! 2007/04/28
c      parameter(MMAX=200000,NMAX=1000000,NDIM=100000)
      parameter(NFR=20)
      dimension asparse(MMAX),ja(MMAX),colsum(3,np_max),ata(NDIM)
      dimension colsumy(3,np_max)  !!! one row

      ! data
      dimension rhs(NMAX)

      ! other
      dimension jssw(11),jlast(12),nssw(3),kssw(4),kountpar(3),dpar(11)
      dimension period(NFR)
c      dimension rmsb(NFR)  ! rms noise in bands
      integer   evlist(500000),klustlist(500000),mevt(500000,2)
      integer   mstat(500000,2)
      dimension evlat(500000),evlon(500000),evdep(500000)
      dimension statlat(500000),statlon(500000),statelev(500000)
      dimension statsum(500000,2),evtsum(500000,2),chiw(100)
      integer   nchi(100)
      integer (kind=8) kountall, kountkerall

      ! various character variables 

      character*72  fname,ident, chbuf ! nchar long
      character*144 fname2, ident_in   ! at least 2*nchar long!
      character*30  dataf2             ! must be <=30 long, as in raymatrix
      character*16 stationcode,stnam,fname3
      character*16 statlist(500000)
      character*8  networkcode,netw
      character*8  netwlist(500000)
      character*3  comp,comp2
      character*32 chcor(11),chcorrection(4)
      character*60 nodesfile
      common/nodes/nodesfile

      ! mask variables (K.S. 2006/11/17)
      logical lmask,reject,ffdatum
      integer mb0,mba,mbt,ntot,nrej
      integer m_bands(NFR),m_ampbands(NFR),m_ttbands(NFR)
      real xcorr_min
      character*256 maskfile 

      
      data chcor/'01:dlnVp',
     &'02:dlnVs',
     &'03:dlnQs',
     &'04:Hypoctr corr (km)',
     &'05:Station corr t_P (s)',
     &'06:Station corr dlnA_P',
     &'07:Origin time corr dTo(s)',
     &'08:Event moment corr dlnA', 
     &'09:Station corr t_S (s)',
     &'10:Station corr dlnA_S',
     &'11:TEST HYPOCENTER'/
      data chcorrection/'01:ellipticity',
     &'02:curst',
     &'03:station elevation',
     &'04:Q-dispersion'/
c     data myrow/1/  !!! one row

! background model (qvp,qvs and qqs are ignored)
      common /modl/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &      Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      write(6,fmt='(a)') 'Running assemblematrix...' 

      iu4 = 4  ! I/O unit for assemblematrix.out.<ident> log file 

      ! Normally jdebug=0. Set -1 if input and output files are ascii
      ! instead of binary. Set 1 for some debugging output.
      jdebug=0

      ! set universal constants at machine precision
      pi=4.0*atan(1.0)
      twopi=2.0*pi
      halfpi=0.5*pi
      r2d=180.0/pi
      sqrt3=sqrt(3.0)

      ! input file with model parametrization
      print *,'Input 3D model grid...'
      call inpgrd
      print *,'3D model has ',np,' grid points'
      npold=np

      ! Read parameters sense switches and uncertainties
      print *,'Input sense switches and prior uncertainties'
      print *,'Parameter switch to be inverted for are set to 1, else 0'
      print *,'[extra: if jssw(9)=-1 then tSstat=sqrt(3)*tPstat,'
      print *,'if jssw(7)=2, then dTo for each cluster separately,'
      print *,'if jssw(8)=2  then dMo for each cluster separately].'
      print *,'  Parameter         jssw dpar'
      do i=1,11
        ! print *, chcor
        write(6, fmt='(/,i2,1x,a32,":",$)') i,chcor(i)
        read(5,*) jssw(i),dpar(i)
        print *, 'on/off', jssw(i), 'value', dpar(i)
      enddo  
      print *,' '
      
      invstat=0
      invevt=0
      klustercor=0
      nallgr=0
      if(jssw(5).eq.1.or.jssw(6).eq.1.or.jssw(9).gt.0) invstat=1
      if(jssw(7).eq.1.or.jssw(8).eq.1) invevt=1
      ! to avoid coupling S and P time corrections
      if(jssw(5).eq.0.and.jssw(9).eq.-1) then 
        jssw(9)=0
        print *,'Warning: no S station correction'
      endif
      if(jssw(7).gt.1.or.jssw(8).gt.1) klustercor=1

      ! determine how many classes of unknowns we have 
      ninvt=0
      do i=1,11
        if(jssw(i).eq.1) ninvt=ninvt+1
        jlast(i)=0
      enddo
      print *,ninvt,' classes of unknowns'
      if(ninvt.eq.0) stop 'error in sense switches'

      ! The column for node j for parameter k is in j+jlast(k)
      ! so jlast(1), for Vp,  remains 0 always
      jlast(2)=0
      if(jssw(1).eq.1) jlast(2)=np           ! jlast for Vp
      jlast(3)=jlast(2)                     
      if(jssw(2).eq.1) jlast(3)=jlast(2)+np  ! jlast for Vs  
      jlast(4)=jlast(3)
      if(jssw(3).eq.1) jlast(4)=jlast(3)+np  ! jlast for Qs

      ! Read correction sense switches (ellipticity, crust,
      ! elevation, Q-dispersion, 1=on, 0=off)
      print *,'Give travel time correction switches: 1=ellipticity,'
      print *,'2=crust,3=station elevation, 4=Q-dispersion'
      print *,'(1=on, 0=off).'
      print *,'[more: If kssw(4)=0, no dispersion correction.'
      print *,'If kssw(4)=1, apply dispersion correction,'
      print *,'reference period=1s.'
      print *,'If kssw(4)=2, correct data to be w.r.t. 1s and'
      print *,'give data-used reference period in the next line.'
      print *,'In this case,data already include dispersion '
      print *,'correction, but NOT w.r.t. 1s.]'
      print *,'e.g. 1 1 1 0 :'
      read *,kssw
c     if(kssw(4).ne.0) 
c    &      print *,'WARNING: Q correction normally not needed',
c    &              ' for cross-correlated delays'
      if(kssw(4).eq.2)  then
        read *, refperiod
c       if (refperiod.gt.10 .or. refperiod.lt.0.1) then
c         print *,'Unreasonable reference period!  kssw(r)=2, '
c         print *,'make sure there is input of reference period.'
c         stop
c       endif
      endif

      ! Get ident (postfix to filenames) for output files
      print *,'Give ident for output matrix file:'
      read(5,fmt='(a)') ident
      print *
      
      ! Check for presence of a mask file assemblematrix.mask.<ident>
      ! (optional file provided by user)
      ! The same mask file gets used for all matrix input files!
      maskfile = 'assemblematrix.mask.'//ident
!      call read_maskfile(maskfile,lmask,m_bands,mb0,xcorr_min,
!     & NFR)
      call read_maskfile(maskfile,lmask,m_bands,mb0,xcorr_min,
     & m_ttbands,mbt,m_ampbands,mba,NFR)


! initialize
      do i=1,NDIM
        ata(i)=0.      ! diagonal of matrix A(Transposed)A
      enddo  
      do i=1,np_max    ! column sums for at most three model parameters
        do j=1,3  
          colsum(j,i)=0.
c         colsumy(j,i)=0.  !!! one row
        enddo
      enddo  
      kklust=0          ! cluster identifier
      kountall=0        ! kountall counts nonzero elements of matrix
      kountkerall=0
      kstat=0           ! station identifier
      lastdemean=0      ! last data that has been demeaned
      lastgroup=0       ! last data that has been grouped
      lastnrow=0        ! last data in previous matrixT/matrixA file
      networkcode=' '
      nevent=0          ! total number of events (hypocenters)
      nfile=0           ! counts number of matrix files assembled
      ngroup=0          ! number of data groups
      nklust=0          ! total number of clusters
      nosigma=0         ! counts number of missing error estimates
      nstat=0           ! total number of station (for corrections)
      stationcode=' '
      nrow=0            ! nrow counts total number of rows (accepted data)
      nrej=0            ! number of rejected data (masked out)
      ntot=0            ! number of total data (nrow+nrej)
      
! open the assembled matrix file for output
      print *,'Assembled matrix file is: ','mat.'//ident
      ! we write the matrix in two steps, using a scratch file first 
      ! (for both matrix file and auxiliary file).
      ! this allows us to write all correction columns in blocks of
      ! the same correction type (and do demean to the data). Only the
      ! hypocentral locations are directly added in the first step.
      if(jdebug.eq.0) then
        open(1,file='dum.'//ident,form='unformatted')
        open(333,file='numrc.'//ident)
        open(444,file='nonzr.'//ident)
        open(555,file='colindx.'//ident)
        open(666,file='colval.'//ident)
      else
        open(1,file='dum.'//ident)
      endif
      open(10,file='data.'//ident) ! output data vector
      print *,'Auxiliary file is: ','aux.'//ident
      open(12,file='dumaux.'//ident)

! open output files
      open(7,file='tex.'//ident)
      open(iu4,file='assemblematrix.out.'//ident)
      write(iu4,fmt='("START_DATA_INPUT")')

      kweighflag=0
      kgroupflag=0

      !--------------------------------------------------------------
      ! Read next matrix file name or group/demean/stop command
      !--------------------------------------------------------------
10    print *; print *,'Give next matrix file name (eg matrixT.P)'
      print *,'or <demean>,<weigh>,<group> or <stop>:'
      print *, fname
      read(5,fmt='(a)') fname

      ! Case 1: Stop -- done reading input
      if(fname.eq.'stop')then
        write(iu4,fmt='("STOP_DATA_INPUT")')
        write(iu4,fmt='(a8,1x,i10)') 'NFILES',nfile
        write(6,*) 'Done reading input.'
        goto 300
      endif

      ! Case 2: demean all data read since last "demean" statement 
      if(fname.eq.'demean') then    
        sum=0
        do i=lastdemean+1,nrow
          sum=sum+rhs(i)
        enddo
        average=sum/(nrow-lastdemean)
        do i=lastdemean+1,nrow
          rhs(i)=rhs(i)-average
        end do
        write(iu4,129) 'DEMEAN', nrow-lastdemean,average
129     format(a8,i10,e12.3)
        lastdemean=nrow
        goto 10  ! get next file name or command
      endif 

      ! Case 3: group together all data since last "group" statement
      !         Each group gets same overall weight (chiw is computed later).
      !         In absence of any "group" statement, chiw=1 for all data.
      if(fname.eq.'group') then
        ngroup=ngroup+1
        kgroupflag=1
        if(kweighflag.gt.0) stop 'Grouping+weighting not allowed'
        if(ngroup.gt.100) stop 'Only 100 data groups allowed'
        nchi(ngroup)=nrow-lastgroup
        nallgr=nallgr+nchi(ngroup)
        chiw(ngroup)=0.
        write(iu4,129) 'GROUP', nrow-lastgroup
        lastgroup=nrow
        goto 10
      endif

      ! Case 4: apply same user-defined weight chiw to all data since last "weigh"
      !         Each datum gets the specified weight.
      if(fname.eq.'weigh') then
        ngroup=ngroup+1
        kweighflag=1
        if(kgroupflag.gt.0) stop 'Grouping+weighting not allowed'
        if(ngroup.gt.100) stop 'Only 100 data groups allowed'
        nchi(ngroup)=nrow-lastgroup
        nallgr=nallgr+nchi(ngroup)
        print *,'Give weight for this set of data:'
        read *,chiw(ngroup)
        write(iu4,129) 'WEIGH', nrow-lastgroup,chiw(ngroup)
        lastgroup=nrow
        goto 10  
      endif  
      
      ! Case 5: Read next matrix file matrixT.* or matrixA.*
      nfile=nfile+1

      kdatatype=0
      if(fname(7:7).eq.'T') kdatatype=1
      if(fname(7:7).eq.'A') kdatatype=2
      if(kdatatype.eq.0) then
        print *, 'Data type not recognized:'
        print *, 'File name must start with matrixX, where X is T or A.'
        stop 
      endif  
      
      print *;
      if(jdebug.eq.0) then      ! unformatted
        print *,'Matrix files must be unformatted (binary)'
        print *,'Opening matrix file ',fname
        open(11,file=fname,status='old',form='unformatted')
        read(11,iostat=ios) linked,dataf2   ! linked file if linked=1
        if(ios.ne.0) stop 'Error in line 1'
        if(linked.eq.1) then
          write(6,'(2a)') 'File: ',fname
          write(6,'(2a)') 'Appears to be linked to: ',dataf2
          write(6,*) "First run linkmatrix to get a difference file"
          stop 'Error in matrix file'
        endif  
        read(11,iostat=ios) np,nssw     ! nr of grid points, par types
        if(ios.ne.0) stop 'Error in line 2'
        read(11,iostat=ios) nband,(period(i),i=1,nband)
        if(ios.ne.0) stop 'Error in line 3 for period'
      else  
        print *,'Matrix files must be formatted (ascii)'
        print *,'jdebug=',jdebug
        print *,'Opening matrix file ',fname
        open(11,file=fname,status='old')
        read(11,fmt='(i2,1x,a)',iostat=ios) linked,dataf2 
        if(ios.ne.0) stop 'Error in line 1 for linked file'
        if(linked.eq.1) then
          write(6,'(2a)') 'File: ',fname
          write(6,'(2a)') 'Appears to be linked to: ',dataf2
          write(6,*) "First run linkmatrix to get a difference file"
          stop 'Error in matrix file'
        endif  
        read(11,*,iostat=ios) np,nssw   ! nr of grid points, par types
        if(ios.ne.0) stop 'Error in line 2 for np, nssw'
        read(11,*,iostat=ios) nband,(period(i),i=1,nband)
        if(ios.ne.0) stop 'Error in line 3 for period'
      endif
      if(np.gt.np_max) stop 'Increase np_max'

      ! check if the file is for the same 3D grid by comparing the
      ! number of nodes
      if(np.ne.npold) stop 'np incompatible with 3D grid file'

      ! Open associated info/mask file: raymatrix.info.T.<ident_in>
      if(kdatatype.eq.1) write(ident_in,118) fname(9:) ! fname = 'matrixT.<ident_in>'
118   format('raymatrix.info.T.',a)
      if(kdatatype.eq.2) write(ident_in,119) fname(9:) ! fname = 'matrixA.<ident_in>'
119   format('raymatrix.info.A.',a)
      
      iu19 = 19
      write(6,fmt='(a)') 'Opening data info file (for mask) ',ident_in      
      open(iu19,file=ident_in,status='unknown')
      read(iu19,*,iostat=ios)  ! one header line
      if(ios.lt.0)then
         write(6,fmt='(a)') '-------------  WARNING!  ---------------'
         write(6,fmt='(a)') 'Data info file non-existent or empty.'
         write(6,fmt='(a)') 'Filename: ', ident_in         
         write(6,fmt='(a)') 'Using ALL data in absence of mask file.'
      endif   
      
      
      if(jdebug.eq.1) write(13,*) 'Now reading ',fname
      write(6,fmt='(a,$)') 'Progress: '   

      ! read one matrix row, rhs element & corrections
      ! data are uniquely identified by combination of ndatum and iband
20    if(jdebug.eq.0) then                    ! unformatted  
        read(11,end=200,iostat=ios) ndatum,kount,iband,dobs,dtheor,s,
     &    ecorr,corcrust,corelev,tstar,qray,ievt,kluster,stationcode,
     &    networkcode,comp,slowness,raz,rdel,angle0,vsrc,krtyp
        if(kount.gt.MMAX) stop 'kount>MMAX, increase dimensions'
        if(ios.ne.0) stop 'Error in matrixfile - header'
        read(11,end=200,iostat=ios) kountpar
        if(ios.ne.0) stop 'Error in matrixfile - kountpar'
        read(11,end=200,iostat=ios) (ja(i),i=1,kount)
        if(ios.ne.0) stop 'Error in matrixfile - ja'
        read(11,end=200,iostat=ios) (asparse(i),i=1,kount)
        if(ios.ne.0) stop 'Error in matrixfile - asparse'
      else  ! formatted
        read(11,*,end=200,iostat=ios) ndatum,kount,iband,dobs,dtheor,s,
     &  ecorr,corcrust,corelev,tstar,qray,ievt,kluster,stationcode,
     &    networkcode,comp,slowness,raz,rdel,angle0,vsrc,krtyp
        write(13,*) 
     &  ndatum,kount,iband,dobs,dtheor,s,
     &  ecorr,corcrust,corelev,tstar,qray,ievt,kluster,stationcode,
     &    networkcode,comp,slowness,raz,rdel,angle0,vsrc,krtyp
        if(ios.ne.0) stop 'Error in matrixfile - header'
        if(kount.gt.MMAX) stop 'kount>MMAX, increase dimensions'
        read(11,*,end=200,iostat=ios) kountpar
        if(ios.ne.0) stop 'Error in matrixfile - kountpar'
        read(11,*,end=200,iostat=ios) (ja(i),i=1,kount)
        if(ios.ne.0) stop 'Error in matrixfile - ja'
        read(11,*,end=200,iostat=ios) (asparse(i),i=1,kount)
        if(ios.ne.0) stop 'Error in matrixfile - asparse'
      endif  

      ntot = ntot+1
      
c...  ! MASKING SECTION    
      ! Read info for masking: one line from raymatrix.info file
      read(iu19,*,iostat=ios) ndatum2,ievt2,itmp,chbuf,iband2,
     & stnam,netw,rlat,rlon,relev,delta,azi,xco,rmsbiband,kunit,
     & comp2,idate,iotime,slat,slon,sdep 

c-- TODO 2007/03/02: reinstate the below as soon as no old raymatrix output files are around any longer
      if(ndatum.ne.ndatum2.or.iband.ne.iband2)then
        print *, 'ERROR: discrepancy between matrix file and'
        print *, 'its corresponding mask file: '
        print *, 'ndatum,ndatum2,iband,iband2'
        print *, ndatum,ndatum2,iband,iband2
        stop
      endif  

      ! Apply mask(s) to this datum
      reject = .false.
      if(lmask)then
        ! Cross-correlation mask
        ! (raytheoretical ttimes have xco==0. Avoid rejecting them.)
        ffdatum = .false.
        if(xco.lt.-1e-6.or.xco.gt.1e-6) ffdatum=.true.
        if(xco.lt.xcorr_min.and.ffdatum) reject=.true.
        ! Frequency band masks
        if(mb0.gt.0)then
          do kk=1,mb0
            if(iband.eq.m_bands(kk)) reject=.true.
          enddo   
        endif
        if(mbt.gt.0)then
         do kk=1,mbt
          if(iband.eq.m_ttbands(kk).and.kdatatype.eq.1) reject=.true.
         enddo   
        endif
        if(mba.gt.0)then
         do kk=1,mba
          if(iband.eq.m_ampbands(kk).and.kdatatype.eq.2) reject=.true.
         enddo   
        endif
        ! Add other masks here...
        ! TODO: mask based on DELTA
        ! TODO: mask based on receiver and source locations (regional inversions)
        ! western_bound, eastern_bound, northern_bound, southern_bound
      endif

      ! if(abs(dobs-dtheor) > 10) then
      !     reject=.true.
      !     print *, dobs-dtheor
      ! endif

      if(reject)then
         nrej = nrej+1
         goto 20    ! skip to next datum
      endif   ! if(lmask)
c...  ! END MASKING SECTION     
       
      nrow=nrow+1   ! data counter (accepted data)
      if(mod(nrow,100).eq.1) write(6,fmt='(a1,$)') 'x' !XXX  
      if(nrow.gt.NMAX) stop 'Too many data, increase NMAX'

      
      ! apply time corrections (added to data predictions)
      if(kdatatype.eq.1) then    ! travel time delays
        if(kssw(1).eq.1) dtheor=dtheor+ecorr    ! ellipticity
        if(kssw(2).eq.1) dtheor=dtheor+corcrust ! crust2.0
        if(kssw(3).eq.1) dtheor=dtheor+corelev  ! station elevation
        ! Apply Q-dispersion correction, reference period is 1 second
        if(kssw(4).eq.1.and.iband.gt.0) 
     &        dtheor=dtheor+tstar*log(period(iband))/pi
        ! NEW 8/13/2007, kssw(4)=2: correct data to be w.r.t. 1 second 
        ! In this case, data already include dispersion correction, 
        ! but NOT w.r.t. 1 s period. We ignore the frequency dependence
        ! of Q even for period<1s because the Amplitude Project uses
        ! a constant Q.
        if(kssw(4).eq.2.and.iband.gt.0) 
     &        dtheor=dtheor+tstar*log(refperiod)/pi
        ddif=dobs-dtheor
      else if(kdatatype.eq.2) then  ! no amp corrections envisaged
        ddif=dobs
      endif
      
      if(s.le.0.) then
        s=1.0
        nosigma=nosigma+1
      endif  
 
      ! always normalize the data vector to unit standard deviation
      rhs(nrow)=ddif/s          

      ! remap ja to column in assembled matrix, normalize to sigma=1
      ! and scale parameters with dpar (so expected value ~ 1)
      ! Note: in raymatrix the columns are always 1 to np for Vp, np+1 to
      ! 2*np for Vs, 2+np=1 to 3*np for Q. In the assembled matrix we
      ! condense these, only count parameters that are inverted for. Thus
      ! if we invert Vs and Q, 1 to np is for Vs, np+1 to 2*np is for Q.
      do i=1,kount
        kpartype=(ja(i)-1)/np+1       ! 1 for Vp, 2 for Vs, 3 for Q 
        if(kpartype.gt.3) stop 'kpartype > 3'
        ja(i)=ja(i)-(kpartype-1)*np+jlast(kpartype)
        asparse(i)=dpar(kpartype)*asparse(i)/s
 
        if(ja(i).gt.np_max) stop 'colsum out-of-bound. Increase np_max!' 

        colsum(kpartype,ja(i))=colsum(kpartype,ja(i))+abs(asparse(i))
c       if(nrow.eq.myrow)  !!! one row 
c    &    colsumy(kpartype,ja(i))=colsumy(kpartype,ja(i))+abs(asparse(i))
      end do  

      ! check stationcode for new station (uncomment is faster)
      call cfstat(stationcode,networkcode,nstat,statlist,netwlist,
     &      kstat)
      if(kstat.eq.0) then             ! add station to the list
        nstat=nstat+1
        kstat=nstat
        if(nstat.gt.500000) stop 'nstat>500000; increase dimensions'
        statlist(kstat)=stationcode
        netwlist(kstat)=networkcode
        statlat(kstat) = rlat
        statlon(kstat) = rlon
        statelev(kstat)= relev
      endif
      statsum(kstat,kdatatype)=statsum(kstat,kdatatype)+rhs(nrow)
      mstat(kstat,kdatatype)=mstat(kstat,kdatatype)+1

      ! find event and cluster index
      ! distinguish clusters for same event but different data files
      if(kluster.gt.99) stop 'kluster>99'
      localkluster=nfile*100+kluster   ! incorporate file # in kluster #
      call cfevt(ievt,nevent,localkluster,evlist,klustlist,
     &      nklust,kklust,klust1)
      if(klustercor.eq.0.and.klust1.gt.0) then
        kklust=klust1 ! if we do not discriminate among clusters
      endif  
      if(kklust.eq.0) then
        nklust=nklust+1
        kklust=nklust
        if(nklust.gt.500000) stop 'nklust>500000; increase dimensions'
        evlist(nklust)=ievt
        klustlist(nklust)=localkluster        ! 0 if klustercor=0
        if(klust1.eq.0)  klust1=kklust    ! first cluster for this event
        evlat(nklust) = slat
        evlon(nklust) = slon
        evdep(nklust) = sdep
      endif
      evtsum(kklust,kdatatype)=evtsum(kklust,kdatatype)+rhs(nrow)
      mevt(kklust,kdatatype)=mevt(kklust,kdatatype)+1

      write(444, *) kount
      write(555, *) (ja(i),i=1,kount)
      write(666, *) (asparse(i),i=1,kount)
      kountkerall = kountkerall + kount

      ! hypocentral parameter corrections
       if(jssw(4).eq.1) then
        ! add hypocentral corrections for dx(East), dy(North), 
        ! dz (negative Depth)
        if(kount+3.gt.NMAX) 
     &      stop 'with hypo added, kount>NMAX; increase dimensions'
        ilast=jlast(4)+(nevent-1)*3    ! locate columns for this event
        if(kdatatype.eq.1) then
          ! First find unit vector in ray direction:
           px=sin(raz)*sin(angle0) ! raz is src azimuth, N over E
           py=cos(raz)*sin(angle0) ! angle0 is ray take-off angle
           pz=cos(angle0)
                                ! (data) shift dT=-(px*dx+py*dy*pz*dz)/vsrc, hence:
           sv=s*vsrc            ! incorporates scaling with error s
           asparse(kount+1)=-dpar(4)*px/sv ! shift East
           asparse(kount+2)=-dpar(4)*py/sv ! shift North
           asparse(kount+3)=-dpar(11)*pz/sv ! shift Up
           
           do i=1,3
              ja(kount+i)=ilast+i ! ja keep track of column numbers
           end do         
           kount=kount+3        ! nonzeros in this row
        endif
        jlast(5)=max(jlast(5),ilast+3)
        ! XXXX make an output file for checking... 
        open(111, file='check_hypocentr_corr')
        write(111, *) kount,dpar(4),px,dpar(4),
     & py,dpar(11),pz,sv
        if(jdebug.eq.1) write(13,*) 'Evcor:',kount,px,py,pz,sv
      else
        jlast(5)=jlast(4)        ! if no corrections, else:
      endif  

      kountall=kountall+kount
      
      ! write raw datum to dumaux.x
      write(12,120) nrow,kount,rhs(nrow),s,ievt,kklust,stationcode,
     &      networkcode,kstat,iband,krtyp,klust1,kdatatype,dtheor
120   format(2i8,e14.5,e12.3,2i8,1x,a16,1x,a8,2i8,3i7,2e14.5)
!     Note that format statement 120 specifies 15 fields, but only
!     14 fields get written to dumaux. The second-to-last field, 
!     ! chiw(igroup) is added only when final aux file is written.
!    
!     K.S. 2010/08/26 format change from the below = Guusts bug fix.
!     &      networkcode,kstat,iband,krtyp,klust1,kdatatype 
! 120   format(2i8,e14.5,e12.3,2i8,1x,a16,1x,a8,2i8,3i7,e14.5)
      ! 2010/04/09 K.S. Format changed to accomodate krtyp>999
      ! (in large data sets). Note: the exact same format change
      ! must be made in mpisolvetomo.f where aux.* file gets read!
! 120   format(2i8,e14.5,e12.3,2i8,1x,a16,1x,a8,2i8,3i4,e14.5)
      if(jdebug.eq.1) write(13,*) nrow,kount,rhs(nrow),s,kstat,klust1,
     &   ' ',stationcode

      ! write one (incomplete) matrix row to scratch file dum.<ident>
      if(jdebug.eq.0) then                    ! unformatted  
        write(1) nrow,kount,iband,ddif,s,kklust,kstat
        write(1) (ja(i),i=1,kount)
        write(1) (asparse(i),i=1,kount)
      else  
        write(1,*) nrow,kount,iband,ddif,s,kklust,kstat
        write(1,*) (ja(i),i=1,kount)
        write(1,*) (asparse(i),i=1,kount)
        write(13,*) 'row ',nrow,' kount=',kount,' last asparse=',
     &        asparse(kount)
      endif  

      ! debug--remove, write out one matirx row
c     if(nrow.eq.myrow) then
c       do k=1,3
c         if(jssw(k).ne.0) then
c           open(17,file='onerow')  !!! write out one matirx row
c           cmaxmy=0.  !!! one row
c           do i=1,np
c             cmaxmy=max(abs(colsumy(k,i)),cmaxmy) !!! one row
c           enddo
c           write(17,*) '3    ',chcor(k),' scale factor:',cmaxmy
c           write(17,*) np
c           do i=1,np
c             amy=0.01*max(log10(colsumy(k,i)/cmaxmy),-3.0) 
c             write(17,346) (points(j,i),j=1,3),amy  
c           enddo
c           close(17)   !!! one row
c         endif
c       enddo
c     endif

      goto 20           ! read next row of matrixfile

200   continue          ! EOF reached
      if(jdebug.gt.0) write(13,*) 'EOF reached'

      write(iu4,128) 'MATXF',nrow-lastnrow,kdatatype,fname
 128  format(a8,i12,i3,2x,a)
      lastnrow = nrow
      ! end Case 5 (read new matrixA/matrixT file)      

      goto 10   ! get next input file name or processing instruction

c...  END OF INPUT
300   continue

      if(jdebug.gt.0) write(13,*) 'End of reading input files'
      write(12,120) 0,0,0.,0.,0,0,'End-of-matrix','xx',0,0,0,0,0

      if(ngroup.eq.0) then
        ngroup=1
        nchi(ngroup)=nrow
        nallgr=nallgr+nchi(ngroup)
      endif  

      if(nallgr.ne.nrow) then
        print *,'Error: not all data in a group.'
        print *,'Cure: end input with a group statement!'
        stop 'Error in data grouping'
      endif  

      ! save scratch files and start reading them back
      close(1)
      close(12)
      if(jdebug.eq.0) then
        open(1,file='dum.'//ident,form='unformatted')
      else
        open(1,file='dum.'//ident)
      endif
      open(12,file='dumaux.'//ident)

      ! Now add the rest of the corrections in blocks
      ! we know already jlast(5) for dtP station corrections
      ! Note: jlast(i+1) -- and not jlast(i) gives the last column 
      ! of the ith correction 
      !
      ! 5: dt_P station corrections
      jlast(6)=jlast(5)         
      if(jssw(5).eq.1) jlast(6)=jlast(5)+nstat
      ! 6: dlnA_P station corrections
      jlast(7)=jlast(6)         
      if(jssw(6).eq.1) jlast(7)=jlast(6)+nstat
      ! 7: dTo event correction (origin time)
      jlast(8)=jlast(7)         
      if(jssw(7).eq.1) jlast(8)=jlast(8)+nevent
      if(jssw(7).eq.2) jlast(8)=jlast(8)+nklust 
      ! 8: dAo event correction (scalar moment Mo)
      jlast(9)=jlast(8)
      if(jssw(8).ge.1) jlast(9)=jlast(9)+nklust      
      ! 9: dt_S station correction
      jlast(10)=jlast(9)        
      if(jssw(9).eq.1) jlast(10)=jlast(9)+nstat
      ! 10: dlnA_S station correction
      jlast(11)=jlast(10)     
      if(jssw(10).eq.1) jlast(11)=jlast(11)+nstat
      ncol=jlast(11)
      if(ncol.gt.NDIM) stop 'ncol>NDIM. Increase NDIM.'

      if(jdebug.eq.0) then
        open(3,file='mat.'//ident,form='unformatted')
      else  
        write(13,*) 'jlast=',jlast
        write(13,*) 'ncol=',ncol
        open(3,file='mat.'//ident)
      endif  
      open(2,file='aux.'//ident)

      ! Output diagnostics to assemblematrix.out.<ident>
      write(iu4,fmt='(/,a)') 'DATA_GROUPS_CHI2_WEIGHTED'
      if(ngroup.le.1) then
        write(iu4,*) 0
        write(iu4,*) 'No data grouping used.'
      else
        write(iu4,*) ngroup
        write(iu4,*) ' grp   nrows    weight'
        do i=1,ngroup
          write(iu4,fmt='(i3,i8,f10.4)') i,nchi(i),
     &          float(nrow)/float(nchi(i)*ngroup) ! called chiw(i) later on
        enddo
      endif  

      ! find data average for all data 
      sum=0
      rms=0
      do i=1,nrow
        sum=sum+rhs(i)
        rms=rms+rhs(i)**2
      enddo
      average=sum/nrow
      rms=sqrt(rms/nrow)

      write(iu4,fmt='(/,a)') 'MATRIX_DIMENSIONS'
      write(iu4,131) nrow,' NROW number of matrix rows'
      write(iu4,131) ncol,' NCOL number of matrix columns'
 131  format(i12,1x,a)

      write(iu4,fmt='(/,e10.4,a)') average, '  Avg_of_ALL_data (univar)'
      write(iu4,fmt='(  e10.4,a)') rms,     '  RMS_of_ALL_data (univar)'

      write(iu4,fmt='(/,a)') 'SOLUTION_VECTOR_COMPONENTS'
      write(iu4,fmt='(i4,1x,a)') 10, ' blocks in model vector'
      write(iu4,fmt='("jssw      ncol    1st_col    last_col",4x,
     &"dpar  parameter_name")')
      do i=1,11
        if(jlast(i+1).gt.jlast(i)) then
          ii = jlast(i+1)-jlast(i)           
          write(iu4,330) jssw(i),ii,jlast(i)+1,jlast(i+1),dpar(i),
     &         chcor(i)
        else
          write(iu4,330) jssw(i),0,0,0,0.0,chcor(i)
        endif  
      enddo
330   format(i5,i10,i10,i10,3x,e10.5,2x,a32)      

      write(iu4,fmt='(/,a)') 'TRAVEL_TIME_CORRECTION_SWITCHES'
      write(iu4,fmt='("kssw      switch_name")')
      do i=1,4
        write(iu4,fmt='(i4,6x,a)') kssw(i),chcorrection(i)
      enddo
      if (kssw(4).eq.2)
     &  write(iu4,fmt='(a,f4.1,a)') 'user defined reference period = ',
     &  refperiod, ' sec'

      write(iu4,fmt='(/,a)') 'SOURCES_AND_RECEIVERS'
      write(iu4,131) nstat, ' NSTAT  stations (receivers)'
      write(iu4,131) nevent,' NEVENT earthquakes'
      write(iu4,131) nklust,' NKLUST data clusters (>=nevent)'
      write(iu4,131) np,    ' NP     points in velocity model'

      ! output station statistics
      write(6,fmt='(//,a)') 'Writing event and station logs...'      
      open(20,file='assemblematrix.stations.'//ident)
      write(20,fmt='(i8,2x,a)') nstat,' stations. Avg before demean.'
      write(20,fmt='("Station",10x,"Netw",7x,"N(dt)",4x,"N(A)",
     &   4x,"dT(av)",2x,"dlnA(av)",4x,"stlat",4x,"stlon",4x,"stelev")')
      write(7,fmt='(//,"Station table entries")')
      do i=1,nstat
        avdt=statsum(i,1)/max(1,mstat(i,1))
        avda=statsum(i,2)/max(1,mstat(i,2))
        write(20,340) statlist(i),netwlist(i),mstat(i,1),mstat(i,2),
     &        avdt,avda,statlat(i),statlon(i),statelev(i)
        write(7,341) statlist(i),netwlist(i),mstat(i,1),mstat(i,2),
     &        avdt,avda
      enddo
340   format(a16,1x,a8,2i8,2f10.3,3f10.3)
341   format(a16,' & ',a8,' & ',i8,' & ',i8,' & ',f10.3,' & ',f10.3,
     &  ' \\\\')
      close(20)

      ! output event statistics
      open(21,file='assemblematrix.events.'//ident)
      write(21,fmt='(i8,2x,a)') nklust, ' nklust (data cluster lines)'
      write(21,fmt='(3x," ievt  klust",3x,"N(dt)",4x,"N(A)",
     &      4x,"dT(av)",2x,"dlnA(av)",4x,"evlat",4x,"evlon",
     &      4x,"evdep")')
      write(7,fmt='(//,"Event table entries")')
      do i=1,nklust
        avdt= evtsum(i,1)/max(1,mevt(i,1))
        avda= evtsum(i,2)/max(1,mevt(i,2))
        write(21,342) evlist(i),klustlist(i),mevt(i,1),mevt(i,2),
     &        avdt,avda,evlat(i),evlon(i),evdep(i)
        write(7,343) evlist(i),klustlist(i),mevt(i,1),mevt(i,2),
     &        avdt,avda
      enddo
342   format(i8,1x,i6,2i8,2f10.3,3f10.3)
343   format(i8,' & ',i3,' & ',i8,' & ',i8,' & ',f10.3,' & ',f10.3,
     &  ' \\\\')
      close(21)
      close(7)

      ! write column density to gmt.xyz file
      write(6,fmt='(/,a)') 'Writing inmapc, columndensity...'      
      open(9,file='runmapc')
      write(9,fmt='(a)') '#! /bin/csh'
      do k=1,3
        if(jssw(k).ne.0) then
          write(fname2,344) k,ident
          open(7,file=fname2)
344       format('columndensity.',i1,'.',a)
          write(fname3,345) k
          open(8,file=fname3)
345       format('inmapc.',i1)
          ! normalize column density
          cmax=0.
          amx=-3.0
          do i=1,np
            cmax=max(abs(colsum(k,i)),cmax)
          enddo  
          write(7,*) '3    ',chcor(k),' scale factor:',cmax
          write(7,*) np
          do i=1,np
            a=0.01*max(log10(colsum(k,i)/cmax),-3.0)    ! log % of max
            amx=max(a,amx)
            write(7,346) (points(j,i),j=1,3),a,colsum(k,i)
346         format(3f8.1,e14.5,e14.5)
          enddo
          close(7)
          ! write input file for map.f
          write(8,*) 0
          write(8,fmt='(a)') fname2
          write(8,fmt='(a)') nodesfile
          write(8,*) -90,90,1
          write(8,*) 0,360,1
          write(8,fmt='(a,i1)') 'coverage',k
          write(8,fmt='(a)') '300','600','1000','1400','1800',
     &          '2200','2600','-1'
          close(8)
          write(9,fmt='("map < ",a)') fname3
          write(9,*) 'set amx =',amx*100.
          write(9,*) 'rm -f dens.cpt'
          write(9,*) 'makecpt -Chot -I -T-3.0/$amx/0.1 -Z > dens.cpt'
          write(9,fmt='("gmtmapc coverage",i1,".",a)') k,'300',k,'600',
     &          k,'1000',k,'1400',k,'1800',k,'2200',k,'2600'
        endif
      enddo  

      kountall=0
      if(jdebug.eq.0) then
        write(3) nrow,ncol
        write(3) jlast
        write(3) jssw
        write(3) dpar
      else
        write(3,*) nrow,ncol
        write(3,*) jlast
        write(3,*) jssw
        write(3,*) dpar
      endif  

      igroup=1
      lastgroup=nchi(1)
      if(chiw(1).eq.0.) chiw(1)=float(nrow)/float(nchi(1)*ngroup)

      if(jdebug.gt.0) 
     &  write(13,*) 'Now re-reading matrix elements, adding corrections'

      ! read next datum from scratch auxiliary file dumaux.<ident>
      ! change the rhs to demeaned value and write to auxiliary file
400   read(12,120) krow0,kount0,rhskrow,s,ievt,kklust0,stationcode,
     &      networkcode,kstat0,iband,krtyp,klust1,kdatatype,dtheor
!     &      networkcode,kstat0,iband,krtyp,klust1,kdatatype ! K.S. before 2010/08/26
      rhskrow0=0.
      if(krow0.gt.lastgroup) then
        igroup=igroup+1
        lastgroup=lastgroup+nchi(igroup)
        if(igroup.gt.ngroup) stop 'igroup error???'   ! debug statement
        if(chiw(igroup).eq.0.) 
     &        chiw(igroup)=float(nrow)/float(nchi(igroup)*ngroup)
      endif
      if(krow0.gt.0) rhskrow0=rhs(krow0)*chiw(igroup)
      if(jdebug.gt.0) write(13,*) krow0,kount0,rhskrow,s,kklust0,
     &      kstat0,' ',stationcode
      if(krow0.eq.0) goto 800   ! last data line

      ! read next matrix row from scratch file
      if(jdebug.eq.0) then                    ! unformatted  
        read(1,end=900) krow,kount,iband,ddif,s,kklust,kstat
        read(1) (ja(i),i=1,kount)
        read(1) (asparse(i),i=1,kount)
      else  
        read(1,*,end=900) krow,kount,iband,ddif,s,kklust,kstat
        write(13,*) krow,kount,iband,ddif,s,kklust,kstat
        read(1,*) (ja(i),i=1,kount)
        read(1,*) (asparse(i),i=1,kount)
        write(13,*) 'krow=',krow,' last asparse=',asparse(kount)
      endif  

      ! write raw datum and scaled, demeaned datum to data.<ident>
      write(10,*) ddif,rhs(krow0)*chiw(igroup),kdatatype

      ! debug output, remove after testing
      if(krow.ne.krow0.or.kount.ne.kount0.or.kklust.ne.kklust0.
     &     or.kstat.ne.kstat0) then
        print *,'BUG! krow,kount,kklust,kstat=',krow,kount,kklust,kstat
        print *,'but in aux file=',krow0,kount0,kklust0,kstat0
        stop
      endif 

      ! Station P time correction (positive for late station delay)?
      if(jssw(5).eq.1.and.kdatatype.eq.1.and.krtyp.eq.1) then
        kount=kount+1
        ja(kount)=jlast(5)+kstat
        asparse(kount)=dpar(5)/s
        if(jdebug.gt.0) write(13,*) chcor(5),asparse(kount)
      endif

      ! Station P amplitude correction?
      if(jssw(6).eq.1.and.kdatatype.eq.2.and.krtyp.eq.1) then
        kount=kount+1
        ja(kount)=jlast(6)+kstat
        asparse(kount)=dpar(6)/s
        if(jdebug.gt.0) write(13,*) chcor(6),asparse(kount)
      endif

      ! Source origin time correction?  !!! kdatatype is always 1!?
      if(jssw(7).gt.0.and.kdatatype.eq.1) then
        kount=kount+1
        kolumn=kklust                   ! if To changes per cluster
        print *, kklust
        if(jssw(1).eq.1) kolumn=klust1  ! if To changes only per event
        ja(kount)=jlast(7)+kolumn
        asparse(kount)=dpar(7)/s
         if(jdebug.gt.0) write(13,*) chcor(7),asparse(kount)
      endif

      ! Source amplitude correction?
      if(jssw(8).gt.0.and.kdatatype.eq.2) then
        kount=kount+1
        ja(kount)=jlast(8)+kklust  
        asparse(kount)=dpar(8)/s
        if(jdebug.gt.0) write(13,*) chcor(8),asparse(kount)
      endif

      ! Station S time correction?
      if(jssw(9).ne.0.and.kdatatype.eq.1.and.krtyp.eq.2) then
        if(jssw(9).eq.1) then
          kount=kount+1
          ja(kount)=jlast(9)+kstat
          asparse(kount)=dpar(9)/s
          if(jdebug.gt.0) write(13,*) chcor(9),asparse(kount)
        else    ! coupling with dtP
          kount=kount+1
          ja(kount)=jlast(5)+kstat
          asparse(kount)=dpar(5)*sqrt3/s
          if(jdebug.gt.0) write(13,*) chcor(9),asparse(kount),' *sqrt3'
        endif  
      endif

      ! Station S amplitude correction?
      if(jssw(10).ne.0.and.kdatatype.eq.2.and.krtyp.eq.2) then
        kount=kount+1
        ja(kount)=jlast(10)+kstat
        asparse(kount)=dpar(10)/s
        if(jdebug.gt.0) write(13,*) chcor(10),asparse(kount)
      endif

      ! K.S. 2010/08/05: stop execution if matrix structure
      ! looks extremely bad, do not ask user
      ! (in practice this seems to 
      ! indicate serious memory problems -- happened with gfortran
      ! compilation, but not with ifort) 
      do i=1,kount
        j=ja(i)
        asparse(i)=asparse(i)*chiw(igroup)
        ata(j)=ata(j)+asparse(i)**2
        if(ata(j).gt.1.0e30) then
          print *,'Large or infinite ata in column',j
          print *,'for row',krow
          print *,'i,asparse,chiw=',i,asparse(i),chiw(igroup)
          print *,kount,iband,ddif,s,kklust,kstat
!          print *,'Type 0 to continue, 1 to stop:'
!          read *,kstop
!          if(kstop.ne.0)then
            !!!!!!! 
            print *,'STOPPED EXECUTION, bad conditioning.'
            print *,'Compiler problem?? Try ifort -mcmodel=large.'
            print *,'K.S. 2010/08/05 experimental.'
            stop
!          endif
        endif  
      enddo  

      kountall=kountall+kount

      if(jdebug.gt.0) write(13,*) 'row',krow,' kount=',kount,
     &  ' kountall=',kountall

      ! write complete matrix row to matrix file mat.<ident>
      if(jdebug.eq.0) then                    ! unformatted  
        write(3) krow,kount,iband,ddif,s,kklust,kstat,chiw(igroup)
        write(3) (ja(i),i=1,kount)
        write(3) (asparse(i),i=1,kount)
      else  
        ! change to free format (3,*) after initial file debugging
        write(3,*) krow,kount,iband,ddif,s,kklust,kstat,chiw(igroup)
        write(3,fmt='(10i7)') (ja(i),i=1,kount)
        write(3,fmt='(10f10.6)') (asparse(i),i=1,kount)
      endif  

      ! write aux.<ident> file with updated value of kount
      write(2,120) krow0,kount,rhskrow0,s,ievt,kklust0,stationcode,
     &      networkcode,kstat0,iband,krtyp,klust1,igroup,
     &      chiw(igroup),dtheor
!    K.S. was until 2010/08/26:
!     &      networkcode,kstat0,iband,krtyp,klust1,igroup,chiw(igroup)

      goto 400

800   continue   ! done evaluating data lines for corrections

      ! aux.<ident>: append station list and event list
      ! aux.<ident> and data.<ident>: finish off with group weights
      write(2,120) 0,0,0.,0.,0,0,'End-of-matrix','xx',0,0,0,0,0,0.,0.
!      write(2,120) 0,0,0.,0.,0,0,'End-of-matrix','xx',0,0,0,0,0,0. ! until 2010/08/26
      write(2,310) nstat,(i,statlist(i),netwlist(i),i=1,nstat)
310   format(i5,/,(i5,1x,a16,1x,a8))
      write(2,320) nklust,(i,evlist(i),klustlist(i),i=1,nklust)
320   format(i5,/,(i5,2i10))

      write(10,*) 'Data groups and weights:'
      do i=1,ngroup
        write(10,*) i,nchi(i),chiw(i)
        write(2,*) i,nchi(i),chiw(i)
      enddo  
      write(2,*) kssw ! inserted 2010/08/26 K.S.
      close(10)
      close(2)
      close(3)

      ! Finish log file assemblematrix.out.<ident>
      write(iu4,fmt='(/,a)') 'DATA_MASKING'     
      write(iu4,fmt='(i10,1x,a)') nrow,' NROW accepted data(rows)' 
      write(iu4,fmt='(i10,1x,a)') nrej,' NREJ rejected data(rows)' 
      if(.not.lmask)then
        write(iu4,fmt='(a,1x,a)') 'NO_SUCH_MASKFILE ',maskfile         
      else  ! info on masking criteria
        write(iu4,fmt='(a12,1x,f7.3)') 'XCORR_MIN',xcorr_min
        do kk=1,mb0
          write(iu4,fmt='(a14,1x,i3)') 'IGNORE_BAND',m_bands(kk)
        enddo
        do kk=1,mbt
          write(iu4,fmt='(a14,1x,i3)') 'IGNORE_TT_BAND',m_ttbands(kk)
        enddo
        do kk=1,mba
          write(iu4,fmt='(a14,1x,i3)') 'IGNORE_AMP_BAND',m_ampbands(kk)
        enddo
      endif  
      print *,'Total number of nonzero matrix elements: ',kountall
      write(iu4,fmt='(/,a)') 'MATRIX_ELEMENTS' 
      write(iu4,*) kountall,' NONZERO_ELEMENTS'

      atamax=0.
      do i=1,ncol
        atamax=max(ata(i),atamax)
      enddo
      write(iu4,*) atamax," ATA_MAX largest element of A'A"
      write(iu4,*) sqrt(atamax)," sqrt(ATA_MAX)"

      write(333, *) nrow 
      write(333, *) np 
      write(333, *) kountkerall

      print *
      print *,"Largest element of A'A=",atamax,", sqrt=",sqrt(atamax)
      print *,'Number of accepted data/rows     : ',nrow  
      print *,'Number of rejected (masked) data : ',nrej
      tmp1 = 100.*float(nrow)/float(nrow+nrej)
      tmp2 = 100.*float(nrej)/float(nrow+nrej)
      write(6,fmt='(f5.1,a)') ,tmp1,'% accepted'
      write(6,fmt='(f5.1,a)') ,tmp2,'% rejected'
      print *

      ! Finish file inmapc.{1,2,3}
      k=0
      do kpartype=1,3
        if(jssw(kpartype).ne.0) then
          write(fname2,844) kpartype,ident
          open(8,file=fname2)
          write(8,*) "3   A'A diag ",chcor(kpartype),atamax,sqrt(atamax)
          write(8,*) np
          do i=1,np
            write(8,346) (points(j,i),j=1,3),ata(k+i)/atamax
          enddo
          close(8)
          k=k+np
        endif
      enddo
844   format('ata.',i1,'.',a)
845   format('atamx.',a)
      write(fname2,845) ident
      open(8,file=fname2)
      write(8,*) atamax,sqrt(atamax)
      close(8)

      !-----------------------------------------
      stop 'Normal end of assemblematrix'
      !-----------------------------------------

900   print *,'Error: scratch file has fewer data than aux file'
      print *,'This occurred at krow=',krow

      end


c========================================================================
      
      subroutine cfstat(stationcode,networkcode,nstat,statlist,
     &      netwlist,kstat)

      ! finds station in existing list of stations or kstat=0
      ! input: stationcode and networkcode=network code for current station
      !        statlist(i) and netwlist(i), i=1,nstat: list of stations
      ! output: kstat=station identifier in the list

      character*16 stationcode,statlist(500000)
      character*8 networkcode,netwlist(500000)
      kstat=0
      do i=1,nstat
        if(statlist(i).eq.stationcode.and.netwlist(i).eq.
     &        networkcode) then
          kstat=i
          return
        endif
      enddo
      return
      end

      subroutine cfevt(ievt,nevent,localkluster,evlist,klustlist,
     &        nklust,kklust,klust1)

      ! finds event and cluster in existing list of events and clusters
      ! input: ievt, localkluster is current event number & cluster #
      !        evlist, klustlist(i), i=1,nklust: list of events
      !        (ievt conforms to event number in file <evlist>)
      ! output: kklust=location of ievt in the list, or 0 if absent.
      !        (every cluster is seen as separate event!)
      !        klust1=first cluster for this event

      ! Note on cluster numbering:
      ! In the data/matrix files, every event has one or more data clusters; 
      ! such data are labeled with variable <kluster>. Because different
      ! data files may contain the same event, we avoid overlap by redefining
      ! <kluster> from file nfile here as <localkluster>=100*nfile+kluster.
      ! The localklusters are stored in klustlist, which has the same
      ! ordering as the list of event numbers in evlist.
      ! Every source-station path then belongs to an event/cluster combination
      ! which is numbered by <kklust>. The first kklust for each ievt is
      ! also stored, to use when we do not discriminate between clusters of the
      ! same event: <klust1>.

      integer evlist(500000),klustlist(500000)
      kklust=0
      klust1=0
      neweve=1
      do i=1,nklust
        if(evlist(i).eq.ievt.and.klust1.eq.0) then
          klust1=i      ! first cluster for this hypocenter
        endif
        if(evlist(i).eq.ievt.and.klustlist(i).eq.localkluster) then
          kklust=i      ! ranking among all clusters for all hypocenters
          return
        endif
        if(evlist(i).eq.ievt) neweve=0
      enddo
      if(neweve.eq.1) nevent=nevent+1   ! new event (hypocenter)
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
        include	'includes/nn.param'
        include 'includes/setdel2.h'

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
        common/nodes/nodesfile

      print *,'Give vertex file name (eg vertices.xyz):'
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
      if(debug)write(*,*)' done nn3d_setup'
 
      if(debug)write(*,fmt='(/"  Number of points = ",i7,
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
	if(debug)write(*,*)' Number of faces on the hull =',nf
	if(debug)write(*,*)' Number of nodes on the hull =',nh
	if(debug)write(*,*)' Counted number of nodes on the hull ='
     &           ,nhc
 	if(writeout)then
	   write(*,*)' vertices of triangles and normals on hull '
	   write(*,*)' (indices are local to hull)'
        end if
        call build_hullneighbour
     &       (nh,hullfaces,nf,nwork2d,
     &        hullneighbour,nhwork1,nhwork2,nhwork3)
	if(debug)then
           write(*,*)' checking hullneighbour array'
           call check_neighbour
     &          (hullneighbour,2,hullfaces,nf,consistent)
           write(*,*)' done hullneighbour array'
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
      write(*,*)' checking neighbour array'
      call check_neighbour(neighbour,3,vertices,nt,consistent)
      write(*,*)' done neighbour array'
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
         if(out)write(*,*)' center of tet',i,' is outside of 
     &the convex hull'
         if(i.ne.loc)then
             write(*,*)' center of',i,' in tetrahedra',loc
             write(*,*)' x:',(xd(j),j=1,3)
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
 
        include 'includes/nn.param'
        include 'includes/setdel2.h'

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

 
	if(loc.eq.0) loc = 1

        call tetloc(xd,points,vertices,neighbour,loc,out)
	if(out)then
           if(extendtetra)then
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
           else
             write(*,*)'points outside hull'
             write(*,*)'xd(1),xd(2),xd(3)',xd(1),xd(2),xd(3)
	     loc = 0
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
           a1 = points(1,v3)-points(1,v2) 
           a2 = points(2,v3)-points(2,v2) 
           a3 = points(3,v3)-points(3,v2) 
           b1 = points(1,v4)-points(1,v2) 
           b2 = points(2,v4)-points(2,v2) 
           b3 = points(3,v4)-points(3,v2) 
           perp(1) = a2*b3-b2*a3
           perp(2) = a3*b1-b3*a1
           perp(3) = a1*b2-b1*a2
           del1 = 0.d0
           do 11 j=1,3
              del1 = del1 + perp(j)*(points(j,v1)-points(j,v2))
  11       continue 
           del2 = 0.d0
           do 12 j=1,3
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
 
        include 'includes/nn.param'
        include 'includes/setdel2.h'

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

      subroutine inpgrd

      include 'includes/nn.param'        ! parameter statements
      include 'includes/setdel3.h'       ! common for np,nt
      include 'includes/setdel2.h'       ! common for points, centres, data etc

      idebug=0
      iwrite=0
      iextend_outside_hull=0
      icall_find_node=1

      call nn_init(idebug,iwrite,icall_find_node,
     *iextend_outside_hull,nt,np)

      print *,'Read model parameterization, nt,np=',nt,np

      return
      end

!     ---------------------------------------------------------------


!     ---------------------------------------------------------------

      subroutine read_maskfile(fnam,lmask,m_bands,mb0,xcorr_min,
     & m_ttbands,mbt,m_ampbands,mba,nbx)
      implicit none
 
      ! Revision log: created 2006/11/17 by Karin Sigloch
      !   Uses function len_trim, which is nonstandard Fortran
      !   but recognized g77 and absoft compilers.
      ! 2007/08/10: no longer writes output file assemblematrix.mask.<ident>
      !   Instead this is now the INPUT file (given in fnam). 
      ! 2007/08/28: added option to selectively discard only amplitudes      
      !   or ttimes in a certain band (IGNORE_AMP_BAND, IGNORE_TT_BAND)  
      !
      ! input:
      character*256 fnam 
      integer nbx               ! max number of frequency bands

      ! output
      ! (second component of 2-vectors contain default value)
      logical lmask      ! whether or not to apply any masking
      integer mb0,mbt,mba
      integer m_bands(nbx)    ! bands in which reject both amp and ttime data
      integer m_ttbands(nbx)  ! bands in which to reject ttime data
      integer m_ampbands(nbx) ! bands in which to reject amp data
      real xcorr_min     ! min. acceptable xcorrelation coefficient
      real snr_min       ! min. acceptable snr
      
      ! local
      integer iu1,iu2,ios,ns,k
      real val
      character*256  line,fnam2
      character*80  nam

      ! -------------------
      ! set default/reference values

      mb0        = 0
      mbt       = 0
      mba       = 0
      xcorr_min = -1.
      snr_min   =  0.
      
      ! --------------------
      
      !   Start reading
      ns = len_trim(fnam)       
      fnam2 = fnam(1:ns)      
      inquire(file=fnam2,exist=lmask)
      if(.not.lmask)then
         write(6    ,fmt='(a,a)') 'NO_SUCH_MASK_FILE ',fnam2(1:ns)
         return
      else   
         write(6,fmt='(a,a)')     'MASK_FILE_CONTENTS ',fnam2(1:ns)
      endif         
      iu1=23      
      open(iu1,file=fnam,status='old',err=999)

10    continue   ! while there are more input lines...
      read(iu1,fmt='(a256)',end=100,err=997) line

       ! skip blank lines and comments
      ns = len_trim(line)
      if(ns.eq.0.or.index(line,'#').eq.1) goto 10

      ! Format is assumed to be two-column: character, numeric_value
      backspace(iu1)
      read(iu1,*,err=996) nam,val
      if(nam=='IGNORE_BAND')then    ! 1) frequency band masking
         mb0 = mb0+1
         m_bands(mb0)= int(val)
      elseif(nam=='XCORR_MIN')then  ! 2) xcorrelation threshold
         xcorr_min = val
      elseif(nam=='IGNORE_AMP_BAND')then
         mba = mba+1
         m_ampbands(mba)= int(val)
      elseif(nam=='IGNORE_TT_BAND')then
         mbt = mbt+1
         m_ttbands(mbt) = int(val)
      ! ... Add more options here...   
      else
         goto 980                   ! default: format error
      endif

      ! copy (valid) line to logfile
      write(6    ,fmt='(a)') line(1:len_trim(line))
      
      goto 10   ! ...read next line of input file

! ----------------

 100  continue
      close(iu1)
      print *
      return

! ------  catch errors  ----------
999   continue
      print *, 'ERROR: could not open maskfile!'
      print *, line
      stop
997   print *, 'ERROR: while reading this line from maskfile:'
      print *, line
      stop
996   print *, 'ERROR while reading char/real pair from maskfile:'
      print *, line
      stop
 980  print *, 'ERROR in read_globals: check format of the'
      print *, 'following line of input file ',len_trim(fnam),': '
      print *, len_trim(line)
      print *, 'Note that comment lines have to start with an #,'
      print *, 'and blank lines must contain blanks only.'
      stop

      end  ! subroutine 

      
