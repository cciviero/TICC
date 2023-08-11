      program mpisolvetomo

      ! matrix rows  for at most MMAX nonzero elements per row, 
      ! NDIM unknowns (including corrections) and NMAX data

      ! compile: mpif77 -o mpisolvetomo mpisolvetomo.f
      ! or with flags: -g -fbounds-check -fno-automatic -w 
      ! (-w suppresses warnings for differing data types in MPI calls)

      ! Change log:
      ! 2010/04/09 K.S. changed format 120, 3i4 --> 3i7 to prevent overflow.
      !   Needed to make the same change in assemblematrix.f and solvetomo.f
      !   but not resolutiontest.f, since free format read there.

      ! a dummy line, test submit

c Files used and their unit numbers:
c     1 vertices file  IN
c     2 aux.xxxx (auxiliary or data file)  IN
c     3 mat.xxxx (matrix file)  IN
c     4 out.solve.xxxx OUT
c     7 tex.solve.xxxx OUT
c     8 atamx.xxxx  IN
c     8 xy.tradeoff.xxxx (GMT file for Lcurve) OUT
c     9 facets file  IN
c     9 solx.xxxx (solution) OUT
c     9 sol.xxxx (full solution, raw chis) OUT
c     10 res.xxxx (residuals file) OUT
c     13 diagnostics.solve OUT
c     21 outl.xxx outliers OUT
      
      include "includes/solvet.h"
ccc      include "includes/mpif.h"  ! Guusts Princeton implementation, 
c                                      ! worked only with mpich
      include "mpif.h"                 ! choose the locally appropriate mpi
      integer status(MPI_STATUS_SIZE)

c solvet.h has default STILL NEEDS TESTING:
c     parameter(MMAX=200000) ! max nonzero in row of a (non-MPI)
c     parameter(MPIX=50000000)  ! max dimension of MPI array segments
c     parameter(NMAX=1000000)! max data vector (incuding damping 0's)
c     parameter(NDIM=200000) ! max model dimension (including corrections)
c     parameter(NXTRY=100)   ! max trials for chi2 convergence
c     parameter(NBR=100)     ! max number of nearest neighbours
c     parameter(NFN=128)     ! max length of file name
c     parameter(NOUT=200000) ! max nr of outliers
c     parameter(MPROC=32)    ! max nr of processors for MPI
c     parameter(NPMAX=100000) ! number of model grid nodes


      dimension rhs(NMAX),sigma(NMAX),ievt(NMAX),v(NDIM)
      REAL, ALLOCATABLE :: colden(:)
      dimension w(NDIM),w1(NDIM),w2(NDIM),w3(NMAX)        ! scratch
      dimension iband(NMAX),kklust(NMAX),u(NMAX),sol(NDIM),kstat(NMAX)
      dimension a(MPIX),ja(MPIX)        ! array storage
      dimension nb(NPMAX),iv(NBR,NPMAX),tvol(NPMAX),tnab(NPMAX)

      ! except during transit from file, ja and a are not used in root
      ! and equivalence statements can save memory for variables not 
      ! used in subprocesses
      equivalence (iv,ja(1))       
      equivalence (w,a(1))
      equivalence (w2,a(NDIM+1))
      equivalence (sol,a(2*NDIM+1))
      equivalence (w3,a(3*NDIM+1))

      dimension jlast(11)
      dimension nv(4),points(3,NPMAX)
      dimension jssw(10),nssw(3),kssw(3),kountpar(3),dpar(10)
      integer evlist(500000),klustlist(500000),mevt(500000,2)
      integer mstat(500000,2)
      INTEGER :: colden_writer, colden_do, nrden
      REAL :: tmptmp, colvalue, xcol, mean_colden
      dimension statsum(500000,2),evtsum(500000,2)
      dimension ytry(NXTRY),chi2try(NXTRY),chi2grp(100),grpw(NMAX)
      integer ngrp(100),ia(NMAX)
      integer knt(0:MPROC),kdisp(0:MPROC)
      integer (kind=8) lasttot, maxperproc, nrow_kntaverage

      ! various character variables. Make sure first equals NFN
      character*128 fname,ident,ident2,fname2,fname3,directory,fname4

      character*30 dataf,dataf1,dataf2
      character*19 chcor(10)
      character*16 stationcodekrow,statlist(500000)
      character*8 phase,networkcode,netwlist(500000)
      character*5 chcor5
      character*3 comp

      data chcor/'dlnVp','dlnVs','dlnQs','Hypoctr corr (km)',
     &'Station corr tP (s)', 'Station corr AP (s)','Origin t corr (s)',
     &'Event corr dlnA-P', 'Station corr tS (s)','Station corr dlnA-S'/

      ! outlier treatment
      logical outl(NMAX)
      data    outl/NMAX*.false./
      real    outlim
      integer klean

c declare: noutl_grp(igrp),chi2_grp(igrp),chi2i_all,chi2g
      dimension noutl_grp(100),chi2_grp(100)
      real      chi2i_all,chi2g
      
      data bigint/2.1e9/        ! largest integer in 32 bit is 2^31

      dimension ytry_lcurve(11)
      !data ytry_lcurve/1.0, 0.01,.05,.07,.1,.15,.2, .3, .4, .6, .8/

cccc   TEST
c      data ytry_lcurve/1.,.0001,.001,.01,.1,.3,.6,.9,.95,.99,.995/
c      data ytry_lcurve/1.,.0001,.001,.003,.005,.007,.01,.02,.03,.05,.1/
      data ytry_lcurve/1.,0.,.0001,.001,.003,.005,.007,.01,.025,.05,.1/
    !   data ytry_lcurve/1.,0.,.0001,.0015,.002,.0025,
    !  &.003,.004,.005,.007,.1/

cccc   TEST
      ! ytry_lcurve must have exactly 11 elements and the first two are
      ! not actually read (the program automatically sets them 
      ! to 1.0 and 0.03, respectively)
      ! The values that used to be hardcoded are:
!      data ytry_lcurve/1.0, 0.03, .1, .2, .3, .4, .5, .6, .7, .8, .9/


c end of declarations


      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr) ! get # of procs
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)  ! get process number
      tcpu0=secnds(0.0)

      klean=0
      ibr=0             ! bailout flag
c end of mpi initialization
c=======================================================================
      if(myid.eq.0) then                ! IF ROOT PROCESS
c=======================================================================

      write(6,fmt='(a)') 'Running mpisolvetomo...' 
      write(6,*) 'nproc=',nproc
         
      itrymax=25        ! max trials for chi2 convergence
      if(itrymax.gt.NXTRY) call bailout(ibr,'increase parameter NXTRY.')

c all user input is here:

c  GUIDE TO USER INPUT PARAMETERS ksmooth/chi2target/epsnorm/epssmooth/epsratio

c  ksmooth = 1 : smoothing matrix S = Identity
c  ksmooth = 1 & chi2target < 0  (single run, abs values of epsnorm and epssmooth
c                as read from input)
c  ksmooth = 1 & chi2target >= 0 (L-curve only (0) or value of chi2>0 specified)
c                input parameter epsratio determines
c                relative weighting of damping and smoothing; input epsnorm and
c                and epssmooth are ignored and internally recomputed as
c                eps = ytry/(1-ytry)  (0<=ytry<=1, i.e. 0<eps<inf)
c                epssmooth=eps*sqratamax     ! scale w.r.t. matrix elements
c                epsnorm=epsratio*epssmooth 
c 
c  ksmooth = 2 : smoothing matrix S = ??
c                epssmooth=1.0 (input value is ignored)
c  ksmooth = 2 & chi2target >= 0 (L-curve only (0) or value of chi2>0 specified)
c                epssmooth is set to 1.0; epsnorm is scaled internally:
c                eps = ytry/(1-ytry)  (0<=ytry<=1, i.e. 0<eps<inf)
c                epsnorm=eps*sqratamax
c  ksmooth = 2 & chi2target < 0 :
c                epsnorm ???? (actual input value, but epssmooth forced to 1.???)
c             
c
c  ksmooth = 3 : smoothing matrix S = ??
c                input epssmooth is actually used;
c  ksmooth = 3 & chi2target >= 0 (L-curve only (0) or value of chi2>0 specified)
c                epsnorm is scaled internally as:
c                epsnorm=eps*sqratamax
c  ksmooth = 3 & chi2target < 0
c                input epsnorm is actually used

c There is a fundamental difference in the way smoothing is effected  between 
c ksmooth=1 and >1:
c ksmooth=1: damping and smoothing is effected by damping the system in  
c the way of (14.15) and (14.24) in Breviary.
c whereby epsnorm is the eps in (14.15) and epssmooth in (14.24) -  
c coupled by the epsratio input.
c
c ksmooth>1: in this case we use (14.26). See subroutine "smooth". 
c The  parameter epssmotth now has a very different meaning, as it simply  
c distributes the weights in case ksmooth=3 [ w=epssmooth for i, (1- epssmooth)/nb(i) 
c for nb(i) neighbours]. If ksmooth=2 epssmooth is  redundant, weights are simply made 
c equal. In both cases epsnorm needs  to be specified to get a norm damping (and L curve).
c
c I was a bit disappointed with the token test I did myself using  ksmooth>1, 
c it did not much speed up convergence (long ago Wim Spakman  used to claim it 
c had better convergence), and it may be wise to stick  to ksmooth=1 if you 
c dont want to lose time, because there may be  unexpected issues arising, 
c whereas ksmooth=1 seems to work fine by now.


      read(5,fmt='(a)') directory       ! directory for mat.*, aux.*
      ldr=length(directory)
      if(directory(ldr:ldr).ne.'/') then
        ldr=ldr+1
        directory(ldr:ldr)='/'
      endif
      read(5,fmt='(a)') ident   ! read ident of matrix file
      read(5,fmt='(a)') ident2  ! and of aux file (normally same)
      read(5,*) chi2target,ksmooth,epsnorm,epssmooth,epsratio
      read(5,*) outlim          ! defines outliers (*sigma)
      read(5,*) itmax1,itmax2   ! # of iters during search and zoomin
      read(5,fmt='(a)') fname   ! vertices file name
      read(5,fmt='(a)') fname2  ! facets file name
      read(5,*) colden_do       ! local smoothing? 

c document correct input:
!      open(13,file=directory(1:ldr)//'diagnostics.solve')
      open(13,file=directory(1:ldr)//'diagnostics.solve.'//ident2)
      write(13,*) 'Input for this run:'
      write(13,fmt='(a)') ident(1:70)    
      write(13,fmt='(a)') ident2(1:70)   
      write(13,*) chi2target,ksmooth,epsnorm,epssmooth,epsratio
      write(13,*) outlim           
      write(13,*) itmax1,itmax2    
      write(13,fmt='(a)') fname(1:70)    
      write(13,fmt='(a)') fname2(1:70)  
      write(6,*) 'Input for this run:'
      write(6,fmt='(a)') ident(1:70)    
      write(6,fmt='(a)') ident2(1:70)   
      write(6,*) chi2target,ksmooth,epsnorm,epssmooth,epsratio
      write(6,*) outlim           
      write(6,*) itmax1,itmax2    
      write(6,fmt='(a)') fname(1:70)    
      write(6,fmt='(a)') fname2(1:70)  

      ! checks on input
      if(ldr.le.0.or.ldr.gt.NFN) then
        write(6,*) 'Directory: ',directory(1:ldr)
        call bailout(ibr,'directory error')
      endif  
      if(ksmooth.lt.1.or.ksmooth.gt.4) then
        write(6,*) 'Bailout, ksmooth error'
        call bailout(ibr,'ksmooth must be 1-4')
      endif  
      write(6,*) 'Opening vertices file'
      open(1,iostat=ios,file=directory(1:ldr)//fname)
      if(ios.ne.0) call bailout(ibr,'Problem with vertices file')
      write(6,*) 'Opening facets file'
      open(9,iostat=ios,file=directory(1:ldr)//fname2)
      if(ios.ne.0) call bailout(ibr,'Problem with facets file')

      ! open damp_density
      open(17,file=directory(1:ldr)//'damp_density')
      ! Skip the first two lines
      read(17,*,iostat=ios) 
      read(17,*,iostat=ios) nrden 
      ALLOCATE(colden(nrden))
      do i=1,nrden
        read(17,*,iostat=ios) tmptmp, tmptmp, tmptmp, tmptmp, colvalue
        colden(i) = colvalue
      enddo
      mean_colden = sum(colden)/(max(1, nrden))
      write(*, *) "Mean of sensitivity kernels: ", mean_colden
      colden_writer = 1
      write(*,*) "Local damping: ", colden_do

! open output files
      write(6,*) 'Opening output file'
      open(4,file=directory(1:ldr)//'out.solve.'//ident2)
      write(4,fmt='("Program mpisolvetomo")')
      write(6,*) 'Time=',secnds(tcpu0)
      write(4,fmt='("Start run at: ",f8.2,/)') secnds(tcpu0)
      write(4,*)

c open matrix file, read first few records
      write(6,*) 'Opening matrix file mat.<ident>'
      open(3,iostat=ios,
     &      file=directory(1:ldr)//'mat.'//ident,
     &      form='unformatted')
      if(ios.ne.0) then
        write(6,*) 'Cannot open mat file'
        call bailout(ibr,'Error opening matrix file')
      endif  
      read(3,iostat=ios) nrow,ncol
      if(ios.ne.0) call bailout(ibr,'error readng matrix file line1')
      write(6,*) 'nrow,ncol=',nrow,ncol
      read(3,iostat=ios) jlast
      if(ios.ne.0) call bailout(ibr,'error readng matrix file jlast')
      write(6,*) 'jlast=',jlast
      read(3,iostat=ios) jssw
      if(ios.ne.0) call bailout(ibr,'error readng matrix file jssw')
      write(6,*) 'jssw=',jssw
      read(3,iostat=ios) dpar
      if(ios.ne.0) call bailout(ibr,'error readng matrix file dpar')
      write(6,*) 'dpar=',dpar
c      npar=jlast(3)
      npar=jlast(4)  ! bug fix acc. to Guust, 2011-05-16
      if(ncol.gt.NDIM) call bailout(ibr,'ncol>NDIM')

c output
      write(4,fmt='(//,"Matrix has ",i8," rows and ",i8," columns.")')
     &  nrow,ncol 
      write(6,fmt='(//,"Matrix has ",i8," rows and ",i8," columns.")')
     &  nrow,ncol 
      write(4,*) 'Chi^2/N target: ',chi2target
      write(4,*) 'Smoothing option ',ksmooth
      if(ksmooth.eq.1) then
        if(chi2target.lt.0.) then
          write(4,*) 'epsnorm=',epsnorm,', epssmooth=',epssmooth
        else 
          write(4,*) 'epsratio=',epsratio
        endif  
      else if(ksmooth.eq.2) then
        epssmooth=1.0
      else if(ksmooth.eq.3) then
        write(4,*) 'epssmooth=',epssmooth
        if(epssmooth.le.0.or.epssmooth.ge.1.0) 
     &          call bailout(ibr,'eps not OK, must be between 0 and 1')
      else if(ksmooth.eq.4) then
        call bailout(ibr, 'ksmooth=4 not yet implemented')
      endif  
      write(4,*) 'Outliers have >',outlim,' st. dev. misfit'
      write(4,*) 'Max ',itmax1,' iterations during L-curve search'
      write(4,*) 'Max ',itmax2,' iterations during rootfinding'
      write(4,fmt='(//,"List of parameters inverted for")')
      write(4,fmt='(//,"Parameter",11x,"jssw      Matrix columns",6x,
     &"dpar")')
      do i=1,10
        if(jlast(i+1).gt.jlast(i)) then
          write(4,70) chcor(i),jssw(i),jlast(i)+1,jlast(i+1),dpar(i)
        endif  
      enddo
70    format(a19,i5,i8,' to ',i8,f10.5)      
      open(7,file=directory(1:ldr)//'tex.solve.'//ident2)
      write(7,*) 'TABLE OF MODEL PARAMETERS'
      write(7,*) 'Parameter & N & $\\sigma_{prior}$ \\\\'
      do i=1,10
        if(jlast(i+1).gt.jlast(i)) then
          write(7,80) chcor(i),jlast(i+1)-jlast(i),dpar(i)
        endif  
      enddo
80    format(a19,' & ',i8,' & ',f10.5,' \\\\')
      write(6,*) 'End of user input '
c end of user in/output

c data statistics - initializing
      chi2start=0.
      ngroup=0
      do i=1,100
        chi2grp(i)=0.
        ngrp(i)=0
      enddo  

c read aux file
      lasttot=0         ! last total of nonzero Aij
      lastrow=0
      open(2,iostat=ios,file=directory(1:ldr)//'aux.'//ident2)
      if(ios.ne.0) then
        write(6,*) 'Error opening aux file'
        call bailout(ibr,'Error opening aux file')
      endif  

      ! The aux file has the data, as well as information on the
      ! sparseness of every row. 

      write(6,*) 'Starting to read aux file' ! read nrow data

100   read(2,120,iostat=ios) krow,kount,rhskrow,sigmakrow,ievtkrow,
     &      kklustkrow,stationcodekrow,networkcode,kstatkrow,ibandkrow,
     &      krtyp,klust1,igrp,grpwkrow,dtheor
120   format(2i8,e14.5,e12.3,2i8,1x,a16,1x,a8,2i8,3i7,2e14.5)
! 120   format(2i8,e14.5,e12.3,2i8,1x,a16,1x,a8,2i8,3i4,e14.5)
      if(ios.ne.0) then
        write(6,*) krow,kount,rhskrow,sigmakrow,ievtkrow,
     &      kklustkrow,stationcodekrow,networkcode,kstatkrow,ibandkrow,
     &      krtyp,klust1,igrp,grpwkrow
        print *,'read error for after row ',lastrow
        call bailout(ibr,'aux file read error')
      endif  
      if(krow.eq.0) goto 140
      if(krow.ne.lastrow+1) then
        write(6,*) 'Previous datum: ',lastrow,' current:',krow
        write(6,*) ievtkrow,stationcodekrow
        call bailout(ibr,'aux: data not in sequence')
      endif
      lastrow=krow
      rhs(krow)=rhskrow         ! datum (group-weighted chi)
      grpw(krow)=grpwkrow       ! weight for this group of data (igrp)
      chi2start=chi2start+(rhskrow/grpwkrow)**2   ! add up chi square
      chi2grp(igrp)=chi2grp(igrp)+(rhskrow/grpwkrow)**2
      ngrp(igrp)=ngrp(igrp)+1   ! count nr of data in each group
      ngroup=max(igrp,ngroup)   ! find number of data groups
      if(ngroup.gt.100) call bailout(ibr,'ngroup exceeds 100')
      sigma(krow)=sigmakrow     ! original standard deviation 
      ievt(krow)=ievtkrow       ! event number for this datum
      kklust(krow)=kklustkrow   ! kluster for this event/datum
      kstat(krow)=kstatkrow     ! station number for this datum
      iband(krow)=ibandkrow     ! frequency band for this datum
      lasttot=lasttot+kount        ! last kount (nonzero Aij)
!      if(float(lasttot).gt.bigint) call bailout(ibr,'lasttot > 2^31')
      goto 100
c end of read data

c now we know how many nonzero Aij we have we can divide them
c equally over all processors

140   nrow=lastrow              ! number of data
      write(6,*) 'nrow=',nrow
      kntaverage=lasttot/nrow+1
      write(6,*) 'kntaverage=',kntaverage
      write(4,*) 'kntaverage=',kntaverage
      write(*,*) 'lasttot=',lasttot
      write(*,*) 'nproc=',nproc
      maxperproc=lasttot/(nproc-1)+kntaverage ! # of Aij per subproc
      write(6,*) 'maxperproc=',maxperproc
      write(4,*) 'maxperproc=',maxperproc
      if(maxperproc.gt.MPIX) then
        write(6,*) 'Number of nonzero Aij=',lasttot
        write(6,*) 'Divided over ',nproc,'-1 processors gives '
        write(6,*) maxperproc,' array elements per subprocess.'
        write(6,*) 'The max allowed is:',MPIX,'. Run this with at'
        write(6,*) 'least ',lasttot/MPIX+1,' processors.'
        write(4,*) 'The max allowed is:',MPIX,'. Run this with at'
        write(4,*) 'least ',lasttot/MPIX+1,' processors.'
        call bailout(ibr,'Not enough memory!')
      else
        write(6,*) 'Number of nonzero Aij=',lasttot
        write(6,*) 'Divided over ',nproc,'-1 processors gives approx.'
        write(6,*) maxperproc,' array elements per subprocess.'
        write(4,*) 'Number of nonzero Aij=',lasttot
        write(4,*) 'Divided over ',nproc,'-1 processors gives approx.'
        write(4,*) maxperproc,' array elements per subprocess.'
      endif  
      if(ibr.eq.1) nrow=-1      ! trick to make all processes bail out
c============================================================
      endif     ! BACK TO ALL
c============================================================

c send data information to all subprocesses
c     write(6,*) 'Now broadcasting'
      call MPI_BCAST(nrow,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(nrow.lt.0) stop 'bailout'
      call MPI_BCAST(ncol,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c     write (6,*) 'Process',myid,'reached barrier 1'
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     write(6,*) 'Process',myid,' passed barrier 1'
      nrowd=nrow+2*ncol         ! nr of rows including damping
      !TODO: check if nrowd can be less if ksmooth>1
      if(nrowd.gt.NMAX) 
     &      call bailout(ibr,'nrowd > NMAX - increase dimensions')
c end of send some data information to subprocessors

c=======================================================================
      if(myid.eq.0) then
c=======================================================================
c read rest of aux file with the station and event information
      write(6,*) 'Now reading station and event lists'
      read(2,142,iostat=ios) nstat,(nn,statlist(i),netwlist(i),
     &      i=1,nstat)
      if(ios.ne.0) call bailout(ibr,'Error reading station list')
      read(2,143,iostat=ios) nklust,(nn,evlist(i),klustlist(i),
     &      i=1,nklust)
      if(ios.ne.0) call bailout(ibr,'Error reading event list')
142   format(i5,/,(i5,1x,a16,1x,a8))
143   format(i5,/,(i5,2i10))
      close(2)
      write(6,*) 'Read lists - nstat,nklust=',nstat,nklust

      chi2start=chi2start/nrow
      print *,'Misfit at start: relative chi2=',chi2start
      if(chi2target.gt.chi2start) then
        print *,'BAILOUT: Even a zero model does better than the'
        print *,'target Chi^2 specified - with outliers still present!'
        call bailout(ibr,'Chi2 target already met by background model')
      endif  


      ! TODO: root directory?
      ! find atamax from atamx.<ident> file
      atamax=0.
      open(8,iostat=ios,file=directory(1:ldr)//'atamx.'//ident)
      if(ios.ne.0) call bailout(ibr,'Could not open atamx file')
      write(6,*) 'Now reading atamx file'
      read(8,*) atamax,sqratamax
      print *,'A^TA max =', atamax
      close(8)
      if(atamax.eq.0.) call bailout(ibr,'Failed to read atamax')

      write(4,145) atamax,sqratamax
145   format('AtA max=',g12.3,' sqrt(AtA max)=',g12.3)
c end of read more data information

c read matrix, divide up among subprocesses. Note that ddif in the
c header is not scaled (though it may be de-meaned), in contrast to
c rhs, which is scaled by the standard error s and groupweight.
c In the inversion, we use rhs (hence its name)
      ip=0
      lasttot=0
      lastrow=0
      kdisp(0)=nrow            ! displacement of u "buffer" at root
      knt(0)=2*ncol            ! nr of elements in u buffer
      kdisp(1)=0               ! subprocess 1 start at u(1) 

      nrow_kntaverage = nrow
      nrow_kntaverage = nrow_kntaverage*kntaverage

      write(6,*) 'Matrix split between sub-processors:'
      write(6,*) 'Proc    offset    #rows   nonzeros      time'
150   ip=ip+1
c     write(6,*) 'ip=',ip,' in root while reading matrix file'
      if(ip.ge.nproc) call bailout(ibr,'Miscalculation ip>nproc')
      ! safeguard agains roundoff:
      if(ip.eq.nproc-1) maxperproc=min(MPIX,maxperproc+nrow_kntaverage)
160   read(3,iostat=ios,end=161) krow,kount,kband,ddif,s,kklustrow,
     &     kstatrow,groupweight
      goto 162
161   krow=0            ! fix matrix files that end with eof io krow=0
      ios=0
162   if(ios.ne.0) then
        write(6,*) 'Read error - last & current row nr:',lastrow,krow
        write(6,*)  krow,kount,kband,ddif,s,kklustrow,kstatrow,
     &        groupweight
        call bailout(ibr,'Matrix file read error')
      endif
      if(krow.gt.0.and.krow.ne.lastrow+1) then
        write(6,*) 'Matrixfile rows not in sequence error'
        call bailout(ibr,'Matrixfile rows not in sequence')
        if(ibr.lt.0) lasttot=-1         ! bailout flag
        goto 164
      endif 
      ! buffer full or last row? then empty it 
      if(lasttot+kount.gt.maxperproc.or.krow.eq.0) then
        write(*,*) '==============================================='
        write(*,*) 'nrow*kntaverage', nrow*kntaverage
        write(*,*) 'nrow_kntaverage', nrow_kntaverage
        write(*,*) 'maxperproc', maxperproc
        write(*,*) 'nrow', nrow
        write(*,*) 'kntaverage', kntaverage
        write(*,*) 'MPIX', MPIX
        write(*,*) '----------------------------------------------'

        backspace 3
        knt(ip)=lastrow-kdisp(ip)
        kdisp(ip+1)=lastrow     ! array u offset at next subprocess
        write(6,163) ip,kdisp(ip),knt(ip),lasttot,secnds(tcpu0)
163     format(i5,3i10,f10.2)
164     ktag=1000+ip
        call MPI_SEND(lasttot,1,MPI_INTEGER,ip,ktag,
     &        MPI_COMM_WORLD,ierr)
        if(ibr.lt.0) stop 'bailout'
        ktag=2000+ip
        call MPI_SEND(ja,lasttot,MPI_INTEGER,ip,ktag,
     &        MPI_COMM_WORLD,ierr)
        ktag=3000+ip
        call MPI_SEND(a,lasttot,MPI_REAL,ip,ktag,
     &        MPI_COMM_WORLD,ierr)
c       write(6,*) 'Root sent a and ja to process ',ip
        lasttot=0
        if(krow.gt.0.and.krow.lt.nrow) goto 150
        write(6,*) 'End of matrix reached krow,nrow=',krow,nrow
      else
        read(3,iostat=ios) (ja(lasttot+i),i=1,kount)
        read(3,iostat=ios) (a(lasttot+i),i=1,kount)
        ia(krow)=kount
        lasttot=lasttot+kount
        lastrow=krow
        goto 160
      endif  
      close(3)

c end of read matrix, divide up and send to subprocessors

c read interior vertices in root
      write(6,*) 'Root reading vertices file'
      read(1,*) nd
      if(nd.ne.3) call bailout(ibr,'vertices file: nd not equal to 3')
      read(1,*) np
      if(np.gt.NPMAX) call bailout(ibr,'np>NPMAX')
      do i=1,np  
        read(1,*) (points(j,i),j=1,3)
      enddo  
      write(6,*) 'Root read vertices file succesfully np=',np
      close(1)

c zero array nb (nr of neighbours of node i)
      do i=1,np
        nb(i)=0
      enddo  

c  read vertices
      read(9,*) nt
      write(6,*) 'Root reading ',nt,' tetrahedra'
      do jtetra=1,nt
        read(9,*) nv
	do i=1,3
	  ii=nv(i)+1    ! qhull numbers from 0, so add 1
	  do j=i+1,4
	    jj=nv(j)+1
	    call adn(ii,jj,nb,iv)	! add to neighbour list
	  end do
	end do
      end do
      write(6,*) 'Succesful read of facets file'
      close(9)

c find approximate volume of tetrahedron
c     write(6,*) 'calling aprvol'
      call aprvol(np,nb,iv,points,tvol,tnab,vscale)
c end of bookkeeping for model grid

c prepare root procesor for next phase: iterating the L curve

      ! GMT file for L curve
      write(6,*) 'Opening ',directory(1:ldr)//'xy.tradeoff.'//ident2
      open(8,file=directory(1:ldr)//'xy.tradeoff.'//ident2) 
      write(6,*) 'Root opened xy.tradeoff file'
      write(8,*) 0.,chi2start
      ! GMT for L curve incl damping
c     open(9,file=directory(1:ldr)//'xy.restot.'//ident2) 
c     write(9,*) 0.,chi2start
      write(4,*)
      write(4,*) 'Next table has parameter values scaled by dpar.'
      write(4,*) 'Trials with reduced iterations in LSQR, itmax=',itmax1
      write(4,190) chi2start
      write(6,*) 'Next table has parameter values scaled by dpar.'
      write(6,*) 'Trials with reduced iterations in LSQR, itmax=',itmax1
      write(6,190) chi2start,0.
190   format(//,'Run  Chi^2(rel)',3x,'Max-par   RMS-par   Max-cor   ',
     & 'RMS-cor       CPU    ytry    RMStot      epsnorm    epssmooth',
     & /,' 1a',f12.3,5(9x,'0'),3x,'1.000',f10.4,'  (with outliers)')

      write(7,*) 'TABLE OF INVERSION RUNS'
      write(7,191) chi2start
191   format('Run & $\\chi^2_{rel}$ & Max-par & RMS-par & Max-cor & ',
     & 'RMS-cor & CPU \\\\',/,'  0 & ',f10.2,5(' & 0'),' \\\\')

      
      ! set parameters for starting model (infinite damping) as itry=1
!      ytry(1)=ytry_lcurve(1)               ! y=eps/(1+eps) so that 0<y<1
      ytry(1)=1.0               ! y=eps/(1+eps) so that 0<y<1  ****
      yup=1.
      chi2up=chi2start

      klean = 0
      ytry(2)=ytry_lcurve(2)   ! y for second run (almost no damping)
!     ytry(2)=0.03              ! y for second run (almost no damping)/CHANGED 2007/10/18
      itmax=itmax1              ! max nr of iterations during trials
      ylow=ytry(2)              ! bracket between ylow/yup or hi/lo chi2
      jzoom=0                   ! rough search first (few iters in lsqr)

      if(ibr.eq.1) itmax=-1

c end of initializations on root

c============================================================
      else              ! SUB-PROCESSES
c============================================================

c root has been sending stuff - absorb it all
c     write(6,*) 'Going into ip loop, myid',myid
      do ip=1,nproc-1
c       write(6,*) 'loop ip=',ip,' for proces',myid,nrow,ncol
        if(myid.eq.ip) then
c         write(6,*) 'process ',myid,' nrow,ncol=',nrow,ncol
          ktag=1000+ip
          call MPI_RECV(lasttot,1,MPI_INTEGER,0,ktag,
     &      MPI_COMM_WORLD,status,ierr)
          if(lasttot.lt.0) stop
          ktag=2000+ip
          call MPI_RECV(ja,lasttot,MPI_INTEGER,0,ktag,
     &      MPI_COMM_WORLD,status,ierr)
          ktag=3000+ip
          call MPI_RECV(a,lasttot,MPI_REAL,0,ktag,
     &      MPI_COMM_WORLD,status,ierr)
c         write(6,*) 'Process ',myid,' received ja and a, t=',
c    &          secnds(tcpu0)
        endif
      enddo  
c============================================================
      endif     ! BACK TO ALL
c============================================================

c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ia,nrow,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kdisp,nproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(knt,nproc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(itmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(itmax.eq.-1) stop 'bailout'
c     write(6,*) 'Process ',myid,' received ia',knt(nproc-1),
c    &  kdisp(nproc-1)
c end of read matrix in subprocessors


c ACTIVATE THIS TEST (set jdebug=1) IF YOU TINKER WITH THE MATRIX
c CODE. MAKE SURE YOU MAKE THE SAME CHANGES IN THE TEST CODE
c jdebug = 1 overwrites epssmooth and epsnorm!
      jdebug=0
      if(jdebug.eq.0) goto 199

c     INTERMEZZO: DO THE CLAERBOUT TEST

c     write(6,*) 'Process ',myid,' now in Claerbout test'
      ! make sure we're all ready to start number crunching
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c=======================================================================
      if(myid.eq.0) then
c=======================================================================
c     write(6,*) 'root passed Barrier'
      
      ! First test smoother S (subprocesses are idle)
      epsnorm=0.
      epssmooth=0.3
      idum=-1
      call randomize(idum,v,ncol)
      call randomize(idum,u,ncol)
      call nullify(w,ncol)
      call nullify(w1,ncol)
      ! w contains Sv
      call smooth(v,w,ncol,epsnorm,epssmooth,nb,iv,ksmooth,jlast)     
      ! w1 contains S'u
      call tsmooth(u,w1,ncol,epsnorm,epssmooth,nb,iv,ksmooth,jlast)    
      call dotp(w,u,ncol,p1)            ! (Sv,u) should equal:
      call dotp(w1,v,ncol,p2)           ! (v,S'u)
      write(13,*) 'Test of smooth and tsmooth:'
      write(13,*) 'epss=0.3, (Sv,u)=',p1,', (v,S^Tu)=',p2
      write(13,*) 'relative error ',abs((p2-p1)/p1)
      write(13,*)
      write(6,*) 'Test of smooth and tsmooth:'
      write(6,*) 'epss=0.3, (Sv,u)=',p1,', (v,S^Tu)=',p2
      write(6,*) 'relative error ',abs((p2-p1)/p1)
      write(6,*)

      ! now test matrix A (involves subprocesses also)
      epsnorm=sqratamax
      epssmooth=sqratamax
      if(ksmooth.gt.1) epssmooth=0.5
      call randomize(idum,v,ncol)       ! New random vectors
      call randomize(idum,u,nrowd)
      do i=1,nrowd              ! save u
        rhs(i)=u(i)
      enddo
      call nullify(w,nrowd)
      call nullify(w1,ncol)
c============================================================
      endif     ! BACK TO ALL
c============================================================

      call MPI_BCAST(u,nrowd,MPI_REAL,0,MPI_COMM_WORLD,ierr)

c     write(6,*) myid,' did bcats of u'


c=======================================================================
      if(myid.eq.0) then
c=======================================================================

c root does the damping part. This must be the same as further on in
c the program, so compare line-by-line if you make changes.

c     Compute A(transpose)*u and store this time in w2
      i=nrow
      do ii=1,ncol              ! norm damping (diagonal matrix)
        i=i+1                   ! true row index
        x=1.0
        if(ii.le.jlast(4)) x=tvol(mod(ii-1,np)+1)
        if(ii.gt.NDIM) then ! debug
          write(6,*) 'ii=',ii,'>NDIM=',NDIM
          write(6,*) 'ncol=',ncol
          stop 'Index ii too large'
        endif  
        w1(ii)=x*epsnorm*u(i)      
      enddo
      if(ksmooth.le.1.and.epssmooth.ne.0) then
        do ii=1,jlast(4)                 ! smoothness damping
          i=i+1                   ! true row index
          kk=mod(ii-1,np)+1
          i0=ii-kk
          w1(ii)=w1(ii)+tnab(kk)*epssmooth*u(i)
          x=tnab(kk)*epssmooth/nb(kk)
          do j=1,nb(kk)
            jj=iv(j,kk)+i0
            if(jj.gt.NDIM) then ! debug
              write(6,*) 'jj=',jj,'>NDIM=',NDIM
              write(6,*) 'ncol=',ncol,'iv,i0=',iv(j,kk),i0
              stop 'Delauney index out of bounds'
            endif  
            w1(jj)=w1(jj)-x*u(i)
          enddo
        enddo
      endif  
      do i=1,ncol       ! needed? 
        w2(i)=0.
      enddo
c     write(6,*) 'Root did A^T*u'
c============================================================
        else              ! SUB-PROCESSES
c=============================================================

c       A(transpose)*u 
        do j=1,ncol
          w1(j)=0.
        enddo 
        j=0
        do irow=kdisp(myid)+1,kdisp(myid)+knt(myid)
          if(.not.outl(irow)) then
            do k=1,ia(irow)
              j=j+1
              jj=ja(j)
              w1(jj)=w1(jj)+a(j)*u(irow)
            enddo
          else
            j=j+ia(irow)
          endif
        enddo  
c       write(6,*) 'Process',myid,' did A^T*u'
c============================================================
      endif     ! BACK TO ALL
c============================================================
c gather the w1's from each subprocess and sum them into w2
c     write(6,*) myid,' going into reduce'
      call MPI_REDUCE(w1,w2,ncol,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
c     write(6,*) 'process',myid,' reduced w1 into w2'


c     broadcast v to subprocesses
      call MPI_BCAST(v,ncol,MPI_REAL,0,MPI_COMM_WORLD,ierr)
c     write(6,*) 'broadcast of v, myid=',myid

      do i=1,nrowd       ! zero u to maintain compatibility in the code
        u(i)=0.          ! (starting u is not 0 in the actual program)
      enddo
      
c now compute A*v and store in u
c=======================================================================
      if(myid.eq.0) then
c=======================================================================
c root does the damping part of Av
      i=nrow              ! skp the regular data segment of A
      do ii=1,ncol        ! norm damping of all paramaters incl corr
        i=i+1                   ! true row index
        x=1.0
        if(ii.le.jlast(4)) x=tvol(mod(ii-1,np)+1)
        u(i)=u(i)+x*epsnorm*v(ii)
      enddo
      ! smoothing constraints within matrix A ?
      if(ksmooth.le.1.and.epssmooth.ne.0.) then
        do ii=1,jlast(4)    ! smoothness damping over model only
          i=i+1             ! true row index
          kk=mod(ii-1,np)+1 ! model node
          i0=ii-kk          ! remember shift in model vector
          u(i)=u(i)+tnab(kk)*epssmooth*v(ii)  ! scale kk using tnab()
          x=tnab(kk)*epssmooth/nb(kk)     
          do j=1,nb(kk)
            u(i)=u(i)-x*v(iv(j,kk)+i0)
          enddo  
        enddo  
      endif
c       write(6,*) 'Process',myid,' did A*v'
c============================================================
        else              ! SUB-PROCESSES
c=============================================================
        j=0
        do irow=kdisp(myid)+1,kdisp(myid)+knt(myid)
          if(outl(irow))then
            u(irow)=0.
            j=j+ia(irow)
          else 
            do k=1,ia(irow)
              j=j+1
              u(irow)=u(irow)+a(j)*v(ja(j))
            end do
          endif  
        end do
c       write(6,*) 'Process',myid,' did A*v'
c============================================================
      endif     ! BACK TO ALL
c============================================================

! Note: all GATHERV statements have been replaced with less
! efficient SEND/RECV until it is clear how to define knt and
! kdisp correctly for GATHERV.
! When that is clear, (1) remove the "do ip" loops, (2) uncomment
! the call the GATHERV, (3) uncomment the mapping from w3 to u in
! root. There are four segments like this one, the last has ktag
! starting at 7000.
! It may be doable to rplace w3 by u and gather "in place" but again the
! mpi documentation was too vague for me to risk that. The SEND/RECV
! version is effectively "in place".

c     call MPI_GATHERV(u(kdisp(myid)+1),knt(myid),MPI_REAL,w3,knt,
c    &      kdisp,MPI_REAL,0,MPI_COMM_WORLD,ierr)
c replacement for gatherv:
        do ip=1,nproc-1
          ktag=4000+ip
          if(myid.eq.ip) then
            call MPI_SEND(u(kdisp(myid)+1),knt(myid),MPI_REAL,0,ktag,
     &        MPI_COMM_WORLD,ierr)
          else if(myid.eq.0) then
            call MPI_RECV(u(kdisp(ip)+1),knt(ip),MPI_REAL,ip,ktag,
     &        MPI_COMM_WORLD,status,ierr)
          endif
        enddo
c     write(6,*) 'Gather myid, offset, knt=',myid,kdisp(myid),knt(myid)

c=======================================================================
      if(myid.eq.0) then
c=======================================================================
! CHECK!
      call dotp(rhs,u,nrowd,p1)           ! (rhs,Av) should equal:
c     call dotp(rhs,w3,nrowd,p1)           ! (rhs,Av) should equal:
      call dotp(w2,v,ncol,p2)           ! (A'rhs,v) 
      write(13,*) 'Test of A*v and A(transpose)*w:'
      write(13,*) 'epsn=',sqratamax
      write(13,*) 'epss=',epssmooth,', (w,Av)=',p1,', (A^Tw,v)=',p2
      write(13,*) 'relative error ',abs((p2-p1)/p1)
      write(6,*) 'Test of A*v and A(transpose)*w:'
      write(6,*) 'epsn=',sqratamax
      write(6,*) 'epss=',epssmooth,', (w,Av)=',p1,', (A^Tw,v)=',p2
      write(6,*) 'relative error ',abs((p2-p1)/p1)
c end of Claerbout's dot product test on root

c============================================================
        endif              ! BACK TO ALL PROCESSES
c============================================================

      call bailout(ibr,'Did Claerbout test only')
c     END OF INTERMEZZO if(jdebug.eq.1)

199   continue

c============================================================

      ! make sure we're all ready to start number crunching
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     write(6,*) myid,' passed barrier before label 200'

c MAIN LOOP  over damping parameter ytry
      itry=1                   ! first try was for 0 change to bg model
      kready=0                 ! Do more inversion runs itry? y/n
                               ! (this is set to 1 once chi2 is OK)
                               
200   itry=itry+1              ! next run (starts with itry=2!)

ccccc  HACK -- delete, is just to shorten Lcruve runs
c      if(itry.gt.3) goto 800    !!!!!!!!!!!!!!!!!!!!!!
ccccc  END HACK -- delete, is just to shorten Lcruve runs


      if(kready.eq.1) goto 800 ! probably not needed
c end of start LSQR iterations
      
c lsqr initializations

c=======================================================================
201   if(myid.eq.0) then       ! IF ROOT PROCESS
c=======================================================================

      if(itry.ge.itrymax) call bailout(ibr,'Maximum trials reached')
      
      eps=ytry(itry)/(1.0-ytry(itry))   ! eps=y/(1-y) ie 0<eps<inf
      write(13,*) 'itry=',itry,' eps=',eps
      if(chi2target.ge.0.) then
        if(ksmooth.eq.1) then  
          epssmooth=eps*sqratamax   ! scale w.r.t. matrix elements
          epsnorm=epsratio*epssmooth  ! all eps are in common damping
        else
          epsnorm=eps*sqratamax
        endif  
      endif  
      
c solve
      ! fill data vector u (was destroyed if lsqr preceded this)
      do i=1,nrow
        u(i)=rhs(i)
      enddo
      write(13,*) 'nrow,nrowd,ncol=',nrow,nrowd,ncol
      do i=nrow+1,nrowd         ! damping part of rhs = 0
        u(i)=0.
      enddo  
      call normlz(nrow,u,beta)  ! nrow i.o. nrowd because rest is 0
      b1=beta 
c============================================================
      endif     ! BACK TO ALL
c============================================================
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      !write(6,*) myid,' passed barrier R12, t=',secnds(tcpu0)
      call MPI_BCAST(u,nrowd,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      !write(6,*) myid,' now has u; starts A(transpose)*u',u(nrow+1)
      !write(6,*) 'at t=',secnds(tcpu0)
 
c=======================================================================
      if(myid.eq.0) then
c=======================================================================
c we do lsqr within the main body of the program, to clarify hierarchy
c if ksmooth equals 1, the smoothing matrix S equals I.

c since damping part of u=0 at start, damping results in 0's for
c A(transose)*u:
      do i=1,ncol 
        sol(i)=0.
        v(i)=0.
        w1(i)=0.    ! w1 is "vtilde" in my notes, i.e. model v=S*w1
      end do

c root process does the damping part of matrix multiplications
c but at initialization the tail of u is zero, so we remain idle
c============================================================
        else              ! SUB-PROCESSES
c=============================================================

c compute A(transpose)*u for part of A and store in w1
c equivalent in lsqr:  call atupv(ksmooth,m,n,u,v,outl)
      do j=1,ncol
        w1(j)=0.
      enddo 
      j=0
      do irow=kdisp(myid)+1,kdisp(myid)+knt(myid)
        if(.not.outl(irow)) then
          do k=1,ia(irow)
            j=j+1
            jj=ja(j)
            w1(jj)=w1(jj)+a(j)*u(irow)
          enddo
        else
          j=j+ia(irow)          ! skip outliers
        endif
      enddo  

c============================================================
      endif     ! BACK TO ALL
c============================================================
c      write(6,*) myid,'passed barrier before R13 t=',secnds(tcpu0) !!XX
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)    !!XX

c gather the w1's from each subprocess and sum them into v
      call MPI_REDUCE(w1,v,ncol,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      !write(6,*) myid,' reduced w1 into v (R13) t=',secnds(tcpu0)

c=======================================================================
      if(myid.eq.0) then
c=======================================================================
c in transpose multiplication, smoothing with S(transpose) comes last
c (no need to add to original v since that v initialized to 0)

      if(ksmooth.gt.1) then
        do i=1,ncol
          w1(i)=v(i)
        enddo
        call tsmooth(w1,v,ncol,epsnorm,epssmooth,nb,iv,ksmooth,jlast)
      endif  

c sum v is now known, normalize and copy into w:
c equivalent statement in lsqr: do i=1,n {w(i)=v(i)}
      call normlz(ncol,v,alfa)
      write(13,*) 'after atupv alfa=',alfa
      rhobar=alfa
      phibar=beta
      ! map back to v
      do i=1,ncol
        w(i)=v(i) 
      end do  
      write(13,202) 0,1.0,beta,secnds(tcpu0)
202   format('iter      xabs   xnrm(%)',9x,'r    phibar      time',/,
     &        i4,22x,3f10.3)
c end of First backprojection

c============================================================
      endif              ! BACK TO ALL PROCESSES
c============================================================

c lsqr iterations


      do iter=1,itmax 

c        write(6,*) myid,' starts iter ',iter,'at t=',secnds(tcpu0) !!XX
        call MPI_BARRIER(MPI_COMM_WORLD,ierr) !!XX
        !write(6,*) 'at t=',secnds(tcpu0)
c=======================================================================
        if(myid.eq.0) then       ! IF ROOT PROCESS
c=======================================================================

c       for AS*v, smooth first if needed
        if(ksmooth.gt.1) then
          call smooth(v,w1,ncol,epsnorm,epssmooth,nb,iv,ksmooth,jlast)

        else

          do i=1,ncol
            w1(i)=v(i)  ! no smoothing, but we work with w1
          enddo

        endif

        aa=-alfa
        do i=1,nrowd 
          u(i)=aa*u(i) 
        end do 

c============================================================
      endif     ! BACK TO ALL
c============================================================
c        write(6,*) myid,' pass barr before R14 t=',secnds(tcpu0) !!XX
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)  !!XX

c       broadcast w1 (smoothed v)  to subprocesses
        call MPI_BCAST(w1,ncol,MPI_REAL,0,MPI_COMM_WORLD,ierr)

c       broadcast new (scaled) u to subprocesses
        call MPI_BCAST(u,nrowd,MPI_REAL,0,MPI_COMM_WORLD,ierr)

c=======================================================================
      if(myid.eq.0) then
c=======================================================================

c root does the damping part of ASv+u
c in lsqr we have here: call avpu(ksmooth,m,n,u,v,outl)

        open(18, file="check_scaling.txt")
        i=nrow              ! skp the regular data segment of A
        do ii=1,ncol        ! norm damping of all paramaters incl corr
          i=i+1                   ! true row index
          x=1.0
          if(ii.le.jlast(4)) x=tvol(mod(ii-1,np)+1) ! grid node
          u(i)=u(i)+x*epsnorm*w1(ii)
        enddo
        ! smoothing constraints within matrix A ?
        if(ksmooth.le.1.and.epssmooth.ne.0.) then
          do ii=1,jlast(4)    ! smoothness damping over model only
            i=i+1             ! true row index
            kk=mod(ii-1,np)+1 ! model grid node

            if (colden_do.eq.1) then
              xcol=1.+colden(mod(ii-1,np)+1)/10.
            else
              xcol=1.
            endif

            if (colden_writer.eq.1) then
              write(18, *) colden(mod(ii-1,np)+1), xcol,
     &                   xcol*x*epsnorm, u(i)
            endif


            i0=ii-kk          ! remember shift in model vector
            u(i)=u(i)+xcol*tnab(kk)*epssmooth*w1(ii)  ! scale kk using tnab()
            x=tnab(kk)*epssmooth/nb(kk)     
            do j=1,nb(kk)
              u(i)=u(i)-xcol*x*w1(iv(j,kk)+i0)
            enddo  
          enddo  
        endif  
        colden_writer = 0
c       !write(6,*) 'root did his damping part of ASv+u'
c============================================================
        else              ! SUB-PROCESSES
c=============================================================

c       ASv+u 
        j=0
        do irow=kdisp(myid)+1,kdisp(myid)+knt(myid)
          if(outl(irow))then
            u(irow)=0.
            j=j+ia(irow)
          else 
            do k=1,ia(irow)
              j=j+1
              u(irow)=u(irow)+a(j)*w1(ja(j))
            end do
          endif  
        end do
        !write(6,*) myid,' did ASv+u t=',secnds(tcpu0)

c============================================================
      endif     ! BACK TO ALL
c============================================================
          
c       write(6,*) myid,' passed barrier after R16 t=',secnds(tcpu0) !!XX
       call MPI_BARRIER(MPI_COMM_WORLD,ierr) !!XX
        
c root has done its job for ASv+u, now get the other segments of u
        call MPI_GATHERV(u(kdisp(myid)+1),knt(myid),MPI_REAL,w3,knt,
     &        kdisp,MPI_REAL,0,MPI_COMM_WORLD,ierr)

c=======================================================================
      if(myid.eq.0) then
c=======================================================================

c transfer in case gatherv was used to transmit u:
        do i=1,nrowd
          u(i)=w3(i)
        enddo  

c       write(6,*) 'root now going into norml for u'
        call normlz(nrowd,u,beta) 

        b=-beta
        do i=1,ncol
          v(i)=b*v(i) 
        end do 
        !write(6,*) 'root now before barrier R17 t=',secnds(tcpu0)

c============================================================
      endif     ! BACK TO ALL
c============================================================

        
c        write(6,*) myid,' passd barr before R17 t=',secnds(tcpu0) !!XX
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)  !!XX

        ! update scaled u at each subrocess
        call MPI_BCAST(u,nrowd,MPI_REAL,0,MPI_COMM_WORLD,ierr)

c       write(6,*) myid,' now has u=',u(1),u(2),u(3),'...'

c       broadcast new v to subprocesses
        call MPI_BCAST(v,ncol,MPI_REAL,0,MPI_COMM_WORLD,ierr)

c       write(6,*) myid,' now has v=',v(1),v(2),v(3),'...'

c=======================================================================
      if(myid.eq.0) then
c=======================================================================
c root does the damping part of A(transpose)*u = w1
        i=nrow
c       write(6,*) 'Root starts damping ', nrow,ncol,np
        do ii=1,ncol              ! norm damping (diagonal matrix)
          i=i+1                   ! true row index
          x=1.0
          if(ii.le.jlast(4)) x=tvol(mod(ii-1,np)+1)     ! grid node
          w1(ii)=x*epsnorm*u(i)      
        enddo

        if(ksmooth.le.1.and.epssmooth.ne.0) then

          do ii=1,jlast(4)                 ! smoothness damping
            i=i+1                   ! true row index
            kk=mod(ii-1,np)+1       ! model grid index

            if (colden_do.eq.1) then
              xcol=1.+colden(mod(ii-1,np)+1)/10.
            else
              xcol=1.
            endif

            i0=ii-kk                ! offset    
            w1(ii)=w1(ii)+xcol*tnab(kk)*epssmooth*u(i)  ! central node
            x=tnab(kk)*epssmooth/nb(kk)
            do j=1,nb(kk)
              jj=iv(j,kk)+i0
              w1(jj)=w1(jj)-xcol*x*u(i)  ! minus neighbours  
            enddo
          enddo
        endif  
c       write(6,*) 'root did its part of A(transpose)*u'

c============================================================
        else              ! SUB-PROCESSES
c=============================================================
c       compute w1=A(transpose)*u (does not yet sum v because 
c       of later smoothing with S)
        do j=1,ncol
          w1(j)=0.
        enddo 
        j=0
        do irow=kdisp(myid)+1,kdisp(myid)+knt(myid)
          if(.not.outl(irow)) then
            do k=1,ia(irow)
              j=j+1
              jj=ja(j)
              w1(jj)=w1(jj)+a(j)*u(irow)
            enddo
          else
            j=j+ia(irow)
          endif
        enddo  
c       write(6,*) myid,' did A(tranpose)*u'

c============================================================
      endif     ! BACK TO ALL
c============================================================
c        write(6,*) myid,' pass bar before R19 t=',secnds(tcpu0) !!XX
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)  !!XX

c collect A(transpose)*u from the subprocesses, sum into w2
        call MPI_REDUCE(w1,w2,ncol,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,
     &        ierr)
c       !write(6,*) myid,' reduced w1 to w2 to get A(transpose)*u'

c=======================================================================
      if(myid.eq.0) then
c=======================================================================
c product with S(transpose), if any, and sum to v
        if(ksmooth.gt.1) then
          call tsmooth(w2,w3,ncol,epsnorm,epssmooth,nb,iv,ksmooth,jlast)
          do i=1,ncol
            v(i)=w3(i)+v(i)
          enddo
        else
          do i=1,ncol
            v(i)=w2(i)+v(i)
          enddo
        endif  
c       write(6,*) 'root now has v'

c       write(6,*) 'root now normalizing v'
        call normlz(ncol,v,alfa)
        rho=sqrt(rhobar*rhobar+beta*beta)
        c=rhobar/rho
        s=beta/rho 
        teta=s*alfa
        rhobar=-c*alfa
        phi=c*phibar
        phibar=s*phibar 
        t1=phi/rho
        t2=-teta/rho
        xabs=0.
        xnrm2=0
        do i=1,ncol
          sol(i)=t1*w(i)+sol(i)
          xabs=xabs+abs(sol(i))         ! compute |x|
          xnrm2=xnrm2+sol(i)*sol(i)     ! and |x^2|
          w(i)=t2*w(i)+v(i)
        end do 
        rr=phibar/b1
        ! file output: diagnostics.solve.<ident2>
        ! xabs  : solution length. xnrm2 = xabs*xabs
        ! phibar: length of the misfit vector. 
        ! rr is phibar scaled by the initial misfit, so that it is 
        ! easier to judge convergence. 
        ! rhobar is not important (it is a measure if the length  
        ! of the next basis function 
        ! in model space and when it gets small it  probably - but not 
        ! certainly - means that  there are no large eigenvalues  left 
        ! to explore)
        !     write(13,203) iter,xabs,sqrt(xnrm2)/ncol,rr,rhobar,phibar
        !
        ! iter      xabs      xnrm         r    rhobar    phibar
        !  0                                  69.100  1421.242
        !  1    49.9         0.000     0.970  -136.996  1379.174
        !  2    69.0         0.000     0.959   121.662  1362.756  
        !
!        write(13,203) iter,xabs,sqrt(xnrm2)/ncol,rr,rhobar,phibar
!203     format(i4,2g12.3,3f10.3)

        write(13,203) iter,xabs,100.*sqrt(xnrm2)/ncol,rr,phibar,
     &        secnds(tcpu0)
203     format(i4,g12.3,4f10.3)

        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
        ! EXPERIMENTAL 2007/07/24: 
c       output of chi2 by group in every iteration
c       1) Starting chi2 is in chi2grp(1:ngroup) and chi2start
c       2) ngrp(1:ngroup) gives number of data in each group.
c          grpw(irow) gives weights.
c       solvetomo reads rhskrow0. assemblematrix has already 
c       multiplied each datum by its weight in inversion: 
c       rhskrow0=rhs(krow0)*chiw(igroup)  -- in assemblematrix
c       But for calculating chi2 we want
c       the absolute numbers in rhs, so solvetomo undoes weigthing
c       by dividing through grpwkrow
c      rhs(krow)=rhskrow         ! datum (normalized)
c      grpw(krow)=grpwkrow       ! weight for this group of data (igrp)
c      chi2start=chi2start+(rhskrow/grpwkrow)**2   ! add up chi square
c      chi2grp(igrp)=chi2grp(igrp)+(rhskrow/grpwkrow)**2
c
c       3) Work down res vector; nrow data in rhs(i) and u(i)
c      do i=1,nrow
c        if(.not.outl(i))then 
c          res=rhs(i)-u(i)
c          chi2=chi2+(res/grpw(i))**2    ! Chi^2 (unweighted residuals)
c        endif  
c      enddo
c      chi2try(itry)=chi2/(nrow-noutl)


c declare: noutl_grp(igrp),chi2_grp(igrp),chi2i_all,chi2g

      if(0.gt.1)then   !!! TESTKEEP OUT
      ioo = 0
      resii=0.
      do igrp=1,ngroup
        iaa=ioo+1
        ioo=ioo+ngrp(igrp)
        noutl_grp(igrp)=0 
        chi2_grp(igrp) =0.
       do i=iaa,ioo
        if(.not.outl(i))then 
          resii=rhs(i)-u(i)
          dchi=(resii/grpw(i))**2
          chi2i_all=chi2i_all+dchi
          chi2_grp(igrp)=chi2_grp(igrp)+dchi  
        else 
          ! count number of outliers in each group
          noutl_grp(igrp)=noutl_grp(igrp)+1 
        endif 
       enddo
       chi2_grp(igrp)=chi2_grp(igrp)/(ngrp(igrp)-noutl_grp(igrp) )
ccc      chi2tryi(itry)=chi2i/(nrow-noutl)
      enddo
      chi2i_all=chi2i_all/(nrow-noutl) 

      ! output to a new diagnostics file
      endif ! if(0.gt.1)then
   

        !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 


c============================================================
      endif     ! BACK TO ALL
c============================================================

      end do    ! end of lsqr iterations
      !write(6,*) myid,' arrived at end of lsqr t=',secnds(tcpu0)

      ! After first iteration (itry=2):
      ! Identify outliers, then redo undamped inversion
      ! without them
      
      if(klean.eq.0) then       ! next segment is only done once

        !write(6,*) myid,' into klean=0 segment'
c=======================================================================
        if(myid.eq.0) then       ! IF ROOT PROCESS
c=======================================================================

          ! apply final smoothing sol=Sv in case of implicit smoothing
          if(ksmooth.gt.1)  then
cc          write(6,*) 'root calls smooth'
           call smooth(sol,v,ncol,epsnorm,epssmooth,nb,iv,ksmooth,jlast)
          else
            do i=1,ncol
              v(i)=sol(i)
            enddo  
          endif  
         
c============================================================
      endif     ! BACK TO ALL
c============================================================

          !write(6,*) myid,' arrived before R20 t=',secnds(tcpu0)
c         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

          !write(6,*) myid,' passed barrier before R20 t=',secnds(tcpu0)
          call MPI_BCAST(v,ncol,MPI_REAL,0,MPI_COMM_WORLD,ierr)

c=======================================================================
      if(myid.eq.0) then
c=======================================================================

          do i=nrow+1,nrowd      ! ignore damping part
            u(i)=0.
          enddo  
c         !write(6,*) 'root set damping part of u to 0'

c============================================================
        else              ! SUB-PROCESSES
c=============================================================
c compute data fit (former subroutine asol)

c compute A*v (formerly asol)

         j=0
         do irow=kdisp(myid)+1,kdisp(myid)+knt(myid)
           u(irow)=0.
           if(outl(irow))then
             j=j+ia(irow)
           else 
             do k=1,ia(irow)
               j=j+1
               u(irow)=u(irow)+a(j)*v(ja(j))
             end do
           endif  
         end do
c        write(6,*) myid,' computed u = A*v'
c============================================================
      endif     ! BACK TO ALL
c============================================================

         klean=1

      ! Gather results
c        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         !write(6,*) myid,' klean=0, passed S21 t=',secnds(tcpu0)
         call MPI_GATHERV(u(kdisp(myid)+1),knt(myid),MPI_REAL,w3,knt,
     &        kdisp,MPI_REAL,0,MPI_COMM_WORLD,ierr)

c=======================================================================
      if(myid.eq.0) then
c=======================================================================

         ! Record outliers: file outl.<ident2>
         fname4(1:5)='outl.'
         fname4(6:NFN)=ident2(1:NFN-5)
         open(21,file=directory(1:ldr)//fname4)  ! file outl.<ident>

c        write(6,*) 'root starts computing chi2start'
         chi2start=0.           ! recompute chi2start
         noutl = 0 

c if gatherv was used, map w3 into u
         do i=1,nrowd
           u(i)=w3(i)
         enddo  
         do i=1,nrow
           if(abs(u(i)-rhs(i)).gt.outlim*grpw(i)) then
             outl(i) = .TRUE.
             rhs(i)  = 0.
             noutl=noutl+1
             write(21,*) 1   ! file outl.<ident>
           else
             write(21,*) 0   ! file outl.<ident>
           endif
           u(i)=rhs(i)
           chi2start=chi2start+(rhs(i)/grpw(i))**2
         enddo
         close(21)  ! file outl.<ident>
         chi2start=chi2start/(nrow-noutl)
         chi2try(1)=chi2start
         if(chi2target.gt.chi2start) then
           print *,'Your target Chi^2/N of ',chi2target
           print *,'Is larger than the null model Chi^2 of ',chi2start
           print *,'after removal of',noutl,' outliers, or',
     &      100.*float(noutl)/float(nrow),'%'
           call bailout(ibr,'Chi^2 target too high')
         endif  
  
         write(4,204) chi2start,0.,noutl,100.*float(noutl)/float(nrow)
         write(6,204) chi2start,0.,noutl,100.*float(noutl)/float(nrow)
204      format(' 1b',f12.3,5(9x,'0'),3x,'1.000',f10.4,'  (removed',i6,
     &      ' outliers or',f6.1,'%)')
c============================================================
      endif     ! BACK TO ALL
c============================================================

         ! broadcast outlier flags to subprocesses
c        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         !write(6,*) myid,' passed barrier R22 t=',secnds(tcpu0)
         call MPI_BCAST(outl,nrow,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
c Compute data fit on root for outlier identification
        
         ! redo computations for minimal damping
         go to 201
      endif             ! end of "klean" loop
c end of end LSQR iterations

c evaluate this solution

c======================================================================
      if(myid.eq.0) then        ! IF ROOT PROCESS
c======================================================================


      ! apply final smoothing v=S*sol in case of implicit smoothing
      if(ksmooth.gt.1)  then   ! smooth sol and store in v
        call smooth(sol,v,ncol,epsnorm,epssmooth,nb,iv,ksmooth,jlast) 
      else
        do i=1,ncol
          v(i)=sol(i)
        enddo  
      endif  

c============================================================
      endif     ! BACK TO ALL
c============================================================
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      !write(6,*) myid,' passed barrier R23 t=',secnds(tcpu0)
      call MPI_BCAST(v,ncol,MPI_REAL,0,MPI_COMM_WORLD,ierr)

c=======================================================================
      if(myid.eq.0) then
c=======================================================================
      do i=nrow+1,nrowd      ! ignore damping part
        u(i)=0.
      enddo  
c============================================================
        else              ! SUB-PROCESSES
c=============================================================
        j=0
        do irow=kdisp(myid)+1,kdisp(myid)+knt(myid)
          u(irow)=0.
          if(outl(irow))then
            j=j+ia(irow)
          else 
            do k=1,ia(irow)
              j=j+1
              u(irow)=u(irow)+a(j)*v(ja(j))
            end do
          endif  
        end do

c============================================================
      endif     ! BACK TO ALL
c============================================================
c       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        !write(6,*) myid,' passed barrier S24 t=',secnds(tcpu0)
        call MPI_GATHERV(u(kdisp(myid)+1),knt(myid),MPI_REAL,w3,knt,
     &        kdisp,MPI_REAL,0,MPI_COMM_WORLD,ierr)

c=======================================================================
      if(myid.eq.0) then
c=======================================================================

      ! compute chi-square on data only (ie ignore damping tail end)
      ! and write residuals to file res.<ident>
      write(fname3,fmt='(i6)') itry+100
      fname3(1:4)='res.'
      fname3(7:7)='.'
      fname3(8:NFN)=ident2(1:NFN-7)
      open(10,file=directory(1:ldr)//fname3)  ! for residuals

c map w3 into u if gatherv was used
      do i=1,nrowd
        u(i)=w3(i)
      enddo  

      chi2=0.
      do i=1,nrow
        if(.not.outl(i))then 
          res=rhs(i)-u(i)
          chi2=chi2+(res/grpw(i))**2    ! Chi^2 (unweighted residuals)
        else
!          print *, 'i,rhs(i),u(i) ',i,rhs(i),u(i)
        endif   
        write(10,*) res 
      enddo
      close(10)  !XX
      chi2try(itry)=chi2/(nrow-noutl)

!      ! write reg.* (regularization part of rhs) 2007/10/23
!      fname3(1:4)='reg.'
!      open(10,file=directory(1:ldr)//fname3)  ! residuals, 2nd part
!      do i=nrow+1,nrowd  ! 2*ncol elements
!        res=rhs(i)-u(i)            
!        write(10,*) res 
!      enddo
!      close(10) 

      ! Evaluate chi2 of this run
      if(itry.eq.2) then
        ! After undamped run: can targeted chi2 be reached?
        chi2low=chi2try(itry) 
        if(chi2low.gt.chi2target.and.chi2target.gt.0.)then
          print *, 'Targeted and minimum chi2:',chi2target,chi2low 
          call bailout(ibr,'STOP: Target chi^2 too low.')
        endif
      else             
        ! After every damped run: is this the last run still
        ! below the targeted chi2, or the first one above?
        ! If yes, save corresponding damping parms ylow, yup. 
        if(chi2target.lt.chi2try(itry)) then
          chi2up=chi2try(itry)
          yup=ytry(itry)
        else
          chi2low=chi2try(itry)
          ylow=ytry(itry)
        endif
        if(jdebug.gt.0) write(13,*) ylow,chi2low,yup,chi2up
      endif  
      
      ! output velocity model solutions solx.<ident2>
      do k=1,3
        if(jssw(k).ne.0) then
          chcor5=chcor(k)(1:5)
          write(fname2,fmt='(i7,".",a5)') itry+100,chcor5
          fname2(1:5)='solx.'
          fname2(14:14)='.'
          fname2(15:NFN)=ident2(1:NFN-14)
c          fname2(15:50)=ident(1:36)
          open(9,file=directory(1:ldr)//fname2)
          write(9,fmt='(i3,4f10.5,1x,a5)') 3,epsnorm,epssmooth,epsratio,
     &        ytry(itry),chcor5
          write(9,*) np
          do i=1,np
            ii=jlast(k)+i
c(bug)      write(9,*) (points(j,i),j=1,3),sol(ii)*dpar(k)
            write(9,*) (points(j,i),j=1,3),v(ii)*dpar(k)
          enddo
          close(9)
        endif  
      enddo  

      ! output full solution sol.<ident2>, in normalized units (chi's)
      write(fname3,fmt='(i6)') itry+100
      fname3(1:4)='sol.'
      fname3(7:7)='.'
      fname3(8:NFN)=ident2(1:NFN-7)
      open(9,file=directory(1:ldr)//fname3) 
      write(9,*) 'Full solution chis. 
     & Next 3 lines: jssw(1:10),dpar(1:10),jlast(1:11)'
      write(9,fmt='(10i12)')   (jssw(k), k=1,10)
      write(9,fmt='(10e12.3)') (dpar(k), k=1,10) 
      write(9,fmt='(11i12)')   (jlast(k),k=1,11)
      do k=1,10  
        do i=jlast(k)+1,jlast(k+1)      
c(bug)    write(9,fmt='(g12.3)')  sol(i)
          write(9,fmt='(g12.3)')  v(i)
        enddo
      enddo
      close(9)

      ! output station/event corrections cor.<ident2> 
      write(fname3,fmt='(i6)') itry+100
      fname3(1:4)='cor.'
      fname3(7:7)='.'
      fname3(8:NFN)=ident2(1:NFN-7)
      open(9,file=directory(1:ldr)//fname3) 
      write(9,*) 'src/rec corrections (7 kinds). 
     & Next 3 lines: jssw(i), dpar(i),n(i); i=4,10'
      write(9,fmt='(7i12)') (jssw(k),k=4,10)
      write(9,fmt='(7e12.3)') (dpar(k),k=4,10) 
      write(9,fmt='(7i12)') (jlast(k+1)-jlast(k),k=4,10)
      write(9,fmt='(i12,a)') ncol, ' ncol total lines in sol(i)'
      write(9,fmt='(i12,a)') ncol-jlast(4),
     &' last (=corrections) lines of sol(i):'
      do k=4,10  
        do i=jlast(k)+1,jlast(k+1)      
c(bug)    write(9,fmt='(g12.3)')  sol(i)*dpar(k)
          write(9,fmt='(g12.3)')  v(i)*dpar(k)
        enddo
      enddo
      close(9)


      ! output diagnostics
      solmax=0.
      solrms=0.
      do i=1,npar                       ! sum over Vp,Vs and Q only
c(bug)  solmax=max(abs(sol(i)),solmax)  ! (keep parameter scaling*dpar)
        solmax=max(abs(v(i)),solmax)  ! (keep parameter scaling*dpar)
        solrms=solrms+sol(i)**2
      enddo
      totrms=solrms
      solrms=sqrt(solrms/npar)
      cormax=0.
      corrms=0.
      do i=npar+1,ncol                  ! sum over all corrections
c(bug)  cormax=max(abs(sol(i)),cormax)
c(bug)  corrms=corrms+sol(i)**2
        cormax=max(abs(v(i)),cormax)
        corrms=corrms+v(i)**2
      enddo
      totrms=totrms+corrms
      if(corrms.gt.0.) corrms=sqrt(corrms/max(1,ncol-npar))
      if(totrms.gt.0.) totrms=sqrt(totrms/ncol)
c     itm=time()        ! nonstandard Fortran
      itm=0
      itm1=0    ! debug fix

      ! Screen and file output: results for this run
      write(6,210) itry,chi2try(itry),solmax,solrms,cormax,corrms,
     &  itm-itm1,ytry(itry),totrms,epsnorm,epssmooth
      write(4,210) itry,chi2try(itry),solmax,solrms,cormax,corrms,
     &  itm-itm1,ytry(itry),totrms,epsnorm,epssmooth
210   format(i3,f12.3,4f10.4,i10,f9.5,f10.4,1x,2g12.3)
      write(7,220) itry,chi2try(itry),solmax,solrms,cormax,corrms,
     &  itm-itm1,ytry(itry)
220   format(i3,' & ',f11.3,4(' $ ',f10.4),' & ',i10,' & ',f9.5,' \\\\')
      write(8,*) solrms,chi2try(itry)
c end of compute chi-square, evaluate solution of this iteration
      
c     estimate next damping value ytry
      if(chi2target.lt.0.) then
        kready=1
        write(6,*) 'End of run for one specific solution'
        goto 700     ! Only one specific solution wanted; done.
      endif  
      if(chi2target.eq.0.and.itry.lt.11) then  ! L-curve: scanning...
        ytry(itry+1)=ytry_lcurve(itry+1)
!!!        ytry(itry+1)=ytry(itry)+.1      ! increment by fixed amount  *** changed, test
!!!        if(itry.eq.2) ytry(3)=0.1                                    *** changed, test
c       write(6,*) 'Scanning next value of L-curve, itry=',itry+1
        goto 700
      else if(chi2target.eq.0.) then    ! L-curve only: done. Break.
        kready=1
        write(6,*) 'Finished computing L-curve only'
        goto 700     ! Only one specific solution wanted; done.
      else if(itry.lt.5) then           ! bisect (chi2target>0) 
        ytry(itry+1)=0.5*(yup+ylow)
c       write(6,*) 'Next try ',itry+1,' by bisection'
        goto 700
      else if(itry.eq.5) then           ! redo ylow with more iters  
        if(itmax.ne.itmax2) then
          write(6,*) 'itmax changed to ',itmax2,' recomputing ylow'
          write(4,*) '(itmax changed to ',itmax2,')'
          itmax=itmax2       ! final trials with more iterations
          ytry(itry+1)=ylow
          goto 700
        else
          ytry(4)=ylow
          chi2try(4)=chi2low
          ytry(5)=yup
          chi2try(5)=chi2up
          write(4,*) '(iterating with starting chi2 brackets:',
     &          chi2low,chi2up,')'
          itry=5
          jzoom=1
          goto 700
        endif
      else if(itry.eq.6) then
        ytry(itry+1)=yup
        write(6,*)  'recomputing yup'
        goto 700
      else  
        dif=abs(chi2target-chi2try(itry))
        if(dif.lt.0.01) then
          kready=1
          write(6,*) 'Exit on dif<0.01'
          goto 700     ! Only one specific solution wanted; done.
        endif  
        if(abs(chi2up-chi2low).lt.0.01) then
          kready=1
          write(6,*) 'Exit on chi2 interval convergence'
          goto 700     ! Only one specific solution wanted; done.
        endif  
        diflast=chi2try(itry-1)-chi2try(itry)
        if(abs(diflast).lt.0.01) then
          write(6,*) 'Exit on chi2 change'
          kready=1
          goto 700     ! Only one specific solution wanted; done.
        endif  
        
c       write(6,*) 'Next try',itry+1,' by Newton interpolation'
        ytry(itry+1)=ylow+(chi2target-chi2low)*(yup-ylow)/(chi2up-
     &        chi2low)
        if(ytry(itry+1).lt.0.0001.or.ytry(itry+1).gt.0.9999) then
          print *, 'WARNING: Lack of convergence'
          kready=1
          goto 700     ! Only one specific solution wanted; done.
        endif  
      endif  
c============================================================
      endif     ! BACK TO ALL
c============================================================

c broadcast the ready flag

c 700   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
700   continue
      !write(6,*) myid,' passed label 700 t=',secnds(tcpu0)
      call MPI_BCAST(kready,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(itmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(kready.eq.0) goto 200  ! next inversion/trial run

      ! END OF INVERSION RUNS

800   continue

c the rest is mostly final model output from root

c=======================================================================
      if(myid.eq.0) then        ! ROOT PROCESS
c=======================================================================
      if(chi2target.eq.0.) goto 900

      
      ! File output:
      ! detailed diagnostics for the final model
      write(4,810) fname2
      write(6,810) fname2
810   format(//,'Final model in file ',a,//,'Parameter diagnostics',
     & ' (unscaled by dpar, statistics over nonzeroes only)',
     & //,'Parameter',15x,'Average',9x,'RMS',9x,'N',4x,'Zeroes')
      write(7,820)
820   format('TABLE OF PARAMETER DIAGNOSTICS',/,'Parameter & ',
     & ' Average & RMS & N & 0 \\\\')
      do k=1,10
        if(jssw(k).ne.0) then
          average=0.
          nzero=0
          rms=0.
          do i=jlast(k)+1,jlast(k+1)
c(bug)      if(abs(sol(i)).lt.1.0e-10) then
            if(abs(v(i)).lt.1.0e-10) then
              nzero=nzero+1
            else  
c(bug)        average=average+sol(i)
c(bug)        rms=rms+sol(i)**2 
              average=average+v(i)
              rms=rms+v(i)**2 
            endif
          enddo
          nn=jlast(k+1)-jlast(k)-nzero
          rms=sqrt(rms/max(1,nn))
          average=average/max(1,nn)
          write(6,830) chcor(k),average*dpar(k),rms*dpar(k),nn,nzero
          write(4,830) chcor(k),average*dpar(k),rms*dpar(k),nn,nzero
830       format(a19,2g12.3,2i10)
          write(7,840) chcor(k),average*dpar(k),rms*dpar(k),nn,nzero
840       format(a19,' & ',g12.3,' & ',g12.3,' & ',i10,' & ',i10,
     &     ' \\\\')
        endif
      enddo

      write(4,850) ngroup
850   format(//,'Number of data groups:',i4,//,'Group',9x,'N',5x,
     &       'Chi^2/N')
      do i=1,ngroup
        write(4,fmt='(i5,i10,g12.3)') i,ngrp(i),chi2grp(i)
      enddo  

      ! Output values of epssmooth and epsnorm for final model
      write(4,890) epsnorm,epssmooth,epsratio,itmax2
890   format(//,'Damping for final model:',/,'epsnorm=',g12.3,
     & /,'epssmooth=',g12.3,/,'epsnorm/epssmooth=',g12.3,/,
     & /,'max nr of iterations for LSQR =',i6)
      
900   print *,'Regular ending of mpisolvetomo'

c============================================================
      endif     ! BACK TO ALL
c============================================================
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
      stop 'Regular stop'

      end


      
      !--------------------------------------------------------------
      
      subroutine adn(ii,jj,nb,iv)

c     add jj to the neighbours of ii
c     nb(),iv() number nodes from 1 (in contrast to qhull)

c     before first call set all nb(i) to 0

c     input: ii,jj
c     output: nb(ii) is number of neighbours of ii
c             iv(j,ii) is j'th vertex connected to ii (j=1,nb(ii))

      include 'includes/solvet.h'
      dimension nb(NPMAX),iv(NBR,NPMAX)

      ni=nb(ii)				! # of neighbours on file
      do i=1,ni
        if(jj.eq.iv(i,ii)) return	! avoid duplications
      end do	
      nb(ii)=nb(ii)+1
      if(nb(ii).gt.NBR) then
        print *,'node ii, nb=',ii,nb(ii)
	print *,'connections:'
	print *,(iv(i,ii),i=1,nb(ii)-1)
	call bailout(ibr,'more than NBR neighbours!')
      endif
      iv(nb(ii),ii)=jj
      
      nb(jj)=nb(jj)+1			! reciprocal neighbour
      if(nb(jj).gt.NBR) then
        print *,'node jj, nb=',jj,nb(jj)
	print *,'connections:'
	print *,(iv(i,jj),i=1,nb(jj)-1)
	call bailout(ibr,'more than NBR neighbours!')
      endif	
      iv(nb(jj),jj)=ii

      return
      end

      subroutine dotp(x,v,n,p)
      dimension x(n),v(n)
      p=0.
      do i=1,n
        p=p+x(i)*v(i)
      enddo
      return
      end

      subroutine randomize(idum,u,m)
      dimension u(m)
      e=0.
      do i=1,m
        u(i)=ran1(idum)
        e=e+u(i)*u(i)
      enddo
      return
      end

      subroutine nullify(u,m)
      dimension u(m)
      do i=1,m
        u(i)=0
      enddo
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

      subroutine aprvol(np,nb,iv,points,tvol,tnab,vscale)
      
      ! compute an approximate tetrahedral volume associated with
      ! each vertex
      ! note: these are not the actual Delauney tetrahedra. Rather,
      ! we assume that each node is in one tetrahedron with facet
      ! length equal to the average of all facets that include the node
      ! We assume this is OK for weighting purposes

      ! input: np - number of grid points
      !        nb,iv (from subroutine adn)
      ! output: tvol(i) for i=1,np - with tvol scaled to average 1
      !         tnab(i) - scaling for nabla square damping
      !         vscale scaled tvol to average 1.

      include 'includes/solvet.h'

      dimension nb(NPMAX),iv(NBR,NPMAX),tvol(NPMAX),tnab(NPMAX)
      dimension points(3,NPMAX)

      if(np.gt.NPMAX) call bailout(ibr,'np too large in aprvol')
      vscale=0.
      do i=1,np
        s=0.
        x=points(1,i)
        y=points(2,i)
        z=points(3,i)
        ! find average facet length
        do j=1,nb(i)
          if(j.gt.NBR) call bailout(ibr,'j error')
          k=iv(j,i)
          if(k.gt.np.or.k.lt.1) call bailout(ibr,'k error')
          xx=points(1,k)
          yy=points(2,k)
          zz=points(3,k)
          d=sqrt((x-xx)**2+(y-yy)**2+(z-zz)**2)
          s=s+d
        enddo
        ! compute tetrahedral volume
        s=s/nb(i)
        tvol(i)=0.11785*s**3    ! dV= sqrt(2)/12 * s^3 if equilateral
        tnab(i)=s*s             ! store squared facet length
        vscale=vscale+tvol(i)   ! total volume of the model
      enddo
      vscale=vscale/np          ! average volume of tetrahedra
      do i=1,np
        tvol(i)=tvol(i)/vscale       ! scale with tetrahedral volume dV
        tnab(i)=tvol(i)/tnab(i)      ! dV * 1/s^2 for 2nd derivative
      enddo  
c     write(6,*) 'return from aprvol with vscale=',vscale

cccc  UNDO the effect of tnab (volume regu for smoothing) -- set to 1 for now
cccc  2007/10/18 Karin
      do i=1,np 
        tnab(i) = 1.
      enddo          
cccc  END UNDO the effect of tnab (volume regu for smoothing) -- set to 1 for now

      return
      end
      

      subroutine smooth(v,w,ncol,epsnorm,epssmooth,nb,iv,ksmooth,jlast)

      ! smoothing routine: w_i = sum_j S_ij v_j
      ! where S smoothes over the nb(i) neighbours of node i

      dimension v(ncol),w(ncol)
      include 'includes/solvet.h'
      dimension nb(NPMAX),iv(NBR,NPMAX),jlast(11)

      ! input: vector v, dimension ncol (model + corrections)
      ! output: vector w, dimension ncol, = v smoothed over model 
      !         parameters but not corrections

      ! smoothing according to ksmooth:
      ! 2: all weights equal to 1/(nb(i)+1) (nb=nr of neighbours of i)
      ! 3: w=epssmooth for i, (1-epssmooth)/nb(i) for neighbours
      ! 4: distance weighting (not yet implemented): w=epssmooth
      !    for node i, (1-epssmooth)/(Dij*sum 1/Dij) for neighbour j
      !    where Dij is distance between nodes i and j.
      ! In all cases the sum of the weights x is 1

      ! Neighbours:
      ! nb(k) is number of neighbours of k
      ! iv(j,k) is j'th vertex connected to k (j=1,nb(k))
      
      do ktype=1,3      ! loop over Vp,Vs,Q
        i0=jlast(ktype)
        i1=i0+1
        i2=jlast(ktype+1)
        do i=i1,i2                      ! matrix row number i
          k=i-jlast(ktype)              ! model grid number
          Sij=epssmooth   
          if(ksmooth.eq.2) Sij=1.0/(nb(k)+1.)  ! all weights equal
          w(i)=Sij*v(i)                 ! target node
          if(ksmooth.eq.3) Sij=(1.-epssmooth)/nb(k)
          do j=1,nb(k)                  ! loop over neighbours j
            kk=iv(j,k)+i0               ! neighbour column number
            w(i)=w(i)+Sij*v(kk)
          enddo
        enddo
      enddo

      do i=jlast(4)+1,jlast(11)
        w(i)=v(i)               ! do not smooth the corrections
      enddo

      return
      end


      subroutine tsmooth(w,v,ncol,epsnorm,epssmooth,nb,iv,ksmooth,jlast)

      ! transpose smoothing routine
      ! v_i = sum_j S_ji w_j
      ! where S smoothes over the nb(i) neighbours of node i

      dimension w(ncol),v(ncol)
      include 'includes/solvet.h'
      dimension nb(NPMAX),iv(NBR,NPMAX),jlast(11)

      ! input: vector w, dimension ncol (model + corrections)
      ! output: vector v, dimension ncol, = w smoothed with
      !         the transpose of the smoothing operator in smooth().

      ! smoothing in smooth() according to ksmooth:
      ! 2: all weights equal to 1/(nb(i)+1) (nb=nr of neighbours of i)
      ! 3: w=epssmooth for i, (1-epssmooth)/nb(i) for neighbours
      ! 4: distance weighting (not yet implemented): w=epssmooth
      !    for node i, (1-epssmooth)/(Dij*sum 1/Dij) for neighbour j
      !    where Dij is distance between nodes i and j.
      ! In all cases the sum of the weights x is 1

      ! Neighbours:
      ! nb(k) is number of neighbours of k
      ! iv(j,k) is j'th vertex connected to k (j=1,nb(k))

      do i=1,jlast(4)
        v(i)=0.               
      enddo
      do i=jlast(4)+1,jlast(11)
        v(i)=w(i)             ! corrections not smoothed
      enddo  
      
      do ktype=1,3      ! loop over Vp,Vs,Q
        i0=jlast(ktype)
        i1=i0+1
        i2=jlast(ktype+1)
        do i=i1,i2              ! row of matrix S
          k=i-i0                ! model grid number
          Sij=epssmooth   
          if(ksmooth.eq.2) Sij=1.0/(nb(k)+1.)       ! equal weights 
          v(i)=v(i)+Sij*w(i)
          if(ksmooth.eq.3) Sij=(1.-epssmooth)/nb(k)
          do j=1,nb(k)          ! loop over neighbours of node k
            kj=iv(j,k)          ! model node # of the neighbour
            kk=kj+i0            ! matrix row number
            if(kk.eq.0) call bailout(ibr,'kk=0 in tsmooth')
            v(kk)=v(kk)+Sij*w(i)
          enddo
        enddo
      enddo

      return
      end

      subroutine bailout(ibr,txt)
      character*(*) txt
      include "mpif.h"      
c      include 'includes/mpif.h'
      write(6,*) 'Bailout!!!'
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)  ! get process number
      write(6,*) 'Bailing out from process ',myid
      write(6,*) txt
c     call MPI_FINALIZE(ierr)
c     stop 'Stop in bailout'
      ibr=1
      return
      end

c     subroutine check(ierr,line)
c     include 'includes/mpif.h'
c     character*(*) line
c     if(ierr.ne.MPI_SUCCESS) then
c       print *, 'MPI error ',ierr
c       print *,line
c       call MPI_FINALIZE(ierr)
c       stop 'Stop in check'
c     endif     
c     return
c     end

      subroutine normlz(n,x,sreturn)
      dimension x(n)
      double precision s,ss,xx    ! against overflow of squares
      
      s=0.d0
      umin=1.0e30
      umax=-1.0e30
      do i=1,n 
        xx=x(i)
        s=s+xx**2 
      end do  
      s=sqrt(s)
      sreturn=s         ! single precision returned
      ss=1./s 
      do i=1,n 
        x(i)=x(i)*ss 
      end do  
      return  
      end 

      function length(ch)
      character*(*) ch
      integer ich
      ich = 1
      do while ((ch(ich:ich) .ne. '  ') .and. (ich .le. len(ch)))    
         ich = ich + 1
      end do
      length = ich - 1
      return
      end
