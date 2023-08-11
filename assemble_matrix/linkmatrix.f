      program linkmatrix

      ! compile: g77 linkmatrix

      ! History:
      ! Created GN 2/18/07

      ! linkmatrix reads a master matrix file, combines it with its
      ! "linked" file into a new matrix. For example, if run on
      ! the PP file, it will find the corrsponding P file, select
      ! common paths, subtract the row in the P matrix from that
      ! in the PP matrix and write a new file.

      ! Data in the master file are already *differences*. The
      ! data in the linked file are ignored.
      
      ! header files for Delauney triangulation
      include '../includes/nn.param'   ! only needed for np_max

      ! matrix rows 
      ! allow for MMAX nonzero elements in one row, MMAX data,
      ! NDIM unknowns (models parameters + corrections)
      parameter(MMAX=200000)
      parameter(NFR=20)
      dimension asparse(MMAX),ja(MMAX)
      dimension asparsel(MMAX),jal(MMAX)
      dimension aa(3*np_max)  !!! one row

      ! other
      dimension kountpar(3),nssw(3),nsswl(3),kountparl(3)

      dimension period(NFR),periodl(NFR)

      ! various character variables 
      character*72 fname,fname2,ident,fname3,none
      character*16 stationcode,stationcodel
      character*8  phase,networkcode,netw
      character*8  phasel,networkcodel,netwl
      character*30 dataf,datafl,dataf2
      character*3  comp,compl

      ! Normally jdebug=0. Set -1 if input and output files are ascii
      ! instead of binary. Set 1 for some debugging output.
      jdebug=0

      write(6,fmt='(a)') 'Running linkmatrix ...' 
      none='None'

      ! Read master matrix file name 
      print *; print *,'Give master matrix file name (eg matrixT.PP)'
      read(5,fmt='(a)') fname

      kdatatype=0
      if(fname(7:7).eq.'T') kdatatype=1
      if(fname(7:7).eq.'A') kdatatype=2
      if(kdatatype.eq.0) stop 'Data type not recognized.'
      
      print *
      if(jdebug.eq.0) then      ! unformatted
        print *,'Matrix files must be unformatted'
        print *,'Opening master file ',fname
        open(11,file=fname,status='old',form='unformatted')
        read(11,iostat=ios) linked,dataf2   ! linked file if linked=1
        if(ios.ne.0) stop 'Error in line 1'
        if(linked.ne.1) stop 'This appears not to have a linked file'
        print *,'File linked to: ',dataf2
        read(11,iostat=ios) np,nssw     ! nr of grid points, par types
        if(ios.ne.0) stop 'Error in line 2'
        read(11,iostat=ios) nband,(period(i),i=1,nband)
        if(ios.ne.0) stop 'Error in line 3 for period'
      else  
        print *,'Matrix files must be formatted, jdebug=',jdebug
        print *,'Opening master file ',fname
        open(11,file=fname,status='old')
        read(11,fmt='(i2,1x,a)',iostat=ios) linked,dataf2 
        if(ios.ne.0) stop 'Error in line 1 for linked file'
        if(linked.ne.1) stop 'This appears not to have a linked file'
        print *,'File linked to: ',dataf2
        read(11,*,iostat=ios) np,nssw   ! nr of grid points, par types
        if(ios.ne.0) stop 'Error in line 2 for np, nssw'
        read(11,*,iostat=ios) nband,(period(i),i=1,nband)
        if(ios.ne.0) stop 'Error in line 3 for period'
      endif
      if(np.gt.np_max) stop 'Increase np_max'

      ! Nown ope linked file
      fname2='matrixT.'//dataf2
      if(kdatatype.eq.2) fname2='matrixA.'//dataf2
      print *,'Linked file name: ',fname2
      if(jdebug.eq.0) then      ! unformatted
        print *,'Opening linked file ',fname2
        open(21,file=fname2,status='old',form='unformatted')
        read(21,iostat=ios) linkl,datafl   ! linked file if linked=1
        if(ios.ne.0) stop 'Error in line 1'
        if(linkl.ne.0) then  ! handle case PPP-PP with a warning
          print *,'WARNING: linked file is itself linked?'
          print *,'linked=',linkl,', to: ',datafl
          print *,'Type <1> to continue'
          read *,kontinue
          if(kontinue.ne.1) stop
        endif
        read(21,iostat=ios) npl,nsswl     ! nr of grid points, par types
        if(ios.ne.0) stop 'Error in line 2'
        read(21,iostat=ios) nbandl,(periodl(i),i=1,nband)
        if(ios.ne.0) stop 'Error in line 3 for period'
      else  
        print *,'Matrix files must be formatted, jdebug=',jdebug
        print *,'Opening matrix file ',fname2
        open(21,file=fname2,status='old')
        read(21,fmt='(i2,1x,a)',iostat=ios) linked,dataf2 
        if(ios.ne.0) stop 'Error in line 1 for linked file'
        if(linked.ne.0) then
          print *,'WARNING: linked file is itself linked?'
          print *,'linked=',linkl,', to: ',datafl
          print *,'Type <1> to continue'
          read *,kontinue
          if(kontinue.ne.1) stop
        endif
        read(21,*,iostat=ios) npl,nsswl   ! nr of grid points, par types
        if(ios.ne.0) stop 'Error in line 2 for np, nssw'
        read(21,*,iostat=ios) nbandl,(periodl(i),i=1,nband)
        if(ios.ne.0) stop 'Error in line 3 for period'
      endif

      ! check if linkage is correct
      if(np.ne.npl) stop 'np incompatible between linked files'
      if(nbandl.gt.nband) stop 'nband for linked file larger?'
      do i=1,3
        if(nssw(i).ne.nsswl(i)) stop 'nssw incompatible between files'
      enddo
      do i=1,nband
        if(abs(period(i)-periodl(i)).gt.0.1) 
     &     stop 'Frequency bands incomaptible between linked files'
      enddo

      print *,'Give ident for the output matrix file (e.g. PP-P):'
      read(5,'(a)') ident
      fname3=fname(1:8)//ident
      if(jdebug.eq.0) then      ! unformatted
        open(31,file=fname3,status='new',form='unformatted')
        write(31) 0,none   ! linked file if linked=1
        write(31) np,nssw     ! nr of grid points, par types
        write(31) nband,(period(i),i=1,nband)
      else  
        open(31,file=fname3,status='new',form='formatted')
        write(31,*) '0 ',none   ! linked file if linked=1
        write(31,*) np,nssw     ! nr of grid points, par types
        write(31,*) nband,(period(i),i=1,nband)
      endif
      open(4,file='out.linkmatrix.'//ident)
      write(4,'(2a)') 'Master file: ',fname
      write(4,'(2a)') 'Linked file: ',fname2
      write(4,'(2a)') 'Difference file: ',fname3
      write(4,10)
10    format(//,9x,'i     ilink      band     event  station',
     &      '      kount  accur(%)')

      if(jdebug.eq.1) write(13,*) 'Now reading ',fname
      write(6,fmt='(a,$)') 'Progress: '   

      ntot=0
      ntotl=0
      npall=0
      do i=1,3
        if(nssw(i).gt.0) npall=npall+np
      enddo  

      ! read one master matrix row, rhs element & corrections
      ! data are uniquely identified by combination of ndatum and iband
20    if(jdebug.eq.0) then                    ! unformatted  
        read(11,end=200,iostat=ios) ndatum,kount,iband,dobs,dtheor,s,
     &    ecorr,corcrust,corelev,tstar,qray,ievt,kluster,stationcode,
     &    networkcode,comp,slowness,raz,rdel,angle0,vsrc,krtyp
        if(ios.ne.0) stop 'Error in master file - header'
        if(kount.gt.MMAX) stop 'kount>MMAX, increase dimensions'
        read(11,end=200,iostat=ios) kountpar
        if(ios.ne.0) stop 'Error in master file - kountpar'
        read(11,end=200,iostat=ios) (ja(i),i=1,kount)
        if(ios.ne.0) stop 'Error in master file - ja'
        read(11,end=200,iostat=ios) (asparse(i),i=1,kount)
        if(ios.ne.0) stop 'Error in master file - asparse'
      else  ! formatted
        read(11,*,end=200,iostat=ios) ndatum,kount,iband,dobs,dtheor,s,
     &  ecorr,corcrust,corelev,tstar,qray,ievt,kluster,stationcode,
     &    networkcode,comp,slowness,raz,rdel,angle0,vsrc,krtyp
        if(jdebug.gt.0) write(13,*) 'Master file:',
     &   ndatum,kount,iband,dobs,dtheor,s,
     &   ecorr,corcrust,corelev,tstar,qray,ievt,kluster,stationcode,
     &   networkcode,comp,slowness,raz,rdel,angle0,vsrc,krtyp
        if(ios.ne.0) stop 'Error in master file - header'
        if(kount.gt.MMAX) stop 'kount>MMAX, increase dimensions'
        read(11,*,end=200,iostat=ios) kountpar
        if(ios.ne.0) stop 'Error in master file - kountpar'
        if(jdebug.gt.0) write(13,*) 'kountpar=',kountpar
        read(11,*,end=200,iostat=ios) (ja(i),i=1,kount)
        if(ios.ne.0) stop 'Error in master file - ja'
        if(jdebug.gt.0) write(13,*) 'ja=',(ja(i),i=1,6),' etc'
        read(11,*,end=200,iostat=ios) (asparse(i),i=1,kount)
        if(ios.ne.0) stop 'Error in master file - asparse'
        if(jdebug.gt.0) write(13,*) 'a=',(asparse(i),i=1,6),' etc'
      endif  

      ntot=ntot+1

      
      ! read next datum in linked file
40    if(jdebug.eq.0) then                    ! unformatted  
        read(21,end=100,iostat=ios) ndl,kountl,ibandl,dobsl,dtheorl,sl,
     &    ecorrl,corcrustl,corelevl,tstarl,qrayl,ievtl,klusterl,
     &    stationcodel,networkcodel,compl,slownessl,razl,rdell,
     &    angle0l,vsrcl,krtypl
        if(ios.ne.0) stop 'Error in linked file - header'
        read(21,end=100,iostat=ios) kountparl
        if(ios.ne.0) stop 'Error in linked file - kountpar'
        read(21,end=100,iostat=ios) (jal(i),i=1,kountl)
        if(ios.ne.0) stop 'Error in linked file - ja'
        read(21,end=100,iostat=ios) (asparsel(i),i=1,kountl)
        if(ios.ne.0) stop 'Error in linked file - asparse'
      else  ! formatted
        read(21,*,end=100,iostat=ios) ndl,kountl,ibandl,dobsl,dtheorl,
     &    sl,ecorrl,corcrustl,corelevl,tstarl,qrayl,ievtl,klusterl,
     &    stationcodel,networkcodel,compl,slownessl,razl,rdell,
     &    angle0l,vsrcl,krtypl
        if(ios.ne.0) stop 'Error in linked file - header'
        if(jdebug.gt.0) write(13,*) 'Linked file:',
     &    ndl,kountl,ibandl,dobsl,dtheorl,sl,
     &    ecorrl,corcrustl,corelevl,tstarl,qrayl,ievtl,klusterl,
     &    stationcodel,networkcodel,compl,slownessl,razl,rdell,
     &    angle0l,vsrcl,krtypl
        read(21,*,end=100,iostat=ios) kountparl
        if(ios.ne.0) stop 'Error in linked file - kountpar'
        read(21,*,end=100,iostat=ios) (jal(i),i=1,kountl)
        if(ios.ne.0) stop 'Error in linked file - ja'
        read(21,*,end=100,iostat=ios) (asparsel(i),i=1,kountl)
        if(ios.ne.0) stop 'Error in linked file - asparse'
      endif  

      ntotl=ntotl+1

      ! skip?
      if(ibandl.ne.iband) goto 40
      if(stationcodel.ne.stationcode) goto 40
      if(ievtl.ne.ievt) goto 40

      ! accept if event, station and frequency band are the same
      do i=1,npall
        aa(i)=0.
      enddo  
      ! subtract the linked rows from the master rows
      do i=1,kount
        j=ja(i)
        aa(j)=aa(j)+asparse(i)
        sum=sum+asparse(i)
      enddo
      do i=1,kountl
        j=jal(i)
        aa(j)=aa(j)-asparsel(i)
      enddo
      kount=0
      ! write out nonzero elements
      sum=0.
      do i=1,npall
        if(aa(i).ne.0.) then
          kount=kount+1
          ja(kount)=i
          asparse(kount)=aa(i)
          sum=sum+asparse(kount)
        endif
      enddo  
      tdif=dtheor-dtheorl     ! theoretical travel time difference
      if(jdebug.eq.0) then  
        write(31,iostat=ios) ndatum,kount,iband,dobs,tdif,s,
     &    ecorr,corcrust,corelev,tstar,qray,ievt,kluster,stationcode,
     &    networkcode,comp,slowness,raz,rdel,angle0,vsrc,krtyp
        write(31) kountpar
        write(31) (ja(i),i=1,kount)
        write(31) (asparse(i),i=1,kount)
      else  
        write(31,60,iostat=ios) ndatum,kount,iband,dobs,tdif,s,
     &    ecorr,corcrust,corelev,tstar,qray,ievt,kluster,stationcode,
     &    networkcode,comp,slowness,raz,rdel,angle0,vsrc,krtyp
60        format(2i8,i4,2f10.2,5f7.2,f10.1,i8,i3,1x,a16,1x,a8,1x,a3,
     &           f10.1,3f10.5,f8.3,i2,2f10.2,i8)
        write(31,*) kountpar
        write(31,*) (ja(i),i=1,kount)
        write(31,*) (asparse(i),i=1,kount)
      endif

      ! write to out file
      acc=100.*abs(sum/tdif)
      write(4,80) ndatum,ndl,iband,ievt,stationcode(1:8),kount,acc
80    format(4i10,2x,a8,i10,f10.1)

      goto 20

100   print *,'End of linked file reached before end of'
      print *,'master file - some data missing or out of'
      print *,'order in the linked file.'
      print *,'Nr of data read from master:',ntot
      print *,'Nr of data read from linked file:',ntotl
      print *,'Problem occurred when looking for event',ievt
      print *,'Station ',stationcode
      write(4,*)'End of linked file reached before end of'
      write(4,*) 'master file - some data missing or out of'
      write(4,*) 'order in the linked file.'
      write(4,*) 'Nr of data read from master:',ntot
      write(4,*) 'Nr of data read from linked file:',ntotl
      write(4,*) 'Problem occurred when looking for event',ievt
      write(4,*) 'Station ',stationcode
      stop 'ERROR'

200   continue
      print *,'End of master file reached.'
      print *,'Nr of data read from master:',ntot
      print *,'Nr of data read from linked file:',ntotl
      write(4,*)
      write(4,*) 'Nr of data read from master:',ntot
      write(4,*) 'Nr of data read from linked file:',ntotl
      stop 'Normal end of linkmatrix'

      end
