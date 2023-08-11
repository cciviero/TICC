      ! compile: g77 raydata.f raytracesubs.f azdel.f

      program raydata

      ! Alpha version. Only limited testing.
      ! Changes 
      ! August 2006: removed relev from elevation corrections when kd<5
      !              added stop for sources above local surface
      !              corrected ptlat/lon for surface source in PPP etc.
      !              changed azimuth in givloc call for deltas > 180
      ! Sept 2006: Depth of source shallower than crust1.0 surface is replaced 
      !            with crust1.0 elev; a list of shallow sources is in src.shallow
      !            Get max. xcorr length in each band and output to xclm.
      !            Output warning message of imprecise hessian against analytical expression
      ! Dec 2006: Replaced Kennetts ellipticity routine with getelcor. 
      ! 2010/08/12 K.S. change name of common block from /mod/ to /modl/
      !  since old name was a fortran standard violation and ifort would
      ! not compile. Also changed in raytracesubs.f
      ! raydata.f(237): error #8038: A referenced intrinsic procedure 
      ! can not have the same name as a common block.   [MOD]
      
      ! traces rays and produces input file for BD kernel
      ! computations, for each source-receiver pair

      ! input:
      ! [3] the raw src/rcv/dT/dA specification file <dataf>

      ! output:
      ! [2] raydata.<dataf>
      ! [4] raydata.out.<dataf>

      ! the program is most efficient if the data are
      ! sorted in segments of lines with equal source depths

      parameter (ityp=36,nla=180,nlo=360)                        ! for rcrust1 arrays
      parameter (NFR=10)                                        ! nr of freq bands

      dimension rseg(20),ktseg(20),kdwn(20)                     ! ray segments
      dimension nfreq(NFR),bndomega(200,NFR),dotm(200,NFR),rmsb(NFR)  ! spectral bands
      dimension tobs(NFR),nbt(NFR),tsig(NFR),corcoeft(NFR), tb_mfi(NFR) ! time data
      dimension xcl(NFR),xcla(NFR),xclm(NFR)                       ! window lengths
      dimension aobs(NFR),nba(NFR),asig(NFR),corcoefa(NFR)          ! amplitude data
      dimension ystart(4),y(4,9000),ytbl(4,4000),jsg(9000)  ! ray vector
      dimension legend(190)                         ! ray leg indexes
      dimension a1(10),a2(10),d1(10),d2(10),trtime(10),indx(10)   ! ray times/dist
      dimension hmf(2,9000),hmb(2,9000),rayvel(9000),rayq(9000) ! Hessian
      dimension ypq(10,9000)                                    ! idem
      dimension rcr2(16),vcr2(16), ecr2(16)                               ! crust1.0
      dimension rref(50),vref(50)                               ! Reference crust

      logical abrt,rythr
      character*72 fname,tablef,chray,line,datacomment
      character*30 dataf,dataf1,dataf2
      character*16 stationcode
      character*8 phase,netw
      character*2 ctype(ityp),types,atype(nlo)
      character*1 yn,star
      character*3 comp

      ! background model
      common /modl/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &     Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      ! the following are filled by subroutine rcrust1:
      common/crust1/amapvp(9,nlo,nla),amaprho(9,nlo,nla),
     &              amapvs(9,nlo,nla),amapthi(8,nlo,nla),
     &              amapele(nlo,nla), am_crust(43621,3),
     &              am_sed(43621,3), am_moho(43621,3)                    ! model crust1.0
      common/glbgrd/types(nlo,nla)                      ! model crust1.0

      data r2d/57.295779513082/,pi/3.14159265359/
      data halfpi/1.570796326795/

      jdebug=0

      ! set some constants
      h=15.     ! integration step size (must equal h in mkrtable)
      itmax=20  ! max nr of iterations for ray "shooting"

c     ! open ellipticity correction table
c     inquire(file='elcordir.tbl',exist=rythr)
c     if(.not. rythr) then
c       print *,'ERROR: cannot find file elcordir.tbl'
c       stop
c     endif  
c     open(21,file='elcordir.tbl',access='direct',recl=80,
c    &      form='formatted')

      ! read crust1.0 data
c      print *,'CRUST1.0'
      print *, fname
      call rcrust1
      print *,'Give model file:'
      read(5,fmt='(a)') fname
      inquire(file=fname,exist=rythr)
      if(.not. rythr) then
          print *, CHAR(27) //'[0;41m[ERROR]'//CHAR(27)//'[0m'
        print *,'cannot find: ',fname
        stop
      else
        print *,'Reading...',fname
      endif  
      call model(fname)
      
      print *,'-------------'
      write(6,*) 'Radius of Earth model:',r(n)
      if(abs(r(n)-6371.).gt.0.1) write(6,*) 'WARNING! Model radius',
     & ' deviates from 6371 - crustal corrections are not OK'
      write(6,*) 'Radius of CMB:',r(noc)
      write(6,*) 'Radius of ICB:',r(nic)
      write(6,*) 'Radius of Moho:',r(nmoh)

      print *,'==========================='
      print *,'Give max number of arrivals to include in computations'
      print *,'of the banana-doughnut kernels, and'
      print *,'their max difference (in sec) w.r.t. first arrival:'
      read *,maxarr,tdifmax
      print *,'Max arrival: ', maxarr
      print *,'Max T diff : ', tdifmax

      print *,'==========================='
      ! open the raw data file for input
      print *,'Give data file name IN THIS DIRECTORY (eg Pshallow):'
      read(5,fmt='(a)') dataf
      inquire(file=dataf,exist=rythr)
      if(.not. rythr) then
          print *, CHAR(27) //'[0;41m[ERROR]'//CHAR(27)//'[0m'
        print *,'cannot find: ',dataf
        stop
      else
        print *,'Opening...',dataf
      endif  
      open(3,file=dataf)
      print *,'-------------'

      ! read filename & twin file name, if any
      read(3,fmt='(a)') dataf1
      if(dataf1.ne.dataf) stop 'Raw data file name does not match'
      read(3,fmt='(a)') dataf2
      
      print *,'Check the data groups...'
      if(dataf2.ne.'None') then
        print *,'Data group: ',dataf1
        print *,'Twinned with: ',dataf2
        inquire(file=dataf2,exist=rythr)
        if(.not. rythr) then
          print *,'WARNING - cannot find twin file: ',dataf2
          print *,'Should I continue (y/n)?'
          read(5,fmt='(a1)') yn
          if(yn.eq.'n') stop
        endif  
      else
          print *, 'Single group (no twin)!'
      endif
      print *,'==========================='

      print *, 'Skip the following lines:'
      ! skip comment lines starting with #
      read(3,fmt='(a72)') line
      do while (line(1:1).eq.'#')
        write(6,fmt='(a72)') line
        read(3,fmt='(a72)') line
      enddo
      backspace(3)
      print *,'==========================='

      ! open the processed datafile for output
      tablef='raydata.'//dataf
      open(2,file=tablef)
      tablef='raydata.out.'//dataf
      open(4,file=tablef)
      write(2,fmt='(a)') fname
      write(2,fmt='(i3,f10.2)') maxarr,tdifmax
      write(2,fmt='(a)') dataf1
      write(2,fmt='(a)') dataf2
      write(4,'(2a)') 'Data group: ',dataf1
      write(4,'(2a)') ' Twinned with: ',dataf2
      write(4,'(2a)') 'Rays computed for model ',fname
      write(4,'(2a)') 'Input from: ',dataf
      write(4,'(2a)') 'Output to: ','raydata.'//dataf

      ! read ray info from the raw data file, write it to the processed file
      call rdray(3,2,chray,phase,rseg,ktseg,kdwn,nseg)
      write(4,'(5a)') 'Phase: ',phase,' (',chray,')'
      
      ! KASRA
      if (phase.eq.'Pdiff') then
          write(*,*) 'Selected phase is...', phase
          write(*,*) 'p is set to (constant) ', 4.439*180.0/pi
      endif

      print *,'==========================='
      print *, 'Read the filters....'
      ! read band info from the raw data file, write it to the processed file
      call rdband(3,2,bndomega,dotm,nfreq,nband)
      ! initialize xclm - maximum xcorr length for each band
      do ib=1,nband
        xclm(ib)=1. ! not 0 to avoid numerical difficulty in raymatrix
      enddo

      print *,'==========================='
      sdepold=-9999.
      ndata=0
      itm=time()        ! g77 - may need change on other compilers
      write(6,*) 'Now processing data for CRUST1 + Indian Ocean:'
c     write(6,fmt='(i8," to ",i8,2x,a26)') ndata,ndata+100,ctime(itm)
      open(7,file='src.shallow.'//dataf) ! sources shallower than crust1.0 surface
      write(7,*) 'Sources shallower than crust1.0 surface:'
      write(7,*) '  idate  srclat   srclon    srcR crust1.0R    dt0',
     & ' kd'
      jshallow=0 ! jshallow counts sources above crust1.0 surface

10    read(3,*,end=900,iostat=ios) idate,iotime,ievt,kluster,
     &     stationcode,netw,comp,
     &     slat,slon,sdep,rlat,rlon,relev,nobst,nobsa,kpole
      if(ios.ne.0) then
        print *, CHAR(27) //'[0;41m[HEADER ERROR]'//CHAR(27)//'[0m'
        write(6,*) 'reading idate=',
     &     idate,', ievt=',
     &     ievt,', station: ', stationcode, netw
        stop 'Data file error in header line'
      endif  
      read(3,*,iostat=ios) kunit,rms0,(rmsb(i),i=1,nband)
      if(ios.ne.0) then
        print *, CHAR(27) //'[0;41m[NOISE ERROR]'//CHAR(27)//'[0m'
        write(6,*) 'reading idate=',idate,', ievt=',
     &     ievt,', station: ',stationcode
        stop 'Data file error in rms line'
      endif  
      do i=1,nobst
        read(3,*,iostat=ios) tobs(i),tsig(i),corcoeft(i),nbt(i),xcl(i),
     &   tb_mfi(i)
        if(nbt(i).gt.0) xclm(nbt(i))=max(xclm(nbt(i)),xcl(i))
        if(ios.ne.0) then
          write(6,*) 'T data error, reading idate=',idate,', ievt=',
     &     ievt,', station: ',stationcode,' tobs nr',i
          stop 'Data file error in tobs line'
        endif  
      enddo  
      do i=1,nobsa
        read(3,*,iostat=ios) aobs(i),asig(i),corcoefa(i),nba(i),xcla(i)
        if(nba(i).gt.0) xclm(nba(i))=max(xclm(nba(i)),xcla(i))
        if(ios.ne.0) then
          write(6,*) 'A data error, reading idate=',idate,', ievt=',
     &     ievt,', station: ',stationcode,' aobs nr',i
          stop 'Data file error in aobs line'
        endif  
      enddo  

      if(jdebug.gt.0) write(13,*) 'datum:',idate,ievt,' ',stationcode,
     &slat,slon,sdep,rlat,rlon,relev

      ! recompute table of i vs phi (angle vs. distance) if new depth
      if(abs(sdep-sdepold).gt.0.01) then
        ang1=90.0+(kdwn(1)-0.5)*0.2
        ang2=179.99
        if(kdwn(1).eq.0) then
          ang2=ang1
          ang1=0.
        endif  
        rseg(1)=r(n)-sdep
        rsrc=rseg(1)
        call velo1(rsrc,iqsrc,1-kdwn(1),ktseg(1),csrc,dc,d2c)
        call mkrtable(rseg,ktseg,kdwn,nseg,ang1,ang2,ytbl,ntbl)
        if(jdebug.gt.0) write(13,*) 'Ray table',ang1,ang2,ntbl,sdep,
     &          sdepold
        sdepold=sdep
      endif

      ! compute epicentral distance del. Note stadis needs geographical
      ! coordinates
      scolat=90.-slat
      rcolat=90.-rlat
      call stadis(scolat,slon,rcolat,rlon,del,az,t0,p0,t2,p1)
      del=360*((kpole+1)/2)+2*(0.5-mod(kpole,2))*del
      if(jdebug.gt.0) write(13,*) 'delta, src azimuth:',del,az
      !write(*,*) 'delta, src azimuth:',del,az
      rdel=del/r2d
      raz=az/r2d

      ! compute ellipticity correction 
c     ! (ADD ecorr to the spherical value to get elliptical Earth time)
c     ! first transform colatitudes to geocentric, radians
c     ! Temporary: do both ellipticity corrections (ellcor, getelcor,
c     ! though getelcor must wait until we know the ray parameter).
      scolatr=geocen((90.-slat)/r2d)
      rcolatr=geocen((90.-rlat)/r2d)    ! needed for getelcor
      slonr=slon/r2d
      rlonr=rlon/r2d
c     call ellref(scolatr)
c     call ellcor(phase,del,sdep,scolatr,az/r2d,ecorr,abrt)
c     if(abrt) ecorr=-99.               ! if no ellipticity correction available

      ! Kasra
      ! Ellipticity correction using BLN Kennetts method (for Pdiff)
      ! az: inside ellip_blnk it will be converted to radians
      call ellip_blnk(phase, del, sdep, scolat, az, tcor)
      ! write(*,*) phase, del, sdep, scolat, az, tcor
      ! END Kasra
      
      ! find all ray arrivals at distance del within tdifmax of first
      call findarr(del,tdifmax,ytbl,ntbl,a1,a2,d1,d2,trtime,narr)
      if(jdebug.eq.1) then
        write(13,*) 'findarr found the following brackets, del=',del
        do i=1,narr
          write(13,fmt='(i2,4f8.2)') i,a1(i),a2(i),d1(i),d2(i)
        enddo
      endif  

      narr=min(maxarr,narr)     ! limit number of arrivals
      ! add del,az and narr to the data line
      ndata=ndata+1
      itm=time()                ! works with g77 but is not standard
      if(mod(ndata,100).eq.0) 
     &  write(6,fmt='(i8," to ",i8,2x,a26)') ndata,ndata+100,ctime(itm)
c      write(6,25) ndata,ievt,stationcode,comp
c25    format(i6," Event nr",i8,1x,"in ",a16,1x,a3)
      write(2,20) idate,iotime,ievt,kluster,stationcode,netw,comp,
     &      slat,slon,sdep,
     &      rlat,rlon,relev,rdel,raz,del,az,narr
      write(2,22) kunit,rms0,(rmsb(i),i=1,nband)
      write(2,*) nobst
      do i=1,nobst
        write(2,21) tobs(i),tsig(i),corcoeft(i),nbt(i),xcl(i)
      enddo  
      write(2,*) nobsa
      do i=1,nobsa
        write(2,21) aobs(i),asig(i),corcoefa(i),nba(i),xcla(i)
      enddo  
20    format(3i8,i4,1x,a16,1x,a8,1x,a3,2f9.3,f7.1,2f9.3,f7.3,
     &       2f9.5,2f9.3,i3)
21    format(3f9.2,i3,f9.1)
22    format(i2,16f10.1)

      ystart(1)=rseg(1)         ! radius (km)
      ystart(3)=0.              ! epicentral distance (radians)
      ystart(4)=0.              ! time (sec)

      ! converge to each arrival in the list
      iconv=0                   ! counts actual arrivals
      do iar=1,narr
        iter=0
        nlegs=0                 ! some dummy values if no convergence
        trtime(iar)=0.
        tstar=0
        qray=0
        p=0
        if(jdebug.eq.1) then
          write(13,*) 'Iterating for iar=',iar
          write(13,*) 'iter        a1        a2        d1        d2'
        endif  
        if(abs(d2(iar)-d1(iar)).lt.0.01) then
          anew=0.5*(a1(iar)+a2(iar))
          goto 80
        endif  
40      anew=(del-d1(iar))*(a2(iar)-a1(iar))/(d2(iar)-d1(iar))+a1(iar)
        if(anew.ge.a2(iar).or.anew.le.a1(iar).or.mod(iter,2).eq.0)
     &       anew=0.5*(a1(iar)+a2(iar))         ! bisect
        ystart(2)=anew/r2d      ! angle i with vertical (radians,0=down)
        iter=iter+1
        call tracer(ystart,rseg,ktseg,kdwn,nseg,y,jsg,rmin,tstar,h,nray)

        ! KASRA
        r_pdiff_compare = 3482.01
        r_pdiff = 3482.01
        if (phase.eq.'Pdiff') then
          do pdiff_i=1,size(jsg)
            if (abs(jsg(pdiff_i)).gt.1.and.abs(jsg(pdiff_i)).lt.6) then
              y(1, pdiff_i) = r_pdiff
              y(2, pdiff_i) = pi/2.+0.1
            endif
            if (abs(jsg(pdiff_i)).eq.1) then
                if (y(1, pdiff_i).lt.r_pdiff_compare) then
                    y(1, pdiff_i) = r_pdiff
                    y(2, pdiff_i) = pi/2.+0.1
                endif
            endif
            if (abs(jsg(pdiff_i)).eq.6) then
                if (y(1, pdiff_i).lt.r_pdiff_compare) then
                    y(1, pdiff_i) = r_pdiff
                    y(2, pdiff_i) = pi/2.+0.1
                    ! y(3, pdiff_i) = (y(3, pdiff_i) + y(3, pdiff_i+1))/2.
                endif
            endif
          enddo
        endif
        ! END KASRA
        
        if(jdebug.gt.0.and.nray.le.0) 
     &    write(13,*) 'WARNING: tracer returned nray=',nray
        if(nray.le.0) goto 99
        qray=1.0e20
        if(tstar.gt.0.) qray=y(4,nray)/tstar
        aold=anew
        delnew=y(3,nray)*r2d
        if(jdebug.gt.0)
     &    write(13,43) iter,a1(iar),a2(iar),d1(iar),d2(iar),anew,delnew
43        format(i5,6f10.3)     
        if(abs(del-delnew).lt.0.01) goto 80
        if(iter.gt.itmax) then
          nray=0
          goto 99      ! lack of convergence
        endif  

        if((delnew-del)*(d1(iar)-del).gt.0.) then
          d1(iar)=delnew
          a1(iar)=anew
        else
          d2(iar)=delnew
          a2(iar)=anew
        endif  
        if((d1(iar)-del)*(d2(iar)-del).gt.0.) then
          nray=0
          goto 99   ! no convergence
        endif
        if(abs(d2(iar)-d1(iar)).lt.0.01) goto 80
        goto 40

        ! converged. Now do dynamic ray tracing
        ! anew is re-computed (but not used...)
80      if(abs(d1(iar)-d2(iar)).gt.0.)
     &   anew=(del-d1(iar))*(a2(iar)-a1(iar))/(d2(iar)-d1(iar))+a1(iar)

        call getelcor(scolatr,slonr,rcolatr,rlonr,y,nray,jsg,rseg,ktseg,
     &      kdwn,p,rayvel,rayq,ecorr)
c    &      kdwn,p,rayvel,rayq,ecorr2)
c       print *,'ecorr,ecorr2=',ecorr,ecorr2
c       write(13,*) 'ecorr,ecorr2=',ecorr,ecorr2
        call hessian(y,ypq,nray,jsg,rseg,ktseg,kdwn,hmf,hmb,p,rayvel,
     &      rayq)
        
        ! Kasra
        ! Comparison between two methods of ellipticity correction
        open(33, file='ellipticity_comparison.'//dataf)
        write(33, *) tcor, ',', ecorr, ',', abs(tcor-ecorr), 
     &       ',', abs(abs(tcor-ecorr)/tcor)*100.0

        ! Pdiff: BLNK method should be used 
        ! Otherwise, they will be compared and if the error is more than
        ! 0.1sec then it will stop the program....
        if(phase.eq.'Pdiff') then
          p = 4.439*180.0/pi
          ecorr = tcor
        else
          if (abs(tcor-ecorr).gt.0.1) then
            write(*,*) 'Different between two methods of ellipticity'
            write(*,*) 'correction is more than 0.1sec'
            stop
          endif
        endif
        ! END Kasra 
        
        call getlegs(y,rayvel,nray,legend,nlegs)
        ! correct travel time for dDelta
        trtime(iar)=y(4,nray)+p*(rdel-y(3,nray)) 
        if(jdebug.gt.0) write(13,*) 'Converged a,t=',anew,trtime(iar)

        iconv=iconv+1

        if(iconv.gt.1) goto 100

        ! Compute the crustal corrections at each reflection point
        ! telev is the elevation correction w.r.t. crust1.0
        ! The elevation correction of crust1.0 w.r.t. background 
        ! model is incorparated in tau
        icor=0
        iray=1
        ptlat=slat
        ptlon=slon          ! geographical location of the point
        rtarget=rsrc
        jup=0               ! sense switch for surface phases like pP
        if(kdwn(1).eq.0) jup=1
        
        ! XXX TEST XXX KASRA - testing crustal corrections 
        ! >>>> uncomment for testing
        !p = 0
        !do ptlat=-90,90
        !do ptlon=-180,180
        ! <<<<< 
        
        !  >>>> comment for testing
        do kseg=1,nseg
        ! <<<<<
        
        ! station elevation correction
        telev=0.
          
          ! >>>> uncomment for testing
          !kd=5
          !if(kseg.gt.2) jup=0
          !rtarget=6371.
          ! <<<<<

          !  >>>> comment for testing
          if(kseg.gt.2) jup=0
          kd=kdwn(kseg)
          rtarget=rseg(kseg)
          ! <<<<<
          
          if(kd.eq.5) iray=nray
          if(jdebug.eq.1) write(13,*) 'kseg,kd,jup,rtarget:',kseg,
     &               kd,jup,rtarget
          tau=0.
          if(kd.eq.2.or.kd.eq.4) goto 90         
          do while(abs(y(1,iray)-rtarget).gt.0.01.and.iray.lt.nray)
            iray=iray+1
          enddo  
          delx=y(3,iray)*r2d  ! distance to epicentre for segment end
          ! calculate azimuth of reflection or end point
          if((mod(del,360.)-180.)*(mod(delx,360.)-180.).ge.0.) then
            azx=az
          else
            azx=mod(az+180.,360.)
          endif
          
          !  >>>> comment for testing
          ! calculate coordinates of reflecion or endpoint:
          if(iray.eq.nray) then
            ptlat=rlat
            ptlon=rlon
          else
            call givloc(scolat,slon,delx,azx,ptlat,ptlon)
          endif
          ! <<<<<
          
          if(jdebug.eq.1) write(13,*) 'iray,delx,ptlat,ptlon=',
     &       iray,delx,ptlat,ptlon
          iray=min(nray,iray+2)  ! avoid problem in next section for PPP
          ! get crust1.0 at start of segment kseg
          ktype=ktseg(kseg)            ! is segment P or S?
          ktype1=ktype
          if(kd.eq.5) ktype=ktseg(kseg-1)
          if(kseg.gt.1) ktype1=ktseg(kseg-1)  ! and incoming?
          ! do departing ray at start of segment
          call getcr2(ptlat,ptlon,ktype,vcr2,rcr2,ncr2,ecr2,rmoho,
     &          vsed)
          
          !print *, '================================'
          !print *, '    vcr2 = ', vcr2(1)
          !print *, '    rcr2 = ', rcr2
          !print *, '    ncr2 = ', ncr2
          !print *, '    ecr2 = ', ecr2
          !print *, '    rmoho = ', rmoho
          !print *, '    vsed = ', vsed
          !print *, '================================'

          vsurf=vsed
          if(vsurf.eq.0)then
              write(*,*)'vsurf sucks: vsurf=', vsurf
          endif
          ! for islands take Vp=6.5, Vs=3.7 in elevation correction:
          if(rcr2(1).le.6370.) vsurf=9.3-2.8*ktype  ! basalt velocity
          ! if source is above crust1.0 surface
          ! replace source radius with crust1.0 surface radius
          ! this assumes r(n)=6371, otherwise the crustal correction
          ! makes no sense.
          if(kseg.eq.1.and.rsrc.gt.rcr2(1)) then
            vsurfref=vp(n)
            if(ktype.eq.2) vsurfref=vs(n)
            cosi=sqrt(1.0-(vsurfref*p/6371.0)**2)
            dt0=(rsrc-rcr2(1))*cosi/vsurfref ! origin time correction
            if(kd.eq.0) dt0=-dt0        ! ray is up from source
            write(7,50) ievt,slat,slon,rsrc,rcr2(1),dt0,kd
            rsrc=rcr2(1) 
            jshallow=jshallow+1
50          format(i8,2f9.2,f8.1,f10.1,f7.3,i3)
          endif
          call inttau(p,rcr2,rsrc,vcr2,ncr2,taucr2,tausrc2)
          tau=0                       ! intercept time
          ! Note: station elevation relev is wrt sealevel (r=6371), 
          ! not surface radius rcr2(1) of crust1.0

          ! find CR2 tau at segment start for kd=0,1 or 3 (i.e. if
          ! it starts from source or reflects at surface)
          if(kd.eq.0) then            ! if ray is up from source
            tau=tausrc2
            cosi=sqrt(1.0-(vsurf*p/6371.0)**2)
c commented out because telev already in crust1.0 tau
c           telev=cosi*(-rcr2(1)+6371.0)/vsurf  
            if(jdebug.gt.0) write(13,*) 'elev,rcr2,telev=',
     &         -rcr2(1)+6371.,rcr2(1),telev
          else if(kd.eq.1) then       ! down from source  
            tau=taucr2-tausrc2
          ! the following assumes surface reflection if r>6360
          else if(kd.eq.3.and.rtarget.gt.6360.) then 
            tau=taucr2
            cosi=sqrt(1.0-(vsurf*p/6371.0)**2)
c commented out because telev already in crust1.0 tau
c           telev=cosi*(-rcr2(1)+6371.0)/vsurf
            if(jdebug.gt.0) write(13,*) 'elev,rcr2,telev=',
     &         -rcr2(1)+6371.,rcr2(1),telev
          endif  
          if(jdebug.eq.1.and.kd.lt.5) 
     &          write(13,*) 'Outgoing tau,tel=',tau,telev

        !  >>>> uncomment for testing
        !  open(55, file='check_crust')
        !  write(55, *) ptlat, ',', ptlon, ',', vcr2(1)
        ! <<<<<

          ! at this point tau is the tau travelled in CR2 from start
          ! of the segment
          ! now add contribution from ray *to* start of segment in
          ! case kd=3 or 5 (i.e. surface reflection or end point)
          if(ktype.ne.ktype1) then      ! if wave type differs
            call getcr2(ptlat,ptlon,ktype1,vcr2,rcr2,ncr2,ecr2,rmoho,
     &          vsed)
            vsurf=vsed
            if(rcr2(1).le.6370.) vsurf=9.3-2.8*ktype1  ! basalt velocity 
            call inttau(p,rcr2,rsrc,vcr2,ncr2,taucr2,tausrc2)
          endif
          ! avoid double counting of pP if jup=1
          if(kd.eq.3.and.rtarget.gt.6360.and.jup.eq.0) then
            tau=tau+taucr2
            cosi=sqrt(1.0-(vsurf*p/6371.0)**2)
c           telev=telev+cosi*(-rcr2(1)+6371.0)/vsurf
            if(jdebug.gt.0) write(13,*) 'elev,rcr2,telev=',
     &         -rcr2(1)+6371.,rcr2(1),telev
          else if(kd.eq.5.and.jup.eq.0) then
            tau=taucr2
            cosi=sqrt(1.0-(vsurf*p/6371.0)**2)
            telev=cosi*(relev-rcr2(1)+6371.0)/vsurf
c            write(*,*)'telev=', telev
c            if(ABS(telev).gt.10000) then
c                write(*,*)'+++++++++++++++++++++++++++++++++++'
c                write(*,*)'relev=', relev
c                write(*,*)'rcr2(1)=', rcr2(1)
c                write(*,*)'relev-rcr2(1)+6371.0=', relev-rcr2(1)+6371.0
c                write(*,*)'vsurf=', vsurf
c                write(*,*)'telev=', telev
c            endif
            if(jdebug.gt.0) write(13,*) 'elev,rcr2,telev=',
     &         relev-rcr2(1)+6371.,rcr2(1),telev
          endif  
          if(jdebug.eq.1.and.kseg.gt.1) 
     &          write(13,*) 'Incoming tau,telev=',taucr2,telev

          ! subtract the reference model tau
          ktype=ktseg(kseg)          ! is departing ray P or S?
          if(kd.eq.5) ktype=ktseg(kseg-1)
          if(kseg.gt.1) ktype1=ktseg(kseg-1)  ! and incoming?
          ! do departing first
          call getref(ktype,vref,rref,nref,rmoho)
          call inttau(p,rref,rsrc,vref,nref,tauref,tausrcr)
          if(kd.eq.0) then            ! if ray is up from source
            tau=tau-tausrcr
          else if(kd.eq.1) then       ! down from source  
            tau=tau-(tauref-tausrcr)
          ! the following assumes surface reflection if r>6360
          else if(kd.eq.3.and.rtarget.gt.6360.) then 
            tau=tau-tauref
          endif  
          if(jdebug.eq.1.and.kd.lt.5) 
     &            write(13,*) 'Outgoing ref tau=',tauref,tausrcr
          ! now add incoming ray
          if(ktype.ne.ktype1) then    ! if wavetype differs
            call getref(ktype1,vref,rref,nref,rmoho)
            call inttau(p,rref,rsrc,vref,nref,tauref,tausrcr)
          endif  
          if(kd.eq.3.and.rtarget.gt.6360.and.jup.eq.0) then
            tau=tau-tauref
          else if(kd.eq.5.and.jup.eq.0) then
            tau=tau-tauref
          endif  
          if(jdebug.eq.1.and.kseg.gt.1) 
     &          write(13,*) 'Incoming ref tau=',tauref

90        if(kd.ne.2.and.kd.ne.4) then
            write(2,95) kseg,kd,ptlat,ptlon,rtarget,tau,telev
            ! Kasra
            open(44, file='ell_ccor.'//dataf)
            write(44, *)  ptlat, ',', ptlon, ',', rtarget, ',', 
     &         ecorr, ',', 
     &         tau, ',', telev 
          else
            write(2,96) kseg,kd,0,0,rtarget,0,0
          endif
95        format(2i3,2f8.2,f8.1,2f8.3)
96        format(2i3,2i8,f8.1,2i8)
          if(jdebug.eq.1)
     &      write(13,*) 'Final crustal correction:',tau,' Elev:',telev
        
        print *, tau
        if (isnan(tau)) print *, tau, ptlat, ptlon
        
        if(kd.ne.2.and.kd.ne.4) then
          if (phase.eq.'Pdiff') then
              trtime_final = tb_mfi(1)
          else
              trtime_final = trtime(iar)
          endif
          open(222, file='info_collector.'//dataf)
          write(222, *) kd, ',', stationcode, ',', netw, ',',
     &          ptlat, ',', ptlon, ',', rtarget, ',', 
     &          slat, ',', slon, ',', ecorr, ',',
     &          tau, ',', telev, ',', tstar, ',',
     &          del, ',', trtime_final, ',', tobs(1), ',', sdep
          endif

        !  >>>> uncomment for testing
        !  enddo
        ! <<<<<

        enddo         ! end of loop over segments for first arriving ray

        goto 100

        ! write one segment line with kseg=0 so raymatrix knows rest is skipped
99      write(2,96) 0,0,0,0,0.,0,0
        nray=0       ! just to be sure


100     if (phase.eq.'Pdiff') then
            trtime_final = tb_mfi(1)
        else
            trtime_final = trtime(iar)
        endif
        
        if(kd.eq.5) then
            open(111, file='synthetic_calculated.'//dataf)
            write(111, *) kd, ',', ievt, ',', idate, ',', 
     &             stationcode, ',', 
     &             netw, ',', ptlat, ',', ptlon, ',', rtarget, ',',
     &             slat, ',', slon, ',', ecorr, ',',
     &             abs(6371-rcr2(1)), ',',
     &             tau, ',', 
     &             ecr2(1)*1000,  ',',
     &             telev, ',', tstar, ',',
     &             del, ',', tb_mfi(1), ',', trtime(iar), ',',
     &             abs(trtime(iar)-tb_mfi(1))/tb_mfi(1)*100, ',',
     &             tobs(1), ',', sdep
        endif

        ! We have already considered the effects of attenuation in
        ! our forward modeling (YSPEC or AXISEM) ---> t*=0.
        ! TODO: it should be more generic! means that it can detect
        ! whether it is yspec or any other code!
!100     write(2,105) iar,nray,nlegs,trtime(iar),ecorr,tstar,qray,p
        write(2,105) iar,nray,nlegs,trtime_final,ecorr,0.,qray,p
105     format(i2,2i5,3f10.2,2e12.3)
        ! END Kasra

        write(2,fmt='(20i5)') (legend(i),i=1,nlegs)
        iq=n
        do i=1,nray
          c=rayvel(i)
          q0=rayq(i)
          h11=hmf(1,i)+hmb(1,i)
          h22=hmf(2,i)+hmb(2,i)
          write(2,120) (y(j,i),j=1,3),c,q0,h11,h22
120       format(f10.3,3f10.5,f10.6,4e14.6)
        enddo

      enddo     ! end of do loop over iar

      goto 10

900   print *, 'End of data file reached'
      print *, 'Total number of raypaths:',ndata
      write(4,*) 'Total number of raypaths (data):',ndata
      close(7)
      open(8,file='xclm.'//dataf) ! maximum xcorr length for each band
      write(8,*) 'band    xclm'
      do ib=1,nband
         write(8,fmt='(i4,f10.3)') ib,xclm(ib)
      enddo
      close(8)
      if(jshallow.gt.0) then
        print *, 'WARNING: The source depths in ',jshallow,' paths with'
        print *, 'sources shallower than crust1.0 surface were replaced'
        print *, 'with source depth at crust1.0 surface.'
      endif

      end
      
      subroutine findarr(del,tdifmax,ytbl,ntbl,a1,a2,d1,d2,trtime,narr)

      ! find all arrivals at distance del and sort them, first
      ! arrival first

      ! input:
      ! del=epicentral distance in degrees
      ! tdifmax=max time differential for later arrivals
      ! ytbl(1,i)=ray angles, (2,i)=distances (deg), (3,i)=times (s)
      ! ntbl=number of rays in table ytbl
      ! ray angles can be ascending or descending but must be monotonic

      ! output
      ! narr=number of arrivals, then for each one:
      ! a1,a2=bracket for ray angle (a1<a2) in degrees
      ! d1,d2=delta's belong to a1,a2 (degr)
      
      dimension a1(10),a2(10),d1(10),d2(10),trtime(10)
      dimension ytbl(4,4000),indx(10)

      jdebug=0
      if(jdebug.gt.0) write(13,*) 'subroutine findarr,del=',del
      if(jdebug.gt.0) write(13,*) '#      i        Delta(deg)   T(s)'

      narr=0
      if(jdebug.gt.0)
     &  write(13,*) 1,ytbl(1,1),ytbl(2,1),ytbl(3,1)
      do i=2,ntbl
        if(jdebug.gt.0)
     &    write(13,*) i,ytbl(1,i),ytbl(2,i),ytbl(3,i)
        if((del-ytbl(2,i-1))*(del-ytbl(2,i)).lt.0.0
     &        .or.del.eq.ytbl(2,i-1)) then
          narr=min(10,narr+1)
          a1(narr)=ytbl(1,i-1)
          a2(narr)=ytbl(1,i)
          d1(narr)=ytbl(2,i-1)
          d2(narr)=ytbl(2,i)
          dd=d2(narr)-d1(narr)
          trtime(narr)=ytbl(3,i-1)
          if(dd.ne.0.) trtime(narr)=ytbl(3,i-1)+(del-d1(narr))*
     &          (ytbl(3,i)-ytbl(3,i-1))/dd
          if(jdebug.gt.0)
     &      write(13,*) 'bracket:',narr,dd,trtime(narr)
          if(narr.ge.10) goto 20
        endif
      enddo

      ! swap if a1>a2
      do i=1,narr
        if(a1(i).gt.a2(i)) then
          aa=a2(i)
          a2(i)=a1(i)
          a1(i)=aa
          aa=d2(i)
          d2(i)=d1(i)
          d1(i)=aa
        endif
      enddo  

20    if(narr.le.1) return

c     write(13,*) 'sort:'
c     do i=1,narr
c       write(13,9) i,trtime(i),a1(i),a2(i),d1(i),d2(i)
c     enddo
      call indexx(narr,trtime,indx)
      if(jdebug.eq.1) write(13,*) 'indx=',(indx(i),i=1,narr)
      call sort(narr,a1,indx)
      call sort(narr,a2,indx)
      call sort(narr,d1,indx)
      call sort(narr,d2,indx)
      call sort(narr,trtime,indx)

c     write(13,*) 'gives:'
c     do i=1,narr
c       write(13,9) i,trtime(i),a1(i),a2(i),d1(i),d2(i)
c     enddo

      nn=narr
      do i=1,nn
        if(trtime(i)-trtime(1).le.tdifmax) narr=i
      enddo  

c     write(13,*) 'tdifmax',tdifmax,' gives narr=',narr
9     format(i3,5f10.3)

      return
      end

      SUBROUTINE INDEXX(N,ARRIN,INDX)
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END

      SUBROUTINE SORT(N,RA,INDX)
      DIMENSION RA(10),WKSP(10),INDX(10)
      DO 11 J=1,N
        WKSP(J)=RA(J)
11    CONTINUE
      DO 12 J=1,N
        RA(J)=WKSP(INDX(J))
12    CONTINUE
      RETURN
      END