      program selectISC

      ! this program selects ISC data from Bob Engdahls
      ! files (eg 98-00RES) and writes them in a form
      ! digestible by raydata.f (specifying 0 bands for
      ! ray theory data), or in the original format.

      parameter(nfmx=50)
      character*80 infile1(nfmx),txtline,outfile,twinfile,gmtfile,
     &      evtfile

      character*8 phasej,selphase
      character*8 kpha
      character*6 sta
      character*1 onset, ahyp, bpsel
      character*3 isol
      character*2 iseq
      logical accept,ppss
      integer eventdate(100000),eventime(100000),kevt(100000)
      integer yr,day,hr,phj,phi,pho,prec,nev,nreject(20)
      real openaz, openaz2
      data nreject/20*0/,ievt0/0/

      elatlo=1.0e10
      elathi=-1.0e10
      slatlo=1.0e10
      slathi=-1.0e10
      elonlo=1.0e10
      elonhi=-1.0e10
      slonlo=1.0e10
      slonhi=-1.0e10
      ipreclo=1000
      iprechi=-1000
      fmblo=1.0e10
      fmbhi=-1.0e10
      deltalo=1.0e10
      deltahi=-1.0e10
      azlo=1.0e10
      azhi=-1.0e10
      ntotlo=1000
      ntothi=-1000
      nlocallo=1000
      nlocalhi=-1000



      i=1
10    print *,'Give input file or <stop>:'
      read(*,'(a)',end=15) infile1(i)
      if(infile1(i).eq.'stop') goto 15
      i=i+1
      goto 10
      

15    kfile=i-1
      klons=0
      klone=0

      print *,'Give output file name (e.g. iscPn1):'
      read (*,'(a)') outfile
      print *,'Give GMT file or <None>:'
      read (*,'(a)') gmtfile
      if(gmtfile.ne.'None') open(3,file=gmtfile)
      print *,'Choose output format:'
      print *,'   [1] Same as Engdahl files'
      print *,'   [2] Suitable for raydata.f:'
      read *,kfm
      twinfile='None'
      if(kfm.eq.2) then
        print *,'Give twinfile for PP, SS etc or <None>:'
        read (*,'(a)') twinfile
        print *,'Do you start a new event file [1] or not [0]?'
        read *,neweventfile
        if(neweventfile.eq.1) then
          print *,'New event file is written'
        else
          print *,'We use existing event file'
        endif  
        print *,'Give event file name:'
        read (*,'(a)') evtfile
        open(4,file=evtfile)
        nnewevt=0       ! running count of events
        neventf=0       ! count of events already known at start
        if(neweventfile.eq.0) then
          print *,'Events not on file will be numbered and added'
          print *,'to ',evtfile
          print *,'Give number to start added event numbering:'
          print *,'(if -1 numbering continuous from existing file)'
          read *,inewevt
          ios=0
          i=0
          do while(ios.eq.0)
            i=i+1
            if(i.gt.100000) stop 'i>100000 - Increase dimensions'
            read(4,*,iostat=ios) kevt(i),eventdate(i),eventime(i)
          enddo
          neventf=i-1
          nnewevt=i-1
          if(inewevt.lt.0) inewevt=i
          print *,i,' events read from file'
        else
          print *,'Give starting number for event numbering:'
          read *,inewevt
        endif  
      endif  
      
      print *,'Give phase to select (e.g. Pn, or ALL):'
      read(*,'(a)') selphase
      if(selphase.eq.'ALL'.and.kfm.eq.2) stop 'ALL not allowed if kfm=2'
      print *,'Give South,North latitude box for event selection:'
      read *,elat1,elat2
      print *,'Normally, we work with 0<longitude<360, but if you'
      print *,'specify  W lon <0 I assume -180<lon<180 works better.'
      print *,'Give West,East box longitude for event selection:'
      read *,elon1,elon2
      if(elon1.lt.0.) klone=1
      print *,'Give South,North latitude for station:'
      read *,slat1,slat2
      print *,'Give West,East longitude for station:'
      read *,slon1,slon2
      if(slon1.lt.0.) klons=1
      print *,'Give min, max magnitude:'
      read *,xm1,xm2
      xm1=xm1-0.001
      xm2=xm2+0.001     ! guard aganis roundoff
      print *,'Give min,max epicentral distance (deg):'
      read *,del1,del2
      print *,'Give max log decimals (e.g. -1):'
      read *,iprecmax
      print *,'Give max open azimuth (e.g. 90):'
      read *,openazmax
      print *,'Give min for local and worldwide stations (eg 3,20):'
      read *,nrloc,nrall
      print *,'Select on bounce point (y/n - y for PP or SS only)?'
      read *,bpsel
      ppss=.false.
      if(bpsel.eq.'y') ppss=.true.
      if(ppss) then
        print *,'Give min, max latitude of bounce point (PP,SS):'
        read *,bplat1,bplat2
        print *,'Give min, max longitude of bounce point (PP,SS):'
        read *,bplon1,bplon2
      endif

      koball=0
      kobaccept=0
      kobphase=0
      ievt=ievt0-1

      open(11,file=outfile)
      if(kfm.eq.2) then
        write(11,'(a)') outfile
        write(11,'(a)') twinfile
        write(11,'(a)') selphase
        write(11,'(a)') selphase

        print *,'Give estimate of standard deviation (s):'
        read *,sigma
      
        print *,'Input of ray segments (one per line):'
        print *,'Give radius (km), kdwn, wavetype (P=1,S=2),'
        print *,'where kdwn=1 for downward start, 0 for upward start,'
        print *,'2 for max turning depth, 3 for reflection, 4 for'
        print *,'transmission, 5 for end of ray'
  
        read(*,*) rseg,kdwn,ktseg
        write(11,40) rseg,kdwn,ktseg
40      format(f8.1,2i5,2x,a4,2x,a1)
        if(kdwn.gt.1) then
          print *,'First ray point should have kdwn 1 or 0'
          stop
        endif  
        do while(kdwn.ne.5)
          read(*,*) rseg,kdwn,ktseg
          write(11,40) rseg,kdwn,ktseg
        end do
        write(11,*) 0     ! zero bands
      endif

      do i=1,kfile

        write(txtline,'(a)') infile1(i)
        write(*,'(a)') infile1(i)
        open(10,file=txtline)

c Example data line:
c 256423  LEQ    82.0  89.5 1998  1  1 1  0 12 38.57  28.615 218.963   8.4 0.0 0.0   35    0      BALM    29.128 217.656  1.300   0.814 231.438        Pg         0   0   0      19.0612   0.814 231.438    0.0        0.000   0.000   0.000  0.000   0.00   0.00          13.52 -2     15.65  -2.13        0.00   0.00   0.00  -2.13 0 0.39        0.000      0.00        0.840  -2.10

20      read(10,30,end=998,iostat=ios)
     &    nev,ahyp,isol,iseq,openaz,openaz2,
     &    iyr,imon,iday,ihold,ihr,imin,sec,
     &    elat,elon,depth,fmb,fms,ntot,ntel,
     &    sta,slat,slon,elev,delta,azim,
     &    onset,phasej,iphj,iphi,ipho,
     &    rdtdd,rdelta,razim,dbot,
     &    gblat,gblon,stadel,bdep,tbath,twater,
     &    obstt,iprec,prett,rawres,
     &    ecor,scor,elcor,resid,iflg,wgt,
     &    tdelta,ttime,
     &    delisc,resisc

30      format(i7,1x,a1,a3,a2,2f6.1,i5,2i3,i2,2i3,
     &    f6.2,2f8.3,f6.1,2f4.1,2i5,5x,
     &    1x,a6,2f8.3,f7.3,2f8.3,5x,
     &    1x,a1,1x,a8,3i4,5x,
     &    f8.4,2f8.3,f7.1,5x,
     &    3f8.3,f7.3,2f7.2,5x,
     &    f10.2,i3,f10.2,f7.2,5x,
     &    4f7.2,i2,f5.2,5x,
     &    f8.3,f10.2,5x,f8.3,f7.2)

       if(ios.ne.0) then
         print *,'Irregular ending of file?'
         goto 998
       endif

       koball=koball+1
       if(selphase.ne.'ALL'.and.phasej.ne.selphase) then
         nreject(1)=nreject(1)+1
         goto 20
       endif  
       if(isol.eq.'XEQ') then
         nreject(2)=nreject(2)+1
         goto 20
       endif  
       kobphase=kobphase+1
 
       ! go from colat to latitude
       slatm=90.0-slat
       elatm=90.0-elat
       bplat=90.0-bplat
 
       ! negative lon, if needed
       slonm=slon
       elonm=elon
       if(klone.eq.1.and.elon.gt.180.) elonm=elon-360.
       if(klons.eq.1.and.slon.gt.180.) slonm=slon-360.
 
 
       ! some input file statistics
       elatlo=min(elatm,elatlo)
       elathi=max(elatm,elathi)
       slatlo=min(slatm,slatlo)
       slathi=max(slatm,slathi)
       elonlo=min(elonm,elonlo)
       elonhi=max(elonm,elonhi)
       slonlo=min(slonm,slonlo)
       slonhi=max(slonm,slonhi)
       ipreclo=min(iprec,ipreclo)
       iprechi=max(iprec,iprechi)
       fmblo=min(fmblo,fmb)
       fmbhi=max(fmbhi,fmb)
       deltalo=min(rdelta,deltalo)
       deltahi=max(rdelta,deltahi)
       azlo=min(azlo,openaz)
       azhi=max(azhi,openaz)
       nlocal=ntot-ntel
       nlocallo=min(nlocallo,nlocal)
       nlocalhi=max(nlocalhi,nlocal)
       ntotlo=min(ntotlo,ntot)
       ntothi=max(ntothi,ntot)
 
       accept=nlocal.ge.nrloc.and.ntot.ge.nrall
       if(.not.accept) then
         nreject(3)=nreject(3)+1
         goto 20
       endif  
       accept=accept.and.elatm.gt.elat1.and.elatm.lt.elat2
       if(.not.accept) then
         nreject(4)=nreject(4)+1
         goto 20
       endif
       accept=accept.and.elonm.gt.elon1.and.elonm.lt.elon2
       if(.not.accept) then
         nreject(5)=nreject(5)+1
         goto 20
       endif
       accept=accept.and.slatm.gt.slat1.and.slatm.lt.slat2
       if(.not.accept) then
         nreject(6)=nreject(6)+1
         goto 20
       endif
       accept=accept.and.slonm.gt.slon1.and.slonm.lt.slon2
       if(.not.accept) then
         nreject(7)=nreject(7)+1
         goto 20
       endif
       accept=accept.and.fmb.ge.xm1.and.fmb.le.xm2
       if(.not.accept) then
         nreject(8)=nreject(8)+1
         goto 20
       endif
       accept=accept.and.rdelta.gt.del1.and.rdelta.lt.del2
       if(.not.accept) then
         nreject(9)=nreject(9)+1
         goto 20
       endif
       accept=accept.and.iprec.le.iprecmax
       if(.not.accept) then
         nreject(10)=nreject(10)+1
         goto 20
       endif
       accept=accept.and.openaz.lt.openazmax
       if(.not.accept) then
         nreject(11)=nreject(11)+1
         goto 20
       endif
       if(ppss) then
         accept=accept.and.gblat.gt.bplat1.and.gblat.lt.bplat2
         if(.not.accept) then
           nreject(12)=nreject(12)+1
           goto 20
         endif
         accept=accept.and.gblon.gt.bplon1.and.gblon.lt.bplon2
         if(.not.accept) then
           nreject(13)=nreject(13)+1
           goto 20
         endif
       endif  
 
       if(.not.accept) goto 20
 
       kobaccept=kobaccept+1    ! accept the datum from this line
 
       if(kfm.eq.2) then  
         idate=julian(iyr,imon,iday)
         iotime=ihr*10000+imin*100+sec
         call findievt(idate,iotime,kevt,eventdate,eventime,
     &           ievt,inewevt,nnewevt)
         write(11,50) idate,iotime,ievt,0,
     &     sta,'ISC','UNK',slatm,slonm,depth,elatm,elonm,elev,1,0,0
50       format(3i8,i4,1x,a8,1x,a3,1x,a3,2f9.3,f7.1,2f9.3,f7.3,3i3)
         write(11,'(3i2)') 0,0,0         ! noise line 
         write(11,'(2f9.2,4i2)',iostat=ios) obstt,sigma,0,0,0,0   ! data line
       else if(kfm.eq.1) then
        write(11,30)
     &    nev,ahyp,isol,iseq,openaz,openaz2,
     &    iyr,imon,iday,ihold,ihr,imin,sec,
     &    elat,elon,depth,fmb,fms,ntot,ntel,
     &    sta,slat,slon,elev,delta,azim,
     &    onset,phasej,iphj,iphi,ipho,
     &    rdtdd,rdelta,razim,dbot,
     &    gblat,gblon,stadel,bdep,tbath,twater,
     &    obstt,iprec,prett,rawres,
     &    ecor,scor,elcor,resid,iflg,wgt,
     &    tdelta,ttime,
     &    delisc,resisc
       endif

       if(gmtfile.ne.'None') then
         write(3,*) slonm,slatm
         write(3,*) elonm,elatm
         write(3,'(a)') '>'
       endif  
    
      goto 20
  998 continue

      enddo

      print *,'elatlo=', elatlo
      print *,'elathi=', elathi
      print *,'slatlo=', slatlo
      print *,'slathi=', slathi
      print *,'elonlo=', elonlo
      print *,'elonhi=', elonhi
      print *,'slonlo=', slonlo
      print *,'slonhi=', slonhi
      print *,'ipreclo=', ipreclo
      print *,'iprechi=', iprechi
      print *,'fmblo=', fmb
      print *,'fmbhi=', fmb
      print *,'deltalo=', deltalo
      print *,'deltahi=', deltahi
      print *,'azlo=', openaz
      print *,'azhi=', openaz

      print *,'Total nr of input lines=',koball
      print *,selphase,' :',kobphase
      print *,'Accepted: ',kobaccept


      print *,' '
      print *,'Rejection statistics:'
      print *,'Rejected because of phase mismatch:',nreject(1)
      print *,'Rejected because of XEQ (bad quality):',nreject(2)
      print *,'Rejected because of few stations:',nreject(3)
      print *,'Rejected because of event latitude:',nreject(4)
      print *,'Rejected because of event longitude:',nreject(5)
      print *,'Rejected because of station latitude:',nreject(6)
      print *,'Rejected because of station longitude:',nreject(7)
      print *,'Rejected because of magnitude:',nreject(8)
      print *,'Rejected because of distance:',nreject(9)
      print *,'Rejected because of low precision:',nreject(10)
      print *,'Rejected because of large open azimuth:',nreject(11)
      print *,'Rejected because of bounce latitude:',nreject(12)
      print *,'Rejected because of bounce longitude:',nreject(13)

      if(nnewevt.ne.neventf) then
        print *,'Event count at start:',neventf
        print *,'Event count at end:',nnewevt
        if(neweventfile.eq.0) then
          print *,'New events were added to the list'
          print *,'Give file name to write extended list of events:'
          read (*,'(a)') evtfile
          close(4)
          open(4,file=evtfile)
        endif
        do i=1,nnewevt
          write(4,*) kevt(i),eventdate(i),eventime(i)
        enddo
      else
        print *,'No new events detected - old event file is still OK'
      endif  
 
      stop
      end

      function julian(iyr,imon,iday)

      ! returns julian data in format YYYYDDD
      ! input year,month,day. 
      ! See code for how I deal with truncated iyr. This won't work
      ! after 2050....

      dimension mday(12)
      data mday/0,31,59,90,120,151,181,212,243,273,304,334/
      
      jj=iyr
      ! is iyr truncated to two digits? Assume 1950<jj<2050
      if(iyr.lt.50) jj=jj+2000
      if(iyr.lt.100) jj=jj+1900

      ! leap year?
      leap=0
      if(mod(jj,4).eq.0) leap=1
      if(mod(jj,100).eq.0) leap=0
      if(mod(jj,400).eq.0) leap=1

      julian=jj*1000+mday(imon)+iday
      if(imon.gt.2) julian=julian+leap
      
      return
      end
         
      subroutine findievt(idate,iotime,kevt,eventdate,eventime,
     &    ievt,inewevt,nnewevt)

      ! compares date and time to event list and returns ievt

      integer kevt(100000),eventdate(100000),eventime(100000)

      ! see if there is a match with the current list
      do i=1,nnewevt
        if(idate.eq.eventdate(i).and.iotime.eq.eventime(i)) then
          ievt=kevt(i)
          return
        endif  
      enddo

      ! in case there is no fit, add event to list:
      nnewevt=nnewevt+1
      eventdate(nnewevt)=idate
      eventime(nnewevt)=iotime
      ievt=inewevt
      inewevt=inewevt+1
      kevt(nnewevt)=ievt
      return
      end
