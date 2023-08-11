        program mat2asc

        ! converts binary file matrixT.xxxx or matrixA.xxxx
        ! into ascii version
        ! see also mat2bin

        character*30 fname
        character stationcode*16,netw*8,comp*3
        dimension jssw(3),kountpar(3)
        dimension asparse(300000),ja(300000)
        dimension period(20)

        print *,'Give matrix file name (eg matrixT.P):'
        read(*,fmt='(a)') fname
        open(11,file=fname,form='unformatted')
        open(12,file='ascii.'//fname)

        read(11) linked,fname
        write(12,*) linked,' ',fname
        read(11) np,jssw
        write(12,fmt='(i8,3i5)') np,jssw
        read(11) nband,(period(i),i=1,nband)
        write(12,fmt='(i5,(20f7.1))') nband,(period(i),i=1,nband)
        print *,'np=',np
        print *,'jssw=',jssw
        print *,'nband=',nband
        if(nband.gt.20) stop 'Increase dimension of period()'
        print *

100     read(11,end=200) ndata,kount,iband,dobs,dtheor,s,ecorr,corcrust,
     &    corelev,tstar,qray,ievt,kluster,stationcode,netw,comp,
     &    slowness,raz,rdel,angle0,vsrc,krtyp
        if(kount.gt.300000) stop 'Increase dimension of asparse,ja'
        read(11) kountpar
        read(11) (ja(i),i=1,kount)
        read(11) (asparse(i),i=1,kount)
        write(12,110) ndata,kount,iband,dobs,dtheor,s,ecorr,corcrust,
     &    corelev,tstar,qray,ievt,kluster,stationcode,netw,comp,
     &    slowness,raz,rdel,angle0,vsrc,krtyp
110     format(2i8,i4,2f10.2,5f7.2,f10.1,i8,i3,1x,a16,1x,a8,1x,a3,
     &         f10.1,3f10.5,f8.3,i2,2f10.2,i8)
        write(12,*) kountpar
        write(12,*) (ja(i),i=1,kount)
        write(12,*) (asparse(i),i=1,kount)
        if(mod(ndata,100).eq.0) write(6,'(a1,$)') 'x'
        goto 100

200     print *,'ndata=',ndata
        stop 'end of file reached'
        end

