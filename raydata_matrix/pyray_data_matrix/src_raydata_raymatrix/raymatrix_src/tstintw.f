      program tstintw

      ! used to test quadrature precision of dintw.f
      ! uses Hung et al GJI 146:289, 2001 eq (13):
      ! |mdot(w)| = (w*tau)/sqrt(2*pi) * exp[-w^2*tau^2/(8*pi^2)]

      parameter (NTM=5000)

      dimension tknl(NTM,2,25),aknl(NTM,2,25),nfreq(25)
      dimension bndomega(NTM,25),dotm(NTM,25),dtbl(25),ntbl(25)

      nband=10
      nw=40  
      tau=10.
      t2=tau**2
      t4=t2*t2
      t6=t2*t4
      pi=3.14159265
      p2=pi*pi
      p4=p2*p2
      wmax=30./tau
      wmin=0
      dw=(wmax-wmin)/(nw-1)
c     nw=27
      open(1,file='mtst.xy')
      open(2,file='Gauss1')
      do j=1,9 !8 ! skip header
        read(2,fmt='(a72)') tmp
      enddo
      do j=1,nband
        read(2,*) tmp
        write(1,*) 'band ',j
        nfreq(j)=nw
c       w=wmin
        do i=1,nw
          read(2,*) bndomega(i,j),dotm(i,j)
          write(1,*) bndomega(i,j),dotm(i,j)
c         bndomega(i,j)=w
c         dotm(i,j)=0.3989*w*tau*exp(-0.01266*(w*tau)**2)
c         write(1,*) w,dotm(i,j)
c         w=w+dw
        enddo
      enddo
      close(1)
      close(2)

      open(1,file='Tth.xy')
      do i=1,250
        a=(i-1)*0.2
        a2=a*a
        f=2563.7*a*exp(-a2*p2/t2)*(a2*a2-5*t2*a2/p2 +3.75*t4/p4)/t6
        write(1,*) a,f
      enddo
      close(1)

      call dintw(tknl,aknl,dtbl,ntbl,dotm,bndomega,nfreq,nband)

      end  
!--------------------------------------------------

        subroutine dintw(tknl,aknl,dtbl,ntbl,dotm,bndomega,nfreq,nband)

        ! computes integrals with sin(w*dt) and cos(w*dt) terms 
        ! for fast interpolation - see Dahlen, eq (78), eg:
        ! tknl = -(1/2*pi) int [w^3 dotm^2 sin Phi]/ int [w^2 dotm^2]
        ! and
        ! aknl = -(1/2*pi) int [w^2 dotm^2 cos Phi]/ int [dotm^2]

        ! input: dotm,bndomega,nfreq and nband
        !        dotm(i,j) is amplitude spectrum; i ranges from 1 to
        !        nfreq(j), j from 1 to nband. bndomega(i,j) is the
        !        circle frequency (rad/s) for dotm(i,j).
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

        parameter(NTM=5000, NFR=25)      ! keep same as in main

        dimension tknl(NTM,2,NFR),aknl(NTM,2,NFR),nfreq(NFR)
c       dimension bndomega(200,NFR),dotm(200,NFR),dtbl(NFR),ntbl(NFR)
        dimension bndomega(NTM,25),dotm(NTM,25),dtbl(25),ntbl(25)
        dimension gt1(NTM),gt2(NTM),ga1(NTM),ga2(NTM)
        dimension gt1old(NTM),gt2old(NTM),ga1old(NTM),ga2old(NTM)

        data twopi/6.2831853/
        data dtmax/50./         ! max detour time 50 seconds

        jdebug=1


        if(nband.gt.NFR) stop 'dintw.f: too many frequency bands'

            open(31,file='TGauss1.xy')
            open(32,file='AGauss1.xy')
        do kb=1,nband           ! we have data for nband freq bands

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
          dwmin=0.01*twopi/dtmax  ! min step size in omega*delay is pi/50
          print *,'dwmin=',dwmin,' for wmax=',wmax

          dtbl(kb)=0.05*twopi/wmax      ! step for dT interpolation
          ntbl(kb)=dtmax/dtbl(kb)+1
          if(ntbl(kb).gt.NTM) then
            print *,'kb,dtmax,wmax=',kb,dtmax,wmax
            print *,'dtbl,ntbl=',dtbl(kb),ntbl(kb)
            stop 'dintw: dtmax too large (nt>NTM)'
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

          if(jdebug.eq.1) then
            write(31,*) ntbl(kb)
            write(32,*) ntbl(kb)
            do i=1,ntbl(kb)
              write(31,*) (i-1)*dtbl(kb),-twopi*tknl(i,1,kb)
              write(32,*) (i-1)*dtbl(kb),-twopi*aknl(i,2,kb)
            enddo
          endif  
  
        end do          ! end loop for frequency bands
            close(31)
            close(32)

        return
        end


