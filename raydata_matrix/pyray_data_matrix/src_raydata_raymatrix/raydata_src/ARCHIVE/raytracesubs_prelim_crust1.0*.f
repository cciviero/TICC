      ! ALPHA version
      ! testing done on a limited set of phases

      ! List of changes:
      ! June 06 - warnings for ray length changed to 8998
      ! July 06 - changed mkraytable to avoid rays trapped in LVZ
      ! Aug 06 - fixed bug in trapezoidal rule in inttau for sources
      !          not on a model interface
      !          fixed nleg bug in getlegs for reflected phases
      !          fixed y() index bug in getlegs for reflected phases
      !          adjusted criterion to recognize reflections in getlegs
      !          fixed bug in continuity of reflected phases at CMB
      ! Sep 06 - in hessian:
      !          Adjusted criterion to avoid taking a turning point as reflection point
      !          Treat the region near source/receiver as a homogeneous
      !                sphere to avoid oscillation of hessian elements.
      ! Nov 06 - Output warning message of imprecise hessian against analytical expression
      !          Add precision test of A53 on backward H11
      ! 2010/08/12 K.S. change name of common block from /mod/ to /modl/
      !  since old name was a fortran standard violation and ifort would
      ! not compile. Also changed in raytracesubs.f
      ! raydata.f(237): error #8038: A referenced intrinsic procedure 
      ! can not have the same name as a common block.   [MOD]


      subroutine hessian(y,ypq,nray,jsg,rseg,ktseg,kdwn,hmf,hmb,p,
     &      rayvel,rayq)

      ! calculate forward and backward hessian matrix (hmf,hmb)
      ! equation numbers refer to Dahlen et al, GJI 141:157-174,2000

      ! input:
      ! y(i,j) - ray geometry from subroutine tracer, j=1,nray
      ! nray - number of ray points
      ! jsg(i) - segment number belonging to ray node i
      ! rseg,ktseg,kdwn - ray description (see subroutine rdray)

      ! output
      ! hmf - M11 and M22 for forward Hessian (A47)
      ! hmb - M11 and M22 for backward Hessian (A48)
      !  (e.g. hmf(2,j) is M22 forward in node j).
      ! slowness p (s/rad)
      ! rayvel - array with model velocity at ray node, <0 if Vs

c     double precision ypq,ynew
      dimension y(4,9000),jsg(9000),rayvel(9000),dv(9000),rayq(9000)
      dimension hmf(2,9000),hmb(2,9000),ypq(10,9000)
      dimension rseg(20),ktseg(20),kdwn(20),ynew(10)
      integer nray

      ! background model
      common /modl/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &      Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      data r2d/57.295779513082/,pi/3.14159265359/
      data halfpi/1.570796326795/,eps/0.0/
      ! TODO: safer to set eps=1e-4 etc, but may run into problem when the turning point
      ! is very close to the 1-D model node
      data r0/20./  ! radius of homogeneous sphere, the same r0 as in raymatrix
      data err0/1.0/  ! error tolerance of hessian against analytical expression

      jdebug=0     ! 0 no debug; 1 check accuracy against analytical expression
                   ! 2 for more debug output
      mray=nray

      ! set initial conditions

      ! forward elements of the partition matrices P and Q
      ypq(1,1)=1.0   ! P1 (A35)
      ypq(2,1)=1.0   ! P2 (A35)
      ypq(3,1)=0.0   ! Q1 (A35)
      ypq(4,1)=0.0   ! Q2 (A35)

      !backward
      ypq(5,1)=0.0   ! P~1 (A37)
      ypq(6,1)=0.0   ! P~2 (A37)
      ypq(7,1)=1.0   ! Q~1 (A37)
      ypq(8,1)=1.0   ! Q~2 (A37)

      ! ray radius and angle with vertical
      ypq(9,1)=y(1,1)  ! r
      ypq(10,1)=y(2,1) ! i (0 for upward)

      iq=n
      iseg=1
      ktype=ktseg(1)

      ! find source velocity (assume source in the entering layer
      ! if at a discontinuity)
      call velo1(rseg(1),iq,1-kdwn(1),ktseg(1),csrc,dc,d2c)
      p=rseg(1)*sin(y(2,1))/csrc     ! slowness
      p2=p*p
      rayvel(1)=csrc
      if(ktseg(1).eq.2) rayvel(1)=-rayvel(1)
      rayq(1)=1.0/Qs(iq)
      dv(1)=dc
      if(jdebug.gt.1) write(13,*) 'hessian p,csrc,i0=',p,csrc,y(2,1)

      if(jdebug.gt.1) then
      write(13,98)
98    format(5x,'i',3x,'jsg',4x,'kd',5x,'phi',9x,'r',6x,'c',3x,'ypq1-8')
      write(13,99) 1,jsg(1),kdwn(iseg),y(3,1)*r2d,y(1,1),rayvel(1),
     &       (ypq(j,1),j=1,8)
      endif

      ! integrate with phi as independent variable, unless a
      ! discontinuity is crossed (when we apply boundary conditions)
      do i=2,mray                       ! step from i-1 to i
        jdown=1
        if(jsg(i).lt.0) jdown=0
        dphi=y(3,i)-y(3,i-1)
        rr=y(1,i-1)
        dr=y(1,i)-rr
        call findiq(iq,jdown,rr)
        drnode=abs(rr-r(iq))    ! distance from node
        
        ! check if node i-1 and i are at same radius of a discontinuity
        if(abs(dr).le.eps.and.drnode.le.eps) then    
          rr2=rr*rr
          csi=cos(y(2,i-1))
          if(jdebug.gt.1) write(13,*) 'Discont i,jsg=',i,jsg(i-1),jsg(i)
          ! c is velocity of incoming ray, c1 of outgoing ray at
          ! this discontinuity
          kupr=jdown    ! case of transmitted ray
          if(jsg(i)*jsg(i-1).lt.0) kupr=1-jdown       ! reflection
          call velo1(rr,iq,kupr,ktype,c,dc,d2c)
          if(jdebug.gt.1)
     &    write(13,*) i-1,' Incoming rr,c=',rr,c,' kupr=',kupr

          ! debug - remove
          if(c.lt.0.2) then
            print *,'rr=',rr,' iq=',iq,' jdown=',jdown
            print *,'ktype=',ktype,' c=',c
            print *,'i=',i,'yi,yi-1=',y(1,i),y(1,i-1)
            stop
          endif  

          ktype1=ktseg(abs(jsg(i)))     ! type may change at point i
          if(jsg(i)*jsg(i-1).gt.0) kupr=1-kupr        ! transmission
          call velo1(rr,iq,kupr,ktype1,c1,dc1,d2c)

          ! debug - remove
          if(c1.lt.0.2) then
            print *,'rr=',rr,' iq=',iq,' jdown=',jdown
            print *,'ktype=',ktype,' c1=',c
            print *,'i=',i,'yi,yi-1=',y(1,i),y(1,i-1)
            stop
          endif  

          rayvel(i)=c1
          if(ktype1.eq.2) rayvel(i)=-rayvel(i)
          rayq(i)=1.0/Qs(iq)
          dv(i)=dc1
          sini=p*c1/rr
          if(jdebug.gt.1)
     &    write(13,*) i,' Outgoing rr,c1=',rr,c1,' kupr=',kupr

          ! The following uses (A32)-(A34). To re-derive subsitute
          ! (A34) first
c         if(abs(sini).ge.1.0.or.kdwn(abs(jsg(i))).eq.3) then ! reflect
c           if(jdebug.gt.1)
c    &      write(13,*) '(reflect, kdwn,sini=',kdwn(abs(jsg(i))),
c    &          sini,')'
c           csi1=cos(y(2,i))
c           ypq(1,i)=-ypq(1,i-1)
c    &           +ypq(3,i-1)*((1.0/(csi*c1)-1.0/(csi1*c))/rr
c    &           -p2*(dc1/csi-dc/csi1)/rr2)
c           ypq(2,i)=ypq(2,i-1)+ypq(4,i-1)*(csi1/c1-csi/c)/rr
c           ypq(3,i)=-ypq(3,i-1)
c           ypq(4,i)=ypq(4,i-1)
c           ypq(5,i)=-ypq(5,i-1)
c    &           +ypq(7,i-1)*((1.0/(csi*c1)-1.0/(csi1*c))/rr
c    &           -p2*(dc1/csi-dc/csi1)/rr2)
c           ypq(6,i)=ypq(6,i-1)+ypq(8,i-1)*(csi1/c1-csi/c)/rr
c           ypq(7,i)=-ypq(7,i-1)
c           ypq(8,i)=ypq(8,i-1)
c           ypq(9,i)=y(1,i)
c           ypq(10,i)=y(2,i)
c         else                                          ! transmission
            if(jdebug.gt.1)
     &      write(13,*) '(refl/transmit, kdwn,sini=',kdwn(abs(jsg(i))),
     &            sini,')'
            csi1=cos(y(2,i))
            ratcs=csi1/csi
            ypq(1,i)=ypq(1,i-1)/ratcs
     &           +ypq(3,i-1)*((1.0/(csi*c1)-1.0/(csi1*c))/rr    ! (A32)
     &           -p2*(dc1/csi-dc/csi1)/rr2)
            ypq(2,i)=ypq(2,i-1)+ypq(4,i-1)*(csi1/c1-csi/c)/rr   ! (A33)
            ypq(3,i)=ypq(3,i-1)*ratcs                           ! (A32)
            ypq(4,i)=ypq(4,i-1)                                 ! (A33)
            ypq(5,i)=ypq(5,i-1)/ratcs                           
     &           +ypq(7,i-1)*((1.0/(csi*c1)-1.0/(csi1*c))/rr
     &           -p2*(dc1/csi-dc/csi1)/rr2)                     ! (A32)
            ypq(6,i)=ypq(6,i-1)+ypq(8,i-1)*(csi1/c1-csi/c)/rr   ! (A33)
            ypq(7,i)=ypq(7,i-1)*ratcs                           ! (A32)
            ypq(8,i)=ypq(8,i-1)                                 ! (A33)
            ypq(9,i)=y(1,i)
            ypq(10,i)=y(2,i)
c         endif
          if(jdebug.gt.1)
     &    write(13,99) i,jsg(i),kdwn(iseg),y(3,i)*r2d,y(1,i),rayvel(i),
     &       (ypq(j,i),j=1,8)
        else
          rr=y(1,i)
          ktype1=ktseg(abs(jsg(i)))     ! type may change at point i
          call velo1(rr,iq,jdown,ktype1,c1,dc1,d2c)
          rayvel(i)=c1
          if(ktype1.eq.2) rayvel(i)=-rayvel(i)
          rayq(i)=1.0/Qs(iq)
          dv(i)=dc1
          call rk4pq(ypq(1,i-1),iq,dphi,ktype,jdown,ynew,p)
          do j=1,10
            ypq(j,i)=ynew(j)
          enddo  

          ! debug
          if(jdebug.gt.1)
     &    write(13,99) i,jsg(i),kdwn(iseg),y(3,i)*r2d,y(1,i),rayvel(i),
     &       (ypq(j,i),j=1,8)
99        format(3i6,f8.2,f10.2,f7.3,8e12.3)

          ypq(9,i)=y(1,i)
          ypq(10,i)=y(2,i)

        endif

        iseg=abs(jsg(i))           ! segment number for next step
        ktype=ktseg(iseg)          ! wave type for next step
        rstop=rseg(iseg)           ! r at end of segment

      end do

      q1tr=ypq(7,mray)
      q2tr=ypq(8,mray)
      q1fr=ypq(3,mray)
      q2fr=ypq(4,mray)
      if(jdebug.gt.1) write(13,*) 'qifr etc.:',q1fr,q2fr,q1tr,q2tr,mray

      ! constants used in the following do loop of raynodes
      y3d=y(3,nray)
      sid=sin(y3d)
      rcosis=y(1,1)*cos(y(2,1))
      rcosisb=-y(1,mray)*cos(y(2,mray))
      xs=y(1,1)*cos(y(3,1))    ! source coordinates in equatorial plane
      ys=y(1,1)*sin(y(3,1))    ! x=r*cos phi,y=r*sin phi
      xr=y(1,mray)*cos(y(3,mray))   ! receiver coordinates in equatorial plane
      yr=y(1,mray)*sin(y(3,mray))    ! x=r*cos phi,y=r*sin phi

      do k=1,mray

c replace ypq(5-8,k) as backward P1,P2,Q1,Q2 see eqs (A43) to (A46)
        ypq(5,k)=ypq(1,k)*q1tr-ypq(5,k)*q1fr         ! P3
        ypq(6,k)=ypq(2,k)*q2tr-ypq(6,k)*q2fr         ! P4
        ypq(7,k)=ypq(3,k)*q1tr-ypq(7,k)*q1fr         ! -Q3
        ypq(8,k)=ypq(4,k)*q2tr-ypq(8,k)*q2fr         ! -Q4

c forward hessian M11=P1/Q1 and M22=P2/Q2  
        if (k.eq.1.or.ypq(3,k).eq.0.0) then          
          hmf(1,k)=sign(999.,ypq(1,k))
        else
          hmf(1,k)=ypq(1,k)/ypq(3,k)        ! M'1=P1/Q1
        endif

        if (k.eq.1.or.ypq(4,k).eq.0.0) then
          hmf(2,k)=sign(999.,ypq(2,k))
        else
          hmf(2,k)=ypq(2,k)/ypq(4,k)       ! M'2=P2/Q2
        endif

        ! get distance between ray node and source
        x0=y(1,k)*cos(y(3,k))    ! ray node coordinates
        y0=y(1,k)*sin(y(3,k))
        diss=sqrt((x0-xs)**2+(y0-ys)**2)
        ! near source. 
        if(k.ne.1.and.diss.lt.r0) then  
          ss=0.  ! ray length from source
          vela=0.  ! average velocity
          do kk=2,k
            ds=abs((y(1,kk-1)-y(1,kk))/cos(y(2,kk-1)))
            ss=ss+ds
            vela=vela+(rayvel(kk)+rayvel(kk-1))/2.*ds
          enddo
          vela=abs(vela/ss)
          hmflim=1./ss/vela
          if(abs(hmf(1,k)).gt.hmflim) hmf(1,k)=hmflim  ! Q1~0 gives trouble
          if(abs(hmf(2,k)).gt.hmflim) hmf(2,k)=hmflim  ! Q2~0 gives trouble
        endif

c backward hessian M11=P3/Q3 and M22=P4/Q4 
        if (k.eq.mray.or.ypq(7,k).eq.0.0) then
          hmb(1,k)=sign(999.,ypq(5,k))
        else
          hmb(1,k)=-ypq(5,k)/ypq(7,k)       ! M"1=P3/Q3
        endif

        if (k.eq.mray.or.ypq(8,k).eq.0.0) then
          hmb(2,k)=sign(999.0,ypq(6,k))
        else
          hmb(2,k)=-ypq(6,k)/ypq(8,k)       ! M"2=P4/Q4
        endif                  

        ! get distance between ray node and receiver
        disr=sqrt((x0-xr)**2+(y0-yr)**2)
        ! near receiver
        if(k.ne.mray.and.disr.lt.r0) then  
          sr=0.  ! ray length from receiver
          vela=0.
          do kk=k,mray-1
            ds=abs((y(1,kk)-y(1,kk+1))/cos(y(2,kk+1)))
            sr=sr+ds
            vela=vela+(rayvel(kk)+rayvel(kk+1))/2.*ds
          enddo
          vela=abs(vela/sr)
          hmblim=1./sr/vela
          if(abs(hmb(1,k)).gt.hmblim) hmb(1,k)=hmblim  ! Q3~0 gives trouble
          if(abs(hmb(2,k)).gt.hmblim) hmb(2,k)=hmblim  ! Q4~0 gives trouble
        endif

        ! debug output - remove
c       write(14,11) k,hmf(1,k)+hmb(1,k),hmf(1,k),hmb(1,k),ypq(1,k),
c    &    ypq(3,k),ypq(5,k),-ypq(7,k),y(3,k)*r2d,y(1,k)
c       write(15,11) k,hmf(2,k)+hmb(2,k),hmf(2,k),hmb(2,k),ypq(2,k),
c    &    ypq(4,k),ypq(6,k),-ypq(8,k),y(3,k)*r2d,y(1,k)
c11      format(i3,7e12.4,f8.3,f7.1)

        ! check for accuracy against analytical expressions
        ! for nodes where Q1,Q2,Q1~,Q2~ are not zero
        if(jdebug.gt.0 .and. k.gt.1.and.k.lt.mray .and.
     &    ypq(3,k).ne.0.0.and.ypq(4,k).ne.0.0 .and.
     &    ypq(7,k).ne.0.0.and.ypq(8,k).ne.0.0) then  
c    &    diss.ge.r0.and.disr.ge.r0) then
          sip=sin(y(3,k))
          sipi=sin(y(3,k)+y(2,k))
          vk=abs(rayvel(k))
          rr=y(1,k)
          sii=sin(y(2,k))
          csi=cos(y(2,k))
          tai=sii/csi
          cti=csi/sii
          a53=rr*csi*ypq(1,k)+(p*p*dv(k)/rr-1.0/vk)*ypq(3,k) ! (A53)
          a53b=-rr*csi*ypq(5,k)-(p*p*dv(k)/rr-1.0/vk)*ypq(7,k) ! (A53) backward
          am2=sipi/(vk*rr*sip)                               ! (A50) 
          ap2=rseg(1)*sipi/(p*vk)                            ! (A49) 
          aq2=rr*rseg(1)*sip/p                               ! (A49) 
          if(jdebug.gt.1) then
            write(13,92) 'node ',k,' radius',rr,' phi',y(3,k)*r2d
            write(13,90) "theor M'2,P2,Q2:",am2,ap2,aq2
            write(13,90) "comp  M'2,P2,Q2:",hmf(2,k),ypq(2,k),ypq(4,k)
90          format(a,3e15.6)
          endif
          err=(hmf(2,k)-am2)/am2*100.
          if(abs(err).gt.err0) write(*,93) "WARNING: M'2 imprecision ",
     &      'for ray node ',k,' err=',err,'%'
          sidp=sin(y3d-y(3,k))
          summ=p*sid/(rr*rr*sip*sidp)                        ! (A51) 
          am2b=summ-am2
          if(jdebug.gt.1)
     &      write(13,91) 'theor M"2:',am2b,', comp:',hmb(2,k)
          err=(hmb(2,k)-am2b)/am2b*100.
          if(abs(err).gt.err0) write(*,93) 'WARNING: M"2 imprecision',
     &      ' for ray node ',k,' err=',err,'%'
          if(jdebug.gt.1)
     &      write(13,91) 'analytical (A53):',a53,', comp:',rcosis
          err=(a53-rcosis)/rcosis*100.
          if(abs(err).gt.err0) write(*,93) "WARNING: A53' imprecision",
     &      ' for ray node ',k,' err=',err,'%'
          err=(a53b-rcosisb)/rcosisb*100.
          if(abs(err).gt.err0) write(*,93) 'WARNING: A53" imprecision',
     &      ' for ray node ',k,' err=',err,'%'
91        format(a,e15.6,a,e15.6)
92        format(a,i5,a,f9.2,a,f9.3)
93        format(2a,i4,a,f6.1,a)
        endif

      end do

      return
      end 

      subroutine tracer(ystart,rseg,ktseg,kdwn,nseg,y,jsg,rmin,
     &      tstar,h,nray)

      ! traces ray specified by starting conditions in ystart and
      ! reflection/transmission/turning behaviour in rseg etc.

      ! input
      ! ystart(1)=r, ystart(2)=angle i, ystart(3)=delta, ystart(4)=time
      !        of the first ray point (usually 3 and 4 are zero)
      ! rseg,ktseg,kdwn, for i=1,...,nseg describe nature of ray (see
      !        main program)
      ! h=integration step in km for Runge Kutta algorithm

      ! output
      ! y(1-4,i), i=1,...,nray: r,angle,delta,time of ray
      ! jsg(i), i=1,...,nray: segment number of node i, with - sign if
      !          the ray travels in the upward direction
      ! nray=number of ray nodes in y (<= 0 if error, eg when
      !      ray depth exceeds turning point depth in rseg)
      ! tstar=traveltime/qray

      ! "error" messages:
      ! nray= 0: ray turns before reaching deeper segment
      !      -1: idem, but at discontinuity (not segment boundary)
      !      -2: reaches deepest segment level before turning
      !      -3: downward transmitted ray becomes evanescent at disc
      !      -4: upward transmitted ray becomes evanescent at discont
      !      -5: upward transmitted ray reflects (not segment boundary)
      !      -6: ray get closer than 10 km to Earth's center

      ! debug output if jdebug=0

      dimension rseg(20),ktseg(20),kdwn(20)
      dimension ystart(4),y(4,9000),dydx(4),jsg(9000)
      ! background model
      common /modl/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &      Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      data r2d/57.295779513082/,pi/3.14159265359/
      data halfpi/1.570796326795/
      data epsdr/0.1/           ! discontinuity tolerance

      jdebug=0
c     modified statement to reduce debug output:
c     if(ystart(2).gt.1.618.and.ystart(2).lt.1.619) jdebug=0      

      if(jdebug.eq.1) write(13,*) 'call to tracer with ystart=',ystart
      do i=1,4
        y(i,1)=ystart(i)
      enddo
      iq=n
      call velo1(rseg(1),iq,1-kdwn(1),ktseg(1),csrc,dc,d2c)
      p=rseg(1)*sin(ystart(2))/csrc     ! slowness
      iseg=1
      nray=1
      iqstp=n
      jdown=kdwn(1)           ! 1 if ray heads downwards, 0 up
      rmin=y(1,1)
      jsg(1)=1
      if(jdown.eq.0) jsg(1)=-1
      epstp=0.                ! makes sure we pass a discontinuity

      if(jdebug.gt.0)
     &write(13,98)
98    format(2x,'nray',9x,'r',9x,'i',7x,'del',9x,'t jd',8x,'dr',5x,'iq')
      if(jdebug.gt.0)
     &write(13,99) nray,y(1,nray),y(2,nray)*r2d,y(3,nray)*r2d,
     &        y(4,nray),jdown,0,iq,jsg(nray)

      ! the following do loop fills y for node nray+1:
      do while(iseg.lt.nseg)
        ! don't trespass rstop as determined by rstp:
        rstart=y(1,nray)-jdown*epstp
        call rstp(rstart,iqstp,jdown,rseg,iseg,rstop,jsegsw)   
        ! step h in the ray direction
        call rk4l(y(1,nray),iq,h,ktseg(iseg),jdown,y(1,nray+1)) 
        dr=y(1,nray+1)-rstop    ! -distance to rstop
        if(jdown.eq.1) dr=-dr
        jsg(nray+1)=iseg
        if(jdown.eq.0) jsg(nray+1)=-iseg
        if(jdebug.gt.0)
     &  write(13,99) nray+1,y(1,nray+1),y(2,nray+1)*r2d,y(3,nray+1)*r2d,
     &        y(4,nray+1),jdown,dr,iq,jsg(nray)
99      format(i6,4f10.3,i3,f10.3,2i7)
        ! debug- remove
!       if(abs(y(1,nray+1)-3484.615).lt.0.003) then
!         print *,'y=',y(1,nray+1)
!       endif  
        if(dr.gt.-epsdr) then       ! near or past rstop, redo to rstop exactly
          epstp=epsdr               ! makes sure it doesn't get stuck here    
          hr=rstop-y(1,nray)        ! adjusted step size    
          if(jdebug.gt.0)
     &    write(13,*) 'call to rk4r,iqstp,rstop,hr=',iqstp,rstop,hr
          ! step hr in r direction
          call rk4r(y(1,nray),iq,hr,ktseg(iseg),jdown,y(1,nray+1)) 
          if(jdebug.gt.0)
     &    write(13,99) nray+1,y(1,nray+1),y(2,nray+1)*r2d,
     &          y(3,nray+1)*r2d,y(4,nray+1),jdown,hr,iq,jsg(nray)

          if(y(1,nray).lt.10.0) then
            nray=-6
            if(jdebug.gt.0) write(13,*) 'Ray hits Earth center'
            return
          endif  
          ! apply continuity conditions

          if(jsegsw.eq.0.) then         ! no segment change

            nray=nray+2
            if(nray.gt.8998) then
              print *,'Ray too long for slowness ',p
              print *,'ystart(2)=',ystart(2),' rad'
              print *,'Last distance ',y(3,nray-1)*r2d
              stop 'nray>8998'
            endif
            y(1,nray)=y(1,nray-1)       ! r continuous
            y(3,nray)=y(3,nray-1)       ! delta continuous
            y(4,nray)=y(4,nray-1)       ! time continuous
            jsg(nray)=iseg
            if(jdown.eq.0) jsg(nray)=-iseg
            call findiq(iq,1-jdown,rstop)  ! find lower node if downwards
            if(jdown.eq.1) then
              sini=p*vps(iq,ktseg(iseg))/rstop
              if(sini.ge.1) then        ! turning point
                iseg=iseg+1
                if(kdwn(iseg).ne.2) then        ! turning unexpected
                  nray=-1
                  return
                else
                  y(2,nray)=pi-y(2,nray-1)
                  jdown=0
                  jsg(nray)=-abs(jsg(nray))
                endif  
              else
                y(2,nray)=pi-asin(sini)   ! i=0 for vertical ray UP
              endif  
            else        ! case of upward ray
              sini=p*vps(iq,ktseg(iseg))/rstop
              if(abs(sini).ge.1.0) then   ! reflects back down!
                nray=-5
                return
              else  
                y(2,nray)=asin(sini)
              endif  
            endif
            if(jdebug.gt.0)
     &      write(13,99) nray,y(1,nray),y(2,nray)*r2d,
     &          y(3,nray)*r2d,y(4,nray),jdown,dr,iq,jsg(nray)

          else                   ! careful, next segment reached

            iseg=iseg+1
            if(jdebug.gt.0)
     &      write(13,*) 'new segment',iseg,rseg(iseg),kdwn(iseg),
     &                  ktseg(iseg)
            if(kdwn(iseg).eq.5) then            ! end of ray
              nray=nray+1
              jsg(nray)=iseg
              if(jdown.eq.0) jsg(nray)=-iseg
              if(nray.gt.8998) then
                print *,'ystart(2)=',ystart(2),' rad'
                print *,'Ray too long for slowness ',p
                print *,'Last distance ',y(3,nray-1)*r2d
                stop 'A nray>8998'
              endif
              if(jdebug.gt.0)
     &        write(13,*) 'end of ray nrry=',nray
              jdebug=0
              call gettstar(y,nray,jsg,ktseg,tstar,qray)
              return
            else if(kdwn(iseg).eq.2) then       ! turning point
              nray=-2
              if(jdebug.gt.0)
     &        write(13,*) 'turning point at disc, nray=-2'
              return
            else if(kdwn(iseg).eq.3) then       ! reflection
              nray=nray+2
              if(nray.gt.8998) then
                print *,'ystart(2)=',ystart(2),' rad'
                print *,'Ray too long for slowness ',p
                print *,'Last distance ',y(3,nray-1)*r2d
                stop 'nray>8998'
              endif
              y(1,nray)=y(1,nray-1)     ! r continuous
              y(3,nray)=y(3,nray-1)     ! delta continuous
              y(4,nray)=y(4,nray-1)     ! time continuous
              jdown=1-jdown             ! reverse direction
              jsg(nray)=iseg
              if(jdown.eq.0) jsg(nray)=-iseg
              if(ktseg(iseg-1).eq.ktseg(iseg)) then   ! if no P/S conversion
                y(2,nray)=pi-y(2,nray-1)
              else
                sinin=sin(y(2,nray-1))
                call findiq(iq,1-jdown,rstop)
                vratio=vps(iq,ktseg(iseg-1))/vps(iq,ktseg(iseg))
                sinout=sinin/vratio
                if(abs(sinout).ge.1.) then        ! evanescent wave
                  nray=-3
                  if(jdebug.gt.0)
     &            write(13,*) 'evanescent wave'
                  return
                endif  
                y(2,nray)=asin(sinout)
                if(jdown.eq.1) y(2,nray)=pi-y(2,nray)
              endif
              if(jdebug.gt.0)
     &        write(13,99) nray,y(1,nray),y(2,nray)*r2d,
     &          y(3,nray)*r2d,y(4,nray),jdown,dr,iq,jsg(nray)
            else                                ! transmission
              nray=nray+2
              if(nray.gt.8998) then
                print *,'ystart(2)=',ystart(2),' rad'
                print *,'Ray too long for slowness ',p
                print *,'Last distance ',y(3,nray-1)*r2d
                stop 'nray>8998'
              endif
              y(1,nray)=y(1,nray-1)     ! r continuous
              y(3,nray)=y(3,nray-1)     ! delta continuous
              y(4,nray)=y(4,nray-1)     ! time continuous
              jsg(nray)=iseg
              if(jdown.eq.0) jsg(nray)=-iseg
              call findiq(iq,1-jdown,rstop)
              if(jdown.eq.1) then
                sini=p*vps(iq,ktseg(iseg))/rstop
                if(abs(sini).ge.1.) then        ! evanescent wave
                  nray=-3
                  if(jdebug.gt.0)
     &            write(13,*) 'evanescent wave'
                  return
                endif  
                y(2,nray)=pi-asin(sini) ! i=0 for vertical ray UP
              else      ! case upward traveling ray
                sini=p*vps(iq,ktseg(iseg))/rstop
                if(abs(sini).ge.1.) then        ! evanescent wave
                  nray=-4
                  if(jdebug.gt.0)
     &            write(13,*) 'evanescent wave'
                  return
                endif  
                y(2,nray)=asin(sini)
                if(jdebug.gt.0)
     &          write(13,99), nray,y(1,nray),y(2,nray)*r2d,
     &            y(3,nray)*r2d,y(4,nray),jdown,dr,iq,jsg(nray)
              endif             ! for (jdown.eq.1) fork
            endif               ! for (kdwn.eq.5) fork
          endif                 ! for (jsegsw.eq.0) fork

        ! we did not cross rstop (fork on dr.ge.0)
        else 
          nray=nray+1
          epstp=0.
          if(nray.gt.8998) then
            print *,'ystart(2)=',ystart(2),' rad'
            print *,'Ray too long for slowness ',p
            print *,'Last distance ',y(3,nray-1)*r2d
            stop 'B nray>8998'
          endif

          ! finally, check if we passed the turning point
          if(jdown.eq.1.and.y(2,nray).lt.halfpi) then
            jdown=0
            iseg=iseg+1
            if(jdebug.gt.0)
     &      write(13,*) 'TP: new segment',iseg,rseg(iseg),kdwn(iseg),
     &                  ktseg(iseg)
            if(kdwn(iseg).ne.2) then  ! ray was not supposed to turn!
              nray=0
              if(jdebug.gt.0)
     &        write(13,*) 'Ray was not supposed to turn, nray=0'
              return
            endif  
          endif                 ! for turning point fork
          jsg(nray)=iseg
          if(jdown.eq.0) jsg(nray)=-iseg
        endif                   ! for (dr.ge.0) fork
        rmin=min(rmin,y(1,nray))
      end do                    ! for (iseg.lt.nseg) loop

      if(jdebug.gt.0)
     &write(13,*) 'Hmmm - you should never have reached this....'
      stop 'software error in tracer'

      return
      end

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
      !           4=transmission, 5=last point of ray) from node i
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

      common /modl/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &      Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      dimension rseg(20),ktseg(20),kdwn(20)
      character*72 chray
      character*1 wtyp(0:2)
      character*4 ktyp(0:5)
      character*16 phase

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

      ierror=0
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
      if(ktseg(i).lt.1.or.ktseg(i).gt.2) stop 'wave type must be 1 or 2'
      if(kdwn(1).gt.1) then
        print *,'First ray point should have kdwn 1 or 0'
        if(kout.ne.0) write(kout,*) 
     &        'First ray point should have kdwn 1 or 0'
        ierror=ierror+1
      endif  
      rmin=rseg(1)
      rmax=rmin
      do while(kdwn(i).ne.5)
        i=i+1
        if(i.gt.20) stop 'Max segments is 20'
        read(kin,*) rseg(i),kdwn(i),ktseg(i)
        if(ktseg(i).eq.0) ktseg(i)=ktseg(i-1)
        if(ktseg(i).lt.1.or.ktseg(i).gt.2) stop 'type must be 1 or 2'
        if(kout.ne.0) 
     &        write(kout,40) rseg(i),kdwn(i),ktseg(i),ktyp(kdwn(i)),
     &        wtyp(ktseg(i))
        if(kdwn(i).lt.2) then
          print *,'kdwn < 2 only allowed for first segment'
          if(kout.ne.0) write(kout,*) 
     &          'kdwn < 2 only allowed for first segment'
          ierror=ierror+1
        endif  
        if(kdwn(i).eq.3.or.kdwn(i).eq.4) then
          call findiq(iq,1,rseg(i))
          if(iq.lt.2.or.iq.gt.n) then
            print *,'Radius outside of model, segment',i
            if(kout.ne.0) write(kout,*) 'Radius outside of model'
            ierror=ierror+1
          else if(r(iq).lt.6370.0.and.abs(r(iq)-r(iq-1)).gt.0.01.or.
     &      abs(rseg(i)-r(iq)).gt.0.01) then
            print *,'Refl/transmission not at discontinuity, segment',i
            if(kout.ne.0) write(kout,*) 
     &            'Refl/transmission not at discontinuity'
            ierror=ierror+1
          endif
        endif  
            
        rmin=min(rmin,rseg(i))
        rmax=max(rmax,rseg(i))
        if(jdebug.eq.1) print *,rseg(i),kdwn(i),ktseg(i),ktyp(kdwn(i)),
     &        wtyp(ktseg(i))
      end do
      if(jdebug.eq.1) print *,'nseg=',i
      nseg=i
      print *,'Rays between r=',rmin,rmax

      if(ierror.gt.0) then
        print *,'Stop because of',ierror,' errors in ray definition'
        stop
      endif  

      return
      end

      subroutine rdband(kin,kout,bndomega,dotm,nfreq,nband)

      ! reads spectral band info
      ! input: I/O unit numbers for input and output kin,kout
      !        (assumed open, no output if kout=0).
      ! output: bndomega=frequencies (rad/s), dotm=mdot(omega), nfreq=# of frequencs

      parameter (NFR=16)
      dimension bndomega(200,NFR),dotm(200,NFR),nfreq(NFR)

      read(kin,*) nband
      if(nband.gt.NFR) stop 'nband > NFR, increase dimensions'
      if(kout.ne.0) write(kout,fmt='(i5)') nband
      do ib=1,nband
        read(kin,*) nfreq(ib)
        if(nfreq(ib).gt.200) stop 'nfreq>200, increase dimensions'
        if(kout.ne.0) write(kout,fmt='(i5)') nfreq(ib)
        do j=1,nfreq(ib)
          read(kin,*) bndomega(j,ib),dotm(j,ib)
          if(kout.ne.0) 
     &          write(kout,fmt='(2f8.4)') bndomega(j,ib),dotm(j,ib)
        end do
      end do

      return
      end


      subroutine rstp(rstart,iq,jdown,rseg,iseg,rstop,jsegsw)   
      dimension rseg(20),ktseg(20),kdwn(20)
      ! background model
      common /modl/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &      Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      data eps/0.01/

      ! determines if discontinuity or segment will be crossed
      ! or if end of segment depth is reached

      ! input: rstart=current radius, 
      !        jdown=direction (1 down, else 0)
      !        iq=guess for layer number belonging to rr,
      !        rseg=segment start radii, iseg=current segment #
      ! output: rstop=r limit
      !        jsegsw=0 if simply a discontinuity (ie transmission)
      !        jsegsw=1 if the r limit is also a segment boundary

      rr=rstart
      call findiq(iq,jdown,rr)  ! r(iq) is directly below (or =) rr
c     if(abs(rr-r(iq)).le.eps) then
c       rr=rr+(2-4*jdown)*eps   ! move 2*eps away from discontinuity
c       call findiq(iq,jdown,rr)
c     endif 
      rs=rseg(iseg+1)           ! start r of next segment
      i=iq
      ! find first discontinuity before r(i) reaches rs:
      if(jdown.eq.1) then       ! downward ray
        do while((i.gt.1).and.r(i)-r(i-1).gt.eps.and.r(i)-rs.gt.0.)
c       do while((i.gt.1).and.r(i)-r(i-1).gt.eps.and.r(i)-rs.gt.eps)
          i=i-1
        enddo  
      else  
        i=min(n,iq+1)
        do while((i.lt.n).and.r(i+1)-r(i).gt.eps.and.rs-r(i).gt.eps)
          i=i+1
        enddo  
      endif  
      jsegsw=0
      rstop=r(i)
      if(jdown.eq.1.and.rstop.lt.rs) rstop=rs
      if(jdown.eq.0.and.rstop.gt.rs) rstop=rs
      if(abs(rstop-rs).lt.eps) jsegsw=1

      return
      end

      subroutine velo1(rr,iq,kupr,ktype,c,dc,d2c)

      ! input: 
      ! rr=radius 
      ! iq is guess for model node directly beneath or at rr
      !    (warning: iq is also output and may be changed, so do 
      !    never use <1> or a common block variable like <n> or <noc> etc!)
      ! kupr=1 if the upper layer of a discontinuity must be selected,
      !      0 if the lower one. kupr is ignored if r is not on a
      !      discontinuity.
      ! ktype=1 for P, 2 for S wave
      ! common block modl/ must be filled (see subroutine model.f)

      ! output: 
      ! iq= node # directly below or at rr 
      !     (is modified if wrong on input!)
      ! c=velocity at radius rr
      ! dc= dc/dr
      ! d2c = d2c/dr2 (second derivative)

      ! method: spline interpolation


      ! background model
      common /modl/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &      Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      data eps/0.01/            ! precision 10 m

      ! find node immediately below r
      call findiq(iq,kupr,rr)

      t=rr-r(iq)                ! distance from node
      if(ktype.eq.1) then
        c=vp(iq)+t*(qvp(1,iq)+t*(qvp(2,iq)+t*qvp(3,iq)))
        dc=qvp(1,iq)+2.*qvp(2,iq)*t+3*qvp(3,iq)*t*t
        d2c=2*qvp(2,iq)+6*t*qvp(3,iq)
      else  
        c=vs(iq)+t*(qvs(1,iq)+t*(qvs(2,iq)+t*qvs(3,iq)))
        if(c.le.0.) then        ! avoid problem in tracer
          c=0.01
          dc=0
          d2c=0
          return
        endif
        dc=qvs(1,iq)+2.*qvs(2,iq)*t+3*qvs(3,iq)*t*t
        d2c=2*qvs(2,iq)+6*t*qvs(3,iq)
      endif  

      return
      end

      subroutine findiq(iq,kupr,rr)

      ! input: iq=guess for layer number
      !        kupr=1 if 'upper' node of discontinuity preferred, else 0
      !        rr=radius
      ! output: iq= layer number at or directly below rr. Ambiguity if rr is at
      !        a discontinuity is decided by kupr


      ! background model
      common /modl/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &      Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      data eps/0.01/

      ! find node immediately below rr
      iq=min(n,max(iq,1))
      do while(iq.gt.1.and.r(iq).gt.rr+eps)
        iq=iq-1
      enddo
      do while(iq.lt.n.and.r(iq+1).lt.rr+eps)
        iq=iq+1
      enddo  

      if(rr-r(iq).ge.eps) return

      ! resolve discontinuity ambiguitiy
      if(iq.lt.n.and.abs(r(iq+1)-r(iq)).lt.eps.and.kupr.eq.1) 
     &          iq=iq+1
      if(iq.gt.1.and.abs(r(iq)-r(iq-1)).lt.eps.and.kupr.eq.0)
     &          iq=iq-1 

      return
      end

      function vps(i,ktype)
      ! returns vp at node i if ktype=1, else vs
      ! background model
      common /modl/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &      Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      if(ktype.eq.1) then
        vps=vp(i)
      else
        vps=vs(i)
      endif
      return
      end

      subroutine model(fname)

      ! Reads a model from file <fname>
      ! Model units are assumed to be MKS, unless the first density
      ! read is less than 100, in which case km and km/s units are
      ! assumed
      ! This is the same file format as used for love.f and rayleigh.f

      ! Format of the model file [unit 1, opened/closed by the routine]:
      ! Line 1: ignored
      ! Line 2: ignored
      ! Line 3: n (# of nodes), nic (top node of inner core), noc (top
      !         node of outer core), nmoh (bottom[!] node of the crust)
      ! Lines 4ff: radius, density, Vp, Vs

      ! The nodes are ordered with INCREASING radius
      ! Discontinuities are indicated by two subsequent nodes
      ! with the same radius. Nodes closer than 0.01 km are collapsed
      ! into a discontinuity.

      ! input: model file name fname (character*72)
      ! output: common block modl/ is filled with model information
      !         including spline coefficients.

      character*72 fname
      dimension wrk(3,500)
      ! background model
      common /modl/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &      Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      jdebug=0

      nmoh=0            ! safeguard against old model files
      open(1,file=fname)
      read(1,*) 
      read(1,*) 
      read(1,*) n,nic,noc,nmoh
      if(nmoh.le.noc.or.nmoh.gt.n) stop 'model: Error in Moho index'
      mks=1
      read(1,*) r(1),rho,vp(1),vs(1)
      if(rho.lt.100.) mks=0
      backspace(1)
      do i=1,n
        read(1,*) r(i),rho,vp(i),vs(i),Qs(i)
        if(mks.eq.1) then
          r(i)=0.001*r(i)
          vp(i)=0.001*vp(i)
          vs(i)=0.001*vs(i)
        endif
        if(jdebug.gt.0) write(13,*) i,r(i),vp(i),vs(i),Qs(i)
      end do  
      close(1)

   60 nsl=n   
      if(vs(nsl).gt.0.) go to 70
   65 nsl=nsl-1
      if(vs(nsl).le.0.) go to 65
   70 nicp1=nic+1
      nocp1=noc+1
      nslp1=nsl+1

      if(nsl.ne.n) stop 'Oceanic reference model not allowed'

c*** spline ***
      rn=r(n)
!!    do 45 i=1,n
!!    if(i.gt.1.and.abs(r(i)-r(i-1)).lt.0.01) r(i)=r(i-1)
! Intel gives the error message: Subscript #1 of the array R has value 0 which is less than the lower bound of 1
! g77 doesn't have this problem
      do 45 i=2,n
      if(abs(r(i)-r(i-1)).lt.0.01) r(i)=r(i-1)
   45 continue
      call drspln(1,n,r,vp,qvp,wrk)
      call drspln(1,n,r,vs,qvs,wrk)
      call drspln(1,n,r,Qs,qqs,wrk)

      if(jdebug.eq.0) return

      write(13,*) ' i    r    qvp'
      do i=1,n
        write(13,*) i,r(i),qvp(1,i),qvp(2,i),qvp(3,i)
      enddo

! added 19/11/06 ----------------------
      if(jdebug.eq.0) return

      ! write GMT files for debugging
      open(41,file='vp.xy')
      open(42,file='dvp.xy')
      open(43,file='ddvp.xy')
      open(44,file='vs.xy')
      open(45,file='dvs.xy')
      open(46,file='ddvs.xy')

      iq=1
      do i=2,n
        step=0.1000001*(r(i)-r(i-1))
        step=max(1.0,step)              ! avoid step=0
        do rr=r(i-1),r(i),step
          call velo1(rr,iq,0,1,c,dc,d2c)        ! P velocity
          write(41,50) rr,c
          write(42,51) rr,dc
          write(43,51) rr,d2c
          if(i.gt.noc) then
            call velo1(rr,iq,0,2,c,dc,d2c)      ! S velocity
            write(44,50) rr,c
            write(45,51) rr,dc
            write(46,51) rr,d2c
          endif
        enddo
      enddo
      close(41)
      close(42)
      close(43)
      close(44)
      close(45)
      close(46)
50    format(f8.1,f8.3)
51    format(f8.1,e12.3)
!-----------------------------------------

      return
      end        

      subroutine drspln(i1,i2,x,y,q,f)

c   rspln computes cubic spline interpolation coefficients
c   for y(x) between grid points i1 and i2 saving them in q.  the
c   interpolation is continuous with continuous first and second
c   derivitives.  it agrees exactly with y at grid points and with the
c   three point first derivitives at both end points (i1 and i2).
c   x must be monotonic but if two successive values of x are equal
c   a discontinuity is assumed and seperate interpolation is done on
c   each strictly monotonic segment.  the arrays must be dimensioned at
c   least - x(i2), y(i2), q(3,i2), and f(3,i2).  f is working storage
c   for rspln.
c                                                     -rpb
      dimension x(i2),y(i2),q(3,i2),f(3,i2),yy(3)
      equivalence (yy(1),y0)
      data yy/3*0./
      j1=i1+1
      y0=0.
c   bail out if there are less than two points total.
      if(i2-i1)13,17,8
 8    a0=x(j1-1)
c   search for discontinuities.
      do 3 i=j1,i2
      b0=a0
      a0=x(i)
      if(a0-b0)3,4,3
 3    continue
 17   j1=j1-1
      j2=i2-2
      go to 5
 4    j1=j1-1
      j2=i-3
c   see if there are enough points to interpolate (at least three).
 5    if(j2+1-j1)9,10,11
c   only two points.  use linear interpolation.
 10   j2=j2+2
      y0=(y(j2)-y(j1))/(x(j2)-x(j1))
      do 15 j=1,3
      q(j,j1)=yy(j)
 15   q(j,j2)=yy(j)
      go to 12
c   more than two points.  do spline interpolation.
 11   a0=0.
      h=x(j1+1)-x(j1)
      h2=x(j1+2)-x(j1)
      y0=h*h2*(h2-h)
      h=h*h
      h2=h2*h2
c   calculate derivitive at near end.
      b0=(y(j1)*(h-h2)+y(j1+1)*h2-y(j1+2)*h)/y0
      b1=b0
c   explicitly reduce banded matrix to an upper banded matrix.
      do 1 i=j1,j2
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h-2.*a0
      h3a=2.*h-3.*a0
      h2b=h2*b0
      q(1,i)=h2/ha
      q(2,i)=-ha/(h2a*h2)
      q(3,i)=-h*h2a/h3a
      f(1,i)=(y0-h*b0)/(h*ha)
      f(2,i)=(h2b-y0*(2.*h-a0))/(h*h2*h2a)
      f(3,i)=-(h2b-3.*y0*ha)/(h*h3a)
      a0=q(3,i)
 1    b0=f(3,i)
c   take care of last two rows.
      i=j2+1
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h*ha
      h2b=h2*b0-y0*(2.*h-a0)
      q(1,i)=h2/ha
      f(1,i)=(y0-h*b0)/h2a
      ha=x(j2)-x(i+1)
      y0=-h*ha*(ha+h)
      ha=ha*ha
c   calculate derivitive at far end.
      y0=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/y0
      q(3,i)=(y0*h2a+h2b)/(h*h2*(h-2.*a0))
      q(2,i)=f(1,i)-q(1,i)*q(3,i)
c   solve upper banded matrix by reverse iteration.
      do 2 j=j1,j2
      k=i-1
      q(1,i)=f(3,k)-q(3,k)*q(2,i)
      q(3,k)=f(2,k)-q(2,k)*q(1,i)
      q(2,k)=f(1,k)-q(1,k)*q(3,k)
 2    i=k
      q(1,i)=b1
c   fill in the last point with a linear extrapolation.
 9    j2=j2+2
      do 14 j=1,3
 14   q(j,j2)=yy(j)
c   see if this discontinuity is the last.
 12   if(j2-i2)6,13,13
c   no.  go back for more.
 6    j1=j2+2
      if(j1-i2)8,8,7
c   there is only one point left after the latest discontinuity.
 7    do 16 j=1,3
 16   q(j,i2)=yy(j)
c   fini.
 13   return
      end

c Eq numbers from Dahlen et al (GJI,2000,appendix)

      subroutine derivs(y,dydx,ir,kdis,kt,c,dc,d2c)

      ! computes derivatives of r,i,phi and t w.r.t. l

      dimension dydx(4),y(4)

      ! input:
      ! y(1)=r, y(2)=i
      ! ir: layer index
      ! kdis: 1 if ray approaches node ir from above, else 0.
      ! kt: wave type (1=P, 2=S)

      ! output
      ! dydx: derivatives of y w.r.t. ray distance dl
      ! c: velocity
      ! dc: dc/dr
      ! dc2: d^2c/dr^2

c calculate velocity and dv/dr at a given r
      call velo1(y(1),ir,kdis,kt,c,dc,d2c)

      sii=sin(y(2))
      coi=cos(y(2))
      ri=1.0/max(0.01,y(1))

      dydx(1)=coi                       ! dr/dl (A5)
      dydx(2)=sii*(dc/c-ri)             ! di/dl
      dydx(3)=sii*ri                    ! dphi/dl
      dydx(4)=1.0/c                     ! dT/dl

      return
      end

      subroutine derivr(y,dydx,ir,kdis,kt,c,dc,d2c)

      ! derivatives of r,i,phi and t w.r.t. r (not l!)

      dimension dydx(4),y(4)

      ! input:
      ! y(1)=r, y(2)=i, y(3)=phi y(4)=t
      ! ir: layer index
      ! kdis: 1 if ray approaches node ir from above, else 0.
      ! kt: wave type (1=P, 2=S)

      ! output
      ! dydx: derivatives of y w.r.t. ray distance dl
      ! c: velocity
      ! dc: dc/dr
      ! dc2: d^2c/dr^2

      ! debug - remove
c     if(y(1).lt.3484.0) then
c       print *, 'deriv is in OC!, y(1)=',y(1)
c     endif  

c calculate velocity and dv/dr at a given r
      rr=y(1)
      call velo1(rr,ir,kdis,kt,c,dc,d2c)

      tgi=tan(y(2))
      coi=cos(y(2))
      ri=1.0/max(0.01,y(1))

      dydx(1)=1.0                       ! dr/dr 
      dydx(2)=tgi*(dc/c-ri)             ! di/dr
      dydx(3)=tgi*ri                    ! dphi/dr
      dydx(4)=1.0/(c*coi)               ! dT/dr

      return
      end

      subroutine rk4l(y,iq,h,ktype,jdown,yout)

      ! integrates r,phi,i,t as function of l

      ! 4th order Runge Kutta solution to y'=f(x,y):
      ! y(n+1)=y(n)+(k1+2k2+2k3+k4)/6

      ! input:
      ! y(1)=r, y(2)=i, y(3)=phi, y(4)=t
      ! iq=guess for node # directly below y(1)
      ! h: step size (km)
      ! ktype=1 for P, 2 for S
      ! jdown=1 if step is decreasing r, 0 if upward ray

      ! output
      ! yout(1)-yout(4): next ray point
      ! iq is adjusted if needed

      dimension dydx(4),y(4),yout(4)
      dimension dym(4),dyt(4),yt(4)
      data halfpi/1.570796326795/

      jdebug = 0
      hh=h*0.5
      h6=h/6.0
      n=4

      kupr=jdown      ! choose proper layer if y(1) at discontinuity
      call derivs(y,dydx,iq,1-kupr,ktype,c,dc,d2c)
      if(jdebug.gt.0) then
        write(13,*) 'Enter rk4l, r,h=',y(1),h
        write(13,97)
        write(13,96) 'k1',y(1),y(2),dydx(1),dydx(2),dydx(3),dydx(4),c,dc
      endif
      do  i=1,n
        yt(i)=y(i)+hh*dydx(i)  
      end do

      call derivs(yt,dyt,iq,kupr,ktype,c,dc,d2c)
      if(jdebug.gt.0)
     &  write(13,96) 'k2',yt(1),yt(2),dyt(1),dyt(2),dyt(3),dyt(4),c,dc
      do i=1,n
        yt(i)=y(i)+hh*dyt(i)    
      enddo
      call derivs(yt,dym,iq,kupr,ktype,c,dc,d2c)
      if(jdebug.gt.0)
     &  write(13,96) 'k3',yt(1),yt(2),dym(1),dym(2),dym(3),dym(4),c,dc

      do i=1,n
        yt(i)=y(i)+h*dym(i)   
        dym(i)=dyt(i)+dym(i)
      end do
      call derivs(yt,dyt,iq,kupr,ktype,c,dc,d2c)
      if(jdebug.gt.0)
     &  write(13,96) 'k4',yt(1),yt(2),dyt(1),dyt(2),dyt(3),dyt(4),c,dc
96    format(a2,8f10.4)
97    format(11x,'r',9x,'i',5x,'dr/dl',5x,'di/dl',5x,'dD/dl',
     &  5x,'dt/dl',9x,'c',5x,'dc/dr')

      do  i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
      enddo

      return
      end

      subroutine rk4r(y,iq,h,ktype,jdown,yout)

      ! integrates l,phi,i,t as function of r

      ! 4th order Runge Kutta solution to y'=f(x,y):
      ! y(n+1)=y(n)+(k1+2k2+2k3+k4)/6

      ! input (note y(1) is not r but l!):
      ! y(1)=l, y(2)=i, y(3)=phi, y(4)=t
      ! iq=guess for node # directly below y(1)
      ! h: step size (km)
      ! ktype=1 for P, 2 for S
      ! jdown=1 if step is decreasing r, 0 if upward ray

      ! output
      ! yout(1)-yout(4): next ray point
      ! iq is adjusted if needed

      dimension dydx(4),y(4),yout(4)
      dimension dym(4),dyt(4),yt(4)

      hh=h*0.5
      h6=h/6.0
      n=4

      kupr=jdown      ! choose proper layer if discontinuity
      call derivr(y,dydx,iq,1-kupr,ktype,c,dc,d2c)
      do  i=1,n
        yt(i)=y(i)+hh*dydx(i)  
      end do

      call derivr(yt,dyt,iq,kupr,ktype,c,dc,d2c)

      do i=1,n
        yt(i)=y(i)+hh*dyt(i)    
      enddo

      call derivr(yt,dym,iq,kupr,ktype,c,dc,d2c)

      do i=1,n
        yt(i)=y(i)+h*dym(i)   
        dym(i)=dyt(i)+dym(i)
      end do

      call derivr(yt,dyt,iq,kupr,ktype,c,dc,d2c)

      do  i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
      enddo

      return
      end

      subroutine mkrtable(rseg,ktseg,kdwn,nseg,a1,a2,ytbl,ntbl)

      ! traces rays and produces ray table ytbl

      ! input:
      ! ray geometry in rseg,ktseg,kdwn,nseg (see subroutine rdray)
      ! a1,a2 start/final ray path angle at source (in degrees,
      !    0<= a1 < a2 < 90 if upward, else 90 < a1 < a2 < 180)
      ! a1 and a2 override ray segment angles. Set (0,90) or (90,180)
      !    to let tracer decide what is the full sweep of angles
      ! model in common modl/ (see subroutine model)

      ! output:
      ! ytbl(1,j)=i, ytbl(2,j)=Delta, ytbl(3,j)=T,  ytbl(4,j)=rturn
      ! for j=1,ntbl, with i and Delta in degrees, T in seconds
      ! a1,a2 - adjusted if out of bounds dictated by the ray

      ! shortcoming:
      ! p and P have to be treated with different ray definition
      ! files

      dimension rseg(20),ktseg(20),kdwn(20),ytbl(4,4000)
      dimension ystart(4),y(4,9000),jsg(9000)
      ! background model
      common /modl/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &      Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      data r2d/57.295779513082/,pi/3.14159265359/
      data halfpi/1.570796326795/, eps/1.0/

      jdebug=0

      h=15.    ! integration step size (keep equal to h in raydata)
      if(a1.gt.a2) stop 'mkrtable: a1>a2'

      iqsrc=n
      ktsrc=ktseg(1)
      rsrc=rseg(1)
      call velo1(rseg(1),iqsrc,1-kdwn(1),ktsrc,csrc,dc,d2c) ! find csrc

      ! check max slowness to avoid that ray turns away from surface
      vmax=vps(iqsrc,ktsrc)
      do i=iqsrc+1,n
        vmax=max(vmax,vps(i,ktsrc))
      enddo
      plim=0.99/vmax
      alim=asin(csrc*plim)
      if(jdebug.eq.1) write(13,*) 'vmax=',vmax,'  alim=',alim*r2d

      if(jdebug.gt.0) write(13,*) 'rsrc, csrc=',rsrc,csrc
      ra1=a1/r2d
      ra2=a2/r2d
      if(kdwn(1).eq.0) then
        amin=max(0.00001,ra1)   ! upward: 0<angle<halfpi
        amax=min(alim,ra2)
      else
        amin=max(pi-alim,ra1)    ! down: halfpi<angle<pi
        amax=min(pi-0.00001,ra2)
      endif  
      if(jdebug.eq.1) write(13,*) 'amin,amax=',amin*r2d,amax*r2d
      ! find slowness range for defined turning points.
      pturn=0.
      ptrns=(rsrc-eps)/csrc
      p=ptrns   ! (defines p for first debug output only)
      do i=1,nseg
        if(kdwn(i).eq.2) then   ! if this is a turning ray
          rturn=rseg(i)+eps
          kturn=ktseg(i)
          call findiq(iqmin,1,rturn)
          cturn=vps(iqmin,kturn)
          p=rturn/cturn
          pturn=max(p,pturn)
          if(jdebug.eq.1) write(13,*) 'cturn,p,pturn=',cturn,p,pturn
          sintrn=pturn*csrc/rsrc
          if(abs(sintrn).gt.1.0) stop 'Unexpected isin>1'
          aturn=asin(sintrn)  ! angle for turning wave
          if(kdwn(1).eq.0) then
            amin=max(amin+0.00001,aturn)        ! up: raise minimum angle
          else
            amax=min(amax-0.00001,pi-aturn)     ! down: lower max angle
          endif  
          if(jdebug.eq.1) write(13,*)'aturn,amin,amax=',aturn,amin,amax
        else if(kdwn(i).eq.3.or.kdwn(i).eq.4) then 
          rtrns=rseg(i)-eps     ! ray needs to reach at least here
          ktrns=ktseg(i-1)      ! to reflect at or transmit this layer 
          call findiq(iqmin,1,rtrns)
          ctrns=vps(iqmin,ktrns)
          p=rtrns/ctrns
          ptrns=min(p,ptrns)
          ktrns=ktseg(i)
          ctrns=vps(iqmin,ktrns)
          p=rtrns/ctrns
          ptrns=min(p,ptrns)
          if(jdebug.eq.1) write(13,*) 'crefl/tr,p,ptrns=',ctrns,p,ptrns
          sintrs=ptrns*csrc/rsrc
          atrns=asin(sintrs)  ! angle for first refl/transmitted wave
          if(kdwn(1).eq.0) then
            amax=min(amax,atrns+0.00001)        ! up: lower max angle
          else
            amin=max(amin,pi-atrns-0.00001)     ! down : raise min angle
          endif  
          if(jdebug.eq.1) write(13,*)'atrns,amin,amax=',atrns,amin,amax
        endif
        if(jdebug.eq.1) write(13,*) i,kdwn(i),p,amin,amax
      enddo  
      ra1=amin
      ra2=amax
      iqmax=n
      if(kdwn(1).eq.1) iqmax=iqsrc

      ddtarget=0.2                      ! target delta spacing (deg)
      drtarget=ddtarget/r2d

      ystart(1)=rseg(1)         ! source radius
      ystart(3)=0.              ! source delta (phi in Dahlen 2000)
      ystart(4)=0.              ! source travel time
      ntbl=0


      ! find dpdd at start of table
      dra=0.0025*min(1.0,ra2-ra1)
10    ystart(2)=ra1
      if(jdebug.gt.0) write(13,*) 'start mkrtable a=',ra1*r2d,' to ',
     & ra2*r2d
      sini=ystart(2)
      if(jdebug.gt.0) write(13,*) 'ntbl angle       delta    time dra'
      call tracer(ystart,rseg,ktseg,kdwn,nseg,y,jsg,rmin,tstar,h,nray)
      if(nray.le.0) then
        ra1=ra1+dra
        if(jdebug.gt.0) write(13,*) 'tracer returns nray=',nray,' i=',
     &     ystart(2)*r2d
        goto 10
      endif  
      ntbl=ntbl+1
      ytbl(1,ntbl)=ystart(2)*r2d        ! store i
      ytbl(2,ntbl)=y(3,nray)*r2d        ! delta
      ytbl(3,ntbl)=y(4,nray)            ! time
      ytbl(4,ntbl)=rmin                 ! turning point radius
      if(jdebug.gt.0)
     &  write(13,*) ntbl,ystart(2)*r2d,y(3,nray)*r2d,y(4,nray),dra*r2d
      alast=ra1
      afirst=ra1
      dlast=y(3,nray)
      dadd=dra/drtarget         ! default in case first ddelta is zero

      do while(ystart(2).lt.ra2-0.0001)
        ystart(2)=min(ra2,ystart(2)+dra)
        call tracer(ystart,rseg,ktseg,kdwn,nseg,y,jsg,rmin,tstar,h,nray)
        if(nray.le.0) then      ! still prblematic at low angles, restart
          ra1=ystart(2)+dra
          ntbl=0
          if(jdebug.gt.0) write(13,*) 'unexpected nray=',nray,' i=',
     &     ystart(2)*r2d
          if(ra1.lt.ra2) goto 10
          if(jdebug.gt.0) write(13,*) 'Unexpected stop of mkrtable'
          return
        endif
        ntbl=ntbl+1
        if(ntbl.gt.4000) stop 'ntbl>4000'
        ytbl(1,ntbl)=ystart(2)*r2d        ! store i
        ytbl(2,ntbl)=y(3,nray)*r2d        ! delta
        ytbl(3,ntbl)=y(4,nray)            ! time
        ytbl(4,ntbl)=rmin                 ! turning point radius
        if(jdebug.gt.0)
     &    write(13,*) ntbl,ystart(2)*r2d,y(3,nray)*r2d,y(4,nray),dra*r2d
        ! step between 0.057 and 0.57 degree (0.001-0.01 rad)
        if(ytbl(4,ntbl).gt.0.) then
          dra=min(0.01,max(0.001,abs(drtarget*dadd)))
        else
          dra=sign(abs(dra),dadd)
        endif
        alast=ystart(2)
        dlast=y(3,nray)
      enddo
      a1=min(afirst,alast)*r2d
      a2=max(afirst,alast)*r2d
      if(jdebug.gt.0) write(13,*) 'Final table between i=',ra1*r2d,
     &   ' and ',alast*r2d

      return

      end



      subroutine derivpq(kupr,ktype,iq,y,dydx,p)

      ! computes derivatives (Dahlen et al., GJI 141:157-174, 2000)

      ! input:
      ! kupr - 1 to select c at upper layer of discontinuity, else 0
      ! ktype - 1 for P, 2 for S
      ! iq - guess for node number directly below rr
      ! y - 10 dimensional vector of P,Q,r and i
      ! p - slowness (s/rad)

      ! vector y is (P1,P2,Q1,Q2,~P1,~P2,~Q1,~Q2,r,i) in
      ! Dahlen's notation

      ! output:
      ! dydx - derivatives of y w.r.t. phi (Dahlen eqs A30, A31)

c     double precision y,dydx
      dimension y(10),dydx(10)

      jdebug=0

      rr=y(9)
      call velo1(rr,iq,kupr,ktype,c,dc,d2c)
      ri=1.0/rr
      r2=rr*rr
      tai=tan(y(10))
      cti=1.0/tai
      cti2=cti*cti
      if(jdebug.gt.0) write(13,*) 'rr,i,d,dc,d2c=',y(9),y(10),c,dc,d2c
c     if(jdebug.gt.0) write(13,*) 'tai,cti,ct2=',tai,cti,cti2
c     if(jdebug.gt.0) write(13,90) '    y=',(y(i),i=1,8)
      c3=c*c*c
      dydx(1)=-y(3)*p*(d2c+dc*cti2*ri)/c        ! (A30)
      dydx(2)=-y(4)*rr*dc/(c3*p)                 ! (A31)
      dydx(3)=y(1)*r2/p                         ! (A30)
      dydx(4)=y(2)*r2/p                         ! (A31)
      dydx(5)=-y(7)*p*(d2c+dc*cti2*ri)/c        ! (A30)
      dydx(6)=-y(8)*rr*dc/(c3*p)                 ! (A31)
      dydx(7)=y(5)*r2/p                         ! (A30)
      dydx(8)=y(6)*r2/p                         ! (A31)
      dydx(9)=rr*cti                             ! (A6)
      dydx(10)=rr*dc/c-1.0                       ! (A6)

      if(jdebug.gt.0) write(13,90) 'deriv:',(dydx(i),i=1,8)
90    format(a6,8e12.3)

      return
      end

      subroutine rk4pq(y,iq,h,ktype,jdown,yout,p)

      ! integrates y w.r.t. phi

      ! 4th order Runge Kutta solution to y'=f(x,y):
      ! y(n+1)=y(n)+(k1+2k2+2k3+k4)/6

      ! input:
      ! y(1-10) (see subroutine hessian ypq)
      ! iq=guess for node # directly below y(1)
      ! h: step size (rad)
      ! ktype=1 for P, 2 for S
      ! jdown=1 if step is decreasing r, 0 if upward ray

      ! output
      ! yout(1)-yout(10): next ray point
      ! iq is adjusted if needed

c     double precision dydx,y,yout,dym,dyt,yt
      dimension dydx(10),y(10),yout(10)
      dimension dym(10),dyt(10),yt(10)
      data halfpi/1.570796326795/

      jdebug=0

      hh=h*0.5
      h6=h/6.0
      n=10

      kupr=jdown      ! choose proper layer if y(1) at discontinuity
      call derivpq(1-kupr,ktype,iq,y,dydx,p)
      do  i=1,n
        yt(i)=y(i)+hh*dydx(i)  
      end do

      call derivpq(kupr,ktype,iq,yt,dyt,p)
      do i=1,n
        yt(i)=y(i)+hh*dyt(i)    
      enddo
      call derivpq(kupr,ktype,iq,yt,dym,p)

      do i=1,n
        yt(i)=y(i)+h*dym(i)   
        dym(i)=dyt(i)+dym(i)
      end do
      call derivpq(kupr,ktype,iq,yt,dyt,p)

      do  i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
      enddo

      if(jdebug.eq.1) write(13,90) yout
90    format('yout=',10e12.3)

      return
      end

      subroutine getcr2(flat,flon,ktype,vcr2,rcr2,ncr2,ecr2,rmoho,vsed)

      ! returns crust1.0 model in a format usable for subroutine inttau.
      ! input: latitude and longitude flat, flon in degrees, ktype (1=P,2=S)
      ! output: arrays vcr2 (velocity) and rcr2 (radius) with ncr2 elements,
      !         and elevation ecr2 in km, rmoho the radius of the 
      !         model bottom, only in case this exceeds the crustal
      !         thickness of the reference model is this the CR2 moho...
      !         vsed is the velocity in the hard sedimentary layer.
      !         rcrust is the radius of upper crust surface
      ! ASSUMES r(n)=6371

      parameter(ityp=36,nla=180,nlo=360)                        ! for rcrust1 arrays

      character*2 types
      dimension vcr2(16),rcr2(16)
      integer :: sed_layer
      real :: sed_thickness

      common/crust1/amapvp(9,nlo,nla),amaprho(9,nlo,nla),
     &              amapvs(9,nlo,nla),amapthi(8,nlo,nla),
     &              amapele(nlo,nla), am_crust(43621,3),
     &              am_sed(43621,3), am_moho(43621,3)                    ! model CRUST1.0
      common/glbgrd/types(nlo,nla)                      ! model CRUST1.0
      ! background model
      common /modl/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &      Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      jdebug=1
      mdebug=0
      
      cola=90.-flat
      cola=min(cola,179.9)              ! fixes S Pole bug
      if(flon.gt.180.)flon=flon-360.
      dx=360.0/nlo
      ilat=min(180,int(cola/dx)+1)
      ilon=min(360,int((flon+180.)/dx)+1)

      ecr2=amapele(ilon,ilat)/1000.             ! converts m to km
c      XXX MT: Hacky. The problem is that crust2 (only 2 sediment layers) ALWAYS
c      has some kind of hard sediment velocity, BUT in crust1 (3 sediment layers)
c      we experience that in some places there is NO hard and/or middle sediment 
c      layer at all, therefore the hacky following lines.
c      ====>>>> CHECK IF THICKNESS IS REQUIRED FOR SEDIMENT LAYER -> EDIT: 
c      in raydata.f we actually only consider the velocity of the sedimentlayer only
c      if the layer is not deeper than 1km.
      vsed=amapvp(5,ilon,ilat)
      sed_layer = 5
      if(amapvp(5,ilon,ilat).eq.0)then
          vsed=amapvp(4,ilon,ilat)
          ! XXX waiting for minakovs answer
          sed_layer = 4 ! I'm intorducing this parameter to later check which layer is getting
          ! what thickness 
          if (amapvp(4,ilon,ilat).eq.0)then
              vsed=amapvp(3,ilon,ilat)
              sed_layer = 3
          endif
      endif
      
      if(ktype.eq.2) vsed=amapvs(5,ilon,ilat) ! XXX this will not work for ktype=2 waves!!!

      a=6371.0+ecr2                             ! r of last solid layer
      ncr2=0
      rr=a
      if(jdebug.eq.1) print *, 'getcr2',flat,flon,ilat,ilon,ktype
      
      ! for the thicknesses we have to deferentiate between the crust1.0 grid and the finer 
      ! minakov gird
      ! corner points of Minakovs map; ulc=upper left corner, etc...
      ulc_lon=40
      ulc_lat=-10
      lrc_lon=80
      lrc_lat=-40
      ! if latitude or longitude is in Minakovs compilation then find corresponding point and get thickness 
      ! of layer
      if (flon .ge. ulc_lon .and. flon .le. lrc_lon .and. 
     &    flat .ge. lrc_lat .and. flat .le. ulc_lat) then
          ! we now the step size from Minakovs files
          if(mdebug.eq.1) print *, '\\\\\\\'
          if(mdebug.eq.1) print *, 'flon', flon, 'flat', flat
          ! -------- LONGITUDE
          flon_remain = flon-int(flon)
          if(mdebug.eq.1) print *, 'flon_remain', flon_remain
          if(flon_remain.eq.0)then
              flon_checked = flon
              if(mdebug.eq.1) print *, 'flon_checked -> 0 =', flon_checked
          else
              flon_checked = flon - modulo(flon,(1/6.))
              if(mdebug.eq.1) print *, 'flon_checked -> !=0 =', flon_checked
          endif
          ilon=min(real(360.),real(flon_checked+180.)+1)
          
          ! -------- LATITUDE
          flat_remain = flat-int(flat)
          if(mdebug.eq.1) print *, 'flat_remain', flat_remain
          if(flat_remain.eq.0)then
              flat_checked = flat
              if(mdebug.eq.1) print *, 'flat_checked -> 0 =', flat_checked
          else
              flat_checked = flat - modulo(flat,(1/6.))
              if(mdebug.eq.1) print *, 'flat_checked -> !=0 =', flat_checked
          endif
          cola=90.-flat_checked
          cola=min(real(cola),real(179.9))   ! fixes S Pole bug
          ilat=min(real(180.),real(cola)+1)
          
          if(jdebug.eq.1) print *, 'ilon', ilon, 'ilat', ilat
          
          ! -------- THICKNESS
          ! print *, am_crust(100,1)
          ! print *, am_crust(100,2)
          ! print *, am_crust(100,3)
          max_rows = 43621
          ! finding the location in the file for the corresponding coordinates
          loc = minloc(abs(am_crust(:,1) - flon), 1)
          loc2 = minloc(abs(am_crust(loc:max_rows, 2) - flat), 1)
          loc_final = loc + loc2 - 1
          
          if(jdebug.eq.1) print *, am_crust(loc_final, 3)
          if(jdebug.eq.1) print *, am_sed(loc_final, 3)
          if(jdebug.eq.1) print *, am_moho(loc_final,3)
          
          ncr2=0

          do i=2,8
              if (i .eq. sed_layer .and. am_sed(loc_final, 3).gt.0) then
                  ncr2=ncr2+1
                  vcr2(ncr2)=amapvp(i,ilon,ilat)
                  sed_thickness = am_sed(loc_final, 3)
                  ! print *, 'thi(',i,')=', sed_tickness, crust_thickness
                  if(jdebug.eq.1) print *, i,ncr2,vcr2(ncr2),rcr2(ncr2)
              
                  ncr2=ncr2+1
                  vcr2(ncr2)=amapvp(i,ilon,ilat)
                  rcr2(ncr2) = rr - sed_tickness - crust_tickness
                  if(jdebug.eq.1) print *, i,ncr2,vcr2(ncr2),rcr2(ncr2)
              
              else if (i .eq. 8 .and. am_crust(loc_final, 3).gt.0) then
                  ncr2=ncr2+1
                  vcr2(ncr2)=amapvp(i,ilon,ilat)
                  crust_thickness = am_crust(loc_final,3)
                  !print *, 'thi(',i,')=', sed_tickness, crust_thickness
                  if(jdebug.eq.1) print *, i,ncr2,vcr2(ncr2),rcr2(ncr2)
              
                  ncr2=ncr2+1
                  vcr2(ncr2)=amapvp(i,ilon,ilat)
                  rcr2(ncr2) = rr - sed_tickness - crust_tickness
                  if(jdebug.eq.1) print *, i,ncr2,vcr2(ncr2),rcr2(ncr2)
              !else 
              !    ncr2=ncr2+1
              !    sed_thickness = 0
              !    crust_thickness = 0
              !    vcr2(ncr2)=0
              end if
          enddo
          r_moho = rr-am_moho(loc_final,3)
          if(rr.gt.r_moho) then
            ncr2=ncr2+1
            rcr2(ncr2)=rr
            vcr2(ncr2)=vp(nmoh-1)           ! Reference model Vp
            if(ktype.eq.2) vcr2(ncr2)=vs(nmoh-1)
            ncr2=ncr2+1
            vcr2(ncr2)=vcr2(ncr2-1)
            rcr2(ncr2)=r(nmoh)
            if(jdebug.eq.1) then
              write(13,*) 'CR2 Moho is at',rr,' vs reference',r(nmoh)
              write(13,*) 'Added:',rr-r(nmoh),' km with v=',
     &           vcr2(ncr2)
            endif
            rr=r(nmoh)
          endif  
          rmoho=rr
          
      else
          ncr2=0
          ! here we are at the original code using either crust2.0 or 1.0
          do i=2,8
            !print *, 'ampthi(',i,')=', amapthi(i,ilon,ilat)
            if(amapthi(i,ilon,ilat).gt.0.) then
c               what basically happens is that we define the upper boundary of the layer 
c               and in the second part of the loop the lower layer: 14 entries without the
c               water layer
                ncr2=ncr2+1
                vcr2(ncr2)=amapvp(i,ilon,ilat)
                if(ktype.eq.2) vcr2(ncr2)=amapvs(i,ilon,ilat)
                rcr2(ncr2)=rr
                !if(jdebug.eq.1) print *, i,ncr2,vcr2(ncr2),rcr2(ncr2)
                rr=rr-amapthi(i,ilon,ilat)
c                  if(i.eq.5) rcrust=rr ! radius of upper crust surface
                ! XXX next lines double - NO, make the figure in the end look right - WHY?
                ncr2=ncr2+1
                vcr2(ncr2)=amapvp(i,ilon,ilat)
                if(ktype.eq.2) vcr2(ncr2)=amapvs(i,ilon,ilat)
                rcr2(ncr2)=rr
                !if(jdebug.eq.1) print *, i,ncr2,vcr2(ncr2),rcr2(ncr2)
c             else
c                 print *, 'i=', i, 'ncr2=', ncr2, 'vcr2(ncr2)=', 
c     &           vcr2(ncr2), 'rcr2(ncr2)=', rcr2(ncr2)
             endif
          enddo
          ! add mantle layer if reference Moho is deeper
          print *, 'nmoh', nmoh, 'r(nmoh)=', r(nmoh)
          if(rr.gt.r(nmoh)) then
            ncr2=ncr2+1
            rcr2(ncr2)=rr
            vcr2(ncr2)=vp(nmoh-1)           ! Reference model Vp
            if(ktype.eq.2) vcr2(ncr2)=vs(nmoh-1)
            ncr2=ncr2+1
            vcr2(ncr2)=vcr2(ncr2-1)
            rcr2(ncr2)=r(nmoh)
            if(jdebug.eq.1) then
              write(13,*) 'CR2 Moho is at',rr,' vs reference',r(nmoh)
              write(13,*) 'Added:',rr-r(nmoh),' km with v=',
     &           vcr2(ncr2)
            endif
            rr=r(nmoh)
          endif  
          rmoho=rr
      endif
      !if(jdebug.eq.1) print *, 'ncr2', ncr2
      !if(jdebug.eq.1) print *, 'rcr2', rcr2(ncr2), 'vcr2', vcr2(ncr2)
      if(jdebug.eq.1) print *, 'Moho at r=',rmoho,' z=',6371-rmoho
      !print *, '======'
      !print *, ilat, ilon
      !print *, 'rcr2=', rcr2
      !print *, 'vcr2=', 
      
      
      !print *, '======'
      !print *, 'ncr2=', ncr2
      !print *, 'rcr2=', rcr2(ncr2)
      !print *, 'vcr2=', vcr2(ncr2)
      !print *, 'rmoho=', rmoho
      
      return
      end
      
      ! if latitude or longitude is in Minakovs compilation then find corresponding point and get thickness 
      ! of layer
      !  if (flon .ge. ulc_lon .and. flon .le. lrc_lon .and. 
      ! &    flat .ge. lrc_lat .and. flat .le. ulc_lat) then
      !      ! we now the step size from Minakovs files
      !      print *, '\\\\\\\'
      !      print *, 'flon', flon, 'flat', flat
      !      ! -------- LONGITUDE
      !      flon_remain = flon-int(flon)
      !      print *, 'flon_remain', flon_remain
      !      if(flon_remain.eq.0)then
      !          flon_checked = flon
      !          print *, 'flon_checked -> 0 =', flon_checked
      !      else
      !          flon_checked = flon - modulo(flon,(1/6.))
      !          print *, 'flon_checked -> !=0 =', flon_checked
      !      endif
      !      ilon=min(real(360.),real(flon_checked+180.)+1)
      !      
      !      ! -------- LATITUDE
      !      flat_remain = flat-int(flat)
      !      print *, 'flat_remain', flat_remain
      !      if(flat_remain.eq.0)then
      !          flat_checked = flat
      !          print *, 'flat_checked -> 0 =', flat_checked
      !      else
      !          flat_checked = flat - modulo(flat,(1/6.))
      !          print *, 'flat_checked -> !=0 =', flat_checked
      !      endif
      !      cola=90.-flat_checked
      !      cola=min(real(cola),real(179.9))   ! fixes S Pole bug
      !      ilat=min(real(180.),real(cola)+1)
      !      
      !      print *, 'ilon', ilon, 'ilat', ilat
      !  else
      !      ! here we are at the original code using either crust2.0 or 1.0
      !      cola=90.-flat
      !      cola=min(cola,179.9)              ! fixes S Pole bug
      !      if(flon.gt.180.)flon=flon-360.
      !      dx=360.0/nlo
      !      ! ilat and ilon are all the coordinates of crust1.0
      !      ilat=min(180,int(cola/dx)+1)
      !      ilon=min(360,int((flon+180.)/dx)+1)
      !  endif

      subroutine getref(ktype,vref,rref,nref,rmoho)

      ! returns reference model in a format usable for subroutine inttau.
      ! input: rmoho (from getcr2), the moho radius in km, ktype (1=P,2=S).
      ! output: arrays vref (velocity) and rref (radius) with nref elements.

      dimension vref(50),rref(50)

      ! background model
      common /modl/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &      Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      nref=1
      vref(1)=vp(n)
      if(ktype.eq.2) vref(1)=vs(n)
      rref(1)=r(n)
      i=n

      ! map reference model in reversed order until Moho is reached
      do while(r(i).gt.rmoho)
        nref=nref+1
        i=i-1
        if(nref.gt.50) stop 'getref: nref>50, increase dimensions'
        if(i.lt.1) stop 'getref: i<1; crustal thickness error?'
        vref(nref)=vp(i)
        if(ktype.eq.2) vref(nref)=vs(i)
        rref(nref)=r(i)
      enddo

      if(abs(rref(nref)-rmoho).lt.0.01) return

      ! interpolate last model node at rmoho
      dr=r(i+1)-r(i)
      if(dr.le.0.) stop 'getref: dr<0'
      vref(nref)=vp(i)+(rmoho-r(i))*(vp(i+1)-vp(i))/dr
      if(ktype.eq.2) vref(nref)=vs(i)+(rmoho-r(i))*(vs(i+1)-vs(i))/dr
      rref(nref)=rmoho

      return
      end

      subroutine inttau(p,r,rsrc,v,n,tau,tausrc)

      ! integrates for the intercept time tau throughout model
      ! input: slowness p (s/rad), radii r (km), velocities (km/s) in n nodes
      ! output: intercept time tau (s) over the whole model and tausrc from
      ! tausrc is only correct if r is in descending order!
      ! source to surface.

      dimension r(50),v(50),f(50)

      jdebug=0
      if(jdebug.eq.1) write(13,98) p
98    format('Enter inttau slowness ',f8.1,/,'  i       r       v',
     & '    sini      dr    dtau     tau')
99    format(i3,f8.1,5f8.3)

      ! integrate tau = integral cos(i)/v dr using trapezium rule
      tau=0
      tausrc=-1
      dtau=0

      ! compute the integrand
      do i=1,n
        sini=p*v(i)/r(i)
        f(i)=sqrt(1.0-sini*sini)/v(i)
      enddo  

      ! trapezoidal quadrature
      do i=2,n
        dr=r(i)-r(i-1)
        dr=abs(dr)                              ! r is descending
        dtau=0.5*dr*(f(i)+f(i-1))
        tau=tau+dtau
        if(jdebug.eq.1) write(13,99) i,r(i),v(i),sini,dr,dtau,tau
        ! the following assumes descending r - are we past the source?
        if(tausrc.eq.-1.and.r(i).le.rsrc) then
          tau0=tau-dtau         ! tau at r(i-1)
          tausrc=tau0+dtau*(rsrc-r(i-1))/(r(i)-r(i-1))
          if(jdebug.eq.1)write(13,99) 0,rsrc,vsrc,sini,abs(dr),0.,tausrc
        endif  
      enddo
      if(rsrc.le.r(n)) tausrc=tau
      if(jdebug.eq.1) write(13,*) 'Exit inttau tau,tausrc=',tau,tausrc

      return
      end

      subroutine rcrust1

c read infos to use to compute scripps crustal corrections
c required by vpcrust1
c 360 key profile: ctype(ityp),fvel(ityp,8),fdep(ityp,8)
c types(ilat,ilon)

c Reads the CRUST1.0 model files CNtype1_key.txt, 'CNtype1.txt and CNelevatio1.txt.
c It first looks for these files in $ directory specified by environment
c variable CRUST10, if it cannot find them it looks in the current directory

c fdep contains thickness of the crustal layers, fdep(ityp,8)
c is the total thickness of the crustal model along that profile

      parameter(ityp=36,nla=180,nlo=360)
c     ---->>>> new stuff
      parameter(np=9)
c     <------

      dimension fvel(ityp,4),fvels(ityp,4),frho(ityp,4),fthi(ityp,4)

c     ---->>>> new stuff
      dimension vp(np,nla,nlo),vs(np,nla,nlo),rho(np,nla,nlo)
      dimension bnd(np,nla,nlo)
c     <<<<----
      
      character dum*1,line*506
      character dum0*5,fname*72
      character*2 ctype(ityp),ctypes(nla,nlo),types,atype(nlo)
      logical ruthere

      common/crust1/amapvp(9,nlo,nla),amaprho(9,nlo,nla),
     &              amapvs(9,nlo,nla),amapthi(8,nlo,nla),
     &              amapele(nlo,nla), am_crust(43621,3),
     &              am_sed(43621,3), am_moho(43621,3)
      common/glbgrd/types(nlo,nla)
      
      ! XXX MINAKOV extra: ATTENION: numbers are hard coded here
      integer :: row, col, max_riws, max_cols

c     ---->>>> new stuff
c read the files for crust1.0 all layer information
      CHARACTER(len=255) :: cwd
      CALL getcwd(cwd)
      WRITE(*,*) 'Working in following directory:', TRIM(cwd)
      
c      inquire(file='crust1.vp',exist=ruthere)
c      if(ruthere) then
c          print *, 'yay'
c      else
c          print *, 'cry me a river'
c      endif

c checking if the files are where they supposed to be... 
c and reading them 

      fname='CNelevatio1.txt'
      write(6,'(2a)') 'Looking for ',fname
      inquire(file=fname,exist=ruthere)
      if(ruthere) then
        open(4,file=fname,status='old')
        fname='CNtype1.txt'
        open(3,file=fname,status='old')
        fname='CNtype1_key.txt'
        open(2,file=fname,status='old')
      else  
        inquire(file='CNelevatio1.txt', exist=ruthere)
        if(.not.ruthere) then
          print *,'Cannot find Crust1.0 files. Either move them to '
          print *,'the current directory or set environment variable'
          print *,'CRUST10 to the directory where these files are.'
          stop 'Error locating Crust1.0 files'
        endif  
        print *,'Opening CRUST1.0 files in current directory'
        open(2,file='CNtype1_key.txt',status='old')
        open(3,file='CNtype1.txt',status='old')
        open(4,file='CNelevatio1.txt',status='old')
      endif  

c reading the crust information files for the crust types

      open(51,file='crust1.vp', status='old')
      open(52,file='crust1.vs', status='old')
      open(53,file='crust1.rho', status='old')
      open(54,file='crust1.bnds', status='old')
      
c here we are reading the Minkov maps
      open(71, file='AMcrustIndian.txt', status='old')
      open(72, file='AMsedIndian.txt', status='old')
      open(73, file='AMmohoIndian.txt', status='old')
      

      print*,' .... read all maps information from CRUST1.0... ' 

      do j=1,nla
         do i=1,nlo
            read(51,*)(vp(k,j,i),k=1,np)
            read(52,*)(vs(k,j,i),k=1,np)
            read(53,*)(rho(k,j,i),k=1,np)
            read(54,*)(bnd(k,j,i),k=1,np)
         enddo
      enddo
      close(51)
      close(52)
      close(53)
      close(54)
      !write(*,*)'shape(bnd)=', shape(bnd)
      !write(*,*)'shape(amapthi)=', shape(amapthi)
      
      ! Minakov extra map part
      max_rows = 43621
      max_cols = 3
      do row = 1, max_rows
          read(71, *) (am_crust(row, col), col=1, max_cols)
          read(72, *) (am_sed(row, col), col=1, max_cols)
          read(73, *) (am_moho(row, col), col=1, max_cols)
      end do
      close(71)
      close(72)
      close(73)
      
      
      read(2,890)dum
      do 101 i=1,ityp
         read(2,899)ctype(i)
         read(2,899)line
         read(line,*)(fvel(i,l),l=1,4)
         read(2,899)line
         read(line,*)(fvels(i,l),l=1,4)
         read(2,899)line
         read(line,*)(frho(i,l),l=1,4)
         read(2,899)line
         read(line,*)(fthi(i,l),l=1,4)
c flip layers
         aux=fvel(i,1)
         fvel(i,1)=fvel(i,2)
         fvel(i,2)=aux
         aux=fvels(i,1)
         fvels(i,1)=fvels(i,2)
         fvels(i,2)=aux
         aux=frho(i,1)
         frho(i,1)=frho(i,2)
         frho(i,2)=aux
         aux=fthi(i,1)
         fthi(i,1)=fthi(i,2)
         fthi(i,2)=aux
 101  continue

c reading the CNtype file but this is just pro forma, beceause we
c don't really need that anymore
c      open(55,file='CNtype1.txt')
c      do j=1,nla
c         read(55,101)(ctypes(j,i),i=1,360)
c      enddo    
c      close(55)
c      write(*,*)'ctypes=', ctypes
c      print*,' .... we came a long way down ... '
c     <<<<----

      write(*,*)'getting crust1.0 information done'

c ---- read in key per cell of global grid (1x1 degree)
c ---- file CNtype1.txt
c      write(*,*) 'ilat=', ilat
c      write(*,*) 'nlo=', nlo
      
c      write(*,*)'vp=', vp(:,90, 90)
      
      read(3,899)line
      read(line,*)flons
      read(4,899)line
      
      do j=1,nla
        read(3,901)ilat,atype
        read(4, *) ilat,(amapele(i,j),i=1,nlo)
        do i=1,nlo
            do l=1,ityp
              if(atype(i).eq.ctype(l))then
                  types(i,j)=ctype(l)
               endif
                  do k=1,9
                      amapvp(k,i,j)=vp(k,j,i)
                      amapvs(k,i,j)=vs(k,j,i)
                      amaprho(k,i,j)=rho(k,j,i)
                  enddo
c      ====>>>> here we have to compute the thickness from depth
                  do k=1,8
                      amapthi(k,i,j)=bnd(k, j, i) - bnd(k+1,j,i)
c                      write(*,*) 'i', i, 'j', j, 'k', k
c                      write(*,*) 'amapthi(k,i,j)', amapthi(k,i,j)
                  enddo
                  goto 10
               
          enddo
  10   enddo
      enddo
c      write(*,*)'vp=', vp(:,90, 90)
c      write(*,*)'amapthi=', amapthi(:,90, 90)
c      write(*,*)'amapvp=', amapvp(5, 180, 180)
c      write(*,*)'thickness=', shape(amapthi(:,:,:))
c      write(*,*)'amapele=', shape(amapele)
c FORMATS
c 101  format(360(1x,a2)) 
c 102  format(a,2i4,2x,a2)
c 103  format(4f7.2)
c 104  format(a,3f7.2)
 
 890  format(////a)
 899  format(a)
 901  format(i4, 1x, 360(1x, a2, 1x)) ! format for the type file
 
c 902  format(i5, 1x, 360(f5.0, 1x)) ! format for the elevation file

      write(*,*)'done read crust1.0'

      return
      end

      subroutine getlegs(y,rayvel,nray,legend,nlegs)


      ! extracts 'legs' of the ray: the segments in between
      ! the end- or reflection points.
      ! transmissions should in principle also be treated
      ! as 'legs', but the error is hopefully small. The
      ! commented lines are ok if you wish to include them.

      ! input:
      ! y(4,i), i=1,...,nray: ray vector
      !         (y1=radius, y2=delta)
      ! rayvel(i), seismic velocity at node i, <0 if Vs

      ! output:
      ! legend(i), i=1,nlegs: end node of each leg
      ! nlegs: number of legs

      dimension y(4,9000),legend(190),rayvel(9000)

      data eps/0.01/    ! must be compatible with eps in rstp
      data halfpi/1.570796327/

      do i=1,190
        legend(i)=0
      enddo  
      nlegs=0
      do i=2,nray-1

c commented lines in case transmissions need to be handled explicitly
c       if(abs(y(1,i)-y(1,i-1)).lt.eps.and.     ! refl/transm
c    &     abs(y(2,i)-y(2,i-1)).lt.0.001.and.
c    &     abs(rayvel(i)-rayvel(i-1)).gt.0.001) then

        ! recognize reflection point
        if(abs(y(1,i)-y(1,i-1)).lt.eps.and.            ! r remains same
     &     (y(2,i)-halfpi)*(y(2,i-1)-halfpi).lt.-eps) then  ! angle flip
          nlegs=nlegs+1
          legend(nlegs)=i-1
          if(nlegs.ge.190) stop 'getlegs: nlegs>190'
        endif

      enddo

      nlegs=nlegs+1
      legend(nlegs)=nray

      return
      end

      subroutine gettstar(y,nray,jsg,ktseg,tstar,qray)

      dimension y(4,9000),jsg(9000),ktseg(20)
      common /modl/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &      Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      ! computes tstar = \int dt/Q = T/Qray from y,jsg and ktseg 
      ! where y is computed by subroutine tracer

      jdebug=0

      iq=n
      tstar=0.
      if(jdebug.gt.0) write(13,*) 'i    rr       qq     dt      tstar'
      kk=1
      do i=2,nray
        dt=y(4,i)-y(4,i-1)       ! travel time in step
        rr=0.5*(y(1,i)+y(1,i-1)) ! average r of step
        call findiq(iq,kk,rr)     ! find nearest node below rr
        t=rr-r(iq)               ! distance to nearest node
        ktype=ktseg(abs(jsg(i-1)))     ! P or S
        qq=Qs(iq)+t*(qqs(1,iq)+t*(qqs(2,iq)+t*qqs(3,iq)))
        if(ktype.eq.1) qq=0.75*(vp(iq)/max(0.01,vs(iq)))**2*qq
        if(vs(iq)*qq.ne.0.) tstar=tstar+dt/qq
        if(jdebug.gt.0) write(13,*) i,rr,qq,dt,tstar,t,iq
      enddo

      qray=1.0e20
      if(tstar.gt.0.) qray=y(4,nray)/tstar

      return
      end

      subroutine getelcor(scolat,slon,rcolat,rlon,y,nray,jsg,rseg,ktseg,
     &      kdwn,p,rayvel,rayq,ecorr)

      ! calculate ellipticity correction

      ! input:
      ! scolat,slon - source colatitude, longitude (rad, geocentric)
      ! rcolat,rlon - receiver coordinates (rad)
      ! y(k,i) - ray geometry from subroutine tracer: k=r,rayangle,phi,t
      !          where r in km, rayangle and phi in rad (t not used).
      ! nray - number of ray nodes i=1,nray
      ! jsg(i) - segment belonging to ray node i, negative if ray upward
      ! rseg(i) - start radius of segment i (km)
      ! ktseg(i) - wave type of segment i (1=P,2=S)
      ! kdwn(i) - ray direction (1=down,0=up,2=turning point,3=refl,
      !           4=transmission, 5=last point of ray)

      ! output:
      ! ecorr - time to be added to spherical earth travel time to
      !         obtain time in elliptical earth
      ! p - ray parameter in rad/s
      ! rayvel - ray velocity at node (>0 if Vp, negative if Vs)
      ! rayq - array with Qs at ray node

      dimension y(4,9000),jsg(9000),rayvel(9000),rayq(9000)
      dimension rseg(20),ktseg(20),kdwn(20)
      dimension e2g(3,3)

      ! background model
      common /modl/r(500),vp(500),vs(500),qvp(3,500),qvs(3,500),
     &      Qs(500),qqs(3,500),nic,noc,nsl,nicp1,nocp1,nslp1,n,nmoh

      data r2d/57.295779513082/, third/0.333333333333333/

      jdebug=0

      ! compute transformation matrix from equatorial->geocentric
      call euler1(scolat,slon,rcolat,rlon,e2g)
      
      ! anchor ellipticity shift eldr to 0 at elliptical surface
      epsa=geteps(r(n))*r(n)
      
      ! get c and ellipticity at source radius
      jdown=kdwn(1)
      iq=n
      call velo1(y(1,1),iq,1-jdown,ktseg(1),c,dc,d2c)
      p=y(1,1)*sin(y(2,1))/c     ! slowness
      if(jdebug.gt.1) write(13,*) 'p=',p
      rayvel(1)=c
      if(ktseg(1).eq.2) rayvel(1)=-c
      eps=geteps(y(1,1))       ! ellipticity at r
      z=e2g(3,1)               ! cartesian z coordinate of source
      rr=y(1,1)
      deltar=eps*rr*(third-z*z)
      dphi=y(3,2)-y(3,1)
      sum=(rr/c)**3*(deltar/rr)*dc*dphi    ! postpone factor 0.5/p
      sumd=deltar*abs(cos(y(2,1)))/c

      if(jdebug.gt.1) then
      write(13,98)
98    format(5x,'i',3x,'jsg',4x,'kd',5x,'phi',9x,'r',8x,'dr',3x,'dtell',
     &       3x,'sumd',4x,'lat',6x,'y2',7x,'c   clast')
      write(13,99) 1,jsg(1),kdwn(1),y(3,1)*r2d,y(1,1),deltar,
     &       0.5*sum/p,sumd,90.-acos(z)*r2d
99    format(3i6,f8.2,2f10.2,2f7.3,f8.1,3f8.2)
      endif

      ! integrate with phi as independent variable, unless a
      ! discontinuity is crossed (when we apply boundary conditions)
      ! we use trapezoidal rule (but omit factor 0.5 while summing)
      do i=2,nray-1                       ! step from i-1 to i
        jdownold=jdown
        jdown=1
        if(jsg(i).lt.0) jdown=0
        dphi=y(3,i+1)-y(3,i-1)
        rr=y(1,i)
        clast=c
        iseg=abs(jsg(i))
        ! choose correct velocity when crossing discontinuity
        kupr=jdown
        dr=y(1,i)-y(1,i-1)
        if(abs(dr).lt.0.0001) kupr=1-jdown
        call velo1(rr,iq,kupr,ktseg(iseg),c,dc,d2c)
        rayvel(i)=c
        if(ktseg(iseg).eq.2) rayvel(i)=-c
        eps=geteps(y(1,i))
        xx=cos(y(3,i))                     ! ray at equator
        yy=sin(y(3,i))
        z=e2g(3,1)*xx+e2g(3,2)*yy           ! z geographic
        deltar=eps*rr*(third-z*z)
        sum=sum+(rr/c)**3*(deltar/rr)*dc*dphi 
        ! test for discontinuity
        if(abs(dr).lt.0.0001) then
          if(jdown.eq.jdownold) then     ! transmission
            sumd=sumd-deltar*(cos(y(2,i))/c-cos(y(2,i-1))/clast)
          else                           ! reflection
            sumd=sumd-2.*deltar*cos(y(2,i))/c ! works for top&bottom
          endif  
        endif
          
        if(jdebug.gt.1) then
          write(13,99) i,jsg(i),kdwn(iseg),y(3,i)*r2d,y(1,i),deltar,
     &       0.5*sum/p,sumd,90.-acos(z)*r2d,y(2,i)*r2d,c,clast
        endif
      enddo  

      ! treat last ray node
      jdown=0
      iq=n
      call velo1(y(1,nray),iq,jdown,ktseg(iseg),c,dc,d2c)
      rayvel(nray)=c
      if(ktseg(iseg).eq.2) rayvel(nray)=-c
      eps=geteps(y(1,nray))
      xx=cos(y(3,nray))                     ! ray at equator
      yy=sin(y(3,nray))
      z=e2g(3,1)*xx+e2g(3,2)*yy           ! z geographic
      rr=y(1,nray)
      deltar=eps*rr*(third-z*z)
      dphi=y(3,nray)-y(3,nray-1)
      sum=sum+(rr/c)**3*(deltar/rr)*dc*dphi    ! postpone factor 0.5/p
      sumd=sumd+deltar*abs(cos(y(2,nray)))/c   
      if(jdebug.gt.1) then
        write(13,99) i,jsg(i),kdwn(iseg),y(3,i)*r2d,y(1,i),deltar,
     &     0.5*sum/p,sumd,90.-acos(z)*r2d
      endif

      ecorr=0.5*sum/p+sumd
      if(jdebug.gt.0) write(13,*) 'ecorr=',ecorr

      return
      end

      subroutine euler1(ts,ps,tr,pr,e2g)

! calculates the Euler matrix e2g to go from
! equatorial coordinates (cartesian) along source-receiver path to 
! geocentric (again Cartesian). Use euler_rot if you also need the
! inverse, but e2g here is the same as in that routine.

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

      dimension e2g(3,3)

      call cart(ts,ps,xi,yi,zi)
      call cart(tr,pr,xr,yr,zr)

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

      subroutine cart(theta,phi,x,y,z)
      ! used in euler_rot and euler1
      ! colatatitude theta and longitude phi -> Cartesian x,y,z
      s=sin(theta)
      x=s*cos(phi)
      y=s*sin(phi)
      z=cos(theta)
      return
      end

      function geteps(x)

      ! returns the ellipticity of the Earth

      data ricb/1229.48/, rcmb/3484.3/

      if(x.gt.5971.) then
        geteps=0.003224+(x-5971.0)*3.2e-7
        return
      else if(x.gt.5701.) then
        geteps=0.003200+(x-5701.)*8.8889e-8
        return
      else if(x.gt.rcmb) then
        geteps=0.002666+(x-rcmb)*2.40899e-7
        return
      else if(x.gt.ricb) then
        geteps=0.002454+(x-ricb)*9.4021e-8
        return
      else
        geteps=0.002454
      endif

      return
      end
