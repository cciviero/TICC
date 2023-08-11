c subroutines to calculate delta, azimuth, and bounce points, etc...
c from Guy Masters

c Changes:
c Aug 06 - adjusted givloc to return p1 in (-180,180)

      subroutine stadis(colat,colon,scolat,scolon,del,az,t0,p0,t1,p1)

c Computes the epicentral distance and azimuth from source to receiver.
c Latitudes are converted to geocentric latitudes prior to performing
c the computations (it is assumed that input latitudes are geographic).
c  input:
c    colat = source colatitude, (degrees)
c    colon =   "    colongitude,  "      
c    scolat    = station colatidue, 
c    scolon    =    "    colongitude, 
c  output:
c    del   = epicentral distance (degrees)
c    az    = azimuth from source to receiver, measure from North (degrees)
      data rpd/1.745329252e-2/,dpr/57.29577951/

c  first do eq coords.
      if(colat.eq.0.) colat=1.0e-5
      t0=geocen(colat*rpd)    
      p0=colon*rpd
      c0= cos(t0)
      s0= sin(t0)
c     print *,t0,p0,c0,s0
c  now do station coords.
      t1=geocen(scolat*rpd)           
      c1= cos(t1)
      s1= sin(t1)
      p1=scolon*rpd
c     print *,t1,c1,s1,p1
c  now calculate distance
      dp=p1-p0
      co=c0*c1+s0*s1*cos(dp)
      si=sqrt(1.-co*co)
c     print *,dp,co,si
      del=atan2(si,co)*dpr
c  now calculate azimuth
      caz=(c1-c0*co)/max(1.0e-30,si*s0)
      dp2=-dp
      saz=-s1*sin(dp2)/si
      az= atan2(saz,caz)*dpr
      if(az.lt.0.0) az=360.0 + az  
c     print *,'stadis returns',del,az
      return
      end
c--------------------------------------------------------------------
      function geogrf(x)
c input:
c   x    = geocentric colatitude (radians)
c output:
c   geogrf = geographic colatitude (radians
c (n.b. fac=(1-f)**2)
c
      data fac/0.993305621334896/
      data pi2/1.570796327/

      geogrf=pi2-atan(cos(x)/(fac*max(1.0e-30,sin(x)))) 
      return                       
      end
c---------------------------------------------------------------------
      subroutine givloc(colat,colon,del,az,t1,p1)
c input:
c   colat,colon = source colatitude and colongitude 
c   del         = distance of point in degrees
c   az          = azimuth of point from source
c output:
c   t1 = point latitude  ( + = N, - = S)
c   p1 = point longitude ( + = E, - = W)
      data rpd/1.745329252e-2/,dpr/57.29577951/
      delr=del*rpd
      azr=az*rpd
      t0=geocen(colat*rpd)    
      ctheta= sin(delr)*sin(t0)*cos(azr) + cos(t0)*cos(delr)
      t1= acos(min(1.,max(-1.,ctheta)))
      if (t0.eq.0.0) then
        p1=az
      elseif (t1.eq.0.0) then
        p1=0.0
      else
        sphi= sin(delr)*sin(azr)/sin(t1)
        cphi=(cos(delr) - cos(t0)*ctheta)/(sin(t0)*sin(t1))
        p1=colon + atan2(sphi,cphi)*dpr
      endif
      t1=90.0-geogrf(t1)*dpr      
      if (p1.gt.360.0) p1 = p1 - 360.0 
      if (p1.gt.180.0) p1 = p1 - 360.0 
      if (p1.lt.-360.0) p1 = p1 + 360.0 
      if (p1.lt.-180.0) p1 = p1 + 360.0
      return
      end
c---------------------------------------------------------------------
      function geocen(arg)
c input:
c   arg    = geographic colatitude (radians)
c output:
c   geocen = geocentric colatitude (radians)
c (n.b. fac=(1-f)**2)
      data fac/0.993305621334896/
      data pi2/1.570796327/
      geocen=pi2-atan(fac*cos(arg)/(max(1.0e-30,sin(arg))))
      return
      end
