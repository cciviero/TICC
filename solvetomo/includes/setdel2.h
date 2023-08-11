      real*8   points  
      real*8   centres 
      real*8   data
      integer  vertices
      integer  neighbour
      integer  nnn
      integer  nflist  
      integer  ntrilist
c     logical  work4   ! not used
c     logical  work5   ! not used
      common/nnsetuparrays/points(3,np_max),
     &                       centres(4,nt_max),
     &                       data(2*np_max),
     &                       vertices(4,nt_max),
     &                       neighbour(4,nt_max),
c    &                       work4(nt_max),  ! not used
c    &                       work5(np_max),  ! not used
     &                       nnn(np_max+1),
     &                       nflist(2,nwork3d),
     &                       ntrilist(nwork3d)
