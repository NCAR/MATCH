c     
c     transform resolution parameters
c
c     -triangle mf = ms, ma = 0
c     -rhombus mf=ms=ma
c     
c     -spectral resolution
      integer mf
!      parameter ( mf = 180 )             ! maximum zonal wave number
      parameter ( mf = plat )             ! maximum zonal wave number

      integer ms
!     parameter ( ms = 180 )             ! maximum polynomial order (for m=0)
      parameter ( ms = plat )             ! maximum polynomial order (for m=0)

      integer ma
      parameter ( ma =  0 )             ! maximum additional polynomial order

c     -2d array sizes (area of pentagon)

c     integer msub
c     parameter ( msub = (mf-ma)*(mf-ma)/2  ) ! area removed from parallelogram
      integer mcoe
c$$$  parameter ( mcoe = (mf+1)*(ms+1)-msub ) ! size of coefficient array
      parameter ( mcoe = (mf+1)*(ms+1)) ! size of coefficient array


      integer nwork
      parameter( nwork=3*plat*(plon+2))
      real(r8) workff(nwork)

      integer ntrig
      parameter( ntrig=3*plon/2 )

      integer idim
      parameter (idim = 18 + 2*(ms+1) + plon)
      integer iwff(idim)

      integer rdim
      parameter (rdim = 2*mcoe*plat + 5*plat + 3*plon/2+1)
      real(r8) rwff(rdim)

      integer wdim
      parameter (wdim = (mf+1)*(ms+2))
      real(r8) wkalff(wdim)

      real(r8) rlff(plat)
      real(r8) piff

      common /fixflcmr/ rwff, wkalff, workff, rlff, piff
      common /fixflcmi/ iwff
