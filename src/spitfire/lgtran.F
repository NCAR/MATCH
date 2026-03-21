      module lgtran

      use precision

      contains
c
c This is the source code for an spherical harmonic transform
c package. All subroutine names begin with LG (for Legendre).
c
c The same form of truncation as in CCM1 is supported by LGTRAN. Truncation 
c may be on an arbitrary pentagon which can reduce to a triangle or rhombus.
c
c The spectral resolution is controlled by 3 parameters:
c       1) MF  - the maximum zonal wave number
c       2) MS  - the maximum polynomial order for zonal wave number zero
c       3) MA  - the maximum additional polynomial order for other waves
c For triangular truncation: MF = MS, MA = 0.
c For rhomboidal truncation: MF = MS = MA.
c Note: MF and MS correspond to the CCM1 truncation parameters M and N, 
c       however MA is K-N in the CCM1 notation.
c
c The grid resolution is controlled by 2 parameters:
c       1) NLON - the number of distinct longitudes
c       2) NLAT - the number of latitudes.
c Although the grid and spectral parameters are specified independantly,
c specific relationships are usually required to maintain the desirable 
c properities of the transform procedure. The actual locations of the 
c latitudes are determined by LGGAUS as the roots of the Legendre 
c polynomial of order NLAT/2. 
c
c There are 2 "permanent" workspaces used by the package. These workspaces 
c are set up by "lginit" and are passed to all of the user callable routines.
c 1. integer workspace ("iwork"), 
c       length: 18 + 2*(ms+1) + nlon
c 2. real(r8) workspace ("rwork"),
c       length: ((mf+1)*(ms+1) - (mf-ma)*(mf-ma+1)/2) * nlat 
c               + 5*nlat + 3*nlon/2+1
c
c A temporary space is required during the initialization of the 
c package and each transform routine requires workspace for the fft's
c 1. initialization real(r8) workspace ("wkalf")
c       length: (mf+1)*(ms+2)
c 2. fft real(r8) workspace ("work")
c       length: 3*(nlon+2)*nlat
c
c The grid to spectral routines expect 3d arrays (nlon,nlat,nlev) as input.
c The spectral to grid routines return 3d arrays (nlon,nlat,nlev) as output.
c The latitudes are assumed to be ordered south to north. 
c nlev is an arbitrary 3rd dimension which may be 1.
c
c Prior to performing any transforms, an initialization routine (LGINIT) must 
c be called once to set the values in the integer and real(r8) work spaces
c
c Subroutines are provided to perform all of the usual transform operations
c found in global models and diagnositcs:
c  1) forward and inverse transform
c  2) gradients
c  3) divergence of a flux
c  4) u,v wind components to vorticity, divergence coefficients and back.
c
c The entry points are listed below by function. the first routine is the
c the normal user entry point and the additonal routines are internal to
c the package, but may be usefull on occasion. The arguments are documented
c in each routine.
c
c workspace size determination
c     subroutine lgdim( nlon, nlat, mf, ms, ma, 
c    $     liwork, lrwork, lwinit, lwfft)
c     
c
c initialization
c     SUBROUTINE LGINIT( nlon, nlat, mf, ms, ma, idim, rdim, wdim
c    $                 , iwork, rwork, wkalf, rlat)
c     SUBROUTINE LGINIT2
c     SUBROUTINE LGGAUS
c     SUBROUTINE LGORD 
c     SUBROUTINE LGNDRE
c     SUBROUTINE LGSTOR
c
c grid point to spectral coefficient transform for a scalar field
c     SUBROUTINE LGGTS( GRID, SPECT, NLEV, iwork, rwork, WORK )
c
c spectral coefficient to grid point transform for a scalar field
c     SUBROUTINE LGSTG( SPECT, GRID, NLEV, iwork, rwork, WORK )
c
c grid point velocity component to spectral vorticity and divergence coefficients
c     SUBROUTINE LGVGTS( UG, VG, ZS, DS, NLEV, A, iwork, rwork, WORK )
c
c spectral vorticity and divergence coefficients to grid point velocity components
c     SUBROUTINE LGVSTG( ZS, DS, UG, VG, NLEV, A, iwork, rwork, WORK )
c
c grid point flux components to spectral coefficients of the flux divergence
c     SUBROUTINE LGDGTS( UG, VG, DS, NLEV, A, iwork, rwork, WORK )
c
c spectral coefficients of a scalar field to grid point values of its gradient
c     SUBROUTINE LGGSTG( SP, GI, GJ, NLEV, A, iwork, rwork, WORK )
c
c
c
c
      subroutine lgdim( nlon, nlat, mf, ms, ma, 
     $     liwork, lrwork, lwinit, lwfft)
c+
c Compute the required dimenstions for the work spaces.
c-

      use precision

      implicit none
      integer mcoe
c+
c Formal parameters, input
c-
      integer
     $     nlon,                ! number of distinct longitudes on grid
     $     nlat,                ! number of latitudes on grid
c+
c 
c-
     $     mf,                  ! maximum zonal wave number
     $     ms,                  ! maximum polynomial order (for m=0)
     $     ma                   ! maximum additional polynomial order
c+
c Formal parameters, output
c-
      integer
     $     liwork,              ! length of the integer work space
     $     lrwork,              ! length of the real(r8) work space
     $     lwinit,              ! length of the initialization work space 
     $     lwfft                ! length of the fft work space 
c+
c Action
c-
      liwork = 18 + 2*(ms+1) + nlon

c$$$      msub = (mf-ma)*(mf-ma+1)/2 ! area removed from parallelogram
c$$$      mcoe = (mf+1)*(ms+1)-msub  ! size of coefficient array
      mcoe = (mf+1)*(ms+1)      ! size of coefficient array
      lrwork = mcoe*nlat + 5*nlat + 3*nlon/2+1

      lwinit = (mf+1)*(ms+2)
      lwfft  = 2*(nlon+2)*nlat
c+
c Return to caller
c-
!      return
      end subroutine lgdim


      subroutine lginit( nlon, nlat, mf, ms, ma, idim, rdim, wdim
     $     , iwork, rwork, wkalf, rlat)
c+
c Lginit initializes the Legendre and Fourier transform arrays in the
c integer and real(r8) work spaces.
c 
c Generate the gauss latitudes and return them to the calling
c program in rlat.
c    - nlat must be even
c    - gauss quadrature is used for integration
c-

      use precision

      implicit none
      integer mcoe
      integer i
c+
c Formal parameters, input
c-
      integer
     $     nlon,                ! number of distinct longitudes on grid
     $     nlat,                ! number of latitudes on grid
     $     mf,                  ! maximum zonal wave number
     $     ms,                  ! maximum polynomial order (for m=0)
     $     ma,                  ! maximum additional polynomial order
     $     idim,                ! length of the integer work space
     $     rdim,                ! length of the real(r8) work space
     $     wdim                 ! length of the temporary space 
c+
c Formal parameters, output
c-
      integer
     $     liwork,              ! length used in the integer work space
     $     lrwork,              ! length used in the real(r8) work space
     $     lwkalf,              ! length used in the temporary space 
     $     iwork(idim)          ! integer work space

      real(r8)
     $     rwork(rdim),         ! real(r8) work space
     $     wkalf(wdim),         ! temporary (initialization) space, trash
     $     rlat(nlat)           ! latitudes in radians

c+
c Local storage
c-
      integer                   ! pointers into workspaces
     $     pjmin, pjmax, pifax, palf, pdlf, pcolat, pcosc, psinc, 
     $     pgwt, pwos2, ptrigs

      logical debug
c+
c Action
c-
c$$$#ifdef DEBUG
c$$$      debug = .true.
c$$$#else
      debug = .false.
c$$$#endif
c+
c Determine whether this is a legitimate truncation
c-
      if (ma .gt. mf) then
         print *, 'lginit: mf must be >= than ma'
         print *, 'mf, ma: ', mf, ma
c$$$         stop
         return
      end if
      if (mf .gt. ms+ma) then
         print *, 'lginit: ms+ma must be >= than mf'
         print *, 'mf, ms, ma, ms+ma: ', mf, ms, ma, ms+ma
c$$$         stop
         return
      end if
c+
c Determine the required workspace lengths and check the dimensions.
c-
      liwork = 18 + 2*(ms+1) + nlon
      if (liwork .gt. idim) then
         print *, 'lginit: integer workspace is too small'
         print *, 'lginit: dimension is ', idim
         print *, 'lginit: required size is ', liwork
c$$$         stop
         return
      end if

c$$$      msub = (mf-ma)*(mf-ma+1)/2 ! area removed from parallelogram
c$$$      mcoe = (mf+1)*(ms+1)-msub  ! size of coefficient array
      mcoe = (mf+1)*(ms+1)     ! size of coefficient array
      lrwork = mcoe*nlat + 5*nlat + 3*nlon/2+1
         
      if (lrwork .gt. rdim) then
         print *, 'lginit: real(r8) workspace is too small'
         print *, 'lginit: dimension is ', rdim
         print *, 'lginit: required size is ', lrwork
c$$$         stop
         return
      end if

      lwkalf = (mf+1)*(ms+2)
      if (lwkalf .gt. wdim) then
         print *, 'lginit: real(r8) workspace is too small'
         print *, 'lginit: dimension is ', wdim
         print *, 'lginit: required size is ', lwkalf
c$$$         stop
         return
      end if
c+
c Set the dimensions and pointers in the integer workspace
c-
      iwork( 1) = nlon
      iwork( 2) = nlat
      iwork( 3) = mf
      iwork( 4) = ms
      iwork( 5) = ma
      iwork( 6) = mcoe
      iwork( 7) = 19                         ! jmin
      iwork( 8) = iwork(7) + ms+1            ! jmax
      iwork( 9) = iwork(8) + ms+1            ! ifax
      iwork(10) = 1                          ! alf
      iwork(11) = iwork(10) + mcoe*nlat/2 ! dlf
      iwork(12) = iwork(11) + mcoe*nlat/2 ! colat
      iwork(13) = iwork(12) + nlat           ! cosc
      iwork(14) = iwork(13) + nlat           ! sinc
      iwork(15) = iwork(14) + nlat           ! gwgt
      iwork(16) = iwork(15) + nlat           ! wos2
      iwork(17) = iwork(16) + nlat           ! trigs
      iwork(18) = 1
      if (debug) then
         print *, 'lginit: iwork(1-18) ', (iwork(i), i=1, 18)
      end if
c+
c Set the pointers from the workspace
c-
c$$$      nlon  = iwork( 1)
c$$$      nlat  = iwork( 2)
c$$$      mf    = iwork( 3)
c$$$      ms    = iwork( 4)
c$$$      ma    = iwork( 5)
c$$$      mcoe  = iwork( 6)
      pjmin = iwork( 7)
      pjmax = iwork( 8)
      pifax = iwork( 9)
      palf  = iwork(10)
      pdlf  = iwork(11)
      pcolat= iwork(12)
      pcosc = iwork(13)
      psinc = iwork(14)
      pgwt  = iwork(15)
      pwos2 = iwork(16)
      ptrigs= iwork(17)
c+
c Set arrays in workspaces
c-
      call lginit2( wkalf, rlat,
     $     nlon, nlat, mf, ms, ma, mcoe,
     $     iwork(pjmin), iwork(pjmax), iwork(pifax),
     $     rwork(palf),  rwork(pdlf),  rwork(pcolat), rwork(pcosc),
     $     rwork(psinc), rwork(pgwt),  rwork(pwos2),  rwork(ptrigs))
c+
c Return to caller
c-
      return
      end subroutine lginit

      subroutine lginit2( wkalf, rlat,
     $     nlon, nlat, mf, ms, ma, mcoe,
     $     jmin, jmax, ifax,
     $     alf, dlf, colat, cosc,
     $     sinc, gwt, wos2, trigs
     $     )
c+
c Lginit2 initializes the Legendre and Fourier transform arrays.
c-

      use precision

      implicit none
      integer j, l, n, nlt
      real(r8) pio2
c+
c Formal parameters, input
c-
      integer
     $     nlon,                ! number of distinct longitudes on grid
     $     nlat,                ! number of latitudes on grid
     $     mf,                  ! maximum zonal wave number
     $     ms,                  ! maximum polynomial order (for m=0)
     $     ma,                  ! maximum additional polynomial order
     $     mcoe                 ! number of spherical harmonic coefficients
c+
c Formal parameters, output
c-
      integer
     $     jmin(0:ms),          ! pointer to min zonal wave number for l=n-m
     $     jmax(0:ms),          ! pointer to max zonal wave number for l=n-m
     $     ifax(nlon)           ! prime factors of nlon

      real(r8)
     $     wkalf((mf+1)*(ms+2)), ! temporary (initialization) space, trash
     $     rlat(nlat),          ! latitudes
     $     alf(mcoe,*),         ! associated legendre functions
     $     dlf(mcoe,*),         ! derivatives of alf
     $     colat(nlat),         ! colatitude
     $     cosc(nlat),          ! cos(colatitude)
     $     sinc(nlat),          ! sin(colatitude)
     $     gwt(nlat),           ! gauss weights
     $     wos2(nlat),          ! gauss weights / sin(colat)**2
     $     trigs(3*nlon/2+1)    ! trig functions for fft.
c+
c Local storage
c-
      logical debug
c+
c Action
c-
c$$$#ifdef DEBUG
c$$$      debug = .true.
c$$$#else
      debug = .false.
c$$$#endif
C
C    -SET UP FACTORS AND TRIG FUNCTIONS FOR FFT
C
c      CALL FAX   ( ifax , NLON, 3 )
c      CALL FFTRIG( trigs, NLON, 3 )
      call set99(trigs, ifax, nlon)
      if (debug) then
         print *, 'lginit2: ifax ', ifax
         print *, 'lginit2: trigs ', trigs
      end if
C
C    -set the index arrays
      JMIN(0) = 1
      JMAX(0) = 1 + MF
      DO 10 L=1,MS+MA-MF
         JMIN(L) = JMIN(L-1) + MF + 1
         JMAX(L) = JMIN(L) + MF
   10 CONTINUE
      N = MF
      DO 20 L=MS+MA-MF+1,MS
         N = N - 1
c$$$         JMIN(L) = JMAX(L-1) + 1
         JMIN(L) = JMIN(L-1) + MF + 1
         JMAX(L) = JMIN(L) + N
   20 CONTINUE
      if (debug) then
         print *, 'lginit2: jmin ', jmin
         print *, 'lginit2: jmax ', jmax
      end if
C
C    -find the Gauss colatitudes (COLAT), weights (GWT), cos(colat) (COSC),
C    -sin(colat) (SINC), and weights over sin(colatitude)**2 (WOS2)
      CALL LGGAUS( NLAT, COSC, GWT, SINC, COLAT, WOS2 )
      pio2 = 0.5*acos(-1.0)
      do j = 1, nlat
         rlat(j) = pio2 - colat(j) 
      end do
      if (debug) then
         print *, 'lginit2: rlat ', rlat
         print *, 'lginit2: cosc ', cosc
         print *, 'lginit2: gwt  ', gwt
         print *, 'lginit2: sinc ', sinc
         print *, 'lginit2: colat ', colat
         print *, 'lginit2: wos2 ', wos2
      end if
C
C
C    -set the Legendre function array for each latitiude
      nlt = nlat/2
      do J=1,NLT
C
C       -evaluate the full parallelogram of associated Legendre functions
         if (debug) then
            print *, 'lginit2: calling legendre for latitude ', j
         end if
         CALL LGNDRE( WKALF, COSC(J), SINC(J), COLAT(J), MS, MF )
C
C       -move the required funtion values to the storage array
         if (debug) then
            print *, 'lginit2: calling lgstor for latitude ', j
         end if
         CALL LGSTOR( WKALF, J,
     $     nlon, nlat, mf, ms, ma, mcoe,
     $     jmin, jmax, alf, dlf, cosc
     $     )
C
      end do
C
!     RETURN
      end subroutine lginit2

      SUBROUTINE LGGAUS( NLAT, COSC, GWT, SINC, COLAT, WOS2 )
C
C  LGGAUS finds the Gaussian latitudes by finding the roots of the
C  ordinary Legendre polynomial of degree NLAT using Newton's
C  iteration method.
C
C  On entry:
C     NLAT - the number of latitudes (degree of the polynomial)
C
C  On exit: for each Gaussian latitude
C     COSC   - cos(colatitude) or sin(latitude)
C     GWT    - the Gaussian weights
C     SINC   - sin(colatitude) or cos(latitude)
C     COLAT  - the colatitudes in radians
C     WOS2   - Gaussian weight over sin**2(colatitude)
C

      use precision

      implicit none

      integer nlat
      integer i, j, nzero
      real(r8) tmp

      real(r8) COSC(NLAT), GWT(NLAT), SINC(NLAT), COLAT(NLAT)
     $        , WOS2(NLAT)

      double precision
     $     xlim, csc, gw, pi,
     $     a, b, c, d, delta, fi, fi1, g, gm, gp, gt
C
C-----------------------------------------------------------------------
C
C    -convergence criterion for iteration of cos latitude
      XLIM  = 1.0D-13
C
C    -the number of zeros between pole and equator
      NZERO = NLAT/2
C
C    -set first guess for cos(colat)
      PI = 3.141592653589793
      DO 10 I=1,NZERO
         COSC(I) = SIN( (I-0.5)*PI/NLAT + PI*0.5 )
   10 CONTINUE
C
C    -constants for determining the derivative of the polynomial
      FI  = NLAT
      FI1 = FI+1.0
      A   = FI*FI1 / SQRT(4.0*FI1*FI1-1.0)
      B   = FI1*FI / SQRT(4.0*FI*FI-1.0)
C
C    -loop over latitudes, iterating the search for each root
      DO 40 I=1,NZERO
         J=0
C
C       -determine the value of the ordinary Legendre polynomial for
C       -the current guess root
         csc = cosc(i)
   30    CALL LGORD( G, CSC, NLAT )
C
C       -determine the derivative of the polynomial at this point
         CALL LGORD( GM, CSC, NLAT-1 )
         CALL LGORD( GP, CSC, NLAT+1 )
         GT = (CSC*CSC-1.0) / (A*GP-B*GM)
C
C       -update the estimate of the root
         DELTA   = G*GT
         CSC = CSC - DELTA
C
C       -if convergence criterion has not been met, keep trying
         J = J+1
         IF( ABS(DELTA).GT.XLIM ) then
            if (j .ge. 10) then
c$$$               stop 'no conv'
               print *, 'gauss latitude iterations not converging'
            else
               GO TO 30
            end if
         end if
C        PRINT*,' LAT NO.',I,J,' ITERATIONS'
C
C       -determine the Gaussian weights
         C  = 2.0 *( 1.0-CSC*CSC )
         CALL LGORD( D, CSC, NLAT-1 )
         D  = D*D*FI*FI
         GW = C *( FI-0.5 ) / D
c+
c Set output arrays, keep all arithmetic in double precision
c-
         cosc(i) = csc
         gwt(i)  = gw
         csc = acos(csc)
         colat(i)= csc
         csc = sin(csc)
         sinc(i) = csc
         wos2(i) = gw / (csc*csc)
   40 CONTINUE
C
C    -if NLAT is odd, set values at the equator
      IF( MOD(NLAT,2) .NE. 0 ) THEN
         I       = NZERO+1
         csc     = 0.0
         C       = 2.0
         CALL LGORD( D, csc, NLAT-1 )
         D       = D*D*FI*FI
         COSC(I) = 0.0
         GWT(I)  = C *( FI-0.5 ) / D
         COLAT(I)= PI*0.5
         SINC(I) = 1.0
         WOS2(I) = GWT(I)
      END IF
C
C    -determine the southern hemisphere values by symmetry
      DO 60 I=NLAT-NZERO+1,NLAT
         COSC(I) =-COSC(NLAT+1-I)
         GWT(I)  = GWT(NLAT+1-I)
         COLAT(I)= PI-COLAT(NLAT+1-I)
         SINC(I) = SINC(NLAT+1-I)
         WOS2(I) = WOS2(NLAT+1-I)
   60 CONTINUE

c     swap all the values so it is ordered south to north
      do i = 1,nlat/2
         j = nlat+1-i
         tmp = cosc(i)
         cosc(i) = cosc(j)
         cosc(j) = tmp
         tmp = wos2(i)
         wos2(i) = wos2(j)
         wos2(j) = tmp
         tmp = sinc(i)
         sinc(i) = sinc(j)
         sinc(j) = tmp
         tmp = colat(i)
         colat(i) = colat(j)
         colat(j) = tmp
         tmp = gwt(i)
         gwt(i) = gwt(j)
         gwt(j) = tmp
      end do
      
c     PRINT*,'NLAT=',NLAT
c     PRINT*,'COLATS'
c     PRINT 101,(I,COLAT(I),COLAT(I)*180./PI,I=1,NLAT)
  101 FORMAT(1X,I3,F16.12,2X,F8.2)
  102 FORMAT(1X,I3,F16.12,2X,F16.12)
c     PRINT*,'COS(COLAT), SIN(COLAT)'
c     PRINT 102,(I,COSC(I),SINC(I),I=1,NLAT)
c     PRINT*,'WEIGHT, GWT/COS**2'
c     PRINT 102,(I,GWT(I),WOS2(I),I=1,NLAT)
C
      RETURN
      END subroutine lggaus

      SUBROUTINE LGORD( F, COSC, N )
C
C  LGORD calculates the value of an ordinary Legendre polynomial at a
C  latitude.
C
C  On entry:
C     COSC - cos(colatitude)
C     N      - the degree of the polynomial
C
C  On exit:
C     F      - the value of the Legendre polynomial of degree N at
C              latitude asin(COSC)
C
C---------------------------------------------------------------------------

      use precision

      implicit double precision (a-h, o-z)
C---------------------------------------------------------------------------
C
C    -determine the colatitude
      COLAT = ACOS(COSC)
C
      C1 = SQRT(2.0)
      DO 20 K=1,N
         C1 = C1 * SQRT( 1.0 - 1.0/(4*K*K) )
   20 CONTINUE
C
      FN = N
      ANG= FN * COLAT
      S1 = 0.0
      C4 = 1.0
      A  =-1.0
      B  = 0.0
      DO 30 K=0,N,2
         IF (K.EQ.N) C4 = 0.5 * C4
         S1 = S1 + C4 * COS(ANG)
         A  = A + 2.0
         B  = B + 1.0
         FK = K
         ANG= COLAT * (FN-FK-2.0)
         C4 = ( A * (FN-B+1.0) / ( B * (FN+FN-A) ) ) * C4
   30 CONTINUE
      F = S1 * C1
C
      RETURN
      END subroutine lgord

      SUBROUTINE LGNDRE ( P, COSC, SINC, COLAT, NS, NF )
C
C  LGNDRE generates the values of the Legendre polynomials for a
C  latitude. A full parallelogram is returned.
C
C  On entry:
C     COSC  - cos(colatitude)
C     SINC  - sin(colatitude)
C     COLAT - colatitude (radians)
C     NS    - maximum polynomial order for zonal wave zero
C     NF    - maximum zonal wave number
C
C  On exit:
C     P     - polynomial values for the requested latitude. Values are
C             returned for polynomials of order n = (m , m+NS+1)
C             for zonal wave number m = (0, NF)
C

      use precision

      implicit none
      integer ns, nf
      real(r8) P(0:NF,0:NS+1)
      real(r8) c1
      integer k, l, m, n
      real(r8) colat, cosc, sinc
      real(r8) c3, c4, c5, c6, c7, c8
      real(r8) a, b, c, d, e, fk, fn
      real(r8) s1, s2
      real(r8) ang
      real(r8) fm, fm1, fm2, fm3
C----------------------------------------------------------------------
C
C    -m=0, n=0
      C1     = SQRT( 2.0 )
      P(0,0) = 1.0 / C1
C
C    -set polynomial values for m=0 and m=1
      DO 20 N=1,NS+2
         FN  = N
         C1  = C1 * SQRT( 1.0 - 1.0 / (4.0*FN*FN) )
         C3  = C1 / SQRT( FN * (FN+1.0) )
         ANG = FN * COLAT
         S1  = 0.0
         S2  = 0.0
         C4  = 1.0
         C5  = FN
         A   =-1.0
         B   = 0.0
         DO 10 K=0,N,2
            S2 = S2 + C5 * SIN(ANG) * C4
            IF (K.EQ.N) C4 = 0.5 * C4
            S1 = S1 + C4 * COS(ANG)
            A  = A + 2.0
            B  = B + 1.0
            FK = K
            ANG= COLAT * (FN-FK-2.0)
            C4 = ( A * (FN-B+1.0) / ( B * (FN+FN-A) ) ) * C4
            C5 = C5 - 2.0
   10    CONTINUE
         IF( N.LE.NS+1 ) P(0,N) = S1 * C1
         P(1,N-1) = S2 * C3
   20 CONTINUE
C
C    -finished if only m=0 and m=1 required
      IF( NF.EQ.1 ) GO TO 60
C
C    -loop over remaining zonal wave numbers
      DO 50 M=2,NF
         FM  = M
         FM1 = FM - 1.0
         FM2 = FM - 2.0
         FM3 = FM - 3.0
         C6  = SQRT( ( 2.0 * FM + 1.0 ) / ( 2.0 * FM ) )
C
C       -set polynomial for n=m
         P(M,0) = C6 * SINC * P(M-1,0)
C
C       -loop over remaining meridional modes
         DO 40 L=1,NS
            FN = M + L
            C7 = (FN*2.0+1.0) / (FN*2.0-1.0)
            C8 = (FM1+FN) / ( (FM+FN) * (FM2+FN) )
            C  = SQRT( (FN*2.0+1.0) / (FN*2.0-3.0) * C8 * (FM3+FN) )
            D  = SQRT( C7 * C8 * (FN-FM1) )
            E  = SQRT( C7 * (FN-FM) / (FN+FM) )
            P(M,L) = C * P(M-2,L  )
     $             - D * P(M-2,L+1) * COSC
     $             + E * P(M  ,L-1) * COSC
   40    CONTINUE
C
C       -set polynomial for n=m+NS+1, used for derivatives (and m+2)
C       -Belousov equation 11
         L  = NS + 1
         FN = M + L
         A  = SQRT( (FN*FN-0.25) / (FN*FN-FM*FM) )
         B  = SQRT( (2.0*FN+1.0) * (FN-FM-1.0) * (FN+FM1)
     $              / ( (2.0*FN-3.0) * (FN-FM) * (FN+FM) ) )
         P(M,L) = 2.0 * A * COSC * P(M,L-1) - B * P(M,L-2)
   50 CONTINUE
 
C
C    -set polynomial values less than machine accuracy to zero
   60 CONTINUE
C     DO 80 M=0,NF
C        DO 70 L=0,NS+1
C           IF( ABS(P(M,L)).LT.1.0E-14 ) P(M,L)=0.0
C  70    CONTINUE
C  80 CONTINUE
C
      RETURN
      END subroutine lgndre

      SUBROUTINE LGSTOR( P, JLAT,
     $     nlon, nlat, mf, ms, ma, mcoe,
     $     jmin, jmax, alf, dlf, cosc
     $     )

      use precision

      implicit none
      integer m
      integer n
      integer jlat
      integer l
      integer j
      integer mm
      integer nn

C
C  LGSTOR moves the required values of the associated Legendre funtions
C  for a latitude from the working array to the storage array in common
C  and computes and stores the derivative functions values. The relation
C  used for the derivative function requires the extra function order in
C  the input array.
C
C  On entry:
C     P     - full parallelogram of function values
C     JLAT  - latitude index (from north pole)
C
C  On exit: common block arrays ALF and DLF set for latitude JLAT
C
      integer
     $     nlon,                ! number of distinct longitudes on grid
     $     nlat,                ! number of latitudes on grid
     $     mf,                  ! maximum zonal wave number
     $     ms,                  ! maximum polynomial order (for m=0)
     $     ma,                  ! maximum additional polynomial order
     $     mcoe                 ! number of spherical harmonic coefficients

      integer
     $     jmin(0:ms),          ! pointer to min zonal wave number for l=n-m
     $     jmax(0:ms)           ! pointer to max zonal wave number for l=n-m

      real(r8)
     $     alf(mcoe,nlat/2),    ! associated legendre functions
     $     dlf(mcoe,nlat/2),    ! derivatives of alf
     $     cosc(nlat)           ! cos(colatitude) or sin(latitude)
C
      real(r8) P(0:MF,0:MS+1)
C
C    -statement function for epsilon (for all n>0)
      real(r8) epsil
      EPSIL(MM,NN) = SQRT( FLOAT(NN*NN-MM*MM) / FLOAT(4*NN*NN-1) )
C
C    -store the required function values
      DO 110 L=0,MS
         M = -1
         DO 100 J=JMIN(L),JMAX(L)
            M = M + 1
            ALF(J,JLAT) = P(M,L)
  100    CONTINUE
  110 CONTINUE
C
C    -compute the derivative functions
C
C    -special section for l=0, (n=m), as n-1 does not exist
      DLF(1,JLAT) = 0.0
      L = 0
         M = 0
         DO 120 J=JMIN(L)+1,JMAX(L)
            M = M + 1
            N = M + L
            DLF(J,JLAT) =
     $                  -   N * COSC(JLAT) * P(M,L)
  120    CONTINUE
C
C    -all other m,n
      DO 140 L=1,MS
         M = -1
         DO 130 J=JMIN(L),JMAX(L)
            M = M + 1
            N = M + L
            DLF(J,JLAT) = (2*N+1) * EPSIL(M,N)  * P(M,L-1)
     $                  -   N * COSC(JLAT) * P(M,L)
  130    CONTINUE
  140 CONTINUE
C
      RETURN
      END subroutine lgstor





      subroutine lggts( grid, spect, nlev, iwork, rwork, work )
C
C  LGGTS transforms a 3d grid array to an array of spherical harmonic
C  coefficients.
C
C  On entry:
C     GRID:   a 3d array of grid values with wrap around points
C     NLEV:   the number of vertical levels
C     WORK:   work space for the Fourier transform, length:2*(NLON+2)*NLAT
C
C  ON EXIT:
C     SPECT:  an array of spherical harmonic coefficients
C-----------------------------------------------------------------------
C

      use precision

      implicit none 
      integer nlon
      integer nlat
      integer ms
      integer mcoe
      integer ifft
      integer nlev

      real(r8) grid(*), work(*), spect(*)

      integer iwork(*)
      real(r8) rwork(*)

      integer                   ! pointers into workspaces
     $     pjmin, pjmax, pifax, palf, pdlf, psinc, 
     $     pgwt, pwos2, ptrigs
C
C-----------------------------------------------------------------------
c+
c Set the pointers from the workspace
c-
      nlon  = iwork( 1)
      nlat  = iwork( 2)
      ms    = iwork( 4)
      mcoe  = iwork( 6)
      pjmin = iwork( 7)
      pjmax = iwork( 8)
      pifax = iwork( 9)
      palf  = iwork(10)
      pdlf  = iwork(11)
      psinc = iwork(14)
      pgwt  = iwork(15)
      pwos2 = iwork(16)
      ptrigs= iwork(17)
c+
c do the transform
c-
      ifft = (nlon+2)*nlat+1
      call lggts2( grid, spect, nlev, nlon, nlat, ms, mcoe,
     $     work(1), work(ifft), iwork(pifax), rwork(ptrigs),
     $     iwork(pjmin), iwork(pjmax), rwork(palf),  rwork(pdlf),  
     $     rwork(psinc), rwork(pgwt),  rwork(pwos2))
c
      return
      end subroutine lggts


      subroutine lggts2( grid, spect, nlev, nlon, nlat, ms, mcoe,
     $     gwrk, work, ifax, trigs,
     $     jmin, jmax, alf, dlf, 
     $     sinc, gwt, wos2)


      use precision

      implicit none
      integer ifax(*)
      integer k
      integer j
      integer i
      integer l
      integer m1
      integer m2
      integer m
      integer jlat
      integer lats
      integer nlev
      integer jlt
      real(r8) c
      real(r8) trigs(*)
      
c
c  lggts2 transforms a grid array to an array of spherical harmonic
c  coefficients.
c
C-----------------------------------------------------------------------
      integer
     $     nlon,                ! number of distinct longitudes on grid
     $     nlat,                ! number of latitudes on grid
     $     ms,                  ! maximum polynomial order (for m=0)
     $     mcoe                 ! number of spherical harmonic coefficients

      integer
     $     jmin(0:ms),          ! pointer to min zonal wave number for l=n-m
     $     jmax(0:ms)           ! pointer to max zonal wave number for l=n-m

      real(r8)
     $     alf(mcoe,*),         ! associated legendre functions
     $     dlf(mcoe,*),         ! derivatives of alf
     $     sinc(nlat),          ! sin(colatitude)
     $     gwt(nlat),           ! gauss weights
     $     wos2(nlat)           ! gauss weights / sin(colat)**2
C
      real(r8)
     $     grid(nlon,nlat,nlev),        ! grid array
     $     spect(2,mcoe,nlev),          ! spherical harmonic coeffs
     $     gwrk(nlon+2,nlat),           ! grid work space for fft
     $     work(nlon+2,nlat)            ! more work space

c      real(r8) trigspjr(3*nlon/2+1)
c      integer ifaxpjr(13)
c      real(r8) gwrkpjr(nlon+2,nlat)
c      real(r8) workpjr(nlon+2,nlat)
      integer isign
C
C-----------------------------------------------------------------------
c+
c zero the coefficient array
c-
c      do j = 1,nlat
c         write (6,*) 'lggts2 d3:', sinc(j), wos2(j), gwt(j)
c      end do

      do k = 1, nlev
         do j = 1, mcoe
            spect(1,j,k) = 0.
            spect(2,j,k) = 0.
         end do
      end do
c+
c loop over outer dimension of variable
c-
      do k = 1, nlev
c+
c copy input grid array to grid work space with 2 extra longitude points
c-
         do j = 1, nlat
            do i = 1, nlon
               gwrk(i,j) = grid(i,j,k)
            end do
         end do
c+
c fourier transform all latitudes, result still in gwrk
c-
#ifndef PJR
c         open(file='/tmp/pjr/fixflx02/binarydata',unit=78,
c     $        form='unformatted')
c         write (78) (ifax(i),i=1,13)
c         write (78) (trigs(i),i=1,3*nlon/2 + 1)
c         write (78) gwrk
c         write (6,*) ' binary data dumped '
#else
c         open(file='/tmp/pjr/fixflx02/binarydata',unit=78,
c     $        form='unformatted') 
c         read (78) (ifax(i),i=1,13)
c         read (78) (trigs(i),i=1,3*nlon/2 + 1)
c         read (78) gwrk
c         write (6,*) ' binary data read '
#endif
            
c         write (77,*) ' lggts2 d1 ', nlon, nlat, k
c         write (77,*) ' ifax '
c         write (77,76) (ifax(i),i=1,13)
c   76    format (10i5)
c         write (77,*) ' trigs '
c         write (77,77) (trigs(i),i=1,3*nlon/2 + 1)
c         write (77,*) ' gwrk '
c         write (77,77) gwrk
c   77    format (1p,10e15.5)
c         do i = 1,13
c            ifaxpjr(i) = ifax(i)
c         end do
c         do i = 1,3*nlon/2+1
c            trigspjr(i) = trigs(i)
c         end do
c         do j = 1,nlat
c            do i = 1,nlon
c               gwrkpjr(i,j) = gwrk(i,j)
c            end do
c         end do

c         call fft991 (gwrk, work, trigs, ifax, 1, nlon+2,
c     $        nlon, nlat, -1)
         isign = -1
c        write (6,*) isign
         call fft991 (gwrk, work, trigs, ifax, 1, nlon+2,
     $        nlon, nlat, isign)
c         write (77,*) ' answers '
c         write (77,77) gwrk
c+
c loop over latitude. Input/output grid data is south to north. 
c Legendre functions are for NH from north to south.
c-
         do lats=1,nlat
            jlat = nlat+1-lats

            if (jlat .le. nlat/2) then
               jlt = jlat
            else
               jlt = nlat +1 -jlat
            end if
c+
c Multiply Fourier coefficient by the Gaussian weight for this
c latitude and the appropriate ALF value and accumulate
c the results. SPECT is complete when the contributions for all
c latitudes has been added. "c" factor accounts for symmetry properties 
c of the associated legendre functions, stored for NH only.
c-
            do l = 0, ms
               m = -1
               if ((jlat .le. nlat/2) .or. (mod(l,2) .eq. 0)) then
                  c = +1.0
               else
                  c = -1.0
               end if
               do j = jmin(l), jmax(l)
                  m  = m + 1
                  m1 = 2*m +1
                  m2 = 2*m +2
                  spect(1,j,k) = spect(1,j,k)
     $                 + c * gwt(jlt) * alf(j,jlt) * gwrk(m1,jlat)
                  spect(2,j,k) = spect(2,j,k)
     $                 + c * gwt(jlt) * alf(j,jlt) * gwrk(m2,jlat)
               end do
            end do
         end do
c         write (6,*) ' lggts2 d2 ', spect(1,1,k), spect(2,1,k),
c     $        spect(1,mcoe/2,k), spect(2,mcoe/2,k)
      end do
C
      return
      end subroutine lggts2



      subroutine lgstg( spect, grid, nlev, iwork, rwork, work )


      use precision

      implicit none
      integer ifft, mcoe, ms, nlat, nlon
      integer nlev
C
C  LGSTG transforms an array of spherical harmonic coefficients to a
C  3d grid array
C
C  On entry:
C     SPECT:  an array of spherical harmonic coefficients.
C     NLEV:   the number of vertical levels
C     WORK:   work space for the Fourier transform, length:2*(nlon+2)*nlat
C
C  On exit:
C     GRID:   a 3d array of grid values
C
C-----------------------------------------------------------------------
C
      real(r8) spect(*), grid(*), work(*)
c
      integer iwork(*)
      real(r8) rwork(*)

      integer                   ! pointers into workspaces
     $     pjmin, pjmax, pifax, palf, pdlf, psinc, 
     $     pgwt, pwos2, ptrigs
c
c-----------------------------------------------------------------------
c+
c Set the pointers from the workspace
c-
      nlon  = iwork( 1)
      nlat  = iwork( 2)
      ms    = iwork( 4)
      mcoe  = iwork( 6)
      pjmin = iwork( 7)
      pjmax = iwork( 8)
      pifax = iwork( 9)
      palf  = iwork(10)
      pdlf  = iwork(11)
      psinc = iwork(14)
      pgwt  = iwork(15)
      pwos2 = iwork(16)
      ptrigs= iwork(17)
c+
c do the transform
c-
      ifft = (nlon+2)*nlat+1
      call lgstg2( spect, grid, nlev, nlon, nlat, ms, mcoe,
     $     work(1), work(ifft), iwork(pifax), rwork(ptrigs),
     $     iwork(pjmin), iwork(pjmax),
     $     rwork(palf),  rwork(pdlf),  
     $     rwork(psinc), rwork(pgwt),  rwork(pwos2))
c
      return
      end subroutine lgstg


      subroutine lgstg2( spect, grid, nlev, nlon, nlat, ms, mcoe,
     $     gwrk, work, ifax, trigs,
     $     jmin, jmax,
     $     alf, dlf, 
     $     sinc, gwt, wos2)


      use precision

      implicit none
      real(r8) c
      integer i, j, jlt, jlat, k, l, lats, m, m1, m2
      integer ifax(*)
      real(r8) a
      integer nlev
      real(r8) trigs(*)
c
c  lgstg2 transform an array of spherical harmonic coefficients to a 
c  3d grid array (no wrap points).
c
c-----------------------------------------------------------------------
      integer
     $     nlon,                ! number of distinct longitudes on grid
     $     nlat,                ! number of latitudes on grid
     $     ms,                  ! maximum polynomial order (for m=0)
     $     mcoe                 ! number of spherical harmonic coefficients

      integer
     $     jmin(0:ms),          ! pointer to min zonal wave number for l=n-m
     $     jmax(0:ms)           ! pointer to max zonal wave number for l=n-m

      real(r8)
     $     alf(mcoe,*),         ! associated legendre functions
     $     dlf(mcoe,*),         ! derivatives of alf
     $     sinc(nlat),          ! sin(colatitude)
     $     gwt(nlat),           ! gauss weights
     $     wos2(nlat)           ! gauss weights / sin(colat)**2
C
      real(r8)
     $     grid(nlon,nlat,nlev),        ! grid array
     $     spect(2,mcoe,nlev),          ! spherical harmonic coeffs
     $     gwrk(nlon+2,nlat),           ! grid work space for fft
     $     work(nlon+2,nlat)            ! more work space

      integer isign

c
c-----------------------------------------------------------------------
c+
c loop over outer dimension of variable
c-
      do k = 1, nlev
c+
c zero the fourier coeficient accumulators
c-
         do j = 1, nlat
            do i = 1, nlon+2
               gwrk(i,j) = 0.0
            end do
         end do
c+
c loop over latitude. Input/output grid data is south to north. 
c Legendre functions are for NH from north to south.
c-
         do lats=1,nlat
            jlat = nlat+1-lats

            if (jlat.le.nlat/2) then
               jlt = jlat
            else
               jlt = nlat +1 -jlat
            end if
c+
c loop over all meridional modes, accumulating fourier coefficients 
c for this latitude
c-
            do l = 0, ms
               m = -1
               if ((jlat .le. nlat/2) .or. (mod(l,2) .eq. 0)) then
                  c = +1.0
               else
                  c = -1.0
               end if
               do j=jmin(l),jmax(l)
                  m  = m + 1
                  m1 = 2*m +1
                  m2 = 2*m +2
                  gwrk(m1,jlat) = gwrk(m1,jlat)
     $                 + c * alf(j,jlt) * spect(1,j,k)
                  gwrk(m2,jlat) = gwrk(m2,jlat)
     $                 + c * alf(j,jlt) * spect(2,j,k)
               end do
            end do

         end do
c+
c fourier transform all latitudes to grid values, still in gwrk 
c-
         isign = 1
         call fft991( gwrk, work, trigs, ifax, 1, nlon+2,
     $        nlon, nlat, isign)
c+
c copy output grid data removing extra points
c-
         do j = 1, nlat
            do i = 1, nlon
               grid(i,j,k) = gwrk(i,j)
            end do
         end do

      end do
c
      return
      end subroutine lgstg2

      subroutine lgvgts( ug, vg, zs, ds, nlev, a, iwork, rwork, work )
      

      use precision

      implicit none
      integer iv, ifft, mcoe, ms
      integer nlon, nlat
      integer nlev
      real(r8) a

c
C  LGVGTS transforms the u and v velocity components (in standard
C  dimensional form without cosine weighting) to the spectral
C  coefficients of vorticity and divergence.
C
C  On entry:
C     UG    - zonal winds
C     VG    - meridional winds
C     NLEV  - the number of vertical levels
C     A     - the radius of the sphere
C     WORK  - work space for the Fourier transform, length:3*(NLON+2)*NLAT
C
C  On exit:
C     ZS    - vorticity spectral coefficients
C     DS    - divergence spectral coefficients
C
C-----------------------------------------------------------------------
C
      real(r8) ug(*), vg(*), work(*), zs(*), ds(*)
c
      integer iwork(*)
      real(r8) rwork(*)

      integer                   ! pointers into workspaces
     $     pjmin, pjmax, pifax, palf, pdlf, psinc, 
     $     pgwt, pwos2, ptrigs
C
C-----------------------------------------------------------------------
c+
c Set the pointers from the workspace
c-
      nlon  = iwork( 1)
      nlat  = iwork( 2)
      ms    = iwork( 4)
      mcoe  = iwork( 6)
      pjmin = iwork( 7)
      pjmax = iwork( 8)
      pifax = iwork( 9)
      palf  = iwork(10)
      pdlf  = iwork(11)
      psinc = iwork(14)
      pgwt  = iwork(15)
      pwos2 = iwork(16)
      ptrigs= iwork(17)
c+
c do the transform
c-
      iv   = (nlon+2)*nlat+1
      ifft = (nlon+2)*nlat+1 + iv
      call lgvgts2( ug, vg, zs, ds, nlev, a, nlon, nlat, ms, mcoe,
     $     work(1), work(iv), work(ifft), iwork(pifax), rwork(ptrigs),
     $     iwork(pjmin), iwork(pjmax), rwork(palf),  rwork(pdlf),  
     $     rwork(psinc), rwork(pgwt),  rwork(pwos2))
c
      return
      end subroutine lgvgts


      subroutine lgvgts2(
     $     ug, vg, zs, ds, nlev, a, nlon, nlat, ms, mcoe,
     $     uwrk, vwrk, work, ifax, trigs,
     $     jmin, jmax, alf, dlf, 
     $     sinc, gwt, wos2)


      use precision

      implicit none
      integer i, j, jlat, jlt, k, l, lats, m, m1, m2
      real(r8) ca, cd
      real(r8) a
      integer nlev
      integer ifax(*)
      real(r8) trigs(*)

c
c  lgvgts2 transforms the u and v velocity components to the
c  spectral coefficients of vorticity and divergence.
c
c-----------------------------------------------------------------------
      integer
     $     nlon,                ! number of distinct longitudes on grid
     $     nlat,                ! number of latitudes on grid
     $     ms,                  ! maximum polynomial order (for m=0)
     $     mcoe                 ! number of spherical harmonic coefficients

      integer
     $     jmin(0:ms),          ! pointer to min zonal wave number for l=n-m
     $     jmax(0:ms)           ! pointer to max zonal wave number for l=n-m

      real(r8)
     $     alf(mcoe,nlat/2),    ! associated legendre functions
     $     dlf(mcoe,nlat/2),    ! derivatives of alf
     $     sinc(nlat),          ! sin(colatitude)
     $     gwt(nlat),           ! gauss weights
     $     wos2(nlat)           ! gauss weights / sin(colat)**2
C
      real(r8)
     $     ug(nlon,nlat,nlev),          ! zonal wind grid array
     $     vg(nlon,nlat,nlev),          ! meridional wind grid array
     $     zs(2,mcoe,nlev),             ! vorticity coeffs
     $     ds(2,mcoe,nlev),             ! divergence coeffs
     $     uwrk(nlon+2,nlat),           ! grid work space for u fft
     $     vwrk(nlon+2,nlat),           ! grid work space for v fft
     $     work(nlon+2,nlat)            ! more work space

      integer isign
C
C-----------------------------------------------------------------------
c+
c zero the coefficient arrays
c-
      do k = 1, nlev
         do j = 1, mcoe
            zs(1,j,k) = 0.
            zs(2,j,k) = 0.
            ds(1,j,k) = 0.
            ds(2,j,k) = 0.
         end do
      end do
c+
c loop over outer dimension of variable
c-
      do k = 1, nlev
c+
c copy input grid array to grid work space with 2 extra longitude points
c weight wind components by cos(latitude) / radius of sphere
c-
         do j = 1, nlat
            do i = 1, nlon
               uwrk(i,j) = ug(i,j,k) * sinc(j) / a
               vwrk(i,j) = vg(i,j,k) * sinc(j) / a
            end do
         end do
c+
c fourier transform all latitudes, result still in gwrk
c-
         isign = -1
         call fft991 (uwrk, work, trigs, ifax, 1, nlon+2,
     $        nlon, nlat, isign)
         call fft991 (vwrk, work, trigs, ifax, 1, nlon+2,
     $        nlon, nlat, isign)
c+
c loop over latitude. Input/output grid data is south to north. 
c Legendre functions are for NH from north to south.
c-
         do lats=1,nlat
            jlat = nlat+1-lats
            if( jlat.le.nlat/2 ) then
               jlt = jlat
            else
               jlt = nlat +1 -jlat
            end if
            do l = 0, ms
               m = -1
               ca = 1.0
               cd = 1.0
               if( jlat.gt.nlat/2 ) then
                  if( mod(l,2).eq.0 ) then
                     cd = -1.0
                  else
                     ca = -1.0
                  end if
               end if
               do j = jmin(l), jmax(l)
                  m  = m + 1
                  m1 = 2*m +1
                  m2 = 2*m +2
                  zs(1,j,k) = zs(1,j,k)
     $                 - ca * wos2(jlt) * vwrk(m2,jlat) * alf(j,jlt) *m
     $                 + cd * wos2(jlt) * uwrk(m1,jlat) * dlf(j,jlt)
                  zs(2,j,k) = zs(2,j,k)
     $                 + ca * wos2(jlt) * vwrk(m1,jlat) * alf(j,jlt) *m
     $                 + cd * wos2(jlt) * uwrk(m2,jlat) * dlf(j,jlt)
                  ds(1,j,k) = ds(1,j,k)
     $                 - ca * wos2(jlt) * uwrk(m2,jlat) * alf(j,jlt) *m
     $                 - cd * wos2(jlt) * vwrk(m1,jlat) * dlf(j,jlt)
                  ds(2,j,k) = ds(2,j,k)
     $                 + ca * wos2(jlt) * uwrk(m1,jlat) * alf(j,jlt) *m
     $                 - cd * wos2(jlt) * vwrk(m2,jlat) * dlf(j,jlt)
               end do
            end do
         end do
      end do
c
      return
      end subroutine lgvgts2


      subroutine lgvstg (zs, ds, ug, vg, nlev, a, iwork, rwork, work)


      use precision

      implicit none
      integer iv, ifft, ms, mcoe, nlat
      real(r8) a
      integer nlev
      integer nlon
c
c  lgvstg transforms the spherical harmonic coefficients of vorticity
c  and divergence to the grid point values of u and v velocities
c
c-----------------------------------------------------------------------
C
C  On entry:
C     ZS    - vorticity spectral coefficients
C     DS    - divergence spectral coefficients
C     NLEV  - the number of vertical levels
C     A     - the radius of the sphere
C     WORK  - work space for the Fourier transform, length:3*(NLON+2)*NLAT
C
C  On exit:
C     UG    - zonal winds
C     VG    - meridional winds
C
C-----------------------------------------------------------------------
C
      real(r8) ug(*), vg(*), work(*), zs(*), ds(*)
c
      integer iwork(*)
      real(r8) rwork(*)

      integer                   ! pointers into workspaces
     $     pjmin, pjmax, pifax, palf, pdlf, psinc, 
     $     pgwt, pwos2, ptrigs
C
C-----------------------------------------------------------------------
c+
c Set the pointers from the workspace
c-
      nlon  = iwork( 1)
      nlat  = iwork( 2)
      ms    = iwork( 4)
      mcoe  = iwork( 6)
      pjmin = iwork( 7)
      pjmax = iwork( 8)
      pifax = iwork( 9)
      palf  = iwork(10)
      pdlf  = iwork(11)
      psinc = iwork(14)
      pgwt  = iwork(15)
      pwos2 = iwork(16)
      ptrigs= iwork(17)
c+
c do the transform
c-
      iv   = (nlon+2)*nlat+1
      ifft = (nlon+2)*nlat+1 + iv
      call lgvstg2( zs, ds, ug, vg, nlev, a, nlon, nlat, ms, mcoe,
     $     work(1), work(iv), work(ifft), iwork(pifax), rwork(ptrigs),
     $     iwork(pjmin), iwork(pjmax), rwork(palf),  rwork(pdlf),  
     $     rwork(psinc), rwork(pgwt),  rwork(pwos2))
c
      return
      end subroutine lgvstg


      subroutine lgvstg2( zs, ds, ug, vg, nlev, a, nlon, nlat, ms, mcoe,
     $     uwrk, vwrk, work, ifax, trigs,
     $     jmin, jmax, alf, dlf, 
     $     sinc, gwt, wos2)


      use precision

      implicit none
      integer i, init, j, lats, jlat, jlt, k, l, m, m1, m2, n
      integer ifax(*)
      real(r8) trigs(*)
      real(r8) a
      integer nlev
      real(r8) ca, cd
c
c  lgvstg2 transforms the spherical harmonic coefficients of vorticity
c  and divergence to the grid values of the u and v velocity components.
c
c-----------------------------------------------------------------------
      integer
     $     nlon,                ! number of distinct longitudes on grid
     $     nlat,                ! number of latitudes on grid
     $     ms,                  ! maximum polynomial order (for m=0)
     $     mcoe                 ! number of spherical harmonic coefficients

      integer
     $     jmin(0:ms),          ! pointer to min zonal wave number for l=n-m
     $     jmax(0:ms)           ! pointer to max zonal wave number for l=n-m

      real(r8)
     $     alf(mcoe,nlat/2),    ! associated legendre functions
     $     dlf(mcoe,nlat/2),    ! derivatives of alf
     $     sinc(nlat),          ! sin(colatitude)
     $     gwt(nlat),           ! gauss weights
     $     wos2(nlat)           ! gauss weights / sin(colat)**2
C
      real(r8)
     $     ug(nlon,nlat,nlev),          ! zonal wind grid array
     $     vg(nlon,nlat,nlev),          ! meridional wind grid array
     $     zs(2,mcoe,nlev),             ! vorticity coeffs
     $     ds(2,mcoe,nlev),             ! divergence coeffs
     $     uwrk(nlon+2,nlat),           ! grid work space for u fft
     $     vwrk(nlon+2,nlat),           ! grid work space for v fft
     $     work(nlon+2,nlat)            ! more work space

      integer isign

C
C-----------------------------------------------------------------------
c+
c loop over outer dimension of variable
c-
      do k = 1, nlev
c+
c zero the fourier coeficient accumulators
c-
         do j = 1, nlat
            do i = 1, nlon+2
               uwrk(i,j) = 0.0
               vwrk(i,j) = 0.0
            end do
         end do
c+
c loop over latitude. Input/output grid data is south to north. 
c Legendre functions are for NH from north to south.
c-
         do lats=1,nlat
            jlat = nlat+1-lats

            if( jlat.le.nlat/2 ) then
               jlt = jlat
            else
               jlt = nlat +1 -jlat
            end if
c+
c loop over all meridional modes, accumulating fourier coefficients 
c for this latitude
c-
            do l = 0, ms
               if( l.eq.0 ) then
                  m = 0
                  init = jmin(l)+1
               else
                  m=-1
                  init = jmin(l)
               end if
               ca = 1.0
               cd = 1.0
               if( jlat.gt.nlat/2 ) then
                  if( mod(l,2).eq.0 ) then
                     cd = -1.0
                  else
                     ca = -1.0
                  end if
               end if
               do j = init, jmax(l)
                  m  = m + 1
                  m1 = 2*m +1
                  m2 = 2*m +2
                  n = m + l
                  uwrk(m1,jlat) = uwrk(m1,jlat)
     $                 +(ca * alf(j,jlt) * ds(2,j,k)*m
     $                 + cd * dlf(j,jlt) * zs(1,j,k) )
     $                 / ( n*(n+1) )
                  uwrk(m2,jlat) = uwrk(m2,jlat)
     $                 -(ca * alf(j,jlt) * ds(1,j,k)*m
     $                 - cd * dlf(j,jlt) * zs(2,j,k) )
     $                 / ( n*(n+1) )
                  vwrk(m1,jlat) = vwrk(m1,jlat)
     $                 +(ca * alf(j,jlt) * zs(2,j,k)*m
     $                 - cd * dlf(j,jlt) * ds(1,j,k) )
     $                 / ( n*(n+1) )
                  vwrk(m2,jlat) = vwrk(m2,jlat)
     $                 -(ca * alf(j,jlt) * zs(1,j,k)*m
     $                 + cd * dlf(j,jlt) * ds(2,j,k) )
     $                 / ( n*(n+1) )
               end do
            end do
         end do
c+
c fourier transform all latitudes to grid values, still in gwrk 
c-
         isign = 1
         call fft991( uwrk, work, trigs, ifax, 1, nlon+2,
     $        nlon, nlat, isign)
         call fft991( vwrk, work, trigs, ifax, 1, nlon+2,
     $        nlon, nlat, isign)
c+
c copy and scale the wind components by radius of sphere * cos(lat)
c-
         do j = 1, nlat
            do i = 1, nlon
               ug(i,j,k) = uwrk(i,j) * a / sinc(j)
               vg(i,j,k) = vwrk(i,j) * a / sinc(j)
            end do
         end do
      end do
c
      return
      end subroutine lgvstg2



      subroutine lgdgts( ug, vg, ds, nlev, a, iwork, rwork, work)
C
C  LGDGTS computes the spectral coefficients of the divergence of a
C  pair of flux components in the longitudinal and latitudinal
C  directions. The input components are in standard dimensional form
C  (without cosine weighting). The input flux may be the u, v
C  horizontal velocity components, in which case the output is the
C  horizontal divergence coefficients.
C
C-----------------------------------------------------------------------
C
C  On entry:
C     UG    - zonal winds
C     VG    - meridional winds
C     NLEV  - the number of vertical levels
C     A     - the radius of the sphere
C     WORK  - work space for the Fourier transform, length:3*(NLON+2)*NLAT
C
C  On exit:
C     DS    - divergence spectral coefficients

C
C-----------------------------------------------------------------------
C

      use precision

      implicit none
      integer ms, mcoe, iv, ifft, nlat, nlon
      integer nlev
      real(r8) a
      real(r8) ug(*), vg(*), work(*), ds(*)
c
      integer iwork(*)
      real(r8) rwork(*)

      integer                   ! pointers into workspaces
     $     pjmin, pjmax, pifax, palf, pdlf, psinc, 
     $     pgwt, pwos2, ptrigs
C
C-----------------------------------------------------------------------
c+
c Set the pointers from the workspace
c-
      nlon  = iwork( 1)
      nlat  = iwork( 2)
      ms    = iwork( 4)
      mcoe  = iwork( 6)
      pjmin = iwork( 7)
      pjmax = iwork( 8)
      pifax = iwork( 9)
      palf  = iwork(10)
      pdlf  = iwork(11)
      psinc = iwork(14)
      pgwt  = iwork(15)
      pwos2 = iwork(16)
      ptrigs= iwork(17)
c+
c do the transform
c-
      iv   = (nlon+2)*nlat+1
      ifft = (nlon+2)*nlat+1 + iv
      call lgdgts2( ug, vg, ds, nlev, a, nlon, nlat, ms, mcoe,
     $     work(1), work(iv), work(ifft), iwork(pifax), rwork(ptrigs),
     $     iwork(pjmin), iwork(pjmax), rwork(palf),  rwork(pdlf),  
     $     rwork(psinc), rwork(pgwt),  rwork(pwos2))
c
      return
      end subroutine lgdgts


      subroutine lgdgts2( ug, vg, ds, nlev, a, nlon, nlat, ms, mcoe,
     $     uwrk, vwrk, work, ifax, trigs,
     $     jmin, jmax, alf, dlf, 
     $     sinc, gwt, wos2)


      use precision

      implicit none
      integer m, m1, m2
      integer i, j, jlt, jlat, k, l, lats
      integer ifax(*)
      real(r8) trigs(*)
      real(r8) a
      integer nlev
      real(r8) ca, cd
c
c
c  lgdgts computes the spectral coefficients of the divergence of a
c  pair of flux components in the longitudinal and latitudinal
c  directions.
c
c-----------------------------------------------------------------------
      integer
     $     nlon,                ! number of distinct longitudes on grid
     $     nlat,                ! number of latitudes on grid
     $     ms,                  ! maximum polynomial order (for m=0)
     $     mcoe                 ! number of spherical harmonic coefficients

      integer
     $     jmin(0:ms),          ! pointer to min zonal wave number for l=n-m
     $     jmax(0:ms)           ! pointer to max zonal wave number for l=n-m

      real(r8)
     $     alf(mcoe,nlat/2),    ! associated legendre functions
     $     dlf(mcoe,nlat/2),    ! derivatives of alf
     $     sinc(nlat),          ! sin(colatitude)
     $     gwt(nlat),           ! gauss weights
     $     wos2(nlat)           ! gauss weights / sin(colat)**2
C
      real(r8)
     $     ug(nlon,nlat,nlev),          ! zonal wind grid array
     $     vg(nlon,nlat,nlev),          ! meridional wind grid array
     $     ds(2,mcoe,nlev),             ! divergence coeffs
     $     uwrk(nlon+2,nlat),           ! grid work space for u fft
     $     vwrk(nlon+2,nlat),           ! grid work space for v fft
     $     work(nlon+2,nlat)            ! more work space

      integer isign
C
C-----------------------------------------------------------------------
c+
c zero the coefficient array
c-
      do k = 1, nlev
         do j = 1, mcoe
            ds(1,j,k) = 0.
            ds(2,j,k) = 0.
         end do
      end do
c+
c loop over outer dimension of variable
c-
      do k = 1, nlev
c+
c copy input grid array to grid work space with 2 extra longitude points
c weight wind components by cos(latitude) / radius of sphere
c-
         do j = 1, nlat
            do i = 1, nlon
               uwrk(i,j) = ug(i,j,k) * sinc(j) / a
               vwrk(i,j) = vg(i,j,k) * sinc(j) / a
            end do
         end do
c+
c fourier transform all latitudes, result still in gwrk
c-
         isign = -1
         call fft991 (uwrk, work, trigs, ifax, 1, nlon+2,
     $        nlon, nlat, isign)
         call fft991 (vwrk, work, trigs, ifax, 1, nlon+2,
     $        nlon, nlat, isign)
c+
c loop over latitude. Input/output grid data is south to north. 
c Legendre functions are for NH from north to south.
c-
         do lats=1,nlat
            jlat = nlat+1-lats
            if( jlat.le.nlat/2 ) then
               jlt = jlat
            else
               jlt = nlat +1 -jlat
            end if
            do l = 0, ms
               m = -1
               ca = 1.0
               cd = 1.0
               if( jlat.gt.nlat/2 ) then
                  if( mod(l,2).eq.0 ) then
                     cd = -1.0
                  else
                     ca = -1.0
                  end if
               end if
               do j = jmin(l), jmax(l)
                  m  = m + 1
                  m1 = 2*m +1
                  m2 = 2*m +2
                  ds(1,j,k) = ds(1,j,k)
     $                 - ca * wos2(jlt) * uwrk(m2,jlat) * alf(j,jlt) *m
     $                 - cd * wos2(jlt) * vwrk(m1,jlat) * dlf(j,jlt)
                  ds(2,j,k) = ds(2,j,k)
     $                 + ca * wos2(jlt) * uwrk(m1,jlat) * alf(j,jlt) *m
     $                 - cd * wos2(jlt) * vwrk(m2,jlat) * dlf(j,jlt)
               end do
            end do
         end do
      end do
c
      return
      end  subroutine lgdgts2



      subroutine lggstg( sp, gi, gj, nlev, a, iwork, rwork, work )


      use precision

      implicit none
      integer iv
      integer ifft
      integer mcoe
      integer ms
      integer nlat
      integer nlon
      integer nlev
      real(r8) a
C
C  LGGSTG transforms the spherical harmonic coefficients of a field to
C  grid values of the longitudal and latitudal components of its gradient.
C
C-----------------------------------------------------------------------
C
C  On entry:
C     SP    - spectral coefficients
C     NLEV  - the number of vertical levels
C     A     - the radius of the sphere
C     WORK  - work space for the Fourier transform, length:3*(NLON+2)*NLAT
C
C  On exit:
C     GI    - longitudinal derivative grid values
C                   1           d
C                _________   _______   F
C                a cos_phi   d_lamda
C     GJ    - latitudinal derivative grid values
C                1     d
C                _   _____   F
C                a   d_phi
C
C-----------------------------------------------------------------------
C
      real(r8) gi(*), gj(*), work(*), sp(*)
c
      integer iwork(*)
      real(r8) rwork(*)

      integer                   ! pointers into workspaces
     $     pjmin, pjmax, pifax, palf, pdlf, psinc, 
     $     pgwt, pwos2, ptrigs
C
C-----------------------------------------------------------------------
c+
c Set the pointers from the workspace
c-
      nlon  = iwork( 1)
      nlat  = iwork( 2)
      ms    = iwork( 4)
      mcoe  = iwork( 6)
      pjmin = iwork( 7)
      pjmax = iwork( 8)
      pifax = iwork( 9)
      palf  = iwork(10)
      pdlf  = iwork(11)
      psinc = iwork(14)
      pgwt  = iwork(15)
      pwos2 = iwork(16)
      ptrigs= iwork(17)
c+
c do the transform
c-
      iv   = (nlon+2)*nlat+1
      ifft = (nlon+2)*nlat+1 + iv
      call lggstg2( sp, gi, gj, nlev, a, nlon, nlat, ms, mcoe,
     $     work(1), work(iv), work(ifft), iwork(pifax), rwork(ptrigs),
     $     iwork(pjmin), iwork(pjmax), rwork(palf),  rwork(pdlf),  
     $     rwork(psinc), rwork(pgwt),  rwork(pwos2))
c
      return
      end subroutine lggstg


      subroutine lggstg2( sp, gi, gj, nlev, a, nlon, nlat, ms, mcoe,
     $     uwrk, vwrk, work, ifax, trigs,
     $     jmin, jmax, alf, dlf, 
     $     sinc, gwt, wos2)


      use precision

      implicit none
      integer i, j, k, lats, jlat, jlt, l
      integer m, m1, m2
      real(r8) ca, cd
      real(r8) trigs(*)
      integer ifax(*)
      integer nlev
      real(r8) a
c
c  lggstg2 transforms the spherical harmonic coefficients of a field to
c  the grid values of its gradient.
c
c-----------------------------------------------------------------------
      integer
     $     nlon,                ! number of distinct longitudes on grid
     $     nlat,                ! number of latitudes on grid
     $     ms,                  ! maximum polynomial order (for m=0)
     $     mcoe                 ! number of spherical harmonic coefficients

      integer
     $     jmin(0:ms),          ! pointer to min zonal wave number for l=n-m
     $     jmax(0:ms)           ! pointer to max zonal wave number for l=n-m

      real(r8)
     $     alf(mcoe,nlat/2),    ! associated legendre functions
     $     dlf(mcoe,nlat/2),    ! derivatives of alf
     $     sinc(nlat),          ! sin(colatitude)
     $     gwt(nlat),           ! gauss weights
     $     wos2(nlat)           ! gauss weights / sin(colat)**2
C
      real(r8)
     $     gi(nlon,nlat,nlev),          ! zonal gradient grid array
     $     gj(nlon,nlat,nlev),          ! meridional gradient grid array
     $     sp(2,mcoe,nlev),             ! spectral coefficients
     $     uwrk(nlon+2,nlat),           ! grid work space for u fft
     $     vwrk(nlon+2,nlat),           ! grid work space for v fft
     $     work(nlon+2,nlat)            ! more work space

      integer isign

C
C-----------------------------------------------------------------------
c+
c loop over outer dimension of variable
c-
      do k = 1, nlev
c+
c zero the fourier coeficient accumulators
c-
         do j = 1, nlat
            do i = 1, nlon+2
               uwrk(i,j) = 0.0
               vwrk(i,j) = 0.0
            end do
         end do
c+
c loop over latitude. Input/output grid data is south to north. 
c Legendre functions are for NH from north to south.
c-
         do lats=1,nlat
            jlat = nlat+1-lats 

            if( jlat.le.nlat/2 ) then
               jlt = jlat
            else
               jlt = nlat +1 -jlat
            end if
c+
c loop over all meridional modes, accumulating fourier coefficients 
c for this latitude
c-
            do l = 0, ms
               m=-1
               ca = 1.0
               cd = 1.0
               if( jlat.gt.nlat/2 ) then
                  if( mod(l,2).eq.0 ) then
                     cd = -1.0
                  else
                     ca = -1.0
                  end if
               end if
               do j = jmin(l), jmax(l)
                  m  = m + 1
                  m1 = 2*m +1
                  m2 = 2*m +2
                  uwrk(m1,jlat) = uwrk(m1,jlat)
     $                 - ca * alf(j,jlt) * sp(2,j,k) * m
                  uwrk(m2,jlat) = uwrk(m2,jlat)
     $                 + ca * alf(j,jlt) * sp(1,j,k) * m
                  vwrk(m1,jlat) = vwrk(m1,jlat)
     $                 + cd * dlf(j,jlt) * sp(1,j,k)
                  vwrk(m2,jlat) = vwrk(m2,jlat)
     $                 + cd * dlf(j,jlt) * sp(2,j,k)
               end do
            end do
         end do
c+
c fourier transform all latitudes to grid values, still in gwrk 
c-
         isign = 1
         call fft991( uwrk, work, trigs, ifax, 1, nlon+2,
     $        nlon, nlat, isign)
         call fft991( vwrk, work, trigs, ifax, 1, nlon+2,
     $        nlon, nlat, isign)
c+
c copy and scale the gradient components
c-
         do j = 1, nlat
            do i = 1, nlon
               gi(i,j,k) = uwrk(i,j) / (a * sinc(j))
               gj(i,j,k) = vwrk(i,j) / (a * sinc(j))
            end do
         end do
      end do
c
      return
      end  subroutine lggstg2

      end module lgtran
